import glob
import logging
import os
import re
import sys
import csv
import yaml
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from argparse import ArgumentParser
from pathlib import Path

import numpy as np
import torch
import uproot
import awkward as ak
import matplotlib.pyplot as plt
from sklearn.metrics import auc, roc_curve
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler
from tqdm import tqdm

from pytorch_lightning import Trainer
from pytorch_lightning.callbacks import EarlyStopping, LearningRateMonitor, ModelCheckpoint
from pytorch_lightning.loggers import TensorBoardLogger

from dataloader import get_dataloader
from model import ABCDLightningModule

MISSING_VALUE = -999.0

def load_config(config_path):
    with open(config_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def resolve_paths(base, paths):
    if isinstance(paths, str):
        paths = [paths]
    resolved = []
    for p in paths:
        full = os.path.expanduser(os.path.join(base, p) if base else p)
        matches = sorted(glob.glob(full))
        if matches:
            resolved.extend(matches)
        else:
            resolved.append(full)
    return resolved

def parse_training_features(training_features_cfg):
    training_features = []
    feature_transforms = {}

    for feat in training_features_cfg:
        if isinstance(feat, dict):
            for name, tf in feat.items():
                training_features.append(name)
                feature_transforms[name] = tf
        else:
            training_features.append(feat)
            feature_transforms[feat] = "none"

    return training_features, feature_transforms


def _extract_expression_variables(expression):
    tokens = re.findall(r"(?<!\.)\b[A-Za-z_][A-Za-z0-9_]*\b", expression)
    excluded = {"np", "abs", "and", "or", "not", "True", "False"}
    return [tok for tok in tokens if tok not in excluded]


def _normalize_derived_vars_cfg(derived_vars_cfg):
    if derived_vars_cfg is None:
        return {}
    if isinstance(derived_vars_cfg, dict):
        return dict(derived_vars_cfg)
    if isinstance(derived_vars_cfg, list):
        normalized = {}
        for item in derived_vars_cfg:
            if not isinstance(item, dict):
                raise ValueError("Each item in 'derived_vars' list must be a mapping of 'new_var: expression'.")
            for key, value in item.items():
                normalized[key] = value
        return normalized
    raise ValueError("'derived_vars' must be either a mapping or a list of single-item mappings.")


def _collect_derived_input_vars(derived_vars_cfg):
    derived_vars_cfg = _normalize_derived_vars_cfg(derived_vars_cfg)

    derived_names = set(derived_vars_cfg.keys())
    needed = []
    for expr in derived_vars_cfg.values():
        if not isinstance(expr, str):
            continue
        for var in _extract_expression_variables(expr):
            if var not in derived_names:
                needed.append(var)
    return list(dict.fromkeys(needed))


def apply_derived_vars(data, derived_vars_cfg):
    derived_vars_cfg = _normalize_derived_vars_cfg(derived_vars_cfg)
    if not derived_vars_cfg:
        return data

    out = {k: np.asarray(v).copy() for k, v in data.items()}
    n_events = _data_length(out)

    for new_var, expression in derived_vars_cfg.items():
        if not isinstance(expression, str):
            raise ValueError(f"derived_vars['{new_var}'] must be a string expression")

        local_dict = {k: np.asarray(v) for k, v in out.items()}
        try:
            values = np.asarray(eval(expression, {"__builtins__": {}, "np": np, "abs": np.abs}, local_dict))
        except NameError as exc:
            raise ValueError(
                f"Failed to evaluate derived var '{new_var}': missing variable in expression '{expression}'"
            ) from exc
        except Exception as exc:
            raise ValueError(
                f"Failed to evaluate derived var '{new_var}' with expression '{expression}': {exc}"
            ) from exc

        if values.ndim == 0:
            values = np.full(n_events, values.item())

        if len(values) != n_events:
            raise ValueError(
                f"Derived var '{new_var}' has length {len(values)} but expected {n_events}"
            )

        out[new_var] = values
        logging.info("Computed derived var '%s' from config expression", new_var)

    return out

def _read_root_frame(path, branches, dataset_idx, sample_idx):
    with uproot.open(path) as root_file:
        arrays = root_file["Events"].arrays(branches, library="np")

    if not isinstance(arrays, dict):
        arrays = {name: arrays[name] for name in arrays.dtype.names}

    columns = {}
    n_events = None
    for branch in branches:
        if branch not in arrays:
            continue

        values = np.asarray(arrays[branch])

        if values.ndim == 1 and values.dtype != object:
            clean_values = values
        else:
            flat_values = [_flatten_awkward_cell(v) for v in values]
            clean_values = np.asarray(flat_values)

        if n_events is None:
            n_events = len(clean_values)
        else:
            n_events = min(n_events, len(clean_values))

        columns[branch] = clean_values

    if n_events is None or n_events == 0:
        return {
            "dataset_idx": np.array([], dtype=np.int32),
            "sample_idx": np.array([], dtype=np.int32),
        }

    trimmed = {name: np.asarray(vals)[:n_events] for name, vals in columns.items()}
    trimmed["dataset_idx"] = np.full(n_events, dataset_idx, dtype=np.int32)
    trimmed["sample_idx"] = np.full(n_events, sample_idx, dtype=np.int32)
    return trimmed


def _flatten_awkward_cell(value):
    if np.isscalar(value):
        return value
    return value[0] if len(value) > 0 else MISSING_VALUE


def _data_length(data):
    if not data:
        return 0
    first_key = next(iter(data))
    return len(data[first_key])


def _apply_mask(data, mask):
    return {k: np.asarray(v)[mask] for k, v in data.items()}


def _concat_chunks(chunks):
    if not chunks:
        return {}

    all_keys = set()
    for chunk in chunks:
        all_keys.update(chunk.keys())

    combined = {}
    for key in all_keys:
        arrays = [np.asarray(chunk[key]) for chunk in chunks if key in chunk]
        combined[key] = np.concatenate(arrays) if arrays else np.array([])
    return combined


def _evaluate_preselection_mask(data, expression):
    local_dict = {k: np.asarray(v) for k, v in data.items()}
    return np.asarray(eval(expression, {"__builtins__": {}, "np": np}, local_dict), dtype=bool)


# Also return fitted min/max so they can be saved for inference-time scaling
def _safe_minmax_scale(values, valid_mask):
    scaled = np.zeros_like(values, dtype=np.float64)
    scaler = MinMaxScaler()
    scaler.fit(values[valid_mask].reshape(-1, 1))
    scaled[valid_mask] = scaler.transform(values[valid_mask].reshape(-1, 1)).ravel()
    return scaled, float(scaler.data_min_[0]), float(scaler.data_max_[0])


def load_data(paths, features, extra_vars, num_workers=1):
    branches = list(dict.fromkeys(features + extra_vars))

    sample_names = []
    for path in paths:
        p = Path(path)
        parent = p.parent.name
        grandparent = p.parent.parent.name if p.parent.parent is not None else ""
        sample_names.append(grandparent if parent.isdigit() and grandparent else parent)

    sample_name_to_idx = {name: idx for idx, name in enumerate(sorted(set(sample_names)))}
    indexed_paths = [
        (dataset_idx, path, sample_name_to_idx[sample_name])
        for dataset_idx, (path, sample_name) in enumerate(zip(paths, sample_names))
    ]
    chunks = [None] * len(indexed_paths)

    if num_workers > 1 and len(indexed_paths) > 1:
        max_workers = min(num_workers, len(indexed_paths))
        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            futures = {
                pool.submit(_read_root_frame, path, branches, dataset_idx, sample_idx): dataset_idx
                for dataset_idx, path, sample_idx in indexed_paths
            }
            for future in tqdm(as_completed(futures), total=len(futures), desc="Loading ROOT files"):
                dataset_idx = futures[future]
                chunks[dataset_idx] = future.result()
    else:
        for dataset_idx, path, sample_idx in tqdm(indexed_paths, total=len(indexed_paths), desc="Loading ROOT files"):
            chunks[dataset_idx] = _read_root_frame(path, branches, dataset_idx, sample_idx)

    chunks = [chunk for chunk in chunks if chunk is not None]

    data = _concat_chunks(chunks)
    logging.info("Loaded %d events from %d files.", _data_length(data), len(paths))
    if "weight" in data:
        data = _apply_mask(data, np.asarray(data["weight"]) > 0)

    return data


# If scaler_output_path is provided, fitted min/max per feature are saved to JSON for use at inference time
def preprocess_data(data, training_features, feature_transforms, constraint_var, scaler_output_path=None):
    out = {k: np.asarray(v).copy() for k, v in data.items()}
    cols_to_clean = list(dict.fromkeys(training_features + [constraint_var]))
    present_cols = [col for col in cols_to_clean if col in out]
    for col in present_cols:
        arr = np.asarray(out[col])
        if arr.dtype == object:
            arr = np.asarray([_flatten_awkward_cell(v) for v in arr], dtype=np.float64)
        else:
            arr = arr.astype(np.float64, copy=False)
        arr = np.where(np.isfinite(arr), arr, np.nan)
        out[col] = np.nan_to_num(arr, nan=MISSING_VALUE, posinf=MISSING_VALUE, neginf=MISSING_VALUE)
    scaler_params = {}
    for feat in training_features:
        feat_arr = np.asarray(out[feat], dtype=np.float64)
        valid = (feat_arr != MISSING_VALUE) & np.isfinite(feat_arr)
        transform = feature_transforms.get(feat, "none")
        if transform == "log":
            positive = valid & (feat_arr > 0)
            feat_arr[valid & ~positive] = MISSING_VALUE
            valid = positive
        if transform == "log":
            feat_arr[valid] = np.log(feat_arr[valid])
        if valid.sum() == 0:
            print(f"WARNING: feature '{feat}' has no valid samples, skipping scaler fit")
            out[feat] = np.zeros_like(feat_arr, dtype=np.float64)
            scaler_params[feat] = {"transform": transform, "min": 0.0, "max": 1.0}
            continue
        out[feat], fmin, fmax = _safe_minmax_scale(feat_arr, valid)
        scaler_params[feat] = {"transform": transform, "min": fmin, "max": fmax}
    if constraint_var not in training_features:
        constraint_arr = np.asarray(out[constraint_var], dtype=np.float64)
        valid_constraint = constraint_arr != MISSING_VALUE
        out[constraint_var], fmin, fmax = _safe_minmax_scale(constraint_arr, valid_constraint)
        scaler_params[constraint_var] = {"transform": "none", "min": fmin, "max": fmax}
    if scaler_output_path is not None:
        import json
        scaler_out = {
            "_training_features": training_features,
            "_constraint_var": constraint_var,
        }
        scaler_out.update(scaler_params)
        with open(scaler_output_path, "w") as f:
            json.dump(scaler_out, f, indent=2)
        logging.info("Saved scaler params to %s", scaler_output_path)
    return out

def normalize_class_weights(sig_data, bkg_data):
    sig = {k: np.asarray(v).copy() for k, v in sig_data.items()}
    bkg = {k: np.asarray(v).copy() for k, v in bkg_data.items()}

    sig_weight_sum = np.sum(sig["weight"])
    bkg_weight_sum = np.sum(bkg["weight"])

    sig["weight"] = sig["weight"] / sig_weight_sum
    bkg["weight"] = bkg["weight"] / bkg_weight_sum
    return sig, bkg


def make_dataloaders(data, training_features, constraint_var, batch_size):
    labels_str = np.asarray(data["label"]).astype(np.int32).astype(str)
    stratify_col = "sample_idx" if "sample_idx" in data else "dataset_idx"
    sample_str = np.asarray(data[stratify_col]).astype(np.int64).astype(str)
    stratify_key = np.char.add(np.char.add(labels_str, "_"), sample_str)
    logging.info("Using stratification column '%s' for train/val split", stratify_col)

    unique_keys, unique_counts = np.unique(stratify_key, return_counts=True)
    rare_keys = unique_keys[unique_counts < 2]
    if len(rare_keys) > 0:
        rare_mask = np.isin(stratify_key, rare_keys)
        stratify_key[rare_mask] = np.char.add(labels_str[rare_mask], "_rare")
        logging.info(
            "Collapsed %d rare strata (<2 events) into label-level rare bins for stable splitting",
            len(rare_keys),
        )

    all_indices = np.arange(_data_length(data))

    train_idx, val_idx = train_test_split(
        all_indices,
        test_size=0.2,
        random_state=42,
        stratify=stratify_key,
    )

    feature_matrix = np.column_stack([data[f] for f in training_features]).astype(np.float32, copy=False)
    constraint_values = np.asarray(data[constraint_var], dtype=np.float32).reshape(-1, 1)
    labels = np.asarray(data["label"], dtype=np.float32)
    weights = np.asarray(data["weight"], dtype=np.float32)

    train_loader = get_dataloader(
        dnn_input_data=torch.from_numpy(feature_matrix[train_idx]),
        constraint_data=torch.from_numpy(constraint_values[train_idx]),
        labels=torch.from_numpy(labels[train_idx]),
        weights=torch.from_numpy(weights[train_idx]),
        batch_size=batch_size,
        use_sampler=False,
    )

    val_loader = get_dataloader(
        dnn_input_data=torch.from_numpy(feature_matrix[val_idx]),
        constraint_data=torch.from_numpy(constraint_values[val_idx]),
        labels=torch.from_numpy(labels[val_idx]),
        weights=torch.from_numpy(weights[val_idx]),
        batch_size=batch_size,
        use_sampler=False,
        is_validation=True,
    )

    return train_loader, val_loader

def run_training(
    model,
    train_loader,
    val_loader,
    output_dir,
    max_epochs,
    flavor,
    devices=[0],
    check_val_every_n_epoch=1,
    early_stopping_patience=20,
    early_stopping_min_delta=1e-4,
):
    run_dir = Path(output_dir)
    run_dir.mkdir(parents=True, exist_ok=True)

    logger = TensorBoardLogger(save_dir=str(run_dir), name=flavor)

    checkpoint = ModelCheckpoint(
        dirpath=str(Path(logger.log_dir) / "checkpoints"),
        filename=f"{model.flavor}-abcdisco-{{epoch:03d}}-{{val_loss:.4f}}",
        monitor="val_loss",
        mode="min",
        save_top_k=5,
        save_weights_only=False,
    )
    lr_monitor = LearningRateMonitor(logging_interval="epoch")

    callbacks = [checkpoint, lr_monitor]
    if early_stopping_patience is not None and early_stopping_patience > 0:
        callbacks.append(
            EarlyStopping(
                monitor="val_loss",
                mode="min",
                patience=early_stopping_patience,
                min_delta=early_stopping_min_delta,
            )
        )

    use_gpu = torch.cuda.is_available()
    trainer = Trainer(
        max_epochs=max_epochs,
        accelerator="gpu" if use_gpu else "cpu",
        devices=devices if use_gpu else None,
        logger=logger,
        check_val_every_n_epoch=check_val_every_n_epoch,
        callbacks=callbacks,
    )

    trainer.fit(model, train_dataloaders=train_loader, val_dataloaders=val_loader)
    return trainer


def _concat_sig_bkg(sig_data, bkg_data):
    column_order = list(sig_data.keys())
    for key in bkg_data.keys():
        if key not in sig_data:
            column_order.append(key)

    combined = {}
    for key in column_order:
        sig_col = np.asarray(sig_data.get(key, np.array([])))
        bkg_col = np.asarray(bkg_data.get(key, np.array([])))
        combined[key] = np.concatenate([sig_col, bkg_col])

    n_sig = _data_length(sig_data)
    n_bkg = _data_length(bkg_data)
    combined["label"] = np.concatenate(
        [
            np.ones(n_sig, dtype=np.int32),
            np.zeros(n_bkg, dtype=np.int32),
        ]
    )
    return combined


def _write_csv_numpy(data, output_csv):
    columns = list(data.keys())
    arrays = [np.asarray(data[col]) for col in columns]

    if not arrays:
        Path(output_csv).write_text("\n", encoding="utf-8")
        return

    lengths = {arr.shape[0] for arr in arrays}
    if len(lengths) != 1:
        raise ValueError(f"All output columns must have the same length, got lengths={sorted(lengths)}")

    all_numeric = all(np.issubdtype(arr.dtype, np.number) for arr in arrays)
    if all_numeric:
        matrix = np.column_stack(arrays)
        np.savetxt(
            output_csv,
            matrix,
            delimiter=",",
            header=",".join(columns),
            comments="",
            fmt="%.10g",
        )
        return

    def _format_cell(val):
        if isinstance(val, (np.floating, float)):
            return f"{float(val):.10g}"
        if isinstance(val, (np.integer, int)):
            return str(int(val))
        if isinstance(val, (np.bool_, bool)):
            return str(int(val))
        if isinstance(val, bytes):
            return val.decode("utf-8")
        return str(val)

    output_path = Path(output_csv)
    with output_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(columns)
        for row in zip(*arrays):
            writer.writerow([_format_cell(v) for v in row])


def _latest_version_dir(output_dir: Path, flavor) -> Path:
    logs_dir = output_dir / flavor
    version_dirs = [p for p in logs_dir.glob("version_*") if p.is_dir()]
    if not version_dirs:
        raise FileNotFoundError(f"No version_* directories found under {logs_dir}")

    def _version_key(path: Path):
        suffix = path.name.replace("version_", "")
        try:
            return (0, int(suffix))
        except ValueError:
            return (1, path.stat().st_mtime)

    return sorted(version_dirs, key=_version_key)[-1]


def _latest_checkpoint_from_config(cfg, flavor) -> Path:
    output_dir = Path(cfg.get("output", "simple_abcdisco_output"))
    version_dir = _latest_version_dir(output_dir, flavor=flavor)
    ckpt_dir = version_dir / "checkpoints"
    ckpts = [p for p in ckpt_dir.glob("*.ckpt") if p.is_file()]
    if not ckpts:
        raise FileNotFoundError(f"No .ckpt files found under {ckpt_dir}")
    return sorted(ckpts, key=lambda p: p.stat().st_mtime)[-1]


def _inference_output_dir_from_checkpoint(checkpoint_path: Path) -> Path:
    if checkpoint_path.parent.name == "checkpoints":
        return checkpoint_path.parent.parent
    return checkpoint_path.parent


def _build_model_from_config(cfg, flavor, input_size, checkpoint_path, device):
    model = ABCDLightningModule(
        input_size=input_size,
        hidden_layers=cfg.get("architecture", [64, 32, 16]),
        learning_rate=cfg.get("learning_rate", 1e-3),
        bce_weight=cfg.get("bce_weight", 1.0),
        disco_lambda=cfg.get("disco_lambda", 0.0),
        flavor=flavor,
        use_batchnorm=cfg.get("use_batchnorm", True),
        dropout=cfg.get("dropout", 0.0),
    )

    checkpoint = torch.load(checkpoint_path, map_location=device)
    if "state_dict" in checkpoint:
        model.load_state_dict(checkpoint["state_dict"], strict=True)
    else:
        model.load_state_dict(checkpoint, strict=True)

    model.to(device)
    model.eval()
    return model


def _batched_scores(model, features_tensor, batch_size, flavor, device):
    scores_0 = []
    scores_1 = []

    with torch.no_grad():
        for start in tqdm(range(0, len(features_tensor), batch_size), desc="Inferring events"):
            batch = features_tensor[start:start + batch_size].to(device)
            logits = model(batch)
            if logits.ndim == 1:
                logits = logits.unsqueeze(-1)
            probs = torch.sigmoid(logits).detach().cpu().numpy()
            scores_0.append(probs[:, 0])
            if flavor == "double":
                scores_1.append(probs[:, 1])

    out_0 = np.concatenate(scores_0) if scores_0 else np.array([], dtype=np.float32)
    out_1 = np.concatenate(scores_1) if scores_1 else np.array([], dtype=np.float32)
    return out_0, out_1


def _plot_roc_curves(data, flavor, output_path):
    labels = np.asarray(data["label"])
    weights = np.asarray(data["weight"]) if "weight" in data else None

    plt.figure(figsize=(8, 7))

    roc_series = [
        ("DNN", "dnn_score")
    ] if flavor == "single" else [
        ("DNN 0", "dnn_0_score"),
        ("DNN 1", "dnn_1_score"),
    ]

    for label_name, score_col in roc_series:
        fpr, tpr, _ = roc_curve(labels, np.asarray(data[score_col]), sample_weight=weights)
        roc_auc = auc(fpr, tpr)
        plt.plot(fpr, tpr, linewidth=2, label=f"{label_name} (AUC={roc_auc:.4f})")

    plt.plot([0, 1], [0, 1], "k--", linewidth=1)
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(f"ROC Curve ({flavor} flavor)")
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


def _plot_score_densities(data, flavor, output_path):
    score_cols = ["dnn_score"] if flavor == "single" else ["dnn_0_score", "dnn_1_score"]
    fig, axes = plt.subplots(1, len(score_cols), figsize=(8 * len(score_cols), 6), squeeze=False)

    labels = np.asarray(data["label"])
    sig_mask = labels == 1
    bkg_mask = labels == 0

    if "weight" in data:
        weights = np.asarray(data["weight"])
        sig_w = weights[sig_mask]
        bkg_w = weights[bkg_mask]
    else:
        sig_w = None
        bkg_w = None

    for idx, col in enumerate(score_cols):
        ax = axes[0, idx]
        values = np.asarray(data[col])
        ax.hist(
            values[sig_mask],
            bins=50,
            range=(0, 1),
            weights=sig_w,
            density=True,
            histtype="step",
            linewidth=2,
            label="Signal",
            color="tab:red",
        )
        ax.hist(
            values[bkg_mask],
            bins=50,
            range=(0, 1),
            weights=bkg_w,
            density=True,
            histtype="step",
            linewidth=2,
            label="Background",
            color="tab:blue",
        )
        ax.set_xlabel(col)
        ax.set_ylabel("Density")
        ax.set_title(f"Score density: {col}")
        ax.legend(loc="best")

    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


def run_inference(args, cfg, flavor, sig_data, bkg_data, training_features, feature_transforms, constraint_var, derived_vars_cfg, is_data=False, inference_data=None):
    logging.info("Starting inference steps...")
    checkpoint_path = Path(args.checkpoint) if args.checkpoint else _latest_checkpoint_from_config(cfg, flavor=flavor)
    logging.info("Using checkpoint: %s", checkpoint_path)

    if is_data:
        full_inference_data = inference_data
    else:
        full_inference_data = _concat_sig_bkg(sig_data, bkg_data)
        
    full_inference_data = apply_derived_vars(full_inference_data, derived_vars_cfg)

    logging.info("Preprocessing inference features...")
    processed = preprocess_data(
        full_inference_data,
        training_features=training_features,
        feature_transforms=feature_transforms,
        constraint_var=constraint_var,
    )

    feature_matrix = np.column_stack([processed[f] for f in training_features]).astype(np.float32, copy=False)
    features_tensor = torch.from_numpy(feature_matrix)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = _build_model_from_config(
        cfg=cfg,
        flavor=flavor,
        input_size=len(training_features),
        checkpoint_path=checkpoint_path,
        device=device,
    )

    score_0, score_1 = _batched_scores(model, features_tensor, cfg.get("batch_size", 8192), flavor, device)

    if flavor == "single":
        full_inference_data["dnn_score"] = score_0
    else:
        full_inference_data["dnn_0_score"] = score_0
        full_inference_data["dnn_1_score"] = score_1

    if is_data:
        default_out = _inference_output_dir_from_checkpoint(checkpoint_path) / f"predictions_{flavor}_data.csv"
    else:
        default_out = _inference_output_dir_from_checkpoint(checkpoint_path) / f"predictions_{flavor}.csv"

    output_csv = Path(args.output_csv) if args.output_csv else default_out
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    logging.info("Writing inference CSV to %s", output_csv)
    _write_csv_numpy(full_inference_data, output_csv)
    logging.info("Done. Wrote %d rows.", _data_length(full_inference_data))

    if not is_data:
        roc_path = output_csv.with_name(f"{output_csv.stem}_roc.png")
        density_path = output_csv.with_name(f"{output_csv.stem}_score_density.png")
        _plot_roc_curves(full_inference_data, flavor=flavor, output_path=roc_path)
        _plot_score_densities(full_inference_data, flavor=flavor, output_path=density_path)
        logging.info("Saved ROC plot to %s", roc_path)
        logging.info("Saved score density plot to %s", density_path)


def main():
    parser = ArgumentParser()
    parser.add_argument("--config", required=True, help="Path to YAML config")
    parser.add_argument("--flavor", choices=["single", "double"], default=None, help="Training flavor: single (one output) or double (two outputs)")
    parser.add_argument("--data", action="store_true", help="Run inference data (without training) using the latest checkpoint from config")
    parser.add_argument("--infer", action="store_true", help="Skip training and run inference only")
    parser.add_argument("--checkpoint", default=None, help="Path to model checkpoint (.ckpt) for inference. If omitted, auto-picks newest checkpoint.")
    parser.add_argument("--output-csv", default=None, help="Output CSV path")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    cfg = load_config(args.config)
    flavor = args.flavor if args.flavor is not None else cfg.get("flavor", "single")
    logging.info("Using flavor=%s", flavor)

    training_features, feature_transforms = parse_training_features(cfg["training_features"])
    constraint_var = cfg["constraint_var"]
    derived_vars_cfg = _normalize_derived_vars_cfg(cfg.get("derived_vars", {}))
    if flavor == "double" and constraint_var not in training_features:
        training_features.append(constraint_var)
        feature_transforms[constraint_var] = cfg.get("constraint_as_feature_transform", "none")
        logging.info(
            "Double flavor: added constraint_var '%s' to training features as regular input",
            constraint_var,
        )

    if cfg.get("auto_include_derived_vars", False):
        auto_tf = cfg.get("auto_include_derived_vars_transform", "none")
        added_count = 0
        for derived_name in derived_vars_cfg.keys():
            if derived_name in training_features:
                continue
            training_features.append(derived_name)
            feature_transforms[derived_name] = auto_tf
            added_count += 1
        logging.info(
            "Auto-include derived vars enabled: added %d derived features with transform='%s'",
            added_count,
            auto_tf,
        )

    derived_var_names = set(derived_vars_cfg.keys())
    derived_input_vars = _collect_derived_input_vars(derived_vars_cfg)

    load_training_features = [f for f in training_features if f not in derived_var_names]
    load_constraint = [] if constraint_var in derived_var_names else [constraint_var]
    load_features = list(dict.fromkeys(load_training_features + load_constraint + derived_input_vars))
    extra_vars = cfg.get("extra_vars", [])
    if "weight" not in extra_vars:
        extra_vars.append("weight")

    sig_base = cfg.get("sig_base_path", cfg.get("input_base_path", ""))
    bkg_base = cfg.get("bkg_base_path", cfg.get("input_base_path", ""))

    sig_paths = resolve_paths(sig_base, cfg["sig_path"])
    bkg_paths = resolve_paths(bkg_base, cfg["bkg_path"])
    io_workers = int(cfg.get("io_workers", min(8, os.cpu_count() or 1)))
    logging.info("Using io_workers=%d for ROOT loading", io_workers)

    if args.data:
        if not args.infer:
            parser.error("--data can only be used with --infer")
        
        data_base = cfg.get("data_base_path", cfg.get("input_base_path", ""))
        if "data_path" not in cfg:
            parser.error("Config must contain 'data_path' when using --data")

        data_paths = resolve_paths(data_base, cfg["data_path"])
        logging.info("Loading real data...")
        real_data = load_data(data_paths, load_features, extra_vars, num_workers=io_workers)
        
        logging.info("Data samples: %d", _data_length(real_data))
        if cfg.get("preselection"):
            data_mask = _evaluate_preselection_mask(real_data, cfg["preselection"])
            real_data = _apply_mask(real_data, data_mask)
        logging.info(
            "After preselection - Data samples: %d",
            _data_length(real_data),
        )

        logging.info("Skipping training. Running inference on data only...")
        run_inference(args, cfg, flavor, None, None, training_features, feature_transforms, constraint_var, derived_vars_cfg, is_data=True, inference_data=real_data)
        return

    logging.info("Loading signal data...")
    sig_data = load_data(sig_paths, load_features, extra_vars, num_workers=io_workers)
    logging.info("Loading background data...")
    bkg_data = load_data(bkg_paths, load_features, extra_vars, num_workers=io_workers)


    logging.info("Signal samples: %d, Background samples: %d", _data_length(sig_data), _data_length(bkg_data))
    if cfg.get("preselection"):
        sig_mask = _evaluate_preselection_mask(sig_data, cfg["preselection"])
        bkg_mask = _evaluate_preselection_mask(bkg_data, cfg["preselection"])
        sig_data = _apply_mask(sig_data, sig_mask)
        bkg_data = _apply_mask(bkg_data, bkg_mask)
    logging.info(
        "After preselection - Signal samples: %d, Background samples: %d",
        _data_length(sig_data),
        _data_length(bkg_data),
    )

    sig_data["label"] = np.ones(_data_length(sig_data), dtype=np.float32)
    bkg_data["label"] = np.zeros(_data_length(bkg_data), dtype=np.float32)
    
    # Keep copies of the original data (with raw weights and features) for inference
    raw_sig_data = {k: np.copy(v) for k, v in sig_data.items()}
    raw_bkg_data = {k: np.copy(v) for k, v in bkg_data.items()}

    sig_data, bkg_data = normalize_class_weights(sig_data, bkg_data)

    logging.info("Preprocessing data...")
    combined_keys = set(sig_data.keys()) | set(bkg_data.keys())
    data = {
        key: np.concatenate([sig_data.get(key, np.array([])), bkg_data.get(key, np.array([]))])
        for key in combined_keys
    }
    data = apply_derived_vars(data, derived_vars_cfg)
    # Save scaler params alongside the output so inference can apply identical scaling
    scaler_path = Path(cfg.get("output", "simple_abcdisco_output")) / "scaler_params.json"
    data = preprocess_data(data, training_features, feature_transforms, constraint_var, scaler_output_path=scaler_path)

    if args.infer:
        logging.info("Skipping training. Running inference only...")
        run_inference(args, cfg, flavor, raw_sig_data, raw_bkg_data, training_features, feature_transforms, constraint_var, derived_vars_cfg)
        return

    logging.info("Creating data loaders...")
    train_loader, val_loader = make_dataloaders(
        data=data,
        training_features=training_features,
        constraint_var=constraint_var,
        batch_size=cfg.get("batch_size", 4096),
    )

    lightning_model = ABCDLightningModule(
        input_size=len(training_features),
        hidden_layers=cfg.get("architecture", [64, 32, 16]),
        learning_rate=cfg.get("learning_rate", 1e-3),
        bce_weight=cfg.get("bce_weight", 1.0),
        disco_lambda=cfg.get("disco_lambda", 0.0),
        flavor=flavor,
        use_batchnorm=cfg.get("use_batchnorm", True),
        dropout=cfg.get("dropout", 0.0),
        weight_decay=cfg.get("weight_decay", 1e-2),
        label_smoothing=cfg.get("label_smoothing", 0.0),
        use_lr_scheduler=cfg.get("use_lr_scheduler", True),
        lr_scheduler_patience=cfg.get("lr_scheduler_patience", 40),
        lr_scheduler_factor=cfg.get("lr_scheduler_factor", 0.5),
        lr_scheduler_min_lr=cfg.get("lr_scheduler_min_lr", 1e-6),
    )

    output_dir = cfg.get("output", "simple_abcdisco_output")

    run_training(
        model=lightning_model,
        train_loader=train_loader,
        val_loader=val_loader,
        output_dir=output_dir,
        max_epochs=cfg.get("n_epochs", 100),
        flavor=flavor,
        devices=cfg.get("devices", [0]),
        check_val_every_n_epoch=cfg.get("check_val_every_n_epoch", 1),
        early_stopping_patience=cfg.get("early_stopping_patience", 40),
        early_stopping_min_delta=cfg.get("early_stopping_min_delta", 1e-4),
    )

    logging.info("Training finished.")

    # Copy scaler params to the versioned checkpoint directory so they stay with the model
    checkpoint_path = _latest_checkpoint_from_config(cfg, flavor=flavor)
    scaler_path_final = _inference_output_dir_from_checkpoint(checkpoint_path) / "scaler_params.json"
    shutil.copy(str(scaler_path), str(scaler_path_final))
    logging.info("Copied scaler params to %s", scaler_path_final)

    logging.info("Running inference...")
    run_inference(args, cfg, flavor, raw_sig_data, raw_bkg_data, training_features, feature_transforms, constraint_var, derived_vars_cfg)
    logging.info("Inference finished.")

if __name__ == "__main__":
    main()
