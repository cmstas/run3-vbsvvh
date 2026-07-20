"""ABCD DNN pipeline entry point.

Orchestrates the three run modes; the heavy lifting lives in the focused modules
(config, root_io, preprocessing, training, inference, plots, checkpoints):

  * train + infer (default): fit the BDT constraint (optional) and the DNN, then
    score sig/bkg MC with the best and last checkpoints.
  * --infer: skip training, score sig/bkg MC with the latest (or given) checkpoint.
  * --infer --data: score real data with the latest (or given) checkpoint.
"""

import logging
from argparse import ArgumentParser

import numpy as np

import bdt as bdt_lib
import checkpoints
import plots
from common import concat_sig_bkg, data_length
from config import RunConfig, load_yaml
from inference import prepare_inference_data, run_inference
from model import ABCDLightningModule
from preprocessing import apply_derived_vars, normalize_class_weights, preprocess_data
from root_io import apply_preselection, load_data
from training import collect_run_artifacts, make_dataloaders, run_training

def parse_args():
    parser = ArgumentParser()
    parser.add_argument("--config", required=True, help="Path to YAML config")
    parser.add_argument("--flavor", choices=["single", "double"], default=None, help="Training flavor: single (one output) or double (two outputs)")
    parser.add_argument("--data", action="store_true", help="Run inference data (without training) using the latest checkpoint from config")
    parser.add_argument("--infer", action="store_true", help="Skip training and run inference only")
    parser.add_argument("--checkpoint", default=None, help="Path to model checkpoint (.ckpt) for inference. If omitted, auto-picks newest checkpoint.")
    parser.add_argument("--output-path", default=None, help="Output predictions path (.parquet)")
    args = parser.parse_args()
    if args.data and not args.infer:
        parser.error("--data can only be used with --infer")
    return args, parser


def infer_on_data(args, parser, run_cfg):
    if "data_path" not in run_cfg.raw:
        parser.error("Config must contain 'data_path' when using --data")

    logging.info("Loading real data...")
    real_data = load_data(
        run_cfg.sample_paths("data"),
        run_cfg.load_features,
        run_cfg.extra_vars_for("data"),
        num_workers=run_cfg.io_workers,
    )
    logging.info("Data samples: %d", data_length(real_data))
    real_data = apply_preselection(real_data, run_cfg.preselection, label="Data")

    logging.info("Skipping training. Running inference on data only...")
    checkpoint_path = args.checkpoint or checkpoints.latest_checkpoint(run_cfg.output_dir, run_cfg.flavor)
    prepared = prepare_inference_data(run_cfg, real_data)
    run_inference(run_cfg, checkpoint_path, prepared, is_data=True, output_path=args.output_path)


def main():
    args, parser = parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    cfg = load_yaml(args.config)
    flavor = args.flavor if args.flavor is not None else cfg.get("flavor", "single")
    logging.info("Using flavor=%s", flavor)

    run_cfg = RunConfig(cfg, flavor)

    if args.data:
        infer_on_data(args, parser, run_cfg)
        return

    logging.info("Loading signal data...")
    sig_data = load_data(
        run_cfg.sample_paths("sig", for_inference=args.infer), run_cfg.load_features, run_cfg.extra_vars_for("sig"),
        num_workers=run_cfg.io_workers,
    )
    logging.info("Loading background data...")
    bkg_data = load_data(
        run_cfg.sample_paths("bkg", for_inference=args.infer), run_cfg.load_features, run_cfg.extra_vars_for("bkg"),
        num_workers=run_cfg.io_workers,
    )

    logging.info("Signal samples: %d, Background samples: %d", data_length(sig_data), data_length(bkg_data))
    sig_data = apply_preselection(sig_data, run_cfg.preselection, label="Signal")
    bkg_data = apply_preselection(bkg_data, run_cfg.preselection, label="Background")

    # Train (or load) the BDT and attach its score so it can serve as the DNN
    # decorrelation constraint. Uses the raw physics weights (class balancing is
    # handled inside the BDT trainer), so this must run before normalize_class_weights.
    if run_cfg.use_bdt:
        if args.infer:
            logging.info("Loading existing BDT (inference mode)...")
            bdt_model = bdt_lib.load_bdt(run_cfg.output_dir)
        else:
            bdt_model = bdt_lib.train_bdt(sig_data, bkg_data, run_cfg.bdt_features, cfg, run_cfg.output_dir)
        sig_data[bdt_lib.BDT_SCORE_NAME] = bdt_lib.predict_bdt(bdt_model, sig_data, run_cfg.bdt_features)
        bkg_data[bdt_lib.BDT_SCORE_NAME] = bdt_lib.predict_bdt(bdt_model, bkg_data, run_cfg.bdt_features)

    # Keep copies of the original data (with raw weights and features) for inference
    raw_sig_data = {k: np.copy(v) for k, v in sig_data.items()}
    raw_bkg_data = {k: np.copy(v) for k, v in bkg_data.items()}

    if args.infer:
        logging.info("Skipping training. Running inference only...")
        checkpoint_path = args.checkpoint or checkpoints.latest_checkpoint(run_cfg.output_dir, run_cfg.flavor)
        prepared = prepare_inference_data(run_cfg, concat_sig_bkg(raw_sig_data, raw_bkg_data))
        run_inference(run_cfg, checkpoint_path, prepared, output_path=args.output_path)
        return

    run_cfg.output_dir.mkdir(parents=True, exist_ok=True)
    plots.plot_weight_distributions(sig_data, bkg_data, run_cfg.output_dir / "input_weight_distributions.png")

    sig_data, bkg_data = normalize_class_weights(sig_data, bkg_data)

    logging.info("Preprocessing data...")
    data = concat_sig_bkg(sig_data, bkg_data)
    data = apply_derived_vars(data, run_cfg.derived_vars)
    # Save scaler params alongside the output so inference can apply identical scaling
    data = preprocess_data(
        data,
        run_cfg.training_features,
        run_cfg.feature_transforms,
        run_cfg.constraint_var,
        scaler_output_path=run_cfg.output_dir / "scaler_params.json",
    )

    logging.info("Creating data loaders...")
    train_loader, val_loader, train_idx, val_idx = make_dataloaders(
        data=data,
        training_features=run_cfg.training_features,
        constraint_var=run_cfg.constraint_var,
        batch_size=run_cfg.batch_size,
    )

    plots.plot_constraint_var_distribution(
        data, run_cfg.constraint_var, train_idx, val_idx,
        run_cfg.output_dir / "constraint_var_distribution.png",
    )

    raw_plot_data = apply_derived_vars(concat_sig_bkg(raw_sig_data, raw_bkg_data), run_cfg.derived_vars)
    plots.plot_input_feature_distributions(
        raw_data=raw_plot_data,
        training_features=run_cfg.training_features,
        train_idx=train_idx,
        val_idx=val_idx,
        output_dir=run_cfg.output_dir,
        feature_transforms=run_cfg.feature_transforms,
    )

    lightning_model = ABCDLightningModule(**run_cfg.model_kwargs(input_size=len(run_cfg.training_features)))

    trainer = run_training(
        model=lightning_model,
        train_loader=train_loader,
        val_loader=val_loader,
        output_dir=run_cfg.output_dir,
        max_epochs=cfg.get("n_epochs", 100),
        flavor=flavor,
        devices=cfg.get("devices", [0]),
        check_val_every_n_epoch=cfg.get("check_val_every_n_epoch", 1),
        early_stopping_patience=cfg.get("early_stopping_patience", 40),
        early_stopping_min_delta=cfg.get("early_stopping_min_delta", 1e-4),
    )

    logging.info("Training finished.")
    version_dir = checkpoints.run_dir_for_checkpoint(
        checkpoints.latest_checkpoint(run_cfg.output_dir, flavor)
    )
    collect_run_artifacts(version_dir, args.config, run_cfg.output_dir)

    logging.info("Running inference...")
    prepared = prepare_inference_data(run_cfg, concat_sig_bkg(raw_sig_data, raw_bkg_data))
    for label, ckpt_path in checkpoints.best_and_last(trainer, run_cfg.output_dir, flavor):
        logging.info("Running inference for %s checkpoint: %s", label, ckpt_path)
        run_inference(run_cfg, ckpt_path, prepared, train_idx=train_idx, val_idx=val_idx,
                      output_path=args.output_path)
    logging.info("Inference finished.")


if __name__ == "__main__":
    main()
