"""Reading event columns from ROOT ntuples into flat numpy dicts."""

import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import uproot
from tqdm import tqdm

from common import (
    apply_mask,
    concat_chunks,
    data_length,
    evaluate_expression,
    flatten_cell,
)
from systematics import SYST_WEIGHTS, ratio_columns


def _syst_ratios(values):
    """Split a length-3 ``[nominal, up, down]`` weight branch into up/down ratios.

    The nominal element is already folded into the event weight upstream, so the
    varied weight is ``weight * element / nominal``. A zero or non-finite nominal
    would make that meaningless, so those events fall back to a ratio of 1 (no
    variation) rather than propagating an inf into the yields.
    """
    triplets = np.asarray([np.asarray(v, dtype=np.float64) for v in values])
    nominal, up, down = triplets[:, 0], triplets[:, 1], triplets[:, 2]

    safe = np.isfinite(nominal) & (nominal != 0)
    ratio_up = np.where(safe, up / np.where(safe, nominal, 1.0), 1.0)
    ratio_dn = np.where(safe, down / np.where(safe, nominal, 1.0), 1.0)
    return ratio_up, ratio_dn


def _read_root_frame(path, branches, dataset_idx, sample_idx):
    with uproot.open(path) as root_file:
        tree = root_file["Events"]
        # Data has no MC weight branches at all, so drop the systematics that are
        # absent rather than letting uproot raise. Only the systematics are
        # forgiven: any other missing branch is a config error and must still
        # fail loudly here rather than silently producing a column-less frame.
        available = tree.keys()
        requested = [b for b in branches if b in available or b not in SYST_WEIGHTS]
        arrays = tree.arrays(requested, library="np")

    if not isinstance(arrays, dict):
        arrays = {name: arrays[name] for name in arrays.dtype.names}

    columns = {}
    n_events = None
    for branch in branches:
        if branch not in arrays:
            continue

        values = np.asarray(arrays[branch])

        # Systematic weights are length-3 vectors; flatten_cell would silently
        # keep only the nominal, so expand them into ratio columns instead.
        if branch in SYST_WEIGHTS:
            up_col, dn_col = ratio_columns(branch)
            ratio_up, ratio_dn = _syst_ratios(values)
            n_events = len(ratio_up) if n_events is None else min(n_events, len(ratio_up))
            columns[up_col] = ratio_up
            columns[dn_col] = ratio_dn
            continue

        if values.ndim == 1 and values.dtype != object:
            clean_values = values
        else:
            flat_values = [flatten_cell(v) for v in values]
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


def load_data(paths, features, extra_vars, num_workers=1):
    branches = list(dict.fromkeys(features + extra_vars))

    # Sample name from the directory layout: <sample>/<chunk>/file.root or
    # <sample>/file.root. Used only as a stable stratification key.
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

    data = concat_chunks(chunks)
    logging.info("Loaded %d events from %d files.", data_length(data), len(paths))
    if "weight" in data:
        weights = np.asarray(data["weight"])
        data = apply_mask(data, (weights > 0) & (weights < 1e4))

    return data


def apply_preselection(data, expression, label=""):
    """Apply the config preselection expression (no-op if expression is falsy)."""
    if not expression:
        return data
    mask = np.asarray(evaluate_expression(data, expression), dtype=bool)
    out = apply_mask(data, mask)
    logging.info("After preselection - %s samples: %d", label or "input", data_length(out))
    return out
