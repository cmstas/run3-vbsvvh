"""Diagnostic plots for training inputs, the ABCD plane, and model performance."""

import logging
import re
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects
import matplotlib.pyplot as plt
from sklearn.metrics import auc, roc_curve
from tensorboard.backend.event_processing.event_accumulator import EventAccumulator

from common import data_length, score_column, to_flat_float_column


def save_tensorboard_plots(log_dir, output_dir):
    ea = EventAccumulator(str(log_dir))
    ea.Reload()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    for tag in ea.Tags().get("scalars", []):
        events = ea.Scalars(tag)
        steps = [e.step for e in events]
        values = [e.value for e in events]
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(steps, values, linewidth=2)
        ax.set_xlabel("Epoch")
        ax.set_ylabel(tag)
        ax.set_title(tag)
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        safe_tag = tag.replace("/", "_")
        plt.savefig(output_dir / f"{safe_tag}.png", dpi=150)
        plt.close()
        logging.info("Saved %s", output_dir / f"{safe_tag}.png")


def plot_input_feature_distributions(
    raw_data,
    training_features,
    train_idx,
    val_idx,
    output_dir,
    feature_transforms=None,
    skip_cols=None,
):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    labels = np.asarray(raw_data["label"])
    weights = np.asarray(raw_data["weight"]) if "weight" in raw_data else np.ones(data_length(raw_data), dtype=np.float64)
    weights = np.where(weights < 1e4, weights, 0.0)  # safeguard against extreme weights that can break the weighted density calculation

    train_mask = np.zeros(data_length(raw_data), dtype=bool)
    val_mask = np.zeros(data_length(raw_data), dtype=bool)
    train_mask[np.asarray(train_idx, dtype=np.int64)] = True
    val_mask[np.asarray(val_idx, dtype=np.int64)] = True

    training_feature_set = set(training_features)
    feature_transforms = feature_transforms or {}
    skip_cols = set(skip_cols or [])
    skip_cols |= {
        "label", "weight", "dataset_idx", "sample_idx", "split",
        "dnn_score", "dnn_0_score", "dnn_1_score",
    }

    def _safe_filename(name):
        return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(name))

    def _weighted_density(values, weights, bins):
        counts, edges = np.histogram(values, bins=bins, weights=weights, density=False)
        total = np.sum(counts)
        if total > 0:
            counts = counts / total
        return counts, edges

    def _numeric_full_column(arr):
        """Full-length float64 column, or None for non-numeric columns."""
        arr = np.asarray(arr)
        if arr.dtype != object and not np.issubdtype(arr.dtype, np.number):
            return None
        try:
            return to_flat_float_column(arr)
        except Exception:
            return None

    color_map = {
        ("bkg", "train"): "tab:blue",
        ("bkg", "val"): "tab:green",
        ("sig", "train"): "tab:red",
        ("sig", "val"): "tab:orange",
    }

    for feat in sorted(k for k in raw_data.keys() if k not in skip_cols):
        full_arr = _numeric_full_column(raw_data[feat])
        if full_arr is None:
            logging.info("Skipping non-numeric feature '%s' for input plotting", feat)
            continue

        valid_mask = np.isfinite(full_arr)
        if valid_mask.sum() < 2:
            logging.info("Skipping feature '%s' because it has fewer than 2 finite values", feat)
            continue

        plot_vals = full_arr[valid_mask]
        plot_labels = labels[valid_mask]
        plot_weights = weights[valid_mask]
        plot_train_mask = train_mask[valid_mask]
        plot_val_mask = val_mask[valid_mask]

        vmin = np.min(plot_vals)
        vmax = np.max(plot_vals)
        if vmin == vmax:
            eps = 0.5 if vmin == 0 else 0.05 * abs(vmin)
            vmin -= eps
            vmax += eps

        bins = np.linspace(vmin, vmax, 51)

        fig, ax = plt.subplots(figsize=(8, 6))

        subsets = [
            ("bkg", "train", (plot_labels == 0) & plot_train_mask, "Background train"),
            ("bkg", "val",   (plot_labels == 0) & plot_val_mask,   "Background val"),
            ("sig", "train", (plot_labels == 1) & plot_train_mask, "Signal train"),
            ("sig", "val",   (plot_labels == 1) & plot_val_mask,   "Signal val"),
        ]

        drew_any = False
        for cls_name, split_name, mask, legend_label in subsets:
            vals = plot_vals[mask]
            wts = plot_weights[mask]
            if len(vals) == 0 or np.sum(wts) <= 0:
                continue

            counts, edges = _weighted_density(vals, wts, bins=bins)
            ax.step(
                edges[:-1],
                counts,
                where="post",
                linewidth=2,
                color=color_map[(cls_name, split_name)],
                label=f"{legend_label} (n={len(vals)})",
            )
            drew_any = True

        if not drew_any:
            plt.close(fig)
            logging.info("Skipping feature '%s' because no drawable subsets were found", feat)
            continue

        ax.set_title(feat)
        ax.set_xlabel(feat)
        ax.set_ylabel("Unit-normalized weighted yield")
        legend = ax.legend(loc="best")
        ax.grid(True, alpha=0.3)

        if feat in training_feature_set:
            transform = feature_transforms.get(feat, "none")
            training_note = "★ used in training"
            if transform != "none":
                training_note += f" (transform: {transform})"

            fig.canvas.draw()
            if legend is not None:
                bbox_disp = legend.get_window_extent(fig.canvas.get_renderer())
                bbox_axes = bbox_disp.transformed(ax.transAxes.inverted())
                x_text = bbox_axes.x0
                y_text = max(0.02, bbox_axes.y0 - 0.06)
            else:
                x_text = 0.02
                y_text = 0.02

            ax.text(
                x_text,
                y_text,
                training_note,
                transform=ax.transAxes,
                fontsize=10,
                color="black",
                ha="left",
                va="top",
                bbox=dict(boxstyle="round,pad=0.25", facecolor="gold", alpha=0.25, edgecolor="goldenrod"),
            )

        output_path = output_dir / f"inputs_{_safe_filename(feat)}.png"
        plt.tight_layout()
        plt.savefig(output_path, dpi=150)
        plt.close()

        logging.info("Saved input feature plot to %s", output_path)


def plot_constraint_var_distribution(data, constraint_var, train_idx, val_idx, output_path):
    labels = np.asarray(data["label"])
    constraint = np.asarray(data[constraint_var])
    weights = np.asarray(data["weight"]) if "weight" in data else None

    fig, axes = plt.subplots(2, 2, figsize=(14, 10), gridspec_kw={"height_ratios": [3, 1]})

    for col, (mask, title) in enumerate([
        (labels == 0, "Background"),
        (labels == 1, "Signal"),
    ]):
        ax_main = axes[0, col]
        ax_ratio = axes[1, col]

        all_idx = np.where(mask)[0]
        train_mask_idx = np.intersect1d(all_idx, train_idx)
        val_mask_idx = np.intersect1d(all_idx, val_idx)

        all_vals = constraint[all_idx]
        train_vals = constraint[train_mask_idx]
        val_vals = constraint[val_mask_idx]

        all_w = weights[all_idx] if weights is not None else None
        train_w = weights[train_mask_idx] if weights is not None else None
        val_w = weights[val_mask_idx] if weights is not None else None

        bins = np.linspace(np.min(all_vals), np.max(all_vals), 51)

        all_counts, _ = np.histogram(all_vals, bins=bins, weights=all_w, density=False)
        train_counts, _ = np.histogram(train_vals, bins=bins, weights=train_w, density=False)
        val_counts, _ = np.histogram(val_vals, bins=bins, weights=val_w, density=False)

        all_counts_norm = all_counts / (all_counts.sum() + 1e-12)
        train_counts_norm = train_counts / (train_counts.sum() + 1e-12)
        val_counts_norm = val_counts / (val_counts.sum() + 1e-12)

        ax_main.step(bins[:-1], all_counts_norm, where="post", linewidth=2,
                     color="black", label=f"All ({len(all_vals)})")
        ax_main.step(bins[:-1], train_counts_norm, where="post", linewidth=2,
                     color="tab:blue", label=f"Train ({len(train_vals)})")
        ax_main.step(bins[:-1], val_counts_norm, where="post", linewidth=2,
                     color="tab:orange", label=f"Val ({len(val_vals)})")
        ax_main.set_ylabel("Density")
        ax_main.set_title(f"{title} {constraint_var} distribution")
        ax_main.legend()
        ax_main.grid(True, alpha=0.3)

        ratio = np.where(train_counts_norm > 1e-12, val_counts_norm / train_counts_norm, np.nan)
        ax_ratio.step(bins[:-1], ratio, where="post", linewidth=2, color="black")
        ax_ratio.axhline(1.0, color="red", linewidth=1, linestyle="--")
        ax_ratio.set_xlabel(constraint_var)
        ax_ratio.set_ylabel("Val / Train")
        ax_ratio.set_title("Val / Train ratio", fontsize=10)
        ax_ratio.set_ylim(0, 2)
        ax_ratio.grid(True, alpha=0.3)

    fig.suptitle(f"{constraint_var} train/val split distribution")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    logging.info("Saved constraint var distribution plot to %s", output_path)


def plot_weight_distributions(sig_data, bkg_data, output_path):
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for ax, data, title, color in [
        (axes[0], sig_data, "Signal weight distribution", "tab:red"),
        (axes[1], bkg_data, "Background weight distribution", "tab:blue"),
    ]:
        ax.hist(np.asarray(data["weight"]), bins=50, histtype="step", linewidth=2, color=color)
        ax.set_xlabel("Weight")
        ax.set_ylabel("Counts")
        ax.set_title(title)
        ax.set_yscale("log")
        ax.grid(True, alpha=0.3)

    fig.suptitle("Input weight distributions")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    logging.info("Saved weight distribution plot to %s", output_path)


# Profile-marker colours for the data plane. Both clear a 3:1 contrast ratio against the
# light end of the Greens ramp (violet 8.2:1, orange 3.1:1), where the old yellow scored
# 2.08:1, and they stay 95 dE apart under simulated protan/deutan/tritan vision.
PROFILE_COLOR = '#4a3aa7'
PROFILE_COLOR_WEIGHTED = '#eb6834'


def _profile_overlay(ax, x, y, bins, xrange, weights=None, color='yellow', label='Profile mean'):
    # Markers sit on a light->dark colormap, so no single colour contrasts with every
    # cell underneath. A white outline separates them from the dark end; the marker
    # colour itself carries the light end.
    ring = [matplotlib.patheffects.withStroke(linewidth=2.5, foreground='white')]

    bin_edges = np.linspace(xrange[0], xrange[1], bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    means, errors = [], []
    for lo, hi in zip(bin_edges[:-1], bin_edges[1:]):
        mask = (x >= lo) & (x < hi)
        vals = y[mask]
        w = weights[mask] if weights is not None else None
        if len(vals) > 0:
            mean = np.average(vals, weights=w)
            variance = np.average((vals - mean) ** 2, weights=w)
            n_eff = (np.sum(w) ** 2 / np.sum(w ** 2)) if w is not None else len(vals)
            errors.append(np.sqrt(variance / n_eff))
            means.append(mean)
        else:
            means.append(np.nan)
            errors.append(np.nan)
    means = np.array(means)
    errors = np.array(errors)
    container = ax.errorbar(bin_centers, means, yerr=errors, color=color,
                            fmt='o', markersize=4, linewidth=1.5, label=label)
    for artist in (container.lines[0], *container.lines[2]):
        if artist is not None:
            artist.set_path_effects(ring)
    ax.legend(fontsize=9)


def plot_abcd_plane(data, flavor, constraint_var, output_path, title_suffix=""):
    labels = np.asarray(data["label"])
    signal = {k: np.asarray(v)[labels == 1] for k, v in data.items()}
    background = {k: np.asarray(v)[labels == 0] for k, v in data.items()}

    score_col = score_column(flavor)

    bins = [50, 50]
    dnn_range = (0, 1)
    constrain_var_range = (0, max(np.asarray(data[constraint_var])))

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for ax, subset, title, cmap in [
        (axes[0], background, 'Background', 'Blues'),
        (axes[1], signal, 'Signal', 'Reds'),
    ]:
        h = ax.hist2d(subset[score_col], subset[constraint_var], bins=bins,
                      range=[dnn_range, constrain_var_range], cmap=cmap,
                      norm=matplotlib.colors.LogNorm())
        plt.colorbar(h[3], ax=ax, label='Counts')
        ax.set_title(title)
        ax.set_xlabel('DNN Score')
        ax.set_ylabel(constraint_var)
        _profile_overlay(ax, subset[score_col], subset[constraint_var], bins[0], dnn_range)
        _profile_overlay(ax, subset[score_col], subset[constraint_var], bins[0], dnn_range,
                         weights=subset["weight"] if "weight" in subset else None,
                         color='orange', label='Profile mean (weighted)')

    fig.suptitle(f'ABCD Plane{" - " + title_suffix if title_suffix else ""}')
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    logging.info("Saved ABCD plane plot to %s", output_path)


def plot_abcd_plane_data(data, flavor, constraint_var, output_path, blind_threshold=0.8, title_suffix=""):
    score = np.asarray(data[score_column(flavor)])
    constraint = np.asarray(data[constraint_var])
    weights = np.asarray(data["weight"]) if "weight" in data else None

    bins = [50, 50]
    dnn_range = (0, 1)
    constrain_var_range = (0, max(constraint))

    counts, xedges, yedges = np.histogram2d(
        score, constraint, bins=bins,
        range=[dnn_range, constrain_var_range], weights=weights,
    )

    # Blank any cell that reaches past the threshold on either axis, not just those fully
    # inside it, so a bin straddling the boundary cannot leak signal-region yield.
    blind_x = xedges[1:] > blind_threshold
    blind_y = yedges[1:] > blind_threshold
    counts[blind_x, :] = np.nan
    counts[:, blind_y] = np.nan

    fig, ax = plt.subplots(figsize=(8, 6))
    mesh = ax.pcolormesh(
        xedges, yedges,
        np.ma.masked_where(~np.isfinite(counts) | (counts <= 0), counts).T,
        cmap='Greens', norm=matplotlib.colors.LogNorm(),
    )
    plt.colorbar(mesh, ax=ax, label='Counts')

    # Profile exactly the cells left visible: the DNN columns left of the boundary, and
    # within them only events below the constraint boundary. Blinding on either axis means
    # even the low-DNN columns are truncated, so averaging every event there would fold in
    # yields this plot is meant to hide. The result is E[constraint | constraint < cut],
    # not the full mean -- a uniform cut across DNN bins, so a trend still reads as
    # correlation, but the values sit below the unblinded signal/background profiles.
    n_open = int(np.argmax(blind_x)) if blind_x.any() else bins[0]
    y_open = yedges[int(np.argmax(blind_y))] if blind_y.any() else yedges[-1]
    if n_open > 0:
        open_range = (xedges[0], xedges[n_open])
        # Cut on both axes rather than leaning on the binning range to drop the high-DNN
        # events, so nothing blinded is handed to the profile in the first place.
        visible = (score < xedges[n_open]) & (constraint < y_open)
        _profile_overlay(ax, score[visible], constraint[visible], n_open, open_range,
                         color=PROFILE_COLOR)
        _profile_overlay(ax, score[visible], constraint[visible], n_open, open_range,
                         weights=weights[visible] if weights is not None else None,
                         color=PROFILE_COLOR_WEIGHTED, label='Profile mean (weighted)')

    # Blinding on either axis leaves an L: the full high-DNN strip, plus the
    # high-constraint strip beside it.
    x0 = xedges[int(np.argmax(blind_x))] if blind_x.any() else xedges[-1]
    for rx, ry, rw, rh in (
        (x0, yedges[0], xedges[-1] - x0, yedges[-1] - yedges[0]),
        (xedges[0], y_open, x0 - xedges[0], yedges[-1] - y_open),
    ):
        if rw > 0 and rh > 0:
            ax.add_patch(plt.Rectangle(
                (rx, ry), rw, rh,
                facecolor='none', edgecolor='grey', hatch='//', linewidth=1.0, zorder=3,
            ))
    if blind_x.any() and blind_y.any():
        ax.text(0.5 * (x0 + xedges[-1]), 0.5 * (y_open + yedges[-1]), 'BLINDED',
                ha='center', va='center', color='grey', fontsize=11, fontweight='bold', zorder=4)

    ax.set_xlim(dnn_range)
    ax.set_ylim(constrain_var_range)
    ax.set_xlabel('DNN Score')
    ax.set_ylabel(constraint_var)
    ax.set_title(f'ABCD Plane - Data{" - " + title_suffix if title_suffix else ""}\n'
                 f'(blinded where either score > {blind_threshold}; profile over visible region only)',
                 fontsize=10)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    logging.info("Saved data ABCD plane plot to %s", output_path)


def plot_decorrelation_check(data, flavor, constraint_var, output_path, title_suffix=""):
    labels = np.asarray(data["label"])
    background = {k: np.asarray(v)[labels == 0] for k, v in data.items()}
    score_col = score_column(flavor)

    score_bins = [0.0, 0.25, 0.5, 0.75, 1.0]

    fig, ax = plt.subplots(figsize=(7, 6))
    for i in range(len(score_bins) - 1):
        mask = (background[score_col] >= score_bins[i]) & (background[score_col] < score_bins[i + 1])
        ax.hist(background[constraint_var][mask], bins=50, density=True,
                histtype='step', label=f'DNN score [{score_bins[i]}, {score_bins[i + 1]}]')

    ax.set_xlabel(constraint_var)
    ax.set_ylabel('Normalised counts')
    ax.legend()
    ax.set_title(f'{constraint_var} in DNN score slices (background){" - " + title_suffix if title_suffix else ""}')
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    logging.info("Saved decorrelation check plot to %s", output_path)


def plot_roc_curves(data, flavor, output_path):
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


def plot_score_densities(data, flavor, output_path):
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
        for mask, w, label, color in [
            (sig_mask, sig_w, "Signal", "tab:red"),
            (bkg_mask, bkg_w, "Background", "tab:blue"),
        ]:
            ax.hist(
                values[mask],
                bins=50,
                range=(0, 1),
                weights=w,
                density=True,
                histtype="step",
                linewidth=2,
                label=label,
                color=color,
            )
        ax.set_xlabel(col)
        ax.set_ylabel("Density")
        ax.set_title(f"Score density: {col}")
        ax.legend(loc="best")

    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


# Plot permutation importance as a horizontal bar chart, sorted by importance.
# Features that hurt the AUC most when shuffled appear at the top.
def plot_permutation_importance(baseline_auc, importances, output_path):
    sorted_feats = sorted(importances, key=importances.get)
    sorted_vals = [importances[f] for f in sorted_feats]

    fig, ax = plt.subplots(figsize=(10, max(6, len(sorted_feats) * 0.25)))
    colors = ["tab:red" if v > 0 else "tab:blue" for v in sorted_vals]
    ax.barh(sorted_feats, sorted_vals, color=colors)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Mean AUC drop when feature is shuffled")
    ax.set_title(f"Permutation Importance (baseline AUC={baseline_auc:.4f})")
    ax.grid(True, alpha=0.3, axis="x")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    logging.info("Saved %s", output_path)
