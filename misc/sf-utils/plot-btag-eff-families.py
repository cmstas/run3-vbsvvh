#!/usr/bin/env python3
"""Plot inter-family b-tag efficiency pulls from raw --btag_eff ROOT outputs."""

import argparse
import importlib.util
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import mplhep as hep
import numpy as np

from btag_eff_families import sample_family


LUMI_FB = {
    "2016preVFP": 19.5, "2016postVFP": 16.8, "2017": 41.5, "2018": 59.8,
    "2022Re-recoBCD": 8.1, "2022Re-recoE+PromptFG": 27.0,
    "2023PromptC": 17.8, "2023PromptD": 9.5, "2024Prompt": 109.95,
}


def converter_module():
    path = Path(__file__).with_name("bEff-convert-to-correction.py")
    spec = importlib.util.spec_from_file_location("btag_converter", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--year", required=True)
    parser.add_argument("--channel", required=True)
    parser.add_argument("--plot-dir", type=Path,
                        help="Persistent diagnostics directory (defaults beside btag_eff.json)")
    parser.add_argument("--lumi", type=float, help="Integrated luminosity in fb^-1")
    parser.add_argument("--sources", nargs="+", help="Only write detail PDFs for these source families")
    parser.add_argument("--skip-matrices", action="store_true")
    parser.add_argument("--skip-pulls", action="store_true")
    return parser.parse_args()


def collect_samples(conv, input_dir, year, channel):
    samples, family_counts, family_vars, family_edges, members = {}, {}, {}, {}, {}
    for sample_dir in sorted(path for path in input_dir.iterdir() if path.is_dir()):
        roots = sorted(sample_dir.glob("*.root"))
        if not roots:
            continue
        sample = sample_dir.name
        family = sample_family(sample)
        counts, variances, edges = conv.read_merged_histograms(roots, year, channel, sample)
        samples[sample] = (family, counts, variances, edges)
        family_counts.setdefault(family, {})
        family_vars.setdefault(family, {})
        conv.add_counts(family_counts[family], counts)
        conv.add_counts(family_vars[family], variances)
        if family in family_edges and not (np.array_equal(edges[0], family_edges[family][0]) and
                                           np.array_equal(edges[1], family_edges[family][1])):
            raise ValueError(f"Binning mismatch in family {family}")
        family_edges[family] = edges
        members.setdefault(family, []).append(sample)
    return samples, family_counts, family_vars, family_edges, members


def efficiencies(conv, counts, variances):
    values, uncertainties, valid = {}, {}, {}
    for flavor in conv.FLAVORS:
        den = counts[(flavor, "den")]
        for wp in conv.WORKING_POINTS:
            num = counts[(flavor, wp)]
            ok = (den > 0) & (num >= 0) & (num <= den)
            values[flavor, wp] = np.divide(num, den, out=np.full_like(den, np.nan), where=ok)
            uncertainties[flavor, wp] = np.divide(
                np.sqrt(np.maximum(variances[(flavor, wp)], 0.0)), den,
                out=np.full_like(den, np.nan), where=ok)
            valid[flavor, wp] = ok
    return values, uncertainties, valid


def pulls(reference, other, flavor, wp):
    values_a, errors_a, valid_a = reference
    values_b, errors_b, valid_b = other
    variance = errors_a[flavor, wp] ** 2 + errors_b[flavor, wp] ** 2
    valid = valid_a[flavor, wp] & valid_b[flavor, wp] & (variance > 0)
    result = np.full(values_a[flavor, wp].shape, np.nan)
    result[valid] = (values_a[flavor, wp][valid] - values_b[flavor, wp][valid]) / np.sqrt(variance[valid])
    return result


def cms_label(ax, lumi, energy):
    rounded_lumi = round(lumi, 2 - int(math.floor(math.log10(abs(lumi))))) if lumi else lumi
    hep.cms.label("Preliminary", data=True, lumi=rounded_lumi, lumi_format="{0:g}",
                  com=energy, ax=ax)


def matrix_plots(groups, results, args, lumi, energy, group_label="family"):
    names = sorted(groups)
    for flavor in ("b", "c", "light"):
        for wp in ("T", "L", "LT", "N"):
            matrix = np.full((len(names), len(names)), np.nan)
            for row, left in enumerate(names):
                for col, right in enumerate(names):
                    if row == col:
                        matrix[row, col] = 0.
                        continue
                    pull = pulls(results[left], results[right], flavor, wp)
                    finite = np.isfinite(pull)
                    if finite.any():
                        matrix[row, col] = np.mean(np.abs(pull[finite]) > 2.)
            fig, ax = plt.subplots(figsize=(max(10, 0.65 * len(names)), max(8, 0.55 * len(names))))
            image = ax.imshow(matrix, vmin=0, vmax=1, cmap="magma", interpolation="nearest")
            ax.set_xticks(range(len(names)), names, rotation=60, ha="right", fontsize=8)
            ax.set_yticks(range(len(names)), names, fontsize=8)
            ax.set_xlabel(f"Comparison {group_label}")
            ax.set_ylabel(f"Comparison {group_label}")
            fig.colorbar(image, ax=ax, label=r"Fraction of valid bins with $|$pull$|>2$")
            cms_label(ax, lumi, energy)
            fig.tight_layout()
            stem = args.plot_dir / f"compatibility_{flavor}_{wp}"
            fig.savefig(stem.with_suffix(".pdf"))
            fig.savefig(stem.with_suffix(".png"), dpi=160)
            plt.close(fig)


def money_plot(groups, results, args, lumi, energy, group_label="family"):
    """One compatibility matrix pooled over every flavor, WP, and kinematic bin."""
    names = sorted(groups)
    matrix = np.full((len(names), len(names)), np.nan)
    for row, left in enumerate(names):
        for col, right in enumerate(names):
            if row == col:
                matrix[row, col] = 0.
                continue
            values = []
            for flavor in ("b", "c", "light"):
                for wp in ("T", "L", "LT", "N"):
                    pull = pulls(results[left], results[right], flavor, wp)
                    values.append(np.abs(pull[np.isfinite(pull)]))
            pooled = np.concatenate(values) if values else np.array([])
            if pooled.size:
                matrix[row, col] = np.mean(pooled > 2.)
    fig, ax = plt.subplots(figsize=(max(10, 0.65 * len(names)), max(8, 0.55 * len(names))))
    image = ax.imshow(matrix, vmin=0, vmax=1, cmap="magma", interpolation="nearest")
    ax.set_xticks(range(len(names)), names, rotation=60, ha="right", fontsize=8)
    ax.set_yticks(range(len(names)), names, fontsize=8)
    ax.set_xlabel(f"Comparison {group_label}")
    ax.set_ylabel(f"Comparison {group_label}")
    for row in range(len(names)):
        for col in range(len(names)):
            value = matrix[row, col]
            if np.isfinite(value) and value < 0.2:
                ax.text(col, row, f"{value:.2f}", color="white", ha="center", va="center", fontsize=6)
    fig.colorbar(image, ax=ax, label=r"Pooled fraction of valid bins with $|$pull$|>2$")
    cms_label(ax, lumi, energy)
    fig.tight_layout()
    stem = args.plot_dir / "compatibility_all_flavours_all_wps"
    fig.savefig(stem.with_suffix(".pdf"))
    fig.savefig(stem.with_suffix(".png"), dpi=180)
    plt.close(fig)


def family_pull_pdfs(groups, results, args, lumi, energy, group_label="family"):
    names = sorted(groups)
    sources = args.sources or names
    unknown = sorted(set(sources) - set(names))
    if unknown:
        raise ValueError(f"Unknown requested source {group_label}s: {unknown}")
    for source in sources:
        others = [name for name in names if name != source]
        with PdfPages(args.plot_dir / f"pulls_{source}.pdf") as pdf:
            for flavor in ("b", "c", "light"):
                fig, axes = plt.subplots(1, 4, figsize=(22, max(7, 0.42 * len(others) + 3)), sharey=True)
                for ax, wp in zip(axes, ("T", "L", "LT", "N")):
                    rows = [np.clip(pulls(results[source], results[other], flavor, wp).reshape(-1), -5, 5)
                            for other in others]
                    image = ax.imshow(rows, aspect="auto", cmap="coolwarm", vmin=-5, vmax=5,
                                      interpolation="nearest")
                    ax.set_xlabel(rf"{wp}: flattened $(p_T, |\eta|)$ bin")
                    ax.set_xticks(np.arange(0, rows[0].size, 8))
                axes[0].set_yticks(range(len(others)), others, fontsize=8)
                axes[0].set_ylabel(f"Comparison {group_label}\n({source} minus other)")
                fig.colorbar(image, ax=axes, label="Pull (clipped at $\pm5$)")
                cms_label(axes[0], lumi, energy)
                fig.subplots_adjust(left=0.28, right=0.90, bottom=0.16, top=0.94, wspace=0.10)
                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)


def summary_json(groups, results, args):
    summary = {}
    for left in sorted(groups):
        for right in sorted(groups):
            if left == right:
                continue
            for flavor in ("b", "c", "light"):
                for wp in ("T", "L", "LT", "N"):
                    pull = pulls(results[left], results[right], flavor, wp)
                    finite = np.abs(pull[np.isfinite(pull)])
                    key = f"{left}__vs__{right}/{flavor}/{wp}"
                    summary[key] = {"valid_bins": int(finite.size),
                                    "fraction_gt2": float(np.mean(finite > 2.)) if finite.size else None,
                                    "fraction_gt3": float(np.mean(finite > 3.)) if finite.size else None,
                                    "max_abs_pull": float(np.max(finite)) if finite.size else None}
    (args.plot_dir / "compatibility_summary.json").write_text(json.dumps(summary, indent=2) + "\n")


def main():
    args = parse_args()
    conv = converter_module()
    if args.plot_dir is None:
        args.plot_dir = (Path(__file__).parents[2] / "preselection" / "corrections" /
                          "scalefactors" / "btagging" / "diagnostics" /
                          f"{args.year}_{args.channel}")
    args.plot_dir.mkdir(parents=True, exist_ok=True)
    _, counts, variances, _, families = collect_samples(conv, args.input_dir, args.year, args.channel)
    results = {}
    for family in families:
        # Sparse signed-weight bins are intentionally masked in comparisons;
        # the converter records any all-MC fallback used for application.
        results[family] = efficiencies(conv, counts[family], variances[family])
    lumi = args.lumi if args.lumi is not None else LUMI_FB.get(args.year)
    energy = 13.6 if args.year.startswith("202") else 13
    plt.style.use(hep.style.CMS)
    if not args.skip_matrices:
        matrix_plots(families, results, args, lumi, energy)
    # Individual pull PDFs are disabled for now; retain the helper for optional future use.
    # if not args.skip_pulls:
    #     family_pull_pdfs(families, results, args, lumi, energy)
    summary_json(families, results, args)


if __name__ == "__main__":
    main()
