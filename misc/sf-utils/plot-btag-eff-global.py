#!/usr/bin/env python3
"""Compare b-tag efficiencies across all channels or all sample families."""

import argparse
import importlib.util
from pathlib import Path

import matplotlib.pyplot as plt
import mplhep as hep

from btag_eff_families import load_config


def load_module(filename, name):
    path = Path(__file__).with_name(filename)
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, required=True,
                        help="Root containing retained YAML source channels and their manifests")
    parser.add_argument("--year", required=True)
    parser.add_argument("--mode", choices=("families", "channels"), required=True,
                        help="families: sum each family over channels; channels: sum all samples per channel")
    parser.add_argument("--plot-dir", type=Path)
    parser.add_argument("--lumi", type=float)
    parser.add_argument("--sources", nargs="+", help="Only write detail PDFs for these groups")
    parser.add_argument("--skip-matrices", action="store_true")
    parser.add_argument("--skip-pulls", action="store_true")
    parser.add_argument("--final", action="store_true",
                        help="Required: compare the final YAML channel/sample groups")
    return parser.parse_args()


def add_all(conv, target, source):
    target.setdefault("counts", {})
    target.setdefault("variances", {})
    conv.add_counts(target["counts"], source[0])
    conv.add_counts(target["variances"], source[1])


def main():
    args = parse_args()
    if not args.final:
        raise ValueError("Global diagnostics are final-only; pass --final")
    conv = load_module("bEff-convert-to-correction.py", "btag_converter")
    plots = load_module("plot-btag-eff-families.py", "btag_plots")
    config = load_config()
    counts, variances, _, _, _, ignored = conv.discover_final_histograms(
        args.input_root, args.year, config)
    grouped = {}
    for (channel, sample), channel_counts in counts.items():
        target = sample if args.mode == "families" else channel
        add_all(conv, grouped.setdefault(target, {}), (channel_counts, variances[channel, sample]))
    results = {name: plots.efficiencies(conv, data["counts"], data["variances"])
               for name, data in grouped.items()}
    if args.plot_dir is None:
        suffix = "final_all_channels_families" if args.mode == "families" else "final_all_samples_channels"
        args.plot_dir = (Path(__file__).parents[2] / "preselection" / "corrections" /
                         "scalefactors" / "btagging" / "diagnostics" / f"{args.year}_{suffix}")
    args.plot_dir.mkdir(parents=True, exist_ok=True)
    args.channel = "all channels" if args.mode == "families" else "all samples"
    args.input_completeness_verified = True
    args.ignored_source_channels = ignored
    lumi = args.lumi if args.lumi is not None else plots.LUMI_FB.get(args.year)
    energy = 13.6 if args.year.startswith("202") else 13
    plt.style.use(hep.style.CMS)
    label = "family" if args.mode == "families" else "channel"
    if not args.skip_matrices:
        plots.matrix_plots(grouped, results, args, lumi, energy, label)
    plots.money_plot(grouped, results, args, lumi, energy, label)
    # Individual pull PDFs are disabled for now; retain the helper for optional future use.
    # if not args.skip_pulls:
    #     plots.family_pull_pdfs(grouped, results, args, lumi, energy, label)
    plots.summary_json(grouped, results, args)


if __name__ == "__main__":
    main()
