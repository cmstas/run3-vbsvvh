#!/usr/bin/env python3
"""Compare b-tag efficiencies across all channels or all sample families."""

import argparse
import importlib.util
from pathlib import Path

import matplotlib.pyplot as plt
import mplhep as hep

from btag_eff_families import DEFAULT_CONFIG, final_channel, final_group


# These selections differ only in trigger paths and can share events.  Their
# covariance is not modelled by the histogram-level diagnostic.
OVERLAPPING_CHANNEL_PAIRS = {
    frozenset(("0lep_1FJ", "0lep_1FJ_met")),
    frozenset(("0lep_2FJ", "0lep_2FJ_met")),
}


def load_module(filename, name):
    path = Path(__file__).with_name(filename)
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--channel-input", action="append", required=True, metavar="CHANNEL=DIR",
                        help="One completed b-tag output directory per channel; repeat this option")
    parser.add_argument("--channel-manifest", action="append", metavar="CHANNEL=PATH",
                        help="Optional Slurm manifest for the matching --channel-input; repeat this option")
    parser.add_argument("--year", required=True)
    parser.add_argument("--mode", choices=("families", "channels"), required=True,
                        help="families: sum each family over channels; channels: sum all samples per channel")
    parser.add_argument("--plot-dir", type=Path)
    parser.add_argument("--lumi", type=float)
    parser.add_argument("--sources", nargs="+", help="Only write detail PDFs for these groups")
    parser.add_argument("--skip-matrices", action="store_true")
    parser.add_argument("--skip-pulls", action="store_true")
    parser.add_argument("--family-config", type=Path, default=DEFAULT_CONFIG,
                        help="Shared preliminary/final grouping YAML")
    parser.add_argument("--final", action="store_true",
                        help="Apply final_merges from the shared YAML before plotting")
    return parser.parse_args()


def channel_inputs(items):
    parsed = {}
    for item in items:
        if "=" not in item:
            raise ValueError(f"Invalid --channel-input {item!r}; use CHANNEL=DIR")
        channel, directory = item.split("=", 1)
        if channel in parsed:
            raise ValueError(f"Duplicate channel input: {channel}")
        parsed[channel] = Path(directory)
    return parsed


def add_all(conv, target, source):
    target.setdefault("counts", {})
    target.setdefault("variances", {})
    conv.add_counts(target["counts"], source[0])
    conv.add_counts(target["variances"], source[1])


def reject_overlapping_channels(channels):
    selected = set(channels)
    overlaps = [sorted(pair) for pair in OVERLAPPING_CHANNEL_PAIRS if pair <= selected]
    if overlaps:
        raise ValueError("Cannot sum overlapping analysis channels without covariance modelling: " +
                         ", ".join(" + ".join(pair) for pair in overlaps))


def main():
    args = parse_args()
    conv = load_module("bEff-convert-to-correction.py", "btag_converter")
    plots = load_module("plot-btag-eff-families.py", "btag_plots")
    inputs = channel_inputs(args.channel_input)
    manifests = channel_inputs(args.channel_manifest or [])
    unknown_manifests = set(manifests) - set(inputs)
    if unknown_manifests:
        raise ValueError(f"Manifest supplied for unknown channel(s): {sorted(unknown_manifests)}")
    if args.mode == "families" and not args.final:
        reject_overlapping_channels(inputs)
    grouped = {}
    completeness = []
    for channel, directory in inputs.items():
        _, family_counts, family_variances, _, families, complete = plots.collect_samples(
            conv, directory, args.year, channel, manifests.get(channel), args.family_config)
        completeness.append(complete is not None)
        sample_groups = {family: final_group("samples", family, args.family_config)
                         for family in families} if args.final else {family: family for family in families}
        channel_group = final_channel(channel, args.family_config) if args.final else channel
        if args.mode == "families":
            for family in families:
                add_all(conv, grouped.setdefault(sample_groups[family], {}),
                        (family_counts[family], family_variances[family]))
        else:
            for family in families:
                add_all(conv, grouped.setdefault(channel_group, {}),
                        (family_counts[family], family_variances[family]))
    results = {name: plots.efficiencies(conv, data["counts"], data["variances"])
               for name, data in grouped.items()}
    if args.plot_dir is None:
        suffix = ("final_all_channels_families" if args.mode == "families" else "final_all_samples_channels") if args.final else ("all_channels_families" if args.mode == "families" else "all_samples_channels")
        args.plot_dir = (Path(__file__).parents[2] / "preselection" / "corrections" /
                         "scalefactors" / "btagging" / "diagnostics" / f"{args.year}_{suffix}")
    args.plot_dir.mkdir(parents=True, exist_ok=True)
    args.channel = "all channels" if args.mode == "families" else "all samples"
    args.input_completeness_verified = bool(completeness) and all(completeness)
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
