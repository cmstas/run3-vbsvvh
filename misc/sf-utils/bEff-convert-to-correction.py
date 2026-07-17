#!/usr/bin/env python3
"""Convert merged weighted b-tag efficiency histograms into correctionlib JSON.

The legacy mode adds one exact MC sample.  Family mode discovers sample
directories, merges raw weighted yields into conservative physics families,
and writes one correction entry per family.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import uproot
import correctionlib.schemav2 as cs

from btag_eff_families import sample_family


FLAVORS = ("b", "c", "light")
WORKING_POINTS = ("T", "L", "LT", "N")
HISTOGRAM_NAMES = {"T": "T", "L": "L", "LT": "LT", "N": "N"}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--input", nargs="+",
                        help="All raw ROOT files for one sample (local paths or root:// URLs)")
    inputs.add_argument("--input-dir", type=Path,
                        help="Directory containing one raw-output directory per exact sample")
    parser.add_argument("--year", required=True, help="Metadata year, e.g. 2024Prompt")
    parser.add_argument("--channel", required=True, help="Preselection channel, e.g. 0lep_0FJ")
    parser.add_argument("--sample", help="Exact RDataFrame sample name from the input JSON (legacy mode)")
    parser.add_argument("--output", required=True, type=Path,
                        help="CorrectionSet JSON to create or update")
    parser.add_argument("--manifest", type=Path,
                        help="Family membership manifest (defaults beside --output in family mode)")
    return parser.parse_args()


def fold_pt_flow(values, path, hist_name):
    """Fold pT flow into edge bins, matching correctionlib's pT clamp flow."""
    # The producer rejects |eta| >= 2.5, so eta flow would indicate a
    # producer/application mismatch rather than information to clamp.
    if np.any(values[:, 0] != 0) or np.any(values[:, -1] != 0):
        raise ValueError(f"{path}:{hist_name} has unexpected eta under/overflow entries")
    central = values[1:-1, 1:-1].astype(float, copy=True)
    central[0, :] += values[0, 1:-1]
    central[-1, :] += values[-1, 1:-1]
    return central


def histogram_arrays_in_application_range(hist, path, hist_name):
    """Return weighted yields and their Sumw2 variances in application bins."""
    values, pt_edges, eta_edges = hist.to_numpy(flow=True)
    variances = hist.variances(flow=True)
    if variances is None:
        raise ValueError(f"{path}:{hist_name} has no Sumw2 information")
    return (fold_pt_flow(values, path, hist_name),
            fold_pt_flow(variances, path, hist_name),
            pt_edges[1:-1], eta_edges[1:-1])


def read_merged_histograms(paths, expected_year, expected_channel, expected_sample):
    """Return weighted yields, Sumw2 variances, and common pT/eta edges."""
    merged = {}
    merged_variances = {}
    edges = None
    for path in paths:
        with uproot.open(path) as root_file:
            if "btag_eff_format" not in root_file:
                raise ValueError(f"{path} has no b-tag efficiency format metadata")
            output_format = root_file["btag_eff_format"].member("fTitle")
            if not output_format.startswith("signed baseweight"):
                raise ValueError(
                    f"{path} contains obsolete unweighted b-tag histograms; regenerate it with the current --btag_eff workflow"
                )
            expected_metadata = {
                "btag_eff_year": expected_year,
                "btag_eff_channel": expected_channel,
                "btag_eff_sample": expected_sample,
            }
            for key, expected_value in expected_metadata.items():
                if key not in root_file:
                    raise ValueError(f"{path} has no {key} metadata")
                actual_value = root_file[key].member("fTitle")
                if actual_value != expected_value:
                    raise ValueError(f"{path} has {key}={actual_value!r}, expected {expected_value!r}")
            for flavor in FLAVORS:
                for state in ("den", *WORKING_POINTS):
                    hist_name = f"btag_{flavor}_{state}"
                    if hist_name not in root_file:
                        raise ValueError(f"{path} does not contain {hist_name}")
                    values, variances, pt_edges, eta_edges = histogram_arrays_in_application_range(
                        root_file[hist_name], path, hist_name)
                    if edges is None:
                        edges = (pt_edges, eta_edges)
                    elif not (np.array_equal(pt_edges, edges[0]) and np.array_equal(eta_edges, edges[1])):
                        raise ValueError(f"Histogram binning in {path} differs from the other inputs")
                    key = (flavor, state)
                    merged[key] = merged.get(key, np.zeros_like(values, dtype=float)) + values
                    merged_variances[key] = merged_variances.get(key, np.zeros_like(variances, dtype=float)) + variances
    return merged, merged_variances, edges


def add_counts(destination, source):
    for key, values in source.items():
        destination[key] = destination.get(key, np.zeros_like(values, dtype=float)) + values


def discover_family_histograms(input_dir, year, channel):
    """Read every exact sample directory and return raw yields grouped by family."""
    if not input_dir.is_dir():
        raise ValueError(f"--input-dir is not a directory: {input_dir}")
    grouped_counts, grouped_variances, grouped_edges, members = {}, {}, {}, {}
    for sample_dir in sorted(path for path in input_dir.iterdir() if path.is_dir()):
        roots = sorted(sample_dir.glob("*.root"))
        if not roots:
            # Allow diagnostics/manifests to live beside the sample directories.
            continue
        sample = sample_dir.name
        family = sample_family(sample)
        counts, variances, edges = read_merged_histograms(roots, year, channel, sample)
        if family in grouped_edges and not (
            np.array_equal(edges[0], grouped_edges[family][0]) and
            np.array_equal(edges[1], grouped_edges[family][1])
        ):
            raise ValueError(f"Histogram binning for {sample} differs within family {family}")
        grouped_edges[family] = edges
        grouped_counts.setdefault(family, {})
        grouped_variances.setdefault(family, {})
        add_counts(grouped_counts[family], counts)
        add_counts(grouped_variances[family], variances)
        members.setdefault(family, []).append(sample)
    if not members:
        raise ValueError(f"No sample directories found in {input_dir}")
    return grouped_counts, grouped_variances, grouped_edges, members


def invalid_count_bins(counts):
    """Return per-flavor masks for bins that cannot form physical efficiencies."""
    masks = {}
    for flavor in FLAVORS:
        denominator = counts[(flavor, "den")]
        tight = counts[(flavor, "T")]
        loose = counts[(flavor, "L")]
        loose_not_tight = counts[(flavor, "LT")]
        untagged = counts[(flavor, "N")]
        identity_scale = np.maximum.reduce([
            np.ones_like(denominator), np.abs(denominator), np.abs(tight),
            np.abs(loose), np.abs(loose_not_tight), np.abs(untagged),
        ])
        identity_tolerance = 1e-10 * identity_scale
        identity_invalid = ((np.abs(loose - (tight + loose_not_tight)) > identity_tolerance) |
                            (np.abs(denominator - (tight + loose_not_tight + untagged)) > identity_tolerance))

        # Signed generator weights can make a small MC bin statistically
        # pathological.  Do not silently turn it into an unweighted result:
        # reject it so the bin can be merged or supplied with more MC.
        tolerance = 1e-10 * identity_scale
        nonempty = np.abs(denominator) > tolerance
        invalid = identity_invalid | ((nonempty & ((denominator <= 0) | (tight < -tolerance) |
                                (loose < -tolerance) | (tight > loose + tolerance) |
                                (loose > denominator + tolerance))) |
                   (~nonempty & ((np.abs(tight) > tolerance) |
                                 (np.abs(loose_not_tight) > tolerance) |
                                 (np.abs(untagged) > tolerance))))
        masks[flavor] = invalid
    return masks


def validate_counts(counts):
    """Validate signed weighted yields before constructing an efficiency map."""
    for flavor, invalid in invalid_count_bins(counts).items():
        if np.any(invalid):
            bad_bins = np.argwhere(invalid).tolist()
            raise ValueError(
                f"Signed-weight b-tag yields are unphysical for flavor {flavor} in bins {bad_bins}; "
                "merge those bins or provide more MC before conversion"
            )


def repair_with_inclusive(counts, variances, inclusive_counts, inclusive_variances):
    """Replace only pathological family bins with the validated all-MC bin."""
    repaired = {key: values.copy() for key, values in counts.items()}
    repaired_variances = {key: values.copy() for key, values in variances.items()}
    replacements = {}
    for flavor, invalid in invalid_count_bins(counts).items():
        if not np.any(invalid):
            continue
        for state in ("den", *WORKING_POINTS):
            repaired[flavor, state][invalid] = inclusive_counts[flavor, state][invalid]
            repaired_variances[flavor, state][invalid] = inclusive_variances[flavor, state][invalid]
        replacements[flavor] = np.argwhere(invalid).tolist()
    validate_counts(repaired)
    return repaired, repaired_variances, replacements


def efficiency(values, denominator):
    output = np.zeros_like(values, dtype=float)
    np.divide(values, denominator, out=output, where=denominator > 0)
    return output


def poisson_efficiency_uncertainty(numerator_variance, denominator):
    """Poisson-style uncertainty sqrt(sumw2_numerator) / sumw_denominator.

    The producer stores signed nominal-weight yields and their Sumw2 values.
    For unit weights this is sqrt(N_tag) / N_all; with event weights, Sumw2
    provides the corresponding weighted Poisson variance.  Keep the raw ROOT
    histograms as well when exact future sample/channel merging is required.
    """
    output = np.zeros_like(numerator_variance, dtype=float)
    np.divide(np.sqrt(np.maximum(numerator_variance, 0.0)), denominator,
              out=output, where=denominator > 0)
    return output


def multibinning(values, pt_edges, eta_edges):
    return {
        "nodetype": "multibinning",
        "inputs": ["pt", "eta"],
        "edges": [pt_edges.tolist(), eta_edges.tolist()],
        "content": values.reshape(-1).tolist(),
        "flow": "clamp",
    }


def sample_category(sample, values, edges):
    pt_edges, eta_edges = edges
    flavor_entries = []
    for flavor in FLAVORS:
        wp_entries = []
        for wp in WORKING_POINTS:
            wp_entries.append({
                "key": wp,
                "value": multibinning(values[(flavor, wp)], pt_edges, eta_edges),
            })
        flavor_entries.append({
            "key": {"b": "B", "c": "C", "light": "L"}[flavor],
            "value": {"nodetype": "category", "input": "WP", "content": wp_entries},
        })
    return {
        "key": sample,
        "value": {"nodetype": "category", "input": "flavor", "content": flavor_entries},
    }


def make_correction(name, sample_entry, description, output_name, output_description):
    return {
        "name": name,
        "description": description,
        "version": 1,
        "inputs": [
            {"name": "sample", "type": "string", "description": "exact RDataFrame sample name"},
            {"name": "flavor", "type": "string", "description": "B/C/L"},
            {"name": "WP", "type": "string", "description": "T/L/LT/N"},
            {"name": "pt", "type": "real", "description": "selected AK4 jet pT"},
            {"name": "eta", "type": "real", "description": "selected AK4 jet eta"},
        ],
        "output": {"name": output_name, "type": "real", "description": output_description},
        "data": {"nodetype": "category", "input": "sample", "content": [sample_entry]},
    }


def update_output(path, correction_specs, replace_entries=False):
    if path.exists():
        payload = json.loads(path.read_text())
    else:
        payload = {"schema_version": 2, "description": "VBS VVH b-tag efficiencies", "corrections": [], "compound_corrections": []}

    expected_inputs = ["sample", "flavor", "WP", "pt", "eta"]
    replaced = set()
    for correction_name, sample_entry, description, output_name, output_description in correction_specs:
        correction = next((item for item in payload["corrections"] if item["name"] == correction_name), None)
        if correction is None:
            payload["corrections"].append(
                make_correction(correction_name, sample_entry, description, output_name, output_description))
            correction = payload["corrections"][-1]
        if ([item["name"] for item in correction["inputs"]] != expected_inputs or
                correction["output"]["name"] != output_name):
            raise ValueError(f"Existing {correction_name} has an incompatible schema; write a new output file")
        entries = correction["data"]["content"]
        if replace_entries and correction_name not in replaced:
            correction["data"]["content"] = []
            replaced.add(correction_name)
        if replace_entries:
            correction["data"]["content"].append(sample_entry)
        else:
            correction["data"]["content"] = [entry for entry in entries if entry["key"] != sample_entry["key"]]
            correction["data"]["content"].append(sample_entry)

    # Validate with correctionlib before touching the output file.
    cs.CorrectionSet.model_validate(payload)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n")


def main():
    args = parse_args()
    prefix = f"btag_{args.year}_{args.channel}"
    if args.input_dir:
        grouped_counts, grouped_variances, grouped_edges, members = discover_family_histograms(
            args.input_dir, args.year, args.channel)
        inclusive_counts, inclusive_variances = {}, {}
        for family in members:
            add_counts(inclusive_counts, grouped_counts[family])
            add_counts(inclusive_variances, grouped_variances[family])
        validate_counts(inclusive_counts)
        specs = []
        fallbacks = {}
        for family in sorted(members):
            counts, variances, fallback_bins = repair_with_inclusive(
                grouped_counts[family], grouped_variances[family],
                inclusive_counts, inclusive_variances)
            if fallback_bins:
                fallbacks[family] = fallback_bins
                print(f"{family}: replaced {sum(len(bins) for bins in fallback_bins.values())} pathological bins with all-MC yields")
            edges = grouped_edges[family]
            efficiency_values = {(flavor, wp): efficiency(counts[(flavor, wp)], counts[(flavor, "den")])
                                 for flavor in FLAVORS for wp in WORKING_POINTS}
            poisson_uncertainties = {(flavor, wp): poisson_efficiency_uncertainty(
                variances[(flavor, wp)], counts[(flavor, "den")])
                                   for flavor in FLAVORS for wp in WORKING_POINTS}
            specs.extend([
                (prefix, sample_category(family, efficiency_values, edges),
                 "Selected-AK4 UParTAK4 b-tag efficiency", "efficiency", "MC tagging efficiency"),
                (f"{prefix}_poisson_unc", sample_category(family, poisson_uncertainties, edges),
                 "Poisson uncertainty of selected-AK4 UParTAK4 b-tag efficiency", "uncertainty",
                 "Poisson-style MC tagging-efficiency uncertainty"),
            ])
        update_output(args.output, specs, replace_entries=True)
        manifest = args.manifest or args.output.with_name(
            f"btag_eff_{args.year}_{args.channel}_families.json")
        manifest.write_text(json.dumps({
            "year": args.year, "channel": args.channel,
            "input_dir": str(args.input_dir), "families": members,
            "inclusive_fallback_bins": fallbacks,
        }, indent=2) + "\n")
        print(f"Wrote {len(members)} family entries and manifest {manifest}")
        return

    if not args.sample:
        raise ValueError("--sample is required with --input")
    counts, variances, edges = read_merged_histograms(args.input, args.year, args.channel, args.sample)
    validate_counts(counts)
    efficiency_values = {(flavor, wp): efficiency(counts[(flavor, wp)], counts[(flavor, "den")])
                         for flavor in FLAVORS for wp in WORKING_POINTS}
    poisson_uncertainties = {(flavor, wp): poisson_efficiency_uncertainty(
        variances[(flavor, wp)], counts[(flavor, "den")])
                           for flavor in FLAVORS for wp in WORKING_POINTS}
    update_output(args.output, [
        (prefix, sample_category(args.sample, efficiency_values, edges),
         "Selected-AK4 UParTAK4 b-tag efficiency", "efficiency", "MC tagging efficiency"),
        (f"{prefix}_poisson_unc", sample_category(args.sample, poisson_uncertainties, edges),
         "Poisson uncertainty of selected-AK4 UParTAK4 b-tag efficiency", "uncertainty",
         "Poisson-style MC tagging-efficiency uncertainty"),
    ])


if __name__ == "__main__":
    main()
