#!/usr/bin/env python3
"""Convert merged raw b-tag efficiency histograms into correctionlib JSON.

Each invocation adds (or replaces) one exact MC sample in one
year/channel correction.  Pass every worker ROOT file for that sample; raw
counts are summed before efficiencies are calculated.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import uproot
import correctionlib.schemav2 as cs


FLAVORS = ("b", "c", "light")
WORKING_POINTS = ("T", "L", "LT", "N")
HISTOGRAM_NAMES = {"T": "T", "L": "L", "LT": "LT", "N": "N"}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", nargs="+", required=True,
                        help="All raw ROOT files for one sample (local paths or root:// URLs)")
    parser.add_argument("--year", required=True, help="Metadata year, e.g. 2024Prompt")
    parser.add_argument("--channel", required=True, help="Preselection channel, e.g. 0lep_0FJ")
    parser.add_argument("--sample", required=True,
                        help="Exact RDataFrame sample name from the input JSON")
    parser.add_argument("--output", required=True, type=Path,
                        help="CorrectionSet JSON to create or update")
    return parser.parse_args()


def read_merged_histograms(paths):
    """Return raw summed arrays and common pT/eta edges for all flavors."""
    merged = {}
    edges = None
    for path in paths:
        with uproot.open(path) as root_file:
            for flavor in FLAVORS:
                for state in ("den", *WORKING_POINTS):
                    hist_name = f"btag_{flavor}_{state}"
                    if hist_name not in root_file:
                        raise ValueError(f"{path} does not contain {hist_name}")
                    values, pt_edges, eta_edges = root_file[hist_name].to_numpy(flow=False)
                    if edges is None:
                        edges = (pt_edges, eta_edges)
                    elif not (np.array_equal(pt_edges, edges[0]) and np.array_equal(eta_edges, edges[1])):
                        raise ValueError(f"Histogram binning in {path} differs from the other inputs")
                    key = (flavor, state)
                    merged[key] = merged.get(key, np.zeros_like(values, dtype=float)) + values
    return merged, edges


def validate_counts(counts):
    """Validate nested-WP identities before constructing an efficiency map."""
    for flavor in FLAVORS:
        denominator = counts[(flavor, "den")]
        tight = counts[(flavor, "T")]
        loose = counts[(flavor, "L")]
        loose_not_tight = counts[(flavor, "LT")]
        untagged = counts[(flavor, "N")]
        if np.any(denominator < 0) or any(np.any(counts[(flavor, state)] < 0)
                                      for state in WORKING_POINTS):
            raise ValueError(f"Negative raw b-tag count found for flavor {flavor}")
        if np.any(tight > loose) or np.any(loose > denominator):
            raise ValueError(f"Nested-WP counts are inconsistent for flavor {flavor}")
        if not np.allclose(loose, tight + loose_not_tight, rtol=0, atol=1e-9):
            raise ValueError(f"L != T + LT for flavor {flavor}")
        if not np.allclose(denominator, tight + loose_not_tight + untagged, rtol=0, atol=1e-9):
            raise ValueError(f"denominator != T + LT + N for flavor {flavor}")


def efficiency(values, denominator):
    output = np.zeros_like(values, dtype=float)
    np.divide(values, denominator, out=output, where=denominator > 0)
    return output


def multibinning(values, pt_edges, eta_edges):
    return {
        "nodetype": "multibinning",
        "inputs": ["pt", "eta"],
        "edges": [pt_edges.tolist(), eta_edges.tolist()],
        "content": values.reshape(-1).tolist(),
        "flow": "clamp",
    }


def sample_category(sample, counts, edges):
    pt_edges, eta_edges = edges
    flavor_entries = []
    for flavor in FLAVORS:
        denominator = counts[(flavor, "den")]
        wp_entries = []
        for wp in WORKING_POINTS:
            wp_entries.append({
                "key": wp,
                "value": multibinning(efficiency(counts[(flavor, wp)], denominator), pt_edges, eta_edges),
            })
        flavor_entries.append({
            "key": {"b": "B", "c": "C", "light": "L"}[flavor],
            "value": {"nodetype": "category", "input": "WP", "content": wp_entries},
        })
    return {
        "key": sample,
        "value": {"nodetype": "category", "input": "flavor", "content": flavor_entries},
    }


def make_correction(name, sample_entry):
    return {
        "name": name,
        "description": "Selected-AK4 UParTAK4 b-tag efficiency",
        "version": 1,
        "inputs": [
            {"name": "sample", "type": "string", "description": "exact RDataFrame sample name"},
            {"name": "flavor", "type": "string", "description": "B/C/L"},
            {"name": "WP", "type": "string", "description": "T/L/LT/N"},
            {"name": "pt", "type": "real", "description": "selected AK4 jet pT"},
            {"name": "eta", "type": "real", "description": "selected AK4 jet eta"},
        ],
        "output": {"name": "efficiency", "type": "real", "description": "MC tagging efficiency"},
        "data": {"nodetype": "category", "input": "sample", "content": [sample_entry]},
    }


def update_output(path, correction_name, sample_entry):
    if path.exists():
        payload = json.loads(path.read_text())
    else:
        payload = {"schema_version": 2, "description": "VBS VVH b-tag efficiencies", "corrections": [], "compound_corrections": []}

    correction = next((item for item in payload["corrections"] if item["name"] == correction_name), None)
    if correction is None:
        payload["corrections"].append(make_correction(correction_name, sample_entry))
    else:
        expected_inputs = ["sample", "flavor", "WP", "pt", "eta"]
        if [item["name"] for item in correction["inputs"]] != expected_inputs:
            raise ValueError(f"Existing {correction_name} has an incompatible schema; write a new output file")
        entries = correction["data"]["content"]
        correction["data"]["content"] = [entry for entry in entries if entry["key"] != sample_entry["key"]]
        correction["data"]["content"].append(sample_entry)

    # Validate with correctionlib before touching the output file.
    cs.CorrectionSet.model_validate(payload)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n")


def main():
    args = parse_args()
    counts, edges = read_merged_histograms(args.input)
    validate_counts(counts)
    update_output(args.output, f"btag_{args.year}_{args.channel}",
                  sample_category(args.sample, counts, edges))


if __name__ == "__main__":
    main()
