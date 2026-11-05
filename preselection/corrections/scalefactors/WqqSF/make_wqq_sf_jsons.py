#!/usr/bin/env python3
"""Write the per-era Wqq SF correctionlib files RDF reads.

Input: the SF provider's single correctionlib file (correction
`particleNet_Wqq_SF`, inputs era / WP (string) / eta / pt / systematic), and
optionally a table of candle MC efficiencies of true W jets in the same layout
without the systematic level.

Output, one file per era, named with the RDF year string:
    Wqq_SF_<year>.json     year in 2016preVFP 2016postVFP 2017 2018 2024Prompt 2025
    Wqq_SF_Run2Comb.json   the combined Run 2 fit, shipped alongside

Each file holds
    Wqq_SF      (WP real, pt real, eta real, systematic)  cumulative SF, score > WP
    Wqq_eff_MC  (WP real, pt real)                        only if --eff is given

Both are piecewise linear in WP between the measured working points, so any WP
in [first, last] is valid and the measured points are reproduced exactly; a WP
outside that range is an error. pT is clamped to the measured bins, plus one
bin above the last edge that repeats the last bin's nominal with its up/down
distances multiplied by HIGH_PT_UNC_SCALE.

Usage:
    python3 make_wqq_sf_jsons.py [--input FILE] [--eff FILE] [--outdir DIR]
"""
import argparse
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_INPUT = os.path.join(HERE, "Wqq_SF_correctionlib_etainclusive.json")
INPUT_CORRECTION = "particleNet_Wqq_SF"

# RDF year string -> era key in the provider's file.
ERAS = {
    "2016preVFP": "2016preVFP",
    "2016postVFP": "2016postVFP",
    "2017": "2017",
    "2018": "2018",
    "2024Prompt": "2024",
    "2025": "2025",
    "Run2Comb": "r2comb",
}
SYSTEMATICS = ("nom", "up", "down")

HIGH_PT_EDGE = 13000.0
HIGH_PT_UNC_SCALE = 3.0
# The WP binning is [lo, hi): nudge the last edge so the last measured WP itself is valid.
LAST_WP_NUDGE = 1e-6


def read_provider_table(path, correction, with_systematic=True):
    """Return (pt_edges, {era: {wp_string: [per-pt-bin value or {syst: value}]}}).

    correction=None takes the file's first correction.
    """
    with open(path) as f:
        payload = json.load(f)
    if correction is None:
        corr = payload["corrections"][0]
    else:
        corr = next(c for c in payload["corrections"] if c["name"] == correction)
    table, pt_edges = {}, None
    for era in corr["data"]["content"]:
        table[era["key"]] = {}
        for wp in era["value"]["content"]:
            binning = wp["value"]
            if binning["nodetype"] != "binning" or binning["input"] != "pt":
                raise ValueError(f"{path}: expected a pt binning under era/WP")
            if pt_edges is None:
                pt_edges = [float(e) for e in binning["edges"]]
            elif [float(e) for e in binning["edges"]] != pt_edges:
                raise ValueError(f"{path}: pt edges differ between eras/WPs")
            cells = []
            for cell in binning["content"]:
                if with_systematic:
                    cells.append({c["key"]: float(c["value"]) for c in cell["content"]})
                else:
                    cells.append(float(cell))
            table[era["key"]][wp["key"]] = cells
    return pt_edges, table


def linear_segment(w0, v0, w1, v1, wp_input):
    slope = (v1 - v0) / (w1 - w0)
    return {
        "nodetype": "formula",
        "expression": "[0]+[1]*x",
        "parser": "TFormula",
        "variables": [wp_input],
        "parameters": [v0 - slope * w0, slope],
    }


def with_high_pt_bin(cells, with_systematic):
    """Append the bin above the last measured edge."""
    last = cells[-1]
    if not with_systematic:
        return cells + [last]
    nom = last["nom"]
    return cells + [{
        "nom": nom,
        "up": nom + HIGH_PT_UNC_SCALE * (last["up"] - nom),
        "down": max(0.0, nom - HIGH_PT_UNC_SCALE * (nom - last["down"])),
    }]


def interpolated_node(era_table, pt_edges, with_systematic):
    """WP binning (piecewise linear) -> pt binning -> [systematic] -> formula."""
    wps = sorted(era_table, key=float)
    wp_values = [float(w) for w in wps]
    wp_edges = wp_values[:-1] + [wp_values[-1] + LAST_WP_NUDGE]
    out_pt_edges = pt_edges + [HIGH_PT_EDGE]
    cells = {w: with_high_pt_bin(era_table[w], with_systematic) for w in wps}

    segments = []
    for i in range(len(wps) - 1):
        lo, hi = wps[i], wps[i + 1]
        pt_content = []
        for ipt in range(len(out_pt_edges) - 1):
            if with_systematic:
                pt_content.append({
                    "nodetype": "category",
                    "input": "systematic",
                    "content": [
                        {"key": s, "value": linear_segment(wp_values[i], cells[lo][ipt][s],
                                                           wp_values[i + 1], cells[hi][ipt][s], "WP")}
                        for s in SYSTEMATICS
                    ],
                })
            else:
                pt_content.append(linear_segment(wp_values[i], cells[lo][ipt],
                                                 wp_values[i + 1], cells[hi][ipt], "WP"))
        segments.append({
            "nodetype": "binning",
            "input": "pt",
            "edges": out_pt_edges,
            "content": pt_content,
            "flow": "clamp",
        })
    return {"nodetype": "binning", "input": "WP", "edges": wp_edges,
            "content": segments, "flow": "error"}, wps


def sf_correction(era_table, pt_edges, provenance):
    data, wps = interpolated_node(era_table, pt_edges, with_systematic=True)
    return {
        "name": "Wqq_SF",
        "version": 1,
        "description": (
            "GloParT-3 Wqq tagging SF for true W jets (FatJet_genMatch 1 or 2), cumulative: "
            "the SF of score > WP, score = (Xcs + Xqq/3) / (Xcs + Xqq/3 + QCD). "
            f"Piecewise linear in WP between the measured points {', '.join(wps)}; a WP outside "
            f"[{wps[0]}, {wps[-1]}] is an error. pt clamped; above {pt_edges[-1]:g} GeV the last "
            f"bin's nominal with {HIGH_PT_UNC_SCALE:g}x its up/down distances. eta is unused "
            f"(inclusive in |eta| < 2.5). systematic: nom/up/down, absolute SF values. "
            f"Source: {provenance}"
        ),
        "inputs": [
            {"name": "WP", "type": "real", "description": "Working point on the score (score > WP)"},
            {"name": "pt", "type": "real", "description": "Fat jet pt (nominal JEC/JER), GeV"},
            {"name": "eta", "type": "real", "description": "Fat jet eta (unused)"},
            {"name": "systematic", "type": "string", "description": "nom, up or down"},
        ],
        "output": {"name": "sf", "type": "real"},
        "data": data,
    }


def eff_correction(era_table, pt_edges, provenance):
    data, wps = interpolated_node(era_table, pt_edges, with_systematic=False)
    return {
        "name": "Wqq_eff_MC",
        "version": 1,
        "description": (
            "Candle MC efficiency of true W jets (tp2) to pass score > WP, for the fail SF "
            "(1 - SF*eff)/(1 - eff) and the multi-WP category weights. Piecewise linear in WP "
            f"between {', '.join(wps)}; pt clamped. Source: {provenance}"
        ),
        "inputs": [
            {"name": "WP", "type": "real", "description": "Working point on the score (score > WP)"},
            {"name": "pt", "type": "real", "description": "Fat jet pt (nominal JEC/JER), GeV"},
        ],
        "output": {"name": "eff", "type": "real"},
        "data": data,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--input", default=DEFAULT_INPUT, help="provider SF file (default: %(default)s)")
    parser.add_argument("--eff", default=None, help="candle MC efficiency table (same layout, no systematic)")
    parser.add_argument("--outdir", default=HERE, help="output directory (default: %(default)s)")
    args = parser.parse_args()

    pt_edges, sf_table = read_provider_table(args.input, INPUT_CORRECTION)
    eff_table = None
    if args.eff:
        eff_pt_edges, eff_table = read_provider_table(args.eff, None, with_systematic=False)
        if eff_pt_edges != pt_edges:
            raise ValueError("efficiency table pt edges differ from the SF table")

    for year, era_key in ERAS.items():
        if era_key not in sf_table:
            raise KeyError(f"era {era_key!r} not in {args.input}")
        provenance = f"{os.path.basename(args.input)}, era {era_key}"
        corrections = [sf_correction(sf_table[era_key], pt_edges, provenance)]
        if eff_table is not None:
            if sorted(eff_table[era_key], key=float) != sorted(sf_table[era_key], key=float):
                raise ValueError(f"era {era_key}: efficiency and SF working points differ")
            corrections.append(eff_correction(eff_table[era_key], pt_edges,
                                              f"{os.path.basename(args.eff)}, era {era_key}"))
        out = {
            "schema_version": 2,
            "description": f"GloParT-3 Wqq tagging scale factors, {year}",
            "corrections": corrections,
        }
        path = os.path.join(args.outdir, f"Wqq_SF_{year}.json")
        with open(path, "w") as f:
            json.dump(out, f, indent=1)
        print(f"wrote {path} ({', '.join(c['name'] for c in corrections)})")


if __name__ == "__main__":
    main()
