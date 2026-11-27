#!/usr/bin/env python3
"""Write the per-era Xbb SF correctionlib files RDF reads.

Input, a folder with the files as the provider delivers them (as in PR #138, commit 837327abd):
    Xbb_SF_<era>.json              correction `Xbb_SF` (WP real, pt real, systematic), one file per era
    Xbb_effMC_correctionlib.json   correction `Xbb_effMC` (era string, WP string, pt real), all eras

Output, one file per era, named with the RDF year string:
    Xbb_SF_<year>.json   year in 2016preVFP 2016postVFP 2017 2018 2022Re-recoBCD 2022Re-recoE+PromptFG
                                 2023PromptC 2023PromptD 2024Prompt 2025

Each file holds
    Xbb_SF      (WP real, pt real, systematic)   cumulative SF, score > WP
    Xbb_eff_MC  (WP real, pt real)               candle MC efficiency of true b jets

Both are piecewise linear in WP between the measured working points, so any WP
in [first, last] is valid and the measured points are reproduced exactly; a WP
outside that range is an error. pT is clamped to the measured bins. The
provider's last pT bin is open-ended; as for Wqq, it is kept up to
HIGH_PT_EDGE, and above it one more bin repeats its nominal with its up/down
distances multiplied by HIGH_PT_UNC_SCALE (the efficiency is repeated as is).

Usage:
    python3 make_xbb_sf_jsons.py --provider DIR [--outdir DIR]
"""
import argparse
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
EFF_FILE = "Xbb_effMC_correctionlib.json"

# RDF year string -> era name in the provider's files (file name Xbb_SF_<era>.json and
# era key of the efficiency). The physical splits are the same: preEE = 2022 C+D,
# postEE = 2022 E+F+G, preBPix = 2023 C, postBPix = 2023 D.
ERAS = {
    "2016preVFP": "2016_PreVFP",
    "2016postVFP": "2016_PostVFP",
    "2017": "2017",
    "2018": "2018",
    "2022Re-recoBCD": "2022_preEE",
    "2022Re-recoE+PromptFG": "2022_postEE",
    "2023PromptC": "2023_preBPix",
    "2023PromptD": "2023_postBPix",
    "2024Prompt": "2024",
    "2025": "2025",
}
SYSTEMATICS = ("nom", "up", "down")

# The provider's last pt bin runs to OPEN_PT_EDGE; it is cut at HIGH_PT_EDGE (as Wqq),
# and the bin above gets HIGH_PT_UNC_SCALE times the uncertainty.
OPEN_PT_EDGE = 13000.0
HIGH_PT_EDGE = 800.0
HIGH_PT_UNC_SCALE = 3.0
# The WP binning is [lo, hi): nudge the last edge so the last measured WP itself is valid.
LAST_WP_NUDGE = 1e-6
WP_DIGITS = 6


def read_provider_sf(path):
    """Return (wps, pt_edges, {wp_string: [per-pt-bin {syst: value}]}) at the measured WPs.

    The provider's SF is already piecewise linear in WP, so the value at a measured WP is
    the segment's line at its edge (the last WP from the last segment's upper end).
    """
    with open(path) as f:
        corr = next(c for c in json.load(f)["corrections"] if c["name"] == "Xbb_SF")
    wp_node = corr["data"]
    if wp_node["nodetype"] != "binning" or wp_node["input"] != "WP":
        raise ValueError(f"{path}: expected a WP binning at the top")
    edges = [float(e) for e in wp_node["edges"]]
    wps = edges[:-1] + [round(edges[-1] - LAST_WP_NUDGE, WP_DIGITS)]
    table, pt_edges = {}, None

    def cells_at(segment, wp):
        if segment["nodetype"] != "binning" or segment["input"] != "pt":
            raise ValueError(f"{path}: expected a pt binning under WP")
        cells = []
        for cat in segment["content"]:
            values = {}
            for c in cat["content"]:
                f = c["value"]
                if f["expression"] != "[0]+[1]*x" or f["variables"] != ["WP"]:
                    raise ValueError(f"{path}: unexpected formula {f['expression']} {f['variables']}")
                p0, p1 = f["parameters"]
                values[c["key"]] = p0 + p1 * wp
            cells.append(values)
        return [float(e) for e in segment["edges"]], cells

    for i, wp in enumerate(wps):
        seg = wp_node["content"][min(i, len(wp_node["content"]) - 1)]
        seg_edges, cells = cells_at(seg, wp)
        if pt_edges is None:
            pt_edges = seg_edges
        elif seg_edges != pt_edges:
            raise ValueError(f"{path}: pt edges differ between WP bins")
        # The segment below must end where this one starts (the SF is continuous in WP)
        if 0 < i < len(wps) - 1:
            _, below = cells_at(wp_node["content"][i - 1], wp)
            for a, b in zip(below, cells):
                if any(abs(a[s] - b[s]) > 1e-9 for s in SYSTEMATICS):
                    raise ValueError(f"{path}: SF not continuous at WP {wp}")
        table[f"{wp:.2f}"] = cells
    return [f"{w:.2f}" for w in wps], pt_edges, table


def read_provider_eff(path):
    """Return {era: (pt_edges, {wp_string: [per-pt-bin value]})}."""
    with open(path) as f:
        corr = next(c for c in json.load(f)["corrections"] if c["name"] == "Xbb_effMC")
    out = {}
    for era in corr["data"]["content"]:
        table, pt_edges = {}, None
        for wp in era["value"]["content"]:
            binning = wp["value"]
            if binning["nodetype"] != "binning" or binning["input"] != "pt":
                raise ValueError(f"{path}: expected a pt binning under era/WP")
            edges = [float(e) for e in binning["edges"]]
            if pt_edges is None:
                pt_edges = edges
            elif edges != pt_edges:
                raise ValueError(f"{path}: era {era['key']}: pt edges differ between WPs")
            table[wp["key"]] = [float(v) for v in binning["content"]]
        out[era["key"]] = (pt_edges, table)
    return out


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
    """Repeat the open-ended last bin above HIGH_PT_EDGE, with the inflated uncertainty for the SF."""
    last = cells[-1]
    if not with_systematic:
        return cells + [last]
    nom = last["nom"]
    return cells + [{
        "nom": nom,
        "up": nom + HIGH_PT_UNC_SCALE * (last["up"] - nom),
        "down": max(0.0, nom - HIGH_PT_UNC_SCALE * (nom - last["down"])),
    }]


def output_pt_edges(pt_edges):
    if pt_edges[-1] != OPEN_PT_EDGE or pt_edges[-2] >= HIGH_PT_EDGE:
        raise ValueError(f"expected an open-ended last pt bin below {HIGH_PT_EDGE:g} GeV, got {pt_edges}")
    return pt_edges[:-1] + [HIGH_PT_EDGE, OPEN_PT_EDGE]


def interpolated_node(era_table, pt_edges, with_systematic):
    """WP binning (piecewise linear) -> pt binning -> [systematic] -> formula."""
    wps = sorted(era_table, key=float)
    wp_values = [float(w) for w in wps]
    wp_edges = wp_values[:-1] + [wp_values[-1] + LAST_WP_NUDGE]
    out_pt_edges = output_pt_edges(pt_edges)
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


def sf_correction(era_table, pt_edges, year, provenance):
    data, wps = interpolated_node(era_table, pt_edges, with_systematic=True)
    return {
        "name": "Xbb_SF",
        "version": 1,
        "description": (
            f"GloParT-3 Xbb (XbbVsQCD) tagging SF, {year}, for true b jets with m_SD in [40, 250] GeV "
            "(mu-tagged method). Cumulative: the SF of score > WP, score = Xbb / (Xbb + QCD). "
            f"Piecewise linear in WP between the measured points {', '.join(wps)}; a WP outside "
            f"[{wps[0]}, {wps[-1]}] is an error. pt clamped; above {HIGH_PT_EDGE:g} GeV the last "
            f"measured bin's nominal with {HIGH_PT_UNC_SCALE:g}x its up/down distances. No eta "
            f"dependence. systematic: nom/up/down, absolute SF values. Source: {provenance}"
        ),
        "inputs": [
            {"name": "WP", "type": "real", "description": "Working point on the score (score > WP)"},
            {"name": "pt", "type": "real", "description": "Fat jet pt (nominal JEC/JER), GeV"},
            {"name": "systematic", "type": "string", "description": "nom, up or down"},
        ],
        "output": {"name": "sf", "type": "real"},
        "data": data,
    }


def eff_correction(era_table, pt_edges, provenance):
    data, wps = interpolated_node(era_table, pt_edges, with_systematic=False)
    return {
        "name": "Xbb_eff_MC",
        "version": 1,
        "description": (
            "Candle MC efficiency of true b jets to pass score > WP, for the fail SF "
            "(1 - SF*eff)/(1 - eff) and the multi-WP category weights. Piecewise linear in WP "
            f"between {', '.join(wps)}; pt clamped, above {HIGH_PT_EDGE:g} GeV the last measured "
            f"bin. No uncertainty. Source: {provenance}"
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
    parser.add_argument("--provider", required=True, help="folder with the files as delivered by the provider")
    parser.add_argument("--outdir", default=HERE, help="output directory (default: %(default)s)")
    args = parser.parse_args()

    eff_tables = read_provider_eff(os.path.join(args.provider, EFF_FILE))
    for year, era in ERAS.items():
        sf_file = f"Xbb_SF_{era}.json"
        wps, pt_edges, sf_table = read_provider_sf(os.path.join(args.provider, sf_file))
        if era not in eff_tables:
            raise KeyError(f"era {era!r} not in {EFF_FILE}")
        eff_pt_edges, eff_table = eff_tables[era]
        if eff_pt_edges != pt_edges:
            raise ValueError(f"era {era}: efficiency pt edges {eff_pt_edges} differ from the SF's {pt_edges}")
        if sorted(eff_table, key=float) != wps:
            raise ValueError(f"era {era}: efficiency working points {sorted(eff_table)} differ from the SF's {wps}")
        corrections = [
            sf_correction(sf_table, pt_edges, year, f"provider file {sf_file}"),
            eff_correction(eff_table, pt_edges, f"provider file {EFF_FILE}, era {era}"),
        ]
        out = {
            "schema_version": 2,
            "description": f"GloParT-3 Xbb tagging scale factors, {year}",
            "corrections": corrections,
        }
        path = os.path.join(args.outdir, f"Xbb_SF_{year}.json")
        with open(path, "w") as f:
            json.dump(out, f, indent=1)
        print(f"wrote {path} ({', '.join(c['name'] for c in corrections)})")


if __name__ == "__main__":
    main()
