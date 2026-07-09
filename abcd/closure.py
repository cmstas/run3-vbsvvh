#!/usr/bin/env python3

import argparse
import os

import numpy as np
import pandas as pd
import yaml

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Per-channel kinematic preselection: how the tagger scores are combined and
# which regions.yaml key holds each cut. Mirrors {1,2,3}fj_analysis.py.
CHANNELS = {
    "1fj": {"combine": "or",  "kin": [("boosted_h_candidate_score", "h_cut"),
                                       ("boosted_v_candidate_score", "v_cut")]},
    "2fj": {"combine": "and", "kin": [("boosted_h_candidate_score", "h_cut"),
                                       ("boosted_v_candidate_score", "v_cut")]},
    "3fj": {"combine": "and", "kin": [("boosted_h_candidate_score", "h_cut"),
                                       ("boosted_v1_candidate_score", "v1_cut"),
                                       ("boosted_v2_candidate_score", "v2_cut")]},
}

DNN_COL = "dnn_score"
BDT_COL = "bdt_score"


def kinematic_mask(df, channel, cuts):
    cfg = CHANNELS[channel]
    parts = [df[col].to_numpy() >= float(cuts[key]) for col, key in cfg["kin"]]
    mask = parts[0].copy()
    for p in parts[1:]:
        mask = (mask & p) if cfg["combine"] == "and" else (mask | p)
    return mask


def region_box(region, dnn_cut, vbs_cut, dnn_vals, bdt_vals):
    """
    Return (enddnn, endbdt): the upper edges of the blind box for the chosen
    region. "Full range" ends are taken from the data (max + tiny epsilon) so
    this works whether bdt_score is a [0,1] probability or a raw Delta-eta.
    """
    dnn_full = (float(dnn_vals.max()) + 1e-6) if dnn_vals.size else dnn_cut
    bdt_full = (float(bdt_vals.max()) + 1e-6) if bdt_vals.size else vbs_cut
    if region == "D":
        enddnn, endbdt = dnn_cut, vbs_cut
    elif region == "CD":
        enddnn, endbdt = dnn_full, vbs_cut
    elif region == "BD":
        enddnn, endbdt = dnn_cut, bdt_full
    else:
        raise ValueError(f"unknown region '{region}'")
    # Blinding guarantee: the box must exclude region A (dnn>=dnn_cut AND bdt>=vbs_cut).
    assert (enddnn <= dnn_cut) or (endbdt <= vbs_cut), (
        f"region '{region}' box does not exclude signal region A"
    )
    return enddnn, endbdt


def _quadrants(H):
    """
    Given a 2D histogram H[dnn_bin, bdt_bin], return cumulative mini-ABCD yields
    over every internal split (i, j):
        A[i,j] = sum over dnn-bin>=i, bdt-bin>=j
        B[i,j] = sum over dnn-bin>=i, bdt-bin<j
        C[i,j] = sum over dnn-bin<i,  bdt-bin>=j
        D[i,j] = sum over dnn-bin<i,  bdt-bin<j
    """
    A = H[::-1, ::-1].cumsum(0).cumsum(1)[::-1, ::-1]
    rowtot = A[:, [0]]   # dnn>=i, all bdt
    coltot = A[[0], :]   # all dnn, bdt>=j
    total = A[0, 0]
    B = rowtot - A
    C = coltot - A
    D = total - rowtot - coltot + A
    return A, B, C, D


def closure_scan(dnn, bdt, w, rangednn, rangebdt):
    """
    Weighted mini-ABCD grids and raw-count grids over the scan of internal
    splits. Counts are used for the min-statistics guard; weighted yields feed
    the closure ratio (identical for unweighted data).
    """
    Hw, _, _ = np.histogram2d(dnn, bdt, bins=[rangednn, rangebdt], weights=w)
    Hc, _, _ = np.histogram2d(dnn, bdt, bins=[rangednn, rangebdt])
    return _quadrants(Hw), _quadrants(Hc)


def plot_hist(errvec, title, path):
    if len(errvec) == 0:
        return
    mean = round(float(np.mean(errvec)), 3)
    std = round(float(np.std(errvec)), 3)
    plt.hist(errvec, 50, density=True, histtype="stepfilled", color="limegreen",
             label=f"mean={mean},\n std={std}")
    plt.legend()
    plt.title(title)
    plt.savefig(path)
    plt.clf()


def run_one(df_kin, dnn_cut, vbs_cut, region, n_scan, min_count, out_dir, tag):
    """Run the closure scan for one (scan, region) and return its summary dict."""
    dnn_all = df_kin[DNN_COL].to_numpy()
    bdt_all = df_kin[BDT_COL].to_numpy()

    enddnn, endbdt = region_box(region, dnn_cut, vbs_cut, dnn_all, bdt_all)

    # Blind box: keep only control-region events strictly inside it.
    box = (dnn_all < enddnn) & (bdt_all < endbdt) & (dnn_all >= 0) & (bdt_all >= 0)
    dnn = dnn_all[box]
    bdt = bdt_all[box]
    w = df_kin["weight"].to_numpy()[box] if "weight" in df_kin.columns else np.ones(box.sum())

    print(f"\n[{tag} | region {region}] box=({DNN_COL}<{enddnn:.4g}, {BDT_COL}<{endbdt:.4g}) "
          f"events={dnn.size}")
    if dnn.size == 0:
        print("  no events in box, skipping")
        return None

    rangednn = np.linspace(0.0, enddnn, n_scan)
    rangebdt = np.linspace(0.0, endbdt, n_scan)
    (Aw, Bw, Cw, Dw), (Ac, Bc, Cc, Dc) = closure_scan(dnn, bdt, w, rangednn, rangebdt)

    valid = (
        (Ac >= min_count) & (Bc >= min_count) & (Cc >= min_count) & (Dc >= min_count)
        & (Dw > 0)
    )
    Apred = np.divide(Bw * Cw, Dw, out=np.zeros_like(Dw), where=Dw > 0)
    valid &= Apred > 0

    diff_avg = np.divide(2.0 * (Aw - Apred), Aw + Apred,
                         out=np.zeros_like(Aw), where=(Aw + Apred) > 0)
    one_minus = np.divide(Aw, Apred, out=np.zeros_like(Aw), where=Apred > 0)
    one_minus = 1.0 - one_minus

    errs = diff_avg[valid]
    errs2 = one_minus[valid]
    abserrs = np.abs(errs)
    abserrs2 = np.abs(errs2)

    if errs.size == 0:
        print(f"  no valid split points (min_count={min_count}); try lowering --min-count")
        return None

    def ms(v):
        return float(np.mean(v)), float(np.std(v))

    m1, s1 = ms(errs)
    m2, s2 = ms(errs2)
    a1m, a1s = ms(abserrs)
    a2m, a2s = ms(abserrs2)
    print(f"  valid split points: {errs.size}")
    print(f"  diff/avg      2(A-Apred)/(A+Apred) : mean={m1:+.4f} std={s1:.4f}")
    print(f"  1 - obs/pred  1 - A/Apred          : mean={m2:+.4f} std={s2:.4f}")
    print(f"  |diff/avg|                          : mean={a1m:.4f} std={a1s:.4f}")
    print(f"  |1 - obs/pred|                      : mean={a2m:.4f} std={a2s:.4f}")

    base = f"closure_{tag}_{region}"
    plot_hist(errs,     f"{tag} {region}: 2(A-Apred)/(A+Apred)", os.path.join(out_dir, base + "_diff_over_avg.png"))
    plot_hist(errs2,    f"{tag} {region}: 1 - A/Apred",          os.path.join(out_dir, base + "_1_minus_ratio.png"))
    plot_hist(abserrs,  f"{tag} {region}: |2(A-Apred)/(A+Apred)|", os.path.join(out_dir, base + "_abs_diff_over_avg.png"))
    plot_hist(abserrs2, f"{tag} {region}: |1 - A/Apred|",        os.path.join(out_dir, base + "_abs_1_minus_ratio.png"))

    return {
        "n_points": int(errs.size),
        "box_events": int(dnn.size),
        "diff_over_avg": {"mean": round(m1, 4), "std": round(s1, 4)},
        "one_minus_ratio": {"mean": round(m2, 4), "std": round(s2, 4)},
        "abs_diff_over_avg": {"mean": round(a1m, 4), "std": round(a1s, 4)},
        "abs_one_minus_ratio": {"mean": round(a2m, 4), "std": round(a2s, 4)},
    }


def load_scan_cuts(args):
    """Return a dict {scan_name: cuts_dict} from --regions or the explicit CLI cuts."""
    if args.regions:
        with open(args.regions) as f:
            doc = yaml.safe_load(f) or {}
        channels = doc.get("channels", {})
        if not channels:
            raise SystemExit(f"No 'channels' found in {args.regions}")
        return channels
    if args.dnn_cut is None or args.vbs_cut is None:
        raise SystemExit("Provide --regions, or both --dnn-cut and --vbs-cut (plus kinematic cuts).")
    cuts = {"dnn_cut": args.dnn_cut, "vbs_cut": args.vbs_cut}
    for col, key in CHANNELS[args.channel]["kin"]:
        val = getattr(args, key.replace("_cut", "") + "_cut", None)
        cuts[key] = 0.0 if val is None else val
    return {"Scan1": cuts}


def main():
    global DNN_COL, BDT_COL
    parser = argparse.ArgumentParser(description="ABCD closure systematic from the data control region")
    parser.add_argument("-i", "--input", required=True, help="Path to the data predictions CSV")
    parser.add_argument("--channel", required=True, choices=list(CHANNELS), help="Analysis channel (sets kinematic preselection)")
    parser.add_argument("--regions", help="regions.yaml with per-scan cuts (as written by the analysis scripts)")
    parser.add_argument("--region", default="D", choices=["D", "CD", "BD", "all"], help="Closure region orientation (default: D)")
    parser.add_argument("-n", "--n-scan", type=int, default=50, help="Number of scan points per axis (default: 50)")
    parser.add_argument("--min-count", type=int, default=20, help="Minimum raw events required in each of A',B',C',D' (default: 1)")
    parser.add_argument("--orthogonal", action="store_true", help="With --regions, remove earlier scans' kinematic regions before each scan (matches the optimization)")
    parser.add_argument("-o", "--outdir", default=None, help="Output directory for plots/summary (default: alongside --input)")
    parser.add_argument("--dnn-col", default=DNN_COL, help=f"DNN score column (default: {DNN_COL})")
    parser.add_argument("--bdt-col", default=BDT_COL, help=f"Second-discriminant column (default: {BDT_COL}; older CSVs use 'vbs_score')")
    # Explicit-cut fallback (when --regions is not given)
    parser.add_argument("--dnn-cut", type=float, default=None)
    parser.add_argument("--vbs-cut", type=float, default=None)
    parser.add_argument("--h-cut", type=float, default=None)
    parser.add_argument("--v-cut", type=float, default=None)
    parser.add_argument("--v1-cut", type=float, default=None)
    parser.add_argument("--v2-cut", type=float, default=None)
    args = parser.parse_args()

    DNN_COL = args.dnn_col
    BDT_COL = args.bdt_col

    out_dir = args.outdir or (os.path.dirname(args.input) or ".")
    os.makedirs(out_dir, exist_ok=True)

    df = pd.read_csv(args.input)
    regions = ["D", "CD", "BD"] if args.region == "all" else [args.region]
    scan_cuts = load_scan_cuts(args)

    prev_kin = np.zeros(len(df), dtype=bool)
    summary = {}
    for scan_name, cuts in scan_cuts.items():
        this_kin = kinematic_mask(df, args.channel, cuts)
        eff_kin = this_kin & ~prev_kin if args.orthogonal else this_kin
        if args.orthogonal:
            prev_kin |= this_kin
        df_kin = df[eff_kin]
        dnn_cut, vbs_cut = float(cuts["dnn_cut"]), float(cuts["vbs_cut"])
        print(f"\n{'='*60}\n{scan_name}: kinematic events={len(df_kin)}  "
              f"dnn_cut={dnn_cut}  vbs_cut={vbs_cut}\n{'='*60}")

        summary[scan_name] = {}
        for region in regions:
            res = run_one(df_kin, dnn_cut, vbs_cut, region, args.n_scan,
                          args.min_count, out_dir, tag=f"{args.channel}_{scan_name}")
            if res is not None:
                summary[scan_name][region] = res

    summary_path = os.path.join(out_dir, f"closure_summary_{args.channel}.yaml")
    with open(summary_path, "w") as f:
        yaml.dump(summary, f, default_flow_style=False, sort_keys=False)
    print(f"\nClosure summary written to {summary_path}")


if __name__ == "__main__":
    main()
