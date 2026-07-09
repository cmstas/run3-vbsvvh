import os

import pandas as pd
import numpy as np
from joblib import Parallel, delayed

def to_threshold_bin_idx(vals, thresholds):
    """
    Pre-bin continuous values into index arrays using thresholds.
    """
    return np.searchsorted(thresholds, vals, side="right") - 1

def threshold_yields_2d_from_idx(d_idx, v_idx, w, n_dnn, n_vbs):
    """
    Computes cumulative 2D yields for >= threshold in both axes (region A).
    """
    if d_idx.size == 0:
        return np.zeros((n_dnn, n_vbs), dtype=float)

    flat_idx = d_idx * n_vbs + v_idx
    H_flat = np.bincount(flat_idx, weights=w, minlength=n_dnn * n_vbs)
    H = H_flat.reshape(n_dnn, n_vbs)

    # cumulative yields for >= threshold in both axes
    return H[::-1, ::-1].cumsum(axis=0).cumsum(axis=1)[::-1, ::-1]


def threshold_ABCD_control_from_idx(d_idx, v_idx, w, n_dnn, n_vbs):
    """
    Cumulative *control-region* yields for every (dnn threshold i, bdt threshold j).

    Splitting the (dnn_score, bdt_score) plane at (cut_disco1[i], cut_disco2[j])::

        C | A          A = dnn >= i AND bdt >= j   (signal region, NOT summed here)
        --+--          B = dnn >= i AND bdt <  j
        D | B          C = dnn <  i AND bdt >= j
                       D = dnn <  i AND bdt <  j

    Only the three control regions B, C, D are ever summed, so region A of the
    input (e.g. data) is never touched -> the signal region stays blind. The
    ABCD background estimate in A is then (B * C) / D.
    """
    z = np.zeros((n_dnn, n_vbs), dtype=float)
    if d_idx.size == 0:
        return z, z.copy(), z.copy()

    flat_idx = d_idx * n_vbs + v_idx
    H = np.bincount(flat_idx, weights=w, minlength=n_dnn * n_vbs).reshape(n_dnn, n_vbs)

    # Suffix/prefix (exclusive) partial sums along the dnn axis.
    suf0 = H[::-1].cumsum(axis=0)[::-1]   # sum over rows a >= i  (dnn >= i)
    pre0 = H.cumsum(axis=0) - H           # sum over rows a <  i  (dnn <  i)

    def suf1(M):                          # sum over cols b >= j  (bdt >= j)
        return M[:, ::-1].cumsum(axis=1)[:, ::-1]

    def pre1(M):                          # sum over cols b <  j  (bdt <  j)
        return M.cumsum(axis=1) - M

    B = pre1(suf0)   # dnn >= i, bdt <  j
    C = suf1(pre0)   # dnn <  i, bdt >= j
    D = pre1(pre0)   # dnn <  i, bdt <  j
    return B, C, D

def optimize_cuts(df_sig, df_bkg,
                  var1_col="boosted_h_candidate_score", var2_col="boosted_v_candidate_score", var3_col=None,
                  disco1_col="dnn_0_score", disco2_col="dnn_1_score",
                  cut_var1_list=None, cut_var2_list=None, cut_var3_list=None,
                  cut_disco1_list=None, cut_disco2_list=None,
                  combine_method="and", n_jobs=16,
                  df_data=None, use_data_driven_bkg=False, min_data_yield=0.0):
    """
    Optimizes 5D cuts to maximize significance.

    The background yield in region A can come from two sources:
      - MC background (default): Y_bkg is the weighted MC yield in region A.
      - Data-driven ABCD prediction (use_data_driven_bkg=True): Y_bkg is
        estimated as (B * C) / D from the *data* control regions. Region A of
        the data is never summed, so the signal region stays blind. df_data
        must be provided. `min_data_yield` requires each of B, C, D to hold at
        least that (weighted) yield for a cell to be eligible, which suppresses
        cuts driven by empty/low-stat control regions.
    """
    if cut_var1_list is None:
        cut_var1_list = np.linspace(0, 1, 21)
    if cut_var2_list is None:
        cut_var2_list = np.linspace(0, 1, 21)
    if cut_var3_list is None:
        if var3_col is not None:
            cut_var3_list = np.linspace(0, 1, 21)
        else:
            cut_var3_list = [None]
    if cut_disco1_list is None:
        cut_disco1_list = np.linspace(0, 1, 101)
    if cut_disco2_list is None:
        cut_disco2_list = np.linspace(0, 1, 101)

    n_dnn = len(cut_disco1_list)
    n_vbs = len(cut_disco2_list)

    if use_data_driven_bkg and df_data is None:
        raise ValueError("use_data_driven_bkg=True requires df_data to be provided")

    # The background-yield source (MC df_bkg or data control regions df_data).
    df_ctrl = df_data if use_data_driven_bkg else df_bkg

    sig_h = df_sig[var1_col].to_numpy()
    sig_v = df_sig[var2_col].to_numpy()
    if var3_col is not None:
        sig_3 = df_sig[var3_col].to_numpy()
    sig_d = df_sig[disco1_col].to_numpy()
    sig_j = df_sig[disco2_col].to_numpy()
    sig_w = df_sig["weight"].to_numpy()

    ctrl_h = df_ctrl[var1_col].to_numpy()
    ctrl_v = df_ctrl[var2_col].to_numpy()
    if var3_col is not None:
        ctrl_3 = df_ctrl[var3_col].to_numpy()
    ctrl_d = df_ctrl[disco1_col].to_numpy()
    ctrl_j = df_ctrl[disco2_col].to_numpy()
    ctrl_w = df_ctrl["weight"].to_numpy()

    # Precompute threshold comparisons once
    sig_h_ge = sig_h[:, None] >= cut_var1_list[None, :]
    sig_v_ge = sig_v[:, None] >= cut_var2_list[None, :]
    if var3_col is not None:
        sig_3_ge = sig_3[:, None] >= cut_var3_list[None, :]

    ctrl_h_ge = ctrl_h[:, None] >= cut_var1_list[None, :]
    ctrl_v_ge = ctrl_v[:, None] >= cut_var2_list[None, :]
    if var3_col is not None:
        ctrl_3_ge = ctrl_3[:, None] >= cut_var3_list[None, :]

    sig_d_idx = to_threshold_bin_idx(sig_d, cut_disco1_list)
    sig_j_idx = to_threshold_bin_idx(sig_j, cut_disco2_list)
    ctrl_d_idx = to_threshold_bin_idx(ctrl_d, cut_disco1_list)
    ctrl_j_idx = to_threshold_bin_idx(ctrl_j, cut_disco2_list)

    sig_valid = (sig_d_idx >= 0) & (sig_j_idx >= 0)
    ctrl_valid = (ctrl_d_idx >= 0) & (ctrl_j_idx >= 0)

    def evaluate_cut_var1(i, cut_var1):
        local_best_significance = -np.inf
        local_best_tuple = None

        sig_mask_h = sig_h_ge[:, i]
        ctrl_mask_h = ctrl_h_ge[:, i]

        for j, cut_var2 in enumerate(cut_var2_list):
            sig_mask_hv = sig_mask_h & sig_v_ge[:, j] if combine_method == "and" else sig_mask_h | sig_v_ge[:, j]
            ctrl_mask_hv = ctrl_mask_h & ctrl_v_ge[:, j] if combine_method == "and" else ctrl_mask_h | ctrl_v_ge[:, j]

            for k, cut_var3 in enumerate(cut_var3_list):
                if var3_col is not None:
                    if combine_method == "and":
                        sig_mask = (sig_mask_hv & sig_3_ge[:, k]) & sig_valid
                        ctrl_mask = (ctrl_mask_hv & ctrl_3_ge[:, k]) & ctrl_valid
                    elif combine_method == "or":
                        sig_mask = (sig_mask_hv | sig_3_ge[:, k]) & sig_valid
                        ctrl_mask = (ctrl_mask_hv | ctrl_3_ge[:, k]) & ctrl_valid
                    else:
                        raise ValueError("combine_method must be 'and' or 'or'")
                else:
                    sig_mask = sig_mask_hv & sig_valid
                    ctrl_mask = ctrl_mask_hv & ctrl_valid

                if not sig_mask.any() and not ctrl_mask.any():
                    continue

                Y_sig = threshold_yields_2d_from_idx(
                    sig_d_idx[sig_mask], sig_j_idx[sig_mask], sig_w[sig_mask], n_dnn, n_vbs
                )

                if use_data_driven_bkg:
                    # Data-driven ABCD estimate of the background in region A,
                    # (B * C) / D, using only the data control regions.
                    B, C, D = threshold_ABCD_control_from_idx(
                        ctrl_d_idx[ctrl_mask], ctrl_j_idx[ctrl_mask], ctrl_w[ctrl_mask], n_dnn, n_vbs
                    )
                    Y_bkg = np.divide(B * C, D, out=np.zeros_like(D), where=D > 0)
                    valid_cells = (
                        (D > 0)
                        & (Y_bkg > 0)
                        & (B >= min_data_yield)
                        & (C >= min_data_yield)
                        & (D >= min_data_yield)
                    )
                else:
                    Y_bkg = threshold_yields_2d_from_idx(
                        ctrl_d_idx[ctrl_mask], ctrl_j_idx[ctrl_mask], ctrl_w[ctrl_mask], n_dnn, n_vbs
                    )
                    valid_cells = None

                s_over_b = np.divide(Y_sig, Y_bkg, out=np.zeros_like(Y_sig), where=Y_bkg > 0)
                term = np.where(
                    Y_bkg > 0,
                    (Y_sig + Y_bkg) * np.log1p(s_over_b) - Y_sig,
                    0.0,
                )
                S = np.sqrt(np.maximum(2.0 * term, 0.0))
                if valid_cells is not None:
                    # Only cells with a usable data-driven background estimate
                    # may be selected; everything else is disqualified.
                    S = np.where(valid_cells, S, -np.inf)

                idx = np.unravel_index(np.argmax(S), S.shape)
                current_significance = S[idx]
                if current_significance > local_best_significance:
                    local_best_significance = current_significance
                    local_best_tuple = (cut_var1, cut_var2, cut_var3, cut_disco1_list[idx[0]], cut_disco2_list[idx[1]])

        return local_best_significance, local_best_tuple

    results = Parallel(n_jobs=n_jobs)(
        delayed(evaluate_cut_var1)(i, cut_var1) for i, cut_var1 in enumerate(cut_var1_list)
    )

    best_significance = -np.inf
    best_tuple = None
    for sig, tup in results:
        if sig > best_significance:
            best_significance = sig
            best_tuple = tup

    if best_tuple is None:
        # Default or fallback if none is found
        return -np.inf, (None, None, None, None, None)

    return best_significance, best_tuple

def get_ABCD_regions(df, cut_var1, cut_var2, cut_disco1, cut_disco2,
                     var1_col="boosted_h_candidate_score", var2_col="boosted_v_candidate_score",
                     disco1_col="dnn_0_score", disco2_col="dnn_1_score",
                     combine_method="and",
                     var3_col=None, cut_var3=None):
    """
    Computes ABCD region yields given the optimized cut thresholds.

    Optionally applies a third tagger cut (var3_col >= cut_var3), combined with
    the first two using the same combine_method.
    """
    # Apply initial kinematic cuts if applicable
    if combine_method == "and":
        kin_mask = (df[var1_col] >= cut_var1) & (df[var2_col] >= cut_var2)
        if var3_col is not None:
            kin_mask = kin_mask & (df[var3_col] >= cut_var3)
    elif combine_method == "or":
        kin_mask = (df[var1_col] >= cut_var1) | (df[var2_col] >= cut_var2)
        if var3_col is not None:
            kin_mask = kin_mask | (df[var3_col] >= cut_var3)
    else:
        raise ValueError("combine_method must be 'and' or 'or'")

    df_filtered = df[kin_mask]
    
    A = df_filtered[(df_filtered[disco1_col] >= cut_disco1) & (df_filtered[disco2_col] >= cut_disco2)].weight.sum()
    B = df_filtered[(df_filtered[disco1_col] >= cut_disco1) & (df_filtered[disco2_col] < cut_disco2)].weight.sum()
    C = df_filtered[(df_filtered[disco1_col] < cut_disco1) & (df_filtered[disco2_col] >= cut_disco2)].weight.sum()
    D = df_filtered[(df_filtered[disco1_col] < cut_disco1) & (df_filtered[disco2_col] < cut_disco2)].weight.sum()

    return A, B, C, D


def plot_background_decorrelation(df_bkg_kin, cut_disco1, cut_disco2, out_dir, scan_idx,
                                  disco1_col="dnn_score", disco2_col="bdt_score",
                                  df_data_kin=None, use_data_driven=False,
                                  profile_name="background_dnn_vs_vbs_profile"):
    """
    Draw the background ABCD 2D map and the mean-<disco2> vs <disco1> profile that
    is used to eyeball whether the two discriminants are decorrelated.

    - MC mode (use_data_driven=False): draws df_bkg_kin (weighted MC), full plane.
    - Data-driven mode (use_data_driven=True): draws df_data_kin and stays BLIND to
      the signal region:
        * the 2D map drops region A (disco1 >= cut AND disco2 >= cut), so the
          top-right quadrant is left blank -- only the control regions B, C, D
          are shown;
        * the profile is computed only in the disco2 < cut_disco2 control band
          (regions B + D), which never contains region A.

    df_*_kin frames must already have the kinematic (tagger) preselection applied.
    Writes: background_abcd_regions_scan_<N>.png and <profile_name>_scan_<N>.png.
    """
    import seaborn as sns
    import matplotlib.pyplot as plt

    if use_data_driven:
        if df_data_kin is None:
            raise ValueError("use_data_driven=True requires df_data_kin")
        in_A = (df_data_kin[disco1_col] >= cut_disco1) & (df_data_kin[disco2_col] >= cut_disco2)
        map_df = df_data_kin[~in_A]                                 # blind: no region A
        prof_df = df_data_kin[df_data_kin[disco2_col] < cut_disco2]  # control band B + D
        weights = None                                              # data is unweighted
        color = "green"
        src_label = "Data"
        map_note = " (region A blinded)"
        prof_note = f" [profile: {disco2_col} < {cut_disco2:.3f} control band]"
    else:
        map_df = df_bkg_kin
        prof_df = df_bkg_kin
        weights = "weight"
        color = "blue"
        src_label = "Background MC"
        map_note = ""
        prof_note = ""

    # 2D ABCD map
    fig, ax = plt.subplots(figsize=(10, 8))
    sns.histplot(data=map_df, x=disco1_col, y=disco2_col, weights=weights,
                 bins=100, color=color, cbar=True, ax=ax)
    ax.axvline(x=cut_disco1, color="red", linestyle="--", linewidth=2, label=f"{disco1_col} = {cut_disco1}")
    ax.axhline(y=cut_disco2, color="red", linestyle="--", linewidth=2, label=f"{disco2_col} = {cut_disco2}")
    ax.set_title(f"{src_label}: ABCD Regions{map_note} (Scan {scan_idx + 1})")
    plt.savefig(os.path.join(out_dir, f"background_abcd_regions_scan_{scan_idx + 1}.png"))
    plt.close()

    # profile of mean disco2 in bins of disco1
    fig, ax = plt.subplots(figsize=(10, 8))
    sns.histplot(data=map_df, x=disco1_col, y=disco2_col, weights=weights,
                 bins=100, color=color, cbar=True, ax=ax)
    bin_means = prof_df.groupby(pd.cut(prof_df[disco1_col], bins=20), observed=False)[disco2_col].mean()
    bin_centers = [interval.mid for interval in bin_means.index.categories]
    ax.plot(bin_centers, bin_means.values, color="red", marker="o", label=f"Mean {disco2_col}")
    ax.set_title(f"{src_label}: {disco1_col} vs {disco2_col} Mean Profile{prof_note} (Scan {scan_idx + 1})")
    plt.savefig(os.path.join(out_dir, f"{profile_name}_scan_{scan_idx + 1}.png"))
    plt.close()
