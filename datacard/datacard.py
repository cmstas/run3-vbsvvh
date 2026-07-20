import os
import sys
import argparse
from pathlib import Path

import yaml
import pandas as pd
import numpy as np
from scipy.stats import gamma

# The systematics definition lives with the code that writes the predictions, so
# the branch list and nuisance naming cannot drift between the two.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "abcd"))
from systematics import ACCEPTANCE_ONLY, SYST_WEIGHTS, nuisance_name, ratio_columns


def get_poisson_uncertainty(n):
    """
    Calculate Poisson uncertainties for integer observation n using Gamma distribution.
    Matches the logic from the old script.
    """
    n = round(n)
    up = round(gamma.ppf(1 - (1 - 0.9973) / 2, n + 1))
    dn = round(gamma.ppf((1 - 0.9973) / 2, n)) if n > 0 else 0
    return up, dn

def get_yields_from_df(df, cuts):
    """
    Get the yield (sum of weights) from a pandas DataFrame applying the given cuts.
    """
    df_cut = df.query(cuts)
    return df_cut["weight"].sum()

def get_data_yields_from_df(df, cuts):
    """
    Get data yields from a pandas DataFrame applying the given cuts.
    Returns integer sum.
    """
    df_cut = df.query(cuts)
    return len(df_cut)

def available_systematics(df):
    """The systematic branches whose up/down ratio columns are present in ``df``."""
    present = []
    for branch in SYST_WEIGHTS:
        up_col, dn_col = ratio_columns(branch)
        if up_col in df.columns and dn_col in df.columns:
            present.append(branch)
    return present


def compute_syst_kappas(sig_df, region_dfs, proc_name, scan_name):
    """lnN kappas for every systematic, per ABCD region, for the signal.

    Only signal is affected: the background rate is a rateParam driven from the
    data control regions, so it carries no MC systematic.

    kappa is the varied-over-nominal yield ratio in the region. For the theory
    variations the inclusive ratio over the whole sample is divided out, leaving
    an acceptance-only effect (the signal cross section is the POI, so its
    normalization must not be constrained here too).

    Returns ``{nuisance: {region: (kappa_dn, kappa_up)}}``, skipping any
    systematic that is identically 1 in every region (e.g. muon SFs in a 0-lep
    channel) so no-op nuisances stay out of the card.
    """
    kappas = {}
    w_all = sig_df["weight"].to_numpy()
    total_nom = w_all.sum()

    for branch in available_systematics(sig_df):
        up_col, dn_col = ratio_columns(branch)

        # Inclusive normalization for the acceptance-only (theory) variations.
        norm_up = norm_dn = 1.0
        if branch in ACCEPTANCE_ONLY and total_nom > 0:
            norm_up = (w_all * sig_df[up_col].to_numpy()).sum() / total_nom
            norm_dn = (w_all * sig_df[dn_col].to_numpy()).sum() / total_nom

        per_region = {}
        for region_id, df_r in region_dfs.items():
            w = df_r["weight"].to_numpy()
            nom = w.sum()
            if nom <= 0 or len(df_r) == 0:
                per_region[region_id] = None
                continue
            k_up = (w * df_r[up_col].to_numpy()).sum() / nom
            k_dn = (w * df_r[dn_col].to_numpy()).sum() / nom
            if norm_up > 0:
                k_up /= norm_up
            if norm_dn > 0:
                k_dn /= norm_dn
            if not (np.isfinite(k_up) and np.isfinite(k_dn) and k_up > 0 and k_dn > 0):
                per_region[region_id] = None
                continue
            per_region[region_id] = (k_dn, k_up)

        values = [v for v in per_region.values() if v is not None]
        if not values:
            continue
        if all(abs(dn - 1) < 5e-5 and abs(up - 1) < 5e-5 for dn, up in values):
            continue

        kappas[nuisance_name(branch, proc_name, scan_name)] = per_region

    return kappas


def build_kin_terms(scan_info):
    """
    Build the list of kinematic score selections (>=) that define the phase space.
    Supports both the single-V phase space (v_cut / boosted_v_candidate_score, e.g. 1FJ/2FJ
    channels) and the two-V phase space (v1_cut, v2_cut / boosted_v{1,2}_candidate_score,
    e.g. the 3FJ channels). The ABCD axes (dnn_score, bdt_score) are handled separately.
    """
    terms = [f"boosted_h_candidate_score >= {scan_info['h_cut']}"]
    if "v1_cut" in scan_info or "v2_cut" in scan_info:
        terms.append(f"boosted_v1_candidate_score >= {scan_info['v1_cut']}")
        terms.append(f"boosted_v2_candidate_score >= {scan_info['v2_cut']}")
    else:
        terms.append(f"boosted_v_candidate_score >= {scan_info['v_cut']}")
    return terms

def create_abcd_datacard_single(process_name, out_name, scan_name, scan_info, orthogonality_cuts, sig_df, data_df, combination="and", unblind=False):
    """
    Create a single datacard for a specific scan.
    """

    yields = {scan_name: {}}
    observations = {scan_name: {}}
    poisson_errs = {scan_name: {}}

    # Base cuts
    dnn_cut = scan_info["dnn_cut"]
    vbs_cut = scan_info["vbs_cut"]

    # Kinematic phase space (h + V candidate scores), shared by all ABCD regions.
    # For single-fatjet channels (combination == "or") an event has either an H or a
    # V candidate (never both), so the candidate-score terms are OR'd; otherwise AND'd.
    kin_combine = " or " if combination == "or" else " and "
    kin_cut = "(" + kin_combine.join(build_kin_terms(scan_info)) + ")"

    A_cut = f"{kin_cut} and (dnn_score > {dnn_cut}) and (bdt_score > {vbs_cut})"
    B_cut = f"{kin_cut} and (dnn_score > {dnn_cut}) and (bdt_score < {vbs_cut})"
    C_cut = f"{kin_cut} and (dnn_score < {dnn_cut}) and (bdt_score > {vbs_cut})"
    D_cut = f"{kin_cut} and (dnn_score < {dnn_cut}) and (bdt_score < {vbs_cut})"
    
    cuts_dict = {"A": A_cut, "B": B_cut, "C": C_cut, "D": D_cut}

    ortho_str = ""
    if orthogonality_cuts:
        ortho_str = " and ".join([f"not ({cut})" for cut in orthogonality_cuts]) + " and "

    region_sig_dfs = {}
    for region_id in ["A", "B", "C", "D"]:
        final_cut = ortho_str + cuts_dict[region_id]
        
        # Data
        obs_val = get_data_yields_from_df(data_df, final_cut)
        if region_id == "A" and not unblind:
            observations[scan_name][region_id] = 1 # blind SR
            poisson_errs[scan_name][region_id] = (1, 1)
        else:
            observations[scan_name][region_id] = int(obs_val)
            poisson_errs[scan_name][region_id] = get_poisson_uncertainty(obs_val)
            
        # Signal
        # compute signal sum weights and sum weights^2 to get stat uncertainty
        df_sig_cut = sig_df.query(final_cut)
        region_sig_dfs[region_id] = df_sig_cut
        sig_sumw = df_sig_cut["weight"].sum()
        sig_sumw2 = (df_sig_cut["weight"] ** 2).sum()
        yields[scan_name][region_id] = sig_sumw
        # multiplicative lnN: 1 + sqrt(sumw2)/sumw (if sumw > 0)
        if sig_sumw > 0:
            sig_err = 1.0 + np.sqrt(sig_sumw2) / sig_sumw
        else:
            sig_err = 1.0
        # store per-region signal stat uncertainty
        poisson_errs[scan_name][f"sig_{region_id}"] = sig_err

    # Write datacard in old format (backgrounds first then signals)

    # if output dir doesn't exist then create it
    out_dir = os.path.dirname(out_name)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir)

    with open(out_name, "w+") as f:
        f.write("imax 4 number of channels\n")
        f.write("jmax 1 number of backgrounds\n")
        f.write("kmax * number of nuisance parameters\n")
        f.write("-" * 150 + "\n")

        # Observation block
        f.write(f"{'bin':<60}{'A':<20}{'B':<20}{'C':<20}{'D':<20}\n")
        f.write(f"{'observation':<60}{observations[scan_name]['A']:<20}{observations[scan_name]['B']:<20}{observations[scan_name]['C']:<20}{observations[scan_name]['D']:<20}\n")
        f.write("-" * 150 + "\n")

        # Process block
        bkg_name = f"TotalBkg_OneLep_{process_name}_{scan_name}"
        f.write(f"{'bin':<60}{'A':<50}{'B':<50}{'C':<50}{'D':<50}{'A':<20}{'B':<20}{'C':<20}{'D':<20}\n")
        f.write(f"{'process':<60}{bkg_name:<50}{bkg_name:<50}{bkg_name:<50}{bkg_name:<50}{'TotalSig':<20}{'TotalSig':<20}{'TotalSig':<20}{'TotalSig':<20}\n")
        f.write(f"{'process':<60}{'1':<50}{'1':<50}{'1':<50}{'1':<50}{'0':<20}{'0':<20}{'0':<20}{'0':<20}\n")
        f.write(f"{'rate':<60}{1:<50}{1:<50}{1:<50}{1:<50}{yields[scan_name]['A']:<20.5f}{yields[scan_name]['B']:<20.5f}{yields[scan_name]['C']:<20.5f}{yields[scan_name]['D']:<20.5f}\n")
        f.write("-" * 150 + "\n")

        # Systematics
        f.write(f"{f'CMS_{process_name}_{scan_name}_control_abcd_syst':<50}{'lnN':<10}{'1.35':<50}{'-':<50}{'-':<50}{'-':<50}{'-':<20}{'-':<20}{'-':<20}{'-':<20}\n")
        f.write(f"{'lumi_13TeV_correlated':<50}{'lnN':<10}{'-':<50}{'-':<50}{'-':<50}{'-':<50}{'1.016':<20}{'1.016':<20}{'1.016':<20}{'1.016':<20}\n")
        
        # Signal stat systematics
        sigA_err = poisson_errs[scan_name]['sig_A']
        sigB_err = poisson_errs[scan_name]['sig_B']
        sigC_err = poisson_errs[scan_name]['sig_C']
        sigD_err = poisson_errs[scan_name]['sig_D']
        f.write(f"{f'CMS_{process_name}_{scan_name}_signal_RegionA':<50}{'lnN':<10}{'-':<50}{'-':<50}{'-':<50}{'-':<50}{sigA_err:<20.5f}{'-':<20}{'-':<20}{'-':<20}\n")
        f.write(f"{f'CMS_{process_name}_{scan_name}_signal_RegionB':<50}{'lnN':<10}{'-':<50}{'-':<50}{'-':<50}{'-':<50}{'-':<20}{sigB_err:<20.5f}{'-':<20}{'-':<20}\n")
        f.write(f"{f'CMS_{process_name}_{scan_name}_signal_RegionC':<50}{'lnN':<10}{'-':<50}{'-':<50}{'-':<50}{'-':<50}{'-':<20}{'-':<20}{sigC_err:<20.5f}{'-':<20}\n")
        f.write(f"{f'CMS_{process_name}_{scan_name}_signal_RegionD':<50}{'lnN':<10}{'-':<50}{'-':<50}{'-':<50}{'-':<50}{'-':<20}{'-':<20}{'-':<20}{sigD_err:<20.5f}\n")

        # Weight-based systematics, signal only (the background rate is a
        # rateParam measured from the data control regions). Asymmetric lnN is
        # written in combine's "down/up" form.
        syst_kappas = compute_syst_kappas(sig_df, region_sig_dfs, process_name, scan_name)
        for nuisance, per_region in sorted(syst_kappas.items()):
            cells = []
            for region_id in ["A", "B", "C", "D"]:
                kv = per_region.get(region_id)
                cells.append("-" if kv is None else f"{kv[0]:.4f}/{kv[1]:.4f}")
            f.write(f"{nuisance:<50}{'lnN':<10}{'-':<50}{'-':<50}{'-':<50}{'-':<50}"
                    f"{cells[0]:<20}{cells[1]:<20}{cells[2]:<20}{cells[3]:<20}\n")
        if not syst_kappas:
            if not available_systematics(sig_df):
                print(f"  WARNING: {scan_name}: signal predictions carry no systematic weight "
                      f"columns (re-run inference to add them); card has only the hardcoded nuisances")
            else:
                print(f"  WARNING: {scan_name}: systematic columns present but no region has a "
                      f"usable signal yield; no weight systematics written")

        f.write("-" * 150 + "\n")

        # RateParams
        f.write(f"A_OneLep_{process_name}_{scan_name} rateParam       A                   {bkg_name}     (@0*@1/@2)\tB_OneLep_{process_name}_{scan_name},C_OneLep_{process_name}_{scan_name},D_OneLep_{process_name}_{scan_name}\n")
        for r in ['B', 'C', 'D']:
            o = observations[scan_name][r]
            up, dn = poisson_errs[scan_name][r]
            f.write(f"{r}_OneLep_{process_name}_{scan_name} rateParam       {r}                   {bkg_name}     {o}\t[{dn},{up}]\n")

    print(f"Datacard written to {out_name}")

def create_abcd_datacards(out_name_prefix, scans, sig_df, data_df, process_name, combination, unblind=False):
    """
    Create individual datacards for all scans using DataFrame inputs and cuts.
    """
    
    orthogonality_cuts = []
    
    for scan_name, scan_info in scans.items():
        out_name = f"{out_name_prefix}_{scan_name}.dat" if out_name_prefix else f"datacard_{scan_name}.dat"
        
        create_abcd_datacard_single(process_name, out_name, scan_name, scan_info, orthogonality_cuts, sig_df, data_df, combination, unblind)

        combine = " or " if combination == "or" else " and "
        ortho_cut = "(" + combine.join(build_kin_terms(scan_info)) + ")"
        orthogonality_cuts.append(ortho_cut)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate ABCD datacard from DataFrames")
    parser.add_argument("--sig", default="signal.parquet", help="Signal predictions file (.parquet)")
    parser.add_argument("--data", default="data.parquet", help="Data predictions file (.parquet)")
    parser.add_argument("--out", default="datacard.txt", help="Output datacard text file")
    parser.add_argument("--config", default="regions.yaml", help="Yaml config containing scans")
    parser.add_argument("--unblind", action="store_true", help="Unblind SR")
    parser.add_argument("--proc", default="vbsvvh1lep", help="Process name for datacard (used in systematics naming)")
    parser.add_argument("--combination", default="and", choices=["and", "or"], help="Use AND or OR combination of orthogonality cuts for multiple scans (default is AND)")
    args = parser.parse_args()
    
    if not os.path.exists(args.config):
        print(f"Config {args.config} not found. Ensure you created it.")
        exit(1)
        
    with open(args.config, 'r') as file:
        config = yaml.safe_load(file)
        
    scans = config.get("channels", {})
    
    # Open trees
    if os.path.exists(args.sig) and os.path.exists(args.data):
        sig_df = pd.read_parquet(args.sig)
        data_df = pd.read_parquet(args.data)

        # If the predictions contain both signal and background (label column), filter label==1
        if 'label' in sig_df.columns:
            sig_df = sig_df[sig_df['label'] == 1].copy()

        # Ensure weight column exists for signal; if not, assume weight=1
        if 'weight' not in sig_df.columns:
            sig_df['weight'] = 1.0

        create_abcd_datacards(args.out, scans, sig_df, data_df, args.proc, args.combination, args.unblind)
    else:
        print("Signal or Data files not found. Provide them via --sig and --data")

