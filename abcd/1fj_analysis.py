import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import os

import numpy as np
import yaml

from analysis import *

def main():
    parser = argparse.ArgumentParser(description="Run 1FJ Analysis")
    parser.add_argument("-i", "--input", required=True, help="Path to the input predictions CSV file")
    parser.add_argument("-d", "--data", required=False, help="Path to the data predictions CSV file (optional)")
    parser.add_argument("-n", "--num-scans", type=int, default=1, help="Number of orthogonal scans to perform")
    args = parser.parse_args()

    out_dir = os.path.dirname(args.input)
    if out_dir == "":
        out_dir = "."

    df = pd.read_csv(args.input)
    df = df.query("boosted_h_candidate_score >= 0.1 | boosted_v_candidate_score >= 0.1")
    
    if args.data:
        df_data = pd.read_csv(args.data)
        df_data = df_data.query("boosted_h_candidate_score >= 0.1 | boosted_v_candidate_score >= 0.1")
    else:
        df_data = None

    previous_cuts = []
    scan_cuts = {}
    
    total_sig_A, total_sig_B, total_sig_C, total_sig_D = 0, 0, 0, 0
    total_bkg_A, total_bkg_B, total_bkg_C, total_bkg_D = 0, 0, 0, 0
    total_data_B, total_data_C, total_data_D = 0, 0, 0

    for scan_idx in range(args.num_scans):
        print(f"\n{'='*50}")
        print(f"Starting Scan {scan_idx + 1}/{args.num_scans}")
        print(f"{'='*50}")

        df_sig = df[df["label"] == 1]
        df_bkg = df[df["label"] == 0]
        
        # Check if we have enough events left
        if len(df_sig) == 0 or len(df_bkg) == 0:
            print("Not enough signal or background events left for further optimization. Stopping.")
            break

        sig, (cut_var1_best, cut_var2_best, _, cut_dnn_best, vbs_cut_best) = optimize_cuts(
            df_sig,
            df_bkg,
            var1_col="boosted_h_candidate_score",
            var2_col="boosted_v_candidate_score",
            disco1_col="dnn_score",
            disco2_col="vbs_detajj",
            cut_var1_list=np.linspace(0.1, 1.0, 46),
            cut_var2_list=np.linspace(0.1, 1.0, 46),
            cut_disco1_list=np.linspace(0., 1.0, 51),
            cut_disco2_list=np.linspace(0., 10.0, 51),
            combine_method="or"
        )
        
        if sig == -np.inf:
            print("No valid cuts could be found. Stopping.")
            break

        print(f"Optimal Cuts (Scan {scan_idx + 1}):")
        if previous_cuts:
            print("  Orthogonality constraints applied:")
            for i, (prev_h, prev_v) in enumerate(previous_cuts):
                print(f"    - NOT (boosted_h_candidate_score >= {prev_h:.4f} OR boosted_v_candidate_score >= {prev_v:.4f})")
        
        print(f"  boosted_h_candidate_score > {cut_var1_best:.4f}")
        print(f"  boosted_v_candidate_score > {cut_var2_best:.4f}")
        print(f"  dnn_score                 > {cut_dnn_best:.4f}")
        print(f"  vbs_detajj                 > {vbs_cut_best:.4f}")
        print("-" * 50)
        print(f"Significance at optimal cuts: {sig:.4f}")
        print("="*50)
        
        a, b, c, d = get_ABCD_regions(df_sig, cut_var1_best, cut_var2_best, cut_dnn_best, vbs_cut_best, disco1_col="dnn_score", disco2_col="vbs_detajj", combine_method="or")
        e, f, g, h = get_ABCD_regions(df_bkg, cut_var1_best, cut_var2_best, cut_dnn_best, vbs_cut_best, disco1_col="dnn_score", disco2_col="vbs_detajj", combine_method="or")
        
        fig, ax = plt.subplots(figsize=(10, 8))
        sns.histplot(data=df_bkg.query("boosted_h_candidate_score >= @cut_var1_best | boosted_v_candidate_score >= @cut_var2_best"), x="dnn_score", y="vbs_detajj", weights="weight", bins=100, color='blue', cbar=True, ax=ax)
        
        # Add lines to mark ABCD regions
        ax.axvline(x=cut_dnn_best, color='red', linestyle='--', linewidth=2, label=f'dnn_score = {cut_dnn_best}')
        ax.axhline(y=vbs_cut_best, color='red', linestyle='--', linewidth=2, label=f'vbs_detajj = {vbs_cut_best}')
        
        # ax.legend()
        ax.set_title(f'Background: ABCD Regions (Scan {scan_idx + 1})')
        plt.savefig(os.path.join(out_dir, f'background_abcd_regions_scan_{scan_idx + 1}.png'))
        plt.close()
        
        # profile plot of bkg in dnn_score vs vbs_detajj, with color representing the mean weight in each bin
        fig, ax = plt.subplots(figsize=(10, 8))
        sns.histplot(data=df_bkg.query("boosted_h_candidate_score >= @cut_var1_best | boosted_v_candidate_score >= @cut_var2_best"), x="dnn_score", y="vbs_detajj", weights="weight", bins=100, color='blue', cbar=True, ax=ax)
        
        # profile of mean of vbs_detajj in bins of dnn_score
        bin_means = df_bkg.query("boosted_h_candidate_score >= @cut_var1_best | boosted_v_candidate_score >= @cut_var2_best").groupby(pd.cut(df_bkg["dnn_score"], bins=20), observed=False)["vbs_detajj"].mean()
        bin_centers = [interval.mid for interval in bin_means.index.categories]
        ax.plot(bin_centers, bin_means.values, color='red', marker='o', label='Mean vbs_detajj')
        ax.set_title(f'Background: dnn_score vs vbs_detajj with Mean Profile (Scan {scan_idx + 1})')
        plt.savefig(os.path.join(out_dir, f'background_dnn_vs_vbs_median_profile_scan_{scan_idx + 1}.png'))
        plt.close()
        
        fig, ax = plt.subplots(figsize=(10, 8))
        sns.histplot(data=df_sig.query("boosted_h_candidate_score >= @cut_var1_best | boosted_v_candidate_score >= @cut_var2_best"), x="dnn_score", y="vbs_detajj", weights="weight", bins=100, color='red', cbar=True, ax=ax)
        
        # Add lines to mark ABCD regions
        ax.axvline(x=cut_dnn_best, color='black', linestyle='--', linewidth=2, label=f'dnn_score = {cut_dnn_best}')
        ax.axhline(y=vbs_cut_best, color='black', linestyle='--', linewidth=2, label=f'vbs_detajj = {vbs_cut_best}')
        
        ax.set_title(f'Signal: ABCD Regions (Scan {scan_idx + 1})')
        plt.savefig(os.path.join(out_dir, f'signal_abcd_regions_scan_{scan_idx + 1}.png'))
        plt.close()
        
        print(f"\nSignal Regions Yields (Scan {scan_idx + 1}):\nA: {a:.4f}, B: {b:.4f}, C: {c:.4f}, D: {d:.4f}")
        print(f"\nBackground Regions Yields (Scan {scan_idx + 1}):\nA: {e:.4f}, B: {f:.4f}, C: {g:.4f}, D: {h:.4f}")
        
        total_sig_A += a; total_sig_B += b; total_sig_C += c; total_sig_D += d
        total_bkg_A += e; total_bkg_B += f; total_bkg_C += g; total_bkg_D += h

        if df_data is not None:
            _, j, k, l = get_ABCD_regions(df_data, cut_var1_best, cut_var2_best, cut_dnn_best, vbs_cut_best, disco1_col="dnn_score", disco2_col="vbs_detajj", combine_method="or")
            _ = -1 # blind
            
            total_data_B += j; total_data_C += k; total_data_D += l

            print(f"\nData Regions Yields (Scan {scan_idx + 1}):\nA: {_:.4f}, B: {j:.4f}, C: {k:.4f}, D: {l:.4f}")

        # Update dataframe for the next scan by removing events that fall in the current combined kinematic region
        # i.e., those that satisfied (boosted_h >= cut_var1_best) OR (boosted_v >= cut_var2_best)
        df = df[~((df["boosted_h_candidate_score"] >= cut_var1_best) | (df["boosted_v_candidate_score"] >= cut_var2_best))]
        if df_data is not None:
            df_data = df_data[~((df_data["boosted_h_candidate_score"] >= cut_var1_best) | (df_data["boosted_v_candidate_score"] >= cut_var2_best))]
        
        previous_cuts.append((cut_var1_best, cut_var2_best))
        scan_cuts[f"Scan{scan_idx + 1}"] = {
            "h_cut": round(float(cut_var1_best), 4),
            "v_cut": round(float(cut_var2_best), 4),
            "dnn_cut": round(float(cut_dnn_best), 4),
            "vbs_cut": round(float(vbs_cut_best), 4),
        }

    regions_path = os.path.join(out_dir, "regions.yaml")
    with open(regions_path, "w") as f:
        yaml.dump({"channels": scan_cuts}, f, default_flow_style=False, sort_keys=False)
    print(f"\nRegions written to {regions_path}")

    print(f"\n{'='*50}")
    print("Summary of Total Yields across all valid scans:")
    print(f"Total Signal Regions Yields:\nA: {total_sig_A:.4f}, B: {total_sig_B:.4f}, C: {total_sig_C:.4f}, D: {total_sig_D:.4f}")
    print(f"\nTotal Background Regions Yields:\nA: {total_bkg_A:.4f}, B: {total_bkg_B:.4f}, C: {total_bkg_C:.4f}, D: {total_bkg_D:.4f}")
    if df_data is not None:
        print(f"\nTotal Data Regions Yields:\nA: -1.0000, B: {total_data_B:.4f}, C: {total_data_C:.4f}, D: {total_data_D:.4f}")
    print("="*50)

if __name__ == "__main__":
    main()