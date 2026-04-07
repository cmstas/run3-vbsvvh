#!/usr/bin/env python3
"""
Script to extract and sum cutflow tables from job stdout files.
"""

import os
import re
import glob
from collections import defaultdict
from pathlib import Path
from argparse import ArgumentParser


def extract_cutflow_table(file_path):
    """
    Extract cutflow table from a job stdout file.
    Returns a list of tuples: [(cut_number, cut_name, sum_w, rel_eff, abs_eff), ...]
    """
    cutflow_data = []
    
    try:
        with open(file_path, 'r') as f:
            content = f.read()
            
        # Find the cutflow section - look for the pattern between the table headers
        # Match lines like: | 0 | All events                                   | 1221175.000 |     1.000 |     1.000 |
        pattern = r'\|\s+(\d+)\s+\|\s+(.+?)\s+\|\s+([\d.]+)\s+\|\s+([\d.]+)\s+\|\s+([\d.]+)\s+\|'
        
        matches = re.findall(pattern, content)
        
        for match in matches:
            cut_num = int(match[0])
            cut_name = match[1].strip()
            sum_w = float(match[2])
            rel_eff = float(match[3])
            abs_eff = float(match[4])
            
            cutflow_data.append({
                'cut_num': cut_num,
                'cut_name': cut_name,
                'sum_w': sum_w,
                'rel_eff': rel_eff,
                'abs_eff': abs_eff,
            })
    
    except (FileNotFoundError, IOError) as e:
        print(f"Warning: Could not read {file_path}: {e}")
    
    return cutflow_data


def aggregate_cutflows(job_stdout_files):
    """
    Aggregate cutflow data from multiple job files.
    Returns a dictionary with cut numbers as keys and aggregated data as values.
    """
    aggregated = defaultdict(lambda: {
        'cut_name': None,
        'sum_w': 0.0,
        'count': 0,
        'rel_eff_sum': 0.0,
        'abs_eff_sum': 0.0,
    })
    
    file_count = 0
    for file_path in job_stdout_files:
        cutflow = extract_cutflow_table(file_path)
        if cutflow:
            file_count += 1
            for entry in cutflow:
                cut_num = entry['cut_num']
                aggregated[cut_num]['cut_name'] = entry['cut_name']
                aggregated[cut_num]['sum_w'] += entry['sum_w']
                aggregated[cut_num]['rel_eff_sum'] += entry['rel_eff']
                aggregated[cut_num]['abs_eff_sum'] += entry['abs_eff']
                aggregated[cut_num]['count'] += 1
    
    return aggregated, file_count


def print_cutflow_table(aggregated, file_count):
    """
    Print the aggregated cutflow table in a nice format.
    """
    if not aggregated:
        print("No cutflow data found!")
        return
    
    # Sort by cut number
    sorted_cuts = sorted(aggregated.items(), key=lambda x: x[0])
    
    # Print header
    print("\n" + "="*120)
    print(f"AGGREGATED CUTFLOW TABLE (from {file_count} job files)")
    print("="*120)
    print(f"{'#':<3} | {'Cut name':<50} | {'Sum(w)':<15} | {'Avg rel. eff.':<15} | {'Avg abs. eff.':<15}")
    print("-"*120)
    
    # Print rows
    for cut_num, data in sorted_cuts:
        avg_rel_eff = data['rel_eff_sum'] / data['count'] if data['count'] > 0 else 0.0
        avg_abs_eff = data['abs_eff_sum'] / data['count'] if data['count'] > 0 else 0.0
        
        print(f"{cut_num:<3} | {data['cut_name']:<50} | {data['sum_w']:<15.1f} | {avg_rel_eff:<15.4f} | {avg_abs_eff:<15.4f}")
    
    print("="*120 + "\n")


def main():
    parser = ArgumentParser(description="Aggregate cutflow tables from multiple job files")
    parser.add_argument("-p", "--path", default="/home/users/aaarora/phys/run3/cmstas-run3-vbsvvh/preselection/condor/jobs/1Lep2FJ_run2-data_1lep_2FJ_r2_2fj_data/", help="Base path to job stdout files")
    args = parser.parse_args()

    # Find all job stdout files
    base_path = args.path
    job_files = glob.glob(f"{base_path}/**/job.*.stdout", recursive=True)
    
    print(f"Found {len(job_files)} job stdout files")
    
    if not job_files:
        print("No job stdout files found!")
        return
    
    # Aggregate cutflows
    aggregated, file_count = aggregate_cutflows(job_files)
    
    print(f"Successfully extracted cutflows from {file_count} files")
    
    # Print aggregated table
    print_cutflow_table(aggregated, file_count)
    
if __name__ == '__main__':
    main()
