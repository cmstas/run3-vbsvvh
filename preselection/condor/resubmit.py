#!/usr/bin/env python3
"""
Resubmit failed or held condor jobs for run3-vbsvvh preselection framework.

Usage:
    # Resubmit all failed jobs
    python condor/resubmit.py --task 0Lep2FJ_run2-sig_0Lep2FJ_v1 --failed

    # Resubmit held jobs
    python condor/resubmit.py --task 0Lep2FJ_run2-sig_0Lep2FJ_v1 --held

    # Resubmit specific jobs by ID
    python condor/resubmit.py --task 0Lep2FJ_run2-sig_0Lep2FJ_v1 --jobs job1 job2

    # Dry run (show what would be resubmitted)
    python condor/resubmit.py --task 0Lep2FJ_run2-sig_0Lep2FJ_v1 --failed --dry-run
"""

import argparse
import json
import os
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

# Import status module for job status determination
from status import (
    get_condor_jobs,
    get_condor_history,
    determine_job_status,
    load_manifest,
    save_manifest
)


CONDOR_OUTPUT_DIR = "jobs"


def extract_cluster_id(condor_output: str) -> Optional[int]:
    """Extract cluster ID from condor_submit output."""
    match = re.search(r'submitted to cluster (\d+)', condor_output)
    if match:
        return int(match.group(1))
    return None


def resubmit_job(job_id: str, job: dict, dry_run: bool = False) -> Optional[int]:
    """
    Resubmit a single job.

    Returns new cluster_id on success, None on failure.
    """
    job_dir = Path(job.get("job_dir", ""))
    submit_file = job_dir / "submit.cmd"

    if not submit_file.exists():
        print(f"  ERROR: Submit file not found: {submit_file}")
        return None

    if dry_run:
        print(f"  [DRY RUN] Would resubmit: {job_id}")
        return -1  # Placeholder

    result = subprocess.run(
        ["condor_submit", str(submit_file)],
        capture_output=True, text=True
    )

    if result.returncode != 0:
        print(f"  ERROR submitting {job_id}: {result.stderr}")
        return None

    cluster_id = extract_cluster_id(result.stdout)
    if cluster_id:
        print(f"  {job_id}: resubmitted to cluster {cluster_id}")
        return cluster_id
    else:
        print(f"  {job_id}: resubmitted (couldn't parse cluster ID)")
        return -1


def main():
    parser = argparse.ArgumentParser(
        description="Resubmit failed condor jobs",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--task", "-t", required=True,
                        help="Task name (directory name under condor/jobs/)")
    parser.add_argument("--failed", "-f", action="store_true",
                        help="Resubmit all failed jobs")
    parser.add_argument("--held", action="store_true",
                        help="Resubmit all held jobs")
    parser.add_argument("--prepared", "-p", action="store_true",
                        help="Submit all prepared jobs (from dry-run)")
    parser.add_argument("--jobs", "-j", nargs="+",
                        help="Specific job IDs to resubmit")
    parser.add_argument("--sample", "-s",
                        help="Resubmit failed jobs for specific sample (regex)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show what would be resubmitted without submitting")
    args = parser.parse_args()

    if not (args.failed or args.held or args.prepared or args.jobs):
        parser.error("Must specify --failed, --held, --prepared, or --jobs")

    # Determine paths
    script_dir = Path(__file__).parent.resolve()
    condor_dir = script_dir
    task_dir = condor_dir / CONDOR_OUTPUT_DIR / args.task

    if not task_dir.exists():
        print(f"ERROR: Task not found: {task_dir}")
        sys.exit(1)

    # Load manifest
    manifest = load_manifest(task_dir)
    if not manifest:
        print(f"ERROR: No manifest.json found in {task_dir}")
        sys.exit(1)

    jobs = manifest.get("jobs", {})
    if not jobs:
        print("No jobs in manifest")
        return

    print(f"Task: {args.task}")
    print(f"Total jobs in manifest: {len(jobs)}")

    # Determine current status of all jobs
    print("\nQuerying condor for current status...")
    condor_jobs, err = get_condor_jobs()
    if err:
        print(f"WARNING: {err}")

    cluster_ids = [
        j["cluster_id"] for j in jobs.values()
        if j.get("cluster_id") is not None
    ]
    missing_ids = [cid for cid in cluster_ids if cid not in condor_jobs]
    condor_history = {}
    if missing_ids:
        condor_history, err = get_condor_history(missing_ids)
        if err:
            print(f"WARNING: {err}")

    # Determine which jobs to resubmit
    jobs_to_resubmit = []

    if args.jobs:
        # Specific jobs
        for job_id in args.jobs:
            if job_id in jobs:
                jobs_to_resubmit.append((job_id, jobs[job_id]))
            else:
                print(f"WARNING: Job ID not found: {job_id}")
    else:
        # By status
        for job_id, job in jobs.items():
            status = determine_job_status(job, condor_jobs, condor_history)

            # Filter by sample if specified
            if args.sample:
                sample = job.get("sample", "")
                if not re.search(args.sample, sample):
                    continue

            if args.failed and status == "failed":
                jobs_to_resubmit.append((job_id, job))
            elif args.held and status == "held":
                jobs_to_resubmit.append((job_id, job))
            elif args.prepared and status == "prepared":
                jobs_to_resubmit.append((job_id, job))

    if not jobs_to_resubmit:
        print("\nNo jobs to resubmit.")
        return

    print(f"\nJobs to resubmit: {len(jobs_to_resubmit)}")

    if args.dry_run:
        print("\n[DRY RUN MODE]")

    # Resubmit jobs
    print("\nResubmitting...")
    success_count = 0
    fail_count = 0

    for job_id, job in jobs_to_resubmit:
        new_cluster_id = resubmit_job(job_id, job, args.dry_run)

        if new_cluster_id is not None:
            success_count += 1
            if not args.dry_run and new_cluster_id > 0:
                # Update manifest
                manifest["jobs"][job_id]["cluster_id"] = new_cluster_id
                manifest["jobs"][job_id]["status"] = "submitted"
                manifest["jobs"][job_id]["submit_time"] = datetime.now().isoformat()
        else:
            fail_count += 1

    # Save updated manifest
    if not args.dry_run:
        manifest["last_resubmit"] = datetime.now().isoformat()
        save_manifest(task_dir, manifest)

    print(f"\n{'='*40}")
    print(f"Resubmitted: {success_count}, Failed: {fail_count}")
    if not args.dry_run:
        print(f"Manifest updated: {task_dir}/manifest.json")


if __name__ == "__main__":
    main()
