#!/usr/bin/env python3
"""
Resubmit failed or stuck SLURM jobs for run3-vbsvvh preselection framework.

Usage:
    # Resubmit all failed jobs
    python slurm/resubmit.py --task <name> --failed

    # Resubmit with more memory (e.g. for OOM failures)
    python slurm/resubmit.py --task <name> --failed --memory 8gb

    # Resubmit specific jobs by ID
    python slurm/resubmit.py --task <name> --jobs job1 job2

    # Resubmit failed jobs for a specific sample
    python slurm/resubmit.py --task <name> --failed --sample "TTToSemiLeptonic"

    # Resubmit stuck jobs (cancels them first, checks node health via SSH)
    python slurm/resubmit.py --task <name> --stuck

    # Dry run (show what would be resubmitted)
    python slurm/resubmit.py --task <name> --failed --dry-run
"""

import argparse
import json
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, Optional

from status import (
    get_slurm_job_states,
    determine_job_status,
    detect_stuck_jobs,
    load_manifest,
    save_manifest
)


SLURM_OUTPUT_DIR = "jobs"


def extract_slurm_job_id(sbatch_output: str) -> Optional[int]:
    """Extract job ID from sbatch output (e.g., 'Submitted batch job 12345')."""
    match = re.search(r'Submitted batch job (\d+)', sbatch_output)
    if match:
        return int(match.group(1))
    return None


def update_sbatch_memory(sbatch_file: Path, memory: str):
    """Update the --mem directive in a .sbatch file."""
    content = sbatch_file.read_text()
    new_content = re.sub(r'#SBATCH --mem=\S+', f'#SBATCH --mem={memory}', content)
    sbatch_file.write_text(new_content)


def resubmit_job(job_id: str, job: dict, memory: Optional[str] = None,
                 dry_run: bool = False) -> Optional[int]:
    """
    Resubmit a single job, optionally with updated memory.

    Returns new slurm_job_id on success, None on failure.
    """
    job_dir = Path(job.get("job_dir", ""))
    sbatch_file = job_dir / "job.sbatch"

    if not sbatch_file.exists():
        print(f"  ERROR: sbatch file not found: {sbatch_file}")
        return None

    # Update memory in sbatch file if requested
    if memory:
        update_sbatch_memory(sbatch_file, memory)

    if dry_run:
        mem_note = f" (mem={memory})" if memory else ""
        print(f"  [DRY RUN] Would resubmit: {job_id}{mem_note}")
        return -1  # Placeholder

    result = subprocess.run(
        ["sbatch", str(sbatch_file)],
        capture_output=True, text=True
    )

    if result.returncode != 0:
        print(f"  ERROR submitting {job_id}: {result.stderr.strip()}")
        return None

    slurm_job_id = extract_slurm_job_id(result.stdout)
    if slurm_job_id:
        mem_note = f", mem={memory}" if memory else ""
        print(f"  {job_id}: resubmitted (job {slurm_job_id}{mem_note})")
        return slurm_job_id
    else:
        print(f"  {job_id}: resubmitted (couldn't parse job ID)")
        return -1


def main():
    parser = argparse.ArgumentParser(
        description="Resubmit failed SLURM jobs",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--task", "-t", required=True,
                        help="Task name (directory name under slurm/jobs/)")
    parser.add_argument("--failed", "-f", action="store_true",
                        help="Resubmit all failed jobs")
    parser.add_argument("--stuck", action="store_true",
                        help="Detect and resubmit stuck jobs (checks node health via SSH, cancels first)")
    parser.add_argument("--prepared", "-p", action="store_true",
                        help="Submit all prepared jobs (from dry-run)")
    parser.add_argument("--jobs", "-j", nargs="+",
                        help="Specific job IDs to resubmit")
    parser.add_argument("--sample", "-s",
                        help="Only resubmit jobs for samples matching this regex")
    parser.add_argument("--memory", "-m", default=None,
                        help="Override memory request (e.g. 8gb, 16gb)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show what would be resubmitted without submitting")
    args = parser.parse_args()

    if not (args.failed or args.prepared or args.jobs or args.stuck):
        parser.error("Must specify --failed, --stuck, --prepared, or --jobs")

    # Determine paths
    script_dir = Path(__file__).parent.resolve()
    task_dir = script_dir / SLURM_OUTPUT_DIR / args.task

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
    print("\nQuerying SLURM for current status...")
    slurm_job_ids = [
        j["slurm_job_id"] for j in jobs.values()
        if j.get("slurm_job_id") is not None
    ]
    slurm_states, err = get_slurm_job_states(slurm_job_ids)
    if err:
        print(f"WARNING: {err}")

    # Determine statuses
    statuses = {}
    for job_id, job in jobs.items():
        statuses[job_id] = determine_job_status(job, slurm_states)

    # Detect stuck jobs via node health check
    stuck_ids = set()
    if args.stuck:
        detected, node_health = detect_stuck_jobs(jobs, statuses, slurm_states)
        stuck_ids = set(detected)

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
        for job_id, job in jobs.items():
            status = statuses[job_id]

            # Filter by sample if specified
            if args.sample:
                sample = job.get("sample", "")
                if not re.search(args.sample, sample):
                    continue

            if args.failed and status == "failed":
                jobs_to_resubmit.append((job_id, job))
            elif args.prepared and status == "prepared":
                jobs_to_resubmit.append((job_id, job))
            elif args.stuck and job_id in stuck_ids:
                jobs_to_resubmit.append((job_id, job))

    if not jobs_to_resubmit:
        print("\nNo jobs to resubmit.")
        return

    print(f"\nJobs to resubmit: {len(jobs_to_resubmit)}")
    if args.memory:
        print(f"Memory override: {args.memory}")

    if args.dry_run:
        print("\n[DRY RUN MODE]")

    # Cancel stuck jobs before resubmitting
    if args.stuck and not args.dry_run:
        slurm_ids_to_cancel = []
        for job_id, job in jobs_to_resubmit:
            sid = job.get("slurm_job_id")
            if sid is not None:
                slurm_ids_to_cancel.append(str(sid))
        if slurm_ids_to_cancel:
            print(f"\nCancelling {len(slurm_ids_to_cancel)} stuck jobs...")
            result = subprocess.run(
                ["scancel"] + slurm_ids_to_cancel,
                capture_output=True, text=True
            )
            if result.returncode != 0:
                print(f"  WARNING: scancel error: {result.stderr.strip()}")
            else:
                print(f"  Cancelled {len(slurm_ids_to_cancel)} jobs")

    # Resubmit jobs
    print("\nResubmitting...")
    success_count = 0
    fail_count = 0

    for job_id, job in jobs_to_resubmit:
        new_slurm_id = resubmit_job(job_id, job, memory=args.memory, dry_run=args.dry_run)

        if new_slurm_id is not None:
            success_count += 1
            if not args.dry_run and new_slurm_id > 0:
                # Update manifest
                manifest["jobs"][job_id]["slurm_job_id"] = new_slurm_id
                manifest["jobs"][job_id]["status"] = "submitted"
                manifest["jobs"][job_id]["submit_time"] = datetime.now().isoformat()
                manifest["jobs"][job_id]["resubmit_count"] = \
                    manifest["jobs"][job_id].get("resubmit_count", 0) + 1
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
