#!/usr/bin/env python3
"""
Check status of condor jobs for run3-vbsvvh preselection framework.

Reads the manifest.json and queries condor to determine job status.
Can also check for output files on XRootD.

Usage:
    # Check status of a task
    python condor/status.py --task 0Lep2FJ_run2-sig_0Lep2FJ_v1

    # List all tasks
    python condor/status.py --list

    # Show detailed status per sample
    python condor/status.py --task 0Lep2FJ_run2-sig_0Lep2FJ_v1 --by-sample

    # Check output files exist
    python condor/status.py --task 0Lep2FJ_run2-sig_0Lep2FJ_v1 --check-output

    # Show timing information
    python condor/status.py --task 0Lep2FJ_run2-sig_0Lep2FJ_v1 --timing

    # Show timing by sample
    python condor/status.py --task 0Lep2FJ_run2-sig_0Lep2FJ_v1 --timing --by-sample

    # Update manifest with status and timing
    python condor/status.py --task 0Lep2FJ_run2-sig_0Lep2FJ_v1 --timing --update
"""

import argparse
import json
import os
import re
import subprocess
import sys
from collections import defaultdict
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple


# Constants
CONDOR_OUTPUT_DIR = "jobs"


def get_condor_jobs() -> tuple[Dict[int, dict], Optional[str]]:
    """
    Query condor_q to get currently running/queued jobs.

    Returns tuple of (dict mapping cluster_id to job info, error message or None).
    """
    try:
        result = subprocess.run(
            ["condor_q", "-json"],
            capture_output=True, text=True, timeout=30
        )
        if result.returncode != 0:
            return {}, f"condor_q failed: {result.stderr.strip()}"

        jobs = json.loads(result.stdout) if result.stdout.strip() else []
        return {job["ClusterId"]: job for job in jobs}, None
    except subprocess.TimeoutExpired:
        return {}, "condor_q timed out after 30 seconds"
    except json.JSONDecodeError as e:
        return {}, f"Failed to parse condor_q output: {e}"
    except FileNotFoundError:
        return {}, "condor_q command not found"


def get_condor_history(cluster_ids: List[int]) -> tuple[Dict[int, dict], Optional[str]]:
    """
    Query condor_history for completed jobs.

    Returns tuple of (dict mapping cluster_id to job info, error message or None).
    """
    if not cluster_ids:
        return {}, None

    try:
        # Query history for specific cluster IDs
        constraint = " || ".join(f"ClusterId=={cid}" for cid in cluster_ids)
        result = subprocess.run(
            ["condor_history", "-json", "-constraint", constraint],
            capture_output=True, text=True, timeout=60
        )
        if result.returncode != 0:
            return {}, f"condor_history failed: {result.stderr.strip()}"

        jobs = json.loads(result.stdout) if result.stdout.strip() else []
        return {job["ClusterId"]: job for job in jobs}, None
    except subprocess.TimeoutExpired:
        return {}, "condor_history timed out after 60 seconds"
    except json.JSONDecodeError as e:
        return {}, f"Failed to parse condor_history output: {e}"
    except FileNotFoundError:
        return {}, "condor_history command not found"


def check_job_log(job_dir: str) -> Optional[dict]:
    """
    Parse condor job log to determine job status.

    Returns dict with status info or None if log not found.
    """
    job_path = Path(job_dir)

    # Find log file
    log_files = list(job_path.glob("job.*.log"))
    if not log_files:
        return None

    log_file = log_files[0]
    try:
        with open(log_file) as f:
            content = f.read()

        info = {"log_file": str(log_file)}

        # Check for job completion
        if "Job terminated" in content:
            info["terminated"] = True
            # Extract return value
            match = re.search(r"return value (\d+)", content)
            if match:
                info["exit_code"] = int(match.group(1))
            # Check for abnormal termination
            if "abnormal" in content.lower():
                info["abnormal"] = True
        else:
            info["terminated"] = False

        # Check for eviction
        if "Job was evicted" in content:
            info["evicted"] = True

        # Check if currently running
        if "Job executing" in content and not info.get("terminated"):
            info["running"] = True

        return info
    except IOError:
        return None


def parse_job_timing(job_dir: str) -> Optional[Dict[str, any]]:
    """
    Parse condor job log to extract timing information.

    Returns dict with start_time, end_time, and duration, or None if not available.
    """
    job_path = Path(job_dir)

    # Find log file (use most recent if multiple exist from resubmissions)
    log_files = sorted(job_path.glob("job.*.log"), key=lambda f: f.stat().st_mtime, reverse=True)
    if not log_files:
        return None

    log_file = log_files[0]
    try:
        with open(log_file) as f:
            content = f.read()

        timing = {}

        # Parse "Job executing" timestamp
        # Format: "001 (12345.000.000) 2024-01-15 10:30:45 Job executing on host..."
        exec_match = re.search(
            r'\d+ \(\d+\.\d+\.\d+\) (\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}) Job executing',
            content
        )
        if exec_match:
            timing["start_time"] = datetime.strptime(exec_match.group(1), "%Y-%m-%d %H:%M:%S")

        # Parse "Job terminated" timestamp
        # Format: "005 (12345.000.000) 2024-01-15 11:45:30 Job terminated."
        term_match = re.search(
            r'\d+ \(\d+\.\d+\.\d+\) (\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}) Job terminated',
            content
        )
        if term_match:
            timing["end_time"] = datetime.strptime(term_match.group(1), "%Y-%m-%d %H:%M:%S")

        # Calculate duration if both times available
        if "start_time" in timing and "end_time" in timing:
            timing["duration"] = timing["end_time"] - timing["start_time"]
        elif "start_time" in timing:
            # Job still running, calculate time since start
            timing["duration"] = datetime.now() - timing["start_time"]
            timing["running"] = True

        return timing if timing else None
    except IOError:
        return None


def format_duration(td: timedelta) -> str:
    """Format a timedelta as a human-readable string."""
    total_seconds = int(td.total_seconds())
    hours, remainder = divmod(total_seconds, 3600)
    minutes, seconds = divmod(remainder, 60)

    if hours > 0:
        return f"{hours}h {minutes}m {seconds}s"
    elif minutes > 0:
        return f"{minutes}m {seconds}s"
    else:
        return f"{seconds}s"


def determine_job_status(
    job: dict,
    condor_jobs: Dict[int, dict],
    condor_history: Dict[int, dict]
) -> str:
    """
    Determine the current status of a job.

    Returns one of: prepared, submitted, running, completed, failed, unknown
    """
    cluster_id = job.get("cluster_id")

    # Not yet submitted
    if cluster_id is None:
        return "prepared"

    # Check if in condor queue
    if cluster_id in condor_jobs:
        condor_status = condor_jobs[cluster_id].get("JobStatus", 0)
        # JobStatus: 1=Idle, 2=Running, 3=Removed, 4=Completed, 5=Held
        if condor_status == 1:
            return "queued"
        elif condor_status == 2:
            return "running"
        elif condor_status == 5:
            return "held"
        else:
            return "submitted"

    # Check history
    if cluster_id in condor_history:
        hist = condor_history[cluster_id]
        exit_code = hist.get("ExitCode", -1)
        if exit_code == 0:
            return "completed"
        else:
            return "failed"

    # Check job log as fallback
    job_dir = job.get("job_dir")
    if job_dir:
        log_info = check_job_log(job_dir)
        if log_info:
            if log_info.get("terminated"):
                if log_info.get("exit_code", -1) == 0:
                    return "completed"
                else:
                    return "failed"
            elif log_info.get("running"):
                return "running"

    # If submitted but can't determine status
    return "unknown"


def load_manifest(task_dir: Path) -> Optional[dict]:
    """Load manifest.json from task directory."""
    manifest_path = task_dir / "manifest.json"
    if not manifest_path.exists():
        return None
    with open(manifest_path) as f:
        return json.load(f)


def save_manifest(task_dir: Path, manifest: dict):
    """Save updated manifest to disk."""
    manifest_path = task_dir / "manifest.json"
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)


def list_tasks(condor_dir: Path) -> List[str]:
    """List all task directories."""
    jobs_dir = condor_dir / CONDOR_OUTPUT_DIR
    if not jobs_dir.exists():
        return []
    return sorted([d.name for d in jobs_dir.iterdir() if d.is_dir()])


def print_status_table(status_counts: Dict[str, int], total: int):
    """Print a formatted status table."""
    # Define status order and colors (ANSI)
    status_order = ["completed", "running", "queued", "submitted", "held", "failed", "prepared", "unknown"]
    colors = {
        "completed": "\033[92m",  # Green
        "running": "\033[94m",    # Blue
        "queued": "\033[93m",     # Yellow
        "submitted": "\033[93m",  # Yellow
        "held": "\033[91m",       # Red
        "failed": "\033[91m",     # Red
        "prepared": "\033[90m",   # Gray
        "unknown": "\033[90m",    # Gray
    }
    reset = "\033[0m"

    print(f"\n{'Status':<12} {'Count':>8} {'Percent':>10}")
    print("-" * 32)

    for status in status_order:
        count = status_counts.get(status, 0)
        if count > 0:
            pct = 100.0 * count / total if total > 0 else 0
            color = colors.get(status, "")
            print(f"{color}{status:<12}{reset} {count:>8} {pct:>9.1f}%")

    print("-" * 32)
    print(f"{'Total':<12} {total:>8}")


def print_sample_summary(jobs: Dict[str, dict], statuses: Dict[str, str]):
    """Print status summary grouped by sample."""
    # Group by sample
    by_sample = defaultdict(lambda: {"statuses": defaultdict(int), "resubmits": 0})
    for job_id, job in jobs.items():
        sample = job.get("sample", "unknown")
        status = statuses.get(job_id, "unknown")
        by_sample[sample]["statuses"][status] += 1
        by_sample[sample]["resubmits"] += job.get("resubmit_count", 0)

    # Print table
    print(f"\n{'Sample':<50} {'Done':>6} {'Run':>6} {'Queue':>6} {'Fail':>6} {'Resub':>6} {'Total':>6}")
    print("-" * 92)

    for sample in sorted(by_sample.keys()):
        counts = by_sample[sample]["statuses"]
        resubmits = by_sample[sample]["resubmits"]
        done = counts.get("completed", 0)
        run = counts.get("running", 0)
        queue = counts.get("queued", 0) + counts.get("submitted", 0)
        fail = counts.get("failed", 0) + counts.get("held", 0)
        total = sum(counts.values())

        # Truncate long sample names
        sample_display = sample[:47] + "..." if len(sample) > 50 else sample
        print(f"{sample_display:<50} {done:>6} {run:>6} {queue:>6} {fail:>6} {resubmits:>6} {total:>6}")


def print_failed_jobs(jobs: Dict[str, dict], statuses: Dict[str, str]):
    """Print details of failed jobs."""
    failed = [(jid, j) for jid, j in jobs.items() if statuses.get(jid) in ("failed", "held")]

    if not failed:
        print("\nNo failed jobs.")
        return

    print(f"\n{'='*60}")
    print(f"Failed Jobs ({len(failed)})")
    print(f"{'='*60}")

    for job_id, job in failed[:20]:  # Limit to first 20
        resubmit_count = job.get('resubmit_count', 0)
        print(f"\n  {job_id}")
        print(f"    Sample: {job.get('sample')}")
        print(f"    Files:  {job.get('n_files', 0)}")
        print(f"    Resubmits: {resubmit_count}")
        print(f"    Dir:    {job.get('job_dir')}")

        # Check for error in stderr (use most recent file)
        job_dir = Path(job.get('job_dir', ''))
        stderr_files = sorted(job_dir.glob("job.*.stderr"), key=lambda f: f.stat().st_mtime, reverse=True)
        if stderr_files:
            try:
                with open(stderr_files[0]) as f:
                    stderr = f.read().strip()
                if stderr:
                    # Show last few lines
                    lines = stderr.split('\n')[-5:]
                    print(f"    Error (last 5 lines):")
                    for line in lines:
                        print(f"      {line}")
            except IOError:
                pass

    if len(failed) > 20:
        print(f"\n  ... and {len(failed) - 20} more failed jobs")


def print_timing_summary(jobs: Dict[str, dict], statuses: Dict[str, str], timings: Dict[str, dict]):
    """Print timing summary for completed and running jobs."""
    # Collect timing data
    completed_times = []
    running_times = []

    for job_id, timing in timings.items():
        if timing and "duration" in timing:
            if timing.get("running"):
                running_times.append((job_id, timing))
            else:
                completed_times.append((job_id, timing))

    if not completed_times and not running_times:
        print("\nNo timing information available.")
        return

    print(f"\n{'='*60}")
    print("Timing Summary")
    print(f"{'='*60}")

    if completed_times:
        durations = [t["duration"] for _, t in completed_times]
        total_duration = sum(durations, timedelta())
        avg_duration = total_duration / len(durations)
        min_duration = min(durations)
        max_duration = max(durations)

        print(f"\nCompleted jobs: {len(completed_times)}")
        print(f"  Average runtime: {format_duration(avg_duration)}")
        print(f"  Min runtime:     {format_duration(min_duration)}")
        print(f"  Max runtime:     {format_duration(max_duration)}")
        print(f"  Total CPU time:  {format_duration(total_duration)}")

    if running_times:
        print(f"\nCurrently running: {len(running_times)}")
        for job_id, timing in sorted(running_times, key=lambda x: x[1]["duration"], reverse=True)[:5]:
            sample = jobs[job_id].get("sample", "unknown")
            print(f"  {sample}: running for {format_duration(timing['duration'])}")
        if len(running_times) > 5:
            print(f"  ... and {len(running_times) - 5} more")


def print_timing_by_sample(jobs: Dict[str, dict], timings: Dict[str, dict]):
    """Print timing summary grouped by sample."""
    # Group by sample
    by_sample = defaultdict(list)
    for job_id, timing in timings.items():
        if timing and "duration" in timing and not timing.get("running"):
            sample = jobs[job_id].get("sample", "unknown")
            by_sample[sample].append(timing["duration"])

    if not by_sample:
        print("\nNo completed job timing available.")
        return

    print(f"\n{'Sample':<50} {'Jobs':>6} {'Avg':>12} {'Total':>12}")
    print("-" * 82)

    for sample in sorted(by_sample.keys()):
        durations = by_sample[sample]
        n_jobs = len(durations)
        total = sum(durations, timedelta())
        avg = total / n_jobs

        sample_display = sample[:47] + "..." if len(sample) > 50 else sample
        print(f"{sample_display:<50} {n_jobs:>6} {format_duration(avg):>12} {format_duration(total):>12}")


def main():
    parser = argparse.ArgumentParser(
        description="Check status of condor jobs",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--task", "-t",
                        help="Task name (directory name under condor/jobs/)")
    parser.add_argument("--list", "-l", action="store_true",
                        help="List all tasks")
    parser.add_argument("--by-sample", "-s", action="store_true",
                        help="Show status grouped by sample")
    parser.add_argument("--failed", "-f", action="store_true",
                        help="Show details of failed jobs")
    parser.add_argument("--update", "-u", action="store_true",
                        help="Update manifest with current status")
    parser.add_argument("--timing", action="store_true",
                        help="Show timing information from job logs")
    parser.add_argument("--no-condor", action="store_true",
                        help="Don't query condor (use manifest/logs only)")
    args = parser.parse_args()

    # Determine paths
    script_dir = Path(__file__).parent.resolve()
    condor_dir = script_dir
    preselection_dir = script_dir.parent

    # List tasks
    if args.list:
        tasks = list_tasks(condor_dir)
        if not tasks:
            print("No tasks found.")
            return
        print("Available tasks:")
        for task in tasks:
            manifest = load_manifest(condor_dir / CONDOR_OUTPUT_DIR / task)
            if manifest:
                n_jobs = len(manifest.get("jobs", {}))
                created = manifest.get("created", "")[:10]
                print(f"  {task} ({n_jobs} jobs, created {created})")
            else:
                print(f"  {task} (no manifest)")
        return

    # Require task for other operations
    if not args.task:
        parser.error("--task is required (use --list to see available tasks)")

    task_dir = condor_dir / CONDOR_OUTPUT_DIR / args.task
    if not task_dir.exists():
        print(f"ERROR: Task not found: {task_dir}")
        print("Use --list to see available tasks")
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
    print(f"Config: {manifest.get('config', 'N/A')}")
    print(f"Created: {manifest.get('created', 'N/A')}")
    print(f"Total jobs: {len(jobs)}")

    # Get condor status
    condor_jobs = {}
    condor_history = {}

    if not args.no_condor:
        print("\nQuerying condor...")
        condor_jobs, err = get_condor_jobs()
        if err:
            print(f"WARNING: {err}")

        # Get cluster IDs from manifest
        cluster_ids = [
            j["cluster_id"] for j in jobs.values()
            if j.get("cluster_id") is not None
        ]

        # Only query history for jobs not in queue
        missing_ids = [cid for cid in cluster_ids if cid not in condor_jobs]
        if missing_ids:
            condor_history, err = get_condor_history(missing_ids)
            if err:
                print(f"WARNING: {err}")

    # Determine status for each job
    statuses = {}
    for job_id, job in jobs.items():
        statuses[job_id] = determine_job_status(job, condor_jobs, condor_history)

    # Count by status and resubmissions
    status_counts = defaultdict(int)
    resubmit_counts = defaultdict(int)  # resubmit_count -> number of jobs
    total_resubmissions = 0
    for job_id, status in statuses.items():
        status_counts[status] += 1
        resubmit_count = jobs[job_id].get("resubmit_count", 0)
        if resubmit_count > 0:
            resubmit_counts[resubmit_count] += 1
            total_resubmissions += resubmit_count

    # Print summary
    print_status_table(dict(status_counts), len(jobs))

    # Print resubmission info if any
    if total_resubmissions > 0:
        print(f"\nResubmissions: {total_resubmissions} total")
        for count in sorted(resubmit_counts.keys()):
            print(f"  {resubmit_counts[count]} job(s) resubmitted {count} time(s)")

    # Print by-sample summary
    if args.by_sample:
        print_sample_summary(jobs, statuses)

    # Print failed jobs
    if args.failed:
        print_failed_jobs(jobs, statuses)

    # Collect and print timing information
    timings = {}
    if args.timing:
        print("\nParsing job logs for timing...")
        for job_id, job in jobs.items():
            job_dir = job.get("job_dir")
            if job_dir:
                timings[job_id] = parse_job_timing(job_dir)

        print_timing_summary(jobs, statuses, timings)
        if args.by_sample:
            print_timing_by_sample(jobs, timings)

    # Update manifest
    if args.update:
        for job_id, status in statuses.items():
            if job_id in manifest["jobs"]:
                manifest["jobs"][job_id]["status"] = status

        # Store timing info if --timing was also specified
        if args.timing and timings:
            for job_id, timing in timings.items():
                if timing and job_id in manifest["jobs"]:
                    timing_info = {}
                    if "start_time" in timing:
                        timing_info["start_time"] = timing["start_time"].isoformat()
                    if "end_time" in timing:
                        timing_info["end_time"] = timing["end_time"].isoformat()
                    if "duration" in timing:
                        timing_info["duration_seconds"] = int(timing["duration"].total_seconds())
                    if timing_info:
                        manifest["jobs"][job_id]["timing"] = timing_info

        manifest["last_checked"] = datetime.now().isoformat()
        save_manifest(task_dir, manifest)
        print(f"\nManifest updated: {task_dir}/manifest.json")


if __name__ == "__main__":
    main()
