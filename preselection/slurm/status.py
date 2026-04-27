#!/usr/bin/env python3
"""
Check status of SLURM jobs for run3-vbsvvh preselection framework.

Reads the manifest.json and queries SLURM to determine job status.

Usage:
    # Check status of a task
    python slurm/status.py --task merged_r2_0lep_0FJ_..._0lep_0FJ

    # List all tasks
    python slurm/status.py --list

    # Show detailed status per sample
    python slurm/status.py --task <name> --by-sample

    # Show failed job details
    python slurm/status.py --task <name> --failed

    # Show timing information
    python slurm/status.py --task <name> --timing

    # Update manifest with current status
    python slurm/status.py --task <name> --update
"""

import argparse
import json
import os
import re
import subprocess
import sys
import time
from collections import defaultdict
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Optional


# Constants
SLURM_OUTPUT_DIR = "jobs"


# ============================================================================
# SLURM Query Functions
# ============================================================================

def get_slurm_job_states(slurm_job_ids: list) -> tuple[Dict, Optional[str]]:
    """
    Query sacct for job states.

    Accepts both integer (individual) and string (array, e.g. "12345_0") job IDs.
    Returns tuple of (dict mapping slurm_job_id to job info, error message or None).
    Keys match the manifest format: int for individual jobs, str for array tasks.
    """
    if not slurm_job_ids:
        return {}, None

    try:
        # Deduplicate base IDs for efficient sacct query
        # "12345_0" -> query "12345" (returns all array tasks)
        query_ids = set()
        for jid in slurm_job_ids:
            jid_str = str(jid)
            if "_" in jid_str:
                query_ids.add(jid_str.split("_")[0])
            else:
                query_ids.add(jid_str)

        job_id_str = ",".join(query_ids)
        result = subprocess.run(
            ["sacct", "-j", job_id_str, "--noheader", "--parsable2",
             "--format=JobID,State,ExitCode,Start,End,Elapsed,NodeList"],
            capture_output=True, text=True, timeout=30
        )
        if result.returncode != 0:
            return {}, f"sacct failed: {result.stderr.strip()}"

        jobs = {}
        for line in result.stdout.strip().split("\n"):
            if not line:
                continue
            fields = line.split("|")
            if len(fields) < 7:
                continue

            job_id_field = fields[0]
            # Skip sub-steps (e.g. 12345.batch, 12345_0.batch)
            if "." in job_id_field:
                continue
            # Expand array summary lines (e.g. 12345_[0-99]) into
            # individual entries so pending tasks are not missed.
            if "[" in job_id_field:
                m = re.match(r"(\d+)_\[(.+)\]", job_id_field)
                if m:
                    base_id = m.group(1)
                    state_info = {
                        "state": fields[1],
                        "exit_code": fields[2],
                        "start": fields[3],
                        "end": fields[4],
                        "elapsed": fields[5],
                        "nodelist": fields[6],
                    }
                    for part in m.group(2).split(","):
                        if "-" in part:
                            lo, hi = part.split("-", 1)
                            for idx in range(int(lo), int(hi) + 1):
                                jobs[f"{base_id}_{idx}"] = state_info
                        else:
                            jobs[f"{base_id}_{part}"] = state_info
                continue

            # Use matching key type: str for array tasks, int for individual
            if "_" in job_id_field:
                key = job_id_field  # e.g. "12345_0"
            else:
                try:
                    key = int(job_id_field)
                except ValueError:
                    continue

            jobs[key] = {
                "state": fields[1],
                "exit_code": fields[2],
                "start": fields[3],
                "end": fields[4],
                "elapsed": fields[5],
                "nodelist": fields[6],
            }

        # sacct doesn't list individual PENDING array tasks, only summary
        # ranges like "12345_[48-213]". Query squeue to fill in pending jobs.
        try:
            sq_result = subprocess.run(
                ["squeue", "-j", job_id_str, "--noheader",
                 "--format=%i|%T|%N"],
                capture_output=True, text=True, timeout=30
            )
            if sq_result.returncode == 0:
                for line in sq_result.stdout.strip().split("\n"):
                    if not line:
                        continue
                    parts = line.strip().split("|")
                    if len(parts) < 3:
                        continue
                    sq_job_id = parts[0].strip()
                    sq_state = parts[1].strip()
                    sq_node = parts[2].strip()
                    if "_" in sq_job_id:
                        key = sq_job_id
                    else:
                        try:
                            key = int(sq_job_id)
                        except ValueError:
                            continue
                    if key not in jobs:
                        jobs[key] = {
                            "state": sq_state,
                            "exit_code": "",
                            "start": "",
                            "end": "",
                            "elapsed": "",
                            "nodelist": sq_node,
                        }
        except (subprocess.TimeoutExpired, FileNotFoundError):
            pass  # squeue not available, proceed with sacct data only

        return jobs, None
    except subprocess.TimeoutExpired:
        return {}, "sacct timed out after 30 seconds"
    except FileNotFoundError:
        return {}, "sacct command not found"


def determine_job_status(job: dict, slurm_states: Dict[int, dict]) -> str:
    """
    Determine the current status of a job.

    Returns one of: prepared, submitted, running, completed, failed, unknown
    """
    slurm_job_id = job.get("slurm_job_id")

    # Not yet submitted
    if slurm_job_id is None:
        return "prepared"

    # Check sacct state
    if slurm_job_id in slurm_states:
        state = slurm_states[slurm_job_id]["state"]
        if state == "COMPLETED":
            return "completed"
        elif state in ("RUNNING", "COMPLETING"):
            return "running"
        elif state == "PENDING":
            return "queued"
        elif state.startswith("CANCELLED") or state in (
                "FAILED", "TIMEOUT", "OUT_OF_MEMORY",
                "NODE_FAIL", "PREEMPTED"):
            return "failed"
        else:
            return "submitted"

    # Check SLURM log files as fallback
    job_dir = job.get("job_dir")
    if job_dir:
        job_path = Path(job_dir)
        out_files = sorted(job_path.glob("slurm_*.out"),
                           key=lambda f: f.stat().st_mtime, reverse=True)
        if out_files:
            try:
                with open(out_files[0]) as f:
                    content = f.read()
                if "Job completed at" in content:
                    return "completed"
                elif "ERROR:" in content:
                    return "failed"
            except IOError:
                pass

    return "unknown"


def parse_elapsed(elapsed_str: str) -> Optional[timedelta]:
    """Parse SLURM elapsed time format (D-HH:MM:SS or HH:MM:SS) to timedelta."""
    try:
        if "-" in elapsed_str:
            days, hms = elapsed_str.split("-")
            days = int(days)
        else:
            days = 0
            hms = elapsed_str

        parts = hms.split(":")
        if len(parts) == 3:
            hours, minutes, seconds = int(parts[0]), int(parts[1]), int(parts[2])
        elif len(parts) == 2:
            hours, minutes, seconds = 0, int(parts[0]), int(parts[1])
        else:
            return None

        return timedelta(days=days, hours=hours, minutes=minutes, seconds=seconds)
    except (ValueError, IndexError):
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


# ============================================================================
# Display Functions
# ============================================================================

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


def list_tasks(slurm_dir: Path) -> List[str]:
    """List all task directories."""
    jobs_dir = slurm_dir / SLURM_OUTPUT_DIR
    if not jobs_dir.exists():
        return []
    return sorted([d.name for d in jobs_dir.iterdir() if d.is_dir()])


def print_status_table(status_counts: Dict[str, int], total: int):
    """Print a formatted status table."""
    status_order = ["completed", "running", "stuck", "queued", "submitted", "failed", "prepared", "unknown"]
    colors = {
        "completed": "\033[92m",  # Green
        "running": "\033[94m",    # Blue
        "stuck": "\033[95m",      # Magenta
        "queued": "\033[93m",     # Yellow
        "submitted": "\033[93m",  # Yellow
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
    by_sample = defaultdict(lambda: defaultdict(int))
    for job_id, job in jobs.items():
        sample = job.get("sample", "unknown")
        status = statuses.get(job_id, "unknown")
        by_sample[sample][status] += 1

    print(f"\n{'Sample':<50} {'Done':>6} {'Run':>6} {'Queue':>6} {'Fail':>6} {'Total':>6}")
    print("-" * 80)

    for sample in sorted(by_sample.keys()):
        counts = by_sample[sample]
        done = counts.get("completed", 0)
        run = counts.get("running", 0)
        queue = counts.get("queued", 0) + counts.get("submitted", 0)
        fail = counts.get("failed", 0)
        total = sum(counts.values())

        sample_display = sample[:47] + "..." if len(sample) > 50 else sample
        print(f"{sample_display:<50} {done:>6} {run:>6} {queue:>6} {fail:>6} {total:>6}")


def check_node_health(node: str, user: str = None,
                      test_path: str = "/blue",
                      timeout_sec: int = 10) -> dict:
    """
    SSH to a compute node and check if jobs are stuck.

    Performs two checks:
    1. Shared filesystem responsiveness (stat /blue)
    2. Process state: if ALL user processes are in D-state (uninterruptible
       sleep) with 0 CPU time, the node is likely hung on I/O.

    Returns dict with keys: healthy (bool), error (str or None).
    """
    if user is None:
        user = os.environ.get("USER", "")

    try:
        # Check filesystem AND process states in one SSH call
        cmd = (
            f"timeout {timeout_sec} stat {test_path} > /dev/null 2>&1"
            f" && echo FS_OK || echo FS_BAD;"
            f" ps -u {user} -o pid=,stat=,time=,comm= 2>/dev/null"
        )
        result = subprocess.run(
            ["ssh", "-o", "ConnectTimeout=5", "-o", "BatchMode=yes",
             node, cmd],
            capture_output=True, text=True, timeout=timeout_sec + 10
        )
        output = result.stdout.strip()
        lines = output.split("\n")

        # Check filesystem
        fs_ok = any("FS_OK" in line for line in lines)
        if not fs_ok:
            return {"healthy": False,
                    "error": f"filesystem unresponsive on {node}"}

        # Check process states: are ALL job processes in D-state?
        # Filter out SSH session processes (sshd, bash, ps) which are
        # spawned by this health check itself.
        ssh_procs = {"sshd", "bash", "ps"}
        proc_lines = [l.strip() for l in lines
                      if l.strip() and "FS_" not in l]
        job_procs = []
        for line in proc_lines:
            fields = line.split()
            if len(fields) >= 4 and fields[3] not in ssh_procs:
                job_procs.append(fields)
        if job_procs:
            all_d_state = all(f[1].startswith("D") for f in job_procs)
            if all_d_state:
                return {"healthy": False,
                        "error": f"all processes in D-state on {node}"}

        return {"healthy": True, "error": None}
    except subprocess.TimeoutExpired:
        return {"healthy": False, "error": f"SSH to {node} timed out"}
    except Exception as e:
        return {"healthy": False, "error": f"SSH to {node} failed: {e}"}


def detect_stuck_jobs(jobs: Dict[str, dict], statuses: Dict[str, str],
                      slurm_states: Dict[int, dict]) -> tuple[List[str], Dict[str, dict]]:
    """
    Detect stuck jobs by checking compute node filesystem health.

    For each node with running jobs, SSHs in and verifies the shared
    filesystem is responsive. Jobs on unresponsive nodes are flagged.

    Returns (list of stuck job_ids, dict of node -> health_info).
    """
    # Group running jobs by node
    node_jobs = defaultdict(list)
    for job_id, job in jobs.items():
        if statuses.get(job_id) != "running":
            continue
        slurm_job_id = job.get("slurm_job_id")
        if slurm_job_id not in slurm_states:
            continue
        node = slurm_states[slurm_job_id].get("nodelist", "")
        if node:
            node_jobs[node].append(job_id)

    if not node_jobs:
        return [], {}

    print(f"\nChecking {len(node_jobs)} node(s) with running jobs...")

    node_health = {}
    stuck = []
    for node, job_ids in node_jobs.items():
        health = check_node_health(node)
        node_health[node] = health
        if health["healthy"]:
            print(f"  {node}: OK ({len(job_ids)} jobs)")
        else:
            print(f"  {node}: \033[91mUNHEALTHY\033[0m - {health['error']}"
                  f" ({len(job_ids)} jobs)")
            stuck.extend(job_ids)

    return stuck, node_health


def print_stuck_jobs(stuck_ids: List[str], jobs: Dict[str, dict],
                     slurm_states: Dict[int, dict]):
    """Print details of jobs detected as stuck on unhealthy nodes."""
    if not stuck_ids:
        return

    yellow = "\033[93m"
    reset = "\033[0m"

    print(f"\n{yellow}{'='*60}")
    print(f"WARNING: Stuck Jobs on Unhealthy Nodes ({len(stuck_ids)})")
    print(f"{'='*60}{reset}")
    print("The shared filesystem is unresponsive on these nodes.")
    print("These jobs will never complete. Cancel and resubmit.\n")

    for job_id in stuck_ids:
        job = jobs[job_id]
        slurm_job_id = job.get("slurm_job_id")
        slurm_info = slurm_states.get(slurm_job_id, {})
        print(f"  {yellow}{job_id}{reset}")
        print(f"    SLURM ID: {slurm_job_id}  Node: {slurm_info.get('nodelist', 'N/A')}"
              f"  Elapsed: {slurm_info.get('elapsed', 'N/A')}")

    ids = " ".join(str(jobs[jid].get("slurm_job_id")) for jid in stuck_ids)
    print(f"\n  To cancel:  scancel {ids}")


def print_failed_jobs(jobs: Dict[str, dict], statuses: Dict[str, str],
                      slurm_states: Dict[int, dict]):
    """Print details of failed jobs."""
    failed = [(jid, j) for jid, j in jobs.items() if statuses.get(jid) == "failed"]

    if not failed:
        print("\nNo failed jobs.")
        return

    red = "\033[91m"
    reset = "\033[0m"

    print(f"\n{'='*60}")
    print(f"Failed Jobs ({len(failed)})")
    print(f"{'='*60}")

    for job_id, job in failed[:20]:
        slurm_job_id = job.get("slurm_job_id")
        slurm_info = slurm_states.get(slurm_job_id, {})

        print(f"\n  {red}{job_id}{reset}")
        print(f"    Sample:    {job.get('sample')}")
        print(f"    Files:     {job.get('n_files', 0)}")
        print(f"    SLURM ID:  {slurm_job_id}")
        print(f"    State:     {slurm_info.get('state', 'N/A')}")
        print(f"    Exit code: {slurm_info.get('exit_code', 'N/A')}")
        print(f"    Node:      {slurm_info.get('nodelist', 'N/A')}")

        # Show tail of stdout and stderr
        job_dir = Path(job.get('job_dir', ''))
        for pattern, label in [("slurm_*.out", "STDOUT"), ("slurm_*.err", "STDERR")]:
            log_files = sorted(job_dir.glob(pattern),
                               key=lambda f: f.stat().st_mtime, reverse=True)
            if log_files:
                try:
                    with open(log_files[0]) as f:
                        content = f.read().strip()
                    if content:
                        lines = content.split('\n')[-5:]
                        print(f"\n    {label} (last {len(lines)} lines):")
                        for line in lines:
                            print(f"      {line}")
                except IOError:
                    pass

    if len(failed) > 20:
        print(f"\n  ... and {len(failed) - 20} more failed jobs")


def print_timing_summary(jobs: Dict[str, dict], slurm_states: Dict[int, dict]):
    """Print timing summary from sacct."""
    completed_times = []
    running_times = []

    for job_id, job in jobs.items():
        slurm_job_id = job.get("slurm_job_id")
        if slurm_job_id not in slurm_states:
            continue
        info = slurm_states[slurm_job_id]
        elapsed = parse_elapsed(info.get("elapsed", ""))
        if elapsed is None:
            continue
        if info["state"] == "COMPLETED":
            completed_times.append((job_id, elapsed))
        elif info["state"] == "RUNNING":
            running_times.append((job_id, elapsed))

    if not completed_times and not running_times:
        print("\nNo timing information available.")
        return

    print(f"\n{'='*60}")
    print("Timing Summary")
    print(f"{'='*60}")

    if completed_times:
        durations = [d for _, d in completed_times]
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
        for job_id, elapsed in sorted(running_times, key=lambda x: x[1], reverse=True)[:5]:
            sample = jobs[job_id].get("sample", "unknown")
            print(f"  {sample}: running for {format_duration(elapsed)}")
        if len(running_times) > 5:
            print(f"  ... and {len(running_times) - 5} more")


def print_timing_by_sample(jobs: Dict[str, dict], slurm_states: Dict[int, dict]):
    """Print timing summary grouped by sample."""
    by_sample = defaultdict(list)
    for job_id, job in jobs.items():
        slurm_job_id = job.get("slurm_job_id")
        if slurm_job_id not in slurm_states:
            continue
        info = slurm_states[slurm_job_id]
        if info["state"] != "COMPLETED":
            continue
        elapsed = parse_elapsed(info.get("elapsed", ""))
        if elapsed is None:
            continue
        sample = job.get("sample", "unknown")
        by_sample[sample].append(elapsed)

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


# ============================================================================
# Main
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Check status of SLURM jobs",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--task", "-t",
                        help="Task name (directory name under slurm/jobs/)")
    parser.add_argument("--list", "-l", action="store_true",
                        help="List all tasks")
    parser.add_argument("--by-sample", "-s", action="store_true",
                        help="Show status grouped by sample")
    parser.add_argument("--failed", "-f", action="store_true",
                        help="Show details of failed jobs")
    parser.add_argument("--update", "-u", action="store_true",
                        help="Update manifest with current status")
    parser.add_argument("--timing", action="store_true",
                        help="Show timing information")
    parser.add_argument("--check-nodes", action="store_true",
                        help="SSH to compute nodes with running jobs and check filesystem health")
    parser.add_argument("--monitor", type=int, nargs="?", const=10, default=None,
                        metavar="MINUTES",
                        help="Monitor mode: check for stuck jobs every N minutes "
                             "(default: 10) and auto-resubmit them. Runs until all "
                             "jobs are completed or failed.")
    args = parser.parse_args()

    # --monitor implies --check-nodes
    if args.monitor is not None:
        args.check_nodes = True

    # Determine paths
    script_dir = Path(__file__).parent.resolve()
    slurm_dir = script_dir

    # List tasks
    if args.list:
        tasks = list_tasks(slurm_dir)
        if not tasks:
            print("No tasks found.")
            return
        print("Available tasks:")
        for task in tasks:
            manifest = load_manifest(slurm_dir / SLURM_OUTPUT_DIR / task)
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

    task_dir = slurm_dir / SLURM_OUTPUT_DIR / args.task
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

    # Query SLURM for job states
    print("\nQuerying SLURM...")
    slurm_job_ids = [
        j["slurm_job_id"] for j in jobs.values()
        if j.get("slurm_job_id") is not None
    ]
    slurm_states, err = get_slurm_job_states(slurm_job_ids)
    if err:
        print(f"WARNING: {err}")

    # Determine status for each job
    statuses = {}
    for job_id, job in jobs.items():
        statuses[job_id] = determine_job_status(job, slurm_states)

    # Count by status
    status_counts = defaultdict(int)
    for status in statuses.values():
        status_counts[status] += 1

    # Print summary
    print_status_table(dict(status_counts), len(jobs))

    # Print by-sample summary
    if args.by_sample:
        print_sample_summary(jobs, statuses)

    # Print failed jobs
    if args.failed:
        print_failed_jobs(jobs, statuses, slurm_states)

    # Check node health for running jobs
    if args.check_nodes:
        stuck_ids, node_health = detect_stuck_jobs(jobs, statuses, slurm_states)
        if stuck_ids:
            # Reclassify stuck jobs
            for job_id in stuck_ids:
                statuses[job_id] = "stuck"
            # Recount
            status_counts = defaultdict(int)
            for status in statuses.values():
                status_counts[status] += 1
            print_status_table(dict(status_counts), len(jobs))
            print_stuck_jobs(stuck_ids, jobs, slurm_states)

    # Monitor mode: auto-resubmit stuck jobs in a loop
    if args.monitor is not None and stuck_ids:
        print(f"\n[Monitor] Resubmitting {len(stuck_ids)} stuck job(s)...")
        resubmit_cmd = [
            sys.executable, str(script_dir / "resubmit.py"),
            "-t", args.task, "--stuck"
        ]
        subprocess.run(resubmit_cmd)

    # Print timing
    if args.timing:
        print_timing_summary(jobs, slurm_states)
        if args.by_sample:
            print_timing_by_sample(jobs, slurm_states)

    # Update manifest
    if args.update:
        for job_id, status in statuses.items():
            if job_id in manifest["jobs"]:
                manifest["jobs"][job_id]["status"] = status

        # Store timing info
        if args.timing:
            for job_id, job in jobs.items():
                slurm_job_id = job.get("slurm_job_id")
                if slurm_job_id in slurm_states:
                    info = slurm_states[slurm_job_id]
                    elapsed = parse_elapsed(info.get("elapsed", ""))
                    if elapsed is not None:
                        manifest["jobs"][job_id]["timing"] = {
                            "start": info.get("start"),
                            "end": info.get("end"),
                            "elapsed": info.get("elapsed"),
                            "duration_seconds": int(elapsed.total_seconds()),
                        }

        manifest["last_checked"] = datetime.now().isoformat()
        save_manifest(task_dir, manifest)
        print(f"\nManifest updated: {task_dir}/manifest.json")

    # Monitor loop: repeat every N minutes until all jobs done
    if args.monitor is not None:
        n_running = status_counts.get("running", 0) + status_counts.get("queued", 0)
        if n_running == 0:
            print(f"\n[Monitor] No running/queued jobs remaining. Done.")
            return

        interval = args.monitor
        print(f"\n[Monitor] {n_running} job(s) still active. "
              f"Next check in {interval} minutes... (Ctrl+C to stop)")
        try:
            time.sleep(interval * 60)
        except KeyboardInterrupt:
            print("\n[Monitor] Stopped.")
            return

        # Re-run main with same args (recursive call via exec to refresh manifest)
        os.execv(sys.executable, [sys.executable] + sys.argv)


if __name__ == "__main__":
    main()
