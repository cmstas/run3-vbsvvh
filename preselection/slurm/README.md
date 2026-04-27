# SLURM Job Submission Framework

Scripts for submitting, monitoring, and resubmitting preselection jobs on SLURM.

All commands should be run from the `preselection/` directory.

## Scripts

| Script | Purpose |
|--------|---------|
| `submit.py` | Split samples into jobs, create sbatch scripts, and submit |
| `status.py` | Check job status, show failures, timing, and node health |
| `resubmit.py` | Resubmit failed or stuck jobs |
| `executable.sh` | The actual job payload (setup CMSSW, run analysis, copy output) |

## Quick Start

### 1. Submit jobs

```bash
python3 slurm/submit.py \
  -c merged_jsons/merged_r2_0lep_0FJ.json \
  -a 0lep_0FJ -r 2 \
  -o /blue/avery/$USER/samples/run3-vbsvvh/r2_0lep_0FJ \
  --files-per-job 1
```

This creates a task directory under `slurm/jobs/<task_name>/` containing:
- `manifest.json` -- tracks all jobs, their inputs, SLURM IDs, and status
- `array.sbatch` -- single SLURM array job covering all sub-jobs
- `job_list.txt` -- maps array task indices to sample/job parameters
- `executable.sh`, `package.tar.gz` -- analysis code shipped to compute nodes
- `<sample_name>/` -- per-job directories with `config.json` and `job.sbatch`

#### Key options

| Flag | Description |
|------|-------------|
| `-c, --config` | Path to merged config JSON (required) |
| `-a, --analysis` | Analysis channel, e.g. `0lep_0FJ` (required) |
| `-r, --run_number` | Run 2 or Run 3 (required) |
| `-o, --output-dir` | Final output directory (required) |
| `--files-per-job N` | Input files per job (default: 10) |
| `--events-per-job N` | Max events per job (overrides `--files-per-job`) |
| `-j, --ncpus N` | CPUs per job (default: 4) |
| `-m, --memory` | Memory per job (default: 8gb) |
| `--sample REGEX` | Only submit samples matching this pattern |
| `--spanet-infer` | Run SPANet inference |
| `--spanet-training` | Generate SPANet training data |
| `-t, --tag` | Optional tag appended to task name |
| `--dry-run` | Prepare everything but don't submit |

### 2. Check status

```bash
# Basic status
python3 slurm/status.py -t <task_name>

# List all tasks
python3 slurm/status.py --list

# Show status grouped by sample
python3 slurm/status.py -t <task_name> --by-sample

# Show failed job details (stderr/stdout tails)
python3 slurm/status.py -t <task_name> --failed

# Show timing statistics
python3 slurm/status.py -t <task_name> --timing

# Check compute node health (SSH to nodes with running jobs)
python3 slurm/status.py -t <task_name> --check-nodes

# Monitor mode: check every 10 min, auto-resubmit stuck jobs
python3 slurm/status.py -t <task_name> --monitor

# Monitor with custom interval (e.g. every 5 minutes)
python3 slurm/status.py -t <task_name> --monitor 5

# Update manifest with current status
python3 slurm/status.py -t <task_name> --update
```

The `--monitor` flag runs in a loop: it checks node health, automatically
cancels and resubmits any stuck jobs, then waits N minutes before repeating.
It stops when all jobs are completed/failed, or on Ctrl+C.

### 3. Resubmit jobs

```bash
# Resubmit all failed jobs
python3 slurm/resubmit.py -t <task_name> --failed

# Resubmit failed jobs with more memory (e.g. for OOM)
python3 slurm/resubmit.py -t <task_name> --failed --memory 8gb

# Resubmit only failed jobs matching a sample pattern
python3 slurm/resubmit.py -t <task_name> --failed --sample "JetHT_Run2018"

# Detect and resubmit stuck jobs (checks node health, cancels, resubmits)
python3 slurm/resubmit.py -t <task_name> --stuck

# Resubmit specific jobs by ID
python3 slurm/resubmit.py -t <task_name> --jobs <job_id_1> <job_id_2>

# Dry run (show what would be resubmitted)
python3 slurm/resubmit.py -t <task_name> --failed --dry-run
```

## Common Failure Modes

### OUT_OF_MEMORY
Jobs killed by SLURM for exceeding the memory limit. Resubmit with `--memory 8gb` (or `16gb`).

### Corrupt input files
`TFile::Init: no keys recovered, file has been made a Zombie` -- the input ROOT file is broken. These need to be re-produced upstream.

### Filesystem transport errors
`Cannot send after transport endpoint shutdown` -- transient `/blue` filesystem glitch. The analysis may have succeeded but the output copy failed. Simply resubmit.

### Stuck jobs (hung node)

Jobs show as RUNNING in SLURM but produce no output. This happens when a
compute node's shared filesystem mount (e.g. `/blue`) hangs. Job processes
enter D-state (uninterruptible sleep) and never make progress, but SLURM
still reports them as RUNNING until they hit the walltime limit.

**How detection works:** `status.py --check-nodes` SSHs to each compute node
with running jobs and checks two things:
1. Whether the shared filesystem responds (`stat /blue`)
2. Whether all job processes are stuck in D-state with 0 CPU time

If both conditions indicate a problem, the jobs on that node are flagged as stuck.

**How to handle stuck jobs:**
```bash
# Detect stuck jobs
python3 slurm/status.py -t <task_name> --check-nodes

# Cancel and resubmit stuck jobs
python3 slurm/resubmit.py -t <task_name> --stuck

# Or use monitor mode to auto-detect and resubmit
python3 slurm/status.py -t <task_name> --monitor
```

## Job Lifecycle

```
submit.py                    status.py              resubmit.py
    |                            |                      |
    v                            v                      v
 prepared --> submitted --> running --> completed       |
                               |                       |
                               +--> failed ----------->+
                               +--> stuck (node issue)->+
```

The `manifest.json` tracks the full history including resubmit counts and timing.
