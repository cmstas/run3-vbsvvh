# Condor Job Submission 

This directory contains scripts for submitting and managing HTCondor batch jobs for the VBS VVH preselection framework.

## Prerequisites

1. Valid grid proxy:
   ```bash
   voms-proxy-init -voms cms -valid 192:00
   ```

2. Config file in `etc/` directory (e.g., `etc/0Lep2FJ_run3-sig.json`)

3. Confirm analysis code compiles locally.


## Quick Start

```bash
cd preselection/

# 1. Submit jobs
python condor/submit.py -c etc/0Lep2FJ_run3-sig.json -a 0Lep2FJ -r 3

# 2. Check status
python condor/status.py --task 0Lep2FJ_run3-sig_0Lep2FJ
```

## Submission

### Basic Usage

```bash
python condor/submit.py -c <config> -a <analysis> -r <run_number> [options]
```

### Required Arguments

| Argument | Description |
|----------|-------------|
| `-c`, `--config` | Path to config JSON file (relative to `preselection/`) |
| `-a`, `--analysis` | Analysis channel, e.g.: `0Lep2FJ`, `1Lep2FJ` |
| `-r`, `--run_number` | Run number: `2` or `3` |

### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `-t`, `--tag` | (none) | Tag for output directory naming |
| `-n`, `--ncpus` | 4 | CPUs per job |
| `-m`, `--memory` | 2G | Memory per job |
| `--files-per-job` | 10 | Input files per job |
| `--sample` | (all) | Regex to filter samples |
| `--dry-run` | false | Prepare jobs without submitting |
| `--monitor` | false | Monitor jobs and auto-resubmit failed ones |
| `--max-resubmits` | 3 | Maximum resubmit attempts per job |
| `--timeout` | (none) | Stop monitoring after N hours |
| `--interval` | 5 | Check interval in minutes |

### Examples

```bash
# Submit all samples in config
python condor/submit.py -c etc/0Lep2FJ_run3-sig.json -a 0Lep2FJ -r 3

# Submit with custom tag and fewer files per job
python condor/submit.py -c etc/0Lep2FJ_run3-sig.json -a 0Lep2FJ -r 3 -t v2 --files-per-job 5

# Submit only specific samples (regex match)
python condor/submit.py -c etc/0Lep2FJ_run3-bkg.json -a 0Lep2FJ -r 3 --sample "QCD"

# Dry run to inspect job setup
python condor/submit.py -c etc/0Lep2FJ_run3-sig.json -a 0Lep2FJ -r 3 --dry-run

# Submit with automatic monitoring (runs in foreground)
python condor/submit.py -c etc/0Lep2FJ_run3-sig.json -a 0Lep2FJ -r 3 --monitor --timeout 12
```

### What Submission Creates

```
condor/jobs/<task_name>/
├── manifest.json          # Job tracking metadata
├── package.tar.gz         # Analysis code tarball
├── executable.sh          # Job executable
├── <sample_name_0>/       # Per-job directories
│   ├── config.json        # Job-specific config
│   ├── submit.cmd         # Condor submit file
│   ├── job.*.log          # Condor log
│   ├── job.*.stdout       # Job stdout
│   └── job.*.stderr       # Job stderr
└── <sample_name_1>/
    └── ...
```

## Checking Status

### Basic Usage

```bash
python condor/status.py --task <task_name> [options]
```

### Examples

```bash
# List all tasks
python condor/status.py --list

# Check status of a task
python condor/status.py --task 0Lep2FJ_run3-sig_0Lep2FJ

# Show status grouped by sample
python condor/status.py --task 0Lep2FJ_run3-sig_0Lep2FJ --by-sample

# Show failed job details (stderr output)
python condor/status.py --task 0Lep2FJ_run3-sig_0Lep2FJ --failed

# Update manifest with current status
python condor/status.py --task 0Lep2FJ_run3-sig_0Lep2FJ --update

# Skip condor queries (use logs only)
python condor/status.py --task 0Lep2FJ_run3-sig_0Lep2FJ --no-condor
```

### Job States

| Status | Description |
|--------|-------------|
| `prepared` | Job created but not submitted |
| `queued` | Waiting in Condor queue |
| `running` | Currently executing |
| `held` | Job held (check logs for error) |
| `completed` | Finished successfully (exit code 0) |
| `failed` | Finished with error (non-zero exit) |
| `unknown` | Status cannot be determined |

## Automatic Monitoring

The `--monitor` flag enables automatic job monitoring and resubmission. This is useful for long-running submissions where you want failed jobs to be automatically retried.

### Usage

Since monitoring runs continuously, use a `screen` or `tmux` session:

```bash
# Start a screen session
screen -S condor_monitor

# Submit with monitoring
python condor/submit.py -c etc/0Lep2FJ_run3-sig.json -a 0Lep2FJ -r 3 --monitor --timeout 12 --max-resubmits 3

# Detach from screen: Ctrl+A, D
# Reconnect later: screen -r condor_monitor
```

### Monitor Options

| Argument | Default | Description |
|----------|---------|-------------|
| `--monitor` | false | Enable automatic monitoring |
| `--max-resubmits` | 3 | Maximum resubmit attempts per job |
| `--timeout` | (none) | Stop monitoring after N hours |
| `--interval` | 5 | Check interval in minutes |

### Behavior

- Checks job status every `--interval` minutes
- Automatically resubmits failed jobs (up to `--max-resubmits` times)
- Only prints output when resubmitting jobs
- Stops when: all jobs complete, timeout reached, or all remaining jobs exceed max resubmits
- Tracks `resubmit_count` per job in manifest

## Resubmitting Failed Jobs (Manual)

For manual resubmission of failed jobs:

```bash
python condor/resubmit.py --task <task_name> [options]
```

### Options

```bash
# Resubmit all failed jobs
python condor/resubmit.py --task 0Lep2FJ_run3-sig_0Lep2FJ --failed

# Resubmit specific samples
python condor/resubmit.py --task 0Lep2FJ_run3-sig_0Lep2FJ --failed --sample "QCD"

# Dry run
python condor/resubmit.py --task 0Lep2FJ_run3-sig_0Lep2FJ --failed --dry-run
```

## Output Location

Job outputs are staged to XRootD:
```
root://redirector.t2.ucsd.edu:1095//store/user/<USER>/vbsvvh/preselection/run3-vbsvvh/<task_name>/<sample_name>/output_<job_idx>.root
```

## Troubleshooting

To troubleshoot, you can check directly individual job logs:
```bash
# View stdout
cat condor/jobs/<task_name>/<sample>/job.*.stdout

# View stderr
cat condor/jobs/<task_name>/<sample>/job.*.stderr

# View condor log
cat condor/jobs/<task_name>/<sample>/job.*.log
```

Use standar Condor commands for manual job control:
```bash
# Remove all jobs for a task
condor_rm <cluster_id>

# Release held jobs
condor_release <cluster_id>

# Check job details
condor_q -l <cluster_id>
```
