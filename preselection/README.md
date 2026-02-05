# Preselection Framework

C++ framework for VBS VVH preselection using RDataFrame. The paths provided in this README assume you are running on the UAF cluster.

## Directory Structure

```
preselection/
├── src/           # C++ source files
├── include/       # Header files
├── etc/           # Config generation scripts and JSON files
├── condor/        # HTCondor batch submission scripts
├── corrections/   # Scale factors and corrections
├── bin/           # Compiled binaries (generated)
└── build/         # Object files (generated)
```

## Quick Start

```bash
# 1. Set up environment
source setup.sh

# 2. Compile
cd preselection/
make -j8

# 3. Generate config
python3 etc/make_config.py --channel 0Lep2FJ_run3 --category sig

# 4. Run locally (for testing)
bin/runAnalysis -i etc/0Lep2FJ_run3-sig.json -a 0Lep2FJ -n 8 --run_number 3

# 5. Or submit to Condor (for production)
python condor/submit.py -c etc/0Lep2FJ_run3-sig.json -a 0Lep2FJ -r 3
```

## Prerequisites

Included in CMSSW 15_0_4:
- ROOT with RDataFrame
- Python 3.x
- Uproot
- Boost (for GoldenJSON)
- Correctionlib

Environment setup script: `setup.sh`

## Compilation
To set up and compile the preselection framework, navigate to the top level directory and run:

```bash
source setup.sh
cd preselection/
make -j8
```
This will compile the source files and place the binary in the bin/ directory.

## Configuration Files
The inputs for the preselection framework are defined using configuration JSON files. The config files define input samples, cross-sections, and metadata in JSON format for RDataFrame's `FromSpec`. Specifics will have to be adapted by the analyzer.

### Generating Configs
To generate a configuration file, run the make-config.py script with the appropriate category and channel name.

```bash
python3 etc/make_config.py --category <category> --channel <channel>
```

For example, to generate a configuration file for signal samples using the skims for the channel `0Lep2FJ_run3`:

```bash
python3 etc/make_config.py --channel 0Lep2FJ_run3 --category sig
```

This creates `0Lep2FJ_run3-sig.json` in the `etc/` directory.

**See [`etc/README.md`](etc/README.md) for detailed documentation.**

## Running Locally
To run the preselection, use the compiled binary and provide the input specification and output file.

For testing or small samples:

```bash
bin/runAnalysis -i <config.json> -a <channel> -n <threads> --run_number <2|3>
```

Example:
```bash
bin/runAnalysis -i etc/0Lep2FJ_run3-sig.json -a 0Lep2FJ -n 8 --run_number 3
```

Full options:
```bash
bin/runAnalysis --help
```

---
## Condor Batch Submission

For production processing, submit jobs to HTCondor:

```bash
# Submit jobs
python condor/submit.py -c <config.json> -a <channel> -r <run_number>

# Check status
python condor/status.py --task <task_name>

# Resubmit failed jobs
python condor/resubmit.py --task <task_name> --failed
```

For automatic monitoring and resubmission:
```bash
screen -S condor_monitor
python condor/submit.py -c <config.json> -a <channel> -r <run_number> --monitor --timeout 12
# Ctrl+A, D to detach
```

**See [`condor/README.md`](condor/README.md) for detailed documentation.**
