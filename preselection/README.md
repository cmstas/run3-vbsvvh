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
# Set up environment
source setup.sh

# Compile
cd preselection/
make -j8

# Run examples in `run_wrapper.sh`, e.g., for running over locally over one signal sample on UAF:
python run_rdf.py etc/input_sample_jsons/run2/sig_c2v1p0_c3_1p0/all_events/2017_VBSWZH_c2v1p0_c3_1p0.json --prefix /ceph/cms/ -o some_dir -n test_small -a all_events -m local -r 2
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
The inputs for the preselection framework are defined using configuration JSON files. The config files define input samples, cross-sections, and metadata in JSON format for RDataFrame's `FromSpec`. The json files for all of the skim sets are provided in the `preselection/etc/input_sample_jsons` directory.  

**See [`etc/README.md`](etc/README.md) for detailed documentation.**

## Notes on running 
To run the preselection, use the compiled binary and provide the input specification and output file. The `run_rdf.py` serves as a wrapper around the `bin/runAnalysis`. This script will:
* **Prepare the json files**: The scrip will prepend a given prefix to the paths in all of the json files (to accommodate running either with direct file paths on UAF, or with xrd). 
* **Merge the given json files**: The underlying analysis code expects to receive a single json, so this script will merge all jsons into one single json. Note that `bin/runAnalysis` expects all samples to be the same kind (either signal, background, or data), so if running locally do not mix multiple kinds of jsons into a single run. The final json is saved locally to the `merged_jsons` dir for reference. 
* **Run the analysis**: It will either directly run the `bin/runAnalysis` (if the mode is local) or will wrap around the condor submission (if the mode is condor).

There are examples in `run_wrapper.sh` for running `run_rdf.py` either locally over single jsons for tests, or for running the full analysis at scale over all skim selections. 

---
## Details on the condor batch submission

For large-scale processing, the HTCondor submission should be used. Some overview information is provided here, with more details in the condor README. 

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
