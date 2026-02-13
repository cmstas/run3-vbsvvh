#!/bin/bash

# Condor job executable for run3-vbsvvh preselection framework
#
# Arguments:
#   $1 - USER
#   $2 - N_CPUS
#   $3 - CONFIG_FILE (just filename, not full path)
#   $4 - JOB_OUTPUT_NAME
#   $5 - ANALYSIS
#   $6 - RUN_NUMBER
#   $7 - SAMPLE_NAME
#   $8 - JOB_IDX
#   $9+ - EXTRA_FLAGS (optional, e.g., --spanet_training, --spanet_infer)

USER=$1
N_CPUS=$2
CONFIG=$3
JOB_OUTPUT_NAME=$4
ANALYSIS=$5
RUN_NUMBER=$6
SAMPLE_NAME=$7
JOB_IDX=$8
shift 8
EXTRA_FLAGS="$@"

# Constants
OUTPUTDIR="output"
OUTPUTFILE="output"
OUTPUT_XRD="root://redirector.t2.ucsd.edu:1095//store/user/$USER/vbsvvh/preselection/"
CMSSW_VERSION='CMSSW_15_0_4'
MAX_RETRIES=5
SLEEP_DURATION="1m"

echo "=========================================="
echo "Job started at $(date)"
echo "=========================================="
echo "Arguments:"
echo "  USER=$USER"
echo "  N_CPUS=$N_CPUS"
echo "  CONFIG=$CONFIG"
echo "  JOB_OUTPUT_NAME=$JOB_OUTPUT_NAME"
echo "  ANALYSIS=$ANALYSIS"
echo "  RUN_NUMBER=$RUN_NUMBER"
echo "  SAMPLE_NAME=$SAMPLE_NAME"
echo "  JOB_IDX=$JOB_IDX"
echo "  EXTRA_FLAGS=$EXTRA_FLAGS"
echo "=========================================="

# Functions
function stageout {
    local COPY_SRC=$1
    local COPY_DEST=$2
    local retries=1
    local COPY_STATUS=1

    until [ $retries -ge $MAX_RETRIES ]; do
        echo "Stageout attempt $((retries+1)): gfal-copy ${COPY_SRC} ${COPY_DEST}"
        env -i X509_USER_PROXY=${X509_USER_PROXY} gfal-copy -p -f -t 7200 --verbose --checksum ADLER32 ${COPY_SRC} ${COPY_DEST}
        COPY_STATUS=$?
        if [ $COPY_STATUS -eq 0 ]; then
            echo "Successful stageout with $retries retries"
            break
        else
            echo "Failed stageout attempt $((retries+1))"
            retries=$((retries+1))
            echo "Sleeping for $SLEEP_DURATION"
            sleep $SLEEP_DURATION
        fi
    done

    if [ $COPY_STATUS -ne 0 ]; then
        echo "ERROR: gfal-copy failed with code $COPY_STATUS after $MAX_RETRIES attempts"
        return 1
    fi
    return 0
}

function source_environment {
    if [ -r "$OSGVO_CMSSW_Path/cmsset_default.sh" ]; then
        echo "Sourcing: $OSGVO_CMSSW_Path/cmsset_default.sh"
        source "$OSGVO_CMSSW_Path/cmsset_default.sh"
    elif [ -r "$OSG_APP/cmssoft/cms/cmsset_default.sh" ]; then
        echo "Sourcing: $OSG_APP/cmssoft/cms/cmsset_default.sh"
        source "$OSG_APP/cmssoft/cms/cmsset_default.sh"
    elif [ -r /cvmfs/cms.cern.ch/cmsset_default.sh ]; then
        echo "Sourcing: /cvmfs/cms.cern.ch/cmsset_default.sh"
        source /cvmfs/cms.cern.ch/cmsset_default.sh
    else
        echo "ERROR: Couldn't find cmsset_default.sh"
        exit 1
    fi
}

function setup_cmssw {
    echo "Setting up CMSSW environment: $CMSSW_VERSION"

    # Use pre-built CMSSW from CVMFS with ONNX Runtime already included
    cd /cvmfs/cms.cern.ch/el8_amd64_gcc12/cms/cmssw/$CMSSW_VERSION
    eval `scramv1 runtime -sh`
    cd -

    echo "CMSSW_BASE=$CMSSW_BASE"
    echo "SCRAM_ARCH=$SCRAM_ARCH"
}

# function setup_cmssw_local {
#     # Old function that creates local CMSSW area and compiles ONNXRuntime
#     # Kept for reference in case pre-built CVMFS version doesn't work
#     echo "Setting up CMSSW environment: $CMSSW_VERSION"
#
#     # Detect architecture
#     if [ -d /cvmfs/cms.cern.ch/el8_amd64_gcc12/cms/cmssw/$CMSSW_VERSION ]; then
#         export SCRAM_ARCH=el8_amd64_gcc12
#     elif [ -d /cvmfs/cms.cern.ch/el9_amd64_gcc12/cms/cmssw/$CMSSW_VERSION ]; then
#         export SCRAM_ARCH=el9_amd64_gcc12
#     else
#         echo "ERROR: Cannot find CMSSW $CMSSW_VERSION"
#         exit 1
#     fi
#
#     # Create local CMSSW area and build ONNXRuntime
#     echo "Creating local CMSSW area..."
#     cmsrel $CMSSW_VERSION
#     cd $CMSSW_VERSION/src
#     eval `scramv1 runtime -sh`
#
#     echo "Adding PhysicsTools/ONNXRuntime package..."
#     # Set temporary git config for cms-addpkg to work
#     git config --global user.name "Condor Job"
#     git config --global user.email "condor@job.local"
#     git config --global user.github "condorjob"
#
#     git cms-addpkg PhysicsTools/ONNXRuntime
#
#     # Verify the package was added
#     if [ ! -d "PhysicsTools/ONNXRuntime" ]; then
#         echo "ERROR: Failed to add PhysicsTools/ONNXRuntime"
#         exit 1
#     fi
#     ls -la PhysicsTools/ONNXRuntime/interface/
#
#     scram b -j $N_CPUS
#
#     cd - > /dev/null
#
#     echo "CMSSW_BASE=$CMSSW_BASE"
#     echo "SCRAM_ARCH=$SCRAM_ARCH"
# }

function run_analysis {
    echo "Running analysis..."

    # Check if --spanet_infer is in EXTRA_FLAGS - if so, disable multithreading and use batch size 518
    if [[ "$EXTRA_FLAGS" == *"--spanet_infer"* ]]; then
        echo "./bin/runAnalysis -b 518 -i $CONFIG -o $OUTPUTFILE --outdir $OUTPUTDIR --ana $ANALYSIS --run_number $RUN_NUMBER $EXTRA_FLAGS"
        ./bin/runAnalysis -b 518 -i $CONFIG -o $OUTPUTFILE --outdir $OUTPUTDIR --ana $ANALYSIS --run_number $RUN_NUMBER $EXTRA_FLAGS
    else
        echo "./bin/runAnalysis -n $N_CPUS -i $CONFIG -o $OUTPUTFILE --outdir $OUTPUTDIR --ana $ANALYSIS --run_number $RUN_NUMBER $EXTRA_FLAGS"
        ./bin/runAnalysis -n $N_CPUS -i $CONFIG -o $OUTPUTFILE --outdir $OUTPUTDIR --ana $ANALYSIS --run_number $RUN_NUMBER $EXTRA_FLAGS
    fi
    return $?
}

function validate_root_file {
    # Validate ROOT file integrity by checking all branches for all events
    # Returns 0 if file is valid, 1 if corrupted
    local FILEPATH=$1
    local TREENAME=${2:-"Events"}

    echo "Validating ROOT file: $FILEPATH (tree: $TREENAME)"

    python3 << EOF
import sys
try:
    import ROOT
    ROOT.gErrorIgnoreLevel = ROOT.kError  # Suppress warnings

    filepath = "${FILEPATH}"
    treename = "${TREENAME}"

    # Check file can be opened
    f = ROOT.TFile.Open(filepath)
    if not f or f.IsZombie():
        print("[VALIDATE] ERROR: Cannot open file or file is zombie")
        sys.exit(1)

    # Check tree exists
    t = f.Get(treename)
    if not t:
        print(f"[VALIDATE] ERROR: Tree '{treename}' not found in file")
        f.Close()
        sys.exit(1)

    nevts = t.GetEntries()
    print(f"[VALIDATE] Found {nevts} entries in tree '{treename}'")

    if nevts == 0:
        print("[VALIDATE] WARNING: Tree has 0 entries (may be expected for some samples)")
        f.Close()
        sys.exit(0)

    # Check all events - GetEntry returns -1 on I/O error
    nbad = 0
    for i in range(nevts):
        if t.GetEntry(i) < 0:
            print(f"[VALIDATE] ERROR: Bad event at index {i}")
            nbad += 1
            if nbad >= 10:
                print("[VALIDATE] ERROR: Too many bad events, aborting check")
                break

    f.Close()

    if nbad > 0:
        print(f"[VALIDATE] FAILED: Found {nbad} corrupted events")
        sys.exit(1)
    else:
        print("[VALIDATE] PASSED: All events readable")
        sys.exit(0)

except Exception as e:
    print(f"[VALIDATE] ERROR: Exception during validation: {e}")
    sys.exit(1)
EOF
    return $?
}

# Main script
echo ""
echo "=== Checking proxy ==="
PROXY_TIMELEFT=$(voms-proxy-info -timeleft 2>&1)
PROXY_STATUS=$?

MIN_PROXY_TIME=$((6 * 3600))  # 6 hours in seconds (21600)

if [ $PROXY_STATUS -ne 0 ] || [ -z "$PROXY_TIMELEFT" ] || [ "$PROXY_TIMELEFT" -le 0 ] 2>/dev/null; then
    echo "ERROR: No valid proxy found or proxy has expired"
    echo "Proxy check output: $PROXY_TIMELEFT"
    echo "Please run: voms-proxy-init -voms cms -valid 192:00"
    exit 1
fi

if [ "$PROXY_TIMELEFT" -lt "$MIN_PROXY_TIME" ]; then
    echo "ERROR: Proxy expires in less than 6 hours (${PROXY_TIMELEFT} seconds remaining)"
    echo "Please renew your proxy: voms-proxy-init -voms cms -valid 192:00"
    exit 1
fi

echo "Proxy time remaining: ${PROXY_TIMELEFT} seconds (~$((PROXY_TIMELEFT / 3600)) hours)"

echo ""
echo "=== Setting up CMS environment ==="
source_environment
setup_cmssw

echo ""
echo "=== Extracting package ==="
## Extract and build inside CMSSW src directory so headers are found
#WORK_DIR=$CMSSW_BASE/src/Analysis

# Extract to local working directory (CMSSW env is set up from CVMFS)
WORK_DIR=$(pwd)/Analysis
mkdir -p $WORK_DIR
tar -xzf package.tar.gz -C $WORK_DIR
cp $CONFIG $WORK_DIR/.
cd $WORK_DIR
ls -la

echo ""
echo "=== Building analysis code ==="
make -j$N_CPUS

if [ $? -ne 0 ]; then
    echo "ERROR: Compilation failed"
    exit 1
fi

echo ""
echo "=== Directory contents before analysis ==="
ls -la
ls -la bin/

echo ""
echo "=== Running analysis ==="
run_analysis
ANALYSIS_STATUS=$?

if [ $ANALYSIS_STATUS -ne 0 ]; then
    echo "ERROR: Analysis failed twice. Aborting."
    exit 1
fi


# Resubmission is now handled via condor/resubmit.py to keep track of resubmits 
# if [ $ANALYSIS_STATUS -ne 0 ]; then
#     echo "Analysis failed with status $ANALYSIS_STATUS, retrying..."
#     run_analysis
#     ANALYSIS_STATUS=$?
#     if [ $ANALYSIS_STATUS -ne 0 ]; then
#         echo "ERROR: Analysis failed twice. Aborting."
#         exit 1
#     fi
# fi

echo ""
echo "=== Directory contents after analysis ==="
ls -la
ls -la $OUTPUTDIR/ 2>/dev/null || echo "Output directory not found"

# Check if output file exists
if [ ! -f "$OUTPUTDIR/$OUTPUTFILE.root" ]; then
    echo "ERROR: Output file not found: $OUTPUTDIR/$OUTPUTFILE.root"
    exit 1
fi

echo ""
echo "=== Validating output ROOT file ==="
validate_root_file "$OUTPUTDIR/$OUTPUTFILE.root" "Events"
if [ $? -ne 0 ]; then
    echo "ERROR: Output ROOT file validation failed"
    echo "Removing corrupted file: $OUTPUTDIR/$OUTPUTFILE.root"
    rm -f "$OUTPUTDIR/$OUTPUTFILE.root"
    exit 1
fi

echo ""
echo "=== Staging out results ==="

# Stage out to: OUTPUT_XRD/JOB_OUTPUT_NAME/SAMPLE_NAME/output_JOB_IDX.root
COPY_SRC="file://$(pwd)/$OUTPUTDIR/$OUTPUTFILE.root"
COPY_DEST="${OUTPUT_XRD}/${JOB_OUTPUT_NAME}/${SAMPLE_NAME}/output_${JOB_IDX}.root"

echo "Copying: ${COPY_SRC} -> ${COPY_DEST}"
stageout "$COPY_SRC" "$COPY_DEST"

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to stage out output file"
    exit 1
fi

# Verify output file exists at the destistation path
echo "=== Verifying output file on XRootD ==="
env -i X509_USER_PROXY=${X509_USER_PROXY} gfal-stat "${COPY_DEST}" > /dev/null 2>&1
if [ $? -ne 0 ]; then
    echo "ERROR: Output verification failed - file not found on XRootD: ${COPY_DEST}"
    exit 1
fi
echo "Output file verified successfully"

# Stage out cutflow file if it exists
if [ -f "$OUTPUTDIR/${OUTPUTFILE}_cutflow.txt" ]; then
    COPY_SRC="file://$(pwd)/$OUTPUTDIR/${OUTPUTFILE}_cutflow.txt"
    COPY_DEST="${OUTPUT_XRD}/${JOB_OUTPUT_NAME}/${SAMPLE_NAME}/output_${JOB_IDX}_cutflow.txt"
    echo "Copying cutflow: ${COPY_SRC} -> ${COPY_DEST}"
    stageout "$COPY_SRC" "$COPY_DEST"
fi

echo ""
echo "=========================================="
echo "Job completed successfully at $(date)"
echo "=========================================="
