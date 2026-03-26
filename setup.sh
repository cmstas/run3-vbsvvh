CMSSW_VERSION=CMSSW_16_1_0_pre2
echo "Setting up CMSSW environment: $CMSSW_VERSION"
# Use pre-built CMSSW from CVMFS with ONNX Runtime already included
source /cvmfs/cms.cern.ch/cmsset_default.sh
. /etc/os-release
if [[ $PLATFORM_ID == *"el9"* ]]; then
    cd /cvmfs/cms.cern.ch/el9_amd64_gcc13/cms/cmssw/$CMSSW_VERSION
elif [[ $PLATFORM_ID == *"el8"* ]]; then
    cd /cvmfs/cms.cern.ch/el8_amd64_gcc13/cms/cmssw/$CMSSW_VERSION
fi
eval `scramv1 runtime -sh`
cd -
echo "CMSSW_BASE=$CMSSW_BASE"
echo "SCRAM_ARCH=$SCRAM_ARCH"
