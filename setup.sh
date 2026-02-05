CMSSW_VERSION=CMSSW_15_0_4
echo "Setting up CMSSW environment: $CMSSW_VERSION"
# Use pre-built CMSSW from CVMFS with ONNX Runtime already included
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /cvmfs/cms.cern.ch/el8_amd64_gcc12/cms/cmssw/$CMSSW_VERSION
eval `scramv1 runtime -sh`
cd -
echo "CMSSW_BASE=$CMSSW_BASE"
echo "SCRAM_ARCH=$SCRAM_ARCH"