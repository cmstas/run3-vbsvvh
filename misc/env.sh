#!/bin/bash

# OBSOLETE: use setup.sh instead for faster setup

CMSSW_VERSION="CMSSW_15_0_4"

source /cvmfs/cms.cern.ch/cmsset_default.sh

if [ ! -d "$CMSSW_VERSION" ]; then
  echo "Creating CMSSW environment: $CMSSW_VERSION"
  cmsrel $CMSSW_VERSION
fi

cd $CMSSW_VERSION/src
eval `scramv1 runtime -sh`

git cms-addpkg PhysicsTools/ONNXRuntime
scram b -j 8
cd -