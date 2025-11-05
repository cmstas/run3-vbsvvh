#ifndef CORRECTIONS_H
#define CORRCTIONS_H

#pragma once

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RVec.hxx"

#include "correction_run3.h"

#include "TRandom3.h"

using correction::CorrectionSet;
using RNode = ROOT::RDF::RNode;
using ROOT::VecOps::RVec;

namespace Run3
{

float looseDFBtagWP(std::string year);
float mediumDFBtagWP(std::string year);
float tightDFBtagWP(std::string year);

// RNode METPhiCorrections(correction::CorrectionSet cset_met_2022, correction::CorrectionSet cset_met_2022EE, RNode df);
// const auto cset_met_2022 = *CorrectionSet::from_file("corrections/met/2022.json");
// const auto cset_met_2022EE = *CorrectionSet::from_file("corrections/met/2022EE.json");

// RNode METUnclusteredCorrections(RNode df, std::string variation);

// // jet mass scale and resolution corrections
// RNode JMR_Corrections(correction::CorrectionSet cset_jet_mass_scale, RNode df, std::string variation); 
// const auto cset_jmr = *CorrectionSet::from_file("corrections/scalefactors/particlenet/jmar.json");

// RNode JMS_Corrections(correction::CorrectionSet cset_jet_mass_scale, RNode df, std::string variation);
// const auto cset_jms = *CorrectionSet::from_file("corrections/scalefactors/particlenet/jmar.json");

// RNode JetEnergyCorrection(correction::CorrectionSet cset_jerc_2022, correction::CorrectionSet cset_jerc_2022EE, RNode df, std::string JEC_type, std::string variation);
// RNode JetEnergyResolution(correction::CorrectionSet cset_jerc_2022, correction::CorrectionSet cset_jerc_2022EE, correction::CorrectionSet cset_jer_smear, RNode df, std::string variation);
// const auto cset_jerc_2022 = *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/2022_Summer22/jet_jerc.json.gz");
// const auto cset_jerc_2022EE = *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/2022_Summer22EE/jet_jerc.json.gz");

// const auto cset_jer_smear = *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/jer_smear.json.gz");

// RNode JetVetoMaps(correction::CorrectionSet cset_jetveto_2022, correction::CorrectionSet cset_jetveto_2022EE, RNode df);
// const auto cset_jetveto_2022 = *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/2022_Summer22/jetvetomaps.json.gz");
// const auto cset_jetveto_2022EE = *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/2022_Summer22EE/jetvetomaps.json.gz");

} // end namespace Run3

#endif