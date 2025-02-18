#ifndef WEIGHTS_H
#define WEIGHTS_H

#pragma once

#include <unordered_map>
#include <string>

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RVec.hxx"

#include "utils.h"
#include "correction.h"

using correction::CorrectionSet;
using RNode = ROOT::RDF::RNode;
using ROOT::VecOps::RVec;


/*
############################################
GOLDEN JSON
############################################
*/
RNode applyGoldenJSONWeight(lumiMask golden, RNode df);
const auto LumiMask = lumiMask::fromJSON("corrections/goldenJson/Cert_Collisions2022_355100_362760_Golden.json");

/*
############################################
PILEUP SFs
############################################
*/
const std::unordered_map<std::string, correction::CorrectionSet> pileupScaleFactors = {
    {"2022Re-recoBCD", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/LUM/2022_Summer22/puWeights.json.gz")},
    {"2022Re-recoE+PromptFG", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/LUM/2022_Summer22EE/puWeights.json.gz")}
};
const std::unordered_map<std::string, std::string> pileupScaleFactors_yearmap = {
    {"2022Re-recoBCD", "Collisions2022_355100_357900_eraBCD_GoldenJson"},
    {"2022Re-recoE+PromptFG", "Collisions2022_359022_362760_eraEFG_GoldenJson"}
};
RNode applyPileupScaleFactors(std::unordered_map<std::string, correction::CorrectionSet> cset_pileup, std::unordered_map<std::string, std::string> year_map, RNode df);

/*
############################################
MUON SFs
############################################
*/
const std::unordered_map<std::string, correction::CorrectionSet> muonScaleFactors = {
    {"2022Re-recoBCD", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2022_Summer22/muon_Z.json.gz")},
    {"2022Re-recoE+PromptFG", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2022_Summer22EE/muon_Z.json.gz")}
};

const std::unordered_map<std::string, std::string> muonIDScaleFactors_yearmap = {
    {"2022Re-recoBCD", "NUM_TightID_DEN_TrackerMuons"},
    {"2022Re-recoE+PromptFG", "NUM_TightID_DEN_TrackerMuons"}
};
RNode applyMuonIDScaleFactors(std::unordered_map<std::string, correction::CorrectionSet> cset_muon, std::unordered_map<std::string, std::string> year_map, RNode df);

const std::unordered_map<std::string, std::string> muonRecoScaleFactors_yearmap = {
    {"2022Re-recoBCD", "NUM_TightPFIso_DEN_TightID"},
    {"2022Re-recoE+PromptFG", "NUM_TightPFIso_DEN_TightID"}
};
RNode applyMuonRecoScaleFactors(std::unordered_map<std::string, correction::CorrectionSet> cset_muon, std::unordered_map<std::string, std::string> year_map, RNode df);

const std::unordered_map<std::string, std::string> muonTriggerScaleFactors_yearmap = {
    {"2022Re-recoBCD", "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight"},
    {"2022Re-recoE+PromptFG", "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight"}
};
RNode applyMuonTriggerScaleFactors(std::unordered_map<std::string, correction::CorrectionSet> cset_muon, std::unordered_map<std::string, std::string> year_map, RNode df);

/*
############################################
ELECTRON SFs
############################################
*/
const std::unordered_map<std::string, correction::CorrectionSet> electronScaleFactors = {
    {"2022Re-recoBCD", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2022_Summer22/electron.json.gz")},
    {"2022Re-recoE+PromptFG", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2022_Summer22EE/electron.json.gz")}
};
const std::unordered_map<std::string, std::string> electronScaleFactors_yearmap = {
    {"2022Re-recoBCD", "Electron-ID-SF"},
    {"2022Re-recoE+PromptFG", "Electron-ID-SF"}
};
RNode applyElectronIDScaleFactors(std::unordered_map<std::string, correction::CorrectionSet> cset_electron, std::unordered_map<std::string, std::string> year_map, RNode df);
RNode applyElectronRecoScaleFactors(std::unordered_map<std::string, correction::CorrectionSet> cset_electron, std::unordered_map<std::string, std::string> year_map, RNode df);

const std::unordered_map<std::string, correction::CorrectionSet> electronTriggerScaleFactors = {
    {"2022Re-recoBCD", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2022_Summer22/electronHlt.json.gz")},
    {"2022Re-recoE+PromptFG", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2022_Summer22EE/electronHlt.json.gz")}
};
const std::unordered_map<std::string, std::string> electronTriggerScaleFactors_yearmap = {
    {"2022Re-recoBCD", "Electron-HLT-SF"},
    {"2022Re-recoE+PromptFG", "Electron-HLT-SF"}
};
RNode applyElectronTriggerScaleFactors(std::unordered_map<std::string, correction::CorrectionSet> cset_electron, std::unordered_map<std::string, std::string> year_map, RNode df);

/*
############################################
PARTICLE NET SFs
############################################
*/

// particle net sfs
// std::unordered_map<std::string, correction::CorrectionSet> PNETWScaleFactors = {
//     {"2022", *CorrectionSet::from_file("corrections/scalefactors/particlenet/pnet.json")},
//     {"2022EE", *CorrectionSet::from_file("corrections/scalefactors/particlenet/pnet.json")}
// };
// RNode applyPNETWScaleFactors(std::unordered_map<std::string, correction::CorrectionSet> cset_pnet_w, RNode df);


/*
############################################
B-TAGGING SFs
############################################
*/
// std::unordered_map<std::string, correction::CorrectionSet> bTaggingScaleFactors = {
//     {"2022", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/2022_Summer22/btagging.json.gz")},
//     {"2022EE", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/2022_Summer22EE/btagging.json.gz")},
//     {"eff", *CorrectionSet::from_file("corrections/scalefactors/btagging/btag_eff.json")}
// };
// RNode applyBTaggingScaleFactors(std::unordered_map<std::string, correction::CorrectionSet> cset_btag, RNode df);

/*
############################################
OTHER SFs
############################################
*/

RNode applyL1PreFiringReweighting(RNode df);
RNode applyPSWeight_FSR(RNode df);
RNode applyPSWeight_ISR(RNode df);
RNode applyLHEScaleWeight_muF(RNode df);
RNode applyLHEScaleWeight_muR(RNode df);
RNode applyLHEWeights_pdf(RNode df);

RNode applyDataWeights(RNode df);
RNode applyMCWeights(RNode df);

#endif