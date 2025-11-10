#ifndef WEIGHTS_RUN3_H
#define WEIGHTS_RUN3_H

#pragma once

#include <unordered_map>
#include <string>
#include <iostream>

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RVec.hxx"

#include "utils.h"
#include "correction.h"

using correction::CorrectionSet;
using RNode = ROOT::RDF::RNode;
using ROOT::VecOps::RVec;

namespace Run3
{

/*
############################################
GOLDEN JSON
############################################
*/
RNode applyGoldenJSONWeight(const lumiMask& golden, RNode df);
const auto LumiMask = lumiMask::fromJSON({"corrections/goldenJson/Cert_Collisions2022_355100_362760_Golden.json", 
    "corrections/goldenJson/Cert_Collisions2023_366442_370790_Golden.json", 
    "corrections/goldenJson/Cert_Collisions2024_378981_386951_Golden.json",
    "corrections/goldenJson/Cert_271036-284044_13TeV_Legacy2016_Collisions16.json",
    "corrections/goldenJson/Cert_294927-306462_13TeV_UL2017_Collisions17_Golden.json",
    "corrections/goldenJson/Cert_314472-325175_13TeV_Legacy2018_Collisions18.json"
});

/*
############################################
PILEUP SFs
############################################
*/
const std::unordered_map<std::string, correction::CorrectionSet> pileupScaleFactors = {
    {"2016preVFP", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/LUM/2016preVFP_UL/puWeights.json.gz")},
    {"2016postVFP", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/LUM/2016postVFP_UL/puWeights.json.gz")},
    {"2017", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/LUM/2017_UL/puWeights.json.gz")},
    {"2018", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/LUM/2018_UL/puWeights.json.gz")},
    {"2022Re-recoBCD", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/LUM/2022_Summer22/puWeights.json.gz")},
    {"2022Re-recoE+PromptFG", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/LUM/2022_Summer22EE/puWeights.json.gz")}
};
const std::unordered_map<std::string, std::string> pileupScaleFactors_yearmap = {
    {"2016preVFP", "Collisions16_UltraLegacy_goldenJSON"},
    {"2016postVFP", "Collisions16_UltraLegacy_goldenJSON"},
    {"2017", "Collisions17_UltraLegacy_goldenJSON"},
    {"2018", "Collisions18_UltraLegacy_goldenJSON"},
    {"2022Re-recoBCD", "Collisions2022_355100_357900_eraBCD_GoldenJson"},
    {"2022Re-recoE+PromptFG", "Collisions2022_359022_362760_eraEFG_GoldenJson"}
};
RNode applyPileupScaleFactors(std::unordered_map<std::string, correction::CorrectionSet> cset_pileup, std::unordered_map<std::string, std::string> year_map, RNode df);

/*
############################################
PILEUP ID SFs
############################################
*/
const std::unordered_map<std::string, correction::CorrectionSet> pileupScaleFactors = {
    {"2016preVFP", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/2016preVFP_UL/jmar.json.gz")},
    {"2016postVFP", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/2016postVFP_UL/jmar.json.gz")},
    {"2017", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/2017_UL/jmar.json.gz")},
    {"2018", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/2018_UL/jmar.json.gz")},
};
const std::unordered_map<std::string, std::string> pileupScaleFactors_yearmap = {
    {"2016preVFP", "PUJetID_eff"},
    {"2016postVFP", "PUJetID_eff"},
    {"2017", "PUJetID_eff"},
    {"2018", "PUJetID_eff"},
};
RNode applyPileupIDScaleFactors(std::unordered_map<std::string, correction::CorrectionSet> cset_pileup, std::unordered_map<std::string, std::string> year_map, RNode df);

/*
############################################
MUON SFs
############################################
*/
const std::unordered_map<std::string, correction::CorrectionSet> muonScaleFactors = {
    {"2016preVFP", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2016preVFP_UL/muon_Z.json.gz")},
    {"2016postVFP", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2016postVFP_UL/muon_Z.json.gz")},
    {"2017", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2017_UL/muon_Z.json.gz")},
    {"2018", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2018_UL/muon_Z.json.gz")},
    {"2022Re-recoBCD", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2022_Summer22/muon_Z.json.gz")},
    {"2022Re-recoE+PromptFG", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2022_Summer22EE/muon_Z.json.gz")},
    // {"2024Prompt", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/2024_Summer24/muon_Z.json.gz")}
};

const std::unordered_map<std::string, std::string> muonIDScaleFactors_yearmap = {
    {"2016preVFP", "NUM_TightID_DEN_TrackerMuons"},
    {"2016postVFP", "NUM_TightID_DEN_TrackerMuons"},
    {"2017", "NUM_TightID_DEN_TrackerMuons"},
    {"2018", "NUM_TightID_DEN_TrackerMuons"},
    {"2022Re-recoBCD", "NUM_TightID_DEN_TrackerMuons"},
    {"2022Re-recoE+PromptFG", "NUM_TightID_DEN_TrackerMuons"},
    // {"2024Prompt", "NUM_TightID_DEN_TrackerMuons"}
};
RNode applyMuonIDScaleFactors(std::unordered_map<std::string, correction::CorrectionSet> cset_muon, std::unordered_map<std::string, std::string> year_map, RNode df);

const std::unordered_map<std::string, std::string> muonRecoScaleFactors_yearmap = {
    {"2022Re-recoBCD", "NUM_TightPFIso_DEN_TightID"},
    {"2022Re-recoE+PromptFG", "NUM_TightPFIso_DEN_TightID"},
    // {"2024Prompt", "NUM_TightPFIso_DEN_TightID"}
};
RNode applyMuonRecoScaleFactors(std::unordered_map<std::string, correction::CorrectionSet> cset_muon, std::unordered_map<std::string, std::string> year_map, RNode df);

const std::unordered_map<std::string, std::string> muonTriggerScaleFactors_yearmap = {
    {"2016preVFP", "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight"},
    {"2016postVFP", "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight"},
    {"2017", "NUM_IsoMu27_DEN_CutBasedIdTight_and_PFIsoTight"},
    {"2018", "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight"},
    {"2022Re-recoBCD", "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight"},
    {"2022Re-recoE+PromptFG", "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight"},
    // {"2024Prompt", "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight"}
};
RNode applyMuonTriggerScaleFactors(std::unordered_map<std::string, correction::CorrectionSet> cset_muon, std::unordered_map<std::string, std::string> year_map, RNode df);

/*
############################################
ELECTRON SFs
############################################
*/
const std::unordered_map<std::string, correction::CorrectionSet> electronScaleFactors = {
    {"2016preVFP", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2016preVFP_UL/electron.json.gz")},
    {"2016postVFP", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2016postVFP_UL/electron.json.gz")},
    {"2017", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2017_UL/electron.json.gz")},
    {"2018", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2018_UL/electron.json.gz")},
    {"2022Re-recoBCD", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2022_Summer22/electron.json.gz")},
    {"2022Re-recoE+PromptFG", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2022_Summer22EE/electron.json.gz")},
    // {"2024Prompt", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/2024_Summer24/electron_v1.json.gz")}
};
const std::unordered_map<std::string, std::string> electronScaleFactors_yearmap = {
    {"2016preVFP", "UL-Electron-ID-SF"},
    {"2016postVFP", "UL-Electron-ID-SF"},
    {"2017", "UL-Electron-ID-SF"},
    {"2018", "UL-Electron-ID-SF"},
    {"2022Re-recoBCD", "Electron-ID-SF"},
    {"2022Re-recoE+PromptFG", "Electron-ID-SF"},
    // {"2024Prompt", "Electron-ID-SF"}
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
std::unordered_map<std::string, correction::CorrectionSet> PNETWScaleFactors = {
    {"2022", *CorrectionSet::from_file("corrections/scalefactors/particlenet/pnet.json")},
    {"2022EE", *CorrectionSet::from_file("corrections/scalefactors/particlenet/pnet.json")}
};
RNode applyPNETWScaleFactors(std::unordered_map<std::string, correction::CorrectionSet> cset_pnet_w, RNode df);

/*
############################################
B-TAGGING SFs
############################################
*/
const std::unordered_map<std::string, correction::CorrectionSet> bTaggingScaleFactors = {
    {"2016preVFP", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/2016preVFP_UL/btagging.json.gz")},
    {"2016postVFP", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/2016postVFP_UL/btagging.json.gz")},
    {"2017", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/2017_UL/btagging.json.gz")},
    {"2018", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/2018_UL/btagging.json.gz")},
    {"2022Re-recoBCD", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/2022_Summer22/btagging.json.gz")},
    {"2022Re-recoE+PromptFG", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/2022_Summer22EE/btagging.json.gz")},
    {"2024Prompt", *CorrectionSet::from_file("/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/2024_Summer24/btagging.json.gz")},
    {"eff", *CorrectionSet::from_file("corrections/scalefactors/btagging/btag_eff.json")}
};

const std::unordered_map<std::string, std::string> bTaggingScaleFactors_HF_corrname = {
    {"2016preVFP", "deepCSV_comb"},
    {"2016postVFP", "deepCSV_comb"},
    {"2017", "deepCSV_comb"},
    {"2018", "deepCSV_comb"}
};

const std::unordered_map<std::string, std::string> bTaggingScaleFactors_LF_corrname = {
    {"2016preVFP", "deepCSV_incl"},
    {"2016postVFP", "deepCSV_incl"},
    {"2017", "deepCSV_incl"},
    {"2018", "deepCSV_incl"}
};

RNode applyBTaggingScaleFactors(std::unordered_map<std::string, correction::CorrectionSet> cset_btag, std::unordered_map<std::string, std::string> corrname_map_HF, std::unordered_map<std::string, std::string> corrname_map_LF, RNode df);

/*
############################################
OTHER SFs
############################################
*/

RNode applyEWKCorrections(CorrectionSet cset_ewk, RNode df);
const auto cset_ewk = *CorrectionSet::from_file("corrections/scalefactors/ewk/EWK.json");

RNode applyL1PreFiringReweighting(RNode df);
RNode applyPSWeight_FSR(RNode df);
RNode applyPSWeight_ISR(RNode df);
RNode applyLHEScaleWeight_muF(RNode df);
RNode applyLHEScaleWeight_muR(RNode df);
RNode applyLHEWeights_pdf(RNode df);

RNode applyDataWeights(RNode df);
RNode applyMCWeights(RNode df);

} // end namespace Run3

#endif //WEIGHTS_RUN3_H