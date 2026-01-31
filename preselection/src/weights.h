#ifndef WEIGHTS
#define WEIGHTS

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

/*
############################################
GOLDEN JSON
############################################
*/
RNode applyGoldenJSONWeight(const lumiMask& golden, RNode df);
const auto LumiMask = lumiMask::fromJSON({
    "corrections/goldenJson/Cert_Collisions2022_355100_362760_Golden.json", 
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
    {"2016preVFP", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/LUM/Run2-2016preVFP-UL-NanoAODv9/latest/puWeights.json.gz")},
    {"2016postVFP", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/LUM/Run2-2016postVFP-UL-NanoAODv9/latest/puWeights.json.gz")},
    {"2017", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/LUM/Run2-2017-UL-NanoAODv9/latest/puWeights.json.gz")},
    {"2018", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/LUM/Run2-2018-UL-NanoAODv9/latest/puWeights.json.gz")},
    {"2022Re-recoBCD", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/LUM/Run3-22CDSep23-Summer22-NanoAODv12/latest/puWeights.json.gz")},
    {"2022Re-recoE+PromptFG", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/LUM/Run3-22EFGSep23-Summer22EE-NanoAODv12/latest/puWeights.json.gz")},
    {"2023PromptC", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/LUM/Run3-23CSep23-Summer23-NanoAODv12/latest/puWeights.json.gz")},
    {"2023PromptD", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/LUM/Run3-23DSep23-Summer23BPix-NanoAODv12/latest/puWeights.json.gz")},
    {"2024Prompt", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/LUM/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/latest/puWeights_BCDEFGHI.json.gz")}
};
const std::unordered_map<std::string, std::string> pileupScaleFactors_yearmap = {
    {"2016preVFP", "Collisions16_UltraLegacy_goldenJSON"},
    {"2016postVFP", "Collisions16_UltraLegacy_goldenJSON"},
    {"2017", "Collisions17_UltraLegacy_goldenJSON"},
    {"2018", "Collisions18_UltraLegacy_goldenJSON"},
    {"2022Re-recoBCD", "Collisions2022_355100_357900_eraBCD_GoldenJson"},
    {"2022Re-recoE+PromptFG", "Collisions2022_359022_362760_eraEFG_GoldenJson"},
    {"2023PromptC", "Collisions2023_366403_369802_eraBC_GoldenJson"},
    {"2023PromptD", "Collisions2023_369803_370790_eraD_GoldenJson"},
    {"2024Prompt", "Collisions24_BCDEFGHI_goldenJSON"}
};
RNode applyPileupScaleFactors(std::unordered_map<std::string, correction::CorrectionSet> cset_pileup, std::unordered_map<std::string, std::string> year_map, RNode df);

// REMOVE PILEUP ID SINCE NO LONGER NEEDED FOR PUPPI

/*
############################################
MUON SFs
############################################
*/
const std::unordered_map<std::string, correction::CorrectionSet> muonScaleFactors = {
    {"2016preVFP", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/MUO/Run2-2016preVFP-UL-NanoAODv9/latest/muon_Z.json.gz")},
    {"2016postVFP", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/MUO/Run2-2016postVFP-UL-NanoAODv9/latest/muon_Z.json.gz")},
    {"2017", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/MUO/Run2-2017-UL-NanoAODv9/latest/muon_Z.json.gz")},
    {"2018", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/MUO/Run2-2018-UL-NanoAODv9/latest/muon_Z.json.gz")},
    {"2022Re-recoBCD", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/MUO/Run3-22CDSep23-Summer22-NanoAODv12/latest/muon_Z.json.gz")},
    {"2022Re-recoE+PromptFG", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/MUO/Run3-22EFGSep23-Summer22EE-NanoAODv12/latest/muon_Z.json.gz")},
    {"2023PromptC", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/MUO/Run3-23CSep23-Summer23-NanoAODv12/latest/muon_Z.json.gz")},
    {"2023PromptD", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/MUO/Run3-23DSep23-Summer23BPix-NanoAODv12/latest/muon_Z.json.gz")},
    {"2024Prompt", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/MUO/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/latest/muon_Z.json.gz")}
};

const std::unordered_map<std::string, std::string> muonIDScaleFactors_yearmap = {
    {"2016preVFP", "NUM_TightID_DEN_TrackerMuons"},
    {"2016postVFP", "NUM_TightID_DEN_TrackerMuons"},
    {"2017", "NUM_TightID_DEN_TrackerMuons"},
    {"2018", "NUM_TightID_DEN_TrackerMuons"},
    {"2022Re-recoBCD", "NUM_TightID_DEN_TrackerMuons"},
    {"2022Re-recoE+PromptFG", "NUM_TightID_DEN_TrackerMuons"},
    {"2023PromptC", "NUM_TightID_DEN_TrackerMuons"},
    {"2023PromptD", "NUM_TightID_DEN_TrackerMuons"},
    {"2024Prompt", "NUM_TightID_DEN_TrackerMuons"}
};
RNode applyMuonIDScaleFactors(std::unordered_map<std::string, correction::CorrectionSet> cset_muon, std::unordered_map<std::string, std::string> year_map, RNode df);

const std::unordered_map<std::string, std::string> muonRecoScaleFactors_yearmap = {
    {"2016preVFP", "NUM_TightPFIso_DEN_TightID"},
    {"2016postVFP", "NUM_TightPFIso_DEN_TightID"},
    {"2017", "NUM_TightPFIso_DEN_TightID"},
    {"2018", "NUM_TightPFIso_DEN_TightID"},
    {"2022Re-recoBCD", "NUM_TightPFIso_DEN_TightID"},
    {"2022Re-recoE+PromptFG", "NUM_TightPFIso_DEN_TightID"},
    {"2023PromptC", "NUM_TightPFIso_DEN_TightID"},
    {"2023PromptD", "NUM_TightPFIso_DEN_TightID"},
    {"2024Prompt", "NUM_TightPFIso_DEN_TightID"}
};
RNode applyMuonRecoScaleFactors(std::unordered_map<std::string, correction::CorrectionSet> cset_muon, std::unordered_map<std::string, std::string> year_map, RNode df);

const std::unordered_map<std::string, std::string> muonTriggerScaleFactors_yearmap = {
    {"2016preVFP", "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight"},
    {"2016postVFP", "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight"},
    {"2017", "NUM_IsoMu27_DEN_CutBasedIdTight_and_PFIsoTight"},
    {"2018", "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight"},
    {"2022Re-recoBCD", "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight"},
    {"2022Re-recoE+PromptFG", "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight"},
    {"2023PromptC", "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight"},
    {"2023PromptD", "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight"},
    {"2024Prompt", "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight"}
};
RNode applyMuonTriggerScaleFactors(std::unordered_map<std::string, correction::CorrectionSet> cset_muon, std::unordered_map<std::string, std::string> year_map, RNode df);

/*
############################################
ELECTRON SFs
############################################
*/
const std::unordered_map<std::string, correction::CorrectionSet> electronScaleFactors = {
    {"2016preVFP", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/Run2-2016preVFP-UL-NanoAODv9/latest/electron.json.gz")},
    {"2016postVFP", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/Run2-2016postVFP-UL-NanoAODv9/latest/electron.json.gz")},
    {"2017", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/Run2-2017-UL-NanoAODv9/latest/electron.json.gz")},
    {"2018", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/Run2-2018-UL-NanoAODv9/latest/electron.json.gz")},
    {"2022Re-recoBCD", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/Run3-22CDSep23-Summer22-NanoAODv12/latest/electron.json.gz")},
    {"2022Re-recoE+PromptFG", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/Run3-22EFGSep23-Summer22EE-NanoAODv12/latest/electron.json.gz")},
    {"2023PromptC", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/Run3-23CSep23-Summer23-NanoAODv12/latest/electron.json.gz")},
    {"2023PromptD", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/Run3-23DSep23-Summer23BPix-NanoAODv12/latest/electron.json.gz")}
    // {"2024Prompt", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/latest/electron.json.gz")}
};
const std::unordered_map<std::string, std::string> electronScaleFactors_yearmap = {
    {"2016preVFP", "UL-Electron-ID-SF"},
    {"2016postVFP", "UL-Electron-ID-SF"},
    {"2017", "UL-Electron-ID-SF"},
    {"2018", "UL-Electron-ID-SF"},
    {"2022Re-recoBCD", "Electron-ID-SF"},
    {"2022Re-recoE+PromptFG", "Electron-ID-SF"},
    {"2023PromptC", "Electron-ID-SF"},
    {"2023PromptD", "Electron-ID-SF"}
    // {"2024Prompt", "Electron-ID-SF"}
};
RNode applyElectronIDScaleFactors(std::unordered_map<std::string, correction::CorrectionSet> cset_electron, std::unordered_map<std::string, std::string> year_map, RNode df);
RNode applyElectronRecoScaleFactors(std::unordered_map<std::string, correction::CorrectionSet> cset_electron, std::unordered_map<std::string, std::string> year_map, RNode df);

const std::unordered_map<std::string, correction::CorrectionSet> electronTriggerScaleFactors = {
    {"2022Re-recoBCD", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/Run3-22CDSep23-Summer22-NanoAODv12/latest/electronHlt.json.gz")},
    {"2022Re-recoE+PromptFG", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/Run3-22EFGSep23-Summer22EE-NanoAODv12/latest/electronHlt.json.gz")},
    {"2023PromptC", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/Run3-23CSep23-Summer23-NanoAODv12/latest/electronHlt.json.gz")},
    {"2023PromptD", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/Run3-23DSep23-Summer23BPix-NanoAODv12/latest/electronHlt.json.gz")}
};
const std::unordered_map<std::string, std::string> electronTriggerScaleFactors_yearmap = {
    {"2022Re-recoBCD", "Electron-HLT-SF"},
    {"2022Re-recoE+PromptFG", "Electron-HLT-SF"},
    {"2023PromptC", "Electron-HLT-SF"},
    {"2023PromptD", "Electron-HLT-SF"}
};
RNode applyElectronTriggerScaleFactors(std::unordered_map<std::string, correction::CorrectionSet> cset_electron, std::unordered_map<std::string, std::string> year_map, RNode df);

/*
############################################
B-TAGGING SFs
############################################
*/
const std::unordered_map<std::string, correction::CorrectionSet> bTaggingScaleFactors = {
    {"2016preVFP", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/BTV/Run2-2016preVFP-UL-NanoAODv9/latest/btagging.json.gz")},
    {"2016postVFP", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/BTV/Run2-2016postVFP-UL-NanoAODv9/latest/btagging.json.gz")},
    {"2017", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/BTV/Run2-2017-UL-NanoAODv9/latest/btagging.json.gz")},
    {"2018", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/BTV/Run2-2018-UL-NanoAODv9/latest/btagging.json.gz")},
    {"2022Re-recoBCD", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/BTV/Run3-22CDSep23-Summer22-NanoAODv12/latest/btagging.json.gz")},
    {"2022Re-recoE+PromptFG", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/BTV/Run3-22EFGSep23-Summer22EE-NanoAODv12/latest/btagging.json.gz")},
    {"2023PromptC", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/BTV/Run3-23CSep23-Summer23-NanoAODv12/latest/btagging.json.gz")},
    {"2023PromptD", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/BTV/Run3-23DSep23-Summer23BPix-NanoAODv12/latest/btagging.json.gz")},
    {"2024Prompt", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/BTV/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/latest/btagging.json.gz")},
    {"eff", *CorrectionSet::from_file("corrections/scalefactors/btagging/btag_eff.json")}
};

// THIS NEEDS TO BE UPDATED FOR PARTICLE TRANSFORMER

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

#endif //WEIGHTS