#ifndef CORRECTIONS_H
#define CORRECTIONS_H

#pragma once

#include <unordered_map>
#include <string>
#include <iostream>
#include <vector>
#include <array>

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RVec.hxx"

#include "correction.h"
#include "TRandom3.h"

using correction::CorrectionSet;
using RNode = ROOT::RDF::RNode;
using ROOT::VecOps::RVec;

/*
############################################
B-TAGGING WORKING POINTS
############################################
*/
RVec<bool> isbTagLoose(std::string year, RVec<float> btag_score);
RVec<bool> isbTagMedium(std::string year, RVec<float> btag_score);
RVec<bool> isbTagTight(std::string year, RVec<float> btag_score);


// Note: Re-using 2024 WPs for 2022-2023
const std::unordered_map <std::string, correction::CorrectionSet> btaggingCorrections = {
    {"2016preVFP", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/BTV/Run2-2016preVFP-UL-NanoAODv15/latest/btagging.json.gz")},
    {"2016postVFP", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/BTV/Run2-2016postVFP-UL-NanoAODv15/latest/btagging.json.gz")},
    {"2017", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/BTV/Run2-2017-UL-NanoAODv15/latest/btagging.json.gz")},
    {"2018", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/BTV/Run2-2018-UL-NanoAODv15/latest/btagging.json.gz")},
    {"2022Re-recoBCD", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/BTV/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/latest/btagging.json.gz")},
    {"2022Re-recoE+PromptFG", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/BTV/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/latest/btagging.json.gz")},
    {"2023PromptC", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/BTV/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/latest/btagging.json.gz")},
    {"2023PromptD", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/BTV/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/latest/btagging.json.gz")},
    {"2024Prompt", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/BTV/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/latest/btagging.json.gz")},
    // Recomendation from: https://btv-wiki.docs.cern.ch/ScaleFactors/Run3Prompt25/
    {"2025", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/BTV/Run3-25Prompt-Summer24-NanoAODv15/latest/btagging.json.gz")}
};

// 2. Map of Numbers: Put the actual numeric thresholds here
static std::unordered_map<std::string, float> btaggingWPMap_Loose = {
    {"2016preVFP",  btaggingCorrections.at("2016preVFP").at("UParTAK4_wp_values")->evaluate({"L"})}, 
    {"2016postVFP", btaggingCorrections.at("2016postVFP").at("UParTAK4_wp_values")->evaluate({"L"})},
    {"2017",        btaggingCorrections.at("2017").at("UParTAK4_wp_values")->evaluate({"L"})},
    {"2018",        btaggingCorrections.at("2018").at("UParTAK4_wp_values")->evaluate({"L"})},
    {"2022Re-recoBCD", btaggingCorrections.at("2022Re-recoBCD").at("UParTAK4_wp_values")->evaluate({"L"})},
    {"2022Re-recoE+PromptFG", btaggingCorrections.at("2022Re-recoE+PromptFG").at("UParTAK4_wp_values")->evaluate({"L"})},
    {"2023PromptC", btaggingCorrections.at("2023PromptC").at("UParTAK4_wp_values")->evaluate({"L"})},
    {"2023PromptD", btaggingCorrections.at("2023PromptD").at("UParTAK4_wp_values")->evaluate({"L"})},
    {"2024Prompt",  btaggingCorrections.at("2024Prompt").at("UParTAK4_wp_values")->evaluate({"L"})},
    {"2025",        btaggingCorrections.at("2025").at("UParTAK4_wp_values")->evaluate({"L"})}
};

static std::unordered_map<std::string, float> btaggingWPMap_Medium = {
    {"2016preVFP",  btaggingCorrections.at("2016preVFP").at("UParTAK4_wp_values")->evaluate({"M"})},
    {"2016postVFP", btaggingCorrections.at("2016postVFP").at("UParTAK4_wp_values")->evaluate({"M"})},
    {"2017",        btaggingCorrections.at("2017").at("UParTAK4_wp_values")->evaluate({"M"})},
    {"2018",        btaggingCorrections.at("2018").at("UParTAK4_wp_values")->evaluate({"M"})},
    {"2022Re-recoBCD", btaggingCorrections.at("2022Re-recoBCD").at("UParTAK4_wp_values")->evaluate({"M"})},
    {"2022Re-recoE+PromptFG", btaggingCorrections.at("2022Re-recoE+PromptFG").at("UParTAK4_wp_values")->evaluate({"M"})},
    {"2023PromptC", btaggingCorrections.at("2023PromptC").at("UParTAK4_wp_values")->evaluate({"M"})},
    {"2023PromptD", btaggingCorrections.at("2023PromptD").at("UParTAK4_wp_values")->evaluate({"M"})},
    {"2024Prompt",  btaggingCorrections.at("2024Prompt").at("UParTAK4_wp_values")->evaluate({"M"})},
    {"2025",        btaggingCorrections.at("2025").at("UParTAK4_wp_values")->evaluate({"M"})}
};

static std::unordered_map<std::string, float> btaggingWPMap_Tight = {
    {"2016preVFP",  btaggingCorrections.at("2016preVFP").at("UParTAK4_wp_values")->evaluate({"T"})},
    {"2016postVFP", btaggingCorrections.at("2016postVFP").at("UParTAK4_wp_values")->evaluate({"T"})},
    {"2017",        btaggingCorrections.at("2017").at("UParTAK4_wp_values")->evaluate({"T"})},
    {"2018",        btaggingCorrections.at("2018").at("UParTAK4_wp_values")->evaluate({"T"})},
    {"2022Re-recoBCD", btaggingCorrections.at("2022Re-recoBCD").at("UParTAK4_wp_values")->evaluate({"T"})},
    {"2022Re-recoE+PromptFG", btaggingCorrections.at("2022Re-recoE+PromptFG").at("UParTAK4_wp_values")->evaluate({"T"})},
    {"2023PromptC", btaggingCorrections.at("2023PromptC").at("UParTAK4_wp_values")->evaluate({"T"})},
    {"2023PromptD", btaggingCorrections.at("2023PromptD").at("UParTAK4_wp_values")->evaluate({"T"})},
    {"2024Prompt",  btaggingCorrections.at("2024Prompt").at("UParTAK4_wp_values")->evaluate({"T"})},
    {"2025",        btaggingCorrections.at("2025").at("UParTAK4_wp_values")->evaluate({"T"})}
};

/*
############################################
MET CORRECTIONS
############################################
*/
// FIXME: met corrections missing for v15
const std::unordered_map<std::string, correction::CorrectionSet> metCorrections = {
    {"2016preVFP", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run2-2016preVFP-UL-NanoAODv9/latest/met.json.gz")},
    {"2016postVFP", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run2-2016postVFP-UL-NanoAODv9/latest/met.json.gz")},
    {"2017", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run2-2017-UL-NanoAODv9/latest/met.json.gz")},
    {"2018", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run2-2018-UL-NanoAODv9/latest/met.json.gz")},
    // {"2022Re-recoBCD", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-22CDSep23-Summer22-NanoAODv12/latest/met_xyCorrections_2022_2022.json.gz")},
    // {"2022Re-recoE+PromptFG", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-22EFGSep23-Summer22EE-NanoAODv12/latest/met_xyCorrections_2022_2022EE.json.gz")},
    // {"2023PromptC", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-23CSep23-Summer23-NanoAODv12/latest/met_xyCorrections_2023_2023.json.gz")},
    // {"2023PromptD", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-23DSep23-Summer23BPix-NanoAODv12/latest/met_xyCorrections_2023_2023BPix.json.gz")}
};
RNode applyMETPhiCorrections(RNode df, bool isData);
RNode applyMETUnclusteredCorrections(RNode df, std::string variation);

/*
############################################
JET MASS SCALE (JMS) AND RESOLUTION (JMR)
############################################

Convention:
  - JMS  : *additive* shift on FatJet_msoftdrop, in GeV. ±1σ variation = jmsr_scale.
           Central = 0.0 GeV (no shift).
  - JMR  : *multiplicative* width factor wrt the gen mass. ±1σ variation =
           jmsr_smear. Central = 1.0 (no widening).

We currently store only central values (= identity) since no calibration has been
derived yet. Replace these with per-era values from the calibration fit when
available; the up/down systematic variants can then be added by passing
variation = "up"/"down" with the derived ±1σ shift/scale.
*/

const std::unordered_map<std::string, float> jetMassScale_central = {
    {"2016preVFP", 0.0f}, {"2016postVFP", 0.0f}, {"2017", 0.0f}, {"2018", 0.0f},
    {"2022Re-recoBCD", 0.0f}, {"2022Re-recoE+PromptFG", 0.0f},
    {"2023PromptC", 0.0f}, {"2023PromptD", 0.0f},
    {"2024Prompt", 0.0f}
};

const std::unordered_map<std::string, float> jetMassResolution_central = {
    {"2016preVFP", 1.0f}, {"2016postVFP", 1.0f}, {"2017", 1.0f}, {"2018", 1.0f},
    {"2022Re-recoBCD", 1.0f}, {"2022Re-recoE+PromptFG", 1.0f},
    {"2023PromptC", 1.0f}, {"2023PromptD", 1.0f},
    {"2024Prompt", 1.0f}
};

// Relative msoftdrop resolution sigma(msd)/msd. Used by the unmatched stochastic branch
// of applyJetMassResolution. Placeholder 1.0 — replace with the per-era value from the same
// calibration fit that derives jetMassResolution_central. With factor f = 1.0 (placeholder)
// the branch is unreachable, so this value has no effect until both maps are updated together.
const std::unordered_map<std::string, float> jetMassResolution_sigmaRel_central = {
    {"2016preVFP", 1.0f}, {"2016postVFP", 1.0f}, {"2017", 1.0f}, {"2018", 1.0f},
    {"2022Re-recoBCD", 1.0f}, {"2022Re-recoE+PromptFG", 1.0f},
    {"2023PromptC", 1.0f}, {"2023PromptD", 1.0f},
    {"2024Prompt", 1.0f}
};

RNode applyJetMassScale(const std::unordered_map<std::string, float>& jms_shift, RNode df);
RNode applyJetMassResolution(const std::unordered_map<std::string, float>& jmr_factor,
                             const std::unordered_map<std::string, float>& jmr_sigma_rel,
                             RNode df);

/*
############################################
JET ENERGY CORRECTIONS
############################################
*/
const std::unordered_map<std::string, correction::CorrectionSet> jetEnergyCorrections = {
    {"2016preVFP", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run2-2016preVFP-UL-NanoAODv15/latest/jet_jerc.json.gz")},
    {"2016postVFP", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run2-2016postVFP-UL-NanoAODv15/latest/jet_jerc.json.gz")},
    {"2017", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run2-2017-UL-NanoAODv15/latest/jet_jerc.json.gz")},
    {"2018", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run2-2018-UL-NanoAODv15/latest/jet_jerc.json.gz")},
    {"2022Re-recoBCD", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-22CDSep23-Summer22-NanoAODv12/latest/jet_jerc.json.gz")},
    {"2022Re-recoE+PromptFG", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-22EFGSep23-Summer22EE-NanoAODv12/latest/jet_jerc.json.gz")},
    {"2023PromptC", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-23CSep23-Summer23-NanoAODv12/latest/jet_jerc.json.gz")},
    {"2023PromptD", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-23DSep23-Summer23BPix-NanoAODv12/latest/jet_jerc.json.gz")},
    {"2024Prompt", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/latest/jet_jerc.json.gz")}
};

// AK8 fat-jet JEC/JER lives in fatJet_jerc.json.gz under the same era directories.
const std::unordered_map<std::string, correction::CorrectionSet> fatJetEnergyCorrections = {
    {"2016preVFP", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run2-2016preVFP-UL-NanoAODv15/latest/fatJet_jerc.json.gz")},
    {"2016postVFP", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run2-2016postVFP-UL-NanoAODv15/latest/fatJet_jerc.json.gz")},
    {"2017", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run2-2017-UL-NanoAODv15/latest/fatJet_jerc.json.gz")},
    {"2018", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run2-2018-UL-NanoAODv15/latest/fatJet_jerc.json.gz")},
    {"2022Re-recoBCD", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-22CDSep23-Summer22-NanoAODv12/latest/fatJet_jerc.json.gz")},
    {"2022Re-recoE+PromptFG", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-22EFGSep23-Summer22EE-NanoAODv12/latest/fatJet_jerc.json.gz")},
    {"2023PromptC", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-23CSep23-Summer23-NanoAODv12/latest/fatJet_jerc.json.gz")},
    {"2023PromptD", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-23DSep23-Summer23BPix-NanoAODv12/latest/fatJet_jerc.json.gz")},
    {"2024Prompt", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/latest/fatJet_jerc.json.gz")}
};

const std::unordered_map<std::string, correction::CorrectionSet> jetEnergyResolution_smear = {
    {"jer_smear", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/JER-Smearing/latest/jer_smear.json.gz")}
};

// JEC tag base — `<TAG>_MC` and `<TAG>_DATA` compounds (`_L1L2L3Res_<algo>`) live in jet_jerc.json.gz.
// The trailing `_MC` here is stripped and replaced with `_DATA` at runtime when running on data.
const std::unordered_map<std::string, std::string> jetEnergyCorrections_JEC_prefix = {
    {"2016preVFP", "Summer20UL16APVNanoV15_V1_MC"},
    {"2016postVFP", "Summer20UL16NanoV15_V1_MC"},
    {"2017", "Summer20UL17NanoV15_V1_MC"},
    {"2018", "Summer20UL18NanoV15_V1_MC"},
    {"2022Re-recoBCD", "Summer22_22Sep2023_V3_MC"},
    {"2022Re-recoE+PromptFG", "Summer22EE_22Sep2023_V3_MC"},
    {"2023PromptC", "Summer23Prompt23_V3_MC"},
    {"2023PromptD", "Summer23BPixPrompt23_V3_MC"},
    {"2024Prompt", "Summer24Prompt24_V2_MC"}
};

const std::unordered_map<std::string, std::string> jetEnergyCorrections_JEC_suffix = {
    {"2016preVFP", "AK4PFPuppi"},
    {"2016postVFP", "AK4PFPuppi"},
    {"2017", "AK4PFPuppi"},
    {"2018", "AK4PFPuppi"},
    {"2022Re-recoBCD", "AK4PFPuppi"},
    {"2022Re-recoE+PromptFG", "AK4PFPuppi"},
    {"2023PromptC", "AK4PFPuppi"},
    {"2023PromptD", "AK4PFPuppi"},
    {"2024Prompt", "AK4PFPuppi"}
};

// AK8 fat-jet variants — same year keys, AK8PFPuppi algo across all eras.
const std::unordered_map<std::string, std::string> fatJetEnergyCorrections_JEC_prefix = {
    {"2016preVFP", "Summer20UL16APVNanoV15_V1_MC"},
    {"2016postVFP", "Summer20UL16NanoV15_V1_MC"},
    {"2017", "Summer20UL17NanoV15_V1_MC"},
    {"2018", "Summer20UL18NanoV15_V1_MC"},
    {"2022Re-recoBCD", "Summer22_22Sep2023_V3_MC"},
    {"2022Re-recoE+PromptFG", "Summer22EE_22Sep2023_V3_MC"},
    {"2023PromptC", "Summer23Prompt23_V3_MC"},
    {"2023PromptD", "Summer23BPixPrompt23_V3_MC"},
    {"2024Prompt", "Summer24Prompt24_V2_MC"}
};

const std::unordered_map<std::string, std::string> fatJetEnergyCorrections_JEC_suffix = {
    {"2016preVFP", "AK8PFPuppi"},
    {"2016postVFP", "AK8PFPuppi"},
    {"2017", "AK8PFPuppi"},
    {"2018", "AK8PFPuppi"},
    {"2022Re-recoBCD", "AK8PFPuppi"},
    {"2022Re-recoE+PromptFG", "AK8PFPuppi"},
    {"2023PromptC", "AK8PFPuppi"},
    {"2023PromptD", "AK8PFPuppi"},
    {"2024Prompt", "AK8PFPuppi"}
};

// Year token embedded in the year-decorrelated Regrouped JEC source names
// (e.g. "Summer20UL18NanoV15_V1_MC_Regrouped_Absolute_2018_AK4PFPuppi"). Confirmed
// against jet_jerc.json.gz / fatJet_jerc.json.gz for every campaign on 2026-05-14.
const std::unordered_map<std::string, std::string> jetEnergyCorrections_yearToken = {
    {"2016preVFP", "2016APV"},
    {"2016postVFP", "2016"},
    {"2017", "2017"},
    {"2018", "2018"},
    {"2022Re-recoBCD", "2022"},
    {"2022Re-recoE+PromptFG", "2022EE"},
    {"2023PromptC", "2023"},
    {"2023PromptD", "2023BPix"},
    {"2024Prompt", "2024"}
};

// JER tag names 
// FIXME: 2024 re-uses the 2023BPix JR until a dedicated 2024 release lands
const std::unordered_map<std::string, std::string> jetEnergyResolution_JER_res_name = {
    {"2016preVFP", "Summer20UL16APV_JRV3_MC_PtResolution_AK4PFPuppi"},
    {"2016postVFP", "Summer20UL16_JRV3_MC_PtResolution_AK4PFPuppi"},
    {"2017", "Summer19UL17_JRV3_MC_PtResolution_AK4PFPuppi"},
    {"2018", "Summer19UL18_JRV2_MC_PtResolution_AK4PFPuppi"},
    {"2022Re-recoBCD", "Summer22_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi"},
    {"2022Re-recoE+PromptFG", "Summer22EE_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi"},
    {"2023PromptC", "Summer23Prompt23_RunCv1234_JRV1_MC_PtResolution_AK4PFPuppi"},
    {"2023PromptD", "Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi"},
    {"2024Prompt", "Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi"}
};

const std::unordered_map<std::string, std::string> jetEnergyResolution_JER_sf_name = {
    {"2016preVFP", "Summer20UL16APV_JRV3_MC_ScaleFactor_AK4PFPuppi"},
    {"2016postVFP", "Summer20UL16_JRV3_MC_ScaleFactor_AK4PFPuppi"},
    {"2017", "Summer19UL17_JRV3_MC_ScaleFactor_AK4PFPuppi"},
    {"2018", "Summer19UL18_JRV2_MC_ScaleFactor_AK4PFPuppi"},
    {"2022Re-recoBCD", "Summer22_22Sep2023_JRV1_MC_ScaleFactor_AK4PFPuppi"},
    {"2022Re-recoE+PromptFG", "Summer22EE_22Sep2023_JRV1_MC_ScaleFactor_AK4PFPuppi"},
    {"2023PromptC", "Summer23Prompt23_RunCv1234_JRV1_MC_ScaleFactor_AK4PFPuppi"},
    {"2023PromptD", "Summer23BPixPrompt23_RunD_JRV1_MC_ScaleFactor_AK4PFPuppi"},
    {"2024Prompt", "Summer23BPixPrompt23_RunD_JRV1_MC_ScaleFactor_AK4PFPuppi"}
};

// AK8 JER tags — same JR releases as AK4 but with AK8PFPuppi algo.
const std::unordered_map<std::string, std::string> fatJetEnergyResolution_JER_res_name = {
    {"2016preVFP", "Summer20UL16APV_JRV3_MC_PtResolution_AK8PFPuppi"},
    {"2016postVFP", "Summer20UL16_JRV3_MC_PtResolution_AK8PFPuppi"},
    {"2017", "Summer19UL17_JRV3_MC_PtResolution_AK8PFPuppi"},
    {"2018", "Summer19UL18_JRV2_MC_PtResolution_AK8PFPuppi"},
    {"2022Re-recoBCD", "Summer22_22Sep2023_JRV1_MC_PtResolution_AK8PFPuppi"},
    {"2022Re-recoE+PromptFG", "Summer22EE_22Sep2023_JRV1_MC_PtResolution_AK8PFPuppi"},
    {"2023PromptC", "Summer23Prompt23_RunCv1234_JRV1_MC_PtResolution_AK8PFPuppi"},
    {"2023PromptD", "Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK8PFPuppi"},
    {"2024Prompt", "Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK8PFPuppi"}
};

const std::unordered_map<std::string, std::string> fatJetEnergyResolution_JER_sf_name = {
    {"2016preVFP", "Summer20UL16APV_JRV3_MC_ScaleFactor_AK8PFPuppi"},
    {"2016postVFP", "Summer20UL16_JRV3_MC_ScaleFactor_AK8PFPuppi"},
    {"2017", "Summer19UL17_JRV3_MC_ScaleFactor_AK8PFPuppi"},
    {"2018", "Summer19UL18_JRV2_MC_ScaleFactor_AK8PFPuppi"},
    {"2022Re-recoBCD", "Summer22_22Sep2023_JRV1_MC_ScaleFactor_AK8PFPuppi"},
    {"2022Re-recoE+PromptFG", "Summer22EE_22Sep2023_JRV1_MC_ScaleFactor_AK8PFPuppi"},
    {"2023PromptC", "Summer23Prompt23_RunCv1234_JRV1_MC_ScaleFactor_AK8PFPuppi"},
    {"2023PromptD", "Summer23BPixPrompt23_RunD_JRV1_MC_ScaleFactor_AK8PFPuppi"},
    {"2024Prompt", "Summer23BPixPrompt23_RunD_JRV1_MC_ScaleFactor_AK8PFPuppi"}
};

// Nominal JEC: removes the JEC stored in NanoAOD (via Jet_rawFactor) and re-applies the
// latest L1FastJet * L2Relative * L3Absolute (* L2L3Residual for data) compound from
// jet_jerc.json.gz. Propagates the change to met_pt / met_phi (Type-I).
RNode applyJetEnergyCorrections(const std::unordered_map<std::string, correction::CorrectionSet>& cset_jerc,
                                const std::unordered_map<std::string, std::string>& jec_prefix_map,
                                const std::unordered_map<std::string, std::string>& jec_suffix_map,
                                RNode df, bool isData);

// JES uncertainty — per-source, suffixed-branch output. Must run after the nominal
// applyJetEnergyCorrections / applyFatJetEnergyCorrections so the shift is applied on top
// of the corrected baseline. Each call writes one set of suffixed columns; the driver
// applyJESVariations loops over the 11 Regrouped V2 sources × {Up, Dn}.
//
// AK4 helper writes Jet_pt_<suffix>, Jet_mass_<suffix>, met_pt_<suffix>, met_phi_<suffix>
// (Type-I MET propagated). AK8 helper writes FatJet_pt_<suffix>, FatJet_mass_<suffix>
// (NOT FatJet_msoftdrop — that has its own JEC recipe) and does NOT propagate to MET
// because Type-I MET is built from AK4 only.
//
// <suffix> = "jes" + column_label + direction, e.g. "jesAbsoluteUp", "jesAbsoluteYearDn".
RNode defineAK4JESVariation(const std::unordered_map<std::string, correction::CorrectionSet>& cset_jerc,
                            const std::unordered_map<std::string, std::string>& jec_prefix_map,
                            const std::unordered_map<std::string, std::string>& jec_suffix_map,
                            const std::unordered_map<std::string, std::string>& year_token_map,
                            RNode df, const std::string& column_label, const std::string& base_source,
                            bool yearDecorrelated, const std::string& direction);

RNode defineFatJetJESVariation(const std::unordered_map<std::string, correction::CorrectionSet>& cset_jerc,
                               const std::unordered_map<std::string, std::string>& jec_prefix_map,
                               const std::unordered_map<std::string, std::string>& jec_suffix_map,
                               const std::unordered_map<std::string, std::string>& year_token_map,
                               RNode df, const std::string& column_label, const std::string& base_source,
                               bool yearDecorrelated, const std::string& direction);

// Driver: emits all 11 Regrouped V2 sources × {Up, Dn} = 22 variations on AK4 + AK8.
RNode applyJESVariations(RNode df);

// When called with false before any analysis, disables all JES/JER variation branches
// (jesVariationSuffixes() returns empty, applyJESVariations becomes a no-op). Use
// --no_systs to activate nominal-only mode.
void setStoreSysts(bool v);

// Public accessor for the 22 variation suffixes (e.g. "jesAbsoluteUp",
// "jesRelativeSampleYearDn"). Returns empty when setStoreSysts(false) has been called.
std::vector<std::string> jesVariationSuffixes();

// JER smearing (MC only). Propagates to met_pt / met_phi as well.
RNode applyJetEnergyResolution(const std::unordered_map<std::string, correction::CorrectionSet>& cset_jerc,
                               const std::unordered_map<std::string, correction::CorrectionSet>& cset_jer_smear,
                               const std::unordered_map<std::string, std::string>& jer_res_map,
                               const std::unordered_map<std::string, std::string>& jer_sf_map,
                               RNode df, std::string variation);

// AK8 (FatJet_*) variants — same recipe as AK4 but reading FatJet_* / GenJetAK8_* branches
// and the AK8PFPuppi compound from fatJet_jerc.json.gz. 
RNode applyFatJetEnergyCorrections(const std::unordered_map<std::string, correction::CorrectionSet>& cset_jerc,
                                   const std::unordered_map<std::string, std::string>& jec_prefix_map,
                                   const std::unordered_map<std::string, std::string>& jec_suffix_map,
                                   RNode df, bool isData);

RNode applyFatJetEnergyResolution(const std::unordered_map<std::string, correction::CorrectionSet>& cset_jerc,
                                  const std::unordered_map<std::string, correction::CorrectionSet>& cset_jer_smear,
                                  const std::unordered_map<std::string, std::string>& jer_res_map,
                                  const std::unordered_map<std::string, std::string>& jer_sf_map,
                                  RNode df, std::string variation);

/*
############################################
JET VETO MAPS
############################################
*/
// Run 2 recommendations: https://cms-jerc.web.cern.ch/Recommendations/#run-2_1
const std::unordered_map<std::string, correction::CorrectionSet> jetVetoMaps = {
    {"2016preVFP", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run2-2016preVFP-UL-NanoAODv9/latest/jetvetomaps.json.gz")},
    {"2016postVFP", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run2-2016postVFP-UL-NanoAODv9/latest/jetvetomaps.json.gz")},
    {"2017", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run2-2017-UL-NanoAODv9/latest/jetvetomaps.json.gz")},
    {"2018", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run2-2018-UL-NanoAODv9/latest/jetvetomaps.json.gz")},
    {"2022Re-recoBCD", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-22CDSep23-Summer22-NanoAODv12/latest/jetvetomaps.json.gz")},
    {"2022Re-recoE+PromptFG", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-22EFGSep23-Summer22EE-NanoAODv12/latest/jetvetomaps.json.gz")},
    {"2023PromptC", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-23CSep23-Summer23-NanoAODv12/latest/jetvetomaps.json.gz")},
    {"2023PromptD", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-23DSep23-Summer23BPix-NanoAODv12/latest/jetvetomaps.json.gz")},
    {"2024Prompt", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/latest/jetvetomaps.json.gz")}
};

const std::unordered_map<std::string, std::string> jetVetoMap_names = {
    {"2016preVFP", "Summer19UL16_V1"},
    {"2016postVFP", "Summer19UL16_V1"},
    {"2017", "Summer19UL17_V1"},
    {"2018", "Summer19UL18_V1"},
    {"2022Re-recoBCD", "Summer22_23Sep2023_RunCD_V1"},
    {"2022Re-recoE+PromptFG", "Summer22EE_23Sep2023_RunEFG_V1"},
    {"2023PromptC", "Summer23Prompt23_RunC_V1"},
    {"2023PromptD", "Summer23BPixPrompt23_RunD_V1"},
    {"2024Prompt", "Summer24Prompt24_RunBCDEFGHI_V1"}
};

RNode applyJetVetoMaps(RNode df);

/*
############################################
ELECTRON SCALE AND SMEARING CORRECTIONS
############################################
*/
const std::unordered_map<std::string, correction::CorrectionSet> electronSSCorrections = {
    {"2016preVFP", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/Run2-2016preVFP-UL-NanoAODv15/latest/electronSS_EtDependent.json.gz")},
    {"2016postVFP", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/Run2-2016postVFP-UL-NanoAODv15/latest/electronSS_EtDependent.json.gz")},
    {"2017", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/Run2-2017-UL-NanoAODv15/latest/electronSS_EtDependent.json.gz")},
    {"2018", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/Run2-2018-UL-NanoAODv15/latest/electronSS_EtDependent.json.gz")},
    {"2022Re-recoBCD", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/Run3-22CDSep23-Summer22-NanoAODv12/latest/electronSS_EtDependent.json.gz")},
    {"2022Re-recoE+PromptFG", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/Run3-22EFGSep23-Summer22EE-NanoAODv12/latest/electronSS_EtDependent.json.gz")},
    {"2023PromptC", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/Run3-23CSep23-Summer23-NanoAODv12/latest/electronSS_EtDependent.json.gz")},
    {"2023PromptD", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/Run3-23DSep23-Summer23BPix-NanoAODv12/latest/electronSS_EtDependent.json.gz")},
    {"2024Prompt", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/EGM/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/latest/electronSS_EtDependent.json.gz")}
};

RNode applyElectronScaleAndSmearing(RNode df, bool isData);

/*
############################################
Others
############################################
*/

RNode HEMCorrection(RNode df, bool isData);

RNode applyDataCorrections(RNode df_);
RNode applyMCCorrections(RNode df_);

#endif // CORRECTIONS_H
