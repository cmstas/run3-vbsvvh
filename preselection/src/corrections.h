#ifndef CORRECTIONS_H
#define CORRECTIONS_H

#pragma once

#include <unordered_map>
#include <string>
#include <iostream>
#include <vector>
#include <array>
#include <stdexcept>

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
JET MASS SCALE (JMS) AND RESOLUTION (JMR) — generic, per-mass-column
############################################

Convention:
  - JMS  : *additive* shift on the target mass column, in GeV. Central = 0.0 GeV.
  - JMR  : *multiplicative* width factor wrt the gen reference mass. Central = 1.0.

applyJetMassScale / applyJetMassResolution (corrections.cpp) are column-agnostic — the
mass branch and (for JMR) the gen mass + gen-index branches are passed in. This
analysis intends to calibrate on the GloParT regressed mass, not FatJet_msoftdrop
(see CORRECTIONS.md § 5-6), and the calibration is not yet derived — the helpers
are currently unwired in applyMCCorrections.

The maps below are identity-valued placeholders kept as a template; the GloParT
calibration will supply its own per-era shift / factor / σ_rel maps. Up/down
systematic variants come from passing ±1σ shift/scale.
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

// Relative mass resolution sigma(m)/m used by the unmatched stochastic branch of
// applyJetMassResolution. Placeholder 1.0 — replace with the per-era value from the same
// calibration fit that derives jetMassResolution_central. With factor f = 1.0 (placeholder)
// the branch is unreachable, so this value has no effect until both maps are updated together.
const std::unordered_map<std::string, float> jetMassResolution_sigmaRel_central = {
    {"2016preVFP", 1.0f}, {"2016postVFP", 1.0f}, {"2017", 1.0f}, {"2018", 1.0f},
    {"2022Re-recoBCD", 1.0f}, {"2022Re-recoE+PromptFG", 1.0f},
    {"2023PromptC", 1.0f}, {"2023PromptD", 1.0f},
    {"2024Prompt", 1.0f}
};

RNode applyJetMassScale(const std::unordered_map<std::string, float>& shift_map,
                     const std::string& mass_col, RNode df);
RNode applyJetMassResolution(const std::unordered_map<std::string, float>& factor_map,
                          const std::unordered_map<std::string, float>& sigma_rel_map,
                          const std::string& mass_col,
                          const std::string& gen_mass_col,
                          const std::string& gen_idx_col,
                          RNode df);

/*
############################################
JET ENERGY CORRECTIONS
############################################

Run 2 + the 2022/2023 Run 3 eras are pinned to the 2026-06-05 JME snapshot (jet_jerc /
fatJet_jerc); the 2024Prompt and 2025 eras are pinned to the newer 2026-07-14 snapshot
(JEC V4/V2, JER JRV2 — see below). JER-Smearing is pinned to 2025-11-03 — NOT the mutable
`latest/` symlink. The tag strings below are the compound/correction keys that exist inside
those pinned files. The "2024Prompt" era uses the Summer24Prompt24 tag (JEC V4 / JER JRV2),
and the "2025" era (2025 data + Summer24 MC) uses the JME-recommended Summer24Prompt25 tag
(JEC V2 / JER JRV2), both pinned to the 2026-07-14 snapshot. The 2026-07-14 update moved the
L2L3Residual corrections to signed-η (from |η|) and refreshed the JER SFs; evalJECCompound
already passes signed eta and resolves inputs by name, so no code change was needed.
correctionlib has no "give me the newest version" API: the code asks
for a tag by name, so the file version and the requested tag string must be bumped
together, deliberately, in one commit. Loading from `latest/` while pinning the tag
string is what previously let the two drift apart (V3 code vs V4 file -> silent no-op;
now a hard error, see applyJetEnergyCorrections in corrections.cpp). To adopt a newer
JME release: pick the new dated dir, find the new recommended tag, and update both
the paths and the tag maps here.

JES
*/
const std::unordered_map<std::string, correction::CorrectionSet> jetEnergyCorrections = {
    {"2016preVFP", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run2-2016preVFP-UL-NanoAODv15/2026-06-05/jet_jerc.json.gz")},
    {"2016postVFP", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run2-2016postVFP-UL-NanoAODv15/2026-06-05/jet_jerc.json.gz")},
    {"2017", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run2-2017-UL-NanoAODv15/2026-06-05/jet_jerc.json.gz")},
    {"2018", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run2-2018-UL-NanoAODv15/2026-06-05/jet_jerc.json.gz")},
    {"2022Re-recoBCD", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-22CDSep23-Summer22-NanoAODv12/2026-06-05/jet_jerc.json.gz")},
    {"2022Re-recoE+PromptFG", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-22EFGSep23-Summer22EE-NanoAODv12/2026-06-05/jet_jerc.json.gz")},
    {"2023PromptC", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-23CSep23-Summer23-NanoAODv12/2026-06-05/jet_jerc.json.gz")},
    {"2023PromptD", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-23DSep23-Summer23BPix-NanoAODv12/2026-06-05/jet_jerc.json.gz")},
    {"2024Prompt", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/2026-07-14/jet_jerc.json.gz")},
    // 2025 data + Summer24 MC — JME-recommended Summer24Prompt25 (JEC + JER + JES self-contained); pinned to 2026-07-03.
    {"2025", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-25Prompt-Summer24-NanoAODv15/2026-07-14/jet_jerc.json.gz")}
};

// AK8 fat-jet JEC/JER lives in fatJet_jerc.json.gz under the same era directories.
const std::unordered_map<std::string, correction::CorrectionSet> fatJetEnergyCorrections = {
    {"2016preVFP", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run2-2016preVFP-UL-NanoAODv15/2026-06-05/fatJet_jerc.json.gz")},
    {"2016postVFP", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run2-2016postVFP-UL-NanoAODv15/2026-06-05/fatJet_jerc.json.gz")},
    {"2017", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run2-2017-UL-NanoAODv15/2026-06-05/fatJet_jerc.json.gz")},
    {"2018", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run2-2018-UL-NanoAODv15/2026-06-05/fatJet_jerc.json.gz")},
    {"2022Re-recoBCD", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-22CDSep23-Summer22-NanoAODv12/2026-06-05/fatJet_jerc.json.gz")},
    {"2022Re-recoE+PromptFG", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-22EFGSep23-Summer22EE-NanoAODv12/2026-06-05/fatJet_jerc.json.gz")},
    {"2023PromptC", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-23CSep23-Summer23-NanoAODv12/2026-06-05/fatJet_jerc.json.gz")},
    {"2023PromptD", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-23DSep23-Summer23BPix-NanoAODv12/2026-06-05/fatJet_jerc.json.gz")},
    {"2024Prompt", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/2026-07-14/fatJet_jerc.json.gz")},
    {"2025", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-25Prompt-Summer24-NanoAODv15/2026-07-14/fatJet_jerc.json.gz")}
};

const std::unordered_map<std::string, correction::CorrectionSet> jetEnergyResolution_smear = {
    {"jer_smear", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/JER-Smearing/2025-11-03/jer_smear.json.gz")}
};

// JEC tag base — `<TAG>_MC` and `<TAG>_DATA` compounds (`_L1L2L3Res_<algo>`) live in jet_jerc.json.gz.
// The trailing `_MC` here is stripped and replaced with `_DATA` at runtime when running on data.
const std::unordered_map<std::string, std::string> jetEnergyCorrections_JEC_prefix = {
    {"2016preVFP", "Summer20UL16APVNanoV15_V1_MC"},
    {"2016postVFP", "Summer20UL16NanoV15_V1_MC"},
    {"2017", "Summer20UL17NanoV15_V1_MC"},
    {"2018", "Summer20UL18NanoV15_V1_MC"},
    {"2022Re-recoBCD", "Summer22_22Sep2023_V4_MC"},
    {"2022Re-recoE+PromptFG", "Summer22EE_22Sep2023_V4_MC"},
    {"2023PromptC", "Summer23Prompt23_V4_MC"},
    {"2023PromptD", "Summer23BPixPrompt23_V4_MC"},
    {"2024Prompt", "Summer24Prompt24_V4_MC"},
    {"2025", "Summer24Prompt25_V2_MC"}
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
    {"2024Prompt", "AK4PFPuppi"},
    {"2025", "AK4PFPuppi"}
};

// AK8 fat-jet variants — same year keys, AK8PFPuppi algo across all eras.
const std::unordered_map<std::string, std::string> fatJetEnergyCorrections_JEC_prefix = {
    {"2016preVFP", "Summer20UL16APVNanoV15_V1_MC"},
    {"2016postVFP", "Summer20UL16NanoV15_V1_MC"},
    {"2017", "Summer20UL17NanoV15_V1_MC"},
    {"2018", "Summer20UL18NanoV15_V1_MC"},
    {"2022Re-recoBCD", "Summer22_22Sep2023_V4_MC"},
    {"2022Re-recoE+PromptFG", "Summer22EE_22Sep2023_V4_MC"},
    {"2023PromptC", "Summer23Prompt23_V4_MC"},
    {"2023PromptD", "Summer23BPixPrompt23_V4_MC"},
    {"2024Prompt", "Summer24Prompt24_V4_MC"},
    {"2025", "Summer24Prompt25_V2_MC"}
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
    {"2024Prompt", "AK8PFPuppi"},
    {"2025", "AK8PFPuppi"}
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
    {"2024Prompt", "2024"},
    {"2025", "2025"}
};

// JER tag names — pinned per era (see JEC section header): Run 2 + 2022/2023 on the
// 2026-06-05 snapshot, 2024Prompt + 2025 on the 2026-07-14 snapshot.
// 2024 ships a dedicated Summer24Prompt24_JRV2 release (no longer reusing 2023BPix).
const std::unordered_map<std::string, std::string> jetEnergyResolution_JER_res_name = {
    {"2016preVFP", "Summer20UL16APV_JRV5_MC_PtResolution_AK4PFPuppi"},
    {"2016postVFP", "Summer20UL16_JRV5_MC_PtResolution_AK4PFPuppi"},
    {"2017", "Summer19UL17_JRV4_MC_PtResolution_AK4PFPuppi"},
    {"2018", "Summer19UL18_JRV3_MC_PtResolution_AK4PFPuppi"},
    {"2022Re-recoBCD", "Summer22_22Sep2023_JRV2_MC_PtResolution_AK4PFPuppi"},
    {"2022Re-recoE+PromptFG", "Summer22EE_22Sep2023_JRV2_MC_PtResolution_AK4PFPuppi"},
    {"2023PromptC", "Summer23Prompt23_RunCv1234_JRV2_MC_PtResolution_AK4PFPuppi"},
    {"2023PromptD", "Summer23BPixPrompt23_RunD_JRV2_MC_PtResolution_AK4PFPuppi"},
    {"2024Prompt", "Summer24Prompt24_JRV2_MC_PtResolution_AK4PFPuppi"},
    {"2025", "Summer24Prompt25_JRV2_MC_PtResolution_AK4PFPuppi"}
};

const std::unordered_map<std::string, std::string> jetEnergyResolution_JER_sf_name = {
    {"2016preVFP", "Summer20UL16APV_JRV5_MC_ScaleFactor_AK4PFPuppi"},
    {"2016postVFP", "Summer20UL16_JRV5_MC_ScaleFactor_AK4PFPuppi"},
    {"2017", "Summer19UL17_JRV4_MC_ScaleFactor_AK4PFPuppi"},
    {"2018", "Summer19UL18_JRV3_MC_ScaleFactor_AK4PFPuppi"},
    {"2022Re-recoBCD", "Summer22_22Sep2023_JRV2_MC_ScaleFactor_AK4PFPuppi"},
    {"2022Re-recoE+PromptFG", "Summer22EE_22Sep2023_JRV2_MC_ScaleFactor_AK4PFPuppi"},
    {"2023PromptC", "Summer23Prompt23_RunCv1234_JRV2_MC_ScaleFactor_AK4PFPuppi"},
    {"2023PromptD", "Summer23BPixPrompt23_RunD_JRV2_MC_ScaleFactor_AK4PFPuppi"},
    {"2024Prompt", "Summer24Prompt24_JRV2_MC_ScaleFactor_AK4PFPuppi"},
    {"2025", "Summer24Prompt25_JRV2_MC_ScaleFactor_AK4PFPuppi"}
};

// AK8 JER tags — same JR releases as AK4 but with AK8PFPuppi algo.
const std::unordered_map<std::string, std::string> fatJetEnergyResolution_JER_res_name = {
    {"2016preVFP", "Summer20UL16APV_JRV5_MC_PtResolution_AK8PFPuppi"},
    {"2016postVFP", "Summer20UL16_JRV5_MC_PtResolution_AK8PFPuppi"},
    {"2017", "Summer19UL17_JRV4_MC_PtResolution_AK8PFPuppi"},
    {"2018", "Summer19UL18_JRV3_MC_PtResolution_AK8PFPuppi"},
    {"2022Re-recoBCD", "Summer22_22Sep2023_JRV2_MC_PtResolution_AK8PFPuppi"},
    {"2022Re-recoE+PromptFG", "Summer22EE_22Sep2023_JRV2_MC_PtResolution_AK8PFPuppi"},
    {"2023PromptC", "Summer23Prompt23_RunCv1234_JRV2_MC_PtResolution_AK8PFPuppi"},
    {"2023PromptD", "Summer23BPixPrompt23_RunD_JRV2_MC_PtResolution_AK8PFPuppi"},
    {"2024Prompt", "Summer24Prompt24_JRV2_MC_PtResolution_AK8PFPuppi"},
    {"2025", "Summer24Prompt25_JRV2_MC_PtResolution_AK8PFPuppi"}
};

const std::unordered_map<std::string, std::string> fatJetEnergyResolution_JER_sf_name = {
    {"2016preVFP", "Summer20UL16APV_JRV5_MC_ScaleFactor_AK8PFPuppi"},
    {"2016postVFP", "Summer20UL16_JRV5_MC_ScaleFactor_AK8PFPuppi"},
    {"2017", "Summer19UL17_JRV4_MC_ScaleFactor_AK8PFPuppi"},
    {"2018", "Summer19UL18_JRV3_MC_ScaleFactor_AK8PFPuppi"},
    {"2022Re-recoBCD", "Summer22_22Sep2023_JRV2_MC_ScaleFactor_AK8PFPuppi"},
    {"2022Re-recoE+PromptFG", "Summer22EE_22Sep2023_JRV2_MC_ScaleFactor_AK8PFPuppi"},
    {"2023PromptC", "Summer23Prompt23_RunCv1234_JRV2_MC_ScaleFactor_AK8PFPuppi"},
    {"2023PromptD", "Summer23BPixPrompt23_RunD_JRV2_MC_ScaleFactor_AK8PFPuppi"},
    {"2024Prompt", "Summer24Prompt24_JRV2_MC_ScaleFactor_AK8PFPuppi"},
    {"2025", "Summer24Prompt25_JRV2_MC_ScaleFactor_AK8PFPuppi"}
};

// Full Type-1 PuppiMET rebuild from RawPuppiMET over the merged AK4 list (Jet + CorrT1METJet)
// with muon subtraction, per the JERC recipe. Writes met_pt / met_phi (nominal) and, on MC with
// systematics enabled, met_pt_<sfx> / met_phi_<sfx> for each JES variation. Must run AFTER
// applyJetEnergyCorrections / applyJetEnergyResolution / applyJESVariations (it consumes
// Jet_jerFactor and the _Jet_shift_<sfx> columns and re-evaluates the JEC compound) and after
// _Jet_rawpt = (1-Jet_rawFactor)*Jet_pt has been defined on the pristine NanoAOD jets.
RNode applyType1MET(RNode df, bool isData);

// Nominal JEC: removes the JEC stored in NanoAOD (via Jet_rawFactor) and re-applies the
// latest L1FastJet * L2Relative * L3Absolute (* L2L3Residual for data) compound from
// jet_jerc.json.gz. MET is NOT propagated here — it is rebuilt in applyType1MET.
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
