#ifndef CORRECTIONS_H
#define CORRECTIONS_H

#pragma once

#include <unordered_map>
#include <string>
#include <iostream>

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
float bTagWPLoose(std::string year);
float bTagWPMedium(std::string year);
float bTagWPTight(std::string year);

/*
############################################
MET CORRECTIONS
############################################
*/
const std::unordered_map<std::string, correction::CorrectionSet> metCorrections = {
    {"2016preVFP", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run2-2016preVFP-UL-NanoAODv9/latest/met.json.gz")},
    {"2016postVFP", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run2-2016postVFP-UL-NanoAODv9/latest/met.json.gz")},
    {"2017", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run2-2017-UL-NanoAODv9/latest/met.json.gz")},
    {"2018", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run2-2018-UL-NanoAODv9/latest/met.json.gz")},
    {"2022Re-recoBCD", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-22CDSep23-Summer22-NanoAODv12/latest/met_xyCorrections_2022_2022.json.gz")},
    {"2022Re-recoE+PromptFG", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-22EFGSep23-Summer22EE-NanoAODv12/latest/met_xyCorrections_2022_2022EE.json.gz")},
    {"2023PromptC", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-23CSep23-Summer23-NanoAODv12/latest/met_xyCorrections_2023_2023.json.gz")},
    {"2023PromptD", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-23DSep23-Summer23BPix-NanoAODv12/latest/met_xyCorrections_2023_2023BPix.json.gz")}
};
RNode applyMETPhiCorrections(std::unordered_map<std::string, correction::CorrectionSet> cset_met, RNode df);
RNode applyMETUnclusteredCorrections(RNode df, std::string variation);

/*
############################################
JET MASS SCALE AND RESOLUTION CORRECTIONS
############################################
*/
const std::unordered_map<std::string, correction::CorrectionSet> jetMassCorrections = {
    {"2016preVFP", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run2-2016preVFP-UL-NanoAODv9/latest/jmar.json.gz")},
    {"2016postVFP", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run2-2016postVFP-UL-NanoAODv9/latest/jmar.json.gz")},
    {"2017", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run2-2017-UL-NanoAODv9/latest/jmar.json.gz")},
    {"2018", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run2-2018-UL-NanoAODv9/latest/jmar.json.gz")}
};
RNode applyJMSCorrections(std::unordered_map<std::string, correction::CorrectionSet> cset_jms, RNode df, std::string variation);
RNode applyJMRCorrections(std::unordered_map<std::string, correction::CorrectionSet> cset_jmr, RNode df, std::string variation);

/*
############################################
JET ENERGY CORRECTIONS
############################################
*/
const std::unordered_map<std::string, correction::CorrectionSet> jetEnergyCorrections = {
    {"2016preVFP", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run2-2016preVFP-UL-NanoAODv9/latest/jet_jerc.json.gz")},
    {"2016postVFP", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run2-2016postVFP-UL-NanoAODv9/latest/jet_jerc.json.gz")},
    {"2017", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run2-2017-UL-NanoAODv9/latest/jet_jerc.json.gz")},
    {"2018", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run2-2018-UL-NanoAODv9/latest/jet_jerc.json.gz")},
    {"2022Re-recoBCD", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-22CDSep23-Summer22-NanoAODv12/latest/jet_jerc.json.gz")},
    {"2022Re-recoE+PromptFG", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-22EFGSep23-Summer22EE-NanoAODv12/latest/jet_jerc.json.gz")},
    {"2023PromptC", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-23CSep23-Summer23-NanoAODv12/latest/jet_jerc.json.gz")},
    {"2023PromptD", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-23DSep23-Summer23BPix-NanoAODv12/latest/jet_jerc.json.gz")},
    {"2024Prompt", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/latest/jet_jerc.json.gz")}
};

const std::unordered_map<std::string, correction::CorrectionSet> jetEnergyResolution_smear = {
    {"jer_smear", *CorrectionSet::from_file("/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/JER-Smearing/latest/jer_smear.json.gz")}
};

const std::unordered_map<std::string, std::string> jetEnergyCorrections_JEC_prefix = {
    {"2016preVFP", "Summer19UL16APV_V7_MC"},
    {"2016postVFP", "Summer19UL16_V7_MC"},
    {"2017", "Summer19UL17_V5_MC"},
    {"2018", "Summer19UL18_V5_MC"},
    {"2022Re-recoBCD", "Summer22_22Sep2023_V2_MC"},
    {"2022Re-recoE+PromptFG", "Summer22EE_22Sep2023_V2_MC"},
    {"2023PromptC", "Summer23Prompt23_V1_MC"},
    {"2023PromptD", "Summer23BPixPrompt23_V1_MC"},
    {"2024Prompt", "Summer24Prompt24_V1_MC"}
};

const std::unordered_map<std::string, std::string> jetEnergyCorrections_JEC_suffix = {
    {"2016preVFP", "AK4PFchs"},
    {"2016postVFP", "AK4PFchs"},
    {"2017", "AK4PFchs"},
    {"2018", "AK4PFchs"},
    {"2022Re-recoBCD", "AK4PFPuppi"},
    {"2022Re-recoE+PromptFG", "AK4PFPuppi"},
    {"2023PromptC", "AK4PFPuppi"},
    {"2023PromptD", "AK4PFPuppi"},
    {"2024Prompt", "AK4PFPuppi"}
};

const std::unordered_map<std::string, std::string> jetEnergyResolution_JER_res_name = {
    {"2016preVFP", "Summer20UL16APV_JRV3_MC_PtResolution_AK4PFchs"},
    {"2016postVFP", "Summer20UL16_JRV3_MC_PtResolution_AK4PFchs"},
    {"2017", "Summer19UL17_JRV2_MC_PtResolution_AK4PFchs"},
    {"2018", "Summer19UL18_JRV2_MC_PtResolution_AK4PFchs"},
    {"2022Re-recoBCD", "Summer22_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi"},
    {"2022Re-recoE+PromptFG", "Summer22EE_22Sep2023_JRV1_MC_PtResolution_AK4PFPuppi"},
    {"2023PromptC", "Summer23Prompt23_JRV1_MC_PtResolution_AK4PFPuppi"},
    {"2023PromptD", "Summer23BPixPrompt23_JRV1_MC_PtResolution_AK4PFPuppi"},
    {"2024Prompt", "Summer24Prompt24_JRV1_MC_PtResolution_AK4PFPuppi"}
};

const std::unordered_map<std::string, std::string> jetEnergyResolution_JER_sf_name = {
    {"2016preVFP", "Summer20UL16APV_JRV3_MC_ScaleFactor_AK4PFchs"},
    {"2016postVFP", "Summer20UL16_JRV3_MC_ScaleFactor_AK4PFchs"},
    {"2017", "Summer19UL17_JRV2_MC_ScaleFactor_AK4PFchs"},
    {"2018", "Summer19UL18_JRV2_MC_ScaleFactor_AK4PFchs"},
    {"2022Re-recoBCD", "Summer22_22Sep2023_JRV1_MC_ScaleFactor_AK4PFPuppi"},
    {"2022Re-recoE+PromptFG", "Summer22EE_22Sep2023_JRV1_MC_ScaleFactor_AK4PFPuppi"},
    {"2023PromptC", "Summer23Prompt23_JRV1_MC_ScaleFactor_AK4PFPuppi"},
    {"2023PromptD", "Summer23BPixPrompt23_JRV1_MC_ScaleFactor_AK4PFPuppi"},
    {"2024Prompt", "Summer24Prompt24_JRV1_MC_ScaleFactor_AK4PFPuppi"}
};

RNode applyJetEnergyCorrections(std::unordered_map<std::string, correction::CorrectionSet> cset_jerc, std::unordered_map<std::string, std::string> jec_prefix_map, std::unordered_map<std::string, std::string> jec_suffix_map, RNode df, std::string JEC_type, std::string variation);
RNode applyJetEnergyResolution(std::unordered_map<std::string, correction::CorrectionSet> cset_jerc, std::unordered_map<std::string, correction::CorrectionSet> cset_jer_smear, std::unordered_map<std::string, std::string> jer_res_map, std::unordered_map<std::string, std::string> jer_sf_map, RNode df, std::string variation);

/*
############################################
JET VETO MAPS
############################################
*/
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
RNode applyJetVetoMaps(RNode df);

#endif // CORRECTIONS_H