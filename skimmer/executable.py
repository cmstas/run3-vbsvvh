#!/usr/bin/python3
# -*- coding: utf-8 -*-

import re
import os
import sys
import glob
import time
import shutil
import subprocess
from argparse import ArgumentParser

import ROOT as r

r.gInterpreter.Declare("""
int find_matching_jet(int target_idx, float target_eta, float target_phi, ROOT::RVec<int> already_matched_jet_indices, ROOT::RVec<int> already_matched_fatjet_indices, ROOT::RVec<float> jet_eta, ROOT::RVec<float> jet_phi, ROOT::RVec<float> fatjet_eta, ROOT::RVec<float> fatjet_phi) {
    int max_jets = 10;
    int max_fatjets = 3;
    const float dR_cut = 0.4f;
    const float fatjet_overlap_cut = 0.8f;
    const float jet_overlap_cut = 0.4f;
                       
    if (target_idx < 0) {
        return -1;
    }
    
    // Calculate dR values for all jets
    ROOT::RVec<float> dR_values;
    for (size_t i = 0; i < jet_eta.size(); ++i) {
        dR_values.push_back(ROOT::VecOps::DeltaR(target_eta, jet_eta[i], target_phi, jet_phi[i]));
    }
    
    auto sorted_indices = ROOT::VecOps::Argsort(dR_values);
    for (int idx : sorted_indices) {
        if (dR_values[idx] >= dR_cut) break; // No more candidates within dR cut
        if (idx >= max_jets) continue; // Skip indices beyond padding limit
        
        // Check if this jet is in the excluded list
        bool is_excluded = false;
        for (int excl_idx : already_matched_jet_indices) {
            if (idx == excl_idx) {
                is_excluded = true;
                break;
            }
        }
        if (is_excluded) continue;
        
        bool overlaps_jet = false;
        for (int matched_j_idx : already_matched_jet_indices) {
            if (matched_j_idx >= 0 && matched_j_idx < max_jets) {
                float dR_jet = ROOT::VecOps::DeltaR(jet_eta[idx], jet_eta[matched_j_idx], jet_phi[idx], jet_phi[matched_j_idx]);
                if (dR_jet < jet_overlap_cut) {
                    overlaps_jet = true;
                    break;
                }
            }
        }

        // Check overlap with matched fatjets only
        bool overlaps_fatjet = false;
        for (int matched_fj_idx : already_matched_fatjet_indices) {
            if (matched_fj_idx >= 0 && matched_fj_idx < max_fatjets) {
                float dR_fatjet = ROOT::VecOps::DeltaR(jet_eta[idx], fatjet_eta[matched_fj_idx], jet_phi[idx], fatjet_phi[matched_fj_idx]);
                if (dR_fatjet < fatjet_overlap_cut) {
                    overlaps_fatjet = true;
                    break;
                }
            }
        }
        if (!overlaps_jet && !overlaps_fatjet) {
            return idx;
        }
    }
    return -1;
}

int find_matching_fatjet(int target_idx, float target_eta, float target_phi, ROOT::RVec<int> already_matched_jet_indices, ROOT::RVec<int> already_matched_fatjet_indices, ROOT::RVec<float> jet_eta, ROOT::RVec<float> jet_phi, ROOT::RVec<float> fatjet_eta, ROOT::RVec<float> fatjet_phi) {
    int max_jets = 10;
    int max_fatjets = 3;
    const float dR_cut = 0.8f;
    const float fatjet_overlap_cut = 0.8f;
    const float jet_overlap_cut = 0.8f;
    
    if (target_idx < 0) {
        return -1;
    }
                       
    // Calculate dR values for all fatjets
    ROOT::RVec<float> dR_values;
    for (size_t i = 0; i < fatjet_eta.size(); ++i) {
        dR_values.push_back(ROOT::VecOps::DeltaR(target_eta, fatjet_eta[i], target_phi, fatjet_phi[i]));
    }
    
    auto sorted_indices = ROOT::VecOps::Argsort(dR_values);
    for (int idx : sorted_indices) {
        if (dR_values[idx] >= dR_cut) break; // No more candidates within dR cut
        if (idx >= max_fatjets) continue; // Skip indices beyond padding limit
        
        // Check if this fatjet is in the excluded list
        bool is_excluded = false;
        for (int excl_idx : already_matched_fatjet_indices) {
            if (idx == excl_idx) {
                is_excluded = true;
                break;
            }
        }
        if (is_excluded) continue;

        // Check overlap with matched jets only
        bool overlaps_jet = false;
        for (int matched_j_idx : already_matched_jet_indices) {
            if (matched_j_idx >= 0 && matched_j_idx < max_jets) {
                float dR_jet = ROOT::VecOps::DeltaR(fatjet_eta[idx], jet_eta[matched_j_idx], fatjet_phi[idx], jet_phi[matched_j_idx]);
                if (dR_jet < jet_overlap_cut) {
                    overlaps_jet = true;
                    break;
                }
            }
        }
        
        bool overlaps_fatjet = false;
        for (int matched_fj_idx : already_matched_fatjet_indices) {
            if (matched_fj_idx >= 0 && matched_fj_idx < max_fatjets) {
                float dR_fatjet = ROOT::VecOps::DeltaR(fatjet_eta[idx], fatjet_eta[matched_fj_idx], fatjet_phi[idx], fatjet_phi[matched_fj_idx]);
                if (dR_fatjet < fatjet_overlap_cut) {
                    overlaps_fatjet = true;
                    break;
                }
            }
        }

        if (!overlaps_jet && !overlaps_fatjet) {
            return idx;
        }
    }
    return -1;
}
                       
ROOT::RVec<int> findHiggsAndDaughters(ROOT::RVec<int>& pdgId, ROOT::RVec<int>& status, ROOT::RVec<short>& motherIdx) {
    ROOT::RVec<int> result = {-1, -1, -1};
    
    int higgs_idx = -1;
    ROOT::RVec<int> hdecay_idx;
    
    // Find intermediate Higgs with status 22 and mother_idx 0
    for (size_t igen = 0; igen < pdgId.size(); ++igen) {
        int part_status = status[igen];
        int part_pdgId = pdgId[igen];
        int mother_idx = motherIdx[igen];
        
        if (mother_idx == 0 && part_status == 22 && part_pdgId == 25) {
            higgs_idx = igen;
            break;
        }
    }
    
    if (higgs_idx == -1) return result;
    
    // Find first daughters of Higgs (non-Higgs particles with Higgs mother)
    for (size_t igen = 0; igen < pdgId.size(); ++igen) {
        int mother_idx = motherIdx[igen];
        if (mother_idx > 0 && mother_idx < pdgId.size()) {
            int mother_pdgId = pdgId[mother_idx];
            int part_pdgId = pdgId[igen];
            
            if (mother_pdgId == 25 && part_pdgId != 25) {
                hdecay_idx.push_back(igen);
            }
        }
    }
    
    result[0] = higgs_idx;
    if (hdecay_idx.size() >= 1) result[1] = hdecay_idx[0];
    if (hdecay_idx.size() >= 2) result[2] = hdecay_idx[1];
    
    return result;
}

int findLastIndex(int current_idx, int current_pdgId, ROOT::RVec<int>& pdgId, ROOT::RVec<short>& motherIdx) {
    int outIdx = current_idx;
    for (size_t igen = 0; igen < pdgId.size(); ++igen) {
        int part_pdgId = pdgId[igen];
        int mother_idx = motherIdx[igen];
        
        if (mother_idx == current_idx) {
            if (part_pdgId == current_pdgId) {
                outIdx = findLastIndex(igen, part_pdgId, pdgId, motherIdx);
            }
        }
    }
    return outIdx;
}

ROOT::RVec<int> findVBosonsAndDaughters(ROOT::RVec<int>& pdgId, ROOT::RVec<int>& status, ROOT::RVec<short>& motherIdx) {
    // Returns: [v1_idx, v1_d1_idx, v1_d2_idx, v2_idx, v2_d1_idx, v2_d2_idx] or [-1, -1, -1, -1, -1, -1]
    ROOT::RVec<int> result = {-1, -1, -1, -1, -1, -1};
    
    ROOT::RVec<int> firstVs_idx;
    
    for (size_t igen = 0; igen < pdgId.size(); ++igen) {
        int part_status = status[igen];
        int part_pdgId = pdgId[igen];
        int mother_idx = motherIdx[igen];
        
        if (mother_idx == 0 && part_status == 22 && (part_pdgId == 23 || abs(part_pdgId) == 24)) {
            firstVs_idx.push_back(igen);
        }
    }
    
    if (firstVs_idx.size() < 2) return result;
    
    ROOT::RVec<int> hadronic_V_indices;
    ROOT::RVec<int> leptonic_V_indices;
    
    for (size_t iV = 0; iV < firstVs_idx.size() && iV < 2; ++iV) {
        int firstV_idx = firstVs_idx[iV];
        int firstV_pdgId = pdgId[firstV_idx];
        
        int lastV_idx = findLastIndex(firstV_idx, firstV_pdgId, pdgId, motherIdx);
        
        int hadronic_daughters = 0;
        int leptonic_daughters = 0;
        
        for (size_t igen = 0; igen < pdgId.size(); ++igen) {
            int mother_idx = motherIdx[igen];
            if (mother_idx == lastV_idx) {
                if (abs(pdgId[igen]) <= 6) {
                    hadronic_daughters++;
                } else {
                    leptonic_daughters++;
                }
            }
        }
        
        if (hadronic_daughters > 0) {
            hadronic_V_indices.push_back(iV);
        } else if (leptonic_daughters > 0) {
            leptonic_V_indices.push_back(iV);
        }
    }
    
    if (hadronic_V_indices.size() > 0) {
        int iV = hadronic_V_indices[0];
        int firstV_idx = firstVs_idx[iV];
        int firstV_pdgId = pdgId[firstV_idx];
        
        int lastV_idx = findLastIndex(firstV_idx, firstV_pdgId, pdgId, motherIdx);
        
        ROOT::RVec<int> vdecays_idx;
        
        for (size_t igen = 0; igen < pdgId.size(); ++igen) {
            int mother_idx = motherIdx[igen];
            if (mother_idx == lastV_idx) {
                if (abs(pdgId[igen]) <= 6) { // quarks only
                   vdecays_idx.push_back(igen);
                }
            }
        }
        
        if (vdecays_idx.size() != 0) {
            result[0] = lastV_idx;
        }
        if (vdecays_idx.size() >= 1) result[1] = vdecays_idx[0];
        if (vdecays_idx.size() >= 2) result[2] = vdecays_idx[1];
    }
    
    if (hadronic_V_indices.size() >= 2) {
        int iV = hadronic_V_indices[1];
        int firstV_idx = firstVs_idx[iV];
        int firstV_pdgId = pdgId[firstV_idx];
        
        int lastV_idx = findLastIndex(firstV_idx, firstV_pdgId, pdgId, motherIdx);
        
        ROOT::RVec<int> vdecays_idx;
        
        for (size_t igen = 0; igen < pdgId.size(); ++igen) {
            int mother_idx = motherIdx[igen];
            if (mother_idx == lastV_idx) {
                if (abs(pdgId[igen]) <= 6) {
                   vdecays_idx.push_back(igen);
                }
            }
        }
        
        if (vdecays_idx.size() != 0) {
            result[3] = lastV_idx;
        }
        if (vdecays_idx.size() >= 1) result[4] = vdecays_idx[0];
        if (vdecays_idx.size() >= 2) result[5] = vdecays_idx[1];
    }
    
    return result;
}

ROOT::RVec<int> findVBSQuarks(ROOT::RVec<int>& pdgId, ROOT::RVec<int>& status, ROOT::RVec<short>& motherIdx) {
    // Returns: [vbs1_idx, vbs2_idx] or [-1, -1]
    ROOT::RVec<int> result = {-1, -1};
    ROOT::RVec<int> vbsquarks_idx;
    
    // Find outgoing quarks with status 23 and mother_idx 0
    for (size_t igen = 0; igen < pdgId.size(); ++igen) {
        int part_status = status[igen];
        int part_pdgId = pdgId[igen];
        int mother_idx = motherIdx[igen];
        
        if (mother_idx == 0 && part_status == 23 && abs(part_pdgId) <= 6 && abs(part_pdgId) >= 1) {
            vbsquarks_idx.push_back(igen);
        }
    }
    
    if (vbsquarks_idx.size() >= 1) result[0] = vbsquarks_idx[0];
    if (vbsquarks_idx.size() >= 2) result[1] = vbsquarks_idx[1];
    
    return result;
}
""")

subprocess.run("python3 -m pip install --user --no-binary=correctionlib correctionlib", shell=True, check=True)
import importlib
correctionlib = importlib.import_module("correctionlib")
correctionlib.register_pyroot_binding()

# Constants
CONDOR_OUTPUT_DIR = "output"
XROOTD_REDIRECTOR = "root://xrootd-cms.infn.it/"
OUTPUT_XRD = "davs://redirector.t2.ucsd.edu:1095//store/user/aaarora/skims"
MAX_RETRIES = 10
SLEEP_DURATION = 60  # 1 minute in seconds

JET_ID_JSONS = {"2024": "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/2024_Summer24/jetid.json.gz"}

class Skimmer():
    def __init__(self, inFiles, outDir, keepDropFile):
        self.inFiles = inFiles
        self.outDir = outDir
        self.keepDropFile = keepDropFile

        self.df = r.RDataFrame("Events", self.inFiles)
        r.RDF.Experimental.AddProgressBar(self.df)
        columns = self.df.GetColumnNames()
        for col in columns:
            if col.startswith("Muon_") or col.startswith("Electron_") or col.startswith("Jet_") or col.startswith("FatJet_"):
                colType = self.df.GetColumnType(col)
                if colType == "ROOT::VecOps::RVec<Float_t>":
                    self.df = self.df.DefaultValueFor(col, r.RVec('Float_t')())
                elif colType == "ROOT::VecOps::RVec<Int_t>":
                    self.df = self.df.DefaultValueFor(col, r.RVec('Int_t')())
                elif colType == "ROOT::VecOps::RVec<UShort_t>":
                    self.df = self.df.DefaultValueFor(col, r.RVec('UShort_t')())
                elif colType == "ROOT::VecOps::RVec<Bool_t>":
                    self.df = self.df.DefaultValueFor(col, r.RVec('Bool_t')())
                elif colType == "ROOT::VecOps::RVec<UChar_t>":
                    self.df = self.df.DefaultValueFor(col, r.RVec('UChar_t')())
                elif colType == "ROOT::VecOps::RVec<Double_t>":
                    self.df = self.df.DefaultValueFor(col, r.RVec('Double_t')())
                elif colType == "ROOT::VecOps::RVec<Long64_t>":
                    self.df = self.df.DefaultValueFor(col, r.RVec('Long64_t')())
                elif colType == "ROOT::VecOps::RVec<Short_t>":
                    self.df = self.df.DefaultValueFor(col, r.RVec('Short_t')())
                else:
                    print(f"Unknown column type: {colType} for column {col}")

        self.df.Display(["Electron_pt"]).Print()

    @staticmethod
    def genSelection(df):
        df = df.Define("ak4_jet_mask", "Jet_pt > 20 && abs(Jet_eta) <= 5 && Jet_jetId > 0") \
            .Define("ak4jet_pt", "Jet_pt[ak4_jet_mask]") \
            .Define("ak4jet_eta", "Jet_eta[ak4_jet_mask]") \
            .Define("ak4jet_phi", "Jet_phi[ak4_jet_mask]") \
            .Define("ak4jet_mass", "Jet_mass[ak4_jet_mask]") \
            .Define("ak4jet_btagDeepFlavB", "Jet_btagDeepFlavB[ak4_jet_mask]")

        # fatjet preselection
        df = df.Define("ak8_jet_mask", "FatJet_pt > 200 && abs(FatJet_eta) <= 2.5 && FatJet_mass > 50 && FatJet_msoftdrop > 40 && FatJet_jetId > 0") \
            .Define("ak8jet_pt", "FatJet_pt[ak8_jet_mask]") \
            .Define("ak8jet_eta", "FatJet_eta[ak8_jet_mask]") \
            .Define("ak8jet_phi", "FatJet_phi[ak8_jet_mask]") \
            .Define("ak8jet_mass", "FatJet_mass[ak8_jet_mask]") \
            .Define("ak8jet_nConstituents", "FatJet_nConstituents[ak8_jet_mask]") \
            .Define("ak8jet_globalParT3_Xbb", "FatJet_globalParT3_Xbb[ak8_jet_mask]") \
            .Define("ak8jet_globalParT3_Xcc", "FatJet_globalParT3_Xcc[ak8_jet_mask]") \
            .Define("ak8jet_globalParT3_Xcs", "FatJet_globalParT3_Xcs[ak8_jet_mask]") \
            .Define("ak8jet_globalParT3_Xqq", "FatJet_globalParT3_Xqq[ak8_jet_mask]") \
            .Define("ak8jet_globalParT3_QCD", "FatJet_globalParT3_QCD[ak8_jet_mask]")
        
        # define bb, qq definition
        df = df.Define("higgs_info", "findHiggsAndDaughters(GenPart_pdgId, GenPart_status, GenPart_genPartIdxMother)") \
            .Define("gen_h_idx", "higgs_info[0]") \
            .Define("Higgs_eta", "gen_h_idx >= 0 ? GenPart_eta[gen_h_idx] : -999.0") \
            .Define("Higgs_phi", "gen_h_idx >= 0 ? GenPart_phi[gen_h_idx] : -999.0") \
            .Define("gen_b1_idx", "higgs_info[1]") \
            .Define("b1_eta", "gen_b1_idx >= 0 ? GenPart_eta[gen_b1_idx] : -999.0") \
            .Define("b1_phi", "gen_b1_idx >= 0 ? GenPart_phi[gen_b1_idx] : -999.0") \
            .Define("gen_b2_idx", "higgs_info[2]") \
            .Define("b2_eta", "gen_b2_idx >= 0 ? GenPart_eta[gen_b2_idx] : -999.0") \
            .Define("b2_phi", "gen_b2_idx >= 0 ? GenPart_phi[gen_b2_idx] : -999.0")

        df = df.Define("vboson_info", "findVBosonsAndDaughters(GenPart_pdgId, GenPart_status, GenPart_genPartIdxMother)") \
            .Define("gen_v1_idx", "vboson_info[0]") \
            .Define("V1_eta", "gen_v1_idx >= 0 ? GenPart_eta[gen_v1_idx] : -999.0") \
            .Define("V1_phi", "gen_v1_idx >= 0 ? GenPart_phi[gen_v1_idx] : -999.0") \
            .Define("gen_v1q1_idx", "vboson_info[1]") \
            .Define("V1_q1_eta", "gen_v1q1_idx >= 0 ? GenPart_eta[gen_v1q1_idx] : -999.0") \
            .Define("V1_q1_phi", "gen_v1q1_idx >= 0 ? GenPart_phi[gen_v1q1_idx] : -999.0") \
            .Define("gen_v1q2_idx", "vboson_info[2]") \
            .Define("V1_q2_eta", "gen_v1q2_idx >= 0 ? GenPart_eta[gen_v1q2_idx] : -999.0") \
            .Define("V1_q2_phi", "gen_v1q2_idx >= 0 ? GenPart_phi[gen_v1q2_idx] : -999.0") \
            .Define("gen_v2_idx", "vboson_info[3]") \
            .Define("V2_eta", "gen_v2_idx >= 0 ? GenPart_eta[gen_v2_idx] : -999.0") \
            .Define("V2_phi", "gen_v2_idx >= 0 ? GenPart_phi[gen_v2_idx] : -999.0") \
            .Define("gen_v2q1_idx", "vboson_info[4]") \
            .Define("V2_q1_eta", "gen_v2q1_idx >= 0 ? GenPart_eta[gen_v2q1_idx] : -999.0") \
            .Define("V2_q1_phi", "gen_v2q1_idx >= 0 ? GenPart_phi[gen_v2q1_idx] : -999.0") \
            .Define("gen_v2q2_idx", "vboson_info[5]") \
            .Define("V2_q2_eta", "gen_v2q2_idx >= 0 ? GenPart_eta[gen_v2q2_idx] : -999.0") \
            .Define("V2_q2_phi", "gen_v2q2_idx >= 0 ? GenPart_phi[gen_v2q2_idx] : -999.0")

        df = df.Define("vbs_info", "findVBSQuarks(GenPart_pdgId, GenPart_status, GenPart_genPartIdxMother)") \
            .Define("gen_vbs_idx1", "vbs_info[0]") \
            .Define("VBSJet1_eta", "vbs_info[0] >= 0 ? GenPart_eta[vbs_info[0]] : -999.0") \
            .Define("VBSJet1_phi", "vbs_info[0] >= 0 ? GenPart_phi[vbs_info[0]] : -999.0") \
            .Define("gen_vbs_idx2", "vbs_info[1]") \
            .Define("VBSJet2_eta", "vbs_info[1] >= 0 ? GenPart_eta[vbs_info[1]] : -999.0") \
            .Define("VBSJet2_phi", "vbs_info[1] >= 0 ? GenPart_phi[vbs_info[1]] : -999.0")

        df = df.Define("vbs1_idx_temp", "find_matching_jet(gen_vbs_idx1, VBSJet1_eta, VBSJet1_phi, ROOT::RVec<int>{}, ROOT::RVec<int>{}, ak4jet_eta, ak4jet_phi, ak8jet_eta, ak8jet_phi)") \
            .Define("vbs2_idx_temp", "find_matching_jet(gen_vbs_idx2, VBSJet2_eta, VBSJet2_phi, ROOT::RVec<int>{vbs1_idx_temp}, ROOT::RVec<int>{}, ak4jet_eta, ak4jet_phi, ak8jet_eta, ak8jet_phi)") \
            .Define("truth_vbs1_idx", "vbs1_idx_temp >= 0 && vbs1_idx_temp < 10 ? vbs1_idx_temp : -1") \
            .Define("truth_vbs2_idx", "vbs2_idx_temp >= 0 && vbs2_idx_temp < 10 ? vbs2_idx_temp : -1")

        df = df.Define("hbb_dR", "ROOT::VecOps::DeltaR(b1_eta, b2_eta, b1_phi, b2_phi)") \
            .Define("hbb_fatjet_idx_temp", "find_matching_fatjet(gen_h_idx, Higgs_eta, Higgs_phi, ROOT::RVec<int>{truth_vbs1_idx, truth_vbs2_idx}, ROOT::RVec<int>{}, ak4jet_eta, ak4jet_phi, ak8jet_eta, ak8jet_phi)") \
            .Define("hbb_fatjet_candidate_b1_dR", "hbb_fatjet_idx_temp >= 0 && hbb_fatjet_idx_temp < ak8jet_eta.size() ? ROOT::VecOps::DeltaR(ak8jet_eta[hbb_fatjet_idx_temp], b1_eta, ak8jet_phi[hbb_fatjet_idx_temp], b1_phi) : 999.0") \
            .Define("hbb_fatjet_candidate_b2_dR", "hbb_fatjet_idx_temp >= 0 && hbb_fatjet_idx_temp < ak8jet_eta.size() ? ROOT::VecOps::DeltaR(ak8jet_eta[hbb_fatjet_idx_temp], b2_eta, ak8jet_phi[hbb_fatjet_idx_temp], b2_phi) : 999.0") \
            .Define("hbb_isBoosted", "hbb_fatjet_idx_temp != -1 && hbb_fatjet_idx_temp < 3 && hbb_dR < 0.8 && hbb_fatjet_candidate_b1_dR < 0.8 && hbb_fatjet_candidate_b2_dR < 0.8") \
            .Define("truth_bh_idx", "hbb_isBoosted ? hbb_fatjet_idx_temp : -1")
        
        df = df.Define("excluded_jet_mask_for_hbb", "ROOT::RVec<int>{truth_vbs1_idx, truth_vbs2_idx}") \
            .Define("matched_fatjet_mask_for_hbb", "ROOT::RVec<int>{}") \
            .Define("b1_idx_temp", "find_matching_jet(gen_b1_idx, b1_eta, b1_phi, excluded_jet_mask_for_hbb, matched_fatjet_mask_for_hbb, ak4jet_eta, ak4jet_phi, ak8jet_eta, ak8jet_phi)") \
            .Define("excluded_jet_mask_for_hbb2", "ROOT::RVec<int>{truth_vbs1_idx, truth_vbs2_idx, b1_idx_temp}") \
            .Define("b2_idx_temp", "find_matching_jet(gen_b2_idx, b2_eta, b2_phi, excluded_jet_mask_for_hbb2, matched_fatjet_mask_for_hbb, ak4jet_eta, ak4jet_phi, ak8jet_eta, ak8jet_phi)") \
            .Define("truth_b1_idx", "b1_idx_temp >= 0 && b1_idx_temp < 10 ? b1_idx_temp : -1") \
            .Define("truth_b2_idx", "b2_idx_temp >= 0 && b2_idx_temp < 10 ? b2_idx_temp : -1")

        df = df.Define("v1qq_dR", "gen_v1_idx != -1 ? ROOT::VecOps::DeltaR(V1_q1_eta, V1_q1_phi, V1_q2_eta, V1_q2_phi) : 999.0") \
            .Define("excluded_fatjet_mask_for_v1", "ROOT::RVec<int>{truth_bh_idx}") \
            .Define("matched_jet_mask_for_v1", "ROOT::RVec<int>{truth_vbs1_idx, truth_vbs2_idx, truth_b1_idx, truth_b2_idx}") \
            .Define("v1qq_fatjet_idx_temp", "find_matching_fatjet(gen_v1_idx, V1_eta, V1_phi, matched_jet_mask_for_v1, excluded_fatjet_mask_for_v1, ak4jet_eta, ak4jet_phi, ak8jet_eta, ak8jet_phi)") \
            .Define("v1qq_fatjet_candidate_q1_dR", "v1qq_fatjet_idx_temp >= 0 && v1qq_fatjet_idx_temp < ak8jet_eta.size() ? ROOT::VecOps::DeltaR(ak8jet_eta[v1qq_fatjet_idx_temp], V1_q1_eta, ak8jet_phi[v1qq_fatjet_idx_temp], V1_q1_phi) : 999.0") \
            .Define("v1qq_fatjet_candidate_q2_dR", "v1qq_fatjet_idx_temp >= 0 && v1qq_fatjet_idx_temp < ak8jet_eta.size() ? ROOT::VecOps::DeltaR(ak8jet_eta[v1qq_fatjet_idx_temp], V1_q2_eta, ak8jet_phi[v1qq_fatjet_idx_temp], V1_q2_phi) : 999.0") \
            .Define("v1qq_isBoosted", "gen_v1_idx != -1 && v1qq_fatjet_idx_temp != -1 && v1qq_fatjet_idx_temp < 3 && v1qq_dR < 0.8 && v1qq_fatjet_candidate_q1_dR < 0.8 && v1qq_fatjet_candidate_q2_dR < 0.8") \
            .Define("truth_bv1_idx", "v1qq_isBoosted ? v1qq_fatjet_idx_temp : -1")

        df = df.Define("v2qq_dR", "gen_v2_idx != -1 ? ROOT::VecOps::DeltaR(V2_q1_eta, V2_q1_phi, V2_q2_eta, V2_q2_phi) : 999.0") \
            .Define("excluded_fatjet_mask_for_v2", "ROOT::RVec<int>{truth_bh_idx, truth_bv1_idx}") \
            .Define("matched_jet_mask_for_v2", "ROOT::RVec<int>{truth_vbs1_idx, truth_vbs2_idx, truth_b1_idx, truth_b2_idx}") \
            .Define("v2qq_fatjet_idx_temp", "find_matching_fatjet(gen_v2_idx, V2_eta, V2_phi, matched_jet_mask_for_v2, excluded_fatjet_mask_for_v2, ak4jet_eta, ak4jet_phi, ak8jet_eta, ak8jet_phi)") \
            .Define("v2qq_fatjet_candidate_q1_dR", "v2qq_fatjet_idx_temp >= 0 && v2qq_fatjet_idx_temp < ak8jet_eta.size() ? ROOT::VecOps::DeltaR(ak8jet_eta[v2qq_fatjet_idx_temp], V2_q1_eta, ak8jet_phi[v2qq_fatjet_idx_temp], V2_q1_phi) : 999.0") \
            .Define("v2qq_fatjet_candidate_q2_dR", "v2qq_fatjet_idx_temp >= 0 && v2qq_fatjet_idx_temp < ak8jet_eta.size() ? ROOT::VecOps::DeltaR(ak8jet_eta[v2qq_fatjet_idx_temp], V2_q2_eta, ak8jet_phi[v2qq_fatjet_idx_temp], V2_q2_phi) : 999.0") \
            .Define("v2qq_isBoosted", "gen_v2_idx != -1 && v2qq_fatjet_idx_temp != -1 && v2qq_fatjet_idx_temp < 3 && v2qq_dR < 0.8 && v2qq_fatjet_candidate_q1_dR < 0.8 && v2qq_fatjet_candidate_q2_dR < 0.8") \
            .Define("truth_bv2_idx", "v2qq_isBoosted ? v2qq_fatjet_idx_temp : -1")
        
        df = df.Define("excluded_jet_mask_for_v1q1", "ROOT::RVec<int>{truth_vbs1_idx, truth_vbs2_idx, truth_b1_idx, truth_b2_idx}") \
            .Define("matched_fatjet_mask_for_v1q", "ROOT::RVec<int>{truth_bh_idx, truth_bv1_idx, truth_bv2_idx}") \
            .Define("v1q1_idx_temp", "find_matching_jet(gen_v1q1_idx, V1_q1_eta, V1_q1_phi, excluded_jet_mask_for_v1q1, matched_fatjet_mask_for_v1q, ak4jet_eta, ak4jet_phi, ak8jet_eta, ak8jet_phi)") \
            .Define("excluded_jet_mask_for_v1q2", "ROOT::RVec<int>{truth_vbs1_idx, truth_vbs2_idx, truth_b1_idx, truth_b2_idx, v1q1_idx_temp}") \
            .Define("v1q2_idx_temp", "find_matching_jet(gen_v1q2_idx, V1_q2_eta, V1_q2_phi, excluded_jet_mask_for_v1q2, matched_fatjet_mask_for_v1q, ak4jet_eta, ak4jet_phi, ak8jet_eta, ak8jet_phi)") \
            .Define("truth_v1q1_idx", "v1q1_idx_temp >= 0 && v1q1_idx_temp < 10 ? v1q1_idx_temp : -1") \
            .Define("truth_v1q2_idx", "v1q2_idx_temp >= 0 && v1q2_idx_temp < 10 ? v1q2_idx_temp : -1")

        df = df.Define("excluded_jet_mask_for_v2q1", "ROOT::RVec<int>{truth_vbs1_idx, truth_vbs2_idx, truth_b1_idx, truth_b2_idx, truth_v1q1_idx, truth_v1q2_idx}") \
            .Define("matched_fatjet_mask_for_v2q", "ROOT::RVec<int>{truth_bh_idx, truth_bv1_idx, truth_bv2_idx}") \
            .Define("v2q1_idx_temp", "find_matching_jet(gen_v2q1_idx, V2_q1_eta, V2_q1_phi, excluded_jet_mask_for_v2q1, matched_fatjet_mask_for_v2q, ak4jet_eta, ak4jet_phi, ak8jet_eta, ak8jet_phi)") \
            .Define("excluded_jet_mask_for_v2q2", "ROOT::RVec<int>{truth_vbs1_idx, truth_vbs2_idx, truth_b1_idx, truth_b2_idx, truth_v1q1_idx, truth_v1q2_idx, v2q1_idx_temp}") \
            .Define("v2q2_idx_temp", "find_matching_jet(gen_v2q2_idx, V2_q2_eta, V2_q2_phi, excluded_jet_mask_for_v2q2, matched_fatjet_mask_for_v2q, ak4jet_eta, ak4jet_phi, ak8jet_eta, ak8jet_phi)") \
            .Define("truth_v2q1_idx", "v2q1_idx_temp >= 0 && v2q1_idx_temp < 10 ? v2q1_idx_temp : -1") \
            .Define("truth_v2q2_idx", "v2q2_idx_temp >= 0 && v2q2_idx_temp < 10 ? v2q2_idx_temp : -1")

        return df
        
    def analyze(self, is_signal):
        self.df = self.df.Define("tight_mu_mask", "Muon_pt > 35. && abs(Muon_eta) < 2.4 && Muon_tightId") \
            .Define("tight_ele_mask", "Electron_pt > 35. && abs(Electron_eta) < 2.5 && Electron_cutBased >= 4") \
            .Filter("Sum(tight_mu_mask) + Sum(tight_ele_mask) < 2") \
            .Define("fatjet_mask", "FatJet_pt > 200 && FatJet_msoftdrop > 10 && FatJet_mass > 10") \
            .Filter("Sum(fatjet_mask) >= 1")

        if self.sample_year in JET_ID_JSONS:
            jet_id_json = JET_ID_JSONS[self.sample_year]

        self.df = self.df.Define("Jet_multiplicity", "Jet_chMultiplicity + Jet_neMultiplicity") \
            .Define("FatJet_multiplicity", "FatJet_chMultiplicity + FatJet_neMultiplicity")

        r.gInterpreter.Declare("""
            #include <ROOT/RVec.hxx>
            using namespace ROOT::VecOps;
            
            RVec<float> evalJetID(const RVec<float>& eta, const RVec<float>& chHEF, const RVec<float>& neHEF,
                                  const RVec<float>& chEmEF, const RVec<float>& neEmEF,
                                  const RVec<float>& muEF, const RVec<int>& chMultiplicity,
                                  const RVec<int>& neMultiplicity, const RVec<int>& multiplicity) {
                auto cset_jetId = correction::CorrectionSet::from_file(\"""" + jet_id_json + """\");
                RVec<float> jetId(eta.size(), 0.0);
                for (size_t i = 0; i < eta.size(); ++i) {
                    jetId[i] += 2 * cset_jetId->at(\"AK4PUPPI_Tight\")->evaluate({eta[i], chHEF[i], neHEF[i], chEmEF[i], neEmEF[i], muEF[i], chMultiplicity[i], neMultiplicity[i], multiplicity[i]});
                    jetId[i] += 4 * cset_jetId->at(\"AK4PUPPI_TightLeptonVeto\")->evaluate({eta[i], chHEF[i], neHEF[i], chEmEF[i], neEmEF[i], muEF[i], chMultiplicity[i], neMultiplicity[i], multiplicity[i]});
                }
                return jetId;
            }

            RVec<float> evalFatJetID(const RVec<float>& eta, const RVec<float>& chHEF, const RVec<float>& neHEF,
                                     const RVec<float>& chEmEF, const RVec<float>& neEmEF,
                                     const RVec<float>& muEF, const RVec<int>& chMultiplicity,
                                     const RVec<int>& neMultiplicity, const RVec<int>& multiplicity) {
                auto cset_fatJetId = correction::CorrectionSet::from_file(\"""" + jet_id_json + """\");
                RVec<float> fatJetId(eta.size(), 0.0);
                for (size_t i = 0; i < eta.size(); ++i) {
                    fatJetId[i] += 2 * cset_fatJetId->at(\"AK8PUPPI_Tight\")->evaluate({eta[i], chHEF[i], neHEF[i], chEmEF[i], neEmEF[i], muEF[i], chMultiplicity[i], neMultiplicity[i], multiplicity[i]});
                    fatJetId[i] += 4 * cset_fatJetId->at(\"AK8PUPPI_TightLeptonVeto\")->evaluate({eta[i], chHEF[i], neHEF[i], chEmEF[i], neEmEF[i], muEF[i], chMultiplicity[i], neMultiplicity[i], multiplicity[i]});
                }
                return fatJetId;
            }

        """)

        self.df = self.df.Define("Jet_jetId", "evalJetID(Jet_eta, Jet_chHEF, Jet_neHEF, Jet_chEmEF, Jet_neEmEF, Jet_muEF, Jet_chMultiplicity, Jet_neMultiplicity, Jet_multiplicity)") \
            .Define("FatJet_jetId", "evalFatJetID(FatJet_eta, FatJet_chHEF, FatJet_neHEF, FatJet_chEmEF, FatJet_neEmEF, FatJet_muEF, FatJet_chMultiplicity, FatJet_neMultiplicity, FatJet_multiplicity)")

        if is_signal:
            self.df = self.genSelection(self.df)

        return self.df.Count().GetValue()

    def Snapshot(self, tag):
        all_cols = [str(col) for col in self.df.GetColumnNames()]
        keep_cols = {col: 0 for col in all_cols}
        comment = re.compile(r"#.*")
        ops = []
        with open(self.keepDropFile, 'r') as f:
            for line in f:
                # convert to python regex
                if len(line) == 0 or line[0] == '#': 
                    continue
                line = re.sub(comment, "", line)
                (op, sel) = line.split()
                if op == "keep":
                    ops.append((sel, 1))
                elif op == "drop":
                    ops.append((sel, 0))
        
        for bre, stat in ops:
            try:
                re.compile(bre)
                for n in all_cols:
                    if re.match(bre, n):
                        keep_cols[n] = stat
            except re.error:
                keep_cols[bre] = stat

        keep_cols = [k for k, v in keep_cols.items() if v == 1]

        self.df.Snapshot("Events", self.outDir + "/" + tag + ".root", keep_cols)
        
        snap_opts = r.RDF.RSnapshotOptions()
        snap_opts.fMode = "UPDATE"

        runs_df = r.RDataFrame("Runs", self.inFiles)
        cols = [str(col) for col in runs_df.GetColumnNames()]
        runs_df.Snapshot("Runs", self.outDir + "/" + tag + ".root", cols, snap_opts)

    @property
    def sample_year(self):
        match = re.search(r'Run3Summer24|RunIII2024Summer24NanoAODv15', self.inFiles[0])
        if match:
            return "2024"
        else:
            return None

def run_skimmer(input_file, output_dir, is_signal):
    print(f"Running skimmer on {input_file}")
    os.makedirs(output_dir, exist_ok=True)
    
    inFiles = [XROOTD_REDIRECTOR + input_file if input_file.startswith('/store') else 'file://' + input_file]
    keepDropFile = "keep_and_drop_skim.txt"
    
    skimmer = Skimmer(inFiles, output_dir, keepDropFile)
    passed = skimmer.analyze(is_signal)
    if passed:
        skimmer.Snapshot("skim")
        return True
    else:
        print("No entries in output")
        return False


def merge_skims(output_dir):
    skim_files = glob.glob(f"{output_dir}/*")
    
    if len(skim_files) == 0:
        print("No output files to merge; exiting...")
        return True
    elif len(skim_files) == 1:
        shutil.move(skim_files[0], f"{output_dir}/output.root")
        return True
    else:
        merge_cmd = ["hadd", f"{output_dir}/output.root"] + skim_files
        print(" ".join(merge_cmd))
        result = subprocess.run(merge_cmd)
        return result.returncode == 0


def determine_output_paths(input_file, is_signal):
    if not is_signal:
        era = input_file.split('/')[3]
        sample_name = input_file.split('/')[4]
        campaign = input_file.split('/')[6]
    else:
        era = input_file.split('/')[6]
        sample_name = input_file.split('/')[7]
        campaign = "private"
        
    output_dir = f"{OUTPUT_XRD}/{era}/{campaign}/{sample_name}"
    return output_dir

def copy_output_file(source, destination):
    print(f"Copying {source} to {destination}")
    
    # Create destination directory
    subprocess.run(["gfal-mkdir", "-p", os.path.dirname(destination)])
    
    # Copy with retries
    for i in range(1, MAX_RETRIES + 1):
        print(f"Attempt {i}")
        result = subprocess.run(["gfal-copy", "-f", source, destination])
        if result.returncode == 0:
            return True
        
        print(f"Failed to copy {source} to {destination}; sleeping for 60s")
        time.sleep(SLEEP_DURATION)
    
    return False

if __name__ == "__main__":
    parser = ArgumentParser(description='Run the NanoAOD skimmer with file transfer.')
    parser.add_argument('proxy', help="Path to the X509 proxy")
    parser.add_argument('input_file', help="Input file path")
    parser.add_argument('job_id', help="Job ID")
    parser.add_argument('is_signal', help='Flag indicating if this is a signal sample', type=int)
    args = parser.parse_args()
    
    # Set up X509 proxy
    os.environ['X509_USER_PROXY'] = args.proxy

    # Run the skimmer
    success = run_skimmer(args.input_file, CONDOR_OUTPUT_DIR, args.is_signal)
    
    # Retry once if failed
    if not success:
        print("Skimmer failed; retrying one more time...")
        success = run_skimmer(args.input_file, CONDOR_OUTPUT_DIR, args.is_signal)
    
    # Merge results
    merge_skims(CONDOR_OUTPUT_DIR)
    
    # Determine output paths
    output_dir = determine_output_paths(args.input_file, args.is_signal)
    
    # Copy the output file
    copy_src = os.path.join(os.getcwd(), f"{CONDOR_OUTPUT_DIR}/output.root")
    copy_dest = f"{output_dir}/output_{args.job_id}.root"
    
    success = copy_output_file(copy_src, copy_dest)
    if not success:
        print(f"Failed to copy output file after {MAX_RETRIES} attempts")
        sys.exit(1)
    
    sys.exit(0)
