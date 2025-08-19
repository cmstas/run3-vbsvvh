from argparse import ArgumentParser

import ROOT as r
r.EnableImplicitMT(64)

import h5py

import awkward as ak
import numpy as np

N_MAX_JETS = 10
N_MAX_FATJETS = 3

BTAGGING_SCORES = {
    "2022": {
        "Loose": 0.0583,
        "Medium": 0.3086,
        "Tight": 0.7183
    },
    "2024": {
        "Loose": 0.0485,
        "Medium": 0.2480,
        "Tight": 0.6708 
    }
}

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

def pad_or_truncate_bool(array, max_len):
    return ak.fill_none(
        ak.pad_none(array, max_len, clip=True), False
    )

def pad_or_truncate(array, max_len):
    return ak.fill_none(
        ak.pad_none(array, max_len, clip=True), 0
    )

def get_array_names(isSignal=False):
    """Return a list of array names needed for the dataset."""
    arrays = [
        # AK4 jets
        "Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass", "Jet_btagDeepFlavB",
        # AK8 jets
        "FatJet_pt", "FatJet_eta", "FatJet_phi", "FatJet_mass", "FatJet_nConstituents", 
        "FatJet_globalParT3_Xbb", "FatJet_globalParT3_Xcc", "FatJet_globalParT3_Xcs", "FatJet_globalParT3_Xqq", "FatJet_globalParT3_QCD",
        # Global
        "PuppiMET_pt", "PuppiMET_phi",
        # Event-level features
        "truth_hbb_idx", "truth_v1qq_idx", "truth_v2qq_idx",
        "truth_hbb1_idx", "truth_hbb2_idx", "truth_v1q1_idx", "truth_v1q2_idx", "truth_v2q1_idx", "truth_v2q2_idx",
        "truth_vbs1_idx", "truth_vbs2_idx"
    ]

    if isSignal:
        arrays.extend(["GenPart_pdgId"])
    
    return arrays

def calculate_met_features(data):
    """Add MET phi-related features."""
    data["MET_cosphi"] = np.cos(data["PuppiMET_phi"])
    data["MET_sinphi"] = np.sin(data["PuppiMET_phi"])
    return data

def calculate_jet_features(data):
    """Add AK4 jet features."""
    data["Jet_cosphi"] = np.cos(data["Jet_phi"])
    data["Jet_sinphi"] = np.sin(data["Jet_phi"])
    data["Jet_isLooseBTag"] = data["Jet_btagDeepFlavB"] > BTAGGING_SCORES["2024"]["Loose"]
    data["Jet_isMediumBTag"] = data["Jet_btagDeepFlavB"] > BTAGGING_SCORES["2024"]["Medium"]
    data["Jet_isTightBTag"] = data["Jet_btagDeepFlavB"] > BTAGGING_SCORES["2024"]["Tight"]
    data["Jet_MASK"] = ak.ones_like(data["Jet_pt"], dtype=bool)
    return data

def calculate_fatjet_features(data):
    """Add AK8 fatjet features."""
    data["FatJet_cosphi"] = np.cos(data["FatJet_phi"])
    data["FatJet_sinphi"] = np.sin(data["FatJet_phi"])
    data["FatJet_MASK"] = ak.ones_like(data["FatJet_pt"], dtype=bool)
    return data

def pad_jet_arrays(data, n_max_jets):
    """Pad or truncate AK4 jet arrays."""
    jet_arrays = [
        "pt", "eta", "mass", "cosphi", "sinphi"
    ]
    
    for arr in jet_arrays:
        data[f"Jet_{arr}"] = pad_or_truncate(data[f"Jet_{arr}"], n_max_jets)
    
    bool_arrays = ["isLooseBTag", "isMediumBTag", "isTightBTag", "MASK"]
    for arr in bool_arrays:
        data[f"Jet_{arr}"] = pad_or_truncate_bool(data[f"Jet_{arr}"], n_max_jets)
    
    return data

def pad_fatjet_arrays(data, n_max_fatjets):
    """Pad or truncate AK8 fatjet arrays."""
    fatjet_arrays = [
        "pt", "eta", "mass", "cosphi", "sinphi", "nConstituents",
        "globalParT3_Xbb", "globalParT3_Xcc", "globalParT3_Xcs", "globalParT3_Xqq", "globalParT3_QCD"
    ]
    
    for arr in fatjet_arrays:
        data[f"FatJet_{arr}"] = pad_or_truncate(data[f"FatJet_{arr}"], n_max_fatjets)
    
    data["FatJet_MASK"] = pad_or_truncate_bool(data["FatJet_MASK"], n_max_fatjets)
    
    return data

def prepare_dataset(input_dataframe, isSignal=False):
    arrays = get_array_names(isSignal)
    
    data = ak.from_rdataframe(input_dataframe, arrays)

    data = calculate_met_features(data)
    data = calculate_jet_features(data)
    data = calculate_fatjet_features(data)
    
    data = pad_jet_arrays(data, N_MAX_JETS)
    data = pad_fatjet_arrays(data, N_MAX_FATJETS)

    return data

def make_dataset(dataframes, output_path):
    sig_df, bkg_df = dataframes

    data_sig = prepare_dataset(sig_df, isSignal=True)
    data_sig["isSignal"] = True

    if bkg_df is not None:
        data_bkg = prepare_dataset(bkg_df)
        data_bkg["isSignal"] = False

        data = ak.concatenate([data_sig, data_bkg], axis=0)
    else:
        data = data_sig

    datasets = {
        # AK4 Jets
        "INPUTS/AK4Jets/MASK": data["Jet_MASK"],
        "INPUTS/AK4Jets/pt": data["Jet_pt"],
        "INPUTS/AK4Jets/eta": data["Jet_eta"],
        "INPUTS/AK4Jets/sinphi": data["Jet_sinphi"],
        "INPUTS/AK4Jets/cosphi": data["Jet_cosphi"],
        "INPUTS/AK4Jets/mass": data["Jet_mass"],
        "INPUTS/AK4Jets/isLooseBTag": data["Jet_isLooseBTag"],
        "INPUTS/AK4Jets/isMediumBTag": data["Jet_isMediumBTag"],
        "INPUTS/AK4Jets/isTightBTag": data["Jet_isTightBTag"],
        
        # AK8 Jets
        "INPUTS/AK8Jets/MASK": data["FatJet_MASK"],
        "INPUTS/AK8Jets/pt": data["FatJet_pt"],
        "INPUTS/AK8Jets/eta": data["FatJet_eta"],
        "INPUTS/AK8Jets/sinphi": data["FatJet_sinphi"],
        "INPUTS/AK8Jets/cosphi": data["FatJet_cosphi"],
        "INPUTS/AK8Jets/mass": data["FatJet_mass"],
        "INPUTS/AK8Jets/nConstituents": data["FatJet_nConstituents"],
        "INPUTS/AK8Jets/parT_Xbb": data["FatJet_globalParT3_Xbb"],
        "INPUTS/AK8Jets/parT_Xqq": data["FatJet_globalParT3_Xqq"],
        "INPUTS/AK8Jets/parT_Xcc": data["FatJet_globalParT3_Xcc"],
        "INPUTS/AK8Jets/parT_Xcs": data["FatJet_globalParT3_Xcs"],
        "INPUTS/AK8Jets/parT_QCD": data["FatJet_globalParT3_QCD"],
        
        # MET
        "INPUTS/MET/pt": data["PuppiMET_pt"],
        "INPUTS/MET/sinphi": data["MET_sinphi"],
        "INPUTS/MET/cosphi": data["MET_cosphi"],

        # Targets - Higgs
        "TARGETS/h/b1": data["truth_hbb1_idx"],
        "TARGETS/h/b2": data["truth_hbb2_idx"],
        "TARGETS/bh/bb": data["truth_hbb_idx"],
        
        # Targets - Vector boson
        "TARGETS/v1/v1q1": data["truth_v1q1_idx"],
        "TARGETS/v1/v1q2": data["truth_v1q2_idx"],
        "TARGETS/v2/v2q1": data["truth_v2q1_idx"],
        "TARGETS/v2/v2q2": data["truth_v2q2_idx"],
        "TARGETS/bv1/qq1": data["truth_v1qq_idx"],
        "TARGETS/bv2/qq2": data["truth_v2qq_idx"],
        
        # Targets - VBS
        "TARGETS/vbs/vbsq1": data["truth_vbs1_idx"],
        "TARGETS/vbs/vbsq2": data["truth_vbs2_idx"],
        
        # Classification
        "CLASSIFICATIONS/EVENT/isSignal": data["isSignal"],
    }

    with h5py.File(output_path, "w") as output:
        for dataset_name, dataset_values in datasets.items():
            output.create_dataset(dataset_name, data=dataset_values.to_numpy())

def genMatching(input_file, isSignal):
    df = r.RDF.Experimental.FromSpec(input_file)
    r.RDF.Experimental.AddProgressBar(df)

    # ak4 jet preselection
    df = df.Define("jets", "Jet_pt > 30 && abs(Jet_eta) <= 5 && Jet_jetId > 0") \
        .Redefine("Jet_pt", "Jet_pt[jets]") \
        .Redefine("Jet_eta", "Jet_eta[jets]") \
        .Redefine("Jet_phi", "Jet_phi[jets]") \
        .Redefine("Jet_mass", "Jet_mass[jets]") \
        .Redefine("Jet_btagDeepFlavB", "Jet_btagDeepFlavB[jets]")

    # fatjet preselection
    df = df.Define("ak8_jets", "FatJet_pt > 250 && abs(FatJet_eta) <= 2.5 && FatJet_mass > 50 && FatJet_msoftdrop > 40 && FatJet_jetId > 0") \
        .Redefine("FatJet_pt", "FatJet_pt[ak8_jets]") \
        .Redefine("FatJet_eta", "FatJet_eta[ak8_jets]") \
        .Redefine("FatJet_phi", "FatJet_phi[ak8_jets]") \
        .Redefine("FatJet_mass", "FatJet_mass[ak8_jets]") \
        .Redefine("FatJet_nConstituents", "FatJet_nConstituents[ak8_jets]") \
        .Redefine("FatJet_globalParT3_Xbb", "FatJet_globalParT3_Xbb[ak8_jets]") \
        .Redefine("FatJet_globalParT3_Xcc", "FatJet_globalParT3_Xcc[ak8_jets]") \
        .Redefine("FatJet_globalParT3_Xcs", "FatJet_globalParT3_Xcs[ak8_jets]") \
        .Redefine("FatJet_globalParT3_Xqq", "FatJet_globalParT3_Xqq[ak8_jets]") \
        .Redefine("FatJet_globalParT3_QCD", "FatJet_globalParT3_QCD[ak8_jets]")
    
    # define bb, qq definition
    if (isSignal):
        df = df.Define("higgs_info", "findHiggsAndDaughters(GenPart_pdgId, GenPart_status, GenPart_genPartIdxMother)") \
            .Define("Higgs_idx", "higgs_info[0]") \
            .Define("Higgs_eta", "Higgs_idx >= 0 ? GenPart_eta[Higgs_idx] : -999.0") \
            .Define("Higgs_phi", "Higgs_idx >= 0 ? GenPart_phi[Higgs_idx] : -999.0") \
            .Define("b1_idx", "higgs_info[1]") \
            .Define("b1_eta", "b1_idx >= 0 ? GenPart_eta[b1_idx] : -999.0") \
            .Define("b1_phi", "b1_idx >= 0 ? GenPart_phi[b1_idx] : -999.0") \
            .Define("b2_idx", "higgs_info[2]") \
            .Define("b2_eta", "b2_idx >= 0 ? GenPart_eta[b2_idx] : -999.0") \
            .Define("b2_phi", "b2_idx >= 0 ? GenPart_phi[b2_idx] : -999.0")

        df = df.Define("vboson_info", "findVBosonsAndDaughters(GenPart_pdgId, GenPart_status, GenPart_genPartIdxMother)") \
            .Define("V1_idx", "vboson_info[0]") \
            .Define("V1_eta", "V1_idx >= 0 ? GenPart_eta[V1_idx] : -999.0") \
            .Define("V1_phi", "V1_idx >= 0 ? GenPart_phi[V1_idx] : -999.0") \
            .Define("V1_q1_idx", "vboson_info[1]") \
            .Define("V1_q1_eta", "V1_q1_idx >= 0 ? GenPart_eta[V1_q1_idx] : -999.0") \
            .Define("V1_q1_phi", "V1_q1_idx >= 0 ? GenPart_phi[V1_q1_idx] : -999.0") \
            .Define("V1_q2_idx", "vboson_info[2]") \
            .Define("V1_q2_eta", "V1_q2_idx >= 0 ? GenPart_eta[V1_q2_idx] : -999.0") \
            .Define("V1_q2_phi", "V1_q2_idx >= 0 ? GenPart_phi[V1_q2_idx] : -999.0") \
            .Define("V2_idx", "vboson_info[3]") \
            .Define("V2_eta", "V2_idx >= 0 ? GenPart_eta[V2_idx] : -999.0") \
            .Define("V2_phi", "V2_idx >= 0 ? GenPart_phi[V2_idx] : -999.0") \
            .Define("V2_q1_idx", "vboson_info[4]") \
            .Define("V2_q1_eta", "V2_q1_idx >= 0 ? GenPart_eta[V2_q1_idx] : -999.0") \
            .Define("V2_q1_phi", "V2_q1_idx >= 0 ? GenPart_phi[V2_q1_idx] : -999.0") \
            .Define("V2_q2_idx", "vboson_info[5]") \
            .Define("V2_q2_eta", "V2_q2_idx >= 0 ? GenPart_eta[V2_q2_idx] : -999.0") \
            .Define("V2_q2_phi", "V2_q2_idx >= 0 ? GenPart_phi[V2_q2_idx] : -999.0")

        df = df.Define("vbs_info", "findVBSQuarks(GenPart_pdgId, GenPart_status, GenPart_genPartIdxMother)") \
            .Define("VBS_idx1", "vbs_info[0]") \
            .Define("VBSJet1_eta", "vbs_info[0] >= 0 ? GenPart_eta[vbs_info[0]] : -999.0") \
            .Define("VBSJet1_phi", "vbs_info[0] >= 0 ? GenPart_phi[vbs_info[0]] : -999.0") \
            .Define("VBS_idx2", "vbs_info[1]") \
            .Define("VBSJet2_eta", "vbs_info[1] >= 0 ? GenPart_eta[vbs_info[1]] : -999.0") \
            .Define("VBSJet2_phi", "vbs_info[1] >= 0 ? GenPart_phi[vbs_info[1]] : -999.0")

        df = df.Filter("Jet_pt.size() > 0") \
            .Define("truth_vbs1_idx_temp", "find_matching_jet(VBS_idx1, VBSJet1_eta, VBSJet1_phi, ROOT::RVec<int>{}, ROOT::RVec<int>{}, Jet_eta, Jet_phi, FatJet_eta, FatJet_phi)") \
            .Define("truth_vbs2_idx_temp", "find_matching_jet(VBS_idx2, VBSJet2_eta, VBSJet2_phi, ROOT::RVec<int>{truth_vbs1_idx_temp}, ROOT::RVec<int>{}, Jet_eta, Jet_phi, FatJet_eta, FatJet_phi)") \
            .Define("truth_vbs1_idx", "truth_vbs1_idx_temp >= 0 && truth_vbs1_idx_temp < 10 ? truth_vbs1_idx_temp : -1") \
            .Define("truth_vbs2_idx", "truth_vbs2_idx_temp >= 0 && truth_vbs2_idx_temp < 10 ? truth_vbs2_idx_temp : -1")

        df = df.Filter("FatJet_pt.size() > 0") \
            .Define("hbb_dR", "ROOT::VecOps::DeltaR(b1_eta, b2_eta, b1_phi, b2_phi)") \
            .Define("hbb_fatjet_idx_temp", "find_matching_fatjet(Higgs_idx, Higgs_eta, Higgs_phi, ROOT::RVec<int>{truth_vbs1_idx, truth_vbs2_idx}, ROOT::RVec<int>{}, Jet_eta, Jet_phi, FatJet_eta, FatJet_phi)") \
            .Define("hbb_fatjet_candidate_b1_dR", "hbb_fatjet_idx_temp >= 0 && hbb_fatjet_idx_temp < FatJet_eta.size() ? ROOT::VecOps::DeltaR(FatJet_eta[hbb_fatjet_idx_temp], b1_eta, FatJet_phi[hbb_fatjet_idx_temp], b1_phi) : 999.0") \
            .Define("hbb_fatjet_candidate_b2_dR", "hbb_fatjet_idx_temp >= 0 && hbb_fatjet_idx_temp < FatJet_eta.size() ? ROOT::VecOps::DeltaR(FatJet_eta[hbb_fatjet_idx_temp], b2_eta, FatJet_phi[hbb_fatjet_idx_temp], b2_phi) : 999.0") \
            .Define("hbb_isBoosted", "hbb_fatjet_idx_temp != -1 && hbb_fatjet_idx_temp < 3 && hbb_dR < 0.8 && hbb_fatjet_candidate_b1_dR < 0.8 && hbb_fatjet_candidate_b2_dR < 0.8") \
            .Define("truth_hbb_idx", "hbb_isBoosted ? hbb_fatjet_idx_temp : -1")
        
        df = df.Define("excluded_jets_for_hbb", "ROOT::RVec<int>{truth_vbs1_idx, truth_vbs2_idx}") \
            .Define("matched_fatjets_for_hbb", "ROOT::RVec<int>{}") \
            .Define("Jet_hbb1_idx_temp", "find_matching_jet(b1_idx, b1_eta, b1_phi, excluded_jets_for_hbb, matched_fatjets_for_hbb, Jet_eta, Jet_phi, FatJet_eta, FatJet_phi)") \
            .Define("excluded_jets_for_hbb2", "ROOT::RVec<int>{truth_vbs1_idx, truth_vbs2_idx, Jet_hbb1_idx_temp}") \
            .Define("truth_hbb2_idx_temp", "find_matching_jet(b2_idx, b2_eta, b2_phi, excluded_jets_for_hbb2, matched_fatjets_for_hbb, Jet_eta, Jet_phi, FatJet_eta, FatJet_phi)") \
            .Define("truth_hbb1_idx", "Jet_hbb1_idx_temp >= 0 && Jet_hbb1_idx_temp < 10 ? Jet_hbb1_idx_temp : -1") \
            .Define("truth_hbb2_idx", "truth_hbb2_idx_temp >= 0 && truth_hbb2_idx_temp < 10 ? truth_hbb2_idx_temp : -1")

        df = df.Define("v1qq_dR", "V1_idx != -1 ? ROOT::VecOps::DeltaR(V1_q1_eta, V1_q1_phi, V1_q2_eta, V1_q2_phi) : 999.0") \
            .Define("excluded_fatjets_for_v1", "ROOT::RVec<int>{truth_hbb_idx}") \
            .Define("matched_jets_for_v1", "ROOT::RVec<int>{truth_vbs1_idx, truth_vbs2_idx, truth_hbb1_idx, truth_hbb2_idx}") \
            .Define("v1qq_fatjet_idx_temp", "find_matching_fatjet(V1_idx, V1_eta, V1_phi, matched_jets_for_v1, excluded_fatjets_for_v1, Jet_eta, Jet_phi, FatJet_eta, FatJet_phi)") \
            .Define("v1qq_fatjet_candidate_q1_dR", "v1qq_fatjet_idx_temp >= 0 && v1qq_fatjet_idx_temp < FatJet_eta.size() ? ROOT::VecOps::DeltaR(FatJet_eta[v1qq_fatjet_idx_temp], V1_q1_eta, FatJet_phi[v1qq_fatjet_idx_temp], V1_q1_phi) : 999.0") \
            .Define("v1qq_fatjet_candidate_q2_dR", "v1qq_fatjet_idx_temp >= 0 && v1qq_fatjet_idx_temp < FatJet_eta.size() ? ROOT::VecOps::DeltaR(FatJet_eta[v1qq_fatjet_idx_temp], V1_q2_eta, FatJet_phi[v1qq_fatjet_idx_temp], V1_q2_phi) : 999.0") \
            .Define("v1qq_isBoosted", "V1_idx != -1 && v1qq_fatjet_idx_temp != -1 && v1qq_fatjet_idx_temp < 3 && v1qq_dR < 0.8 && v1qq_fatjet_candidate_q1_dR < 0.8 && v1qq_fatjet_candidate_q2_dR < 0.8") \
            .Define("truth_v1qq_idx", "v1qq_isBoosted ? v1qq_fatjet_idx_temp : -1")

        df = df.Define("v2qq_dR", "V2_idx != -1 ? ROOT::VecOps::DeltaR(V2_q1_eta, V2_q1_phi, V2_q2_eta, V2_q2_phi) : 999.0") \
            .Define("excluded_fatjets_for_v2", "ROOT::RVec<int>{truth_hbb_idx, truth_v1qq_idx}") \
            .Define("matched_jets_for_v2", "ROOT::RVec<int>{truth_vbs1_idx, truth_vbs2_idx, truth_hbb1_idx, truth_hbb2_idx}") \
            .Define("v2qq_fatjet_idx_temp", "find_matching_fatjet(V2_idx, V2_eta, V2_phi, matched_jets_for_v2, excluded_fatjets_for_v2, Jet_eta, Jet_phi, FatJet_eta, FatJet_phi)") \
            .Define("v2qq_fatjet_candidate_q1_dR", "v2qq_fatjet_idx_temp >= 0 && v2qq_fatjet_idx_temp < FatJet_eta.size() ? ROOT::VecOps::DeltaR(FatJet_eta[v2qq_fatjet_idx_temp], V2_q1_eta, FatJet_phi[v2qq_fatjet_idx_temp], V2_q1_phi) : 999.0") \
            .Define("v2qq_fatjet_candidate_q2_dR", "v2qq_fatjet_idx_temp >= 0 && v2qq_fatjet_idx_temp < FatJet_eta.size() ? ROOT::VecOps::DeltaR(FatJet_eta[v2qq_fatjet_idx_temp], V2_q2_eta, FatJet_phi[v2qq_fatjet_idx_temp], V2_q2_phi) : 999.0") \
            .Define("v2qq_isBoosted", "V2_idx != -1 && v2qq_fatjet_idx_temp != -1 && v2qq_fatjet_idx_temp < 3 && v2qq_dR < 0.8 && v2qq_fatjet_candidate_q1_dR < 0.8 && v2qq_fatjet_candidate_q2_dR < 0.8") \
            .Define("truth_v2qq_idx", "v2qq_isBoosted ? v2qq_fatjet_idx_temp : -1")
        
        df = df.Define("excluded_jets_for_v1q1", "ROOT::RVec<int>{truth_vbs1_idx, truth_vbs2_idx, truth_hbb1_idx, truth_hbb2_idx}") \
            .Define("matched_fatjets_for_v1q", "ROOT::RVec<int>{truth_hbb_idx, truth_v1qq_idx, truth_v2qq_idx}") \
            .Define("truth_v1q1_idx_temp", "find_matching_jet(V1_q1_idx, V1_q1_eta, V1_q1_phi, excluded_jets_for_v1q1, matched_fatjets_for_v1q, Jet_eta, Jet_phi, FatJet_eta, FatJet_phi)") \
            .Define("excluded_jets_for_v1q2", "ROOT::RVec<int>{truth_vbs1_idx, truth_vbs2_idx, truth_hbb1_idx, truth_hbb2_idx, truth_v1q1_idx_temp}") \
            .Define("truth_v1q2_idx_temp", "find_matching_jet(V1_q2_idx, V1_q2_eta, V1_q2_phi, excluded_jets_for_v1q2, matched_fatjets_for_v1q, Jet_eta, Jet_phi, FatJet_eta, FatJet_phi)") \
            .Define("truth_v1q1_idx", "truth_v1q1_idx_temp >= 0 && truth_v1q1_idx_temp < 10 ? truth_v1q1_idx_temp : -1") \
            .Define("truth_v1q2_idx", "truth_v1q2_idx_temp >= 0 && truth_v1q2_idx_temp < 10 ? truth_v1q2_idx_temp : -1")

        df = df.Define("excluded_jets_for_v2q1", "ROOT::RVec<int>{truth_vbs1_idx, truth_vbs2_idx, truth_hbb1_idx, truth_hbb2_idx, truth_v1q1_idx, truth_v1q2_idx}") \
            .Define("matched_fatjets_for_v2q", "ROOT::RVec<int>{truth_hbb_idx, truth_v1qq_idx, truth_v2qq_idx}") \
            .Define("truth_v2q1_idx_temp", "find_matching_jet(V2_q1_idx, V2_q1_eta, V2_q1_phi, excluded_jets_for_v2q1, matched_fatjets_for_v2q, Jet_eta, Jet_phi, FatJet_eta, FatJet_phi)") \
            .Define("excluded_jets_for_v2q2", "ROOT::RVec<int>{truth_vbs1_idx, truth_vbs2_idx, truth_hbb1_idx, truth_hbb2_idx, truth_v1q1_idx, truth_v1q2_idx, truth_v2q1_idx_temp}") \
            .Define("truth_v2q2_idx_temp", "find_matching_jet(V2_q2_idx, V2_q2_eta, V2_q2_phi, excluded_jets_for_v2q2, matched_fatjets_for_v2q, Jet_eta, Jet_phi, FatJet_eta, FatJet_phi)") \
            .Define("truth_v2q1_idx", "truth_v2q1_idx_temp >= 0 && truth_v2q1_idx_temp < 10 ? truth_v2q1_idx_temp : -1") \
            .Define("truth_v2q2_idx", "truth_v2q2_idx_temp >= 0 && truth_v2q2_idx_temp < 10 ? truth_v2q2_idx_temp : -1")

    else:
        df = df.Define("truth_hbb_idx", "-1") \
            .Define("truth_hbb1_idx", "-1") \
            .Define("truth_hbb2_idx", "-1") \
            .Define("truth_v1qq_idx", "-1") \
            .Define("truth_v2qq_idx", "-1") \
            .Define("truth_v1q1_idx", "-1") \
            .Define("truth_v1q2_idx", "-1") \
            .Define("truth_v2q1_idx", "-1") \
            .Define("truth_v2q2_idx", "-1") \
            .Define("truth_vbs1_idx", "-1") \
            .Define("truth_vbs2_idx", "-1") \
            .Define("hbb_isBoosted", "false") \
            .Define("v1qq_isBoosted", "false") \
            .Define("v2qq_isBoosted", "false") \
            .Filter("event % 5 == 0")
    return df

if __name__ == "__main__":
    df_bkg = None
    df_sig = genMatching(f"/home/users/aaarora/phys/run3/run3-vbsvvh/preselection/etc/1Lep2FJ-sig.json", isSignal=True)
    # df_bkg = genMatching(f"/home/users/aaarora/phys/run3/SPANet/input_proc/1Lep2FJ1.5-bkg.json", isSignal=False)

    df_sig.Snapshot("Events", "out.root")
    print("Number of evts in sig: ", df_sig.Count().GetValue())
    # print("Number of evts in bkg: ", df_bkg.Count().GetValue())
    
    make_dataset((df_sig, df_bkg), output_path=f"../data/output.h5")