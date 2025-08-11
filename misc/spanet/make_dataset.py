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
ROOT::RVec<float> get_dR(float eta1, float phi1, ROOT::RVec<float> eta2, ROOT::RVec<float> phi2) {
    ROOT::RVec<float> dR;
    for (size_t i = 0; i < eta2.size(); ++i) {
        dR.push_back(ROOT::VecOps::DeltaR(eta1, eta2[i], phi1, phi2[i]));
    }
    return dR;
}

int find_matching_jet(ROOT::RVec<float> dR_values, ROOT::RVec<int> excluded_indices) {
    int max_jets = 10;
    const float dR_cut = 0.4f;
    auto sorted_indices = ROOT::VecOps::Argsort(dR_values);
    for (int idx : sorted_indices) {
        if (dR_values[idx] >= dR_cut) break; // No more candidates within dR cut
        if (idx >= max_jets) continue; // Skip indices beyond padding limit
        bool is_excluded = false;
        for (int excl_idx : excluded_indices) {
            if (idx == excl_idx) {
                is_excluded = true;
                break;
            }
        }
        if (!is_excluded) {
            return idx;
        }
    }
    return -1;
}

int find_matching_fatjet(ROOT::RVec<float> dR_values, ROOT::RVec<int> excluded_indices) {
    int max_fatjets = 3;
    const float dR_cut = 0.8f;
    auto sorted_indices = ROOT::VecOps::Argsort(dR_values);
    for (int idx : sorted_indices) {
        if (dR_values[idx] >= dR_cut) break; // No more candidates within dR cut
        if (idx >= max_fatjets) continue; // Skip indices beyond padding limit
        bool is_excluded = false;
        for (int excl_idx : excluded_indices) {
            if (idx == excl_idx) {
                is_excluded = true;
                break;
            }
        }
        if (!is_excluded) {
            return idx;
        }
    }
    return -1;
}
                       
int num_hadronic_gauge_bosons(ROOT::RVec<int> pdgId, ROOT::RVec<short> motherIdx) {
    // particles at idx 2 and 3 are hadronic gauge bosons
    auto daughterIdx = [&] (int idx) {
        ROOT::RVec<int> daughters;
        for (size_t i = 0; i < pdgId.size(); ++i) {
            if (motherIdx[i] == idx) {
                daughters.push_back(i);
            }
        }
        return daughters;
    };
                       
    std::function<bool(int)> is_hadronic = [&] (int idx) -> bool {
        if (idx < 0 || idx >= pdgId.size()) return false;
        if (pdgId[idx] == 23 || pdgId[idx] == 24 || pdgId[idx] == -24) {
            auto daughters = daughterIdx(idx);
            for (int daughter : daughters) {
                if (pdgId[daughter] == 23 || pdgId[daughter] == 24 || pdgId[daughter] == -24) {
                    if (is_hadronic(daughter)) {
                        return true;
                    }
                } else if (abs(pdgId[daughter]) < 6) {
                    return true;
                } else if (abs(pdgId[daughter]) == 11 || abs(pdgId[daughter]) == 13 || abs(pdgId[daughter]) == 15) {
                    return false;
                }
            }
        }
        return false;
    };

    // return 0 if both are leptonic, 1 if only genpart 2 is hadronic, 2 if only genpart 3 is hadronic, 3 if both are hadronic
    int count = 0;
    if (is_hadronic(2)) count += 1;
    if (is_hadronic(3)) count += 2;
    return count;
}

int get_higgs_boson_idx(ROOT::RVec<int>& pdgId, ROOT::RVec<short>& motherIdx) {
    const int bId = 5;
    const int hId = 25;
    ROOT::RVec<int> bIndices, bBarIndices;

    for (size_t i = 0; i < pdgId.size(); ++i) {
        if (pdgId[i] == bId) bIndices.push_back(i);
        else if (pdgId[i] == -bId) bBarIndices.push_back(i);
    }

    auto traceToTopHiggs = [&](int idx) {
        while (idx >= 0 && pdgId[idx] == hId && motherIdx[idx] >= 0 && pdgId[motherIdx[idx]] == hId) {
            idx = motherIdx[idx];
        }
        return (pdgId[idx] == hId) ? idx : -1;
    };

    for (int bIdx : bIndices) {
        for (int bBarIdx : bBarIndices) {
            int motherB = motherIdx[bIdx];
            int motherBBar = motherIdx[bBarIdx];

            if (motherB == -1 || motherBBar == -1) continue;
            if (motherB != motherBBar) continue;
            if (pdgId[motherB] != hId) continue;

            int topHiggsIdx = traceToTopHiggs(motherB);
            if (topHiggsIdx == 4) {
                return motherB;
            }
        }
    }

    return -1;
}                    
                       
"""
)

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
        "FatJet_hbb_idx", "FatJet_v1qq_idx", "FatJet_v2qq_idx",
        "Jet_hbb1_idx", "Jet_hbb2_idx", "Jet_v1q1_idx", "Jet_v1q2_idx", "Jet_v2q1_idx", "Jet_v2q2_idx",
        "Jet_vbs1_idx", "Jet_vbs2_idx",
        "hbb_isBoosted", "hbb_isResolved", "v1qq_isBoosted", "v1qq_isResolved", "v2qq_isBoosted", "v2qq_isResolved",
    ]

    if isSignal:
        arrays.extend(["GenPart_pdgId", "GenPart_motherPdgId"])
    
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
        "TARGETS/h/b1": data["Jet_hbb1_idx"],
        "TARGETS/h/b2": data["Jet_hbb2_idx"],
        "TARGETS/bh/bb": data["FatJet_hbb_idx"],
        
        # Targets - Vector boson
        "TARGETS/v1/v1q1": data["Jet_v1q1_idx"],
        "TARGETS/v1/v1q2": data["Jet_v1q2_idx"],
        "TARGETS/v2/v2q1": data["Jet_v2q1_idx"],
        "TARGETS/v2/v2q2": data["Jet_v2q2_idx"],
        "TARGETS/bv1/qq1": data["FatJet_v1qq_idx"],
        "TARGETS/bv2/qq2": data["FatJet_v2qq_idx"],
        
        # Targets - VBS
        "TARGETS/vbs/vbsq1": data["Jet_vbs1_idx"],
        "TARGETS/vbs/vbsq2": data["Jet_vbs2_idx"],
        
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
        df = df.Define("GenPart_motherPdgId", "Take(GenPart_pdgId, GenPart_genPartIdxMother)")
        df = df.Define("Higgs_idx", "get_higgs_boson_idx(GenPart_pdgId, GenPart_genPartIdxMother)") \
            .Define("Higgs_eta", "GenPart_eta[Higgs_idx]") \
            .Define("Higgs_phi", "GenPart_phi[Higgs_idx]") \
            .Define("VBSJet1_eta", "GenPart_eta[5]") \
            .Define("VBSJet1_phi", "GenPart_phi[5]") \
            .Define("VBSJet2_eta", "GenPart_eta[6]") \
            .Define("VBSJet2_phi", "GenPart_phi[6]") \
            .Define("num_v", "num_hadronic_gauge_bosons(GenPart_pdgId, GenPart_genPartIdxMother)") \
            .Define("V1_idx", "num_v == 1 ? 2 : num_v == 2 ? 3 : num_v == 3 ? 2 : -1") \
            .Define("V2_idx", "num_v == 3 ? 3 : -1") \
            .Define("V1_eta", "V1_idx != -1 ? GenPart_eta[V1_idx] : -999") \
            .Define("V1_phi", "V1_idx != -1 ? GenPart_phi[V1_idx] : -999") \
            .Define("V2_eta", "V2_idx != -1 ? GenPart_eta[V2_idx] : -999") \
            .Define("V2_phi", "V2_idx != -1 ? GenPart_phi[V2_idx] : -999") \
            .Define("qs_from_v1", "abs(GenPart_pdgId) <= 5 && GenPart_genPartIdxMother == V1_idx") \
            .Define("qs_from_v2", "abs(GenPart_pdgId) <= 5 && GenPart_genPartIdxMother == V2_idx") \
            .Define("V1_q1_eta", "GenPart_eta[qs_from_v1][0]") \
            .Define("V1_q1_phi", "GenPart_phi[qs_from_v1][0]") \
            .Define("V1_q2_eta", "GenPart_eta[qs_from_v1][1]") \
            .Define("V1_q2_phi", "GenPart_phi[qs_from_v1][1]") \
            .Define("V2_q1_eta", "GenPart_eta[qs_from_v2][0]") \
            .Define("V2_q1_phi", "GenPart_phi[qs_from_v2][0]") \
            .Define("V2_q2_eta", "GenPart_eta[qs_from_v2][1]") \
            .Define("V2_q2_phi", "GenPart_phi[qs_from_v2][1]") \
            .Define("bs_from_higgs", "abs(GenPart_pdgId) == 5 && GenPart_genPartIdxMother == Higgs_idx") \
            .Define("b1_eta", "GenPart_eta[bs_from_higgs][0]") \
            .Define("b1_phi", "GenPart_phi[bs_from_higgs][0]") \
            .Define("b2_eta", "GenPart_eta[bs_from_higgs][1]") \
            .Define("b2_phi", "GenPart_phi[bs_from_higgs][1]")

        # VBS matching - Use argsort approach with bounds checking
        df = df.Filter("Jet_pt.size() > 0") \
            .Define("vbs1_dR", "get_dR(VBSJet1_eta, VBSJet1_phi, Jet_eta, Jet_phi)") \
            .Define("vbs2_dR", "get_dR(VBSJet2_eta, VBSJet2_phi, Jet_eta, Jet_phi)") \
            .Define("Jet_vbs1_idx_temp", "find_matching_jet(vbs1_dR, ROOT::RVec<int>{})") \
            .Define("Jet_vbs2_idx_temp", "find_matching_jet(vbs2_dR, ROOT::RVec<int>{Jet_vbs1_idx_temp})") \
            .Define("Jet_vbs1_idx", "Jet_vbs1_idx_temp >= 0 && Jet_vbs1_idx_temp < 10 ? Jet_vbs1_idx_temp : -1") \
            .Define("Jet_vbs2_idx", "Jet_vbs2_idx_temp >= 0 && Jet_vbs2_idx_temp < 10 ? Jet_vbs2_idx_temp : -1")

        # Boosted Hbb matching - Use argsort approach with bounds checking
        df = df.Filter("FatJet_pt.size() > 0") \
            .Define("hbb_fatjet_dR", "get_dR(Higgs_eta, Higgs_phi, FatJet_eta, FatJet_phi)") \
            .Define("hbb_dR", "ROOT::VecOps::DeltaR(b1_eta, b1_phi, b2_eta, b2_phi)") \
            .Define("hbb_fatjet_idx_temp", "Higgs_idx != -1 && hbb_fatjet_dR.size() > 0 ? find_matching_fatjet(hbb_fatjet_dR, ROOT::RVec<int>{}) : -1") \
            .Define("hbb_fatjet_candidate_b1_dR", "hbb_fatjet_idx_temp >= 0 && hbb_fatjet_idx_temp < FatJet_eta.size() ? ROOT::VecOps::DeltaR(FatJet_eta[hbb_fatjet_idx_temp], FatJet_phi[hbb_fatjet_idx_temp], b1_eta, b1_phi) : 999.0") \
            .Define("hbb_fatjet_candidate_b2_dR", "hbb_fatjet_idx_temp >= 0 && hbb_fatjet_idx_temp < FatJet_eta.size() ? ROOT::VecOps::DeltaR(FatJet_eta[hbb_fatjet_idx_temp], FatJet_phi[hbb_fatjet_idx_temp], b2_eta, b2_phi) : 999.0") \
            .Define("hbb_isBoosted", "hbb_fatjet_idx_temp != -1 && hbb_fatjet_idx_temp < 3 && hbb_dR < 0.8 && hbb_fatjet_candidate_b1_dR < 0.8 && hbb_fatjet_candidate_b2_dR < 0.8") \
            .Define("FatJet_hbb_idx", "hbb_isBoosted ? hbb_fatjet_idx_temp : -1")
        
        # Resolved Hbb matching - Use argsort approach with bounds checking
        df = df.Define("b1_jet_dR", "get_dR(b1_eta, b1_phi, Jet_eta, Jet_phi)") \
            .Define("b2_jet_dR", "get_dR(b2_eta, b2_phi, Jet_eta, Jet_phi)") \
            .Define("excluded_jets_for_hbb", "ROOT::RVec<int>{Jet_vbs1_idx, Jet_vbs2_idx}") \
            .Define("Jet_hbb1_idx_temp", "Higgs_idx != -1 && b1_jet_dR.size() > 0 ? find_matching_jet(b1_jet_dR, excluded_jets_for_hbb) : -1") \
            .Define("excluded_jets_for_hbb2", "ROOT::RVec<int>{Jet_vbs1_idx, Jet_vbs2_idx, Jet_hbb1_idx_temp}") \
            .Define("Jet_hbb2_idx_temp", "Higgs_idx != -1 && b2_jet_dR.size() > 0 ? find_matching_jet(b2_jet_dR, excluded_jets_for_hbb2) : -1") \
            .Define("Jet_hbb1_idx", "Jet_hbb1_idx_temp >= 0 && Jet_hbb1_idx_temp < 10 ? Jet_hbb1_idx_temp : -1") \
            .Define("Jet_hbb2_idx", "Jet_hbb2_idx_temp >= 0 && Jet_hbb2_idx_temp < 10 ? Jet_hbb2_idx_temp : -1") \
            .Define("hbb_isResolved", "!hbb_isBoosted && Jet_hbb1_idx != -1 && Jet_hbb2_idx != -1")

        # Boosted V1qq matching - Use argsort approach with bounds checking
        df = df.Define("v1qq_fatjet_dR", "V1_idx != -1 ? get_dR(V1_eta, V1_phi, FatJet_eta, FatJet_phi) : ROOT::RVec<float>()") \
            .Define("v1qq_dR", "V1_idx != -1 ? ROOT::VecOps::DeltaR(V1_q1_eta, V1_q1_phi, V1_q2_eta, V1_q2_phi) : 999.0") \
            .Define("excluded_fatjets_for_v1", "ROOT::RVec<int>{FatJet_hbb_idx}") \
            .Define("v1qq_fatjet_idx_temp", "V1_idx != -1 && v1qq_fatjet_dR.size() > 0 ? find_matching_fatjet(v1qq_fatjet_dR, excluded_fatjets_for_v1) : -1") \
            .Define("v1qq_fatjet_candidate_q1_dR", "v1qq_fatjet_idx_temp >= 0 && v1qq_fatjet_idx_temp < FatJet_eta.size() ? ROOT::VecOps::DeltaR(FatJet_eta[v1qq_fatjet_idx_temp], FatJet_phi[v1qq_fatjet_idx_temp], V1_q1_eta, V1_q1_phi) : 999.0") \
            .Define("v1qq_fatjet_candidate_q2_dR", "v1qq_fatjet_idx_temp >= 0 && v1qq_fatjet_idx_temp < FatJet_eta.size() ? ROOT::VecOps::DeltaR(FatJet_eta[v1qq_fatjet_idx_temp], FatJet_phi[v1qq_fatjet_idx_temp], V1_q2_eta, V1_q2_phi) : 999.0") \
            .Define("v1qq_isBoosted", "V1_idx != -1 && v1qq_fatjet_idx_temp != -1 && v1qq_fatjet_idx_temp < 3 && v1qq_dR < 0.8 && v1qq_fatjet_candidate_q1_dR < 0.8 && v1qq_fatjet_candidate_q2_dR < 0.8") \
            .Define("FatJet_v1qq_idx", "v1qq_isBoosted ? v1qq_fatjet_idx_temp : -1")

        # Boosted V2qq matching - Use argsort approach with bounds checking
        df = df.Define("v2qq_fatjet_dR", "V2_idx != -1 ? get_dR(V2_eta, V2_phi, FatJet_eta, FatJet_phi) : ROOT::RVec<float>()") \
            .Define("v2qq_dR", "V2_idx != -1 ? ROOT::VecOps::DeltaR(V2_q1_eta, V2_q1_phi, V2_q2_eta, V2_q2_phi) : 999.0") \
            .Define("excluded_fatjets_for_v2", "ROOT::RVec<int>{FatJet_hbb_idx, FatJet_v1qq_idx}") \
            .Define("v2qq_fatjet_idx_temp", "V2_idx != -1 && v2qq_fatjet_dR.size() > 0 ? find_matching_fatjet(v2qq_fatjet_dR, excluded_fatjets_for_v2) : -1") \
            .Define("v2qq_fatjet_candidate_q1_dR", "v2qq_fatjet_idx_temp >= 0 && v2qq_fatjet_idx_temp < FatJet_eta.size() ? ROOT::VecOps::DeltaR(FatJet_eta[v2qq_fatjet_idx_temp], FatJet_phi[v2qq_fatjet_idx_temp], V2_q1_eta, V2_q1_phi) : 999.0") \
            .Define("v2qq_fatjet_candidate_q2_dR", "v2qq_fatjet_idx_temp >= 0 && v2qq_fatjet_idx_temp < FatJet_eta.size() ? ROOT::VecOps::DeltaR(FatJet_eta[v2qq_fatjet_idx_temp], FatJet_phi[v2qq_fatjet_idx_temp], V2_q2_eta, V2_q2_phi) : 999.0") \
            .Define("v2qq_isBoosted", "V2_idx != -1 && v2qq_fatjet_idx_temp != -1 && v2qq_fatjet_idx_temp < 3 && v2qq_dR < 0.8 && v2qq_fatjet_candidate_q1_dR < 0.8 && v2qq_fatjet_candidate_q2_dR < 0.8") \
            .Define("FatJet_v2qq_idx", "v2qq_isBoosted ? v2qq_fatjet_idx_temp : -1")
        
        # Resolved V1qq matching - Use argsort approach with bounds checking
        df = df.Define("v1q1_jet_dR", "V1_idx != -1 ? get_dR(V1_q1_eta, V1_q1_phi, Jet_eta, Jet_phi) : ROOT::RVec<float>()") \
            .Define("v1q2_jet_dR", "V1_idx != -1 ? get_dR(V1_q2_eta, V1_q2_phi, Jet_eta, Jet_phi) : ROOT::RVec<float>()") \
            .Define("excluded_jets_for_v1q1", "ROOT::RVec<int>{Jet_vbs1_idx, Jet_vbs2_idx, Jet_hbb1_idx, Jet_hbb2_idx}") \
            .Define("Jet_v1q1_idx_temp", "V1_idx != -1 && v1q1_jet_dR.size() > 0 ? find_matching_jet(v1q1_jet_dR, excluded_jets_for_v1q1) : -1") \
            .Define("excluded_jets_for_v1q2", "ROOT::RVec<int>{Jet_vbs1_idx, Jet_vbs2_idx, Jet_hbb1_idx, Jet_hbb2_idx, Jet_v1q1_idx_temp}") \
            .Define("Jet_v1q2_idx_temp", "V1_idx != -1 && v1q2_jet_dR.size() > 0 ? find_matching_jet(v1q2_jet_dR, excluded_jets_for_v1q2) : -1") \
            .Define("Jet_v1q1_idx", "Jet_v1q1_idx_temp >= 0 && Jet_v1q1_idx_temp < 10 ? Jet_v1q1_idx_temp : -1") \
            .Define("Jet_v1q2_idx", "Jet_v1q2_idx_temp >= 0 && Jet_v1q2_idx_temp < 10 ? Jet_v1q2_idx_temp : -1") \
            .Define("v1qq_isResolved", "V1_idx != -1 && !v1qq_isBoosted && Jet_v1q1_idx != -1 && Jet_v1q2_idx != -1")

        # Resolved V2qq matching - Use argsort approach with bounds checking
        df = df.Define("v2q1_jet_dR", "V2_idx != -1 ? get_dR(V2_q1_eta, V2_q1_phi, Jet_eta, Jet_phi) : ROOT::RVec<float>()") \
            .Define("v2q2_jet_dR", "V2_idx != -1 ? get_dR(V2_q2_eta, V2_q2_phi, Jet_eta, Jet_phi) : ROOT::RVec<float>()") \
            .Define("excluded_jets_for_v2q1", "ROOT::RVec<int>{Jet_vbs1_idx, Jet_vbs2_idx, Jet_hbb1_idx, Jet_hbb2_idx, Jet_v1q1_idx, Jet_v1q2_idx}") \
            .Define("Jet_v2q1_idx_temp", "V2_idx != -1 && v2q1_jet_dR.size() > 0 ? find_matching_jet(v2q1_jet_dR, excluded_jets_for_v2q1) : -1") \
            .Define("excluded_jets_for_v2q2", "ROOT::RVec<int>{Jet_vbs1_idx, Jet_vbs2_idx, Jet_hbb1_idx, Jet_hbb2_idx, Jet_v1q1_idx, Jet_v1q2_idx, Jet_v2q1_idx_temp}") \
            .Define("Jet_v2q2_idx_temp", "V2_idx != -1 && v2q2_jet_dR.size() > 0 ? find_matching_jet(v2q2_jet_dR, excluded_jets_for_v2q2) : -1") \
            .Define("Jet_v2q1_idx", "Jet_v2q1_idx_temp >= 0 && Jet_v2q1_idx_temp < 10 ? Jet_v2q1_idx_temp : -1") \
            .Define("Jet_v2q2_idx", "Jet_v2q2_idx_temp >= 0 && Jet_v2q2_idx_temp < 10 ? Jet_v2q2_idx_temp : -1") \
            .Define("v2qq_isResolved", "V2_idx != -1 && !v2qq_isBoosted && Jet_v2q1_idx != -1 && Jet_v2q2_idx != -1")

    else:
        df = df.Define("FatJet_hbb_idx", "-1") \
            .Define("Jet_hbb1_idx", "-1") \
            .Define("Jet_hbb2_idx", "-1") \
            .Define("FatJet_v1qq_idx", "-1") \
            .Define("FatJet_v2qq_idx", "-1") \
            .Define("Jet_v1q1_idx", "-1") \
            .Define("Jet_v1q2_idx", "-1") \
            .Define("Jet_v2q1_idx", "-1") \
            .Define("Jet_v2q2_idx", "-1") \
            .Define("Jet_vbs1_idx", "-1") \
            .Define("Jet_vbs2_idx", "-1") \
            .Define("hbb_isBoosted", "false") \
            .Define("hbb_isResolved", "false") \
            .Define("v1qq_isBoosted", "false") \
            .Define("v1qq_isResolved", "false") \
            .Define("v2qq_isBoosted", "false") \
            .Define("v2qq_isResolved", "false") \
            .Filter("event % 5 == 0")
    return df

if __name__ == "__main__":
    df_bkg = None
    df_sig = genMatching(f"/home/users/aaarora/phys/run3/run3-vbsvvh/preselection/etc/1Lep2FJ-sig.json", isSignal=True)
    # df_bkg = genMatching(f"/home/users/aaarora/phys/run3/SPANet/input_proc/1Lep2FJ1.5-bkg.json", isSignal=False)

    # df_sig.Snapshot("Events", "out.root", get_array_names(isSignal=True))
    print("Number of evts in sig: ", df_sig.Count().GetValue())
    # print("Number of evts in bkg: ", df_bkg.Count().GetValue())
    
    make_dataset((df_sig, df_bkg), output_path=f"../data/output.h5")