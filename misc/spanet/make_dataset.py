from argparse import ArgumentParser

import ROOT as r
r.EnableImplicitMT(64)

import h5py

import uproot
import awkward as ak
import numpy as np

N_MAX_JETS = 10
N_MAX_FATJETS = 3

r.gInterpreter.Declare("""
ROOT::RVec<float> get_dR(float eta1, float phi1, ROOT::RVec<float> eta2, ROOT::RVec<float> phi2) {
    ROOT::RVec<float> dR;
    for (size_t i = 0; i < eta2.size(); ++i) {
        dR.push_back(ROOT::VecOps::DeltaR(eta1, eta2[i], phi1, phi2[i]));
    }
    return dR;
}
                       
int get_hadronic_gauge_boson_idx(ROOT::RVec<int> pdgId, ROOT::RVec<int> motherIdx) {
    int final_choice = -1;
    for (size_t i = 0; i < pdgId.size(); ++i) {
        if (abs(pdgId[i]) <= 5 && (abs(pdgId[motherIdx[i]]) == 24 || abs(pdgId[motherIdx[i]]) == 23)) {
            int current_idx = motherIdx[i];
            while (current_idx != 2 && current_idx != 3) { 
                if (abs(pdgId[motherIdx[current_idx]]) == 24 || abs(pdgId[motherIdx[current_idx]]) == 23) {
                    current_idx = motherIdx[current_idx];
                } else {
                    break;
                }
            }
            final_choice = current_idx;
            break;
        }
    }
    return final_choice;
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
        "FatJet_pt", "FatJet_eta", "FatJet_phi", "FatJet_mass", 
        "FatJet_nConstituents", "FatJet_particleNet_XbbVsQCD", "FatJet_particleNet_XqqVsQCD",
        # Global
        "MET_pt", "MET_phi",
        # Leptons
        "lepton_pt", "lepton_eta", "lepton_phi", "lepton_mass",
        # Event-level features
        "FatJet_hbb_idx", "FatJet_vqq_idx",
        "Jet_hbb1_idx", "Jet_hbb2_idx", "Jet_vqq1_idx", "Jet_vqq2_idx",
        "Jet_vbs1_idx", "Jet_vbs2_idx",
        "hbb_isBoosted", "hbb_isResolved", "vqq_isBoosted", "vqq_isResolved",
    ]

    if isSignal:
        arrays.extend(["GenPart_pdgId", "GenPart_motherPdgId"])
    
    return arrays

def filter_signal_events(data):
    """Filter signal events that have exactly one lepton from W decay."""
    lepton_mask = (
        (abs(data.GenPart_pdgId) == 11) | (abs(data.GenPart_pdgId) == 13)
    ) & (abs(data.GenPart_motherPdgId) == 24)

    event_mask = ak.sum(ak.fill_none(lepton_mask, False), axis=1) == 1
    return data[event_mask]

def calculate_met_features(data):
    """Add MET phi-related features."""
    data["MET_cosphi"] = np.cos(data["MET_phi"])
    data["MET_sinphi"] = np.sin(data["MET_phi"])
    return data

def calculate_lepton_features(data):
    """Add lepton phi-related features."""
    data["lepton_cosphi"] = np.cos(data["lepton_phi"])
    data["lepton_sinphi"] = np.sin(data["lepton_phi"])
    return data

def calculate_jet_features(data):
    """Add AK4 jet features."""
    data["Jet_cosphi"] = np.cos(data["Jet_phi"])
    data["Jet_sinphi"] = np.sin(data["Jet_phi"])
    data["Jet_isLooseBTag"] = data["Jet_btagDeepFlavB"] > 0.0583
    data["Jet_isMediumBTag"] = data["Jet_btagDeepFlavB"] > 0.3086
    data["Jet_isTightBTag"] = data["Jet_btagDeepFlavB"] > 0.7183
    data["Jet_MASK"] = ak.ones_like(data["Jet_pt"], dtype=bool)
    return data

def calculate_fatjet_features(data):
    """Add AK8 fatjet features."""
    data["FatJet_cosphi"] = np.cos(data["FatJet_phi"])
    data["FatJet_sinphi"] = np.sin(data["FatJet_phi"])
    data["FatJet_MASK"] = ak.ones_like(data["FatJet_pt"], dtype=bool)
    return data

def electronSelections(df):
    return df.Define("Electron_SC_eta", "Electron_eta + Electron_deltaEtaSC") \
        .Define("_looseElectrons", 
            "Electron_pt > 7 &&"
            "abs(Electron_SC_eta) < 2.5 && "
            "((abs(Electron_SC_eta) <= 1.479 && abs(Electron_dxy) <= 0.05 && abs(Electron_dz) < 0.1) || (abs(Electron_dxy) <= 0.1 && abs(Electron_dz) < 0.2)) && "
            "abs(Electron_sip3d) < 8 && "
            "Electron_cutBased >= 2 && "
            "Electron_pfRelIso03_all < 0.4 && "
            "Electron_lostHits <= 1") \
        .Define("nElectron_Loose", "nElectron == 0 ? 0 : Sum(_looseElectrons)") \
        .Define("_tightElectrons", "_looseElectrons &&" 
            "Electron_pt > 30 && "
            "Electron_cutBased >= 4 && "
            "Electron_pfRelIso03_all < 0.15 && "
            "Electron_hoe < 0.1 && "
            "Electron_eInvMinusPInv > -0.04 && "
            "((abs(Electron_SC_eta) <= 1.479 && Electron_sieie < 0.011) || Electron_sieie <= 0.030) && "
            "Electron_convVeto == true && "
            "Electron_tightCharge == 2 && "
            "Electron_lostHits == 0") \
        .Define("nElectron_Tight", "nElectron_Loose == 0 ? 0 : Sum(_tightElectrons)") \
        .Redefine("Electron_pt", "Electron_pt[_tightElectrons]") \
        .Redefine("Electron_eta", "Electron_eta[_tightElectrons]") \
        .Redefine("Electron_SC_eta", "Electron_SC_eta[_tightElectrons]") \
        .Redefine("Electron_phi", "Electron_phi[_tightElectrons]") \
        .Redefine("Electron_mass", "Electron_mass[_tightElectrons]")

def muonSelections(df):
    return df.Define("_looseMuons", 
            "Muon_pt > 5 && "
            "Muon_pfIsoId >= 2 && "
            "abs(Muon_eta) < 2.4 && "
            "abs(Muon_dxy) < 0.2 && "
            "abs(Muon_dz) < 0.5 && "
            "abs(Muon_sip3d) < 8 && "
            "Muon_looseId == 1") \
        .Define("nMuon_Loose", "nMuon == 0 ? 0 : Sum(_looseMuons)") \
        .Define("_tightMuons", "_looseMuons && "
            "Muon_pt > 30 && "
            "Muon_pfIsoId > 4 && "
            "Muon_tightCharge == 2 && "
            "Muon_highPurity && "
            "Muon_tightId") \
        .Define("nMuon_Tight", "nMuon_Loose == 0 ? 0 : Sum(_tightMuons)") \
        .Redefine("Muon_pt", "Muon_pt[_tightMuons]") \
        .Redefine("Muon_eta", "Muon_eta[_tightMuons]") \
        .Redefine("Muon_phi", "Muon_phi[_tightMuons]") \
        .Redefine("Muon_mass", "Muon_mass[_tightMuons]")
    
def leptonSelections(df):
    df = electronSelections(df)
    df = muonSelections(df)
    return df.Define("_isLepton", "(nElectron_Loose + nMuon_Loose == 1) && (nElectron_Tight + nMuon_Tight == 1)") \
        .Define("_isElectron", "nElectron_Loose == 1 && nMuon_Loose == 0 && nElectron_Tight == 1 && nMuon_Tight == 0") \
        .Define("_isMuon", "nElectron_Loose == 0 && nMuon_Loose == 1 && nElectron_Tight == 0 && nMuon_Tight == 1") \
        .Define("lepton_pt", "_isElectron ? Electron_pt[0] : Muon_pt[0]") \
        .Define("lepton_eta", "_isElectron ? Electron_eta[0] : Muon_eta[0]") \
        .Define("lepton_phi", "_isElectron ? Electron_phi[0] : Muon_phi[0]") \
        .Define("lepton_mass", "_isElectron ? Electron_mass[0] : Muon_mass[0]")

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
        "particleNet_XbbVsQCD", "particleNet_XqqVsQCD"
    ]
    
    for arr in fatjet_arrays:
        data[f"FatJet_{arr}"] = pad_or_truncate(data[f"FatJet_{arr}"], n_max_fatjets)
    
    data["FatJet_MASK"] = pad_or_truncate_bool(data["FatJet_MASK"], n_max_fatjets)
    
    return data

def prepare_dataset(input_dataframe, isSignal=False):
    arrays = get_array_names(isSignal)
    
    data = ak.from_rdataframe(input_dataframe, arrays)

    if isSignal:
        data = filter_signal_events(data)

    data = calculate_met_features(data)
    data = calculate_jet_features(data)
    data = calculate_fatjet_features(data)
    data = calculate_lepton_features(data)
    
    data = pad_jet_arrays(data, N_MAX_JETS)
    data = pad_fatjet_arrays(data, N_MAX_FATJETS)

    return data

def make_dataset(dataframes, output_path):
    sig_df, bkg_df = dataframes

    data_sig = prepare_dataset(sig_df, isSignal=True)
    data_bkg = prepare_dataset(bkg_df)

    data_sig["isSignal"] = True
    data_bkg["isSignal"] = False

    # data = ak.concatenate([data_sig, data_bkg], axis=0)
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
        "INPUTS/AK8Jets/particleNet_XbbVsQCD": data["FatJet_particleNet_XbbVsQCD"],
        "INPUTS/AK8Jets/particleNet_XqqVsQCD": data["FatJet_particleNet_XqqVsQCD"],
        
        # MET
        "INPUTS/MET/pt": data["MET_pt"],
        "INPUTS/MET/sinphi": data["MET_sinphi"],
        "INPUTS/MET/cosphi": data["MET_cosphi"],

        # Leptons
        "INPUTS/Lepton/pt": data["lepton_pt"],
        "INPUTS/Lepton/eta": data["lepton_eta"],
        "INPUTS/Lepton/sinphi": data["lepton_sinphi"],
        "INPUTS/Lepton/cosphi": data["lepton_cosphi"],
        "INPUTS/Lepton/mass": data["lepton_mass"],
        
        # Targets - Higgs
        "TARGETS/h/b1": data["Jet_hbb1_idx"],
        "TARGETS/h/b2": data["Jet_hbb2_idx"],
        "TARGETS/bh/bb": data["FatJet_hbb_idx"],
        
        # Targets - Vector boson
        "TARGETS/v/vq1": data["Jet_vqq1_idx"],
        "TARGETS/v/vq2": data["Jet_vqq2_idx"],
        "TARGETS/bv/qq": data["FatJet_vqq_idx"],
        
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
        .Redefine("FatJet_particleNet_XbbVsQCD", "FatJet_particleNet_XbbVsQCD[ak8_jets]") \
        .Redefine("FatJet_particleNet_XqqVsQCD", "FatJet_particleNet_XqqVsQCD[ak8_jets]")
    
    df = leptonSelections(df)
    df = df.Filter("_isLepton && Sum(jets) >= 2 && Sum(ak8_jets) >= 1")

    # define bb, qq definition
    if (isSignal):
        df = df.Define("GenPart_motherPdgId", "Take(GenPart_pdgId, GenPart_genPartIdxMother)") \
            .Filter("Sum((GenPart_pdgId == 11 || GenPart_pdgId == 13) && (abs(GenPart_motherPdgId) == 24)) == 1")
        
        df = df.Define("Higgs_eta", "GenPart_eta[4]") \
            .Define("Higgs_phi", "GenPart_phi[4]") \
            .Define("VBSJet1_eta", "GenPart_eta[5]") \
            .Define("VBSJet1_phi", "GenPart_phi[5]") \
            .Define("VBSJet2_eta", "GenPart_eta[6]") \
            .Define("VBSJet2_phi", "GenPart_phi[6]") \
            .Define("V_idx", "get_hadronic_gauge_boson_idx(GenPart_pdgId, GenPart_genPartIdxMother)") \
            .Filter("V_idx != -1") \
            .Define("V_eta", "GenPart_eta[V_idx]") \
            .Define("V_phi", "GenPart_phi[V_idx]") \
            .Define("qs_from_v", "abs(GenPart_pdgId) <= 5 && (abs(GenPart_motherPdgId == 24) || abs(GenPart_motherPdgId == 23))") \
            .Define("V_q1_eta", "GenPart_eta[qs_from_v][0]") \
            .Define("V_q1_phi", "GenPart_phi[qs_from_v][0]") \
            .Define("V_q2_eta", "GenPart_eta[qs_from_v][1]") \
            .Define("V_q2_phi", "GenPart_phi[qs_from_v][1]") \
            .Define("bs_from_higgs", "abs(GenPart_pdgId) == 5 && abs(GenPart_motherPdgId) == 25") \
            .Define("b1_eta", "GenPart_eta[bs_from_higgs][0]") \
            .Define("b1_phi", "GenPart_phi[bs_from_higgs][0]") \
            .Define("b2_eta", "GenPart_eta[bs_from_higgs][1]") \
            .Define("b2_phi", "GenPart_phi[bs_from_higgs][1]")

        # VBS matching        
        df = df.Define("vbs1_dR", "get_dR(VBSJet1_eta, VBSJet1_phi, Jet_eta, Jet_phi)") \
            .Define("vbs2_dR", "get_dR(VBSJet2_eta, VBSJet2_phi, Jet_eta, Jet_phi)") \
            .Define("vbs1_closest_jet", "ROOT::VecOps::Min(vbs1_dR)") \
            .Define("vbs2_closest_jet", "ROOT::VecOps::Min(vbs2_dR)") \
            .Define("Jet_vbs1_idx", "(vbs1_closest_jet < 0.4) ? (int)ROOT::VecOps::ArgMin(vbs1_dR) : -1") \
            .Define("Jet_vbs2_idx", "(vbs2_closest_jet < 0.4) ? (int)ROOT::VecOps::ArgMin(vbs2_dR) : -1") \
            .Filter("(Jet_vbs1_idx != Jet_vbs2_idx)")

        # Boosted Hbb matching
        df = df.Define("hbb_fatjet_dR", "get_dR(Higgs_eta, Higgs_phi, FatJet_eta, FatJet_phi)") \
            .Define("hbb_dR", "ROOT::VecOps::DeltaR(b1_eta, b1_phi, b2_eta, b2_phi)") \
            .Define("hbb_fatjet_idx_temp", "(int)ROOT::VecOps::ArgMin(hbb_fatjet_dR)") \
            .Define("hbb_fatjet_candidate_b1_dR", "ROOT::VecOps::DeltaR(FatJet_eta[hbb_fatjet_idx_temp], FatJet_phi[hbb_fatjet_idx_temp], b1_eta, b1_phi)") \
            .Define("hbb_fatjet_candidate_b2_dR", "ROOT::VecOps::DeltaR(FatJet_eta[hbb_fatjet_idx_temp], FatJet_phi[hbb_fatjet_idx_temp], b2_eta, b2_phi)") \
            .Define("hbb_isBoosted", "Min(hbb_fatjet_dR) < 0.8 && hbb_dR < 0.8 && hbb_fatjet_candidate_b1_dR < 0.8 && hbb_fatjet_candidate_b2_dR < 0.8") \
            .Define("FatJet_hbb_idx", "hbb_isBoosted ? (int)ROOT::VecOps::ArgMin(hbb_fatjet_dR) : -1")
        
        # Resolved Hbb matching
        df = df.Define("b1_jet_dR", "get_dR(b1_eta, b1_phi, Jet_eta, Jet_phi)") \
            .Define("b2_jet_dR", "get_dR(b2_eta, b2_phi, Jet_eta, Jet_phi)") \
            .Define("b1_closest_jet", "ROOT::VecOps::Min(b1_jet_dR)") \
            .Define("b2_closest_jet", "ROOT::VecOps::Min(b2_jet_dR)") \
            .Define("b1_closest_jet_idx", "(int)ROOT::VecOps::ArgMin(b1_jet_dR)") \
            .Define("b2_closest_jet_idx", "(int)ROOT::VecOps::ArgMin(b2_jet_dR)") \
            .Define("hbb_isResolved", "b1_closest_jet < 0.4 && b2_closest_jet < 0.4 && (b1_closest_jet_idx != b2_closest_jet_idx)") \
            .Define("Jet_hbb1_idx", "(hbb_isResolved && b1_closest_jet_idx != Jet_vbs1_idx && b1_closest_jet_idx != Jet_vbs2_idx) ? b1_closest_jet_idx : -1") \
            .Define("Jet_hbb2_idx", "(hbb_isResolved && b2_closest_jet_idx != Jet_vbs1_idx && b2_closest_jet_idx != Jet_vbs2_idx) ? b2_closest_jet_idx : -1")

        # Boosted Vqq matching
        df = df.Define("vqq_fatjet_dR", "get_dR(V_eta, V_phi, FatJet_eta, FatJet_phi)") \
            .Define("vqq_dR", "ROOT::VecOps::DeltaR(V_q1_eta, V_q1_phi, V_q2_eta, V_q2_phi)") \
            .Define("vqq_fatjet_idx_temp", "(int)ROOT::VecOps::ArgMin(vqq_fatjet_dR)") \
            .Define("vqq_fatjet_candidate_q1_dR", "ROOT::VecOps::DeltaR(FatJet_eta[vqq_fatjet_idx_temp], FatJet_phi[vqq_fatjet_idx_temp], V_q1_eta, V_q1_phi)") \
            .Define("vqq_fatjet_candidate_q2_dR", "ROOT::VecOps::DeltaR(FatJet_eta[vqq_fatjet_idx_temp], FatJet_phi[vqq_fatjet_idx_temp], V_q2_eta, V_q2_phi)") \
            .Define("vqq_isBoosted", "Min(vqq_fatjet_dR) < 0.8 && vqq_dR < 0.8 && vqq_fatjet_candidate_q1_dR < 0.8 && vqq_fatjet_candidate_q2_dR < 0.8") \
            .Define("FatJet_vqq_idx", "(vqq_isBoosted && FatJet_hbb_idx != vqq_fatjet_idx_temp) ? (int)ROOT::VecOps::ArgMin(vqq_fatjet_dR) : -1") \
        
        # Resolved Vqq matching
        df = df.Define("q1_jet_dR", "get_dR(V_q1_eta, V_q1_phi, Jet_eta, Jet_phi)") \
            .Define("q2_jet_dR", "get_dR(V_q2_eta, V_q2_phi, Jet_eta, Jet_phi)") \
            .Define("q1_closest_jet", "ROOT::VecOps::Min(q1_jet_dR)") \
            .Define("q2_closest_jet", "ROOT::VecOps::Min(q2_jet_dR)") \
            .Define("q1_closest_jet_idx", "(int)ROOT::VecOps::ArgMin(q1_jet_dR)") \
            .Define("q2_closest_jet_idx", "(int)ROOT::VecOps::ArgMin(q2_jet_dR)") \
            .Define("vqq_isResolved", "q1_closest_jet < 0.4 && q2_closest_jet < 0.4 && (q1_closest_jet_idx != q2_closest_jet_idx)") \
            .Define("Jet_vqq1_idx", "(vqq_isResolved && q1_closest_jet_idx != Jet_hbb1_idx && q1_closest_jet_idx != Jet_hbb2_idx && q1_closest_jet_idx != Jet_vbs1_idx && q1_closest_jet_idx != Jet_vbs2_idx) ? q1_closest_jet_idx : -1") \
            .Define("Jet_vqq2_idx", "(vqq_isResolved && q2_closest_jet_idx != Jet_hbb1_idx && q2_closest_jet_idx != Jet_hbb2_idx && q2_closest_jet_idx != Jet_vbs1_idx && q2_closest_jet_idx != Jet_vbs2_idx) ? q2_closest_jet_idx : -1")

    else:
        df = df.Define("FatJet_hbb_idx", "-1") \
            .Define("Jet_hbb1_idx", "-1") \
            .Define("Jet_hbb2_idx", "-1") \
            .Define("FatJet_vqq_idx", "-1") \
            .Define("Jet_vqq1_idx", "-1") \
            .Define("Jet_vqq2_idx", "-1") \
            .Define("Jet_vbs1_idx", "-1") \
            .Define("Jet_vbs2_idx", "-1") \
            .Define("hbb_isBoosted", "false") \
            .Define("hbb_isResolved", "false") \
            .Define("vqq_isBoosted", "false") \
            .Define("vqq_isResolved", "false") \
            .Filter("event % 10 == 0")
    return df

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--train", action="store_true", help="Is this a training sample?")
    args = parser.parse_args()

    # Open the input file
    if args.train:
        suffix = "_train"
    else:
        suffix = "_test"

    df_sig = genMatching(f"OneLep2FJ-sig{suffix}.json", isSignal=True)
    df_bkg = genMatching(f"OneLep2FJ-bkg{suffix}.json", isSignal=False)

    print("Number of evts in sig: ", df_sig.Count().GetValue())
    print("Number of evts in bkg: ", df_bkg.Count().GetValue())
    
    make_dataset((df_sig, df_bkg), output_path=f"output{suffix}.h5")