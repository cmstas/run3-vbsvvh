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

if __name__ == "__main__":
    df_bkg = None
    df_sig = genMatching(f"/home/users/aaarora/phys/run3/run3-vbsvvh/preselection/etc/1Lep2FJ-sig.json", isSignal=True)
    # df_bkg = genMatching(f"/home/users/aaarora/phys/run3/SPANet/input_proc/1Lep2FJ1.5-bkg.json", isSignal=False)

    df_sig.Snapshot("Events", "out.root")
    print("Number of evts in sig: ", df_sig.Count().GetValue())
    # print("Number of evts in bkg: ", df_bkg.Count().GetValue())
    
    make_dataset((df_sig, df_bkg), output_path=f"../data/output.h5")