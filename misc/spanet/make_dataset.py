from argparse import ArgumentParser

import ROOT as r
r.EnableImplicitMT(64)

import h5py

import awkward as ak
import numpy as np

N_MAX_JETS = 10
N_MAX_FATJETS = 3

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
        "Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass", "Jet_isLooseBTag", "Jet_isMediumBTag", "Jet_isTightBTag",
        # AK8 jets
        "FatJet_pt", "FatJet_eta", "FatJet_phi", "FatJet_mass", "FatJet_nConstituents", 
        "FatJet_globalParT3_Xbb", "FatJet_globalParT3_Xcc", "FatJet_globalParT3_Xcs", "FatJet_globalParT3_Xqq", "FatJet_globalParT3_QCD",
        # Global
        "PuppiMET_pt", "PuppiMET_phi",
        # Lepton
        "Lepton_pt", "Lepton_eta", "Lepton_phi", "Lepton_mass", "Lepton_charge",
    ]

    if isSignal:
        arrays.extend([
            "truth_h_idx", "truth_b1_idx", "truth_b2_idx",
            "truth_v1_idx", "truth_v1q1_idx", "truth_v1q2_idx", 
            "truth_v2_idx", "truth_v2q1_idx", "truth_v2q2_idx",
            "truth_vbs1_idx", "truth_vbs2_idx"
        ])
    
    return arrays

def calculate_met_features(data):
    """Add MET phi-related features."""
    data["MET_cosphi"] = np.cos(data["PuppiMET_phi"])
    data["MET_sinphi"] = np.sin(data["PuppiMET_phi"])
    data["MET_MASK"] = ak.ones_like(data["PuppiMET_pt"], dtype=bool)
    return data

def calculate_jet_features(data):
    """Add AK4 jet features."""
    data["Jet_cosphi"] = np.cos(data["Jet_phi"])
    data["Jet_sinphi"] = np.sin(data["Jet_phi"])
    data["Jet_isLooseBTag"] = data["Jet_isLooseBTag"]
    data["Jet_isMediumBTag"] = data["Jet_isMediumBTag"]
    data["Jet_isTightBTag"] = data["Jet_isTightBTag"]
    data["Jet_MASK"] = ak.ones_like(data["Jet_pt"], dtype=bool)
    return data

def calculate_fatjet_features(data):
    """Add AK8 fatjet features."""
    data["FatJet_cosphi"] = np.cos(data["FatJet_phi"])
    data["FatJet_sinphi"] = np.sin(data["FatJet_phi"])
    data["FatJet_MASK"] = ak.ones_like(data["FatJet_pt"], dtype=bool)
    return data

def calculate_lepton_features(data):
    """Add lepton features."""
    data["Lepton1_pt"] = data["Lepton_pt"][:, 0]
    data["Lepton1_eta"] = data["Lepton_eta"][:, 0]
    data["Lepton1_phi"] = data["Lepton_phi"][:, 0]
    data["Lepton1_mass"] = data["Lepton_mass"][:, 0]
    data["Lepton1_charge"] = data["Lepton_charge"][:, 0]
    data["Lepton1_cosphi"] = np.cos(data["Lepton1_phi"])
    data["Lepton1_sinphi"] = np.sin(data["Lepton1_phi"])
    data["Lepton1_MASK"] = (data["nLeptons"] >= 1)

    data["Lepton2_pt"] = data["Lepton_pt"][:, 1]
    data["Lepton2_eta"] = data["Lepton_eta"][:, 1]
    data["Lepton2_phi"] = data["Lepton_phi"][:, 1]
    data["Lepton2_mass"] = data["Lepton_mass"][:, 1]
    data["Lepton2_charge"] = data["Lepton_charge"][:, 1]
    data["Lepton2_cosphi"] = np.cos(data["Lepton2_phi"])
    data["Lepton2_sinphi"] = np.sin(data["Lepton2_phi"])
    data["Lepton2_MASK"] = (data["nLeptons"] >= 2)
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

def pad_lepton_arrays(data):
    """Pad or truncate lepton arrays if needed."""
    lepton_arrays = [
        "pt", "eta", "mass", "phi", "charge"
    ]

    data["nLeptons"] = ak.num(data["Lepton_pt"])

    for arr in lepton_arrays:
        data[f"Lepton_{arr}"] = pad_or_truncate(data[f"Lepton_{arr}"], 2)

    return data

def prepare_dataset(input_dataframe, isSignal=False):
    arrays = get_array_names(isSignal)
    
    data = ak.from_rdataframe(input_dataframe, arrays)

    data = calculate_met_features(data)
    data = calculate_jet_features(data)
    data = calculate_fatjet_features(data)

    data = pad_lepton_arrays(data)
    data = calculate_lepton_features(data)

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
        "INPUTS/MET/MASK": data["MET_MASK"],
        "INPUTS/MET/pt": data["PuppiMET_pt"],
        "INPUTS/MET/sinphi": data["MET_sinphi"],
        "INPUTS/MET/cosphi": data["MET_cosphi"],

        # Lepton

        "INPUTS/Lepton1/MASK": data["Lepton1_MASK"],
        "INPUTS/Lepton1/pt": data["Lepton1_pt"],
        "INPUTS/Lepton1/eta": data["Lepton1_eta"],
        "INPUTS/Lepton1/sinphi": data["Lepton1_sinphi"],
        "INPUTS/Lepton1/cosphi": data["Lepton1_cosphi"],
        "INPUTS/Lepton1/mass": data["Lepton1_mass"],
        "INPUTS/Lepton1/charge": data["Lepton1_charge"],

        "INPUTS/Lepton2/MASK": data["Lepton2_MASK"],
        "INPUTS/Lepton2/pt": data["Lepton2_pt"],
        "INPUTS/Lepton2/eta": data["Lepton2_eta"],
        "INPUTS/Lepton2/sinphi": data["Lepton2_sinphi"],
        "INPUTS/Lepton2/cosphi": data["Lepton2_cosphi"],
        "INPUTS/Lepton2/mass": data["Lepton2_mass"],
        "INPUTS/Lepton2/charge": data["Lepton2_charge"],

        # Targets - Higgs
        "TARGETS/h/b1": data["truth_b1_idx"],
        "TARGETS/h/b2": data["truth_b2_idx"],
        "TARGETS/bh/bb": data["truth_h_idx"],
        
        # Targets - Vector boson
        "TARGETS/v1/v1q1": data["truth_v1q1_idx"],
        "TARGETS/v1/v1q2": data["truth_v1q2_idx"],
        "TARGETS/v2/v2q1": data["truth_v2q1_idx"],
        "TARGETS/v2/v2q2": data["truth_v2q2_idx"],
        "TARGETS/bv1/qq1": data["truth_v1_idx"],
        "TARGETS/bv2/qq2": data["truth_v2_idx"],
        
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
    df_sig = r.RDataFrame("Events", "../data/sig_spanet_training_data.root")
    
    make_dataset((df_sig, df_bkg), output_path=f"../data/output.h5")