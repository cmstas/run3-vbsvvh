import ROOT as r
r.EnableImplicitMT(32)

from glob import glob

import mplhep as hep
import matplotlib.pyplot as plt
import numpy as np

r.gInterpreter.Declare("""
using namespace ROOT::VecOps;

RVec<float> VVdR(const RVec<float>& vec_eta1, const RVec<float>& vec_phi1, const RVec<float>& vec_eta2, const RVec<float>& vec_phi2) {
    if (vec_eta1.empty()) {
        return RVec<float>();
    }
    if (vec_eta2.empty()) {
        return RVec<float>(vec_eta1.size(), 999.0f);
    }
    RVec<float> out(vec_eta1.size());
    for (size_t i = 0; i < vec_eta1.size(); i++) {
        float mindR = 999.;
        for (size_t j = 0; j < vec_eta2.size(); j++) {
            float dR = ROOT::VecOps::DeltaR(vec_eta1[i], vec_eta2[j], vec_phi1[i], vec_phi2[j]);
            if (dR < mindR) {
                mindR = dR;
            }
        }
        out[i] = mindR;
    }
    return out;
}
""")

r3_data = glob("/ceph/cms/store/user/aaarora/VBSVVH_skim_v30/Run3_Data_v15_v30_1Lep1FJ/Muon*2024*/*.root")
# r2_data = glob("/ceph/cms/store/user/aaarora/VBSVVH_skim_v30/Run2_Data_v15_v30_1Lep1FJ/Muon*/*")

df = r.RDataFrame("Events", r3_data)

# MET filters (Run 3)
df = df.Filter(
    "Flag_goodVertices && "
    "Flag_globalSuperTightHalo2016Filter && "
    "Flag_EcalDeadCellTriggerPrimitiveFilter && "
    "Flag_BadPFMuonFilter && "
    "Flag_BadPFMuonDzFilter && "
    "Flag_hfNoisyHitsFilter && "
    "Flag_eeBadScFilter && "
    "Flag_ecalBadCalibFilter",
    "MET filters"
)

# Trigger selection (Run 3)
df = df.Filter("HLT_IsoMu24")

# Electron selections
df = (df.Define("Electron_SC_eta", "Electron_eta + Electron_deltaEtaSC")
    .Define("_vetoElectrons",
        "Electron_pt > 10 && "
        "abs(Electron_SC_eta) < 2.5 && "
        "((abs(Electron_SC_eta) <= 1.479 && abs(Electron_dxy) <= 0.05 && abs(Electron_dz) < 0.1) || "
        "((abs(Electron_SC_eta) > 1.479) && abs(Electron_dxy) <= 0.1 && abs(Electron_dz) < 0.2)) && "
        "Electron_cutBased >= 1")
    .Define("_looseElectrons", "_vetoElectrons && Electron_cutBased >= 2")
    .Define("nElectron_Veto", "nElectron == 0 ? 0 : Sum(_vetoElectrons)")
    .Define("nElectron_Loose", "nElectron_Veto == 0 ? 0 : Sum(_looseElectrons)")
    .Define("electron_pt", "Electron_pt[_vetoElectrons]")
    .Define("electron_eta", "Electron_eta[_vetoElectrons]")
    .Define("electron_phi", "Electron_phi[_vetoElectrons]")
    .Define("electron_mass", "Electron_mass[_vetoElectrons]")
    .Define("electron_charge", "Electron_charge[_vetoElectrons]")
    .Define("electron_cutBased", "Electron_cutBased[_vetoElectrons]")
    .Define("electron_isLoose", "electron_cutBased >= 2")
    .Define("electron_isMedium", "electron_cutBased >= 3")
    .Define("electron_isTight", "electron_cutBased >= 4")
    .Define("nElectron_Tight", "Sum(electron_isTight)"))

# Muon selections
df = (df.Define("_looseMuons",
        "Muon_pt > 10 && "
        "Muon_pfIsoId >= 2 && "
        "abs(Muon_eta) < 2.4 && "
        "abs(Muon_dxy) < 0.2 && "
        "abs(Muon_dz) < 0.5 && "
        "abs(Muon_sip3d) < 8 && "
        "Muon_looseId")
    .Define("_mediumMuons", "_looseMuons && Muon_pfIsoId >= 3 && Muon_mediumId")
    .Define("nMuon_Loose", "nMuon == 0 ? 0 : Sum(_looseMuons)")
    .Define("nMuon_Medium", "nMuon_Loose == 0 ? 0 : Sum(_mediumMuons)")
    .Define("muon_pt", "Muon_pt[_looseMuons]")
    .Define("muon_eta", "Muon_eta[_looseMuons]")
    .Define("muon_phi", "Muon_phi[_looseMuons]")
    .Define("muon_mass", "Muon_mass[_looseMuons]")
    .Define("muon_charge", "Muon_charge[_looseMuons]")
    .Define("muon_mediumId", "Muon_mediumId[_looseMuons]")
    .Define("muon_pfIsoId", "Muon_pfIsoId[_looseMuons]")
    .Define("muon_isMedium", "muon_mediumId && (muon_pfIsoId >= 3)")
    .Define("muon_isTight", "muon_mediumId && (muon_pfIsoId >= 4)")
    .Define("muon_isVeryTight", "muon_mediumId && (muon_pfIsoId >= 5)")
    .Define("nMuon_Tight", "Sum(muon_isTight)"))

# Combined lepton collection (pt-sorted)
df = (df.Define("lepton_pt", "Concatenate(electron_pt, muon_pt)")
    .Define("_leptonSorted", "Argsort(-lepton_pt)")
    .Redefine("lepton_pt", "Take(lepton_pt, _leptonSorted)")
    .Define("lepton_eta", "Take(Concatenate(electron_eta, muon_eta), _leptonSorted)")
    .Define("lepton_phi", "Take(Concatenate(electron_phi, muon_phi), _leptonSorted)")
    .Define("lepton_mass", "Take(Concatenate(electron_mass, muon_mass), _leptonSorted)")
    .Define("lepton_charge", "Take(Concatenate(electron_charge, muon_charge), _leptonSorted)"))

# AK8 fat jet selection (cleaned against leptons)
df = (df.Define("_dR_ak8_lep", "VVdR(FatJet_eta, FatJet_phi, lepton_eta, lepton_phi)")
    .Define("_good_ak8jets",
        "_dR_ak8_lep > 0.8 && "
        "FatJet_pt > 250 && "
        "abs(FatJet_eta) <= 2.5 && "
        "FatJet_msoftdrop > 40 && "
        "FatJet_jetId > 0")
    .Define("nFatJets", "Sum(_good_ak8jets)")
    .Define("FatJet_HvsQCD", "FatJet_globalParT3_Xbb / (FatJet_globalParT3_Xbb + FatJet_globalParT3_QCD)")
    .Define("FatJet_VvsQCD", "(FatJet_globalParT3_Xqq/3 + FatJet_globalParT3_Xcs) / (FatJet_globalParT3_Xqq/3 + FatJet_globalParT3_Xcs + FatJet_globalParT3_QCD)")
    .Define("nfatjet", "Sum(_good_ak8jets)")
    .Define("fatjet_pt", "FatJet_pt[_good_ak8jets]")
    .Define("fatjet_eta", "FatJet_eta[_good_ak8jets]")
    .Define("fatjet_phi", "FatJet_phi[_good_ak8jets]")
    .Define("fatjet_mass", "FatJet_mass[_good_ak8jets]")
    .Define("fatjet_msoftdrop", "FatJet_msoftdrop[_good_ak8jets]")
    .Define("fatjet_tau1", "FatJet_tau1[_good_ak8jets]")
    .Define("fatjet_tau2", "FatJet_tau2[_good_ak8jets]")
    .Define("fatjet_HvsQCD", "FatJet_HvsQCD[_good_ak8jets]")
    .Define("fatjet_VvsQCD", "FatJet_VvsQCD[_good_ak8jets]")
    .Define("ht_fatjets", "Sum(fatjet_pt)"))

# AK4 jet selection (Run 3; cleaned against leptons and fat jets)

df = (df.Define("_dR_ak4_lep", "VVdR(Jet_eta, Jet_phi, lepton_eta, lepton_phi)")
    .Define("_dR_ak4_fatjet", "VVdR(Jet_eta, Jet_phi, fatjet_eta, fatjet_phi)")
    .Define("_good_ak4jets",
        "_dR_ak4_lep > 0.4 && "
        "((Jet_pt > 30 && (abs(Jet_eta) <= 2.5 || abs(Jet_eta) >= 3.0) && abs(Jet_eta) < 2.5) || "
        "(Jet_pt > 50 && abs(Jet_eta) > 2.5 && abs(Jet_eta) < 3.0)) && "
        "Jet_jetId >= 2")
    .Define("njet", "Sum(_good_ak4jets)")
    .Define("jet_pt", "Jet_pt[_good_ak4jets]")
    .Define("jet_eta", "Jet_eta[_good_ak4jets]")
    .Define("jet_phi", "Jet_phi[_good_ak4jets]")
    .Define("jet_mass", "Jet_mass[_good_ak4jets]")
    .Define("ht_jet", "Sum(jet_pt)"))

# 1lep_2FJ event preselection (no trigger applied)
df = df.Filter(
    "((nMuon_Loose == 1 && nMuon_Tight == 1 && nElectron_Veto == 0 && nElectron_Loose == 0 && nElectron_Tight == 0) || "
    "(nMuon_Loose == 0 && nMuon_Tight == 0 && nElectron_Veto == 1 && nElectron_Loose == 1 && nElectron_Tight == 1)) && "
    "(lepton_pt[0] > 30)",
    "1-lepton selection"
)

df = df.Define("ht_new", "ht_jet + lepton_pt[0]")

df = df.Filter("nfatjet >= 2", "at-least 2 fat jets")
df = df.Filter("njet >= 2", "at-least 2 jets")

df.Count()

df_trigger = df.Filter("HLT_PFHT1050")

h_trigger = df_trigger.Histo1D(("HT_trigger", "HT distribution after trigger selection", 100, 600, 1200), "ht_jet")
h_before = df.Histo1D(("HT_before", "HT distribution before trigger selection", 100, 600, 1200), "ht_jet")

hb = h_before.GetValue()
ht = h_trigger.GetValue()

fig, (ax, ax_ratio) = plt.subplots(
    2, 1, sharex=True, figsize=(8, 6),
    gridspec_kw={"height_ratios": [3, 1]}
)

hep.histplot(hb, label="Before trigger selection", histtype="step", color="blue", ax=ax)
hep.histplot(ht, label="After trigger selection", histtype="step", color="red", ax=ax)

nbins = hb.GetNbinsX()
edges = np.array([hb.GetBinLowEdge(i) for i in range(1, nbins + 2)])
before_vals = np.array([hb.GetBinContent(i) for i in range(1, nbins + 1)])
trigger_vals = np.array([ht.GetBinContent(i) for i in range(1, nbins + 1)])

ratio = np.divide(trigger_vals, before_vals, out=np.zeros_like(trigger_vals, dtype=float), where=before_vals > 0)
ax_ratio.step(edges[:-1], ratio, where="post", color="black")
ax_ratio.set_ylabel("After/Before")
ax_ratio.set_xlabel("HT [GeV]")
ax.set_ylabel("Events")
ax.legend()

plt.tight_layout()


