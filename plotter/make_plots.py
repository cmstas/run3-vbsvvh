#!/usr/bin/env python3

from plotter import *

hists = [
        Hist1D("muon_pt", r"Muon $p_T$ (GeV)", (40, 40, 200)),
        Hist1D("muon_eta", r"Muon $\eta$", (40, -2.5, 2.5)),
        Hist1D("muon_phi", r"Muon $\phi$", (40, -3.14, 3.14)),
        Hist1D("electron_pt", r"Electron $p_T$ (GeV)", (40, 40, 200)),
        Hist1D("electron_sceta", r"Electron $\eta$", (40, -2.5, 2.5)),
        Hist1D("electron_phi", r"Electron $\phi$", (40, -3.14, 3.14)),
        Hist1D("lepton_pt", r"Lepton $p_T$ (GeV)", (40, 40, 200)),
        Hist1D("lepton_eta", r"Lepton $\eta$", (40, -2.5, 2.5)),
        Hist1D("lepton_phi", r"Lepton $\phi$", (40, -3.14, 3.14)),
        Hist1D("vbs1_pt", r"VBS jet 1 $p_T$ (GeV)", (40, 40, 200)),
        Hist1D("vbs1_eta", r"VBS jet 1 $\eta$", (40, -2.5, 2.5)),
        Hist1D("vbs1_phi", r"VBS jet 1 $\phi$", (40, -3.14, 3.14)),
        Hist1D("vbs2_pt", r"VBS jet 2 $p_T$ (GeV)", (40, 40, 200)),
        Hist1D("vbs2_eta", r"VBS jet 2 $\eta$", (40, -2.5, 2.5)),
        Hist1D("vbs2_phi", r"VBS jet 2 $\phi$", (40, -3.14, 3.14)),
        Hist1D("vbs_detajj", r"$\Delta \eta_{jj}$", (40, 0, 10)),
        Hist1D("vbs_dphijj", r"$\Delta \phi_{jj}$", (40, 0, 3.14)),
        Hist1D("vbs_mjj", r"$m_{jj}$ (GeV)", (40, 0, 2000)),
        Hist1D("vbs_ptjj", r"$p_T^{jj}$ (GeV)", (40, 0, 1000)),
        Hist1D("vbs_score", r"VBS Tagger score", (40, 0, 1)),
        Hist1D("hbb_pt", r"$p_T^{bb}$ (GeV)", (40, 0, 1000)),
        Hist1D("hbb_score", r"Hbb ParticleNet score", (40, 0, 1)),
        Hist1D("hbb_eta", r"$\eta^{bb}$", (40, -2.5, 2.5)),
        Hist1D("hbb_phi", r"$\phi^{bb}$", (40, -3.14, 3.14)),
        Hist1D("wqq_pt", r"$p_T^{qq}$ (GeV)", (40, 0, 1000)),
        Hist1D("wqq_eta", r"$\eta^{qq}$", (40, -2.5, 2.5)),
        Hist1D("wqq_phi", r"$\phi^{qq}$", (40, -3.14, 3.14)),
        Hist1D("wqq_score", r"Wqq ParticleNet score", (40, 0, 1)),
    ]

bkg_samples_labels = {
        "DY": "Drell-Yan",
        "TTbar": r"$t\bar t$",
        "WJets": "W + Jets",
        "QCD": "QCD",
        "Boson": "Boson"
}

BASE_PATH = "/data/userdata/aaarora/vbsvvhAnalysis/preselection/OneLep2FJ/"

# plotter = Plotter(sig=[BASE_PATH + "sig_1_5.root", BASE_PATH + "sig.root", BASE_PATH + "sig_2_0.root"], sig_samples_labels=["1.5", "1.7", "2.0"])
plotter = Plotter(bkg=[BASE_PATH + "bkg.root"], bkg_samples_labels=bkg_samples_labels, data=[BASE_PATH + "data.root"])
plotter.make_plots(hists, save=True)