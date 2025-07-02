#!/usr/bin/env python3

from plotter import *

NUM_BINS = 20

hists = [
        Hist1D("vbs1_pt", r"VBS jet 1 $p_T$ (GeV)", (NUM_BINS, 40, 200)),
        Hist1D("vbs1_eta", r"VBS jet 1 $\eta$", (NUM_BINS, -2.5, 2.5)),
        Hist1D("vbs1_phi", r"VBS jet 1 $\phi$", (NUM_BINS, -3.14, 3.14)),
        Hist1D("vbs2_pt", r"VBS jet 2 $p_T$ (GeV)", (NUM_BINS, 40, 200)),
        Hist1D("vbs2_eta", r"VBS jet 2 $\eta$", (NUM_BINS, -2.5, 2.5)),
        Hist1D("vbs2_phi", r"VBS jet 2 $\phi$", (NUM_BINS, -3.14, 3.14)),
        Hist1D("vbs_score", r"VBS Tagger score", (NUM_BINS, 0, 1)),
        Hist1D("vbs_detajj", r"$VBS \Delta \eta_{jj}$", (NUM_BINS, 0, 10)),
        Hist1D("vbs_mjj", r"VBS $m_{jj}$ (GeV)", (NUM_BINS, 0, 2000)),
        Hist1D("hbb_pt", r"$p_T^{bb}$ (GeV)", (NUM_BINS, 0, 1000)),
        Hist1D("bh_score", r"Hbb ParticleNet score", (NUM_BINS, 0, 1)),
        Hist1D("hbb_eta", r"$\eta^{bb}$", (NUM_BINS, -2.5, 2.5)),
        Hist1D("hbb_phi", r"$\phi^{bb}$", (NUM_BINS, -3.14, 3.14)),
        Hist1D("wqq_pt", r"$p_T^{qq}$ (GeV)", (NUM_BINS, 0, 1000)),
        Hist1D("wqq_eta", r"$\eta^{qq}$", (NUM_BINS, -2.5, 2.5)),
        Hist1D("wqq_phi", r"$\phi^{qq}$", (NUM_BINS, -3.14, 3.14)),
        Hist1D("bv_score", r"Wqq ParticleNet score", (NUM_BINS, 0, 1)),
        Hist1D("hbb_score", r"$x_{bb}$ score", (NUM_BINS, 0, 1)),
        Hist1D("wqq_score", r"$x_{qq}$ score", (NUM_BINS, 0, 1)),
    ]

bkg_samples_labels = {
        "TTbar": r"$t\bar t$",
        "Boson": "Boson",
        "DY": "Drell-Yan",
        "WJets": "W + Jets",
        "QCD": "QCD",
}

BASE_PATH = "/data/userdata/aaarora/vbsvvhAnalysis/preselection/OneLep2FJ/"

# plotter = Plotter(sig=[BASE_PATH + "sig_1_5.root", BASE_PATH + "sig.root", BASE_PATH + "sig_2_0.root"], sig_samples_labels=["1.5", "1.7", "2.0"])
plotter = Plotter(sig=BASE_PATH + "sig.root", bkg=[BASE_PATH + "bkg.root"], bkg_samples_labels=bkg_samples_labels)
plotter.make_plots(hists, save=True, density=True)
