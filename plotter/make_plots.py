#!/usr/bin/env python3

from plotter import *

NUM_BINS = 20

hists = [
        Hist1D("spanet_vbs_detajj", r"$VBS \Delta \eta_{jj}$", (NUM_BINS, 0, 10)),
        Hist1D("spanet_vbs_mjj", r"VBS $m_{jj}$ (GeV)", (NUM_BINS, 0, 2000)),
        Hist1D("spanet_bh_score", "bH tagger score", (NUM_BINS, 0, 1)),
        Hist1D("spanet_bv1_w_score", "bV W tagger score", (NUM_BINS, 0, 1)),
    ]

bkg_samples_labels = {
        "TTbar": r"$t\bar t$",
        "Boson": "Boson",
        "DY": "Drell-Yan",
        "WJets": "W + Jets",
        "Other": "Other",
}

BASE_PATH = "/data/userdata/aaarora/vbsvvhAnalysis/preselection/1Lep2FJ/"

plotter = Plotter(sig=BASE_PATH + "sig.root", bkg=[BASE_PATH + "bkg.root"], bkg_samples_labels=bkg_samples_labels, cut="spanet_bh_score > 0", year=2024)
plotter.make_plots(hists, save=True, density=True)
