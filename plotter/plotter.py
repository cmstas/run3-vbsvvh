#!/usr/bin/env python
from dataclasses import dataclass, field

import ROOT as r
r.EnableImplicitMT()

import matplotlib.pyplot as plt
import mplhep as hep
hep.style.use(hep.style.CMS)

from pathlib import Path

@dataclass
class Hist1D:
    var: str
    xlabel: str
    binning: tuple
    hist_data: list = field(default_factory=list)
    hist_bkg: list = field(default_factory=list)
    hist_sig: list = field(default_factory=list)

@dataclass
class Hist2D:
    xvar: str
    yvar: str
    xlabel: str
    ylabel: str
    xbinning: tuple
    ybinning: tuple
    hist_data: list = field(default_factory=list)
    hist_bkg: list = field(default_factory=list)
    hist_sig: list = field(default_factory=list)

class Plotter:
    def __init__(self, sig=None, bkg=None, data=None, bkg_samples_labels=None, cut=None, year=None):
        self.year = year
        if data:
            self.df_data = r.RDataFrame("Events", data).Filter(cut) if cut else r.RDataFrame("Events", data)
        else:
            self.df_data = None
    
        if sig:
            self.df_sig = r.RDataFrame("Events", sig).Filter(cut) if cut else r.RDataFrame("Events", sig)
        else:
            self.df_sig = None
        
        if bkg:
            self.df_bkg = r.RDataFrame("Events", bkg).Filter(cut) if cut else r.RDataFrame("Events", bkg)
        else:
            self.df_bkg = None

        self.bkg_samples_labels = bkg_samples_labels
        if self.bkg_samples_labels is None:
            print("No background labels provided, will use single histogram for background")

    def plot1D(self, histogram, save=True):
        try:
            if self.df_data and self.df_bkg and histogram.hist_data and histogram.hist_bkg:
                hist_ratio = histogram.hist_data[0].GetValue().Clone()
                hist_bkg_total = histogram.hist_bkg[0].GetValue().Clone()
                [hist_bkg_total.Add(hist.GetValue()) for hist in histogram.hist_bkg[1:]]
                hist_ratio.Divide(hist_bkg_total)
            else:
                hist_ratio = None

            if hist_ratio:
                fig, ax = plt.subplots(2, 1, gridspec_kw={"height_ratios": (4, 1)})
                ax_main = ax[0]
                ax_ratio = ax[1]
            else:
                fig, ax = plt.subplots(1, 1)
                ax_main = ax
                ax_ratio = None

            if histogram.hist_sig:
                histogram.hist_sig[0].Scale(1000)
                hep.histplot(histogram.hist_sig, label="Signal x 1000", ax=ax_main, histtype="step", color="red", linewidth=2, yerr=False)

            if histogram.hist_bkg and self.bkg_samples_labels is None:
                hep.histplot(histogram.hist_bkg, ax=ax_main, histtype="fill", label="Background")
            elif histogram.hist_bkg and self.bkg_samples_labels:
                hep.histplot(histogram.hist_bkg, ax=ax_main, histtype="fill", stack=True, label=self.bkg_samples_labels.values())

            if histogram.hist_data:
                hep.histplot(histogram.hist_data, label="Data", ax=ax_main, histtype="errorbar", color="black")

            if hist_ratio:
                hep.histplot(hist_ratio, color="black", ax=ax_ratio, histtype="errorbar")
            
            year = self.year if self.year else "Run3"
            hep.cms.label("Preliminary", data=True, year=year, ax=ax_main)
            
            if hist_ratio:
                ax_main.legend()
                ax_main.set_xlabel("")
                ax_main.set_ylabel("Events")
                ax_ratio.set_xlabel(histogram.xlabel)
                ax_ratio.set_ylabel("Data / MC")
                ax_ratio.set_ylim(0.8, 1.2)
                ax_ratio.axhline(1, color="black", linestyle="--")
            else:
                ax_main.legend()
                ax_main.set_xlabel("")
                ax_main.set_ylabel("Events")
                ax_main.set_xlabel(histogram.xlabel)
            
            if save:
                Path("plots").mkdir(parents=True, exist_ok=True)
                fig.tight_layout()
                plt.savefig(f"plots/{histogram.var}.png")

            plt.close(fig)

        except Exception as e:
            print(f"Error {e} in", histogram.var)

    def plot2D(self, histogram):
        # TODO: Implement the 2D plotting functionality
        raise NotImplementedError

    def make_plots(self, hists, save=True):
        for histogram in hists:
            if isinstance(histogram, Hist1D):
                if self.df_data:
                    histogram.hist_data = [self.df_data.Histo1D((histogram.var, histogram.var, *histogram.binning), histogram.var, "weight")]
                if self.df_sig:
                    histogram.hist_sig = [self.df_sig.Histo1D((histogram.var, histogram.var, *histogram.binning), histogram.var, "weight")]
                if self.df_bkg:
                    histogram.hist_bkg = []
                    if self.bkg_samples_labels is None:
                        histogram.hist_bkg.append(
                            self.df_bkg.Histo1D((histogram.var, histogram.var, *histogram.binning), histogram.var, "weight")
                        )
                    else:
                        for sample in self.bkg_samples_labels.keys():
                            histogram.hist_bkg.append(
                                self.df_bkg.Filter(f"sample_type == \"{sample}\"").Histo1D((histogram.var, histogram.var, *histogram.binning), histogram.var, "weight")
                            )
            elif isinstance(histogram, Hist2D):
                pass

        for hist in hists:
            if isinstance(hist, Hist1D):
                self.plot1D(hist, save)
            elif isinstance(hist, Hist2D):
                self.plot2D(hist, save)