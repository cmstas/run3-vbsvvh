#!/usr/env/python3

import ROOT as r
r.EnableImplicitMT()

import numpy as np

import argparse

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description="Reweight Signal")
    argparser.add_argument("--input", type=str, required=True, help="Input file")
    argparser.add_argument("--output", type=str, required=True, help="Output file")
    argparser.add_argument("--coupling", type=float, required=True, help="Coupling value")
    args = argparser.parse_args()

    df = r.RDataFrame("Events", args.input)
    r.RDF.Experimental.AddProgressBar(df)

    couplings = [round(i, 1) for i in np.linspace(-2, 2, 41)]
    couplings.remove(1.7)

    coupling_idx = lambda x: couplings.index(x)

    reweight = f"weight * LHEReweightingWeight[{coupling_idx(args.coupling)}]"
    df = df.Redefine("weight", reweight)

    df.Snapshot("Events", args.output)