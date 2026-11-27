#ifndef XBB_SF_H
#define XBB_SF_H

#pragma once

#include <string>
#include <vector>

#include "ROOT/RDataFrame.hxx"

using RNode = ROOT::RDF::RNode;

/*
############################################
XBB TAGGING SCALE FACTORS (GloParT-3)
############################################

Per-fat-jet Xbb tagging SFs, for the working points each channel lists in
applyxbbsf.yaml (same format and behaviour as the Wqq SFs, see wqq_sf.h):

  2lep_1FJ:
    flags: [[0.50]]

  flags  one branch triple per WP set, for working points used as a hard cut or as
         flags (e.g. NN inputs). The weight is the nested-category weight of
         cumulativeWPCategoryWeight (category_weight.h): a one-WP set gives the pass
         SF / fail SF (1 - SF*eff)/(1 - eff); a multi-WP set gives the weight of the
         exclusive score interval the jet falls in. With a hard cut at the WP the
         failing jets are gone, so only the pass SF enters. Needs the candle MC
         efficiency (Xbb_eff_MC) in the era file.

Branches, on the full FatJet_ collection (the fatjet_ aliases follow):
  FatJet_XbbSF_flags0p50_{nom,up,down}

The score is FatJet_Hbb = Xbb / (Xbb + QCD). Only true b jets (FatJet_genMatch 5,
the jet contains both b quarks of an H->bb) get SFs; every other jet, and every jet
in data, gets 1. The SFs are evaluated at the nominal FatJet_pt and are the same in
every JES/JER variation. up/down shift the SFs of every WP of a set together.

Payloads: corrections/scalefactors/XbbSF/Xbb_SF_<year>.json (make_xbb_sf_jsons.py),
corrections Xbb_SF (WP, pt, systematic) and Xbb_eff_MC (WP, pt). Unlike Wqq_SF, Xbb_SF
has no eta input.
*/
RNode applyXbbScaleFactors(RNode df, const std::string &channel, bool isData);

#endif
