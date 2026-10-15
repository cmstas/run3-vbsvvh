#ifndef WQQ_SF_H
#define WQQ_SF_H

#pragma once

#include <string>
#include <vector>

#include "ROOT/RDataFrame.hxx"

using RNode = ROOT::RDF::RNode;

/*
############################################
WQQ TAGGING SCALE FACTORS (GloParT-3)
############################################

Per-fat-jet Wqq tagging SFs, for the working points each channel lists in
applywqqsf.yaml:

  0lep_1FJ_met:
    cut:   [0.80]
    flags: [[0.60], [0.40, 0.60, 0.80]]

  cut    one branch triple per WP, for a hard cut score > WP: a passing jet gets the
         cumulative SF, a failing jet 1 (its event is cut away). No efficiency enters.
  flags  one branch triple per WP set, for working points used as flags (e.g. NN
         inputs), where passing and failing jets both stay. The weight is the
         nested-category weight of cumulativeWPCategoryWeight (category_weight.h):
         a one-WP set gives the pass SF / fail SF (1 - SF*eff)/(1 - eff); a multi-WP
         set gives the weight of the exclusive score interval the jet falls in. Needs
         the candle MC efficiency (Wqq_eff_MC) in the era file.

Branches, on the full FatJet_ collection (the fatjet_ aliases follow):
  FatJet_WqqSF_cut0p80_{nom,up,down}
  FatJet_WqqSF_flags0p60_{nom,up,down}
  FatJet_WqqSF_flags0p40_0p60_0p80_{nom,up,down}

Only true W jets (FatJet_genMatch 1 prompt W or 2 W from top) get SFs; every other
jet, and every jet in data, gets 1. The SFs are evaluated at the nominal FatJet_pt
and are the same in every JES/JER variation. up/down shift the SFs of every WP of a
set together.

Payloads: corrections/scalefactors/WqqSF/Wqq_SF_<year>.json (make_wqq_sf_jsons.py),
correction Wqq_SF (WP, pt, eta, systematic) and, for flags, Wqq_eff_MC (WP, pt).
*/
RNode applyWqqScaleFactors(RNode df, const std::string &channel, bool isData);

#endif
