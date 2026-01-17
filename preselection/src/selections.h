#ifndef SELECTIONS_RUN3_H
#define SELECTIONS_RUN3_H

#pragma once

#include <unordered_map>

#include "ROOT/RDataFrame.hxx"

#include "utils.h"
#include "corrections.h"

using RNode = ROOT::RDF::RNode;

const std::unordered_map<std::string, std::string> TriggerMap = {
    {"1Lep2FJ", "HLT_Ele30_WPTight_Gsf || HLT_IsoMu24"},
    {"1Lep1FJ", "HLT_Ele30_WPTight_Gsf || HLT_IsoMu24"},
    {"0Lep3FJ", "HLT_PFHT1050"},
    // {"0Lep2FJ", "HLT_BLAH"},
    // {"0Lep2FJMET", "HLT_BLAH"},
};

RNode TriggerSelections(RNode df_, std::string channel, const std::unordered_map<std::string, std::string> &trigger_map);
RNode ElectronSelections(RNode df);
RNode MuonSelections(RNode df);
RNode AK8JetsSelection(RNode df, std::string run_number);
RNode AK4JetsSelection(RNode df, std::string run_number);

RNode runPreselection(RNode df_, std::string channel, std::string run_number, bool noCut);

#endif // SELECTIONS_H