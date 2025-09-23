#ifndef SELECTIONS_RUN2_H
#define SELECTIONS_RUN2_H

#pragma once

#include <unordered_map>

#include "ROOT/RDataFrame.hxx"

#include "utils.h"
#include "corrections.h"

using RNode = ROOT::RDF::RNode;

namespace Run2
{
    // FIXME
    const std::unordered_map<std::string, std::string> TriggerMap = {
        {"1Lep2FJ", "HLT_Ele30_WPTight_Gsf || HLT_IsoMu24"},
        {"1Lep1FJ", "HLT_Ele30_WPTight_Gsf || HLT_IsoMu24"},
        {"0Lep3FJ", "HLT_PFHT1050"},
        // {"0Lep2FJ", "HLT_BLAH"},
        // {"0Lep2FJMET", "HLT_BLAH"},
    };

    RNode runPreselection(RNode df_, std::string channel, bool noCut);

    RNode defineCorrectedCols(RNode df);
    RNode redefineGenColumns(RNode df);
    RNode flagSelections(RNode df);
    RNode TriggerSelections(RNode df_, std::string channel, const std::unordered_map<std::string, std::string> &trigger_map);
    RNode AK8JetsPreselection(RNode df_);
    RNode ElectronPreselection(RNode df_);
    RNode MuonPreselection(RNode df_);
    RNode AK4JetsPreselection(RNode df);
    RNode VBSJetsPreselection(RNode df);

    RNode ElectronSelections(RNode df);
    RNode MuonSelections(RNode df);

} // end namespace

#endif // SELECTIONS_H