#ifndef SELECTIONS_RUN3_H
#define SELECTIONS_RUN3_H

#pragma once

#include <unordered_map>
#include <regex>

#include "ROOT/RDataFrame.hxx"

#include "utils.h"
#include "corrections.h"

using RNode = ROOT::RDF::RNode;

const std::unordered_map<std::string, std::string> TriggerMap = {
	{"1lep_2FJ", 
        "(is2016 && (HLT_Ele27_eta2p1_WPTight_Gsf == true || HLT_IsoMu24 == true || HLT_IsoTkMu24 == true)) || "
        "(is2017 && (HLT_Ele32_WPTight_Gsf_L1DoubleEG == true || HLT_IsoMu27 == true)) || "
        "(is2018 && (HLT_Ele32_WPTight_Gsf == true || HLT_IsoMu24 == true)) || "
        "(isRun3 && (HLT_Ele30_WPTight_Gsf == true || HLT_IsoMu24 == true))"
    },
    {"1lep_1FJ", 
        "(is2016 && (HLT_Ele27_eta2p1_WPTight_Gsf == true || HLT_IsoMu24 == true || HLT_IsoTkMu24 == true)) || "
        "(is2017 && (HLT_Ele32_WPTight_Gsf_L1DoubleEG == true || HLT_IsoMu27 == true)) || "
        "(is2018 && (HLT_Ele32_WPTight_Gsf == true || HLT_IsoMu24 == true)) || "
        "(isRun3 && (HLT_Ele30_WPTight_Gsf == true || HLT_IsoMu24 == true))"
    },
    {"0lep_3FJ", "HLT_PFHT1050"},
};

RNode ElectronSelections(RNode df);
RNode MuonSelections(RNode df);
RNode AK8JetsSelection(RNode df);
RNode AK4JetsSelection(RNode df);

RNode runPreselection(RNode df_);

#endif // SELECTIONS_H
