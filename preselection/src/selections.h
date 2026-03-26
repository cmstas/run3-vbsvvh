#ifndef SELECTIONS_RUN3_H
#define SELECTIONS_RUN3_H

#pragma once

#include <unordered_map>
#include <regex>

#include "ROOT/RDataFrame.hxx"

#include "utils.h"
#include "corrections.h"

using RNode = ROOT::RDF::RNode;

RNode ElectronSelections(RNode df);
RNode MuonSelections(RNode df);
RNode AK8JetsSelection(RNode df);
RNode AK4JetsSelection(RNode df);

RNode runPreselection(RNode df_, std::string channel, bool isData, bool noCut);

// Trigger selections
RNode TriggerSelections(RNode df_, std::string trigger_logic_string);
// Trigger logic sringr, paste this in from the output of the script for building this logic
const std::string trigger_pass_no_overlap_string = "(is2016 && ( (((shortname==\"SingleElectron\") || !isData) && ((HLT_Ele27_WPTight_Gsf == true) || (HLT_Ele25_eta2p1_WPTight_Gsf == true) || (HLT_Ele27_eta2p1_WPLoose_Gsf == true))) && !(((HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ == true) || (HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL == true) || (HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL == true) || (HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ == true) || (HLT_TripleMu_12_10_5 == true) || (HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL == true) || (HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ == true) || (HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ == true) || (HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL == true) || (HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ == true) || (HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL == true) || (HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ == true) || (HLT_Mu8_DiEle12_CaloIdL_TrackIdL == true) || (HLT_DiMu9_Ele9_CaloIdL_TrackIdL == true) || (HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL == true) || (HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ == true) || (HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL == true) || (HLT_IsoMu24 == true) || (HLT_IsoTkMu24 == true) || (HLT_IsoMu22_eta2p1 == true) || (HLT_IsoTkMu22_eta2p1 == true) || (HLT_IsoMu22 == true) || (HLT_IsoTkMu22 == true) || (HLT_IsoMu27 == true)) && isData) )) || (is2017 && ( (((shortname==\"SingleElectron\") || !isData) && ((HLT_Ele32_WPTight_Gsf == true) || (HLT_Ele35_WPTight_Gsf == true))) && !(((HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ == true) || (HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 == true) || (HLT_TripleMu_12_10_5 == true) || (HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL == true) || (HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ == true) || (HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ == true) || (HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ == true) || (HLT_Mu8_DiEle12_CaloIdL_TrackIdL == true) || (HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ == true) || (HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ == true) || (HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL == true) || (HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ == true) || (HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL == true) || (HLT_IsoMu24 == true) || (HLT_IsoMu27 == true)) && isData) )) || (is2018 && ( (((shortname==\"EGamma\") || !isData) && ((HLT_Ele32_WPTight_Gsf == true) || (HLT_Ele35_WPTight_Gsf == true) || (HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL == true) || (HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ == true) || (HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL == true))) && !(((HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ == true) || (HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 == true) || (HLT_TripleMu_12_10_5 == true) || (HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL == true) || (HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ == true) || (HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ == true) || (HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ == true) || (HLT_Mu8_DiEle12_CaloIdL_TrackIdL == true) || (HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ == true) || (HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ == true) || (HLT_IsoMu24 == true) || (HLT_IsoMu27 == true)) && isData) ))";

#endif // SELECTIONS_H
