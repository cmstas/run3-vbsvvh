#ifndef WEIGHTS_H
#define WEIGHTS_H

#pragma once

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RVec.hxx"

#include "utils.h"
#include "corrections_run2.h"
#include "correction.h"

using correction::CorrectionSet;
using RNode = ROOT::RDF::RNode;
using ROOT::VecOps::RVec;

namespace Run2
{

// golden json
RNode goodRun(lumiMask golden, RNode df);
const auto LumiMask = lumiMask::fromJSON("corrections/goldenJson/Cert_271036-325175_13TeV_allRun2_JSON.json");

// l1 prefiring weights
RNode L1PreFiringWeight(RNode df);

RNode EWKCorrections(correction::CorrectionSet cset_ewk, RNode df);
const auto cset_ewk = *CorrectionSet::from_file("corrections/scalefactors/ewk/EWK.json");
// pileup sfs
RNode pileupScaleFactors(correction::CorrectionSet cset_pileup_2016preVFP, correction::CorrectionSet cset_pileup_2016postVFP, correction::CorrectionSet cset_pileup_2017, correction::CorrectionSet cset_pileup_2018, RNode df);
const auto cset_pileup_2016preVFP = *CorrectionSet::from_file("corrections/pileup/2016preVFP_puWeights.json");
const auto cset_pileup_2016postVFP = *CorrectionSet::from_file("corrections/pileup/2016postVFP_puWeights.json");
const auto cset_pileup_2017 = *CorrectionSet::from_file("corrections/pileup/2017_puWeights.json");
const auto cset_pileup_2018 = *CorrectionSet::from_file("corrections/pileup/2018_puWeights.json");

RNode pileupIDScaleFactors(correction::CorrectionSet cset_pileup_2016preVFP, correction::CorrectionSet cset_pileup_2016postVFP, correction::CorrectionSet cset_pileup_2017, correction::CorrectionSet cset_pileup_2018, RNode df);
const auto cset_pileupID_2016preVFP = *CorrectionSet::from_file("corrections/pileup/2016preVFP_puID.json");
const auto cset_pileupID_2016postVFP = *CorrectionSet::from_file("corrections/pileup/2016postVFP_puID.json");
const auto cset_pileupID_2017 = *CorrectionSet::from_file("corrections/pileup/2017_puID.json");
const auto cset_pileupID_2018 = *CorrectionSet::from_file("corrections/pileup/2018_puID.json");

// muon sfs
RNode muonScaleFactors_trigger(correction::CorrectionSet cset_muon_2016preVFP, correction::CorrectionSet cset_muon_2016postVFP, correction::CorrectionSet cset_muon_2017, correction::CorrectionSet cset_muon_2018, RNode df);
RNode muonScaleFactors_ID(correction::CorrectionSet cset_muon_2016preVFP, correction::CorrectionSet cset_muon_2016postVFP, correction::CorrectionSet cset_muon_2017, correction::CorrectionSet cset_muon_2018, RNode df);
RNode muonScaleFactors_ttHID(correction::CorrectionSet cset_muon_tth, RNode df);
RNode muonScaleFactors_ttHISO(correction::CorrectionSet cset_muon_tth, RNode df);
const auto cset_muon_ID_2016preVFP = *CorrectionSet::from_file("corrections/scalefactors/muon/2016preVFP_muon_ID.json");
const auto cset_muon_ID_2016postVFP = *CorrectionSet::from_file("corrections/scalefactors/muon/2016postVFP_muon_ID.json");
const auto cset_muon_ID_2017 = *CorrectionSet::from_file("corrections/scalefactors/muon/2017_muon_ID.json");
const auto cset_muon_ID_2018 = *CorrectionSet::from_file("corrections/scalefactors/muon/2018_muon_ID.json");
const auto cset_muon_ttH = *CorrectionSet::from_file("corrections/scalefactors/muon/ttH_Muon_SF.json");

// electron sfs
RNode electronScaleFactors_Reco(correction::CorrectionSet cset_electron_2016preVFP, correction::CorrectionSet cset_electron_2016postVFP, correction::CorrectionSet cset_electron_2017, correction::CorrectionSet cset_electron_2018, RNode df);
RNode electronScaleFactors_ID(correction::CorrectionSet cset_electron_2016preVFP, correction::CorrectionSet cset_electron_2016postVFP, correction::CorrectionSet cset_electron_2017, correction::CorrectionSet cset_electron_2018, RNode df);
RNode electronScaleFactors_ttHID(correction::CorrectionSet cset_electron_tth, RNode df);
RNode electronScaleFactors_ttHISO(correction::CorrectionSet cset_electron_tth, RNode df);
RNode electronScaleFactors_Trigger(correction::CorrectionSet cset_electron_trigger, RNode df);
const auto cset_electron_ID_2016preVFP = *CorrectionSet::from_file("corrections/scalefactors/electron/2016preVFP_electron_ID.json");
const auto cset_electron_ID_2016postVFP = *CorrectionSet::from_file("corrections/scalefactors/electron/2016postVFP_electron_ID.json");
const auto cset_electron_ID_2017 = *CorrectionSet::from_file("corrections/scalefactors/electron/2017_electron_ID.json");
const auto cset_electron_ID_2018 = *CorrectionSet::from_file("corrections/scalefactors/electron/2018_electron_ID.json");
const auto cset_electron_ttH = *CorrectionSet::from_file("corrections/scalefactors/electron/ttH_Electron_SF.json");
const auto cset_electron_trigger = *CorrectionSet::from_file("corrections/scalefactors/electron/trigger.json");

// particle net sfs
RNode PNET_W_ScaleFactors_2016preVFP(correction::CorrectionSet cset_pnet_w, RNode df);
RNode PNET_W_ScaleFactors_2016postVFP(correction::CorrectionSet cset_pnet_w, RNode df);
RNode PNET_W_ScaleFactors_2017(correction::CorrectionSet cset_pnet_w, RNode df);
RNode PNET_W_ScaleFactors_2018(correction::CorrectionSet cset_pnet_w, RNode df);
RNode PNET_H_ScaleFactors_2016preVFP(correction::CorrectionSet cset_pnet_h, RNode df);
RNode PNET_H_ScaleFactors_2016postVFP(correction::CorrectionSet cset_pnet_h, RNode df);
RNode PNET_H_ScaleFactors_2017(correction::CorrectionSet cset_pnet_h, RNode df);
RNode PNET_H_ScaleFactors_2018(correction::CorrectionSet cset_pnet_h, RNode df);
const auto cset_pnet_w = *CorrectionSet::from_file("corrections/scalefactors/particlenet/pnet.json");
const auto cset_pnet_h = *CorrectionSet::from_file("corrections/scalefactors/particlenet/pnet.json");

// btagging sfs
RNode bTaggingScaleFactors_LF(correction::CorrectionSet cset_btag_2016preVFP, correction::CorrectionSet cset_btag_2016postVFP, correction::CorrectionSet cset_btag_2017, correction::CorrectionSet cset_btag_2018, correction::CorrectionSet cset_btag_eff, RNode df);
RNode bTaggingScaleFactors_HF(correction::CorrectionSet cset_btag_2016preVFP, correction::CorrectionSet cset_btag_2016postVFP, correction::CorrectionSet cset_btag_2017, correction::CorrectionSet cset_btag_2018, correction::CorrectionSet cset_btag_eff, RNode df);
const auto cset_btag_2016preVFP = *CorrectionSet::from_file("corrections/scalefactors/btagging/2016preVFP.json");
const auto cset_btag_2016postVFP = *CorrectionSet::from_file("corrections/scalefactors/btagging/2016postVFP.json");
const auto cset_btag_2017 = *CorrectionSet::from_file("corrections/scalefactors/btagging/2017.json");
const auto cset_btag_2018 = *CorrectionSet::from_file("corrections/scalefactors/btagging/2018.json");
const auto cset_btag_eff = *CorrectionSet::from_file("corrections/scalefactors/btagging/btag_eff.json");

RNode PSWeight_FSR(RNode df);
RNode PSWeight_ISR(RNode df);
RNode LHEScaleWeight_muF(RNode df);
RNode LHEScaleWeight_muR(RNode df);
RNode LHEWeights_pdf(RNode df);

RNode applyDataWeights(RNode df);
RNode applyMCWeights(RNode df);

} // end namespace Run2

#endif