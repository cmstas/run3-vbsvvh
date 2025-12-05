#ifndef CORRECTIONS_H
#define CORRCTIONS_H

#pragma once

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RVec.hxx"

#include "correction.h"

#include "TRandom3.h"

using correction::CorrectionSet;
using RNode = ROOT::RDF::RNode;
using ROOT::VecOps::RVec;

namespace Run2
{

float looseDFBtagWP(std::string year);
float mediumDFBtagWP(std::string year);
float tightDFBtagWP(std::string year);

RNode defineCorrectedCols(RNode df);

RNode HEMCorrection(RNode df);

RNode METUnclusteredCorrections(RNode df, std::string variation);

// jet mass scale and resolution corrections
RNode JMR_Corrections(correction::CorrectionSet cset_jet_mass_scale, RNode df, std::string variation); 
const auto cset_jmr = *CorrectionSet::from_file("corrections/scalefactors/particlenet/jmar.json");

RNode JMS_Corrections(correction::CorrectionSet cset_jet_mass_scale, RNode df, std::string variation);
const auto cset_jms = *CorrectionSet::from_file("corrections/scalefactors/particlenet/jmar.json");

RNode JetEnergyCorrection(correction::CorrectionSet cset_jerc_2016preVFP, correction::CorrectionSet cset_jerc_2016postVFP, correction::CorrectionSet cset_jerc_2017, correction::CorrectionSet cset_jerc_2018, RNode df, std::string JEC_type, std::string variation);
RNode JetEnergyResolution(correction::CorrectionSet cset_jerc_2016preVFP, correction::CorrectionSet cset_jerc_2016postVFP, correction::CorrectionSet cset_jerc_2017, correction::CorrectionSet cset_jerc_2018, correction::CorrectionSet cset_jer_smear, RNode df, std::string variation);
const auto cset_jerc_2016preVFP = *CorrectionSet::from_file("corrections/jets/2016preVFP_jet_jerc.json");
const auto cset_jerc_2016postVFP = *CorrectionSet::from_file("corrections/jets/2016postVFP_jet_jerc.json");
const auto cset_jerc_2017 = *CorrectionSet::from_file("corrections/jets/2017_jet_jerc.json");
const auto cset_jerc_2018 = *CorrectionSet::from_file("corrections/jets/2018_jet_jerc.json");
const auto cset_jer_smear = *CorrectionSet::from_file("corrections/jets/jer_smear.json");

} // end namespace Run2
#endif