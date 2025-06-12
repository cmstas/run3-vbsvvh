#ifndef COMMONSELECTIONS_H
#define COMMONSELECTIONS_H

#include "ROOT/RDataFrame.hxx"

#include "inference.hpp"
#include "TMVA/RBDT.hxx"

#include "utils.h"
#include "corrections.h"

using RNode = ROOT::RDF::RNode;

namespace CommonSelections {
    RNode EventFilters(RNode df);
    RNode ElectronSelections(RNode df);
    RNode MuonSelections(RNode df);
    RNode AK8JetsSelection(RNode df);
    RNode AK4JetsSelection(RNode df);
    RNode VBSJetSelections(RNode df_, TMVA::Experimental::RBDT &vbstagger);
}
#endif // COMMONSELECTIONS_H