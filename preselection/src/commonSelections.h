#ifndef COMMONSELECTIONS_H
#define COMMONSELECTIONS_H

#include "ROOT/RDataFrame.hxx"

#include "utils.h"
#include "corrections.h"

using RNode = ROOT::RDF::RNode;

RNode CommonSelections(RNode df);
RNode EventFilters(RNode df);
RNode ElectronSelections(RNode df);
RNode MuonSelections(RNode df);
RNode AK8JetsSelection(RNode df);
RNode AK4JetsSelection(RNode df);
RNode VBSJetsSelection(RNode df);


#endif // COMMONSELECTIONS_H