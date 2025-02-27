#ifndef SELECTIONS_H
#define SELECTIONS_H

#include "ROOT/RDataFrame.hxx"

#include "utils.h"

using RNode = ROOT::RDF::RNode;

namespace OneLep2FJ {

  RNode TriggerSelections(RNode df);
  RNode ElectronSelections(RNode df);
  RNode MuonSelections(RNode df);
  RNode LeptonSelections(RNode df);
  RNode HbbSelections(RNode df);
  RNode VqqSelections(RNode df);
  RNode AK4Selections(RNode df);
  RNode VBSJetsSelections(RNode df);

  RNode ObjectSelections(RNode df);
  RNode EventSelections(RNode df);
  RNode runAnalysis(RNode df_);

} // OneLep2FJ

#endif // SELECTIONS_H