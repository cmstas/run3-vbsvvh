#ifndef SELECTIONS_ZEROLEP_H
#define SELECTIONS_ZEROLEP_H

#include "ROOT/RDataFrame.hxx"

#include "utils.h"

using RNode = ROOT::RDF::RNode;

namespace ZeroLep3FJ {

  RNode TriggerSelections(RNode df);
  RNode LeptonSelections(RNode df);
  RNode EventSelections(RNode df);
  RNode runPreselection(RNode df_);

}

#endif // SELECTIONS_H