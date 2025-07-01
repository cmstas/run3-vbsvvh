#ifndef SELECTIONS_ONELEP2FJ_H
#define SELECTIONS_ONELEP2FJ_H

#include "ROOT/RDataFrame.hxx"
#include "spanet.hpp"
#include "inference.hpp"

#include "commonSelections.h"
#include "utils.h"

using RNode = ROOT::RDF::RNode;

namespace OneLep2FJ {
  RNode TriggerSelections(RNode df);
  RNode runPreselection(RNode df_);
  RNode LeptonSelections(RNode df_);
} // OneLep2FJ

#endif // SELECTIONS_H