#ifndef GENSELECTIONS_H
#define GENSELECTIONS_H

#pragma once

#include "ROOT/RDataFrame.hxx"

using RNode = ROOT::RDF::RNode;

int find_matching_jet(int target_idx, float target_eta, float target_phi, ROOT::RVec<int> already_matched_jet_indices, ROOT::RVec<int> already_matched_fatjet_indices, ROOT::RVec<float> jet_eta, ROOT::RVec<float> jet_phi, ROOT::RVec<float> fatjet_eta, ROOT::RVec<float> fatjet_phi);
int find_matching_fatjet(int target_idx, float target_eta, float target_phi, ROOT::RVec<int> already_matched_jet_indices, ROOT::RVec<int> already_matched_fatjet_indices, ROOT::RVec<float> jet_eta, ROOT::RVec<float> jet_phi, ROOT::RVec<float> fatjet_eta, ROOT::RVec<float> fatjet_phi);

RNode GenSelections(RNode df_);

#endif // GENSELECTIONS_H