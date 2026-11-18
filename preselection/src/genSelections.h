#ifndef GENSELECTIONS_H
#define GENSELECTIONS_H

#pragma once

#include "ROOT/RDataFrame.hxx"

using RNode = ROOT::RDF::RNode;

int find_matching_jet(int target_idx, float target_eta, float target_phi, ROOT::RVec<int> already_matched_jet_indices, ROOT::RVec<int> already_matched_fatjet_indices, ROOT::RVec<float> jet_eta, ROOT::RVec<float> jet_phi, ROOT::RVec<float> fatjet_eta, ROOT::RVec<float> fatjet_phi);
int find_matching_fatjet(int target_idx, float target_eta, float target_phi, ROOT::RVec<int> already_matched_jet_indices, ROOT::RVec<int> already_matched_fatjet_indices, ROOT::RVec<float> jet_eta, ROOT::RVec<float> jet_phi, ROOT::RVec<float> fatjet_eta, ROOT::RVec<float> fatjet_phi);

RNode GenSelections(RNode df_);

// Per-fat-jet gen truth for MC, on the full FatJet_ collection (any channel):
//   FatJet_genMatch      0 nothing, 1 prompt W, 2 W from top, 3 top, 4 prompt Z, 5 H->bb
//   FatJet_genDrVqMax    max dR(jet, q) over the candidate boson's two quarks
//   FatJet_genDrTopB     dR(jet, b) for the b from the candidate W's top
//   FatJet_genVqFlavour  0 light, 1 c and no b, 2 at least one b, -1 no candidate
// Needs GenPart_{pdgId,genPartIdxMother,statusFlags,eta,phi}. See fatjet_gen_truth in
// genSelections.cpp for the matching rules.
RNode FatJetGenTruth(RNode df_);

#endif // GENSELECTIONS_H