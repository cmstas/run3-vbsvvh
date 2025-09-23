#include "truth_selections.h"
#include <cassert>

RNode addTruthVariables(RNode df){
  auto df_out = df.Define("dR_Hd1d2", delta_R, {"truthH_daught1_pt", "truthH_daught1_eta", "truthH_daught1_phi", "truthH_daught1_mass", "truthH_daught2_pt", "truthH_daught2_eta", "truthH_daught2_phi", "truthH_daught2_mass"})
                    .Define("dR_V1d1d2", delta_R, {"truthV1_daught1_pt", "truthV1_daught1_eta", "truthV1_daught1_phi", "truthV1_daught1_mass", "truthV1_daught2_pt", "truthV1_daught2_eta", "truthV1_daught2_phi", "truthV1_daught2_mass"})
                    .Define("dR_V2d1d2", delta_R, {"truthV2_daught1_pt", "truthV2_daught1_eta", "truthV2_daught1_phi", "truthV2_daught1_mass", "truthV2_daught2_pt", "truthV2_daught2_eta", "truthV2_daught2_phi", "truthV2_daught2_mass"})

                    .Define("dR_HV1", delta_R, {"truthH_pt", "truthH_eta", "truthH_phi", "truthH_mass", "truthV1_pt", "truthV1_eta", "truthV1_phi", "truthV1_mass"})
                    .Define("dR_HV2", delta_R, {"truthH_pt", "truthH_eta", "truthH_phi", "truthH_mass", "truthV2_pt", "truthV2_eta", "truthV2_phi", "truthV2_mass"})
                    .Define("dR_V1V2", delta_R, {"truthV1_pt", "truthV1_eta", "truthV1_phi", "truthV1_mass", "truthV2_pt", "truthV2_eta", "truthV2_phi", "truthV2_mass"});

  return df_out;
}

// Given the list of boson daughters (partons or matched GenJe), find the matching jet in the collection. 
//     - Matching priority: H (leading daughter -> trailing daughter) -> Leading V -> Trailing V
//     - A jet index cannot be reused
std::vector<int> truthMatch_GenJets(
    float truthH_daught1_pt, float truthH_daught1_eta, float truthH_daught1_phi, float truthH_daught1_mass, int truthH_daught1_pdgId,
    float truthH_daught2_pt, float truthH_daught2_eta, float truthH_daught2_phi, float truthH_daught2_mass, int truthH_daught2_pdgId,
    float truthV1_daught1_pt, float truthV1_daught1_eta, float truthV1_daught1_phi, float truthV1_daught1_mass, int truthV1_daught1_pdgId,
    float truthV1_daught2_pt, float truthV1_daught2_eta, float truthV1_daught2_phi, float truthV1_daught2_mass, int truthV1_daught2_pdgId,
    float truthV2_daught1_pt, float truthV2_daught1_eta, float truthV2_daught1_phi, float truthV2_daught1_mass, int truthV2_daught1_pdgId,
    float truthV2_daught2_pt, float truthV2_daught2_eta, float truthV2_daught2_phi, float truthV2_daught2_mass, int truthV2_daught2_pdgId,
    float truthVBSq1_ptOrd_pt, float truthVBSq1_ptOrd_eta, float truthVBSq1_ptOrd_phi, float truthVBSq1_ptOrd_mass, int truthVBSq1_ptOrd_pdgId,
    float truthVBSq2_ptOrd_pt, float truthVBSq2_ptOrd_eta, float truthVBSq2_ptOrd_phi, float truthVBSq2_ptOrd_mass, int truthVBSq2_ptOrd_pdgId,
    const ROOT::RVecF GenJet_pt, const ROOT::RVecF GenJet_eta, const ROOT::RVecF GenJet_phi, const ROOT::RVecF GenJet_mass)
{

  std::vector<int> GenJet_matched_indices {};

  // match H
  int genjet_hdaught1_idx = truthMatchedAk4Idx(truthH_daught1_pt, truthH_daught1_eta, truthH_daught1_phi, truthH_daught1_mass, 
  GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass, GenJet_matched_indices);
  GenJet_matched_indices.push_back(genjet_hdaught1_idx);

  int genjet_hdaught2_idx = truthMatchedAk4Idx(truthH_daught2_pt, truthH_daught2_eta, truthH_daught2_phi, truthH_daught2_mass, GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass, GenJet_matched_indices);
  GenJet_matched_indices.push_back(genjet_hdaught2_idx);

  // match V1
  int genjet_v1daught1_idx = truthMatchedAk4Idx(truthV1_daught1_pt, truthV1_daught1_eta, truthV1_daught1_phi, truthV1_daught1_mass, GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass, GenJet_matched_indices);
  GenJet_matched_indices.push_back(genjet_v1daught1_idx);

  int genjet_v1daught2_idx = truthMatchedAk4Idx(truthV1_daught2_pt, truthV1_daught2_eta, truthV1_daught2_phi, truthV1_daught2_mass, GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass, GenJet_matched_indices);
  GenJet_matched_indices.push_back(genjet_v1daught2_idx);

  // match V2
  int genjet_v2daught1_idx = truthMatchedAk4Idx(truthV2_daught1_pt, truthV2_daught1_eta, truthV2_daught1_phi, truthV2_daught1_mass, GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass, GenJet_matched_indices);
  GenJet_matched_indices.push_back(genjet_v2daught1_idx);

  int genjet_v2daught2_idx = truthMatchedAk4Idx(truthV2_daught2_pt, truthV2_daught2_eta, truthV2_daught2_phi, truthV2_daught2_mass, GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass, GenJet_matched_indices);
  GenJet_matched_indices.push_back(genjet_v2daught2_idx);

  // match leading vbs
  float smallest_dR = 0.4;
  int genjet_vbs1_idx = truthMatchedAk4Idx(truthVBSq1_ptOrd_pt, truthVBSq1_ptOrd_eta, truthVBSq1_ptOrd_phi, truthVBSq1_ptOrd_mass,
                                           GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass, GenJet_matched_indices, smallest_dR = smallest_dR);
  GenJet_matched_indices.push_back(genjet_vbs1_idx);

  // match trailing vbs
  int genjet_vbs2_idx = truthMatchedAk4Idx(truthVBSq2_ptOrd_pt, truthVBSq2_ptOrd_eta, truthVBSq2_ptOrd_phi, truthVBSq2_ptOrd_mass, 
  GenJet_pt, GenJet_eta, GenJet_phi, GenJet_mass, GenJet_matched_indices, smallest_dR = smallest_dR);

  GenJet_matched_indices.push_back(genjet_vbs2_idx);

  return GenJet_matched_indices;

}


// Given the list of bosons (partons or matched GenJetAK8), boson daughters (partons only), find the matching jet in the collection.
//     - Matching priority: H -> Leading V -> Trailing V
//     - A jet index cannot be reused  
std::vector<int> truthMatch_GenJetsAK8_HV1V2(
    float truthH_pt, float truthH_eta, float truthH_phi, float truthH_mass,
    float truthH_daught1_pt, float truthH_daught1_eta, float truthH_daught1_phi, float truthH_daught1_mass, int truthH_daught1_pdgId,
    float truthH_daught2_pt, float truthH_daught2_eta, float truthH_daught2_phi, float truthH_daught2_mass, int truthH_daught2_pdgId,
    float truthV1_pt, float truthV1_eta, float truthV1_phi, float truthV1_mass,
    float truthV1_daught1_pt, float truthV1_daught1_eta, float truthV1_daught1_phi, float truthV1_daught1_mass, int truthV1_daught1_pdgId,
    float truthV1_daught2_pt, float truthV1_daught2_eta, float truthV1_daught2_phi, float truthV1_daught2_mass, int truthV1_daught2_pdgId,
    float truthV2_pt, float truthV2_eta, float truthV2_phi, float truthV2_mass,
    float truthV2_daught1_pt, float truthV2_daught1_eta, float truthV2_daught1_phi, float truthV2_daught1_mass, int truthV2_daught1_pdgId,
    float truthV2_daught2_pt, float truthV2_daught2_eta, float truthV2_daught2_phi, float truthV2_daught2_mass, int truthV2_daught2_pdgId,
    const ROOT::RVecF GenJetAK8_pt, const ROOT::RVecF GenJetAK8_eta, const ROOT::RVecF GenJetAK8_phi, const ROOT::RVecF GenJetAK8_mass)
{
  std::vector<int> GenJetAK8_matched_indices {};

  // match H
  int genjetak8_h_idx = -1;
  if (isQuarkOrGluon(truthH_daught1_pdgId) && isQuarkOrGluon(truthH_daught2_pdgId)) {
    genjetak8_h_idx = truthMatchedAk8Idx_daughtContainement(truthH_pt, truthH_eta, truthH_phi, truthH_mass,
                                           truthH_daught1_pt, truthH_daught1_eta, truthH_daught1_phi, truthH_daught1_mass,
                                           truthH_daught2_pt, truthH_daught2_eta, truthH_daught2_phi, truthH_daught2_mass,
                                           GenJetAK8_pt, GenJetAK8_eta, GenJetAK8_phi, GenJetAK8_mass, GenJetAK8_matched_indices);
  }                                          
  GenJetAK8_matched_indices.push_back(genjetak8_h_idx);

  // match V1
  int genjetak8_v1_idx = -1;
  if (isQuarkOrGluon(truthV1_daught1_pdgId) && isQuarkOrGluon(truthV1_daught2_pdgId))
  {
    genjetak8_v1_idx = truthMatchedAk8Idx_daughtContainement(truthV1_pt, truthV1_eta, truthV1_phi, truthV1_mass,
                                          truthV1_daught1_pt, truthV1_daught1_eta, truthV1_daught1_phi, truthV1_daught1_mass,
                                          truthV1_daught2_pt, truthV1_daught2_eta, truthV1_daught2_phi, truthV1_daught2_mass,
                                          GenJetAK8_pt, GenJetAK8_eta, GenJetAK8_phi, GenJetAK8_mass, GenJetAK8_matched_indices);
  }
  GenJetAK8_matched_indices.push_back(genjetak8_v1_idx);

  // match V2
  int genjetak8_v2_idx = -1;
  if (isQuarkOrGluon(truthV2_daught1_pdgId) && isQuarkOrGluon(truthV2_daught2_pdgId))
  {
    genjetak8_v2_idx = truthMatchedAk8Idx_daughtContainement(truthV2_pt, truthV2_eta, truthV2_phi, truthV2_mass,
                                            truthV2_daught1_pt, truthV2_daught1_eta, truthV2_daught1_phi, truthV2_daught1_mass,
                                            truthV2_daught2_pt, truthV2_daught2_eta, truthV2_daught2_phi, truthV2_daught2_mass,
                                            GenJetAK8_pt, GenJetAK8_eta, GenJetAK8_phi, GenJetAK8_mass, GenJetAK8_matched_indices);
  }
  GenJetAK8_matched_indices.push_back(genjetak8_v2_idx);

  return GenJetAK8_matched_indices;
}

// Find the AK8 jet matched to the given GenJet. Return indices of matched Ak8 jets to H, V1, V2 in order.
std::vector<int> truthMatch_RecoGenAK8_HV1V2(
    float genH_pt, float genH_eta, float genH_phi, float genH_mass,
    float genV1_pt, float genV1_eta, float genV1_phi, float genV1_mass,
    float genV2_pt, float genV2_eta, float genV2_phi, float genV2_mass,
    const ROOT::RVecF AK8Jet_pt, const ROOT::RVecF AK8Jet_eta, const ROOT::RVecF AK8Jet_phi, const ROOT::RVecF AK8Jet_mass)
{
  std::vector<int> AK8Jet_matched_indices{};
  float matching_dR = 0.8;

  // match H

  int h_idx = -1;
  if (genH_pt>=0) {
    h_idx = truthMatchedJetArrIdx(genH_pt, genH_eta, genH_phi, genH_mass,
                               AK8Jet_pt, AK8Jet_eta, AK8Jet_phi, AK8Jet_mass,
                               AK8Jet_matched_indices, matching_dR);
  }
  AK8Jet_matched_indices.push_back(h_idx);

  // match V1
  int v1_idx = -1;
  if (genV1_pt>=0) {
    v1_idx = truthMatchedJetArrIdx(genV1_pt, genV1_eta, genV1_phi, genV1_mass,
                                AK8Jet_pt, AK8Jet_eta, AK8Jet_phi, AK8Jet_mass,
                                AK8Jet_matched_indices, matching_dR);
  }
  AK8Jet_matched_indices.push_back(v1_idx);

  // match V2
  int v2_idx = -1;
  if (genV2_pt >= 0){
    v2_idx = truthMatchedJetArrIdx(genV2_pt, genV2_eta, genV2_phi, genV2_mass,
                                AK8Jet_pt, AK8Jet_eta, AK8Jet_phi, AK8Jet_mass,
                                AK8Jet_matched_indices, matching_dR);
  }
  AK8Jet_matched_indices.push_back(v2_idx);

  return AK8Jet_matched_indices;
}

RNode reconstructTruthEvent(RNode df)
{
  std::cout << " -> Run reconstructTruthEvent()" << std::endl;

  // Order truth VBS jets by pt and store new variables for bookkeeping
  auto df_out = df.Define("VBSJetsAlreadyPtOrdered", [](float pt1, float pt2)
                          {
                              if (pt1 >= pt2) {
                                return true;
                              }
                              return false; },
                          {"truthVBSq1_pt", "truthVBSq2_pt"})
                    .Define("truthVBSq1_ptOrd_pt", "VBSJetsAlreadyPtOrdered ? truthVBSq1_pt : truthVBSq2_pt")
                    .Define("truthVBSq1_ptOrd_eta", "VBSJetsAlreadyPtOrdered ? truthVBSq1_eta : truthVBSq2_eta")
                    .Define("truthVBSq1_ptOrd_phi", "VBSJetsAlreadyPtOrdered ? truthVBSq1_phi : truthVBSq2_phi")
                    .Define("truthVBSq1_ptOrd_mass", "VBSJetsAlreadyPtOrdered ? truthVBSq1_mass : truthVBSq2_mass")
                    .Define("truthVBSq1_ptOrd_pdgId", "VBSJetsAlreadyPtOrdered ? truthVBSq1_pdgId : truthVBSq2_pdgId")
                    .Define("truthVBSq2_ptOrd_pt", "VBSJetsAlreadyPtOrdered ? truthVBSq2_pt : truthVBSq1_pt")
                    .Define("truthVBSq2_ptOrd_eta", "VBSJetsAlreadyPtOrdered ? truthVBSq2_eta : truthVBSq1_eta")
                    .Define("truthVBSq2_ptOrd_phi", "VBSJetsAlreadyPtOrdered ? truthVBSq2_phi : truthVBSq1_phi")
                    .Define("truthVBSq2_ptOrd_mass", "VBSJetsAlreadyPtOrdered ? truthVBSq2_mass : truthVBSq1_mass")
                    .Define("truthVBSq2_ptOrd_pdgId", "VBSJetsAlreadyPtOrdered ? truthVBSq2_pdgId : truthVBSq1_pdgId");

  // Find GenJet matched to H/V1/V2 daughters and VBS jets
  // - Matching priority: H daughters (leading->trailing) -> V1 daughters -> V2 daughters -> vbs1 -> vbs2
  df_out = df_out.Define("GenJet_matchedIndices", truthMatch_GenJets,
                         {"truthH_daught1_pt", "truthH_daught1_eta", "truthH_daught1_phi", "truthH_daught1_mass", "truthH_daught1_pdgId",
                          "truthH_daught2_pt", "truthH_daught2_eta", "truthH_daught2_phi", "truthH_daught2_mass", "truthH_daught2_pdgId",
                          "truthV1_daught1_pt", "truthV1_daught1_eta", "truthV1_daught1_phi", "truthV1_daught1_mass", "truthV1_daught1_pdgId",
                          "truthV1_daught2_pt", "truthV1_daught2_eta", "truthV1_daught2_phi", "truthV1_daught2_mass", "truthV1_daught2_pdgId",
                          "truthV2_daught1_pt", "truthV2_daught1_eta", "truthV2_daught1_phi", "truthV2_daught1_mass", "truthV2_daught1_pdgId",
                          "truthV2_daught2_pt", "truthV2_daught2_eta", "truthV2_daught2_phi", "truthV2_daught2_mass", "truthV2_daught2_pdgId",
                          "truthVBSq1_ptOrd_pt", "truthVBSq1_ptOrd_eta", "truthVBSq1_ptOrd_phi", "truthVBSq1_ptOrd_mass", "truthVBSq1_ptOrd_pdgId",
                          "truthVBSq2_ptOrd_pt", "truthVBSq2_ptOrd_eta", "truthVBSq2_ptOrd_phi", "truthVBSq2_ptOrd_mass", "truthVBSq2_ptOrd_pdgId",
                          "GenJet_pt", "GenJet_eta", "GenJet_phi", "GenJet_mass"});

  // Set matched GenJet variables
  df_out = df_out.Define("GenJet_Hdaught1_pt", "GenJet_matchedIndices[0] > -1 ? GenJet_pt[GenJet_matchedIndices[0]] : -99")
               .Define("GenJet_Hdaught1_eta", "GenJet_matchedIndices[0] > -1 ? GenJet_eta[GenJet_matchedIndices[0]] : -99")
               .Define("GenJet_Hdaught1_phi", "GenJet_matchedIndices[0] > -1 ? GenJet_phi[GenJet_matchedIndices[0]] : -99")
               .Define("GenJet_Hdaught1_mass", "GenJet_matchedIndices[0] > -1 ? GenJet_mass[GenJet_matchedIndices[0]] : -99")

               .Define("GenJet_Hdaught2_pt", "GenJet_matchedIndices[1] > -1 ? GenJet_pt[GenJet_matchedIndices[1]] : -99")
               .Define("GenJet_Hdaught2_eta", "GenJet_matchedIndices[1] > -1 ? GenJet_eta[GenJet_matchedIndices[1]] : -99")
               .Define("GenJet_Hdaught2_phi", "GenJet_matchedIndices[1] > -1 ? GenJet_phi[GenJet_matchedIndices[1]] : -99")
               .Define("GenJet_Hdaught2_mass", "GenJet_matchedIndices[1] > -1 ? GenJet_mass[GenJet_matchedIndices[1]] : -99")

               .Define("GenJet_V1daught1_pt", "GenJet_matchedIndices[2] > -1 ? GenJet_pt[GenJet_matchedIndices[2]] : -99")
               .Define("GenJet_V1daught1_eta", "GenJet_matchedIndices[2] > -1 ? GenJet_eta[GenJet_matchedIndices[2]] : -99")
               .Define("GenJet_V1daught1_phi", "GenJet_matchedIndices[2] > -1 ? GenJet_phi[GenJet_matchedIndices[2]] : -99")
               .Define("GenJet_V1daught1_mass", "GenJet_matchedIndices[2] > -1 ? GenJet_mass[GenJet_matchedIndices[2]] : -99")

               .Define("GenJet_V1daught2_pt", "GenJet_matchedIndices[3] > -1 ? GenJet_pt[GenJet_matchedIndices[3]] : -99")
               .Define("GenJet_V1daught2_eta", "GenJet_matchedIndices[3] > -1 ? GenJet_eta[GenJet_matchedIndices[3]] : -99")
               .Define("GenJet_V1daught2_phi", "GenJet_matchedIndices[3] > -1 ? GenJet_phi[GenJet_matchedIndices[3]] : -99")
               .Define("GenJet_V1daught2_mass", "GenJet_matchedIndices[3] > -1 ? GenJet_mass[GenJet_matchedIndices[3]] : -99")

               .Define("GenJet_V2daught1_pt", "GenJet_matchedIndices[4] > -1 ? GenJet_pt[GenJet_matchedIndices[4]] : -99")
               .Define("GenJet_V2daught1_eta", "GenJet_matchedIndices[4] > -1 ? GenJet_eta[GenJet_matchedIndices[4]] : -99")
               .Define("GenJet_V2daught1_phi", "GenJet_matchedIndices[4] > -1 ? GenJet_phi[GenJet_matchedIndices[4]] : -99")
               .Define("GenJet_V2daught1_mass", "GenJet_matchedIndices[4] > -1 ? GenJet_mass[GenJet_matchedIndices[4]] : -99")

               .Define("GenJet_V2daught2_pt", "GenJet_matchedIndices[5] > -1 ? GenJet_pt[GenJet_matchedIndices[5]] : -99")
               .Define("GenJet_V2daught2_eta", "GenJet_matchedIndices[5] > -1 ? GenJet_eta[GenJet_matchedIndices[5]] : -99")
               .Define("GenJet_V2daught2_phi", "GenJet_matchedIndices[5] > -1 ? GenJet_phi[GenJet_matchedIndices[5]] : -99")
               .Define("GenJet_V2daught2_mass", "GenJet_matchedIndices[5] > -1 ? GenJet_mass[GenJet_matchedIndices[5]] : -99")

               .Define("GenJet_vbs1_pt", "GenJet_matchedIndices[6] > -1 ? GenJet_pt[GenJet_matchedIndices[6]] : -99")
               .Define("GenJet_vbs1_eta", "GenJet_matchedIndices[6] > -1 ? GenJet_eta[GenJet_matchedIndices[6]] : -99")
               .Define("GenJet_vbs1_phi", "GenJet_matchedIndices[6] > -1 ? GenJet_phi[GenJet_matchedIndices[6]] : -99")
               .Define("GenJet_vbs1_mass", "GenJet_matchedIndices[6] > -1 ? GenJet_mass[GenJet_matchedIndices[6]] : -99")

               .Define("GenJet_vbs2_pt", "GenJet_matchedIndices[7] > -1 ? GenJet_pt[GenJet_matchedIndices[7]] : -99")
               .Define("GenJet_vbs2_eta", "GenJet_matchedIndices[7] > -1 ? GenJet_eta[GenJet_matchedIndices[7]] : -99")
               .Define("GenJet_vbs2_phi", "GenJet_matchedIndices[7] > -1 ? GenJet_phi[GenJet_matchedIndices[7]] : -99")
               .Define("GenJet_vbs2_mass", "GenJet_matchedIndices[7] > -1 ? GenJet_mass[GenJet_matchedIndices[7]] : -99");
  

  // Find GenJetAK8 matched to H/V1/V2, requiring containment of truth daughters
  df_out = df_out.Define("GenJetAK8_matchedIndices", truthMatch_GenJetsAK8_HV1V2,
                         {"truthH_pt", "truthH_eta", "truthH_phi", "truthH_mass",
                          "truthH_daught1_pt", "truthH_daught1_eta", "truthH_daught1_phi", "truthH_daught1_mass", "truthH_daught1_pdgId",
                          "truthH_daught2_pt", "truthH_daught2_eta", "truthH_daught2_phi", "truthH_daught2_mass", "truthH_daught2_pdgId",
                          "truthV1_pt", "truthV1_eta", "truthV1_phi", "truthV1_mass",
                          "truthV1_daught1_pt", "truthV1_daught1_eta", "truthV1_daught1_phi", "truthV1_daught1_mass", "truthV1_daught1_pdgId",
                          "truthV1_daught2_pt", "truthV1_daught2_eta", "truthV1_daught2_phi", "truthV1_daught2_mass", "truthV1_daught2_pdgId",
                          "truthV2_pt", "truthV2_eta", "truthV2_phi", "truthV2_mass",
                          "truthV2_daught1_pt", "truthV2_daught1_eta", "truthV2_daught1_phi", "truthV2_daught1_mass", "truthV2_daught1_pdgId",
                          "truthV2_daught2_pt", "truthV2_daught2_eta", "truthV2_daught2_phi", "truthV2_daught2_mass", "truthV2_daught2_pdgId",
                          "GenJetAK8_pt", "GenJetAK8_eta", "GenJetAK8_phi", "GenJetAK8_mass"});

  // Set matched GenJetAK8 variables
  df_out = df_out.Define("GenJetAK8_H_pt", "GenJetAK8_matchedIndices[0] > -1 ? GenJetAK8_pt[GenJetAK8_matchedIndices[0]] : -99")
               .Define("GenJetAK8_H_eta", "GenJetAK8_matchedIndices[0] > -1 ? GenJetAK8_eta[GenJetAK8_matchedIndices[0]] : -99")
               .Define("GenJetAK8_H_phi", "GenJetAK8_matchedIndices[0] > -1 ? GenJetAK8_phi[GenJetAK8_matchedIndices[0]] : -99")
               .Define("GenJetAK8_H_mass", "GenJetAK8_matchedIndices[0] > -1 ? GenJetAK8_mass[GenJetAK8_matchedIndices[0]] : -99")

               .Define("GenJetAK8_V1_pt", "GenJetAK8_matchedIndices[1] > -1 ? GenJetAK8_pt[GenJetAK8_matchedIndices[1]] : -99")
               .Define("GenJetAK8_V1_eta", "GenJetAK8_matchedIndices[1] > -1 ? GenJetAK8_eta[GenJetAK8_matchedIndices[1]] : -99")
               .Define("GenJetAK8_V1_phi", "GenJetAK8_matchedIndices[1] > -1 ? GenJetAK8_phi[GenJetAK8_matchedIndices[1]] : -99")
               .Define("GenJetAK8_V1_mass", "GenJetAK8_matchedIndices[1] > -1 ? GenJetAK8_mass[GenJetAK8_matchedIndices[1]] : -99")

               .Define("GenJetAK8_V2_pt", "GenJetAK8_matchedIndices[2] > -1 ? GenJetAK8_pt[GenJetAK8_matchedIndices[2]] : -99")
               .Define("GenJetAK8_V2_eta", "GenJetAK8_matchedIndices[2] > -1 ? GenJetAK8_eta[GenJetAK8_matchedIndices[2]] : -99")
               .Define("GenJetAK8_V2_phi", "GenJetAK8_matchedIndices[2] > -1 ? GenJetAK8_phi[GenJetAK8_matchedIndices[2]] : -99")
               .Define("GenJetAK8_V2_mass", "GenJetAK8_matchedIndices[2] > -1 ? GenJetAK8_mass[GenJetAK8_matchedIndices[2]] : -99");


  // Add some variables to check matching
  df_out = df_out.Define("dR_Hd1_truth_matchedGenJet", delta_R, {"truthH_daught1_pt", "truthH_daught1_eta", "truthH_daught1_phi", "truthH_daught1_mass", "GenJet_Hdaught1_pt", "GenJet_Hdaught1_eta", "GenJet_Hdaught1_phi", "GenJet_Hdaught1_mass"})
               .Define("dR_Hd2_truth_matchedGenJet", delta_R, {"truthH_daught2_pt", "truthH_daught2_eta", "truthH_daught2_phi", "truthH_daught2_mass", "GenJet_Hdaught2_pt", "GenJet_Hdaught2_eta", "GenJet_Hdaught2_phi", "GenJet_Hdaught2_mass"})
               .Define("dR_V1d1_truth_matchedGenJet", delta_R, {"truthV1_daught1_pt", "truthV1_daught1_eta", "truthV1_daught1_phi", "truthV1_daught1_mass", "GenJet_V1daught1_pt", "GenJet_V1daught1_eta", "GenJet_V1daught1_phi", "GenJet_V1daught1_mass"})
               .Define("dR_V1d2_truth_matchedGenJet", delta_R, {"truthV1_daught2_pt", "truthV1_daught2_eta", "truthV1_daught2_phi", "truthV1_daught2_mass", "GenJet_V1daught2_pt", "GenJet_V1daught2_eta", "GenJet_V1daught2_phi", "GenJet_V1daught2_mass"})
               .Define("dR_V2d1_truth_matchedGenJet", delta_R, {"truthV2_daught1_pt", "truthV2_daught1_eta", "truthV2_daught1_phi", "truthV2_daught1_mass", "GenJet_V2daught1_pt", "GenJet_V2daught1_eta", "GenJet_V2daught1_phi", "GenJet_V2daught1_mass"})
               .Define("dR_V2d2_truth_matchedGenJet", delta_R, {"truthV2_daught2_pt", "truthV2_daught2_eta", "truthV2_daught2_phi", "truthV2_daught2_mass", "GenJet_V2daught2_pt", "GenJet_V2daught2_eta", "GenJet_V2daught2_phi", "GenJet_V2daught2_mass"})
               .Define("dR_Vbs1_truth_matchedGenJet", delta_R, {"truthVBSq1_ptOrd_pt", "truthVBSq1_ptOrd_eta", "truthVBSq1_ptOrd_phi", "truthVBSq1_ptOrd_mass", "GenJet_vbs1_pt", "GenJet_vbs1_eta", "GenJet_vbs1_phi", "GenJet_vbs1_mass"})
               .Define("dR_Vbs2_truth_matchedGenJet", delta_R, {"truthVBSq2_ptOrd_pt", "truthVBSq2_ptOrd_eta", "truthVBSq2_ptOrd_phi", "truthVBSq2_ptOrd_mass", "GenJet_vbs2_pt", "GenJet_vbs2_eta", "GenJet_vbs2_phi", "GenJet_vbs2_mass"})
               .Define("dR_H_truth_matchedGenJet", delta_R, {"truthH_pt", "truthH_eta", "truthH_phi", "truthH_mass", "GenJetAK8_H_pt", "GenJetAK8_H_eta", "GenJetAK8_H_phi", "GenJetAK8_H_mass"})
               .Define("dR_V1_truth_matchedGenJet", delta_R, {"truthV1_pt", "truthV1_eta", "truthV1_phi", "truthV1_mass", "GenJetAK8_V1_pt", "GenJetAK8_V1_eta", "GenJetAK8_V1_phi", "GenJetAK8_V1_mass"})
               .Define("dR_V2_truth_matchedGenJet", delta_R, {"truthV2_pt", "truthV2_eta", "truthV2_phi", "truthV2_mass", "GenJetAK8_V2_pt", "GenJetAK8_V2_eta", "GenJetAK8_V2_phi", "GenJetAK8_V2_mass"});

  return df_out;
}

RNode reconstructRecoEventWithTruthInfo_Resolved(RNode df, std::string recoJetsCollection, std::string recoJetsCollectionAK8)
{
  std::cout << " -> Run reconstructRecoEventWithTruthInfo_Resolved(recoJetsCollection = " << recoJetsCollection << ")" << std::endl;

  // Find reco ak4 jets matched to GenJets matched to bosons daughters
  // matching priority: Higgs (leading daughter -> trailing daughter) -> leading V -> trailing V

  std::vector recoJetsVars = { recoJetsCollection + "_pt", // order matters
                               recoJetsCollection + "_eta",
                               recoJetsCollection + "_phi",
                               recoJetsCollection + "_mass" };

  std::vector recoJetsVarsAK8 = {recoJetsCollectionAK8 + "_pt", // order matters
                              recoJetsCollectionAK8 + "_eta",
                              recoJetsCollectionAK8 + "_phi",
                              recoJetsCollectionAK8 + "_mass"};
 
  auto df_out = df.Define(recoJetsCollection + "_matchedIndices", truthMatchGenJetToRecoWithoutAK8overlap,
                          {"GenJet_matchedIndices", "GenJet_pt", "GenJet_eta", "GenJet_phi", "GenJet_mass",
                           recoJetsVars.at(0), recoJetsVars.at(1), recoJetsVars.at(2), recoJetsVars.at(3),
                           recoJetsCollectionAK8 + "_matchedIndices",
                           recoJetsVarsAK8.at(0), recoJetsVarsAK8.at(1), recoJetsVarsAK8.at(2), recoJetsVarsAK8.at(3)});

  df_out = df_out.Define("hq1_idx_" + recoJetsCollection, recoJetsCollection + "_matchedIndices[0]")
               .Define("hq2_idx_" + recoJetsCollection, recoJetsCollection + "_matchedIndices[1]")
               .Define("v1q1_idx_" + recoJetsCollection, recoJetsCollection + "_matchedIndices[2]")
               .Define("v1q2_idx_" + recoJetsCollection, recoJetsCollection + "_matchedIndices[3]")
               .Define("v2q1_idx_" + recoJetsCollection, recoJetsCollection + "_matchedIndices[4]")
               .Define("v2q2_idx_" + recoJetsCollection, recoJetsCollection + "_matchedIndices[5]")
               .Define("vbs1_idx_" + recoJetsCollection, recoJetsCollection + "_matchedIndices[6]")
               .Define("vbs2_idx_" + recoJetsCollection, recoJetsCollection + "_matchedIndices[7]");


  // Old matching: match truth partons directly to AK4 jets (keeping for reference, can delete eventually)
  // >>
  df_out = df_out.Define(recoJetsCollection + "_oldMatchedIndices", truthMatch_GenJets,
                         {"truthH_daught1_pt", "truthH_daught1_eta", "truthH_daught1_phi", "truthH_daught1_mass", "truthH_daught1_pdgId",
                          "truthH_daught2_pt", "truthH_daught2_eta", "truthH_daught2_phi", "truthH_daught2_mass", "truthH_daught2_pdgId",
                          "truthV1_daught1_pt", "truthV1_daught1_eta", "truthV1_daught1_phi", "truthV1_daught1_mass", "truthV1_daught1_pdgId",
                          "truthV1_daught2_pt", "truthV1_daught2_eta", "truthV1_daught2_phi", "truthV1_daught2_mass", "truthV1_daught2_pdgId",
                          "truthV2_daught1_pt", "truthV2_daught1_eta", "truthV2_daught1_phi", "truthV2_daught1_mass", "truthV2_daught1_pdgId",
                          "truthV2_daught2_pt", "truthV2_daught2_eta", "truthV2_daught2_phi", "truthV2_daught2_mass", "truthV2_daught2_pdgId",
                          "truthVBSq1_ptOrd_pt", "truthVBSq1_ptOrd_eta", "truthVBSq1_ptOrd_phi", "truthVBSq1_ptOrd_mass", "truthVBSq1_ptOrd_pdgId",
                          "truthVBSq2_ptOrd_pt", "truthVBSq2_ptOrd_eta", "truthVBSq2_ptOrd_phi", "truthVBSq2_ptOrd_mass", "truthVBSq2_ptOrd_pdgId",
                          recoJetsVars.at(0), recoJetsVars.at(1), recoJetsVars.at(2), recoJetsVars.at(3)});

  df_out = df_out.Define("hq1_idx_oldMatching_" + recoJetsCollection, recoJetsCollection + "_oldMatchedIndices[0]")
               .Define("hq2_idx_oldMatching_" + recoJetsCollection, recoJetsCollection + "_oldMatchedIndices[1]")
               .Define("v1q1_idx_oldMatching_" + recoJetsCollection, recoJetsCollection + "_oldMatchedIndices[2]")
               .Define("v1q2_idx_oldMatching_" + recoJetsCollection, recoJetsCollection + "_oldMatchedIndices[3]")
               .Define("v2q1_idx_oldMatching_" + recoJetsCollection, recoJetsCollection + "_oldMatchedIndices[4]")
               .Define("v2q2_idx_oldMatching_" + recoJetsCollection, recoJetsCollection + "_oldMatchedIndices[5]")
               .Define("vbs1_idx_oldMatching_" + recoJetsCollection, recoJetsCollection + "_oldMatchedIndices[6]")
               .Define("vbs2_idx_oldMatching_" + recoJetsCollection, recoJetsCollection + "_oldMatchedIndices[7]");
  // << 

  // Define some variables to check the matching

  // Matching deltaR
  df_out = df_out.Define("dR_Hd1_truth_matched"+recoJetsCollection, dR_particle_ak4AtIdx, {"truthH_daught1_pt", "truthH_daught1_eta", "truthH_daught1_phi", "truthH_daught1_mass", recoJetsVars.at(0), recoJetsVars.at(1), recoJetsVars.at(2), recoJetsVars.at(3), "hq1_idx_"+recoJetsCollection})
               .Define("dR_Hd2_truth_matched"+recoJetsCollection, dR_particle_ak4AtIdx, {"truthH_daught2_pt", "truthH_daught2_eta", "truthH_daught2_phi", "truthH_daught2_mass", recoJetsVars.at(0), recoJetsVars.at(1), recoJetsVars.at(2), recoJetsVars.at(3), "hq2_idx_"+recoJetsCollection})
               .Define("dR_V1d1_truth_matched"+recoJetsCollection, dR_particle_ak4AtIdx, {"truthV1_daught1_pt", "truthV1_daught1_eta", "truthV1_daught1_phi", "truthV1_daught1_mass", recoJetsVars.at(0), recoJetsVars.at(1), recoJetsVars.at(2), recoJetsVars.at(3), "v1q1_idx_"+recoJetsCollection})
               .Define("dR_V1d2_truth_matched"+recoJetsCollection, dR_particle_ak4AtIdx, {"truthV1_daught2_pt", "truthV1_daught2_eta", "truthV1_daught2_phi", "truthV1_daught2_mass", recoJetsVars.at(0), recoJetsVars.at(1), recoJetsVars.at(2), recoJetsVars.at(3), "v1q2_idx_"+recoJetsCollection})
               .Define("dR_V2d1_truth_matched"+recoJetsCollection, dR_particle_ak4AtIdx, {"truthV2_daught1_pt", "truthV2_daught1_eta", "truthV2_daught1_phi", "truthV2_daught1_mass", recoJetsVars.at(0), recoJetsVars.at(1), recoJetsVars.at(2), recoJetsVars.at(3), "v2q1_idx_"+recoJetsCollection})
               .Define("dR_V2d2_truth_matched"+recoJetsCollection, dR_particle_ak4AtIdx, {"truthV2_daught2_pt", "truthV2_daught2_eta", "truthV2_daught2_phi", "truthV2_daught2_mass", recoJetsVars.at(0), recoJetsVars.at(1), recoJetsVars.at(2), recoJetsVars.at(3), "v2q2_idx_"+recoJetsCollection})
               .Define("dR_Vbs1_truth_matched"+recoJetsCollection, dR_particle_ak4AtIdx, {"truthVBSq1_ptOrd_pt", "truthVBSq1_ptOrd_eta", "truthVBSq1_ptOrd_phi", "truthVBSq1_ptOrd_mass", recoJetsVars.at(0), recoJetsVars.at(1), recoJetsVars.at(2), recoJetsVars.at(3), "vbs1_idx_"+recoJetsCollection})
               .Define("dR_Vbs2_truth_matched"+recoJetsCollection, dR_particle_ak4AtIdx, {"truthVBSq2_ptOrd_pt", "truthVBSq2_ptOrd_eta", "truthVBSq2_ptOrd_phi", "truthVBSq2_ptOrd_mass", recoJetsVars.at(0), recoJetsVars.at(1), recoJetsVars.at(2), recoJetsVars.at(3), "vbs2_idx_"+recoJetsCollection})

               .Define("dR_Hd1_matchedGenJet_matched" + recoJetsCollection, dR_particle_ak4AtIdx, {"GenJet_Hdaught1_pt", "GenJet_Hdaught1_eta", "GenJet_Hdaught1_phi", "GenJet_Hdaught1_mass", recoJetsVars.at(0), recoJetsVars.at(1), recoJetsVars.at(2), recoJetsVars.at(3), "hq1_idx_"+recoJetsCollection})
               .Define("dR_Hd2_matchedGenJet_matched" + recoJetsCollection, dR_particle_ak4AtIdx, {"GenJet_Hdaught2_pt", "GenJet_Hdaught2_eta", "GenJet_Hdaught2_phi", "GenJet_Hdaught2_mass", recoJetsVars.at(0), recoJetsVars.at(1), recoJetsVars.at(2), recoJetsVars.at(3), "hq2_idx_"+recoJetsCollection})
               .Define("dR_V1d1_matchedGenJet_matched" + recoJetsCollection, dR_particle_ak4AtIdx, {"GenJet_V1daught1_pt", "GenJet_V1daught1_eta", "GenJet_V1daught1_phi", "GenJet_V1daught1_mass", recoJetsVars.at(0), recoJetsVars.at(1), recoJetsVars.at(2), recoJetsVars.at(3), "v1q1_idx_"+recoJetsCollection})
               .Define("dR_V1d2_matchedGenJet_matched" + recoJetsCollection, dR_particle_ak4AtIdx, {"GenJet_V1daught2_pt", "GenJet_V1daught2_eta", "GenJet_V1daught2_phi", "GenJet_V1daught2_mass", recoJetsVars.at(0), recoJetsVars.at(1), recoJetsVars.at(2), recoJetsVars.at(3), "v1q2_idx_"+recoJetsCollection})
               .Define("dR_V2d1_matchedGenJet_matched" + recoJetsCollection, dR_particle_ak4AtIdx, {"GenJet_V2daught1_pt", "GenJet_V2daught1_eta", "GenJet_V2daught1_phi", "GenJet_V2daught1_mass", recoJetsVars.at(0), recoJetsVars.at(1), recoJetsVars.at(2), recoJetsVars.at(3), "v2q1_idx_"+recoJetsCollection})
               .Define("dR_V2d2_matchedGenJet_matched" + recoJetsCollection, dR_particle_ak4AtIdx, {"GenJet_V2daught2_pt", "GenJet_V2daught2_eta", "GenJet_V2daught2_phi", "GenJet_V2daught2_mass", recoJetsVars.at(0), recoJetsVars.at(1), recoJetsVars.at(2), recoJetsVars.at(3), "v2q2_idx_"+recoJetsCollection})
               .Define("dR_Vbs1_matchedGenJet_matched" + recoJetsCollection, dR_particle_ak4AtIdx, {"truthVBSq1_ptOrd_pt", "truthVBSq1_ptOrd_eta", "truthVBSq1_ptOrd_phi", "truthVBSq1_ptOrd_mass", recoJetsVars.at(0), recoJetsVars.at(1), recoJetsVars.at(2), recoJetsVars.at(3), "vbs1_idx_"+recoJetsCollection})
               .Define("dR_Vbs2_matchedGenJet_matched" + recoJetsCollection, dR_particle_ak4AtIdx, {"truthVBSq2_ptOrd_pt", "truthVBSq2_ptOrd_eta", "truthVBSq2_ptOrd_phi", "truthVBSq2_ptOrd_mass", recoJetsVars.at(0), recoJetsVars.at(1), recoJetsVars.at(2), recoJetsVars.at(3), "vbs2_idx_"+recoJetsCollection});

  // Reconstructed dijet invariant mass
  df_out = df_out.Define("H_recoDijet_" + recoJetsCollection + "_m", dijetInvMass_fromIdx, {recoJetsVars.at(0), recoJetsVars.at(1), recoJetsVars.at(2), recoJetsVars.at(3), "hq1_idx_" + recoJetsCollection, "hq2_idx_" + recoJetsCollection})
               .Define("V1_recoDijet_" + recoJetsCollection + "_m", dijetInvMass_fromIdx, {recoJetsVars.at(0), recoJetsVars.at(1), recoJetsVars.at(2), recoJetsVars.at(3), "v1q1_idx_" + recoJetsCollection, "v1q2_idx_" + recoJetsCollection})
               .Define("V2_recoDijet_" + recoJetsCollection + "_m", dijetInvMass_fromIdx, {recoJetsVars.at(0), recoJetsVars.at(1), recoJetsVars.at(2), recoJetsVars.at(3), "v2q1_idx_" + recoJetsCollection, "v2q2_idx_" + recoJetsCollection});

  // Count number of quarks matched to a reco jet
  df_out = df_out.Define("nTruthMatchedBosonDaught_" + recoJetsCollection, [](const int &hq1, const int &hq2, const int &v1q1, const int &v1q2, const int &v2q1, const int &v2q2) -> int {
                          std::vector<int> values = {hq1, hq2, v1q1, v1q2, v2q1, v2q2};
                          // return number of indices not equal to -1
                          return std::count_if(values.begin(), values.end(), [](double x) { return x != -1; }); },
                         {"hq1_idx_" + recoJetsCollection, "hq2_idx_" + recoJetsCollection,
                          "v1q1_idx_" + recoJetsCollection, "v1q2_idx_" + recoJetsCollection, "v2q1_idx_" + recoJetsCollection, "v2q2_idx_" + recoJetsCollection});

  return df_out;
}

RNode reconstructRecoEventWithTruthInfo_Boosted(RNode df, std::string recoJetsCollection)
{
  std::cout << " -> Run reconstructRecoEventWithTruthInfo_Boosted(recoJetsCollection = " << recoJetsCollection << ")" << std::endl;
  //
  // Matching priority: Higgs, Leading V, Trailing V
  // Matching requirements: 
  //    1. minimum dR < 0.8 
  //    2. The two boson daughters are contained in the matched large-R jet 
  //

  std::vector recoJetsVars = { recoJetsCollection + "_pt", // order matters
                               recoJetsCollection + "_eta",
                               recoJetsCollection + "_phi",
                               recoJetsCollection + "_mass" };

  auto df_out = df.Define(recoJetsCollection + "_matchedIndices", truthMatch_RecoGenAK8_HV1V2,
                          {"GenJetAK8_H_pt", "GenJetAK8_H_eta", "GenJetAK8_H_phi", "GenJetAK8_H_mass",
                           "GenJetAK8_V1_pt", "GenJetAK8_V1_eta", "GenJetAK8_V1_phi", "GenJetAK8_V1_mass",
                           "GenJetAK8_V2_pt", "GenJetAK8_V2_eta", "GenJetAK8_V2_phi", "GenJetAK8_V2_mass",
                           recoJetsVars.at(0), recoJetsVars.at(1), recoJetsVars.at(2), recoJetsVars.at(3)});

  df_out = df_out.Define("h_idx_" + recoJetsCollection, recoJetsCollection + "_matchedIndices[0]")
               .Define("v1_idx_" + recoJetsCollection, recoJetsCollection + "_matchedIndices[1]")
               .Define("v2_idx_" + recoJetsCollection, recoJetsCollection + "_matchedIndices[2]");

  // Old matching: match truth partons directly to AK8 jets (keeping for reference, can delete eventually)
  // >>
  df_out = df_out.Define(recoJetsCollection + "_oldMatchedIndices_HV1V2", truthMatch_GenJetsAK8_HV1V2,
                         {"truthH_pt", "truthH_eta", "truthH_phi", "truthH_mass",
                          "truthH_daught1_pt", "truthH_daught1_eta", "truthH_daught1_phi", "truthH_daught1_mass", "truthH_daught1_pdgId",
                          "truthH_daught2_pt", "truthH_daught2_eta", "truthH_daught2_phi", "truthH_daught2_mass", "truthH_daught2_pdgId",
                          "truthV1_pt", "truthV1_eta", "truthV1_phi", "truthV1_mass",
                          "truthV1_daught1_pt", "truthV1_daught1_eta", "truthV1_daught1_phi", "truthV1_daught1_mass", "truthV1_daught1_pdgId",
                          "truthV1_daught2_pt", "truthV1_daught2_eta", "truthV1_daught2_phi", "truthV1_daught2_mass", "truthV1_daught2_pdgId",
                          "truthV2_pt", "truthV2_eta", "truthV2_phi", "truthV2_mass",
                          "truthV2_daught1_pt", "truthV2_daught1_eta", "truthV2_daught1_phi", "truthV2_daught1_mass", "truthV2_daught1_pdgId",
                          "truthV2_daught2_pt", "truthV2_daught2_eta", "truthV2_daught2_phi", "truthV2_daught2_mass", "truthV2_daught2_pdgId",
                          recoJetsVars.at(0), recoJetsVars.at(1), recoJetsVars.at(2), recoJetsVars.at(3)});

  df_out = df_out.Define("h_idx_oldMatching_" + recoJetsCollection, recoJetsCollection + "_oldMatchedIndices_HV1V2[0]")
               .Define("v1_idx_oldMatching_" + recoJetsCollection, recoJetsCollection + "_oldMatchedIndices_HV1V2[1]")
               .Define("v2_idx_oldMatching_" + recoJetsCollection, recoJetsCollection + "_oldMatchedIndices_HV1V2[2]");
  // <<
  
  // Count number of bosons matched to a large-R jet
  df_out = df_out.Define("nTruthMatchedBosons_" + recoJetsCollection, [](const int &h_idx, const int &v1_idx, const int &v2_idx) -> int {
                          std::vector<int> values = {h_idx, v1_idx, v2_idx};
                          // return number of indices not equal to -1
                          return std::count_if(values.begin(), values.end(), [](double x) { return x != -1; }); },
                         {"h_idx_" + recoJetsCollection, "v1_idx_" + recoJetsCollection,
                          "v2_idx_" + recoJetsCollection});

  df_out = df_out.Define("dR_H_truth_matched" + recoJetsCollection, dR_particle_ak4AtIdx, {"truthH_pt", "truthH_eta", "truthH_phi", "truthH_mass", recoJetsVars.at(0), recoJetsVars.at(1), recoJetsVars.at(2), recoJetsVars.at(3), "h_idx_" + recoJetsCollection})
               .Define("dR_V1_truth_matched" + recoJetsCollection, dR_particle_ak4AtIdx, {"truthV1_pt", "truthV1_eta", "truthV1_phi", "truthV1_mass", recoJetsVars.at(0), recoJetsVars.at(1), recoJetsVars.at(2), recoJetsVars.at(3), "v1_idx_" + recoJetsCollection})
               .Define("dR_V2_truth_matched" + recoJetsCollection, dR_particle_ak4AtIdx, {"truthV2_pt", "truthV2_eta", "truthV2_phi", "truthV2_mass", recoJetsVars.at(0), recoJetsVars.at(1), recoJetsVars.at(2), recoJetsVars.at(3), "v2_idx_" + recoJetsCollection})

               .Define("dR_H_matchedGenJet_matched" + recoJetsCollection, dR_particle_ak4AtIdx, {"GenJetAK8_H_pt", "GenJetAK8_H_eta", "GenJetAK8_H_phi", "GenJetAK8_H_mass", recoJetsVars.at(0), recoJetsVars.at(1), recoJetsVars.at(2), recoJetsVars.at(3), "h_idx_" + recoJetsCollection})
               .Define("dR_V1_matchedGenJet_matched" + recoJetsCollection, dR_particle_ak4AtIdx, {"GenJetAK8_V1_pt", "GenJetAK8_V1_eta", "GenJetAK8_V1_phi", "GenJetAK8_V1_mass", recoJetsVars.at(0), recoJetsVars.at(1), recoJetsVars.at(2), recoJetsVars.at(3), "v1_idx_" + recoJetsCollection})
               .Define("dR_V2_matchedGenJet_matched" + recoJetsCollection, dR_particle_ak4AtIdx, {"GenJetAK8_V2_pt", "GenJetAK8_V2_eta", "GenJetAK8_V2_phi", "GenJetAK8_V2_mass", recoJetsVars.at(0), recoJetsVars.at(1), recoJetsVars.at(2), recoJetsVars.at(3), "v2_idx_" + recoJetsCollection});

  df_out = df_out.Define("H_matched_" + recoJetsCollection + "_m", [](const ROOT::RVecF &jm, int matched_idx)
                         {
                            float mass = -99;
                            if(matched_idx != -1){
                              mass = jm.at(matched_idx);
                            }
                            return mass; },
                         {recoJetsVars.at(3), "h_idx_" + recoJetsCollection})
               .Define("V1_matched_" + recoJetsCollection + "_m", [](const ROOT::RVecF &jm, int matched_idx)
                       {
                            float mass = -99;
                            if(matched_idx != -1){
                              mass = jm.at(matched_idx);
                            }
                            return mass; },
                       {recoJetsVars.at(3), "v1_idx_" + recoJetsCollection})
               .Define("V2_matched_" + recoJetsCollection + "_m", [](const ROOT::RVecF &jm, int matched_idx)
                       {
                            float mass = -99;
                            if(matched_idx != -1){
                              mass = jm.at(matched_idx);
                            }
                            return mass; },
                       {recoJetsVars.at(3), "v2_idx_" + recoJetsCollection});
  return df_out;
}



//
// AK4 matching with final state quarks 
//
bool isQuark(int pdgId)
{
    std::vector<int> possible_pdgids = {1, 2, 3, 4, 5, 6, -1, -2, -3, -4, -5, -6};

    return std::find(possible_pdgids.begin(), possible_pdgids.end(), pdgId) != possible_pdgids.end();
}

bool isQuarkOrGluon(int pdgId)
{
  std::vector<int> possible_pdgids = {1, 2, 3, 4, 5, 6, -1, -2, -3, -4, -5, -6, 21};

  return std::find(possible_pdgids.begin(), possible_pdgids.end(), pdgId) != possible_pdgids.end();
}

std::vector<int> truthMatchGenJetToRecoWithoutAK8overlap(
    std::vector<int> indices_to_match,                                                                           // indices to be used to retrieve jets in inj collection
    const ROOT::RVecF &injpt, const ROOT::RVecF &injeta, const ROOT::RVecF &injphi, const ROOT::RVecF &injm,     // collection to retrieve jets from using indices_to_match
    const ROOT::RVecF &outjpt, const ROOT::RVecF &outjeta, const ROOT::RVecF &outjphi, const ROOT::RVecF &outjm, // collection to match jets in inj collection to
    std::vector<int> indices_to_notoverlap,                                                                      // indices to be used to retrieve jets in ak8 collection for overlap removal
    const ROOT::RVecF &jak8pt, const ROOT::RVecF &jak8eta, const ROOT::RVecF &jak8phi, const ROOT::RVecF &jak8m  // ak8 collection to retrieve jets from using indices_to_notoverlap
)
{
  double matchingDR = 0.3;
  // Build AK8 jets not to overlap with
  std::vector<ROOT::Math::PtEtaPhiMVector> ak8OverlapRemoval_p4;
  for (unsigned int i = 0; i < indices_to_notoverlap.size(); i++){
    int ak8_index = indices_to_notoverlap.at(i);
    if (ak8_index > -1)
    {
      ROOT::Math::PtEtaPhiMVector vec(jak8pt.at(ak8_index), jak8eta.at(ak8_index), jak8phi.at(ak8_index), jak8m.at(ak8_index));
      ak8OverlapRemoval_p4.push_back(vec);
    }
    else{
      ROOT::Math::PtEtaPhiMVector vec(0,0,0,0);
      ak8OverlapRemoval_p4.push_back(ROOT::Math::PtEtaPhiMVector(0, 0, 0, 0));
    }
  }

  int n_injets = indices_to_match.size();
  int n_outjets = outjpt.size();
  int n_ORjets = ak8OverlapRemoval_p4.size();

  assert(n_ORjets == 3); // code assumes the full list is passed. If a boson was not matched, it handles the case of idx = -1. 

  std::vector<int> used_indices{};
  std::vector<int> matched_indices{};

  for (int i = 0; i < n_injets; i++)
  {
    double smallestDR = matchingDR;
    int in_idx = indices_to_match.at(i);

    // if the jet to be matched was already not matched, leave it unmatched (-1)
    if (in_idx < 0)
    {
      matched_indices.push_back(in_idx);
      continue;
    }

    // otherwise, build its 4-vector
    ROOT::Math::PtEtaPhiMVector in_p4(injpt.at(in_idx), injeta.at(in_idx), injphi.at(in_idx), injm.at(in_idx));
    int candidate_out_idx = -1; // default if no match is found

    for (int j = 0; j < n_outjets; j++)
    {
      ROOT::Math::PtEtaPhiMVector out_p4(outjpt.at(j), outjeta.at(j), outjphi.at(j), outjm.at(j));
      double dR = ROOT::Math::VectorUtil::DeltaR(in_p4, out_p4);

      if (dR < smallestDR)
      {
        
        // check reco ak4 jet index has not already been matched 
        if (std::find(used_indices.begin(), used_indices.end(), j) == used_indices.end())
        {

          // check the reco ak4 jet is not contained in any overlap removal jet 
          bool doesNotOverlap = true;
          for( int k=0; k < n_ORjets; k++){
            // if i is a boson daughter (i<=5), only apply overlap with ak8 jets matched to bosons that are not the daughter's mother
            if ((indices_to_notoverlap.at(k) > -1) && ((i / 2 != k) || (i > 5)) )
              {
                ROOT::Math::PtEtaPhiMVector ak8p4 = ak8OverlapRemoval_p4.at(k);
                float dR_overlap = ROOT::Math::VectorUtil::DeltaR(out_p4, ak8p4);
                if (dR_overlap <= 0.8)
                {
                  //std::cout << "Overlapping between quark " << i << " and boson " << k << " with dR = " << dR_overlap << std::endl;
                  doesNotOverlap = false;
                  break;
                }
              }
          }
          if (doesNotOverlap)
          {
            // if arrived here, jet has not already been used and does not overlap with specified ak8jets
            candidate_out_idx = j;
            smallestDR = dR;
          }
        }
      }
    }

    matched_indices.push_back(candidate_out_idx);
    used_indices.push_back(candidate_out_idx); // Do not reuse the index
  }

  return matched_indices;
}

int truthMatchedAk4Idx(float ppt, float peta, float pphi, float pm, const ROOT::RVecF &jpt, const ROOT::RVecF &jeta, const ROOT::RVecF &jphi, const ROOT::RVecF &jm, const std::vector<int> &used_idxs, float smallest_dR)
{
  int candidate_idx = -1;
  int size = jpt.size();
  // if (!isQuarkOrGluon(pdgId))
  // {
  //   return candidate_idx; //-1
  // }
  ROOT::Math::PtEtaPhiMVector parent_p4(ppt, peta, pphi, pm);
  for (int i = 0; i < size; ++i)
  {
    ROOT::Math::PtEtaPhiMVector jet_p4(jpt.at(i), jeta.at(i), jphi.at(i), jm.at(i));
    double dR = ROOT::Math::VectorUtil::DeltaR(parent_p4, jet_p4);
    if (dR < smallest_dR)
    {
      // check reco ak4 jet index has not already been matched to another parton
      if (std::find(used_idxs.begin(), used_idxs.end(), i) == used_idxs.end()){
        candidate_idx = i;
        smallest_dR = dR;
      }
    }
  }
  return candidate_idx;
}

int truthMatchedJetArrIdx(float ppt, float peta, float pphi, float pm,  
    const ROOT::RVecF &jpt, const ROOT::RVecF &jeta, const ROOT::RVecF &jphi, const ROOT::RVecF &jm, 
    const std::vector<int> &used_idxs, float smallest_dR)
{
  int candidate_idx = -1;
  //double smallest_dR = 0.8;
  int njets = jpt.size();

  ROOT::Math::PtEtaPhiMVector parent_p4(ppt, peta, pphi, pm);
  for (int i = 0; i < njets; ++i)
  {
    // check reco jet index has not already been matched to another parton
    if (std::find(used_idxs.begin(), used_idxs.end(), i) != used_idxs.end())
    {
      continue;
    }
    ROOT::Math::PtEtaPhiMVector jet_p4(jpt.at(i), jeta.at(i), jphi.at(i), jm.at(i));
    double dR_jet_parent = ROOT::Math::VectorUtil::DeltaR(parent_p4, jet_p4);

    // Find jet closest to parent
    if (dR_jet_parent < smallest_dR)
    {
      candidate_idx = i;
      smallest_dR = dR_jet_parent;
    }
  }
  return candidate_idx;
}

int truthMatchedAk8Idx_daughtContainement(float ppt, float peta, float pphi, float pm,  
    float d1_pt, float d1_eta, float d1_phi, float d1_m,
    float d2_pt, float d2_eta, float d2_phi, float d2_m,
    const ROOT::RVecF &jpt, const ROOT::RVecF &jeta, const ROOT::RVecF &jphi, const ROOT::RVecF &jm, 
    const std::vector<int> &used_idxs)
  {
  int candidate_idx = -1;
  double smallest_dR = 0.8;
  int nak8 = jpt.size();

  // if (!(isQuarkOrGluon(d1_pdgId) && isQuarkOrGluon(d2_pdgId)))
  // {
  //   //std::cout << "It's not a quark, returning. " << std::endl;
  //   return candidate_idx; //-1
  // }

  ROOT::Math::PtEtaPhiMVector parent_p4(ppt, peta, pphi, pm);
  ROOT::Math::PtEtaPhiMVector daught1_p4(d1_pt, d1_eta, d1_phi, d1_m);
  ROOT::Math::PtEtaPhiMVector daught2_p4(d2_pt, d2_eta, d2_phi, d2_m);
  for (int i = 0; i < nak8; ++i)
  {
    // check reco jet index has not already been matched to another parton
    if (std::find(used_idxs.begin(), used_idxs.end(), i) != used_idxs.end()){
      continue;
    }
    ROOT::Math::PtEtaPhiMVector jet_p4(jpt.at(i), jeta.at(i), jphi.at(i), jm.at(i));
    double dR_jet_parent = ROOT::Math::VectorUtil::DeltaR(parent_p4, jet_p4);
    //std::cout << "dR_jet_parent = " << dR_jet_parent << std::endl;

    // Find jet closest to parent
    if (dR_jet_parent < smallest_dR)
    {
      //std::cout << "Passed dR matching" << std::endl;
      // Only match the parent if the daughters are contained in the large-R jet
      double dR_jet_daught1 = ROOT::Math::VectorUtil::DeltaR(daught1_p4, jet_p4);
      double dR_jet_daught2 = ROOT::Math::VectorUtil::DeltaR(daught2_p4, jet_p4);

      //std::cout << "dR_jet_daught1 = " << dR_jet_daught1 << ", dR_jet_daught2 = " << dR_jet_daught2 << std::endl;
      if ( (dR_jet_daught1 < 0.8) && (dR_jet_daught2 < 0.8)) {
        //std::cout << "Passed containmenet" << std::endl;
        candidate_idx = i;
        smallest_dR = dR_jet_parent;    
      }
    }
  }
  //std::cout << "Returning candidate_idx = " << candidate_idx << std::endl;
  return candidate_idx;
}

float dR_particle_ak4AtIdx(float &ppt, float &peta, float &pphi, float &pm, const ROOT::RVecF &jpt, const ROOT::RVecF &jeta, const ROOT::RVecF &jphi, const ROOT::RVecF &jm, int j_idx)
{
  float dR = -99;
  ROOT::Math::PtEtaPhiMVector truth_p4(ppt, peta, pphi, pm);
  if (j_idx != -1)
  {
    ROOT::Math::PtEtaPhiMVector jet_p4(jpt.at(j_idx), jeta.at(j_idx), jphi.at(j_idx), jm.at(j_idx));
    dR = ROOT::Math::VectorUtil::DeltaR(truth_p4, jet_p4);
  }
  return dR;
}

float delta_R(float pt1, float eta1, float phi1, float m1, float pt2, float eta2, float phi2, float m2)
{
  ROOT::Math::PtEtaPhiMVector first_p4(pt1, eta1, phi1, m1);
  ROOT::Math::PtEtaPhiMVector second_p4(pt2, eta2, phi2, m2);
  return ROOT::Math::VectorUtil::DeltaR(first_p4, second_p4);
}

float dijetInvMass(float pt1, float eta1, float phi1, float m1, float pt2, float eta2, float phi2, float m2)
{
  ROOT::Math::PtEtaPhiMVector v1_p4(pt1, eta1, phi1, m1);
  ROOT::Math::PtEtaPhiMVector v2_p4(pt2, eta2, phi2, m2);
  ROOT::Math::PtEtaPhiMVector v_sum = v1_p4 + v2_p4;
  return v_sum.M();
}

float dijetInvMass_fromIdx(const ROOT::RVecF &vec_pt, const ROOT::RVecF &vec_eta, const ROOT::RVecF &vec_phi, const ROOT::RVecF &vec_m, int idx1, int idx2)
{
  float m = -99;
  if (idx1 != -1 && idx2 != -1)
  {
    ROOT::Math::PtEtaPhiMVector j1(vec_pt.at(idx1), vec_eta.at(idx1), vec_phi.at(idx1), vec_m.at(idx1));
    ROOT::Math::PtEtaPhiMVector j2(vec_pt.at(idx2), vec_eta.at(idx2), vec_phi.at(idx2), vec_m.at(idx2));
    ROOT::Math::PtEtaPhiMVector j12 = j1 + j2;
    m = j12.M();
  }
  return m;
}