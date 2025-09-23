#ifndef TRUTHSELECTIONS_H
#define TRUTHSELECTIONS_H

#include "ROOT/RDataFrame.hxx"
#include "Math/VectorUtil.h"
#include "utils.h"

RNode reconstructRecoEventWithTruthInfo_Resolved(RNode df_, std::string recoJetsCollection = "goodAK4Jets", std::string recoJetsCollectionAK8= "goodAK8Jets");
RNode reconstructRecoEventWithTruthInfo_Boosted(RNode df, std::string recoJetsCollection = "goodAK8Jets");
RNode findTruthMatchedVBSJets(RNode df, std::string recoAK4JetsCollection, std::string recoAK8JetsCollection);
RNode addTruthVariables(RNode df);
int truthMatchedAk4Idx(float ppt, float peta, float pphi, float pm,
                       const ROOT::RVecF &jpt, const ROOT::RVecF &jeta, const ROOT::RVecF &jphi, const ROOT::RVecF &jm,
                       const std::vector<int> &used_idxs, float smallest_dR = 0.3);
int truthMatchedAk8Idx_daughtContainement(float ppt, float peta, float pphi, float pm,
                       float d1_pt, float d1_eta, float d1_phi, float d1_m,
                       float d2_pt, float d2_eta, float d2_phi, float d2_m,
                       const ROOT::RVecF &jpt, const ROOT::RVecF &jeta, const ROOT::RVecF &jphi, const ROOT::RVecF &jm,
                       const std::vector<int> &used_idxs);
int truthMatchedJetArrIdx(float ppt, float peta, float pphi, float pm,
                          const ROOT::RVecF &jpt, const ROOT::RVecF &jeta, const ROOT::RVecF &jphi, const ROOT::RVecF &jm,
                          const std::vector<int> &used_idxs, float smallest_dR);

std::vector<int> truthMatchedAK4VBS(float tvbs1_pt, float tvbs1_eta, float tvbs1_phi, float tvbs1_m,
                                    float tvbs2_pt, float tvbs2_eta, float tvbs2_phi, float tvbs2_m,
                                    const ROOT::RVecF &jpt, const ROOT::RVecF &jeta, const ROOT::RVecF &jphi, const ROOT::RVecF &jm,
                                    const ROOT::RVecF &jak8pt, const ROOT::RVecF &jak8eta, const ROOT::RVecF &jak8phi, const ROOT::RVecF &jak8m,
                                    const int &h_ak8idx, const int &v1_ak8idx, const int &v2_ak8idx,
                                    const int &h_q1_ak4idx, const int &h_q2_ak4idx, const int &v1_q1_ak4idx,
                                    const int &v1_q2_ak4idx, const int &v2_q1_ak4idx, const int &v2_q2_ak4idx);
bool isQuark(int pdgId);
bool isQuarkOrGluon(int pdgId);
float ak4_reco_dr_idx(const ROOT::RVecF &ppt, const ROOT::RVecF &peta, const ROOT::RVecF &pphi, const ROOT::RVecF &pm, const ROOT::RVecF &jpt, const ROOT::RVecF &jeta, const ROOT::RVecF &jphi, const ROOT::RVecF &jm, int j_idx);
float dR_particle_ak4AtIdx(float &ppt, float &peta, float &pphi, float &pm, const ROOT::RVecF &jpt, const ROOT::RVecF &jeta, const ROOT::RVecF &jphi, const ROOT::RVecF &jm, int j_idx);
float delta_R(float pt1, float eta1, float phi1, float m1, float pt2, float eta2, float phi2, float m2);
float dijetInvMass_fromIdx(const ROOT::RVecF &vec_pt, const ROOT::RVecF &vec_eta, const ROOT::RVecF &vec_phi, const ROOT::RVecF &vec_m, int idx1, int idx2);
float dijetInvMass(float pt1, float eta1, float phi1, float m1, float pt2, float eta2, float phi2, float m2);

std::vector<int> truthMatch_GenJets_HV1V2daughters(
    float truthH_daught1_pt, float truthH_daught1_eta, float truthH_daught1_phi, float truthH_daught1_mass,
    float truthH_daught2_pt, float truthH_daught2_eta, float truthH_daught2_phi, float truthH_daught2_mass,
    float truthV1_daught1_pt, float truthV1_daught1_eta, float truthV1_daught1_phi, float truthV1_daught1_mass,
    float truthV1_daught2_pt, float truthV1_daught2_eta, float truthV1_daught2_phi, float truthV1_daught2_mass,
    float truthV2_daught1_pt, float truthV2_daught1_eta, float truthV2_daught1_phi, float truthV2_daught1_mass,
    float truthV2_daught2_pt, float truthV2_daught2_eta, float truthV2_daught2_phi, float truthV2_daught2_mass,
    const ROOT::RVecF GenJet_pt, const ROOT::RVecF GenJet_eta, const ROOT::RVecF GenJet_phi, const ROOT::RVecF GenJet_mass);

std::vector<int> truthMatch_GenJets(
    float truthH_daught1_pt, float truthH_daught1_eta, float truthH_daught1_phi, float truthH_daught1_mass, int truthH_daught1_pdgId,
    float truthH_daught2_pt, float truthH_daught2_eta, float truthH_daught2_phi, float truthH_daught2_mass, int truthH_daught2_pdgId,
    float truthV1_daught1_pt, float truthV1_daught1_eta, float truthV1_daught1_phi, float truthV1_daught1_mass, int truthV1_daught1_pdgId,
    float truthV1_daught2_pt, float truthV1_daught2_eta, float truthV1_daught2_phi, float truthV1_daught2_mass, int truthV1_daught2_pdgId,
    float truthV2_daught1_pt, float truthV2_daught1_eta, float truthV2_daught1_phi, float truthV2_daught1_mass, int truthV2_daught1_pdgId,
    float truthV2_daught2_pt, float truthV2_daught2_eta, float truthV2_daught2_phi, float truthV2_daught2_mass, int truthV2_daught2_pdgId,
    float truthVBSq1_ptOrd_pt, float truthVBSq1_ptOrd_eta, float truthVBSq1_ptOrd_phi, float truthVBSq1_ptOrd_mass, int truthVBSq1_ptOrd_pdgId,
    float truthVBSq2_ptOrd_pt, float truthVBSq2_ptOrd_eta, float truthVBSq2_ptOrd_phi, float truthVBSq2_ptOrd_mass, int truthVBSq2_ptOrd_pdgId,
    const ROOT::RVecF GenJet_pt, const ROOT::RVecF GenJet_eta, const ROOT::RVecF GenJet_phi, const ROOT::RVecF GenJet_mass);

RNode reconstructTruthEvent(RNode df);

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
    const ROOT::RVecF GenJetAK8_pt, const ROOT::RVecF GenJetAK8_eta, const ROOT::RVecF GenJetAK8_phi, const ROOT::RVecF GenJetAK8_mass);

std::vector<int> truthMatch_RecoGenAK8_HV1V2(
    float genH_pt, float genH_eta, float genH_phi, float genH_mass,
    float genV1_pt, float genV1_eta, float genV1_phi, float genV1_mass,
    float genV2_pt, float genV2_eta, float genV2_phi, float genV2_mass,
    const ROOT::RVecF AK8Jet_pt, const ROOT::RVecF AK8Jet_eta, const ROOT::RVecF AK8Jet_phi, const ROOT::RVecF AK8Jet_mass);
    
std::vector<int> truthMatch_GenJets_VBS(float truthVBSq1_ptOrd_pt, float truthVBSq1_ptOrd_eta, float truthVBSq1_ptOrd_phi, float truthVBSq1_ptOrd_mass,
                                        float truthVBSq2_ptOrd_pt, float truthVBSq2_ptOrd_eta, float truthVBSq2_ptOrd_phi, float truthVBSq2_ptOrd_mass,
                                        const ROOT::RVecF GenJet_pt, const ROOT::RVecF GenJet_eta, const ROOT::RVecF GenJet_phi, const ROOT::RVecF GenJet_mass);

std::vector<int> truthMatchGenJetToRecoWithoutAK8overlap(
    std::vector<int> indices_to_match,                                                                           // indices to be used to retrieve jets in inj collection
    const ROOT::RVecF &injpt, const ROOT::RVecF &injeta, const ROOT::RVecF &injphi, const ROOT::RVecF &injm,     // collection to retrieve jets from using indices
    const ROOT::RVecF &outjpt, const ROOT::RVecF &outjeta, const ROOT::RVecF &outjphi, const ROOT::RVecF &outjm, // collection to match jets to
    std::vector<int> indices_to_notoverlap,                                                                      // indices to be used to retrieve jets in kak8 collection
    const ROOT::RVecF &jak8pt, const ROOT::RVecF &jak8eta, const ROOT::RVecF &jak8phi, const ROOT::RVecF &jak8m // collection to avoid overlap with
    );


#endif // TRUTHSELECTIONS_H