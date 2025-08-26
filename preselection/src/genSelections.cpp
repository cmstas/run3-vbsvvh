#include "genSelections.h"

int find_matching_jet(int target_idx, float target_eta, float target_phi, ROOT::RVec<int> already_matched_jet_indices, ROOT::RVec<int> already_matched_fatjet_indices, ROOT::RVec<float> jet_eta, ROOT::RVec<float> jet_phi, ROOT::RVec<float> fatjet_eta, ROOT::RVec<float> fatjet_phi) {
    int max_jets = 10;
    int max_fatjets = 3;
    const float dR_cut = 0.4f;
    const float fatjet_overlap_cut = 0.8f;
    const float jet_overlap_cut = 0.4f;
                       
    if (target_idx < 0) {
        return -1;
    }
    
    // Calculate dR values for all jets
    ROOT::RVec<float> dR_values;
    for (size_t i = 0; i < jet_eta.size(); ++i) {
        dR_values.push_back(ROOT::VecOps::DeltaR(target_eta, jet_eta[i], target_phi, jet_phi[i]));
    }
    
    auto sorted_indices = ROOT::VecOps::Argsort(dR_values);
    for (int idx : sorted_indices) {
        if (dR_values[idx] >= dR_cut) break; // No more candidates within dR cut
        if (idx >= max_jets) continue; // Skip indices beyond padding limit
        
        // Check if this jet is in the excluded list
        bool is_excluded = false;
        for (int excl_idx : already_matched_jet_indices) {
            if (idx == excl_idx) {
                is_excluded = true;
                break;
            }
        }
        if (is_excluded) continue;
        
        bool overlaps_jet = false;
        for (int matched_j_idx : already_matched_jet_indices) {
            if (matched_j_idx >= 0 && matched_j_idx < max_jets) {
                float dR_jet = ROOT::VecOps::DeltaR(jet_eta[idx], jet_eta[matched_j_idx], jet_phi[idx], jet_phi[matched_j_idx]);
                if (dR_jet < jet_overlap_cut) {
                    overlaps_jet = true;
                    break;
                }
            }
        }

        // Check overlap with matched fatjets only
        bool overlaps_fatjet = false;
        for (int matched_fj_idx : already_matched_fatjet_indices) {
            if (matched_fj_idx >= 0 && matched_fj_idx < max_fatjets) {
                float dR_fatjet = ROOT::VecOps::DeltaR(jet_eta[idx], fatjet_eta[matched_fj_idx], jet_phi[idx], fatjet_phi[matched_fj_idx]);
                if (dR_fatjet < fatjet_overlap_cut) {
                    overlaps_fatjet = true;
                    break;
                }
            }
        }
        if (!overlaps_jet && !overlaps_fatjet) {
            return idx;
        }
    }
    return -1;
}

int find_matching_fatjet(int target_idx, float target_eta, float target_phi, ROOT::RVec<int> already_matched_jet_indices, ROOT::RVec<int> already_matched_fatjet_indices, ROOT::RVec<float> jet_eta, ROOT::RVec<float> jet_phi, ROOT::RVec<float> fatjet_eta, ROOT::RVec<float> fatjet_phi) {
    int max_jets = 10;
    int max_fatjets = 3;
    const float dR_cut = 0.8f;
    const float fatjet_overlap_cut = 0.8f;
    const float jet_overlap_cut = 0.8f;
    
    if (target_idx < 0) {
        return -1;
    }
                       
    // Calculate dR values for all fatjets
    ROOT::RVec<float> dR_values;
    for (size_t i = 0; i < fatjet_eta.size(); ++i) {
        dR_values.push_back(ROOT::VecOps::DeltaR(target_eta, fatjet_eta[i], target_phi, fatjet_phi[i]));
    }
    
    auto sorted_indices = ROOT::VecOps::Argsort(dR_values);
    for (int idx : sorted_indices) {
        if (dR_values[idx] >= dR_cut) break; // No more candidates within dR cut
        if (idx >= max_fatjets) continue; // Skip indices beyond padding limit
        
        // Check if this fatjet is in the excluded list
        bool is_excluded = false;
        for (int excl_idx : already_matched_fatjet_indices) {
            if (idx == excl_idx) {
                is_excluded = true;
                break;
            }
        }
        if (is_excluded) continue;

        // Check overlap with matched jets only
        bool overlaps_jet = false;
        for (int matched_j_idx : already_matched_jet_indices) {
            if (matched_j_idx >= 0 && matched_j_idx < max_jets) {
                float dR_jet = ROOT::VecOps::DeltaR(fatjet_eta[idx], jet_eta[matched_j_idx], fatjet_phi[idx], jet_phi[matched_j_idx]);
                if (dR_jet < jet_overlap_cut) {
                    overlaps_jet = true;
                    break;
                }
            }
        }
        
        bool overlaps_fatjet = false;
        for (int matched_fj_idx : already_matched_fatjet_indices) {
            if (matched_fj_idx >= 0 && matched_fj_idx < max_fatjets) {
                float dR_fatjet = ROOT::VecOps::DeltaR(fatjet_eta[idx], fatjet_eta[matched_fj_idx], fatjet_phi[idx], fatjet_phi[matched_fj_idx]);
                if (dR_fatjet < fatjet_overlap_cut) {
                    overlaps_fatjet = true;
                    break;
                }
            }
        }

        if (!overlaps_jet && !overlaps_fatjet) {
            return idx;
        }
    }
    return -1;
}
                       
RNode GenSelections(RNode df_) {
    auto df = df_.Define("gen_vbs1_eta", "gen_vbs1_idx >= 0 ? GenPart_eta[gen_vbs1_idx] : -999.0f")
        .Define("gen_vbs1_phi", "gen_vbs1_idx >= 0 ? GenPart_phi[gen_vbs1_idx] : -999.0f")
        .Define("gen_vbs2_eta", "gen_vbs2_idx >= 0 ? GenPart_eta[gen_vbs2_idx] : -999.0f")
        .Define("gen_vbs2_phi", "gen_vbs2_idx >= 0 ? GenPart_phi[gen_vbs2_idx] : -999.0f")
        .Define("gen_h_eta", "gen_h_idx >= 0 ? GenPart_eta[gen_h_idx] : -999.0f")
        .Define("gen_h_phi", "gen_h_idx >= 0 ? GenPart_phi[gen_h_idx] : -999.0f")
        .Define("gen_b1_eta", "gen_b1_idx >= 0 ? GenPart_eta[gen_b1_idx] : -999.0f")
        .Define("gen_b1_phi", "gen_b1_idx >= 0 ? GenPart_phi[gen_b1_idx] : -999.0f")
        .Define("gen_b2_eta", "gen_b2_idx >= 0 ? GenPart_eta[gen_b2_idx] : -999.0f")
        .Define("gen_b2_phi", "gen_b2_idx >= 0 ? GenPart_phi[gen_b2_idx] : -999.0f")
        .Define("gen_v1_eta", "gen_v1_idx >= 0 ? GenPart_eta[gen_v1_idx] : -999.0f")
        .Define("gen_v1_phi", "gen_v1_idx >= 0 ? GenPart_phi[gen_v1_idx] : -999.0f")
        .Define("gen_v2_eta", "gen_v2_idx >= 0 ? GenPart_eta[gen_v2_idx] : -999.0f")
        .Define("gen_v2_phi", "gen_v2_idx >= 0 ? GenPart_phi[gen_v2_idx] : -999.0f")
        .Define("gen_v1q1_eta", "gen_v1q1_idx >= 0 ? GenPart_eta[gen_v1q1_idx] : -999.0f")
        .Define("gen_v1q1_phi", "gen_v1q1_idx >= 0 ? GenPart_phi[gen_v1q1_idx] : -999.0f")
        .Define("gen_v1q2_eta", "gen_v1q2_idx >= 0 ? GenPart_eta[gen_v1q2_idx] : -999.0f")
        .Define("gen_v1q2_phi", "gen_v1q2_idx >= 0 ? GenPart_phi[gen_v1q2_idx] : -999.0f")
        .Define("gen_v2q1_eta", "gen_v2q1_idx >= 0 ? GenPart_eta[gen_v2q1_idx] : -999.0f")
        .Define("gen_v2q1_phi", "gen_v2q1_idx >= 0 ? GenPart_phi[gen_v2q1_idx] : -999.0f")
        .Define("gen_v2q2_eta", "gen_v2q2_idx >= 0 ? GenPart_eta[gen_v2q2_idx] : -999.0f")
        .Define("gen_v2q2_phi", "gen_v2q2_idx >= 0 ? GenPart_phi[gen_v2q2_idx] : -999.0f");
    
    df = df.Define("empty_jet_exclusions", "ROOT::RVec<int>{}");

    df = df.Define("vbs1_idx_temp", find_matching_jet, {"gen_vbs1_idx", "gen_vbs1_eta", "gen_vbs1_phi", "empty_jet_exclusions", "empty_jet_exclusions", "Jet_eta", "Jet_phi", "FatJet_eta", "FatJet_phi"})
        .Define("vbs1_exclusions", "ROOT::RVec<int>{vbs1_idx_temp}")
        .Define("vbs2_idx_temp", find_matching_jet, {"gen_vbs2_idx", "gen_vbs2_eta", "gen_vbs2_phi", "vbs1_exclusions", "empty_jet_exclusions", "Jet_eta", "Jet_phi", "FatJet_eta", "FatJet_phi"})
        .Define("truth_vbs1_idx", "vbs1_idx_temp >= 0 && vbs1_idx_temp < 10 ? vbs1_idx_temp : -1")
        .Define("truth_vbs2_idx", "vbs2_idx_temp >= 0 && vbs2_idx_temp < 10 ? vbs2_idx_temp : -1");

    df = df.Define("excluded_jets_for_bh", "ROOT::RVec<int>{truth_vbs1_idx, truth_vbs2_idx}")
        .Define("hbb_dR", "ROOT::VecOps::DeltaR(gen_b1_eta, gen_b2_eta, gen_b1_phi, gen_b2_phi)")
        .Define("hbb_fatjet_idx_temp", find_matching_fatjet, {"gen_h_idx", "gen_h_eta", "gen_h_phi", "excluded_jets_for_bh", "empty_jet_exclusions", "Jet_eta", "Jet_phi", "FatJet_eta", "FatJet_phi"})
        .Define("hbb_fatjet_candidate_b1_dR", "hbb_fatjet_idx_temp >= 0 && hbb_fatjet_idx_temp < FatJet_eta.size() ? ROOT::VecOps::DeltaR(FatJet_eta[hbb_fatjet_idx_temp], gen_b1_eta, FatJet_phi[hbb_fatjet_idx_temp], gen_b1_phi) : 999.0")
        .Define("hbb_fatjet_candidate_b2_dR", "hbb_fatjet_idx_temp >= 0 && hbb_fatjet_idx_temp < FatJet_eta.size() ? ROOT::VecOps::DeltaR(FatJet_eta[hbb_fatjet_idx_temp], gen_b2_eta, FatJet_phi[hbb_fatjet_idx_temp], gen_b2_phi) : 999.0")
        .Define("hbb_isBoosted", "hbb_fatjet_idx_temp != -1 && hbb_fatjet_idx_temp < 3 && hbb_dR < 0.8 && hbb_fatjet_candidate_b1_dR < 0.8 && hbb_fatjet_candidate_b2_dR < 0.8")
        .Define("truth_h_idx", "hbb_isBoosted ? hbb_fatjet_idx_temp : -1");
    
    df = df.Define("excluded_jet_mask_for_hbb", "ROOT::RVec<int>{truth_vbs1_idx, truth_vbs2_idx}")
        .Define("b1_idx_temp", find_matching_jet, {"gen_b1_idx", "gen_b1_eta", "gen_b1_phi", "excluded_jet_mask_for_hbb", "empty_jet_exclusions", "Jet_eta", "Jet_phi", "FatJet_eta", "FatJet_phi"})
        .Define("excluded_jet_mask_for_hbb2", "ROOT::RVec<int>{truth_vbs1_idx, truth_vbs2_idx, b1_idx_temp}")
        .Define("b2_idx_temp", find_matching_jet, {"gen_b2_idx", "gen_b2_eta", "gen_b2_phi", "excluded_jet_mask_for_hbb2", "empty_jet_exclusions", "Jet_eta", "Jet_phi", "FatJet_eta", "FatJet_phi"})
        .Define("truth_b1_idx", "b1_idx_temp >= 0 && b1_idx_temp < 10 ? b1_idx_temp : -1")
        .Define("truth_b2_idx", "b2_idx_temp >= 0 && b2_idx_temp < 10 ? b2_idx_temp : -1");

    df = df.Define("v1qq_dR", "gen_v1_idx != -1 ? ROOT::VecOps::DeltaR(gen_v1q1_eta, gen_v1q1_phi, gen_v1q2_eta, gen_v1q2_phi) : 999.0")
        .Define("matched_jet_mask_for_v1", "ROOT::RVec<int>{truth_vbs1_idx, truth_vbs2_idx, truth_b1_idx, truth_b2_idx}")
        .Define("excluded_fatjet_mask_for_v1", "ROOT::RVec<int>{truth_h_idx}")
        .Define("v1qq_fatjet_idx_temp", find_matching_fatjet, {"gen_v1_idx", "gen_v1_eta", "gen_v1_phi", "matched_jet_mask_for_v1", "excluded_fatjet_mask_for_v1", "Jet_eta", "Jet_phi", "FatJet_eta", "FatJet_phi"})
        .Define("v1qq_fatjet_candidate_q1_dR", "v1qq_fatjet_idx_temp >= 0 && v1qq_fatjet_idx_temp < FatJet_eta.size() ? ROOT::VecOps::DeltaR(FatJet_eta[v1qq_fatjet_idx_temp], gen_v1q1_eta, FatJet_phi[v1qq_fatjet_idx_temp], gen_v1q1_phi) : 999.0")
        .Define("v1qq_fatjet_candidate_q2_dR", "v1qq_fatjet_idx_temp >= 0 && v1qq_fatjet_idx_temp < FatJet_eta.size() ? ROOT::VecOps::DeltaR(FatJet_eta[v1qq_fatjet_idx_temp], gen_v1q2_eta, FatJet_phi[v1qq_fatjet_idx_temp], gen_v1q2_phi) : 999.0")
        .Define("v1qq_isBoosted", "gen_v1_idx != -1 && v1qq_fatjet_idx_temp != -1 && v1qq_fatjet_idx_temp < 3 && v1qq_dR < 0.8 && v1qq_fatjet_candidate_q1_dR < 0.8 && v1qq_fatjet_candidate_q2_dR < 0.8")
        .Define("truth_v1_idx", "v1qq_isBoosted ? v1qq_fatjet_idx_temp : -1");

    df = df.Define("v2qq_dR", "gen_v2_idx != -1 ? ROOT::VecOps::DeltaR(gen_v2q1_eta, gen_v2q1_phi, gen_v2q2_eta, gen_v2q2_phi) : 999.0")
        .Define("excluded_fatjet_mask_for_v2", "ROOT::RVec<int>{truth_h_idx, truth_v1_idx}")
        .Define("matched_jet_mask_for_v2", "ROOT::RVec<int>{truth_vbs1_idx, truth_vbs2_idx, truth_b1_idx, truth_b2_idx}")
        .Define("v2qq_fatjet_idx_temp", find_matching_fatjet, {"gen_v2_idx", "gen_v2_eta", "gen_v2_phi", "matched_jet_mask_for_v2", "excluded_fatjet_mask_for_v2", "Jet_eta", "Jet_phi", "FatJet_eta", "FatJet_phi"})
        .Define("v2qq_fatjet_candidate_q1_dR", "v2qq_fatjet_idx_temp >= 0 && v2qq_fatjet_idx_temp < FatJet_eta.size() ? ROOT::VecOps::DeltaR(FatJet_eta[v2qq_fatjet_idx_temp], gen_v2q1_eta, FatJet_phi[v2qq_fatjet_idx_temp], gen_v2q1_phi) : 999.0")
        .Define("v2qq_fatjet_candidate_q2_dR", "v2qq_fatjet_idx_temp >= 0 && v2qq_fatjet_idx_temp < FatJet_eta.size() ? ROOT::VecOps::DeltaR(FatJet_eta[v2qq_fatjet_idx_temp], gen_v2q2_eta, FatJet_phi[v2qq_fatjet_idx_temp], gen_v2q2_phi) : 999.0")
        .Define("v2qq_isBoosted", "gen_v2_idx != -1 && v2qq_fatjet_idx_temp != -1 && v2qq_fatjet_idx_temp < 3 && v2qq_dR < 0.8 && v2qq_fatjet_candidate_q1_dR < 0.8 && v2qq_fatjet_candidate_q2_dR < 0.8")
        .Define("truth_v2_idx", "v2qq_isBoosted ? v2qq_fatjet_idx_temp : -1");
    
    df = df.Define("excluded_jet_mask_for_v1q1", "ROOT::RVec<int>{truth_vbs1_idx, truth_vbs2_idx, truth_b1_idx, truth_b2_idx}")
        .Define("matched_fatjet_mask_for_v1q", "ROOT::RVec<int>{truth_h_idx, truth_v2_idx}")
        .Define("v1q1_idx_temp", find_matching_jet, {"gen_v1q1_idx", "gen_v1q1_eta", "gen_v1q1_phi", "excluded_jet_mask_for_v1q1", "matched_fatjet_mask_for_v1q", "Jet_eta", "Jet_phi", "FatJet_eta", "FatJet_phi"})
        .Define("excluded_jet_mask_for_v1q2", "ROOT::RVec<int>{truth_vbs1_idx, truth_vbs2_idx, truth_b1_idx, truth_b2_idx, v1q1_idx_temp}")
        .Define("v1q2_idx_temp", find_matching_jet, {"gen_v1q2_idx", "gen_v1q2_eta", "gen_v1q2_phi", "excluded_jet_mask_for_v1q2", "matched_fatjet_mask_for_v1q", "Jet_eta", "Jet_phi", "FatJet_eta", "FatJet_phi"})
        .Define("truth_v1q1_idx", "v1q1_idx_temp >= 0 && v1q1_idx_temp < 10 ? v1q1_idx_temp : -1")
        .Define("truth_v1q2_idx", "v1q2_idx_temp >= 0 && v1q2_idx_temp < 10 ? v1q2_idx_temp : -1");

    df = df.Define("excluded_jet_mask_for_v2q1", "ROOT::RVec<int>{truth_vbs1_idx, truth_vbs2_idx, truth_b1_idx, truth_b2_idx, truth_v1q1_idx, truth_v1q2_idx}")
        .Define("matched_fatjet_mask_for_v2q", "ROOT::RVec<int>{truth_h_idx, truth_v1_idx}")
        .Define("v2q1_idx_temp", find_matching_jet, {"gen_v2q1_idx", "gen_v2q1_eta", "gen_v2q1_phi", "excluded_jet_mask_for_v2q1", "matched_fatjet_mask_for_v2q", "Jet_eta", "Jet_phi", "FatJet_eta", "FatJet_phi"})
        .Define("excluded_jet_mask_for_v2q2", "ROOT::RVec<int>{truth_vbs1_idx, truth_vbs2_idx, truth_b1_idx, truth_b2_idx, truth_v1q1_idx, truth_v1q2_idx, v2q1_idx_temp}")
        .Define("v2q2_idx_temp", find_matching_jet, {"gen_v2q2_idx", "gen_v2q2_eta", "gen_v2q2_phi", "excluded_jet_mask_for_v2q2", "matched_fatjet_mask_for_v2q", "Jet_eta", "Jet_phi", "FatJet_eta", "FatJet_phi"})
        .Define("truth_v2q1_idx", "v2q1_idx_temp >= 0 && v2q1_idx_temp < 10 ? v2q1_idx_temp : -1")
        .Define("truth_v2q2_idx", "v2q2_idx_temp >= 0 && v2q2_idx_temp < 10 ? v2q2_idx_temp : -1");

    return df;
}