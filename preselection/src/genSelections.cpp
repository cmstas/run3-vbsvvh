#include "genSelections.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <vector>

// Match a gen quark to an AK4 jet.
//
// Overlap removal is deliberately limited to two things:
//   - an index already assigned to an earlier object is unavailable;
//   - a jet lying inside an AK8 jet that a boosted boson already claimed is unavailable.
// There is NO AK4-AK4 dR cleaning: anti-kt with R=0.4 already decided that two entries in
// the collection are two jets, so vetoing a candidate for being close to an already-matched
// jet would only throw away correct matches.
//
// NOTE: the AK4 collection passed in ("jet") is already cleaned against every good fat jet
// at dR > 0.8 in AK4JetsSelection, so the containment check below cannot currently fire. It
// is kept because it is exactly the check needed if this is ever pointed at the uncleaned
// "jetNoFJClean" collection, where the resolved/boosted split would be decided against the
// fat jets actually matched to a boson rather than against every good fat jet in the event.
int find_matching_jet(int target_idx, float target_eta, float target_phi, ROOT::RVec<int> already_matched_jet_indices, ROOT::RVec<int> already_matched_fatjet_indices, ROOT::RVec<float> jet_eta, ROOT::RVec<float> jet_phi, ROOT::RVec<float> fatjet_eta, ROOT::RVec<float> fatjet_phi) {
    int max_jets = 10;
    int max_fatjets = 3;
    const float dR_cut = 0.4f;
    const float fatjet_containment_cut = 0.8f;

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

        // Already assigned to an earlier object
        bool is_excluded = false;
        for (int excl_idx : already_matched_jet_indices) {
            if (idx == excl_idx) {
                is_excluded = true;
                break;
            }
        }
        if (is_excluded) continue;

        // Inside an AK8 jet that a boosted boson already claimed
        bool inside_matched_fatjet = false;
        for (int matched_fj_idx : already_matched_fatjet_indices) {
            if (matched_fj_idx >= 0 && matched_fj_idx < max_fatjets) {
                float dR_fatjet = ROOT::VecOps::DeltaR(jet_eta[idx], fatjet_eta[matched_fj_idx], jet_phi[idx], fatjet_phi[matched_fj_idx]);
                if (dR_fatjet < fatjet_containment_cut) {
                    inside_matched_fatjet = true;
                    break;
                }
            }
        }
        if (inside_matched_fatjet) continue;

        return idx;
    }
    return -1;
}

// Match a gen boson to an AK8 jet.
//
// As above, no AK8-AK8 dR cleaning: an AK8 index already claimed by an earlier boson is
// removed from the pool, but a distinct fat jet is a distinct fat jet. The "contains an AK4
// jet already assigned to a VBS quark" veto below is inert against the FJ-cleaned AK4
// collection, and is kept for the same reason as the containment check above.
int find_matching_fatjet(int target_idx, float target_eta, float target_phi, ROOT::RVec<int> already_matched_jet_indices, ROOT::RVec<int> already_matched_fatjet_indices, ROOT::RVec<float> jet_eta, ROOT::RVec<float> jet_phi, ROOT::RVec<float> fatjet_eta, ROOT::RVec<float> fatjet_phi) {
    int max_jets = 10;
    int max_fatjets = 3;
    const float dR_cut = 0.8f;
    const float jet_containment_cut = 0.8f;

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

        // Already claimed by an earlier boson
        bool is_excluded = false;
        for (int excl_idx : already_matched_fatjet_indices) {
            if (idx == excl_idx) {
                is_excluded = true;
                break;
            }
        }
        if (is_excluded) continue;

        // Contains an AK4 jet already assigned to a VBS quark
        bool contains_matched_jet = false;
        for (int matched_j_idx : already_matched_jet_indices) {
            if (matched_j_idx >= 0 && matched_j_idx < max_jets) {
                float dR_jet = ROOT::VecOps::DeltaR(fatjet_eta[idx], jet_eta[matched_j_idx], fatjet_phi[idx], jet_phi[matched_j_idx]);
                if (dR_jet < jet_containment_cut) {
                    contains_matched_jet = true;
                    break;
                }
            }
        }
        if (contains_matched_jet) continue;

        return idx;
    }
    return -1;
}

// |pdgId| in [1,6] (a quark) or 21 (a gluon). Anything else -- a charged lepton, a
// neutrino, or the 0 used when an index is unfilled -- fails.
static std::string isQuarkOrGluon(const std::string &col) {
    return "((std::abs(" + col + ") >= 1 && std::abs(" + col + ") <= 6) || std::abs(" + col + ") == 21)";
}

RNode GenSelections(RNode df_) {
    // AK4 collection used for the resolved matching. "jet" is cleaned against every good
    // fat jet at dR > 0.8 upstream, which for a 1-fat-jet channel costs a resolved match
    // only when that fat jet is matched to no truth boson. Switching these to
    // "jetNoFJClean_eta"/"jetNoFJClean_phi" moves the decision into find_matching_jet,
    // where it is made against the fat jets actually matched to a boson.
    const std::string jet_eta_col = "jet_eta";
    const std::string jet_phi_col = "jet_phi";

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

    // Require the boson's daughters to be quarks or gluons before matching anything to it.
    // The skim fills gen_<v>q<n>_idx with whatever the V decayed to (leptons and neutrinos
    // included), and gen_b1/b2_idx follow any H decay, not just H -> bb. The matching is purely
    // geometric, so without this gate a hadronic tau or a neutrino could be matched to a jet and
    // labelled as a quark. Gating the TARGET index is enough: find_matching_jet and
    // find_matching_fatjet return -1 for a negative target, so both the resolved labels and
    // the boosted flag fall through to "unmatched".
    df = df.Define("_gen_b1_pdgId",   "gen_b1_idx   >= 0 ? GenPart_pdgId[gen_b1_idx]   : 0")
        .Define("_gen_b2_pdgId",   "gen_b2_idx   >= 0 ? GenPart_pdgId[gen_b2_idx]   : 0")
        .Define("_gen_v1q1_pdgId", "gen_v1q1_idx >= 0 ? GenPart_pdgId[gen_v1q1_idx] : 0")
        .Define("_gen_v1q2_pdgId", "gen_v1q2_idx >= 0 ? GenPart_pdgId[gen_v1q2_idx] : 0")
        .Define("_gen_v2q1_pdgId", "gen_v2q1_idx >= 0 ? GenPart_pdgId[gen_v2q1_idx] : 0")
        .Define("_gen_v2q2_pdgId", "gen_v2q2_idx >= 0 ? GenPart_pdgId[gen_v2q2_idx] : 0")
        .Define("_h_isHadronic",  isQuarkOrGluon("_gen_b1_pdgId")   + " && " + isQuarkOrGluon("_gen_b2_pdgId"))
        .Define("_v1_isHadronic", isQuarkOrGluon("_gen_v1q1_pdgId") + " && " + isQuarkOrGluon("_gen_v1q2_pdgId"))
        .Define("_v2_isHadronic", isQuarkOrGluon("_gen_v2q1_pdgId") + " && " + isQuarkOrGluon("_gen_v2q2_pdgId"))
        .Define("_gen_h_idx_had",    "_h_isHadronic  ? gen_h_idx    : -1")
        .Define("_gen_b1_idx_had",   "_h_isHadronic  ? gen_b1_idx   : -1")
        .Define("_gen_b2_idx_had",   "_h_isHadronic  ? gen_b2_idx   : -1")
        .Define("_gen_v1_idx_had",   "_v1_isHadronic ? gen_v1_idx   : -1")
        .Define("_gen_v1q1_idx_had", "_v1_isHadronic ? gen_v1q1_idx : -1")
        .Define("_gen_v1q2_idx_had", "_v1_isHadronic ? gen_v1q2_idx : -1")
        .Define("_gen_v2_idx_had",   "_v2_isHadronic ? gen_v2_idx   : -1")
        .Define("_gen_v2q1_idx_had", "_v2_isHadronic ? gen_v2q1_idx : -1")
        .Define("_gen_v2q2_idx_had", "_v2_isHadronic ? gen_v2q2_idx : -1");

    df = df.Define("_empty_exclusions", "ROOT::RVec<int>{}");

    // 1. VBS quarks -> AK4. Nothing has been claimed yet.
    df = df.Define("_vbs1_idx_temp", find_matching_jet, {"gen_vbs1_idx", "gen_vbs1_eta", "gen_vbs1_phi", "_empty_exclusions", "_empty_exclusions", jet_eta_col, jet_phi_col, "fatjet_eta", "fatjet_phi"})
        .Define("_vbs1_exclusions", "ROOT::RVec<int>{_vbs1_idx_temp}")
        .Define("_vbs2_idx_temp", find_matching_jet, {"gen_vbs2_idx", "gen_vbs2_eta", "gen_vbs2_phi", "_vbs1_exclusions", "_empty_exclusions", jet_eta_col, jet_phi_col, "fatjet_eta", "fatjet_phi"})
        .Define("truth_vbs1_idx", "_vbs1_idx_temp >= 0 && _vbs1_idx_temp < 10 ? _vbs1_idx_temp : -1")
        .Define("truth_vbs2_idx", "_vbs2_idx_temp >= 0 && _vbs2_idx_temp < 10 ? _vbs2_idx_temp : -1")
        .Define("_matched_vbs_jets", "ROOT::RVec<int>{truth_vbs1_idx, truth_vbs2_idx}");

    // 2. Boosted H -> AK8. Veto fat jets containing a VBS jet.
    df = df.Define("_hbb_dR", "ROOT::VecOps::DeltaR(gen_b1_eta, gen_b2_eta, gen_b1_phi, gen_b2_phi)")
        .Define("_hbb_fatjet_idx_temp", find_matching_fatjet, {"_gen_h_idx_had", "gen_h_eta", "gen_h_phi", "_matched_vbs_jets", "_empty_exclusions", jet_eta_col, jet_phi_col, "fatjet_eta", "fatjet_phi"})
        .Define("_hbb_fatjet_candidate_b1_dR", "_hbb_fatjet_idx_temp >= 0 ? ROOT::VecOps::DeltaR(fatjet_eta[_hbb_fatjet_idx_temp], gen_b1_eta, fatjet_phi[_hbb_fatjet_idx_temp], gen_b1_phi) : 999.0")
        .Define("_hbb_fatjet_candidate_b2_dR", "_hbb_fatjet_idx_temp >= 0 ? ROOT::VecOps::DeltaR(fatjet_eta[_hbb_fatjet_idx_temp], gen_b2_eta, fatjet_phi[_hbb_fatjet_idx_temp], gen_b2_phi) : 999.0")
        .Define("_hbb_isBoosted", "_h_isHadronic && _hbb_fatjet_idx_temp >= 0 && _hbb_dR < 0.8 && _hbb_fatjet_candidate_b1_dR < 0.8 && _hbb_fatjet_candidate_b2_dR < 0.8")
        .Define("truth_h_idx", "_hbb_isBoosted ? _hbb_fatjet_idx_temp : -1");

    // 3. Boosted V1 -> AK8. H's fat jet is out of the pool.
    df = df.Define("_v1qq_dR", "gen_v1_idx >= 0 ? ROOT::VecOps::DeltaR(gen_v1q1_eta, gen_v1q2_eta, gen_v1q1_phi, gen_v1q2_phi) : 999.0")
        .Define("_matched_fatjets_for_v1", "ROOT::RVec<int>{truth_h_idx}")
        .Define("_v1qq_fatjet_idx_temp", find_matching_fatjet, {"_gen_v1_idx_had", "gen_v1_eta", "gen_v1_phi", "_matched_vbs_jets", "_matched_fatjets_for_v1", jet_eta_col, jet_phi_col, "fatjet_eta", "fatjet_phi"})
        .Define("_v1qq_fatjet_candidate_q1_dR", "_v1qq_fatjet_idx_temp >= 0 ? ROOT::VecOps::DeltaR(fatjet_eta[_v1qq_fatjet_idx_temp], gen_v1q1_eta, fatjet_phi[_v1qq_fatjet_idx_temp], gen_v1q1_phi) : 999.0")
        .Define("_v1qq_fatjet_candidate_q2_dR", "_v1qq_fatjet_idx_temp >= 0 ? ROOT::VecOps::DeltaR(fatjet_eta[_v1qq_fatjet_idx_temp], gen_v1q2_eta, fatjet_phi[_v1qq_fatjet_idx_temp], gen_v1q2_phi) : 999.0")
        .Define("_v1qq_isBoosted", "_v1_isHadronic && _v1qq_fatjet_idx_temp >= 0 && _v1qq_dR < 0.8 && _v1qq_fatjet_candidate_q1_dR < 0.8 && _v1qq_fatjet_candidate_q2_dR < 0.8")
        .Define("truth_v1_idx", "_v1qq_isBoosted ? _v1qq_fatjet_idx_temp : -1");

    // 4. Boosted V2 -> AK8. H's and V1's fat jets are out of the pool.
    df = df.Define("_v2qq_dR", "gen_v2_idx >= 0 ? ROOT::VecOps::DeltaR(gen_v2q1_eta, gen_v2q2_eta, gen_v2q1_phi, gen_v2q2_phi) : 999.0")
        .Define("_matched_fatjets_for_v2", "ROOT::RVec<int>{truth_h_idx, truth_v1_idx}")
        .Define("_v2qq_fatjet_idx_temp", find_matching_fatjet, {"_gen_v2_idx_had", "gen_v2_eta", "gen_v2_phi", "_matched_vbs_jets", "_matched_fatjets_for_v2", jet_eta_col, jet_phi_col, "fatjet_eta", "fatjet_phi"})
        .Define("_v2qq_fatjet_candidate_q1_dR", "_v2qq_fatjet_idx_temp >= 0 ? ROOT::VecOps::DeltaR(fatjet_eta[_v2qq_fatjet_idx_temp], gen_v2q1_eta, fatjet_phi[_v2qq_fatjet_idx_temp], gen_v2q1_phi) : 999.0")
        .Define("_v2qq_fatjet_candidate_q2_dR", "_v2qq_fatjet_idx_temp >= 0 ? ROOT::VecOps::DeltaR(fatjet_eta[_v2qq_fatjet_idx_temp], gen_v2q2_eta, fatjet_phi[_v2qq_fatjet_idx_temp], gen_v2q2_phi) : 999.0")
        .Define("_v2qq_isBoosted", "_v2_isHadronic && _v2qq_fatjet_idx_temp >= 0 && _v2qq_dR < 0.8 && _v2qq_fatjet_candidate_q1_dR < 0.8 && _v2qq_fatjet_candidate_q2_dR < 0.8")
        .Define("truth_v2_idx", "_v2qq_isBoosted ? _v2qq_fatjet_idx_temp : -1");

    // Every fat jet that a boson actually claimed; AK4 jets inside these are unavailable
    // to the resolved matching below.
    df = df.Define("_matched_fatjets", "ROOT::RVec<int>{truth_h_idx, truth_v1_idx, truth_v2_idx}");

    // 5. Resolved H -> AK4.
    df = df.Define("_b1_idx_temp", find_matching_jet, {"_gen_b1_idx_had", "gen_b1_eta", "gen_b1_phi", "_matched_vbs_jets", "_matched_fatjets", jet_eta_col, jet_phi_col, "fatjet_eta", "fatjet_phi"})
        .Define("_excluded_jets_for_b2", "ROOT::RVec<int>{truth_vbs1_idx, truth_vbs2_idx, _b1_idx_temp}")
        .Define("_b2_idx_temp", find_matching_jet, {"_gen_b2_idx_had", "gen_b2_eta", "gen_b2_phi", "_excluded_jets_for_b2", "_matched_fatjets", jet_eta_col, jet_phi_col, "fatjet_eta", "fatjet_phi"})
        .Define("truth_b1_idx", "_b1_idx_temp >= 0 && _b1_idx_temp < 10 ? _b1_idx_temp : -1")
        .Define("truth_b2_idx", "_b2_idx_temp >= 0 && _b2_idx_temp < 10 ? _b2_idx_temp : -1");

    // 6. Resolved V1 -> AK4.
    df = df.Define("_excluded_jets_for_v1q1", "ROOT::RVec<int>{truth_vbs1_idx, truth_vbs2_idx, truth_b1_idx, truth_b2_idx}")
        .Define("_v1q1_idx_temp", find_matching_jet, {"_gen_v1q1_idx_had", "gen_v1q1_eta", "gen_v1q1_phi", "_excluded_jets_for_v1q1", "_matched_fatjets", jet_eta_col, jet_phi_col, "fatjet_eta", "fatjet_phi"})
        .Define("_excluded_jets_for_v1q2", "ROOT::RVec<int>{truth_vbs1_idx, truth_vbs2_idx, truth_b1_idx, truth_b2_idx, _v1q1_idx_temp}")
        .Define("_v1q2_idx_temp", find_matching_jet, {"_gen_v1q2_idx_had", "gen_v1q2_eta", "gen_v1q2_phi", "_excluded_jets_for_v1q2", "_matched_fatjets", jet_eta_col, jet_phi_col, "fatjet_eta", "fatjet_phi"})
        .Define("truth_v1q1_idx", "_v1q1_idx_temp >= 0 && _v1q1_idx_temp < 10 ? _v1q1_idx_temp : -1")
        .Define("truth_v1q2_idx", "_v1q2_idx_temp >= 0 && _v1q2_idx_temp < 10 ? _v1q2_idx_temp : -1");

    // 7. Resolved V2 -> AK4.
    df = df.Define("_excluded_jets_for_v2q1", "ROOT::RVec<int>{truth_vbs1_idx, truth_vbs2_idx, truth_b1_idx, truth_b2_idx, truth_v1q1_idx, truth_v1q2_idx}")
        .Define("_v2q1_idx_temp", find_matching_jet, {"_gen_v2q1_idx_had", "gen_v2q1_eta", "gen_v2q1_phi", "_excluded_jets_for_v2q1", "_matched_fatjets", jet_eta_col, jet_phi_col, "fatjet_eta", "fatjet_phi"})
        .Define("_excluded_jets_for_v2q2", "ROOT::RVec<int>{truth_vbs1_idx, truth_vbs2_idx, truth_b1_idx, truth_b2_idx, truth_v1q1_idx, truth_v1q2_idx, _v2q1_idx_temp}")
        .Define("_v2q2_idx_temp", find_matching_jet, {"_gen_v2q2_idx_had", "gen_v2q2_eta", "gen_v2q2_phi", "_excluded_jets_for_v2q2", "_matched_fatjets", jet_eta_col, jet_phi_col, "fatjet_eta", "fatjet_phi"})
        .Define("truth_v2q1_idx", "_v2q1_idx_temp >= 0 && _v2q1_idx_temp < 10 ? _v2q1_idx_temp : -1")
        .Define("truth_v2q2_idx", "_v2q2_idx_temp >= 0 && _v2q2_idx_temp < 10 ? _v2q2_idx_temp : -1");

    return df;
}

namespace {

using ROOT::VecOps::RVec;

// ---------------------------------------------------------------------------
// Per-fat-jet gen-truth (MC only), computed for every FatJet so any channel can
// read it and the per-variation FatJet_isGood_<sfx> masks select it like any other
// FatJet_ attribute. It depends on eta/phi only, so it is the same in every JES/JER
// variation.
//
// Strategy: JET-FIRST. Ask what truth object, if any, is contained in this jet --
// never start from the sample's expected topology and check containment afterwards.
// The second direction forces an a priori choice exactly where it is least safe: in
// tW->4Q, TTW-WtoQQ, TTWW and TTWZ a top-W and a prompt W are both hadronic in the
// same event, and any fixed priority between them is applied blind.
//
// IDENTITY BEFORE CONTAINMENT. The candidate is the hadronic W or Z, or the H->bb,
// whose two quarks best match the jet (smallest max dR); what it is then decides
// which containment test applies. The answer is one code:
//
//   FatJet_genMatch   0 nothing   1 prompt W   2 W from top   3 top   4 prompt Z
//                     5 H->bb
//
// 0 covers both "no candidate boson in the record" and "one exists but is not in
// the cone"; FatJet_genDrVqMax (999 vs finite-but->=kMatchDR) separates them.
//
//   FatJet_genDrVqMax    max dR(jet, q) over the candidate boson's two quarks
//   FatJet_genDrTopB     dR(jet, b) for the b from the candidate W's top (999 if the
//                        W is not from a top, or the top has no b daughter)
//   FatJet_genVqFlavour  0 light, 1 at least one c and no b, 2 at least one b,
//                        -1 no candidate boson in the event
//
// Only H->bb is a Higgs candidate; other Higgs decays never match.
//
// The W quarks are contained within kMatchDR. The top's b is tested against the full
// AK8 radius, kTopBDR, so the top / W-from-top split leaves no b partly in the cone:
//   top          both W quarks within kMatchDR AND dR(jet, b) <  kTopBDR
//   W from top   both W quarks within kMatchDR AND dR(jet, b) >= kTopBDR
// Only a b counts as the top's partner: a CKM-suppressed t->sW / t->dW is labelled by
// its W alone.
//
// NanoAOD v15 types, verified against the v30 skims:
//   GenPart_pdgId            Int_t    -> int
//   GenPart_genPartIdxMother Short_t  -> short
//   GenPart_statusFlags      UShort_t -> unsigned short
// isLastCopy is bit 13 of statusFlags.
namespace gentruth {

// Containment radius for the candidate boson's two quarks.
constexpr float kMatchDR = 0.6f;

// The top's b counts as inside the jet anywhere within the AK8 radius.
constexpr float kTopBDR = 0.8f;

constexpr float kNone = 999.f;

enum : int { kMatchNone = 0, kMatchPromptW = 1, kMatchTopW = 2, kMatchTop = 3, kMatchPromptZ = 4, kMatchHbb = 5 };

enum : int { kFlavNoV = -1, kFlavLight = 0, kFlavCharm = 1, kFlavBottom = 2 };

struct FatJetTruth {
    RVec<int>   match;
    RVec<float> drVqMax;
    RVec<float> drTopB;
    RVec<int>   vqFlavour;
};

enum class Boson { W, Z, H };

// A last-copy W or Z decaying to quarks, or H->bb: its two quarks, and for a W from
// a top, the top's b.
struct HadronicV {
    int q1, q2;
    Boson kind;
    bool fromTop;
    int topB;
};

} // namespace gentruth

gentruth::FatJetTruth fatjet_gen_truth(
    const RVec<float>& fj_eta, const RVec<float>& fj_phi,
    const RVec<int>& pdg,
    const RVec<short>& mom,
    const RVec<unsigned short>& flags,
    const RVec<float>& geta,
    const RVec<float>& gphi)
{
    using namespace gentruth;
    const int n = static_cast<int>(pdg.size());
    auto isLastCopy = [&](int i) { return (flags[i] >> 13) & 1; };

    // Walk the mother chain from i up to a |pdg|==6 top; return its index or -1.
    auto topAncestor = [&](int i) -> int {
        int cur = static_cast<int>(mom[i]);
        for (int guard = 0; cur >= 0 && cur < n && guard < 50; ++guard) {
            if (std::abs(pdg[cur]) == 6) return cur;
            cur = static_cast<int>(mom[cur]);
        }
        return -1;
    };

    // The quark daughters of the boson at index v. A last-copy W/Z/H decays either to
    // quarks or to leptons, never both -- downstream radiation attaches to the quarks,
    // not back to the boson -- so no lepton guard is needed. (Measured on 240k
    // last-copy W in tW: the quark multiplicity is exactly 0 or 2, never 1 or 3.) A
    // count != 2 is a record we do not understand, so it is rejected rather than
    // silently truncated to the first two.
    auto quarkDaughters = [&](int v, int& q1, int& q2) -> bool {
        q1 = -1; q2 = -1; int nq = 0;
        for (int j = 0; j < n; ++j) {
            if (static_cast<int>(mom[j]) != v) continue;
            const int ad = std::abs(pdg[j]);
            if (ad < 1 || ad > 5) continue;
            ++nq;
            if      (q1 < 0) q1 = j;
            else if (q2 < 0) q2 = j;
        }
        return nq == 2;
    };

    // Candidates are found once per event; every jet is then compared to all of them.
    std::vector<HadronicV> cands;
    for (int i = 0; i < n; ++i) {
        const int ad = std::abs(pdg[i]);
        if ((ad != 24 && ad != 23 && ad != 25) || !isLastCopy(i)) continue;
        const Boson kind = ad == 24 ? Boson::W : ad == 23 ? Boson::Z : Boson::H;
        HadronicV v{-1, -1, kind, false, -1};
        if (!quarkDaughters(i, v.q1, v.q2)) continue;   // no quark decay, or rejected
        if (kind == Boson::H && (std::abs(pdg[v.q1]) != 5 || std::abs(pdg[v.q2]) != 5)) continue;
        if (kind == Boson::W) {
            const int top = topAncestor(i);
            if (top >= 0) {
                v.fromTop = true;
                for (int j = 0; j < n; ++j) {
                    if (static_cast<int>(mom[j]) == top && std::abs(pdg[j]) == 5) { v.topB = j; break; }
                }
            }
        }
        cands.push_back(v);
    }

    const size_t nfj = fj_eta.size();
    FatJetTruth out{RVec<int>(nfj, kMatchNone), RVec<float>(nfj, kNone),
                    RVec<float>(nfj, kNone), RVec<int>(nfj, kFlavNoV)};
    if (cands.empty()) return out;

    for (size_t k = 0; k < nfj; ++k) {
        auto dRjet = [&](int i) { return ROOT::VecOps::DeltaR(fj_eta[k], geta[i], fj_phi[k], gphi[i]); };

        const HadronicV* best = nullptr;
        float best_dr = 1e9f;
        for (const auto& v : cands) {
            const float dr = std::max(dRjet(v.q1), dRjet(v.q2));
            if (dr < best_dr) { best_dr = dr; best = &v; }
        }

        const int a1 = std::abs(pdg[best->q1]), a2 = std::abs(pdg[best->q2]);
        const int flav = (a1 == 5 || a2 == 5) ? kFlavBottom
                       : (a1 == 4 || a2 == 4) ? kFlavCharm
                       :                        kFlavLight;
        const float dr_top_b = (best->fromTop && best->topB >= 0) ? dRjet(best->topB) : kNone;

        int match;
        if      (best_dr >= kMatchDR)        match = kMatchNone;
        else if (best->kind == Boson::H)     match = kMatchHbb;
        else if (best->kind == Boson::Z)     match = kMatchPromptZ;
        else if (!best->fromTop)             match = kMatchPromptW;
        else if (dr_top_b < kTopBDR)         match = kMatchTop;
        else                                 match = kMatchTopW;

        out.match[k]     = match;
        out.drVqMax[k]   = best_dr;
        out.drTopB[k]    = dr_top_b;
        out.vqFlavour[k] = flav;
    }
    return out;
}

} // anonymous namespace

RNode FatJetGenTruth(RNode df_) {
    return df_.Define("_fatjet_gen_truth", fatjet_gen_truth,
                      {"FatJet_eta", "FatJet_phi", "GenPart_pdgId", "GenPart_genPartIdxMother",
                       "GenPart_statusFlags", "GenPart_eta", "GenPart_phi"})
              .Define("FatJet_genMatch",     [](const gentruth::FatJetTruth& t) { return t.match; },     {"_fatjet_gen_truth"})
              .Define("FatJet_genDrVqMax",   [](const gentruth::FatJetTruth& t) { return t.drVqMax; },   {"_fatjet_gen_truth"})
              .Define("FatJet_genDrTopB",    [](const gentruth::FatJetTruth& t) { return t.drTopB; },    {"_fatjet_gen_truth"})
              .Define("FatJet_genVqFlavour", [](const gentruth::FatJetTruth& t) { return t.vqFlavour; }, {"_fatjet_gen_truth"});
}
