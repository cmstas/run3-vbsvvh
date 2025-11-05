#include "spanet_inference_run2.h"
#include <vector>
#include <algorithm>
#include <utility>
#include <set>
#include <limits>

namespace SPANet {

// FIXME this is ugly! Data members of the derived class (e.g. INPUT_NAMES) are constructed before the 
// base class constructor body, but after the base class initializer list arguments are evaluated. 
// This means INPUT_NAMES has to already exist when you passed it into the base class constructor.
const std::vector<const char*> SPANetInferenceRun2::INPUT_NAMES = {
    "AK4Jets_data", "AK4Jets_mask", "AK8Jets_data", "AK8Jets_mask", 
    "HadronicActivity_data", "HadronicActivity_mask"
};

SPANetInferenceRun2::SPANetInferenceRun2(const std::string &model_path, size_t batch_size)
    : SPANetInferenceBase(model_path, batch_size, INPUT_NAMES, MAX_AK4_JETS, AK4_FEATURES, MAX_AK8_JETS, AK8_FEATURES, MET_FEATURES, MAX_LEPTONS, LEPTON_FEATURES) {}

std::vector<std::shared_ptr<EventDataBase>> SPANetInferenceRun2::extractEventsFromDataFrame(RNode df) {
    std::cout << "--> SPANetInferenceRun2::extractEventsFromDataFrame()" << std::endl;

    auto ak4_pt_vec = df.Take<RVec<float>>("goodAK4Jets_pt").GetValue();
    auto ak4_eta_vec = df.Take<RVec<float>>("goodAK4Jets_eta").GetValue();
    auto ak4_phi_vec = df.Take<RVec<float>>("goodAK4Jets_phi").GetValue();
    auto ak4_mass_vec = df.Take<RVec<float>>("goodAK4Jets_mass").GetValue();
    auto ak4_isTightBTag_vec = df.Take<RVec<int>>("goodAK4Jets_isTightBTag").GetValue();
    auto ak4_isMediumBTag_vec = df.Take<RVec<int>>("goodAK4Jets_isMediumBTag").GetValue();
    auto ak4_isLooseBTag_vec = df.Take<RVec<int>>("goodAK4Jets_isLooseBTag").GetValue();

    auto ak8_pt_vec = df.Take<RVec<float>>("goodAK8Jets_pt").GetValue();
    auto ak8_eta_vec = df.Take<RVec<float>>("goodAK8Jets_eta").GetValue();
    auto ak8_phi_vec = df.Take<RVec<float>>("goodAK8Jets_phi").GetValue();
    auto ak8_mass_vec = df.Take<RVec<float>>("goodAK8Jets_msoftdrop").GetValue();

    auto ak8_nConstituents_vec = df.Take<RVec<unsigned char>>("goodAK8Jets_nConstituents").GetValue();
    auto ak8_HbbScore_vec = df.Take<RVec<float>>("goodAK8Jets_HbbScore").GetValue();
    auto ak8_WqqScore_vec = df.Take<RVec<float>>("goodAK8Jets_WqqScore").GetValue();

    auto met_pt_vec = df.Take<float>("MET_pt").GetValue();
    auto ht_ak4_vec = df.Take<float>("ht_goodAK4Jets").GetValue();
    auto ht_ak8_vec = df.Take<float>("ht_goodAK8Jets").GetValue();

    auto rdf_entry_vec = df.Take<ULong64_t>("rdfentry_").GetValue();
    auto rdf_slot_vec = df.Take<unsigned int>("rdfslot_").GetValue();

    std::vector<std::shared_ptr<EventDataBase>> events;
    size_t n_events = ak4_pt_vec.size();
    events.reserve(n_events);

    for (size_t i = 0; i < n_events; ++i) {
        auto event = std::make_shared<EventData>();

        event->ak4_pt.assign(ak4_pt_vec[i].begin(), ak4_pt_vec[i].end());
        event->ak4_eta.assign(ak4_eta_vec[i].begin(), ak4_eta_vec[i].end());
        event->ak4_phi.assign(ak4_phi_vec[i].begin(), ak4_phi_vec[i].end());
        event->ak4_mass.assign(ak4_mass_vec[i].begin(), ak4_mass_vec[i].end());
        event->ak4_isTightBTag.assign(ak4_isTightBTag_vec[i].begin(), ak4_isTightBTag_vec[i].end());
        event->ak4_isMediumBTag.assign(ak4_isMediumBTag_vec[i].begin(), ak4_isMediumBTag_vec[i].end());
        event->ak4_isLooseBTag.assign(ak4_isLooseBTag_vec[i].begin(), ak4_isLooseBTag_vec[i].end());

        event->ak8_pt.assign(ak8_pt_vec[i].begin(), ak8_pt_vec[i].end());
        event->ak8_eta.assign(ak8_eta_vec[i].begin(), ak8_eta_vec[i].end());
        event->ak8_phi.assign(ak8_phi_vec[i].begin(), ak8_phi_vec[i].end());
        event->ak8_mass.assign(ak8_mass_vec[i].begin(), ak8_mass_vec[i].end());
        event->ak8_nConstituents.assign(ak8_nConstituents_vec[i].begin(), ak8_nConstituents_vec[i].end());
        event->ak8_HbbScore.assign(ak8_HbbScore_vec[i].begin(), ak8_HbbScore_vec[i].end());
        event->ak8_WqqScore.assign(ak8_WqqScore_vec[i].begin(), ak8_WqqScore_vec[i].end());

        event->met_pt = met_pt_vec[i];
        event->ht_ak4 = ht_ak4_vec[i];
        event->ht_ak8 = ht_ak8_vec[i];

        event->rdf_entry = rdf_entry_vec[i];
        event->rdf_slot = rdf_slot_vec[i];

        events.push_back(event);
    }

    return events;
}

void SPANetInferenceRun2::fillBatchTensors(const std::vector<std::shared_ptr<EventDataBase>>& events_base, size_t actual_batch_size) {
    //std::cout << " --> SPANetInferenceRun2::fillBatchTensors() " << std::endl;

    std::fill(ak4_flat_jets_.begin(), ak4_flat_jets_.begin() + actual_batch_size * MAX_AK4_JETS * AK4_FEATURES, 0.0f);
    std::fill(ak8_flat_jets_.begin(), ak8_flat_jets_.begin() + actual_batch_size * MAX_AK8_JETS * AK8_FEATURES, 0.0f);
    std::fill(ak4_mask_char_.begin(), ak4_mask_char_.begin() + actual_batch_size * MAX_AK4_JETS, 0);
    std::fill(ak8_mask_char_.begin(), ak8_mask_char_.begin() + actual_batch_size * MAX_AK8_JETS, 0);

    for (size_t batch_idx = 0; batch_idx < actual_batch_size; ++batch_idx) {
        auto event_base = events_base[batch_idx];
        auto event = std::dynamic_pointer_cast<EventData>(event_base);

        const size_t max_ak4 = std::min<size_t>(MAX_AK4_JETS, event->ak4_pt.size());
        for (size_t i = 0; i < max_ak4; ++i) {
            ak4_mask_char_[batch_idx * MAX_AK4_JETS + i] = 1;

            const size_t base_idx = batch_idx * MAX_AK4_JETS * AK4_FEATURES + i * AK4_FEATURES;
            ak4_flat_jets_[base_idx]     = std::log(event->ak4_mass[i] + 1.0f);
            ak4_flat_jets_[base_idx + 1] = std::log(event->ak4_pt[i] + 1.0f);
            ak4_flat_jets_[base_idx + 2] = event->ak4_eta[i];
            ak4_flat_jets_[base_idx + 3] = std::sin(event->ak4_phi[i]);
            ak4_flat_jets_[base_idx + 4] = std::cos(event->ak4_phi[i]);
            ak4_flat_jets_[base_idx + 5] = static_cast<float>(event->ak4_isTightBTag[i]);
            ak4_flat_jets_[base_idx + 6] = static_cast<float>(event->ak4_isMediumBTag[i]);
            ak4_flat_jets_[base_idx + 7] = static_cast<float>(event->ak4_isLooseBTag[i]);
        }

        const size_t max_ak8 = std::min<size_t>(MAX_AK8_JETS, event->ak8_pt.size());
        for (size_t i = 0; i < max_ak8; ++i) {
            ak8_mask_char_[batch_idx * MAX_AK8_JETS + i] = 1;

            const size_t base_idx = batch_idx * MAX_AK8_JETS * AK8_FEATURES + i * AK8_FEATURES;
            ak8_flat_jets_[base_idx]     = std::log(event->ak8_mass[i] + 1.0f);
            ak8_flat_jets_[base_idx + 1] = std::log(event->ak8_pt[i] + 1.0f);
            ak8_flat_jets_[base_idx + 2] = event->ak8_eta[i];
            ak8_flat_jets_[base_idx + 3] = std::sin(event->ak8_phi[i]);
            ak8_flat_jets_[base_idx + 4] = std::cos(event->ak8_phi[i]);
            ak8_flat_jets_[base_idx + 5] = static_cast<float>(event->ak8_nConstituents[i]);
            ak8_flat_jets_[base_idx + 6] = event->ak8_HbbScore[i];
            ak8_flat_jets_[base_idx + 7] = event->ak8_WqqScore[i];
        }

        const size_t event_base_idx = batch_idx * MET_FEATURES;
        met_inputs_[event_base_idx]     = std::log(event->met_pt + 1.0f);
        met_inputs_[event_base_idx + 1] = std::log(event->ht_ak4 + 1.0f);
        met_inputs_[event_base_idx + 2] = std::log(event->ht_ak8 + 1.0f);

        met_mask_char_[batch_idx] = 1;
    }
    
}



std::vector<int> SPANetInferenceRun2::assign_all_objects_maxprob(
    std::vector<std::vector<float>> vbs_assignment,
    std::vector<std::vector<float>> h_assignment,
    std::vector<std::vector<float>> bh_assignment,
    std::vector<std::vector<float>> v1_assignment,
    std::vector<std::vector<float>> v2_assignment,
    std::vector<std::vector<float>> bv1_assignment,
    std::vector<std::vector<float>> bv2_assignment,
    float vbs_detection,
    float h_detection,
    float bh_detection,
    float v1_detection,
    float v2_detection,
    float bv1_detection,
    float bv2_detection,
    RVec<float> Jet_eta,
    RVec<float> Jet_phi,
    RVec<float> FatJet_eta,
    RVec<float> FatJet_phi
) {
    bool DEBUG = true;
    if (DEBUG) {
        std::cout << "Event processing" << std::endl;
        std::cout << " Boosted assignments:" << std::endl;
        std::cout << " bh = {";
        for (const auto& v : bh_assignment) {
            if (!v.empty()) std::cout << "[" << v[0] << ", " << static_cast<int>(v[0]) << "], ";
        }
        std::cout << "}" << std::endl;
        std::cout << " bv1 = {";
        for (const auto& v : bv1_assignment) {
            if (!v.empty()) std::cout << "[" << v[0] << ", " << static_cast<int>(v[0]) << "], ";
        }
        std::cout << "}" << std::endl;
        std::cout << " bv2 = {";
        for (const auto& v : bv2_assignment) {
            if (!v.empty()) std::cout << "[" << v[0] << ", " << static_cast<int>(v[0]) << "], ";
        }
        std::cout << "}" << std::endl;
    }

    std::vector<int> result(11, -1);
    std::set<int> used_jets;
    std::set<int> used_fatjets;

   
    // Helper to get sorted candidates for resolved (descending prob)
    auto get_sorted_candidates = [](const std::vector<std::vector<float>>& matrix) -> std::vector<std::vector<float>> {
        std::vector<std::vector<float>> cands;
        size_t n = matrix.size();
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i + 1; j < n; ++j) {
                float prob = matrix[i][j];
                cands.push_back({prob, static_cast<float>(i), static_cast<float>(j)});
            }
        }
        std::sort(cands.begin(), cands.end(), [](const std::vector<float>& a, const std::vector<float>& b) {
            return a[0] > b[0];
        });
        return cands;
    };

    // Helper to get sorted candidates for boosted (descending prob)
    auto get_sorted_boosted = [](const std::vector<std::vector<float>>& vec) -> std::vector<std::vector<float>> {
        std::vector<std::vector<float>> cands;
        size_t n = vec.size();
        for (size_t k = 0; k < n; ++k) {
            if (!vec[k].empty()) {
                float prob = vec[k][0];
                cands.push_back({prob, static_cast<float>(k)});
            }
        }
        std::sort(cands.begin(), cands.end(), [](const std::vector<float>& a, const std::vector<float>& b) {
            return a[0] > b[0];
        });
        return cands;
    };

    // Extract sorted candidates
    auto cand_vbs = get_sorted_candidates(vbs_assignment);
    auto cand_h = get_sorted_candidates(h_assignment);
    auto cand_v1 = get_sorted_candidates(v1_assignment);
    auto cand_v2 = get_sorted_candidates(v2_assignment);
    auto cand_bh = get_sorted_boosted(bh_assignment);
    auto cand_bv1 = get_sorted_boosted(bv1_assignment);
    auto cand_bv2 = get_sorted_boosted(bv2_assignment);

    // Filter for boosted
    auto filter_candidates = [](const std::vector<std::vector<float>>& cand, float th) {
        std::vector<std::vector<float>> filtered;
        for (const auto& c : cand) {
            if (c[0] > th) {
                filtered.push_back(c);
            }
        }
        return filtered;
    };

    float threshold = 0.5f;
    std::vector<std::vector<std::vector<float>>> boosted_cands = {
        filter_candidates(cand_bh, threshold),
        filter_candidates(cand_bv1, threshold),
        filter_candidates(cand_bv2, threshold)
    };

    // Helper to find best fatjet
    auto get_best_fatjet = [](const std::vector<std::vector<float>>& candidates, const std::set<int>& used) -> std::pair<int, float> {
        for (const auto& cand : candidates) {
            int fj = static_cast<int>(cand[1]);
            if (used.find(fj) == used.end()) {
                return {fj, cand[0]}; // assumes candidates are probability ordered
            }
        }
        return {-1, -std::numeric_limits<float>::infinity()};
    };

    // Helper to find best pair
    auto get_best_pair = [](const std::vector<std::vector<float>>& candidates, const std::set<int>& used) -> std::tuple<int, int, float> {
        for (const auto& cand : candidates) {
            int j1 = static_cast<int>(cand[1]);
            int j2 = static_cast<int>(cand[2]);
            if (j1 != j2 && used.find(j1) == used.end() && used.find(j2) == used.end()) {
                return {j1, j2, cand[0]}; // assumes candidates are probability ordered
            }
        }
        return {-1, -1, -std::numeric_limits<float>::infinity()};
    };

    // Boosted assignment logic
    std::vector<int> boosted_result_pos = {4, 9, 10}; // bh, bv1, bv2 positions

    std::vector<float> top_probs(3, -std::numeric_limits<float>::infinity());
    for (int i = 0; i < 3; ++i) {
        if (!boosted_cands[i].empty()) {
            top_probs[i] = boosted_cands[i][0][0];
        }
    }

    int first_idx = -1;
    float max_prob = -std::numeric_limits<float>::infinity();
    for (int i = 0; i < 3; ++i) {
        if (top_probs[i] > max_prob) {
            max_prob = top_probs[i];
            first_idx = i;
        }
    }

    if (first_idx != -1) {
        auto [fj, prob] = get_best_fatjet(boosted_cands[first_idx], used_fatjets);
        if (fj != -1) {
            if (DEBUG) {
                std::cout << "First matched boosted boson: jet idx = " << fj << ", original idx = " << first_idx << std::endl;
            }
            result[boosted_result_pos[first_idx]] = fj;
            used_fatjets.insert(fj);
        }

        std::vector<int> remaining;
        for (int i = 0; i < 3; ++i) {
            if (i != first_idx && !boosted_cands[i].empty()) {
                remaining.push_back(i);
            }
        }

        if (!remaining.empty()) {
            if (remaining.size() == 1) {
                int idx = remaining[0];
                auto [fj, prob] = get_best_fatjet(boosted_cands[idx], used_fatjets);
                if (fj != -1) {
                    if (DEBUG) {
                        std::cout << "Found second matched boosted boson: idx = " << fj << ", original idx = " << idx << std::endl;
                    }
                    result[boosted_result_pos[idx]] = fj;
                    used_fatjets.insert(fj);
                }
            } else {
               std::vector<std::tuple<int, float, int>> best_probs;
                for (int idx : remaining) {
                    auto [fj, prob] = get_best_fatjet(boosted_cands[idx], used_fatjets);
                    best_probs.push_back({idx, prob, fj});
                }

                if (DEBUG) {
                    std::cout << "Potential second boosted bosons probs: {";
                    for (const auto& [idx, prob, fj] : best_probs) {
                        std::cout << prob << ", ";
                    }
                    std::cout << "}" << std::endl;
                }

                // Sort descending by prob
                std::sort(best_probs.begin(), best_probs.end(),
                        [](const std::tuple<int, float, int>& a, const std::tuple<int, float, int>& b) {
                            return std::get<1>(a) > std::get<1>(b);
                        });

                // Assign the best one
                int best_idx = std::get<0>(best_probs[0]);
                int best_fj = std::get<2>(best_probs[0]);
                float best_p = std::get<1>(best_probs[0]);
                if (best_fj != -1) {
                    if (DEBUG) {
                        std::cout << "Found second matched boosted boson: idx = " << best_fj << ", original idx = " << best_idx << std::endl;
                    }
                    result[boosted_result_pos[best_idx]] = best_fj;
                    used_fatjets.insert(best_fj);
                } else if (DEBUG) {
                    std::cout << "No jets left after overlap removal for either boson. Returning." << std::endl;
                }

                // Last one: pick the remaining index
                int last_idx = remaining[best_idx == remaining[0] ? 1 : 0];
                auto [last_fj, last_prob] = get_best_fatjet(boosted_cands[last_idx], used_fatjets);
                if (last_fj != -1) {
                    if (DEBUG) {
                        std::cout << "Found third matched boosted boson: idx = " << last_fj << ", original idx = " << last_idx << std::endl;
                    }
                    result[boosted_result_pos[last_idx]] = last_fj;
                    used_fatjets.insert(last_fj);
                }
            }
        } else if (DEBUG) {
            std::cout << "No boosted boson left to match. Returning." << std::endl;
        }           
    }

    // Resolved
    std::vector<bool> is_boosted(3, false);
    is_boosted[0] = (result[4] != -1);
    is_boosted[1] = (result[9] != -1);
    is_boosted[2] = (result[10] != -1);

    std::vector<std::vector<std::vector<float>>> resolved_cands;
    std::vector<std::pair<int, int>> resolved_pos;

    if (!is_boosted[0]) {
        resolved_cands.push_back(cand_h);
        resolved_pos.emplace_back(2, 3);
    }
    if (!is_boosted[1]) {
        resolved_cands.push_back(cand_v1);
        resolved_pos.emplace_back(5, 6);
    }
    if (!is_boosted[2]) {
        resolved_cands.push_back(cand_v2);
        resolved_pos.emplace_back(7, 8);
    }

    if (DEBUG) {
        std::cout << " Resolved assignments:" << std::endl;
        for (const auto& label : {"h", "v1", "v2"}) {
            const auto& cand = (label == std::string("h") ? cand_h : (label == std::string("v1") ? cand_v1 : cand_v2));
            std::cout << " " << label << " = {";
            for (const auto& c : cand) {
                std::cout << "[" << c[0] << ", " << static_cast<int>(c[1]) << ", " << static_cast<int>(c[2]) << "], ";
            }
            std::cout << "}" << std::endl;
        }
        std::cout << "Resolved bosons to match: {";
        for (const auto& [j1, j2] : resolved_pos) {
            std::string label = (j1 == 2 ? "h" : (j1 == 5 ? "v1" : "v2"));
            std::cout << label << ", ";
        }
        std::cout << "}" << std::endl;
    }

    if (!resolved_cands.empty()) {
        std::vector<float> res_top_probs(resolved_cands.size(), -std::numeric_limits<float>::infinity());
        for (size_t i = 0; i < resolved_cands.size(); ++i) {
            if (!resolved_cands[i].empty()) {
                res_top_probs[i] = resolved_cands[i][0][0]; // candidate i, highest probability jet pair, probability value
            }
        }

        int res_first_idx = -1;
        float res_max_prob = -std::numeric_limits<float>::infinity();
        for (size_t i = 0; i < resolved_cands.size(); ++i) {
            if (res_top_probs[i] > res_max_prob) {
                res_max_prob = res_top_probs[i];
                res_first_idx = static_cast<int>(i);
            }
        }

        if (DEBUG && res_first_idx != -1) {
            std::cout << "Boson with max reconstruction prob: index = " << res_first_idx << std::endl;
        }

        if (res_first_idx != -1) {
            auto [j1, j2, prob] = get_best_pair(resolved_cands[res_first_idx], used_jets);
            if (j1 != -1) {
                if (DEBUG) {
                    std::cout << "Selected pair: (" << j1 << ", " << j2 << ")" << std::endl;
                }
                result[resolved_pos[res_first_idx].first] = j1;
                result[resolved_pos[res_first_idx].second] = j2;
                used_jets.insert(j1);
                used_jets.insert(j2);
            }

            std::vector<int> res_remaining;
            for (size_t i = 0; i < resolved_cands.size(); ++i) {
                if (static_cast<int>(i) != res_first_idx && !resolved_cands[i].empty()) {
                    res_remaining.push_back(static_cast<int>(i));
                }
            }

            if (!res_remaining.empty()) {
                if (res_remaining.size() == 1) {
                    if (DEBUG) {
                        std::cout << "Only one remaining boson." << std::endl;
                    }
                    int idx = res_remaining[0];
                    auto [j1, j2, prob] = get_best_pair(resolved_cands[idx], used_jets);
                    if (j1 != -1) {
                        result[resolved_pos[idx].first] = j1;
                        result[resolved_pos[idx].second] = j2;
                        used_jets.insert(j1);
                        used_jets.insert(j2);
                    }
                }
            } else {
                std::vector<std::tuple<int, float, int, int>> res_best_probs;
                for (int idx : res_remaining) {
                    auto [j1, j2, prob] = get_best_pair(resolved_cands[idx], used_jets);
                    res_best_probs.push_back({idx, prob, j1, j2});
                }
                if (DEBUG) {
                    std::cout << "\nAfter first selection:" << std::endl;
                    std::cout << "Best probs for remaining: {";
                    for (const auto& [idx, prob, j1, j2] : res_best_probs) {
                        std::cout << prob << ", ";
                    }
                    std::cout << "}" << std::endl;
                }

                // Sort descending by prob
                std::sort(res_best_probs.begin(), res_best_probs.end(),
                          [](const std::tuple<int, float, int, int> &a, const std::tuple<int, float, int, int> &b) {
                              return std::get<1>(a) > std::get<1>(b);
                          });

                int best_idx = std::get<0>(res_best_probs[0]);
                int best_j1 = std::get<2>(res_best_probs[0]);
                int best_j2 = std::get<3>(res_best_probs[0]);
                float best_p = std::get<1>(res_best_probs[0]);
                if (best_j1 != -1) {
                    result[resolved_pos[best_idx].first] = best_j1;
                    result[resolved_pos[best_idx].second] = best_j2;
                    used_jets.insert(best_j1);
                    used_jets.insert(best_j2);
                } else if (DEBUG) {
                        std::cout << "No valid pairs for remaining bosons." << std::endl;
                }

                int last_idx = res_remaining[best_idx == res_remaining[0] ? 1 : 0];
                auto [last_j1, last_j2, last_prob] = get_best_pair(resolved_cands[last_idx], used_jets);
                if (last_j1 != -1) {
                    if (DEBUG) {
                        std::cout << "\nAfter second selection:" << std::endl;
                        std::cout << "Last pair: (" << last_j1 << ", " << last_j2 << ")" << std::endl;
                    }
                    result[resolved_pos[last_idx].first] = last_j1;
                    result[resolved_pos[last_idx].second] = last_j2;
                    used_jets.insert(last_j1);
                    used_jets.insert(last_j2);
                }    
            }
        } 
    } else if (DEBUG) {
        std::cout << "No remaining bosons, returning." << std::endl;
    }

    // VBS
    if (DEBUG) {
        std::cout << "\nStarting VBS Jets Selection" << std::endl;
        std::cout << "vbs = {";
        for (const auto& c : cand_vbs) {
            std::cout << "[" << c[0] << ", " << static_cast<int>(c[1]) << ", " << static_cast<int>(c[2]) << "], ";
        }
        std::cout << "}" << std::endl;
    }
    auto [vbs_j1, vbs_j2, vbs_prob] = get_best_pair(cand_vbs, used_jets);
    if (vbs_j1 != -1 && vbs_prob >= 0.0f) {
        if (DEBUG) {
            std::cout << "vbs_best_pair = (" << vbs_j1 << ", " << vbs_j2 << "), vbs_best_p = " << vbs_prob << std::endl;
        }
        result[0] = vbs_j1;
        result[1] = vbs_j2;
    }

    if (DEBUG) {
        std::cout << "matched_bosons (boosted) = {";
        std::vector<std::string> bb_labels = {"bh", "bv1", "bv2"};
        for (int i = 0; i < 3; ++i) {
            int pos = boosted_result_pos[i];
            if (result[pos] != -1) {
                std::cout << bb_labels[i] << ": " << result[pos] << ", ";
            }
        }
        std::cout << "}" << std::endl;
        std::cout << "Selected resolved jet pairs: {";
        for (size_t i = 0; i < resolved_pos.size(); ++i) {
            int j1 = result[resolved_pos[i].first];
            int j2 = result[resolved_pos[i].second];
            if (j1 != -1 && j2 != -1) {
                std::string label = (resolved_pos[i].first == 2 ? "h" : (resolved_pos[i].first == 5 ? "v1" : "v2"));
                std::cout << label << ": (" << j1 << ", " << j2 << "), ";
            }
        }
        std::cout << "}" << std::endl;
        if (result[0] != -1 && result[1] != -1) {
            std::cout << "matched_bosons (vbs) = (" << result[0] << ", " << result[1] << ")" << std::endl;
        }
    }

    return result;
}

RNode SPANetInferenceRun2::ParseSpanetInference(RNode df_) {
    std::cout << "--> SPANetInferenceRun2::ParseSpanetInference()" << std::endl;

    auto df = df_.Define("_all_assignments", assign_all_objects_maxprob, {
        "spanet_vbs_assignment", "spanet_h_assignment", "spanet_bh_assignment",
        "spanet_v1_assignment", "spanet_v2_assignment", "spanet_bv1_assignment", "spanet_bv2_assignment",
        "spanet_vbs_detection", "spanet_h_detection", "spanet_bh_detection",
        "spanet_v1_detection", "spanet_v2_detection", "spanet_bv1_detection", "spanet_bv2_detection",
        "goodAK4Jets_eta", "goodAK4Jets_phi", "goodAK8Jets_eta", "goodAK8Jets_phi"
    })
    .Define("vbs1_idx", "_all_assignments[0]")
    .Define("vbs2_idx", "_all_assignments[1]")
    .Define("h1_idx", "_all_assignments[2]")
    .Define("h2_idx", "_all_assignments[3]")
    .Define("bh_idx", "_all_assignments[4]")
    .Define("v1_j1_idx", "_all_assignments[5]")
    .Define("v1_j2_idx", "_all_assignments[6]")
    .Define("v2_j1_idx", "_all_assignments[7]")
    .Define("v2_j2_idx", "_all_assignments[8]")
    .Define("bv1_idx", "_all_assignments[9]")
    .Define("bv2_idx", "_all_assignments[10]");

    // VBS jet variables
    df = df.Define("spanet_vbs1_pt", "vbs1_idx >= 0 ? goodAK4Jets_pt[vbs1_idx] : -999.0f")
           .Define("spanet_vbs1_eta", "vbs1_idx >= 0 ? goodAK4Jets_eta[vbs1_idx] : -999.0f")
           .Define("spanet_vbs1_phi", "vbs1_idx >= 0 ? goodAK4Jets_phi[vbs1_idx] : -999.0f")
           .Define("spanet_vbs1_mass", "vbs1_idx >= 0 ? goodAK4Jets_mass[vbs1_idx] : -999.0f")
           .Define("spanet_vbs2_pt", "vbs2_idx >= 0 ? goodAK4Jets_pt[vbs2_idx] : -999.0f")
           .Define("spanet_vbs2_eta", "vbs2_idx >= 0 ? goodAK4Jets_eta[vbs2_idx] : -999.0f")
           .Define("spanet_vbs2_phi", "vbs2_idx >= 0 ? goodAK4Jets_phi[vbs2_idx] : -999.0f")
           .Define("spanet_vbs2_mass", "vbs2_idx >= 0 ? goodAK4Jets_mass[vbs2_idx] : -999.0f")
           .Define("spanet_vbs_detajj", "vbs1_idx >= 0 && vbs2_idx >= 0 ? abs(spanet_vbs1_eta - spanet_vbs2_eta) : -999.0f")
           .Define("spanet_vbs_mjj", "vbs1_idx >= 0 && vbs2_idx >= 0 ? (ROOT::Math::PtEtaPhiMVector(spanet_vbs1_pt, spanet_vbs1_eta, spanet_vbs1_phi, spanet_vbs1_mass) + "
                                     "ROOT::Math::PtEtaPhiMVector(spanet_vbs2_pt, spanet_vbs2_eta, spanet_vbs2_phi, spanet_vbs2_mass)).M() : -999.0f");

    // Resolved Higgs variables
    df = df.Define("spanet_h1_pt", "h1_idx >= 0 ? goodAK4Jets_pt[h1_idx] : -999.0f")
           .Define("spanet_h1_eta", "h1_idx >= 0 ? goodAK4Jets_eta[h1_idx] : -999.0f")
           .Define("spanet_h1_phi", "h1_idx >= 0 ? goodAK4Jets_phi[h1_idx] : -999.0f")
           .Define("spanet_h1_mass", "h1_idx >= 0 ? goodAK4Jets_mass[h1_idx] : -999.0f")
           .Define("spanet_h2_pt", "h2_idx >= 0 ? goodAK4Jets_pt[h2_idx] : -999.0f")
           .Define("spanet_h2_eta", "h2_idx >= 0 ? goodAK4Jets_eta[h2_idx] : -999.0f")
           .Define("spanet_h2_phi", "h2_idx >= 0 ? goodAK4Jets_phi[h2_idx] : -999.0f")
           .Define("spanet_h2_mass", "h2_idx >= 0 ? goodAK4Jets_mass[h2_idx] : -999.0f")
           .Define("spanet_h_mjj", "h1_idx >= 0 && h2_idx >= 0 ? (ROOT::Math::PtEtaPhiMVector(spanet_h1_pt, spanet_h1_eta, spanet_h1_phi, spanet_h1_mass) + "
                                   "ROOT::Math::PtEtaPhiMVector(spanet_h2_pt, spanet_h2_eta, spanet_h2_phi, spanet_h2_mass)).M() : -999.0f");

    // Boosted Higgs variables
    df = df.Define("spanet_bh_eta", "bh_idx >= 0 ? goodAK8Jets_eta[bh_idx] : -999.0f")
             .Define("spanet_bh_phi", "bh_idx >= 0 ? goodAK8Jets_phi[bh_idx] : -999.0f")
             .Define("spanet_bh_msoftdrop", "bh_idx >= 0 ? goodAK8Jets_msoftdrop[bh_idx] : -999.0f")
             .Define("spanet_bh_pt", "bh_idx >= 0 ? goodAK8Jets_pt[bh_idx] : -999.0f")
             .Define("spanet_bh_HbbScore", "bh_idx >= 0 ? goodAK8Jets_HbbScore[bh_idx] : -999.0f")
             .Define("spanet_bh_WqqScore", "bh_idx >= 0 ? goodAK8Jets_WqqScore[bh_idx] : -999.0f");

    // Resolved V1 variables
    df = df.Define("spanet_v1_j1_pt", "v1_j1_idx >= 0 ? goodAK4Jets_pt[v1_j1_idx] : -999.0f")
           .Define("spanet_v1_j1_eta", "v1_j1_idx >= 0 ? goodAK4Jets_eta[v1_j1_idx] : -999.0f")
           .Define("spanet_v1_j1_phi", "v1_j1_idx >= 0 ? goodAK4Jets_phi[v1_j1_idx] : -999.0f")
           .Define("spanet_v1_j1_mass", "v1_j1_idx >= 0 ? goodAK4Jets_mass[v1_j1_idx] : -999.0f")
           .Define("spanet_v1_j2_pt", "v1_j2_idx >= 0 ? goodAK4Jets_pt[v1_j2_idx] : -999.0f")
           .Define("spanet_v1_j2_eta", "v1_j2_idx >= 0 ? goodAK4Jets_eta[v1_j2_idx] : -999.0f")
           .Define("spanet_v1_j2_phi", "v1_j2_idx >= 0 ? goodAK4Jets_phi[v1_j2_idx] : -999.0f")
           .Define("spanet_v1_j2_mass", "v1_j2_idx >= 0 ? goodAK4Jets_mass[v1_j2_idx] : -999.0f")
           .Define("spanet_v1_mjj", "v1_j1_idx >= 0 && v1_j2_idx >= 0 ? (ROOT::Math::PtEtaPhiMVector(spanet_v1_j1_pt, spanet_v1_j1_eta, spanet_v1_j1_phi, spanet_v1_j1_mass) + "
                                    "ROOT::Math::PtEtaPhiMVector(spanet_v1_j2_pt, spanet_v1_j2_eta, spanet_v1_j2_phi, spanet_v1_j2_mass)).M() : -999.0f");

    // Resolved V2 variables
    df = df.Define("spanet_v2_j1_pt", "v2_j1_idx >= 0 ? goodAK4Jets_pt[v2_j1_idx] : -999.0f")
           .Define("spanet_v2_j1_eta", "v2_j1_idx >= 0 ? goodAK4Jets_eta[v2_j1_idx] : -999.0f")
           .Define("spanet_v2_j1_phi", "v2_j1_idx >= 0 ? goodAK4Jets_phi[v2_j1_idx] : -999.0f")
           .Define("spanet_v2_j1_mass", "v2_j1_idx >= 0 ? goodAK4Jets_mass[v2_j1_idx] : -999.0f")
           .Define("spanet_v2_j2_pt", "v2_j2_idx >= 0 ? goodAK4Jets_pt[v2_j2_idx] : -999.0f")
           .Define("spanet_v2_j2_eta", "v2_j2_idx >= 0 ? goodAK4Jets_eta[v2_j2_idx] : -999.0f")
           .Define("spanet_v2_j2_phi", "v2_j2_idx >= 0 ? goodAK4Jets_phi[v2_j2_idx] : -999.0f")
           .Define("spanet_v2_j2_mass", "v2_j2_idx >= 0 ? goodAK4Jets_mass[v2_j2_idx] : -999.0f")
           .Define("spanet_v2_mjj", "v2_j1_idx >= 0 && v2_j2_idx >= 0 ? (ROOT::Math::PtEtaPhiMVector(spanet_v2_j1_pt, spanet_v2_j1_eta, spanet_v2_j1_phi, spanet_v2_j1_mass) + "
                                    "ROOT::Math::PtEtaPhiMVector(spanet_v2_j2_pt, spanet_v2_j2_eta, spanet_v2_j2_phi, spanet_v2_j2_mass)).M() : -999.0f");

    // Boosted V1 variables
    df = df.Define("spanet_bv1_eta", "bv1_idx >= 0 ? goodAK8Jets_eta[bv1_idx] : -999.0f")
             .Define("spanet_bv1_phi", "bv1_idx >= 0 ? goodAK8Jets_phi[bv1_idx] : -999.0f")
             .Define("spanet_bv1_msoftdrop", "bv1_idx >= 0 ? goodAK8Jets_msoftdrop[bv1_idx] : -999.0f")
             .Define("spanet_bv1_pt", "bv1_idx >= 0 ? goodAK8Jets_pt[bv1_idx] : -999.0f")
             .Define("spanet_bv1_HbbScore", "bv1_idx >= 0 ? goodAK8Jets_HbbScore[bv1_idx] : -999.0f")
             .Define("spanet_bv1_WqqScore", "bv1_idx >= 0 ? goodAK8Jets_WqqScore[bv1_idx] : -999.0f");

    // Boosted V2 variables
    df = df.Define("spanet_bv2_eta", "bv2_idx >= 0 ? goodAK8Jets_eta[bv2_idx] : -999.0f")
             .Define("spanet_bv2_phi", "bv2_idx >= 0 ? goodAK8Jets_phi[bv2_idx] : -999.0f")
             .Define("spanet_bv2_msoftdrop", "bv2_idx >= 0 ? goodAK8Jets_msoftdrop[bv2_idx] : -999.0f")
             .Define("spanet_bv2_pt", "bv2_idx >= 0 ? goodAK8Jets_pt[bv2_idx] : -999.0f")
             .Define("spanet_bv2_HbbScore", "bv1_idx >= 0 ? goodAK8Jets_HbbScore[bv1_idx] : -999.0f")
             .Define("spanet_bv2_WqqScore", "bv1_idx >= 0 ? goodAK8Jets_WqqScore[bv1_idx] : -999.0f");

    return df;
}

}