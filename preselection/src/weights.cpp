#include "weights.h"

#include <atomic>
#include <cmath>
#include <stdexcept>

namespace {
constexpr double kBTagDenominatorEpsilon = 1.e-8;
std::atomic<unsigned long long> g_btag_negative_intermediate{0};
std::atomic<unsigned long long> g_btag_tiny_denominator{0};
std::atomic<unsigned long long> g_btag_invalid_probability{0};

double unityForInvalidBTagWeight(std::atomic<unsigned long long> &counter) {
    ++counter;
    return 1.;
}

std::string bTagEfficiencyFamily(const std::string &sample) {
    const auto contains = [&sample](const std::string &needle) {
        return sample.find(needle) != std::string::npos;
    };
    // Keep this ordered consistently with misc/sf-utils/btag_eff_families.py.
    if (contains("VBSWWH_") || contains("VBSWZH_") || contains("VBSZZH_")) return "VBS_VVH";
    if (contains("VBS-SSWW")) return "VBS_SSWW";
    if (contains("WWJJ")) return "WWJJ";
    if (contains("ZZJJ")) return "ZZJJ";
    if (contains("GluGluToContinto2Z") || contains("GluGlutoContinto2Z")) return "ggZZ_continuum";
    if (contains("GluGluH-")) return "ggH";
    if (contains("GluGluZH") || contains("WplusH") || contains("WminusH") || contains("ZH-")) return "VH";
    if (contains("TTH-")) return "ttH";
    if (contains("TTWW") || contains("TTWZ") || contains("TTW-") || contains("TZQB")) return "ttV_rare";
    if (contains("TTto") || contains("TTLL") || contains("TTLNu")) return "ttbar";
    if (contains("TWminus") || contains("TbarWplus")) return "tW";
    if (contains("TBbarQ") || contains("TbarBQ") || contains("TBbarto") || contains("TbarBto")) return "single_top";
    if (contains("WWW-") || contains("WWZ-") || contains("WZZ-") || contains("ZZZ-")) return "triboson";
    if (contains("WZ_") || contains("WZto")) return "WZ";
    if (contains("WW_") || contains("WWto")) return "WW";
    if (contains("ZZ_") || contains("ZZto")) return "ZZ";
    if (contains("DYto")) return "DY";
    if (contains("WtoLNu")) return "W_leptonic";
    if (contains("Wto2Q")) return "W_hadronic";
    if (contains("Zto2Q")) return "Z_hadronic";
    if (contains("QCD-4Jets")) return "QCD_HT";
    if (contains("QCD_Bin")) return "QCD_PT";
    throw std::runtime_error("No b-tag efficiency family is configured for sample " + sample);
}
} // namespace

void resetBTagDiagnostics() {
    g_btag_negative_intermediate = 0;
    g_btag_tiny_denominator = 0;
    g_btag_invalid_probability = 0;
}

void printBTagDiagnostics(std::ostream &out) {
    const auto negative = g_btag_negative_intermediate.load();
    const auto tiny_denominator = g_btag_tiny_denominator.load();
    const auto invalid_probability = g_btag_invalid_probability.load();
    if (negative == 0 && tiny_denominator == 0 && invalid_probability == 0) return;
    out << "[BTag SF diagnostics] negative intermediate probabilities=" << negative
        << ", tiny denominators=" << tiny_denominator
        << ", invalid efficiencies/probabilities=" << invalid_probability << '\n';
}

correction::CorrectionSet loadBTagEfficiencyCorrectionSet() {
    return *CorrectionSet::from_file("corrections/scalefactors/btagging/btag_eff.json");
}

/*
############################################
GOLDEN JSON
############################################
*/

RNode applyGoldenJSONWeight(const lumiMask& golden, RNode df){
    auto goldenjson = [&golden](unsigned int &run, unsigned int &luminosityBlock){ return golden.accept(run, luminosityBlock); };
    return df.Filter(goldenjson, {"run", "luminosityBlock"});
}

/*
############################################
PILEUP SFs
############################################
*/

RNode applyPileupScaleFactors(std::unordered_map<std::string, correction::CorrectionSet> cset_pileup, std::unordered_map<std::string, std::string> year_map, RNode df) {
    auto eval_correction = [cset_pileup, year_map] (std::string year, float ntrueint) {
        RVec<double> pileup_weights;
        if (cset_pileup.find(year) == cset_pileup.end()) {
            static std::unordered_set<std::string> warned_years;
            if (warned_years.find(year) == warned_years.end()) {
                std::cout << "Warning: Pileup correction set for year " << year << " not found. Setting pileup weights to 1." << std::endl;
                warned_years.insert(year);
            }
            pileup_weights.push_back(1.0);
            pileup_weights.push_back(1.0);
            pileup_weights.push_back(1.0);
            return pileup_weights;
        }
        auto correctionset = cset_pileup.at(year).at(year_map.at(year));
        pileup_weights.push_back(correctionset->evaluate({ntrueint, "nominal"}));
        pileup_weights.push_back(correctionset->evaluate({ntrueint, "up"}));
        pileup_weights.push_back(correctionset->evaluate({ntrueint, "down"}));
        return pileup_weights;
    };
    return df.Define("weight_pileup", eval_correction, {"year", "Pileup_nTrueInt"});
}

/*
############################################
MUON SFs
############################################
*/

RNode applyMuonIDScaleFactors(std::unordered_map<std::string, correction::CorrectionSet> cset_muon, std::unordered_map<std::string, std::string> year_map, RNode df) {
    auto eval_correction = [cset_muon, year_map] (std::string year, const RVec<float> eta, const RVec<float> pt) {
        RVec<double> muon_sf_weights = {1., 1., 1.};
        if (eta.empty()) {
            return muon_sf_weights;
        }
        if (cset_muon.find(year) == cset_muon.end()) {
            static std::unordered_set<std::string> warned_years;
            if (warned_years.find(year) == warned_years.end()) {
                std::cout << "Warning: Muon ID correction set for year " << year << " not found. Setting muon ID weights to 1." << std::endl;
                warned_years.insert(year);
            }
            return muon_sf_weights;
        }
        auto correctionset = cset_muon.at(year).at(year_map.at(year));
        for (size_t i = 0; i < eta.size(); i++) {
            float pt_min = 15.1;
            float pt_to_pass = std::max(pt[i],pt_min);
            muon_sf_weights[0] *= correctionset->evaluate({abs(eta[i]), pt_to_pass, "nominal"});
            muon_sf_weights[1] *= correctionset->evaluate({abs(eta[i]), pt_to_pass, "systup"});
            muon_sf_weights[2] *= correctionset->evaluate({abs(eta[i]), pt_to_pass, "systdown"});
        }
        return muon_sf_weights;
    };
    return df.Define("weight_muonid", eval_correction, {"year", "muon_eta", "muon_pt"});
}

RNode applyMuonRecoScaleFactors(std::unordered_map<std::string, correction::CorrectionSet> cset_muon, std::unordered_map<std::string, std::string> year_map, RNode df) {
    auto eval_correction = [cset_muon, year_map] (std::string year, const RVec<float> eta, const RVec<float> pt) {
        RVec<double> muon_sf_weights = {1., 1., 1.};
        if (eta.empty()) {
            return muon_sf_weights;
        }
        if (cset_muon.find(year) == cset_muon.end()) {
            static std::unordered_set<std::string> warned_years;
            if (warned_years.find(year) == warned_years.end()) {
                std::cout << "Warning: Muon Reco correction set for year " << year << " not found. Setting muon reco weights to 1." << std::endl;
                warned_years.insert(year);
            }
            return muon_sf_weights;
        }
        auto correctionset = cset_muon.at(year).at(year_map.at(year));
        for (size_t i = 0; i < eta.size(); i++) {
            float pt_max = 15.1;
            float pt_to_pass = std::max(pt[i],pt_max);
            muon_sf_weights[0] *= correctionset->evaluate({abs(eta[i]), pt_to_pass, "nominal"});
            muon_sf_weights[1] *= correctionset->evaluate({abs(eta[i]), pt_to_pass, "systup"});
            muon_sf_weights[2] *= correctionset->evaluate({abs(eta[i]), pt_to_pass, "systdown"});
        }
        return muon_sf_weights;
    };
    return df.Define("weight_muonreco", eval_correction, {"year", "muon_eta", "muon_pt"});
}

RNode applyMuonTriggerScaleFactors(std::unordered_map<std::string, correction::CorrectionSet> cset_muon, std::unordered_map<std::string, std::string> year_map, RNode df) {
    auto eval_correction = [cset_muon, year_map] (std::string year, const RVec<float> eta, const RVec<float> pt) {
        RVec<double> muon_sf_weights = {1., 1., 1.};
        if (eta.empty()) {
            return muon_sf_weights;
        }
        if (cset_muon.find(year) == cset_muon.end()) {
            static std::unordered_set<std::string> warned_years;
            if (warned_years.find(year) == warned_years.end()) {
                std::cout << "Warning: Muon Trigger correction set for year " << year << " not found. Setting muon trigger weights to 1." << std::endl;
                warned_years.insert(year);
            }
            return muon_sf_weights;
        }
        auto correctionset = cset_muon.at(year).at(year_map.at(year));
        for (size_t i = 0; i < eta.size(); i++) {
            float pt_min = 29.1;
            float pt_to_pass = std::max(pt[i],pt_min);
            muon_sf_weights[0] *= correctionset->evaluate({abs(eta[i]), pt_to_pass, "nominal"});
            muon_sf_weights[1] *= correctionset->evaluate({abs(eta[i]), pt_to_pass, "systup"});
            muon_sf_weights[2] *= correctionset->evaluate({abs(eta[i]), pt_to_pass, "systdown"});
        }
        return muon_sf_weights;
    };
    return df.Define("weight_muontrigger", eval_correction, {"year", "muon_eta", "muon_pt"});
}

/*
############################################
ELECTRON SFs
############################################
*/

RNode applyElectronIDScaleFactors(std::unordered_map<std::string, correction::CorrectionSet> cset_electron, std::unordered_map<std::string, std::string> year_map, RNode df) {
    auto eval_correction = [cset_electron, year_map] (std::string year, const RVec<float> eta, const RVec<float> pt) {
        RVec<double> electron_sf_weights = {1., 1., 1.};
        if (eta.empty()) {
            return electron_sf_weights;
        }
        if (cset_electron.find(year) == cset_electron.end()) {
            static std::unordered_set<std::string> warned_years;
            if (warned_years.find(year) == warned_years.end()) {
                std::cout << "Warning: Electron ID correction set for year " << year << " not found. Setting electron ID weights to 1." << std::endl;
                warned_years.insert(year);
            }
            return electron_sf_weights;
        }
        auto correctionset = cset_electron.at(year).at(year_map.at(year));
        for (size_t i = 0; i < eta.size(); i++) {
            electron_sf_weights[0] *= correctionset->evaluate({year, "sf", "Tight", eta[i], pt[i]});
            electron_sf_weights[1] *= correctionset->evaluate({year, "sfup", "Tight", eta[i], pt[i]});
            electron_sf_weights[2] *= correctionset->evaluate({year, "sfdown", "Tight", eta[i], pt[i]});
        }
        return electron_sf_weights;
    };
    return df.Define("weight_electronid", eval_correction, {"year", "electron_SC_eta", "electron_pt"});
}

RNode applyElectronRecoScaleFactors(std::unordered_map<std::string, correction::CorrectionSet> cset_electron, std::unordered_map<std::string, std::string> year_map, RNode df) {
    auto eval_correction = [cset_electron, year_map] (std::string year, const RVec<float> eta, const RVec<float> pt) {
        RVec<double> electron_sf_weights = {1., 1., 1.};
        if (eta.empty()) {
            return electron_sf_weights;
        }
        if (cset_electron.find(year) == cset_electron.end()) {
            static std::unordered_set<std::string> warned_years;
            if (warned_years.find(year) == warned_years.end()) {
                std::cout << "Warning: Electron Reco correction set for year " << year << " not found. Setting electron reco weights to 1." << std::endl;
                warned_years.insert(year);
            }
            return electron_sf_weights;
        }
        auto correctionset = cset_electron.at(year).at(year_map.at(year));

        // Simple way to check if the year refers to Run 2
        bool is_run2 = (year.find("2016") != std::string::npos ||
                        year.find("2017") != std::string::npos ||
                        year.find("2018") != std::string::npos);

        for (size_t i = 0; i < eta.size(); i++) {
            if (is_run2) {
                // Run 2 Scale Factors
                if (pt[i] >= 20) {
                    electron_sf_weights[0] *= correctionset->evaluate({year, "sf", "RecoAbove20", eta[i], pt[i]});
                    electron_sf_weights[1] *= correctionset->evaluate({year, "sfup", "RecoAbove20", eta[i], pt[i]});
                    electron_sf_weights[2] *= correctionset->evaluate({year, "sfdown", "RecoAbove20", eta[i], pt[i]});
                } else {
                    electron_sf_weights[0] *= correctionset->evaluate({year, "sf", "RecoBelow20", eta[i], pt[i]});
                    electron_sf_weights[1] *= correctionset->evaluate({year, "sfup", "RecoBelow20", eta[i], pt[i]});
                    electron_sf_weights[2] *= correctionset->evaluate({year, "sfdown", "RecoBelow20", eta[i], pt[i]});
                }
            } else {
                // Run 3 Scale Factors
                if (pt[i] >= 20 && pt[i] < 75) {
                    electron_sf_weights[0] *= correctionset->evaluate({year, "sf", "Reco20to75", eta[i], pt[i]});
                    electron_sf_weights[1] *= correctionset->evaluate({year, "sfup", "Reco20to75", eta[i], pt[i]});
                    electron_sf_weights[2] *= correctionset->evaluate({year, "sfdown", "Reco20to75", eta[i], pt[i]});
                } else if (pt[i] >= 75) {
                    electron_sf_weights[0] *= correctionset->evaluate({year, "sf", "RecoAbove75", eta[i], pt[i]});
                    electron_sf_weights[1] *= correctionset->evaluate({year, "sfup", "RecoAbove75", eta[i], pt[i]});
                    electron_sf_weights[2] *= correctionset->evaluate({year, "sfdown", "RecoAbove75", eta[i], pt[i]});
                } else {
                    electron_sf_weights[0] *= correctionset->evaluate({year, "sf", "RecoBelow20", eta[i], pt[i]});
                    electron_sf_weights[1] *= correctionset->evaluate({year, "sfup", "RecoBelow20", eta[i], pt[i]});
                    electron_sf_weights[2] *= correctionset->evaluate({year, "sfdown", "RecoBelow20", eta[i], pt[i]});
                }
            }
        }
        return electron_sf_weights;
    };
    return df.Define("weight_electronreco", eval_correction, {"year", "electron_SC_eta", "electron_pt"});
}

RNode applyElectronTriggerScaleFactors(std::unordered_map<std::string, correction::CorrectionSet> cset_electron, std::unordered_map<std::string, std::string> year_map, RNode df) {
    auto eval_correction = [cset_electron, year_map] (std::string year, const RVec<float> eta, const RVec<float> pt) {
        RVec<double> electron_sf_weights = {1., 1., 1.};
        if (eta.empty()) {
            return electron_sf_weights;
        }
        if (cset_electron.find(year) == cset_electron.end()) {
            static std::unordered_set<std::string> warned_years;
            if (warned_years.find(year) == warned_years.end()) {
                std::cout << "Warning: Electron Trigger correction set for year " << year << " not found. Setting electron trigger weights to 1." << std::endl;
                warned_years.insert(year);
            }
            return electron_sf_weights;
        }
        auto correctionset = cset_electron.at(year).at(year_map.at(year));
        for (size_t i = 0; i < eta.size(); i++) {
            electron_sf_weights[0] *= correctionset->evaluate({year, "sf", "HLT_SF_Ele30_TightID", eta[i], pt[i]});
            electron_sf_weights[1] *= correctionset->evaluate({year, "sfup", "HLT_SF_Ele30_TightID", eta[i], pt[i]});
            electron_sf_weights[2] *= correctionset->evaluate({year, "sfdown", "HLT_SF_Ele30_TightID", eta[i], pt[i]});
        }
        return electron_sf_weights;
    };
    return df.Define("weight_electrontrigger", eval_correction, {"year", "electron_SC_eta", "electron_pt"});
}

/*
############################################
BTAG SFs
############################################
*/

RNode applyBTaggingScaleFactors(std::unordered_map<std::string, correction::CorrectionSet> cset_btag,
                                 std::unordered_map<std::string, std::string> corrname_map_HF,
                                 std::unordered_map<std::string, std::string> corrname_map_LF,
                                 const std::string &channel, RNode df) {

    auto calc_btag_sf = [cset_btag, channel] (const std::string &year, const std::string &sample,
                                               const RVec<float> &eta, const RVec<float> &pt,
                                               const RVec<unsigned char> &jetflavor, const RVec<bool> &is_tight,
                                               const RVec<bool> &is_loose, const std::string &correction_name,
                                               bool is_heavy_flavor) {
        RVec<double> btag_sf_weights = {1., 1., 1., 1., 1.};
        if (eta.empty()) return btag_sf_weights;
        if (eta.size() != pt.size() || eta.size() != jetflavor.size() ||
            eta.size() != is_tight.size() || eta.size() != is_loose.size())
            throw std::runtime_error("B-tag input collections have inconsistent sizes");

        const auto cset_it = cset_btag.find(year);
        const auto cset_eff_it = cset_btag.find("eff");
        if (cset_it == cset_btag.end() || cset_eff_it == cset_btag.end())
            throw std::runtime_error("B-tag SF or efficiency correction set is unavailable for year " + year);

        const std::string efficiency_name = "btag_" + year + "_" + channel;
        const std::string legacy_efficiency_name = "btag_" + year;
        bool use_legacy_efficiency = false;
        std::string efficiency_sample = sample;
        decltype(cset_eff_it->second.at(efficiency_name)) efficiency;
        try {
            efficiency = cset_eff_it->second.at(efficiency_name);
        } catch (const std::exception &) {
            try {
                efficiency = cset_eff_it->second.at(legacy_efficiency_name);
                use_legacy_efficiency = true;
            } catch (const std::exception &) {
                throw std::runtime_error("B-tag efficiency map " + efficiency_name +
                                         " is unavailable; run --btag_eff and convert the matching MC sample first");
            }
        }
        if (!use_legacy_efficiency) {
            std::string family_key;
            try {
                family_key = bTagEfficiencyFamily(sample);
            } catch (const std::exception &) {
                // Exact-sample payloads are valid even when no family has been
                // configured.  Do not substitute an inclusive default.
            }
            const auto has_entry = [&](const std::string &key) {
                try {
                    // Any in-range point tests the sample category; the maps
                    // clamp kinematics and all entries share the same binning.
                    (void) efficiency->evaluate({key, "B", "T", 30., 0.});
                    return true;
                } catch (const std::exception &) {
                    return false;
                }
            };
            if (!family_key.empty() && has_entry(family_key)) {
                efficiency_sample = family_key;
            } else if (has_entry(sample)) {
                efficiency_sample = sample;
            } else {
                throw std::runtime_error("B-tag efficiency entries are unavailable for " +
                                         year + ":" + channel + ":" + sample +
                                         "; attempted family key " +
                                         (family_key.empty() ? std::string("<unconfigured>") : family_key) +
                                         " and exact key " + sample);
            }
        }

        const auto &sf_correction = cset_it->second.at(correction_name);
        auto jet_weight = [&](double sf_tight, double sf_loose, double eff_tight,
                              double eff_loose, bool tight, bool loose) {
            if (!std::isfinite(sf_tight) || !std::isfinite(sf_loose) ||
                !std::isfinite(eff_tight) || !std::isfinite(eff_loose) ||
                !(0. <= eff_tight && eff_tight <= eff_loose && eff_loose <= 1.))
                return unityForInvalidBTagWeight(g_btag_invalid_probability);
            const double q_tight = sf_tight * eff_tight;
            const double q_loose = sf_loose * eff_loose;
            if (!std::isfinite(q_tight) || !std::isfinite(q_loose) ||
                !(0. <= q_tight && q_tight <= q_loose && q_loose <= 1.))
                return unityForInvalidBTagWeight(g_btag_invalid_probability);
            if (tight) return sf_tight;
            if (loose) {
                const double denominator = eff_loose - eff_tight;
                const double numerator = q_loose - q_tight;
                if (numerator < 0.) return unityForInvalidBTagWeight(g_btag_negative_intermediate);
                if (std::abs(denominator) < kBTagDenominatorEpsilon)
                    return unityForInvalidBTagWeight(g_btag_tiny_denominator);
                return numerator / denominator;
            }
            const double denominator = 1. - eff_loose;
            const double numerator = 1. - q_loose;
            if (numerator < 0.) return unityForInvalidBTagWeight(g_btag_invalid_probability);
            if (std::abs(denominator) < kBTagDenominatorEpsilon)
                return unityForInvalidBTagWeight(g_btag_tiny_denominator);
            return numerator / denominator;
        };

        for (std::size_t i = 0; i < eta.size(); ++i) {
            bool process_jet = false;
            const char *flavor_label = nullptr;
            const int flavor = std::abs(jetflavor[i]);
            if (is_heavy_flavor && flavor == 5) {
                process_jet = true;
                flavor_label = "B";
            } else if (is_heavy_flavor && flavor == 4) {
                process_jet = true;
                flavor_label = "C";
            } else if (!is_heavy_flavor && flavor != 5 && flavor != 4) {
                process_jet = true;
                flavor_label = "L";
            }
            if (!process_jet || std::abs(eta[i]) >= 2.5f) continue;
            if (is_tight[i] && !is_loose[i]) {
                unityForInvalidBTagWeight(g_btag_invalid_probability);
                continue;
            }

            double eff_tight = 0.;
            double eff_loose = 0.;
            try {
                if (use_legacy_efficiency) {
                    eff_tight = efficiency->evaluate({flavor_label, "T", pt[i], eta[i]});
                    eff_loose = efficiency->evaluate({flavor_label, "L", pt[i], eta[i]});
                } else {
                    eff_tight = efficiency->evaluate({efficiency_sample, flavor_label, "T", pt[i], eta[i]});
                    eff_loose = efficiency->evaluate({efficiency_sample, flavor_label, "L", pt[i], eta[i]});
                }
            } catch (const std::exception &) {
                throw std::runtime_error("B-tag efficiency entries are unavailable for " +
                                         year + ":" + channel + ":" + sample);
            }
            const double abs_eta = std::abs(eta[i]);
            const auto evaluate_sf = [&](const char *variation, const char *wp) {
                return sf_correction->evaluate({variation, wp, flavor, abs_eta, pt[i]});
            };

            btag_sf_weights[0] *= jet_weight(evaluate_sf("central", "T"), evaluate_sf("central", "L"),
                                              eff_tight, eff_loose, is_tight[i], is_loose[i]);
            btag_sf_weights[1] *= jet_weight(evaluate_sf("up_uncorrelated", "T"), evaluate_sf("up_uncorrelated", "L"),
                                              eff_tight, eff_loose, is_tight[i], is_loose[i]);
            btag_sf_weights[2] *= jet_weight(evaluate_sf("down_uncorrelated", "T"), evaluate_sf("down_uncorrelated", "L"),
                                              eff_tight, eff_loose, is_tight[i], is_loose[i]);
            btag_sf_weights[3] *= jet_weight(evaluate_sf("up_correlated", "T"), evaluate_sf("up_correlated", "L"),
                                              eff_tight, eff_loose, is_tight[i], is_loose[i]);
            btag_sf_weights[4] *= jet_weight(evaluate_sf("down_correlated", "T"), evaluate_sf("down_correlated", "L"),
                                              eff_tight, eff_loose, is_tight[i], is_loose[i]);
        }
        return btag_sf_weights;
    };

    auto eval_HF = [cset_btag, corrname_map_HF, calc_btag_sf] (const std::string& year, const std::string &sample,
                                                                 const RVec<float>& eta, const RVec<float>& pt,
                                                                 const RVec<unsigned char>& jetflavor, const RVec<bool> &is_tight,
                                                                 const RVec<bool> &is_loose) {
        auto corrname_it = corrname_map_HF.find(year);
        if (corrname_it == corrname_map_HF.end())
            throw std::runtime_error("No UParTAK4 HF b-tag SF payload is configured for year " + year);
        return calc_btag_sf(year, sample, eta, pt, jetflavor, is_tight, is_loose, corrname_it->second, true);
    };

    auto eval_LF = [cset_btag, corrname_map_LF, calc_btag_sf] (const std::string& year, const std::string &sample,
                                                                 const RVec<float>& eta, const RVec<float>& pt,
                                                                 const RVec<unsigned char>& jetflavor, const RVec<bool> &is_tight,
                                                                 const RVec<bool> &is_loose) {
        auto corrname_it = corrname_map_LF.find(year);
        if (corrname_it == corrname_map_LF.end())
            throw std::runtime_error("No UParTAK4 LF b-tag SF payload is configured for year " + year);
        return calc_btag_sf(year, sample, eta, pt, jetflavor, is_tight, is_loose, corrname_it->second, false);
    };

    auto select_uncorrelated = [] (const RVec<double> &weights) {
        return RVec<double>{weights[0], weights[1], weights[2]};
    };
    auto select_correlated = [] (const RVec<double> &weights) {
        return RVec<double>{weights[0], weights[3], weights[4]};
    };
    return df.Define("btagging_sf_HF_all", eval_HF,
                     {"year", "name", "jet_eta", "jet_pt", "jet_hadronFlavour", "jet_isTightBTag", "jet_isLooseBTag"})
             .Define("btagging_sf_LF_all", eval_LF,
                     {"year", "name", "jet_eta", "jet_pt", "jet_hadronFlavour", "jet_isTightBTag", "jet_isLooseBTag"})
             .Define("weight_btagging_sf_HF_uncorrelated", select_uncorrelated, {"btagging_sf_HF_all"})
             .Define("weight_btagging_sf_HF_correlated", select_correlated, {"btagging_sf_HF_all"})
             .Define("weight_btagging_sf_LF_uncorrelated", select_uncorrelated, {"btagging_sf_LF_all"})
             .Define("weight_btagging_sf_LF_correlated", select_correlated, {"btagging_sf_LF_all"})
             // Keep the historical names as aliases for the uncorrelated variation.
             .Define("weight_btagging_sf_HF", select_uncorrelated, {"btagging_sf_HF_all"})
             .Define("weight_btagging_sf_LF", select_uncorrelated, {"btagging_sf_LF_all"});
}

/*
############################################
OTHER SFs
############################################
*/

// See https://github.com/cmstas/run3-vbsvvh/pull/28#issuecomment-3820814039
RNode applyEWKCorrections(correction::CorrectionSet cset_ewk, RNode df){
    auto eval_correction = [cset_ewk] (RVec<float> LHEPart_pt, RVec<float> LHEPart_eta, RVec<float> LHEPart_phi, RVec<float> LHEPart_mass, RVec<int> LHEPart_pdgId, int do_ewk_corr) {
        if(do_ewk_corr == 0) return 1.;
        else{
            TLorentzVector TEWKq1, TEWKq2, TEWKlep, TEWKnu;
            TEWKq1.SetPtEtaPhiM(LHEPart_pt[4],LHEPart_eta[4],LHEPart_phi[4],LHEPart_mass[4]);
            TEWKq2.SetPtEtaPhiM(LHEPart_pt[5],LHEPart_eta[5],LHEPart_phi[5],LHEPart_mass[5]);
            TEWKlep.SetPtEtaPhiM(LHEPart_pt[2],LHEPart_eta[2],LHEPart_phi[2],LHEPart_mass[2]);
            TEWKnu.SetPtEtaPhiM(LHEPart_pt[3],LHEPart_eta[3],LHEPart_phi[3],LHEPart_mass[3]);
            int chargequark[7] = {0,-1,2,-1,2,-1,2};
            int EWKpdgq1 = LHEPart_pdgId[4];
            int EWKpdgq2 = LHEPart_pdgId[5];
            int EWKsignq1 = (EWKpdgq1 > 0) - (EWKpdgq1 < 0);
            int EWKsignq2 = (EWKpdgq2 > 0) - (EWKpdgq2 < 0);
            double EWKMass_q12 = (TEWKq1 + TEWKq2).M();
            double EWKMass_lnu = (TEWKlep + TEWKnu).M();
            double fabscharge=(fabs((double)(EWKsignq1 * chargequark[abs(EWKpdgq1)] + (EWKsignq2 * chargequark[abs(EWKpdgq2)]))))/3;
            double EWKbjet_pt = -999;
            if(fabscharge ==1){
                if( EWKMass_q12 >= 70 && EWKMass_q12 < 90  && 
                    EWKMass_lnu >= 70 && EWKMass_lnu < 90){
                    return 0.;
                }
            }
            if(EWKMass_q12 >= 95){
                if( abs(EWKpdgq1) == 5 && abs(EWKpdgq2) == 5){
                    if(TEWKq1.Pt() > TEWKq2.Pt())  EWKbjet_pt = TEWKq1.Pt();
                    else                           EWKbjet_pt = TEWKq2.Pt();
                }else if(abs(EWKpdgq1) == 5){
                    EWKbjet_pt = TEWKq1.Pt();
                }else if(abs(EWKpdgq2) == 5){
                    EWKbjet_pt = TEWKq2.Pt();
                }
            }
            if(EWKbjet_pt > -998){
                return cset_ewk.at("EWK")->evaluate({EWKbjet_pt});
            }
            else return 1.;
        }
    };
    return df.Define("weight_ewk", eval_correction, {"LHEPart_pt", "LHEPart_eta", "LHEPart_phi", "LHEPart_mass", "LHEPart_pdgId", "do_ewk_corr"});
}

RNode applyL1PreFiringReweighting(RNode df){
    auto eval_correction = [] (float L1prefire, float L1prefireup, float L1prefiredown) {
        return RVec<float>{L1prefire, L1prefireup, L1prefiredown};
    };
    // TODO: check what this is in v15
    // return df.Define("weight_l1prefiring", eval_correction, {"L1PreFiringWeight_Nom", "L1PreFiringWeight_Up", "L1PreFiringWeight_Dn"});
    return df;
}

RNode applyPSWeight_FSR(RNode df) {
    auto eval_correction = [] (const RVec<float> PSWeight) {
        return RVec<float>{1., PSWeight[1], PSWeight[3]};
    };
    return df.Define("weight_PSFSR", eval_correction, {"PSWeight"});
}

RNode applyPSWeight_ISR(RNode df) {
    auto eval_correction = [] (const RVec<float> PSWeight) {
        return RVec<float>{1., PSWeight[0], PSWeight[2]};
    };
    return df.Define("weight_PSISR", eval_correction, {"PSWeight"});
}

RNode applyLHEScaleWeight_muF(RNode df) {
    auto eval_correction = [] (const RVec<float> LHEScaleWeight) {
        return RVec<float>{1., LHEScaleWeight[5], LHEScaleWeight[3]};
    };
    return df.Define("weight_muF", eval_correction, {"LHEScaleWeight"});
}

RNode applyLHEScaleWeight_muR(RNode df) {
    auto eval_correction = [] (const RVec<float> LHEScaleWeight) {
        return RVec<float>{1., LHEScaleWeight[7], LHEScaleWeight[1]};
    };
    return df.Define("weight_muR", eval_correction, {"LHEScaleWeight"});
}

RNode applyDataWeights(RNode df_) {
    return applyGoldenJSONWeight(LumiMask, df_);
}

RNode applyMCWeights(RNode df_, const std::string &channel, bool apply_btag_sf) {
    // Check for LHE branches (not present in all samples, e.g. QCD)
    auto colNames = df_.GetColumnNames();
    auto hasColumn = [&colNames](const std::string& name) {
        return std::find(colNames.begin(), colNames.end(), name) != colNames.end();
    };
    bool hasLHEScale = hasColumn("LHEScaleWeight");
    bool hasLHEPart = hasColumn("LHEPart_pt");

    auto df = applyPileupScaleFactors(pileupScaleFactors, pileupScaleFactors_yearmap, df_);
    df = applyMuonIDScaleFactors(muonScaleFactors, muonIDScaleFactors_yearmap, df);
    df = applyMuonRecoScaleFactors(muonScaleFactors, muonRecoScaleFactors_yearmap, df);
    df = applyMuonTriggerScaleFactors(muonScaleFactors, muonTriggerScaleFactors_yearmap, df);

    df = applyElectronIDScaleFactors(electronScaleFactors, electronScaleFactors_yearmap, df);
    df = applyElectronRecoScaleFactors(electronScaleFactors, electronScaleFactors_yearmap, df);
    df = applyElectronTriggerScaleFactors(electronTriggerScaleFactors, electronTriggerScaleFactors_yearmap, df);

    if (apply_btag_sf) {
        auto btag_corrections = bTaggingScaleFactors;
        btag_corrections.emplace("eff", loadBTagEfficiencyCorrectionSet());
        resetBTagDiagnostics();
        df = applyBTaggingScaleFactors(std::move(btag_corrections), bTaggingScaleFactors_HF_corrname, bTaggingScaleFactors_LF_corrname, channel, df);
    } else {
        auto unit_btag_weight = [] () { return RVec<double>{1., 1., 1.}; };
        df = df.Define("weight_btagging_sf_HF", unit_btag_weight)
               .Define("weight_btagging_sf_LF", unit_btag_weight)
               .Define("weight_btagging_sf_HF_uncorrelated", unit_btag_weight)
               .Define("weight_btagging_sf_HF_correlated", unit_btag_weight)
               .Define("weight_btagging_sf_LF_uncorrelated", unit_btag_weight)
               .Define("weight_btagging_sf_LF_correlated", unit_btag_weight);
    }

    if (hasLHEPart) {
        df = applyEWKCorrections(cset_ewk, df);
    } else {
        df = df.Define("weight_ewk", [] () { return 1.; }, {});
    }

    df = applyL1PreFiringReweighting(df);
    df = applyPSWeight_FSR(df);
    df = applyPSWeight_ISR(df);

    if (hasLHEScale) {
        df = applyLHEScaleWeight_muF(df);
        df = applyLHEScaleWeight_muR(df);
    } else {
        df = df.Define("weight_muF", [] () { return RVec<float>{1.f, 1.f, 1.f}; }, {});
        df = df.Define("weight_muR", [] () { return RVec<float>{1.f, 1.f, 1.f}; }, {});
    }

    return df.Redefine("weight",
        "weight *"
        "weight_pileup[0] * "
        "weight_muonid[0] * "
        "weight_muonreco[0] * "
        "weight_muontrigger[0] * "
        "weight_electronid[0] * "
        "weight_electronreco[0] * "
        "weight_electrontrigger[0] * "
        "weight_btagging_sf_HF[0] * "
        "weight_btagging_sf_LF[0] * "
        "weight_ewk * "
        // "weight_l1prefiring[0] * "
        "weight_PSISR[0] * "
        "weight_PSFSR[0] * "
        "weight_muF[0] * "
        "weight_muR[0]");
}
