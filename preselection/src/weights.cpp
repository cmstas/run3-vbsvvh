#include "weights.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>

namespace {
constexpr double kBTagDenominatorEpsilon = 1.e-8;
std::atomic<unsigned long long> g_btag_negative_intermediate{0};
std::atomic<unsigned long long> g_btag_tiny_denominator{0};
std::atomic<unsigned long long> g_btag_invalid_probability{0};

constexpr const char *kBTagFamilyConfig = "corrections/scalefactors/btagging/btag_eff_families.yaml";

struct BTagFamilyConfig {
    std::vector<std::pair<std::string, std::vector<std::string>>> preliminary;
    std::map<std::string, std::vector<std::string>> final_samples;
    std::map<std::string, std::vector<std::string>> final_channels;
    std::vector<std::string> excluded_channels;
};

std::string trim(std::string value) {
    const auto begin = value.find_first_not_of(" \t");
    if (begin == std::string::npos) return "";
    const auto end = value.find_last_not_of(" \t");
    return value.substr(begin, end - begin + 1);
}

const BTagFamilyConfig &bTagFamilyConfig() {
    static const BTagFamilyConfig config = [] {
        std::ifstream input(kBTagFamilyConfig);
        if (!input) throw std::runtime_error("Cannot read b-tag family configuration " + std::string(kBTagFamilyConfig));
        BTagFamilyConfig parsed;
        std::set<std::string> preliminary_names, excluded_names;
        std::string line, section, final_kind, current_group;
        const auto fail = [](const std::string &message) {
            throw std::runtime_error("Invalid canonical b-tag family YAML: " + message);
        };
        while (std::getline(input, line)) {
            if (line.find('\t') != std::string::npos) fail("tabs are not supported");
            const auto comment = line.find('#');
            if (comment != std::string::npos) line.erase(comment);
            const auto content = trim(line);
            if (content.empty()) continue;
            const auto indent = line.find_first_not_of(" \t");
            if (indent == 0) {
                if (content != "preliminary_families:" && content != "final_merges:" &&
                    content != "excluded_source_channels:") fail("unknown top-level key " + content);
                section = content.substr(0, content.size() - 1);
                final_kind.clear();
                current_group.clear();
                continue;
            }
            if (section == "preliminary_families") {
                if (indent == 2 && content.back() == ':') {
                    current_group = content.substr(0, content.size() - 1);
                    if (current_group.empty() || !preliminary_names.insert(current_group).second)
                        fail("duplicate or empty preliminary family");
                    parsed.preliminary.emplace_back(current_group, std::vector<std::string>{});
                } else if (indent == 4 && content.rfind("- ", 0) == 0 && !parsed.preliminary.empty()) {
                    parsed.preliminary.back().second.push_back(trim(content.substr(2)));
                } else fail("invalid preliminary_families indentation or syntax");
                continue;
            }
            if (section == "excluded_source_channels") {
                if (indent != 2 || content.rfind("- ", 0) != 0) fail("invalid excluded_source_channels entry");
                const auto channel = trim(content.substr(2));
                if (channel.empty() || !excluded_names.insert(channel).second) fail("duplicate or empty excluded channel");
                parsed.excluded_channels.push_back(channel);
                continue;
            }
            if (section != "final_merges") fail("content outside a supported YAML section");
            if (indent == 2 && content.back() == ':') {
                final_kind = content.substr(0, content.size() - 1);
                if (final_kind != "samples" && final_kind != "channels") fail("unknown final_merges kind " + final_kind);
                current_group.clear();
            } else if (indent == 4 && content.back() == ':' && !final_kind.empty()) {
                current_group = trim(content.substr(0, content.size() - 1));
                auto &groups = final_kind == "samples" ? parsed.final_samples : parsed.final_channels;
                if (current_group.empty() || !groups.emplace(current_group, std::vector<std::string>{}).second)
                    fail("duplicate or empty final group");
            } else if (indent == 6 && content.rfind("- ", 0) == 0 && !current_group.empty()) {
                auto &groups = final_kind == "samples" ? parsed.final_samples : parsed.final_channels;
                const auto member = trim(content.substr(2));
                if (member.empty()) fail("empty final group member");
                groups.at(current_group).push_back(member);
            } else fail("invalid final_merges indentation or syntax");
        }
        if (parsed.preliminary.empty() || parsed.final_samples.empty() || parsed.final_channels.empty() ||
            parsed.excluded_channels.empty()) fail("missing required non-empty mapping");
        for (const auto &[family, needles] : parsed.preliminary)
            if (needles.empty() || std::any_of(needles.begin(), needles.end(), [](const auto &x) { return x.empty(); }))
                fail("empty preliminary family " + family);
        const auto validate_groups = [&fail](const auto &groups, const std::set<std::string> &expected, const char *kind) {
            std::map<std::string, unsigned int> membership;
            for (const auto &[group, members] : groups) {
                if (members.empty()) fail(std::string("empty final ") + kind + " group " + group);
                for (const auto &member : members) ++membership[member];
            }
            std::set<std::string> observed;
            for (const auto &[member, count] : membership) {
                if (count != 1) fail(std::string("duplicate final ") + kind + " member " + member);
                observed.insert(member);
            }
            if (observed != expected) fail(std::string("final ") + kind + " membership does not match canonical sources");
        };
        validate_groups(parsed.final_samples, preliminary_names, "sample");
        std::set<std::string> retained_channels;
        for (const auto &[_, members] : parsed.final_channels)
            retained_channels.insert(members.begin(), members.end());
        validate_groups(parsed.final_channels, retained_channels, "channel");
        for (const auto &channel : parsed.excluded_channels)
            if (retained_channels.count(channel)) fail("excluded channel is also retained: " + channel);
        return parsed;
    }();
    return config;
}

double unityForInvalidBTagWeight(std::atomic<unsigned long long> &counter) {
    ++counter;
    return 1.;
}

std::string finalGroup(const std::map<std::string, std::vector<std::string>> &groups,
                       const std::string &name, const std::string &kind) {
    std::string match;
    for (const auto &[group, members] : groups) {
        if (std::find(members.begin(), members.end(), name) == members.end()) continue;
        if (!match.empty()) throw std::runtime_error(name + " occurs in multiple final b-tag " + kind + " groups");
        match = group;
    }
    if (match.empty()) throw std::runtime_error(name + " is not assigned to a final b-tag " + kind + " group");
    return match;
}

std::string bTagEfficiencyFamily(const std::string &sample) {
    const auto &config = bTagFamilyConfig();
    for (const auto &[family, needles] : config.preliminary) {
        for (const auto &needle : needles) {
            if (sample.find(needle) != std::string::npos)
                return finalGroup(config.final_samples, family, "sample");
        }
    }
    throw std::runtime_error("No preliminary b-tag efficiency family is configured for sample " + sample);
}

std::string bTagEfficiencyChannel(const std::string &channel) {
    return finalGroup(bTagFamilyConfig().final_channels, channel, "channel");
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

        const bool is_legacy_year = (year == "2016preVFP" || year == "2016postVFP" ||
                                     year == "2017" || year == "2018");
        std::string efficiency_channel;
        if (!is_legacy_year) {
            try {
                efficiency_channel = bTagEfficiencyChannel(channel);
            } catch (const std::exception &error) {
                throw std::runtime_error("B-tag final channel mapping is unavailable for year=" + year +
                                         ", requested_channel=" + channel + ": " + error.what());
            }
        }
        const std::string efficiency_name = "btag_" + year + "_" + efficiency_channel;
        const std::string legacy_efficiency_name = "btag_" + year;
        const bool use_legacy_efficiency = is_legacy_year;
        std::string efficiency_sample;
        if (!use_legacy_efficiency) {
            try {
                efficiency_sample = bTagEfficiencyFamily(sample);
            } catch (const std::exception &error) {
                throw std::runtime_error("B-tag final sample mapping is unavailable for year=" + year +
                                         ", requested_channel=" + channel + ", sample=" + sample +
                                         ": " + error.what());
            }
        }
        decltype(cset_eff_it->second.at(efficiency_name)) efficiency;
        try {
            efficiency = cset_eff_it->second.at(use_legacy_efficiency ? legacy_efficiency_name : efficiency_name);
        } catch (const std::exception &) {
            throw std::runtime_error("B-tag efficiency correction is unavailable for year=" + year +
                                     ", requested_channel=" + channel + ", final_channel=" + efficiency_channel +
                                     ", requested_correction=" + (use_legacy_efficiency ? legacy_efficiency_name : efficiency_name));
        }
        if (!use_legacy_efficiency) {
            try {
                // Strict final-payload contract: no preliminary/exact/inclusive fallback.
                (void) efficiency->evaluate({efficiency_sample, "B", "T", 30., 0.});
            } catch (const std::exception &) {
                throw std::runtime_error("B-tag efficiency sample key is unavailable for year=" + year +
                                         ", requested_channel=" + channel + ", final_channel=" + efficiency_channel +
                                         ", sample=" + sample + ", final_sample=" + efficiency_sample +
                                         ", correction=" + efficiency_name);
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
