#include "wqq_sf.h"
#include "category_weight.h"

#include <algorithm>
#include <array>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>

#include "ROOT/RVec.hxx"
#include "correction.h"

using ROOT::VecOps::RVec;

namespace {

constexpr const char *kConfigPath = "applywqqsf.yaml";
constexpr const char *kPayloadDir = "corrections/scalefactors/WqqSF/";
const std::array<std::string, 6> kYears = {"2016preVFP", "2016postVFP", "2017", "2018", "2024Prompt", "2025"};
const std::array<std::string, 3> kSystematics = {"nom", "up", "down"};
constexpr double kDenominatorEpsilon = 1.e-8;
// Probe pt for checking at graph construction that every requested WP is in range.
constexpr double kProbePt = 300.;

struct WqqSFConfig {
    std::vector<double> cut;
    std::vector<std::vector<double>> flags;
    bool empty() const { return cut.empty() && flags.empty(); }
};

struct WqqSFPayload {
    correction::Correction::Ref sf;
    correction::Correction::Ref eff;  // null when the file has no Wqq_eff_MC
};

std::string trim(std::string value) {
    const auto begin = value.find_first_not_of(" \t");
    if (begin == std::string::npos) return "";
    const auto end = value.find_last_not_of(" \t");
    return value.substr(begin, end - begin + 1);
}

// "[0.4, 0.6]" -> {0.4, 0.6}; "[]" -> {}.
std::vector<double> parseNumberList(const std::string &text, const std::string &line) {
    const auto value = trim(text);
    if (value.size() < 2 || value.front() != '[' || value.back() != ']')
        throw std::runtime_error("Invalid applywqqsf.yaml list: " + line);
    std::vector<double> out;
    std::stringstream items(value.substr(1, value.size() - 2));
    std::string item;
    while (std::getline(items, item, ',')) {
        item = trim(item);
        if (item.empty()) continue;
        std::size_t used = 0;
        double wp = 0.;
        try { wp = std::stod(item, &used); } catch (const std::exception &) { used = 0; }
        if (used != item.size()) throw std::runtime_error("Invalid working point '" + item + "' in applywqqsf.yaml: " + line);
        out.push_back(wp);
    }
    return out;
}

// "[[0.6], [0.4, 0.6, 0.8]]" -> {{0.6}, {0.4, 0.6, 0.8}}; "[]" -> {}.
std::vector<std::vector<double>> parseNestedList(const std::string &text, const std::string &line) {
    const auto value = trim(text);
    if (value.size() < 2 || value.front() != '[' || value.back() != ']')
        throw std::runtime_error("Invalid applywqqsf.yaml list: " + line);
    std::vector<std::vector<double>> out;
    const auto inner = value.substr(1, value.size() - 2);
    std::size_t pos = 0;
    while (true) {
        const auto open = inner.find('[', pos);
        if (open == std::string::npos) {
            if (!trim(inner.substr(pos)).empty() && trim(inner.substr(pos)) != ",")
                throw std::runtime_error("Invalid applywqqsf.yaml list: " + line);
            break;
        }
        const auto between = trim(inner.substr(pos, open - pos));
        if (!(between.empty() || (between == "," && !out.empty())))
            throw std::runtime_error("Invalid applywqqsf.yaml list: " + line);
        const auto close = inner.find(']', open);
        if (close == std::string::npos) throw std::runtime_error("Invalid applywqqsf.yaml list: " + line);
        out.push_back(parseNumberList(inner.substr(open, close - open + 1), line));
        pos = close + 1;
    }
    return out;
}

void validateWorkingPoints(const std::vector<double> &wps, const std::string &what, const std::string &channel) {
    if (wps.empty()) throw std::runtime_error("Empty " + what + " working-point set for channel " + channel + " in applywqqsf.yaml");
    for (std::size_t i = 1; i < wps.size(); ++i)
        if (!(wps[i - 1] < wps[i]))
            throw std::runtime_error(what + " working points must be strictly increasing for channel " + channel + " in applywqqsf.yaml");
}

const std::map<std::string, WqqSFConfig> &applyWqqSFConfig() {
    static const auto config = [] {
        std::ifstream input(kConfigPath);
        if (!input) throw std::runtime_error("Cannot read Wqq SF configuration " + std::string(kConfigPath));
        std::map<std::string, WqqSFConfig> parsed;
        std::string channel, raw;
        while (std::getline(input, raw)) {
            auto line = raw;
            const auto comment = line.find('#');
            if (comment != std::string::npos) line.erase(comment);
            if (trim(line).empty()) continue;
            const bool indented = line.front() == ' ' || line.front() == '\t';
            line = trim(line);
            const auto separator = line.find(':');
            if (separator == std::string::npos) throw std::runtime_error("Invalid applywqqsf.yaml entry: " + raw);
            const auto key = trim(line.substr(0, separator));
            const auto value = trim(line.substr(separator + 1));
            if (!indented) {
                if (key.empty() || !value.empty()) throw std::runtime_error("Invalid applywqqsf.yaml channel line: " + raw);
                if (!parsed.emplace(key, WqqSFConfig{}).second)
                    throw std::runtime_error("Duplicate channel in applywqqsf.yaml: " + key);
                channel = key;
                continue;
            }
            if (channel.empty()) throw std::runtime_error("applywqqsf.yaml entry before any channel: " + raw);
            auto &entry = parsed.at(channel);
            if (key == "cut") {
                if (!entry.cut.empty()) throw std::runtime_error("Duplicate cut entry for channel " + channel + " in applywqqsf.yaml");
                entry.cut = parseNumberList(value, raw);
                auto sorted = entry.cut;
                std::sort(sorted.begin(), sorted.end());
                if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end())
                    throw std::runtime_error("Duplicate cut working point for channel " + channel + " in applywqqsf.yaml");
            } else if (key == "flags") {
                if (!entry.flags.empty()) throw std::runtime_error("Duplicate flags entry for channel " + channel + " in applywqqsf.yaml");
                entry.flags = parseNestedList(value, raw);
                for (const auto &set : entry.flags) validateWorkingPoints(set, "flags", channel);
                for (std::size_t i = 0; i < entry.flags.size(); ++i)
                    for (std::size_t j = i + 1; j < entry.flags.size(); ++j)
                        if (entry.flags[i] == entry.flags[j])
                            throw std::runtime_error("Duplicate flags set for channel " + channel + " in applywqqsf.yaml");
            } else {
                throw std::runtime_error("Unknown key '" + key + "' for channel " + channel + " in applywqqsf.yaml (expected cut or flags)");
            }
        }
        if (parsed.empty()) throw std::runtime_error("applywqqsf.yaml contains no channels");
        return parsed;
    }();
    return config;
}

const WqqSFConfig &wqqSFConfigForChannel(const std::string &channel) {
    const auto &config = applyWqqSFConfig();
    const auto entry = config.find(channel);
    if (entry == config.end()) throw std::runtime_error("Channel " + channel + " is missing from applywqqsf.yaml");
    return entry->second;
}

// 0.6 -> "0p60", 0.785 -> "0p785": at least two decimals, trailing zeros beyond that dropped.
std::string wpToken(double wp) {
    char buffer[32];
    std::snprintf(buffer, sizeof(buffer), "%.4f", wp);
    std::string text(buffer);
    while (text.size() > 4 && text.back() == '0' && text[text.size() - 3] != '.') text.pop_back();
    std::replace(text.begin(), text.end(), '.', 'p');
    return text;
}

std::map<std::string, WqqSFPayload> loadPayloads() {
    std::map<std::string, WqqSFPayload> out;
    for (const auto &year : kYears) {
        const std::string path = std::string(kPayloadDir) + "Wqq_SF_" + year + ".json";
        if (!std::filesystem::exists(path)) continue;
        const auto cset = correction::CorrectionSet::from_file(path);
        WqqSFPayload payload;
        payload.sf = cset->at("Wqq_SF");
        for (const auto &[name, corr] : *cset)
            if (name == "Wqq_eff_MC") payload.eff = corr;
        out.emplace(year, payload);
    }
    return out;
}

bool isTrueW(int gen_match) { return gen_match == 1 || gen_match == 2; }

} // namespace

RNode applyWqqScaleFactors(RNode df, const std::string &channel, bool isData) {
    const auto &config = wqqSFConfigForChannel(channel);
    if (config.empty()) return df;

    // One branch triple per cut WP, then one per flags set; output index = 3*triple + systematic.
    std::vector<std::string> stems;
    for (const double wp : config.cut) stems.push_back("FatJet_WqqSF_cut" + wpToken(wp));
    for (const auto &set : config.flags) {
        std::string stem = "FatJet_WqqSF_flags";
        for (std::size_t i = 0; i < set.size(); ++i) stem += (i ? "_" : "") + wpToken(set[i]);
        stems.push_back(stem);
    }

    if (isData) {
        for (const auto &stem : stems)
            for (const auto &syst : kSystematics)
                df = df.Define(stem + "_" + syst,
                               [](const RVec<float> &pt) { return RVec<float>(pt.size(), 1.f); }, {"FatJet_pt"});
        return df;
    }

    const auto payloads = std::make_shared<const std::map<std::string, WqqSFPayload>>(loadPayloads());
    if (payloads->empty()) throw std::runtime_error("No Wqq SF payloads found in " + std::string(kPayloadDir));
    for (const auto &[year, payload] : *payloads) {
        if (!config.flags.empty() && !payload.eff)
            throw std::runtime_error("Channel " + channel + " requests Wqq SF flags, but Wqq_SF_" + year +
                                     ".json has no Wqq_eff_MC (the candle MC efficiency the fail SF needs)");
        std::vector<double> requested = config.cut;
        for (const auto &set : config.flags) requested.insert(requested.end(), set.begin(), set.end());
        for (const double wp : requested) {
            try {
                payload.sf->evaluate({wp, kProbePt, 0., std::string("nom")});
                if (!config.flags.empty()) payload.eff->evaluate({wp, kProbePt});
            } catch (const std::exception &) {
                throw std::runtime_error("Wqq working point " + std::to_string(wp) + " requested by channel " + channel +
                                         " is outside the measured range of Wqq_SF_" + year + ".json");
            }
        }
    }

    const std::vector<double> cut_wps = config.cut;
    const std::vector<std::vector<double>> flag_sets = config.flags;
    auto evaluate = [payloads, cut_wps, flag_sets](const std::string &year, const RVec<float> &pt, const RVec<float> &eta,
                                                   const RVec<float> &score, const RVec<int> &gen_match) {
        const std::size_t n_triples = cut_wps.size() + flag_sets.size();
        std::vector<RVec<float>> out(3 * n_triples, RVec<float>(pt.size(), 1.f));
        const auto entry = payloads->find(year);
        for (std::size_t jet = 0; jet < pt.size(); ++jet) {
            if (!isTrueW(gen_match[jet])) continue;
            if (entry == payloads->end())
                throw std::runtime_error("No Wqq SF payload Wqq_SF_" + year + ".json for a true W jet");
            const auto &payload = entry->second;
            const double jet_pt = pt[jet], jet_eta = eta[jet];
            const auto sf = [&](double wp, const std::string &syst) {
                return payload.sf->evaluate({wp, jet_pt, jet_eta, syst});
            };

            std::size_t triple = 0;
            for (const double wp : cut_wps) {
                if (score[jet] > wp)
                    for (std::size_t s = 0; s < kSystematics.size(); ++s)
                        out[3 * triple + s][jet] = static_cast<float>(sf(wp, kSystematics[s]));
                ++triple;
            }
            for (const auto &set : flag_sets) {
                std::vector<double> eff;
                std::vector<bool> passed;
                for (const double wp : set) {
                    eff.push_back(payload.eff->evaluate({wp, jet_pt}));
                    passed.push_back(score[jet] > wp);
                }
                for (std::size_t s = 0; s < kSystematics.size(); ++s) {
                    std::vector<double> sfs;
                    for (const double wp : set) sfs.push_back(sf(wp, kSystematics[s]));
                    const auto result = cumulativeWPCategoryWeight(sfs, eff, passed, kDenominatorEpsilon);
                    if (result.status != CategoryWeightStatus::Ok) {
                        std::ostringstream message;
                        message << "Invalid Wqq SF category weight (status " << static_cast<int>(result.status)
                                << ") for year " << year << ", pt " << jet_pt << ", score " << score[jet]
                                << ", systematic " << kSystematics[s];
                        throw std::runtime_error(message.str());
                    }
                    out[3 * triple + s][jet] = static_cast<float>(result.value);
                }
                ++triple;
            }
        }
        return out;
    };

    df = df.Define("_wqq_sf", evaluate, {"year", "FatJet_pt", "FatJet_eta", "FatJet_Wqq", "FatJet_genMatch"});
    for (std::size_t triple = 0; triple < stems.size(); ++triple)
        for (std::size_t s = 0; s < kSystematics.size(); ++s)
            df = df.Define(stems[triple] + "_" + kSystematics[s],
                           "_wqq_sf[" + std::to_string(3 * triple + s) + "]");
    return df;
}
