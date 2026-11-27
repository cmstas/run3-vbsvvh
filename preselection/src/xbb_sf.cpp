#include "xbb_sf.h"
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

constexpr const char *kConfigPath = "applyxbbsf.yaml";
constexpr const char *kPayloadDir = "corrections/scalefactors/XbbSF/";
const std::array<std::string, 10> kYears = {"2016preVFP", "2016postVFP", "2017", "2018",
                                            "2022Re-recoBCD", "2022Re-recoE+PromptFG", "2023PromptC", "2023PromptD",
                                            "2024Prompt", "2025"};
const std::array<std::string, 3> kSystematics = {"nom", "up", "down"};
constexpr double kDenominatorEpsilon = 1.e-8;
// Probe pt for checking at graph construction that every requested WP is in range.
constexpr double kProbePt = 300.;

// The flags sets of one channel.
using XbbSFConfig = std::vector<std::vector<double>>;

struct XbbSFPayload {
    correction::Correction::Ref sf;
    correction::Correction::Ref eff;  // null when the file has no Xbb_eff_MC
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
        throw std::runtime_error("Invalid applyxbbsf.yaml list: " + line);
    std::vector<double> out;
    std::stringstream items(value.substr(1, value.size() - 2));
    std::string item;
    while (std::getline(items, item, ',')) {
        item = trim(item);
        if (item.empty()) continue;
        std::size_t used = 0;
        double wp = 0.;
        try { wp = std::stod(item, &used); } catch (const std::exception &) { used = 0; }
        if (used != item.size()) throw std::runtime_error("Invalid working point '" + item + "' in applyxbbsf.yaml: " + line);
        out.push_back(wp);
    }
    return out;
}

// "[[0.6], [0.4, 0.6, 0.8]]" -> {{0.6}, {0.4, 0.6, 0.8}}; "[]" -> {}.
std::vector<std::vector<double>> parseNestedList(const std::string &text, const std::string &line) {
    const auto value = trim(text);
    if (value.size() < 2 || value.front() != '[' || value.back() != ']')
        throw std::runtime_error("Invalid applyxbbsf.yaml list: " + line);
    std::vector<std::vector<double>> out;
    const auto inner = value.substr(1, value.size() - 2);
    std::size_t pos = 0;
    while (true) {
        const auto open = inner.find('[', pos);
        if (open == std::string::npos) {
            if (!trim(inner.substr(pos)).empty() && trim(inner.substr(pos)) != ",")
                throw std::runtime_error("Invalid applyxbbsf.yaml list: " + line);
            break;
        }
        const auto between = trim(inner.substr(pos, open - pos));
        if (!(between.empty() || (between == "," && !out.empty())))
            throw std::runtime_error("Invalid applyxbbsf.yaml list: " + line);
        const auto close = inner.find(']', open);
        if (close == std::string::npos) throw std::runtime_error("Invalid applyxbbsf.yaml list: " + line);
        out.push_back(parseNumberList(inner.substr(open, close - open + 1), line));
        pos = close + 1;
    }
    return out;
}

void validateWorkingPoints(const std::vector<double> &wps, const std::string &channel) {
    if (wps.empty()) throw std::runtime_error("Empty flags working-point set for channel " + channel + " in applyxbbsf.yaml");
    for (std::size_t i = 1; i < wps.size(); ++i)
        if (!(wps[i - 1] < wps[i]))
            throw std::runtime_error("flags working points must be strictly increasing for channel " + channel + " in applyxbbsf.yaml");
}

const std::map<std::string, XbbSFConfig> &applyXbbSFConfig() {
    static const auto config = [] {
        std::ifstream input(kConfigPath);
        if (!input) throw std::runtime_error("Cannot read Xbb SF configuration " + std::string(kConfigPath));
        std::map<std::string, XbbSFConfig> parsed;
        std::string channel, raw;
        while (std::getline(input, raw)) {
            auto line = raw;
            const auto comment = line.find('#');
            if (comment != std::string::npos) line.erase(comment);
            if (trim(line).empty()) continue;
            const bool indented = line.front() == ' ' || line.front() == '\t';
            line = trim(line);
            const auto separator = line.find(':');
            if (separator == std::string::npos) throw std::runtime_error("Invalid applyxbbsf.yaml entry: " + raw);
            const auto key = trim(line.substr(0, separator));
            const auto value = trim(line.substr(separator + 1));
            if (!indented) {
                if (key.empty() || !value.empty()) throw std::runtime_error("Invalid applyxbbsf.yaml channel line: " + raw);
                if (!parsed.emplace(key, XbbSFConfig{}).second)
                    throw std::runtime_error("Duplicate channel in applyxbbsf.yaml: " + key);
                channel = key;
                continue;
            }
            if (channel.empty()) throw std::runtime_error("applyxbbsf.yaml entry before any channel: " + raw);
            if (key != "flags")
                throw std::runtime_error("Unknown key '" + key + "' for channel " + channel + " in applyxbbsf.yaml (expected flags)");
            auto &sets = parsed.at(channel);
            if (!sets.empty()) throw std::runtime_error("Duplicate flags entry for channel " + channel + " in applyxbbsf.yaml");
            sets = parseNestedList(value, raw);
            for (const auto &set : sets) validateWorkingPoints(set, channel);
            for (std::size_t i = 0; i < sets.size(); ++i)
                for (std::size_t j = i + 1; j < sets.size(); ++j)
                    if (sets[i] == sets[j])
                        throw std::runtime_error("Duplicate flags set for channel " + channel + " in applyxbbsf.yaml");
        }
        if (parsed.empty()) throw std::runtime_error("applyxbbsf.yaml contains no channels");
        return parsed;
    }();
    return config;
}

const XbbSFConfig &xbbSFConfigForChannel(const std::string &channel) {
    const auto &config = applyXbbSFConfig();
    const auto entry = config.find(channel);
    if (entry == config.end()) throw std::runtime_error("Channel " + channel + " is missing from applyxbbsf.yaml");
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

std::map<std::string, XbbSFPayload> loadPayloads() {
    std::map<std::string, XbbSFPayload> out;
    for (const auto &year : kYears) {
        const std::string path = std::string(kPayloadDir) + "Xbb_SF_" + year + ".json";
        if (!std::filesystem::exists(path)) continue;
        const auto cset = correction::CorrectionSet::from_file(path);
        XbbSFPayload payload;
        payload.sf = cset->at("Xbb_SF");
        for (const auto &[name, corr] : *cset)
            if (name == "Xbb_eff_MC") payload.eff = corr;
        out.emplace(year, payload);
    }
    return out;
}

// True b jet: the jet contains both b quarks of an H->bb (FatJet_genMatch 5, see genSelections.h)
bool isTrueBB(int gen_match) { return gen_match == 5; }

} // namespace

RNode applyXbbScaleFactors(RNode df, const std::string &channel, bool isData) {
    const XbbSFConfig flag_sets = xbbSFConfigForChannel(channel);
    if (flag_sets.empty()) return df;

    // One branch triple per flags set; output index = 3*triple + systematic.
    std::vector<std::string> stems;
    for (const auto &set : flag_sets) {
        std::string stem = "FatJet_XbbSF_flags";
        for (std::size_t i = 0; i < set.size(); ++i) stem += (i ? "_" : "") + wpToken(set[i]);
        stems.push_back(stem);
    }

    // Checked before the data branch: one configuration is valid for data and MC alike.
    const auto payloads = std::make_shared<const std::map<std::string, XbbSFPayload>>(loadPayloads());
    if (payloads->empty()) throw std::runtime_error("No Xbb SF payloads found in " + std::string(kPayloadDir));
    for (const auto &[year, payload] : *payloads) {
        if (!payload.eff)
            throw std::runtime_error("Channel " + channel + " requests Xbb SF flags, but Xbb_SF_" + year +
                                     ".json has no Xbb_eff_MC (the candle MC efficiency the fail SF needs)");
        for (const auto &set : flag_sets) {
            for (const double wp : set) {
                try {
                    payload.sf->evaluate({wp, kProbePt, std::string("nom")});
                    payload.eff->evaluate({wp, kProbePt});
                } catch (const std::exception &) {
                    throw std::runtime_error("Xbb working point " + std::to_string(wp) + " requested by channel " + channel +
                                             " is outside the measured range of Xbb_SF_" + year + ".json");
                }
            }
        }
    }

    if (isData) {
        for (const auto &stem : stems)
            for (const auto &syst : kSystematics)
                df = df.Define(stem + "_" + syst,
                               [](const RVec<float> &pt) { return RVec<float>(pt.size(), 1.f); }, {"FatJet_pt"});
        return df;
    }

    auto evaluate = [payloads, flag_sets](const std::string &year, const RVec<float> &pt,
                                          const RVec<float> &score, const RVec<int> &gen_match) {
        std::vector<RVec<float>> out(3 * flag_sets.size(), RVec<float>(pt.size(), 1.f));
        const auto entry = payloads->find(year);
        for (std::size_t jet = 0; jet < pt.size(); ++jet) {
            if (!isTrueBB(gen_match[jet])) continue;
            if (entry == payloads->end())
                throw std::runtime_error("No Xbb SF payload Xbb_SF_" + year + ".json for a true b jet");
            const auto &payload = entry->second;
            const double jet_pt = pt[jet];

            std::size_t triple = 0;
            for (const auto &set : flag_sets) {
                std::vector<double> eff;
                std::vector<bool> passed;
                for (const double wp : set) {
                    eff.push_back(payload.eff->evaluate({wp, jet_pt}));
                    passed.push_back(score[jet] > wp);
                }
                for (std::size_t s = 0; s < kSystematics.size(); ++s) {
                    std::vector<double> sfs;
                    for (const double wp : set) sfs.push_back(payload.sf->evaluate({wp, jet_pt, kSystematics[s]}));
                    const auto result = cumulativeWPCategoryWeight(sfs, eff, passed, kDenominatorEpsilon);
                    if (result.status != CategoryWeightStatus::Ok) {
                        std::ostringstream message;
                        message << "Invalid Xbb SF category weight (status " << static_cast<int>(result.status)
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

    df = df.Define("_xbb_sf", evaluate, {"year", "FatJet_pt", "FatJet_Hbb", "FatJet_genMatch"});
    for (std::size_t triple = 0; triple < stems.size(); ++triple)
        for (std::size_t s = 0; s < kSystematics.size(); ++s)
            df = df.Define(stems[triple] + "_" + kSystematics[s],
                           "_xbb_sf[" + std::to_string(3 * triple + s) + "]");
    return df;
}
