#include "btag_settings.h"

#include <fstream>
#include <algorithm>
#include <map>
#include <set>
#include <stdexcept>

namespace {

std::string trim(std::string value) {
    const auto begin = value.find_first_not_of(" \t");
    if (begin == std::string::npos) return "";
    const auto end = value.find_last_not_of(" \t");
    return value.substr(begin, end - begin + 1);
}

const std::map<std::string, std::vector<std::string>> &applyBTagConfig() {
    static const auto config = [] {
        constexpr const char *path = "applybtag.yaml";
        std::ifstream input(path);
        if (!input) throw std::runtime_error("Cannot read b-tag application configuration " + std::string(path));

        std::map<std::string, std::vector<std::string>> parsed;
        std::string line;
        while (std::getline(input, line)) {
            const auto comment = line.find('#');
            if (comment != std::string::npos) line.erase(comment);
            line = trim(line);
            if (line.empty()) continue;

            const auto separator = line.find(':');
            if (separator == std::string::npos || line.find(':', separator + 1) != std::string::npos)
                throw std::runtime_error("Invalid applybtag.yaml entry: " + line);
            const auto channel = trim(line.substr(0, separator));
            auto value = trim(line.substr(separator + 1));
            if (channel.empty() || value.size() < 2 || value.front() != '[' || value.back() != ']')
                throw std::runtime_error("Invalid applybtag.yaml entry: " + line);
            value = trim(value.substr(1, value.size() - 2));
            std::vector<std::string> working_points;
            if (!value.empty()) {
                std::size_t start = 0;
                while (start <= value.size()) {
                    const auto comma = value.find(',', start);
                    auto wp = trim(value.substr(start, comma == std::string::npos ? std::string::npos : comma - start));
                    if (std::find(kBTagInclusiveWorkingPoints.begin(), kBTagInclusiveWorkingPoints.end(), wp) == kBTagInclusiveWorkingPoints.end())
                        throw std::runtime_error("Unknown working point in applybtag.yaml entry: " + wp);
                    working_points.push_back(wp);
                    if (comma == std::string::npos) break;
                    start = comma + 1;
                }
            }
            if (std::adjacent_find(working_points.begin(), working_points.end()) != working_points.end())
                throw std::runtime_error("Duplicate working point in applybtag.yaml entry: " + line);
            for (std::size_t i = 1; i < working_points.size(); ++i)
                if (std::find(kBTagInclusiveWorkingPoints.begin(), kBTagInclusiveWorkingPoints.end(), working_points[i - 1]) >=
                    std::find(kBTagInclusiveWorkingPoints.begin(), kBTagInclusiveWorkingPoints.end(), working_points[i]))
                    throw std::runtime_error("Working points must be ordered L,M,T,XT,XXT: " + line);
            if (!parsed.emplace(channel, std::move(working_points)).second)
                throw std::runtime_error("Duplicate channel in applybtag.yaml: " + channel);
        }
        if (parsed.empty()) throw std::runtime_error("applybtag.yaml contains no channels");
        return parsed;
    }();
    return config;
}

std::map<std::string, BTagEraConfig> loadBTagTaggerConfig(const std::string &path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot read b-tagger configuration " + path);
    const auto fail = [&path](const std::string &message) {
        throw std::runtime_error("Invalid b-tagger configuration " + path + ": " + message);
    };

    // Two sections (taggers, eras), each mapping a name to `key: value` fields.
    std::map<std::string, std::map<std::string, std::map<std::string, std::string>>> sections;
    std::string line, section, name;
    while (std::getline(input, line)) {
        if (line.find('\t') != std::string::npos) fail("tabs are not supported");
        const auto comment = line.find('#');
        if (comment != std::string::npos) line.erase(comment);
        const auto content = trim(line);
        if (content.empty()) continue;
        const auto indent = line.find_first_not_of(' ');
        const auto separator = content.find(':');
        if (separator == std::string::npos) fail("missing ':' in " + content);
        const auto key = trim(content.substr(0, separator));
        const auto value = trim(content.substr(separator + 1));
        if (indent == 0 && value.empty() && (key == "taggers" || key == "eras")) {
            section = key;
            name.clear();
        } else if (indent == 2 && value.empty() && !section.empty() && !key.empty()) {
            name = key;
            if (!sections[section].emplace(name, std::map<std::string, std::string>{}).second)
                fail("duplicate " + section + " entry " + name);
        } else if (indent == 4 && !name.empty() && !value.empty()) {
            if (!sections[section][name].emplace(key, value).second) fail("duplicate key " + key + " in " + name);
        } else fail("invalid indentation or syntax: " + content);
    }

    const auto field = [&fail](const std::map<std::string, std::string> &fields, const std::string &key,
                               const std::string &owner) {
        const auto it = fields.find(key);
        if (it == fields.end()) fail(owner + " has no " + key);
        return it->second;
    };
    const auto check_keys = [&fail](const std::map<std::string, std::string> &fields,
                                    const std::set<std::string> &allowed, const std::string &owner) {
        for (const auto &[key, _] : fields)
            if (!allowed.count(key)) fail("unknown key " + key + " in " + owner);
    };
    const auto &taggers = sections["taggers"];
    std::map<std::string, BTagEraConfig> eras;
    for (const auto &[era, fields] : sections["eras"]) {
        check_keys(fields, {"tagger", "payload", "sf_max_abseta"}, era);
        BTagEraConfig config;
        config.tagger = field(fields, "tagger", era);
        const auto tagger = taggers.find(config.tagger);
        if (tagger == taggers.end()) fail(era + " uses unknown tagger " + config.tagger);
        check_keys(tagger->second, {"score_branch", "payload_prefix"}, config.tagger);
        config.score_branch = field(tagger->second, "score_branch", config.tagger);
        config.payload_prefix = field(tagger->second, "payload_prefix", config.tagger);
        config.payload = field(fields, "payload", era);
        try { config.sf_max_abseta = std::stod(field(fields, "sf_max_abseta", era)); }
        catch (const std::invalid_argument &) { fail("sf_max_abseta of " + era + " is not a number"); }
        if (!(config.sf_max_abseta > 0.)) fail("sf_max_abseta of " + era + " must be positive");
        eras.emplace(era, config);
    }
    if (eras.empty()) fail("no eras are configured");
    return eras;
}

} // namespace

const std::map<std::string, BTagEraConfig> &bTagEraConfigs() {
    static const auto config = loadBTagTaggerConfig("corrections/scalefactors/btagging/btag_taggers.yaml");
    return config;
}

const BTagEraConfig &bTagEraConfig(const std::string &year) {
    const auto &config = bTagEraConfigs();
    const auto entry = config.find(year);
    if (entry == config.end())
        throw std::runtime_error("No b-tagger is configured for year " + year + " in btag_taggers.yaml");
    return entry->second;
}

std::vector<std::string> bTagScoreBranches() {
    std::set<std::string> branches;
    for (const auto &[_, config] : bTagEraConfigs()) branches.insert(config.score_branch);
    return {branches.begin(), branches.end()};
}

const correction::CorrectionSet &bTagPayload(const std::string &year) {
    static const auto payloads = [] {
        std::map<std::string, correction::CorrectionSet> out;
        for (const auto &[era, config] : bTagEraConfigs())
            out.emplace(era, *correction::CorrectionSet::from_file(config.payload));
        return out;
    }();
    const auto entry = payloads.find(year);
    if (entry == payloads.end())
        throw std::runtime_error("No b-tag SF payload is configured for year " + year + " in btag_taggers.yaml");
    return entry->second;
}

std::vector<std::string> bTagWorkingPointsForChannel(const std::string &channel) {
    const auto &config = applyBTagConfig();
    const auto entry = config.find(channel);
    if (entry == config.end())
        throw std::runtime_error("Channel " + channel + " is missing from applybtag.yaml");
    return entry->second;
}

bool bTagScaleFactorsEnabled(const std::string &channel) {
    return !bTagWorkingPointsForChannel(channel).empty();
}
