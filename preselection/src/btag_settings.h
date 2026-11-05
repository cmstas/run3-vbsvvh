#ifndef BTAG_SETTINGS_H
#define BTAG_SETTINGS_H

#pragma once

#include <array>
#include <string>
#include <string_view>
#include <cmath>
#include <algorithm>
#include <map>
#include <vector>

#include "correction.h"

// Canonical b-tag working-point order, from loosest to tightest.  Keep this
// definition shared by efficiency production and runtime SF application.
inline constexpr std::array<std::string_view, 5> kBTagInclusiveWorkingPoints = {
    "L", "M", "T", "XT", "XXT"
};

// Mutually exclusive observed categories corresponding to the ordered WPs.
inline constexpr std::array<std::string_view, 6> kBTagExclusiveCategories = {
    "N", "LnotM", "MnotT", "TnotXT", "XTnotXXT", "XXT"
};

// AK4 b-tagger of an era, from corrections/scalefactors/btagging/btag_taggers.yaml.
struct BTagEraConfig {
    std::string tagger;
    std::string score_branch;   // NanoAOD discriminant, e.g. Jet_btagUParTAK4B
    std::string payload_prefix; // BTV corrections <prefix>_comb, _light and _wp_values
    std::string payload;        // BTV btagging.json.gz
    double sf_max_abseta;
};

const std::map<std::string, BTagEraConfig> &bTagEraConfigs();
const BTagEraConfig &bTagEraConfig(const std::string &year);
// Distinct discriminant branches of the configured taggers.
std::vector<std::string> bTagScoreBranches();
const correction::CorrectionSet &bTagPayload(const std::string &year);

inline double bTagMaxAbsEta(const std::string &year) {
    return (year == "2016preVFP" || year == "2016postVFP") ? 2.4 : 2.5;
}

// Some SF payloads (e.g. Run-2 UParTAK4) are tabulated below |eta|=2.4, while
// the selected-jet acceptance remains |eta|<2.5 outside 2016. Clamp only the SF
// lookup coordinate to the valid payload domain; do not discard those jets.
inline double bTagSFAbsEta(const std::string &year, double eta) {
    return std::min(std::abs(eta), std::nextafter(bTagEraConfig(year).sf_max_abseta, 0.0));
}

inline std::string bTagSafeYearToken(std::string year) {
    if (year == "2022Re-recoBCD") return "2022Re_recoBCD";
    if (year == "2022Re-recoE+PromptFG") return "2022Re_recoE_PromptFG";
    return year;
}

std::vector<std::string> bTagWorkingPointsForChannel(const std::string &channel);
bool bTagScaleFactorsEnabled(const std::string &channel);

#endif
