#ifndef BTAG_SETTINGS_H
#define BTAG_SETTINGS_H

#pragma once

#include <string>

inline double bTagMaxAbsEta(const std::string &year) {
    // UParTAK4_comb fixed-WP payloads in Run-2 NanoAODv15 have abseta
    // edges [0, 2.4]; the Summer24 payload extends to 2.5.
    return (year == "2016preVFP" || year == "2016postVFP" ||
            year == "2017" || year == "2018") ? 2.4 : 2.5;
}

inline std::string bTagSafeYearToken(std::string year) {
    if (year == "2022Re-recoBCD") return "2022Re_recoBCD";
    if (year == "2022Re-recoE+PromptFG") return "2022Re_recoE_PromptFG";
    return year;
}

#endif
