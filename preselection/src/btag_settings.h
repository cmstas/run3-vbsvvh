#ifndef BTAG_SETTINGS_H
#define BTAG_SETTINGS_H

#pragma once

#include <string>

inline double bTagMaxAbsEta(const std::string &year) {
    return (year == "2016preVFP" || year == "2016postVFP") ? 2.4 : 2.5;
}

inline std::string bTagSafeYearToken(std::string year) {
    if (year == "2022Re-recoBCD") return "2022Re_recoBCD";
    if (year == "2022Re-recoE+PromptFG") return "2022Re_recoE_PromptFG";
    return year;
}

#endif
