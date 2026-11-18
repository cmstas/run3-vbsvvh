#ifndef CATEGORY_WEIGHT_H
#define CATEGORY_WEIGHT_H

#pragma once

#include <cmath>
#include <cstddef>
#include <vector>

// Per-object weight for a set of nested cumulative working points (loosest first),
// given the cumulative data/MC scale factors sf[i] and the MC efficiencies eff[i] of
// passing working point i, and which working points the object passes.
//
//   fails the loosest      (1 - sf[0]*eff[0]) / (1 - eff[0])
//   passes i, fails i+1    (sf[i]*eff[i] - sf[i+1]*eff[i+1]) / (eff[i] - eff[i+1])
//   passes the tightest    sf[N-1]
//
// With a single working point this is the usual pass SF / fail SF pair. Used by the
// b-tag weights and the Wqq tagging SFs.
enum class CategoryWeightStatus {
    Ok,
    InvalidProbability,
    NonnestedWorkingPoints,
    NegativeIntermediate,
    TinyDenominator,
};

struct CategoryWeightResult {
    double value;
    CategoryWeightStatus status;
};

inline CategoryWeightResult cumulativeWPCategoryWeight(const std::vector<double> &sf,
                                                       const std::vector<double> &eff,
                                                       const std::vector<bool> &passed,
                                                       double denominator_epsilon) {
    const auto failed = [](CategoryWeightStatus status) { return CategoryWeightResult{1., status}; };
    for (std::size_t index = 0; index < passed.size(); ++index) {
        if (!std::isfinite(sf[index]) || !std::isfinite(eff[index]))
            return failed(CategoryWeightStatus::InvalidProbability);
        if (index > 0 && passed[index] && !passed[index - 1])
            return failed(CategoryWeightStatus::NonnestedWorkingPoints);
        if (!(0. <= eff[index] && eff[index] <= 1.))
            return failed(CategoryWeightStatus::InvalidProbability);
    }
    for (std::size_t index = 1; index < eff.size(); ++index)
        if (eff[index - 1] < eff[index])
            return failed(CategoryWeightStatus::InvalidProbability);

    std::vector<double> q(sf.size());
    for (std::size_t index = 0; index < q.size(); ++index) {
        q[index] = sf[index] * eff[index];
        if (!std::isfinite(q[index]) || q[index] < 0. || q[index] > 1.)
            return failed(CategoryWeightStatus::InvalidProbability);
        if (index > 0 && q[index - 1] < q[index])
            return failed(CategoryWeightStatus::InvalidProbability);
    }

    int passed_index = -1;
    for (int index = static_cast<int>(passed.size()) - 1; index >= 0; --index)
        if (passed[static_cast<std::size_t>(index)]) { passed_index = index; break; }
    if (passed_index < 0) {
        if (std::abs(1. - eff[0]) < denominator_epsilon)
            return failed(CategoryWeightStatus::TinyDenominator);
        return {(1. - q[0]) / (1. - eff[0]), CategoryWeightStatus::Ok};
    }
    if (passed_index == static_cast<int>(passed.size()) - 1) return {sf.back(), CategoryWeightStatus::Ok};
    const std::size_t current = static_cast<std::size_t>(passed_index);
    const std::size_t next = current + 1;
    if (q[current] < q[next])
        return failed(CategoryWeightStatus::NegativeIntermediate);
    if (std::abs(eff[current] - eff[next]) < denominator_epsilon)
        return failed(CategoryWeightStatus::TinyDenominator);
    return {(q[current] - q[next]) / (eff[current] - eff[next]), CategoryWeightStatus::Ok};
}

#endif
