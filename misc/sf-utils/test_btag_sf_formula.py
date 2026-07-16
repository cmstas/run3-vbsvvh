#!/usr/bin/env python3
"""Small regression test for the two-nested-WP event reweighting formula.

Run with ``python3 misc/sf-utils/test_btag_sf_formula.py``.  The same three
branches are implemented in preselection/src/weights.cpp.
"""

import math


EPS = 1e-8


def jet_weight(sf_tight, sf_loose, eff_tight, eff_loose, is_tight, is_loose):
    if not (0.0 <= eff_tight <= eff_loose <= 1.0):
        return 1.0
    q_tight = sf_tight * eff_tight
    q_loose = sf_loose * eff_loose
    if not (0.0 <= q_tight <= q_loose <= 1.0):
        return 1.0
    if is_tight:
        return sf_tight
    if is_loose:
        numerator = q_loose - q_tight
        denominator = eff_loose - eff_tight
        if numerator < 0.0 or abs(denominator) < EPS:
            return 1.0
        return numerator / denominator
    numerator = 1.0 - q_loose
    denominator = 1.0 - eff_loose
    if numerator < 0.0 or abs(denominator) < EPS:
        return 1.0
    return numerator / denominator


def close(actual, expected):
    assert math.isclose(actual, expected, rel_tol=0.0, abs_tol=1e-12), (actual, expected)


def main():
    # Unit SFs preserve every observed category.
    close(jet_weight(1.0, 1.0, 0.4, 0.7, True, True), 1.0)
    close(jet_weight(1.0, 1.0, 0.4, 0.7, False, True), 1.0)
    close(jet_weight(1.0, 1.0, 0.4, 0.7, False, False), 1.0)

    sf_tight, sf_loose, eff_tight, eff_loose = 0.9, 1.1, 0.4, 0.7
    close(jet_weight(sf_tight, sf_loose, eff_tight, eff_loose, True, True), sf_tight)
    close(jet_weight(sf_tight, sf_loose, eff_tight, eff_loose, False, True),
          (sf_loose * eff_loose - sf_tight * eff_tight) / (eff_loose - eff_tight))
    close(jet_weight(sf_tight, sf_loose, eff_tight, eff_loose, False, False),
          (1.0 - sf_loose * eff_loose) / (1.0 - eff_loose))

    # Required safety fallbacks.
    close(jet_weight(1.1, 0.5, 0.4, 0.7, False, True), 1.0)  # negative LT data probability
    close(jet_weight(0.9, 1.1, 0.7, 0.7, False, True), 1.0)  # tiny LT denominator
    close(jet_weight(0.9, 1.1, 0.4, 1.0, False, False), 1.0)  # tiny N denominator
    close(jet_weight(0.9, 1.1, 0.8, 0.7, False, True), 1.0)  # invalid nesting

    # Invalid scaled probabilities must be rejected for every observed category.
    for is_tight, is_loose in ((True, True), (False, True), (False, False)):
        close(jet_weight(1.3, 1.2, 0.8, 0.9, is_tight, is_loose), 1.0)

    # Event factors are products of observed per-jet categories.
    event = (jet_weight(sf_tight, sf_loose, eff_tight, eff_loose, True, True) *
             jet_weight(sf_tight, sf_loose, eff_tight, eff_loose, False, False))
    close(event, sf_tight * (1.0 - sf_loose * eff_loose) / (1.0 - eff_loose))
    print("two-WP b-tag SF formula checks passed")


if __name__ == "__main__":
    main()
