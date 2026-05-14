#include "corrections.h"

/*
############################################
B-TAGGING WORKING POINTS
############################################
*/

RVec<bool> isbTagLoose(std::string year, RVec<float> btag_score) {
    return btag_score > btaggingWPMap_Loose.at(year);
}

RVec<bool> isbTagMedium(std::string year, RVec<float> btag_score) {
    return btag_score > btaggingWPMap_Medium.at(year);
}

RVec<bool> isbTagTight(std::string year, RVec<float> btag_score) {
    return btag_score > btaggingWPMap_Tight.at(year);
}

/*
############################################
JET MASS SCALE (JMS) AND RESOLUTION (JMR)
############################################

Parametrization is applied on `FatJet_msoftdrop`:

  JMS  : msd' = max(0, msd + Δshift[year])               // additive shift in GeV
  JMR  : matched   →  msd' = max(0, m_gen + f * (msd - m_gen))   
         unmatched →  msd' = max(0, msd * (1 + √max(f²-1,0) * σrel * N(0,1)))
                      // σrel = 1.0 is a placeholder; replace with the per-era msd
                      // resolution when JMR is derived. Branch is dead while f == 1.0.

Current factors are placeholders (Δshift = 0 GeV and f = 1.0), both reduce to the identity.
Replace with the calibration outputs once derived. 
TO DO: add systematic up/down, add a `variation` argument and pass ±1σ shift/scale instead of the central.

*/

RNode applyJetMassScale(const std::unordered_map<std::string, float>& jms_shift, RNode df) {
    auto eval = [jms_shift](std::string year, RVec<float> msd) {
        auto it = jms_shift.find(year);
        if (it == jms_shift.end()) {
            static std::unordered_set<std::string> warned_years;
            if (warned_years.insert(year).second) {
                std::cout << "Warning: JMS shift not defined for year " << year
                          << ". Skipping." << std::endl;
            }
            return msd;
        }
        const float shift = it->second;
        if (shift == 0.0f || msd.empty()) return msd;
        // apply jet mass cale shift to every fat jet's soft-drop mass in the event
        for (auto& m : msd) m = std::max(0.0f, m + shift);
        return msd;
    };
    return df.Redefine("FatJet_msoftdrop", eval, {"year", "FatJet_msoftdrop"});
}

RNode applyJetMassResolution(const std::unordered_map<std::string, float>& jmr_factor,
                             const std::unordered_map<std::string, float>& jmr_sigma_rel,
                             RNode df) {
    auto eval = [jmr_factor, jmr_sigma_rel](std::string year,
                             RVec<float> msd,
                             RVec<short> genJet_idx, RVec<float> genJet_mass,
                             unsigned int lumi, unsigned long long event) {
        auto it = jmr_factor.find(year);
        if (it == jmr_factor.end()) {
            static std::unordered_set<std::string> warned_years;
            if (warned_years.insert(year).second) {
                std::cout << "Warning: JMR factor not defined for year " << year
                          << ". Skipping." << std::endl;
            }
            return msd;
        }
        const float f = it->second;
        if (f == 1.0f || msd.empty()) return msd;

        const float widen = std::sqrt(std::max(f * f - 1.0f, 0.0f));
        // sigma_rel = relative msoftdrop resolution σ(msd)/msd, used by the unmatched
        // stochastic branch. Per-era placeholder maps jetMassResolution_sigmaRel_central
        // in corrections.h; replace with the value from the same calibration fit as
        // jmr_factor. Falls back to 1.0 if the year is missing from the map.
        auto sr_it = jmr_sigma_rel.find(year);
        const float sigma_rel = (sr_it != jmr_sigma_rel.end()) ? sr_it->second : 1.0f;
        TRandom3 rng((static_cast<unsigned long long>(lumi) << 10) ^ event);
        for (size_t i = 0; i < msd.size(); ++i) {
            const bool matched = (genJet_idx[i] >= 0
                                  && genJet_idx[i] < static_cast<int>(genJet_mass.size()));
            if (matched) {
                const float m_gen = genJet_mass[genJet_idx[i]];
                msd[i] = std::max(0.0f, m_gen + f * (msd[i] - m_gen));
            } else {
                const float kick = static_cast<float>(rng.Gaus(0.0, 1.0));
                msd[i] = std::max(0.0f, msd[i] * (1.0f + widen * sigma_rel * kick));
            }
        }
        return msd;
    };
    return df.Redefine("FatJet_msoftdrop", eval,
                       {"year", "FatJet_msoftdrop",
                        "FatJet_genJetAK8Idx", "GenJetAK8_mass",
                        "luminosityBlock", "event"});
}

/*
############################################
HEM Corrections
############################################
*/

RNode HEMCorrection(RNode df, bool isData) {
    auto HEMCorrections = [isData](unsigned int run, unsigned long long event, std::string sample_year, RVec<float> pt, RVec<float> eta, RVec<float> phi, RVec<float> jet_id) {
        RVec<bool> jet_mask;   
        if (sample_year == "2018" && ((isData && run >= 319077) || (!isData && event % 100 < 64))) {
            jet_mask = (jet_id >= 2 && pt > 15.0); // NanoAOD jetID convention https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD#Jets
                                                    // should still work for current skimmer, which sets jetId==3 : "pass tight ID, fail tightLepVeto", jetId==7 : "pass tight and tightLepVeto ID"
            auto eta_ = eta[jet_mask];
            auto phi_ = phi[jet_mask];
            for (size_t i = 0; i < eta_.size(); i++) {
                if (eta_[i] > -3.2 && eta_[i] < -1.3 && phi_[i] > -1.57 && phi_[i] < -0.87) {
                    return false;
                }
            }
        }
        return true;
    };

    return df.Define("HEMVeto", HEMCorrections, {"run", "event", "year", "Jet_pt", "Jet_eta", "Jet_phi", "Jet_jetId"}).Filter("HEMVeto");
}

/*
############################################
MET PHI CORRECTIONS
############################################
*/

RNode applyMETPhiCorrections(RNode df, bool isData) {
    auto eval_correction = [isData] (std::string year, float pt, float phi, unsigned char npvs, unsigned int run) {
        double pt_corr = pt;
        double phi_corr = phi;
        
        if (metCorrections.find(year) == metCorrections.end()) {
            static std::unordered_set<std::string> warned_years;
            if (warned_years.find(year) == warned_years.end()) {
                std::cout << "Warning: MET correction set for year " << year << " not found. Skipping MET phi corrections." << std::endl;
                warned_years.insert(year);
            }
            return std::make_pair(static_cast<float>(pt_corr), static_cast<float>(phi_corr));
        }

        std::string pt_corr_name = isData ? "pt_metphicorr_puppimet_data" : "pt_metphicorr_puppimet_mc";
        std::string phi_corr_name = isData ? "phi_metphicorr_puppimet_data" : "phi_metphicorr_puppimet_mc";

        // Note the correction json is only defined for up to 6500, so do not pass larger values
        float pt_max = 6499.9;
        float pt_to_pass = std::min(pt,pt_max);
        pt_corr = metCorrections.at(year).at(pt_corr_name)->evaluate({pt_to_pass, phi, static_cast<double>(npvs), static_cast<double>(run)});
        phi_corr = metCorrections.at(year).at(phi_corr_name)->evaluate({pt_to_pass, phi, static_cast<double>(npvs), static_cast<double>(run)});
        
        return std::make_pair(static_cast<float>(pt_corr), static_cast<float>(phi_corr));
    };
    return df.Define("_MET_phicorr", eval_correction, {"year", "PuppiMET_pt", "PuppiMET_phi", "PV_npvs", "run"})
            .Define("met_pt", "_MET_phicorr.first")            
            .Define("met_phi", "_MET_phicorr.second");
}

/*
############################################
MET UNCLUSTERED CORRECTIONS
############################################
*/

RNode applyMETUnclusteredCorrections(RNode df, std::string variation) {
    if (variation == "up") {
        return df.Define("_MET_uncert_dx", "met_pt * TMath::Cos(met_phi) + MET_MetUnclustEnUpDeltaX")
                .Define("_MET_uncert_dy", "met_pt * TMath::Sin(met_phi) + MET_MetUnclustEnUpDeltaY")
                .Redefine("met_pt", "(float)TMath::Sqrt(_MET_uncert_dx * _MET_uncert_dx + _MET_uncert_dy * _MET_uncert_dy)");
    }

    else if (variation == "down") {
        return df.Define("_MET_uncert_dx", "met_pt * TMath::Cos(met_phi) - MET_MetUnclustEnUpDeltaX")
                .Define("_MET_uncert_dy", "met_pt * TMath::Sin(met_phi) - MET_MetUnclustEnUpDeltaY")
                .Redefine("met_pt", "(float)TMath::Sqrt(_MET_uncert_dx * _MET_uncert_dx + _MET_uncert_dy * _MET_uncert_dy)");
    }
    return df;
}

/*
############################################
JET ENERGY CORRECTIONS — nominal 
############################################

Recipe (per AK4 jet, year-aware):
  pt_raw   = (1 - Jet_rawFactor) * Jet_pt     // undo NanoAOD JEC
  mass_raw = (1 - Jet_rawFactor) * Jet_mass

  factor   = compound(Jet_area, Jet_eta, pt_raw, rho [, run])
             // compound = L1FastJet * L2Relative * L3Absolute (* L2L3Residual for DATA),
             // with `inputs_update=[JetPt]` so each stage feeds the next its updated pt.
  pt_new   = pt_raw   * factor
  mass_new = mass_raw * factor

Type-I MET is rebuilt on top of the NanoAOD baseline:
  met_x_new = met_x_old + sum_jets (pt_NanoAOD - pt_new) * cos(phi)
  met_y_new = met_y_old + sum_jets (pt_NanoAOD - pt_new) * sin(phi)
where the sum runs over jets with pt_new > 15 GeV and (neEmEF + chEmEF) < 0.9.
FIXME: this needs to be updated to the new recommendation: https://indico.cern.ch/event/1644923/contributions/6916115/attachments/3211593/5720863/260202_JMEGeneral_Type1METWithNano_Nurfikri.pdf

NanoAOD branches consumed:
  Jet_pt, Jet_mass, Jet_eta, Jet_phi, Jet_area, Jet_rawFactor, Jet_neEmEF, Jet_chEmEF,
  Rho_fixedGridRhoFastjetAll, run, year
*/

namespace {
    inline std::string makeCompoundJECName(const std::string& mc_prefix,
                                           const std::string& algo,
                                           bool isData) {
        // jec_prefix_map values end in "_MC" by convention; replace with "_DATA" for data.
        std::string base = mc_prefix;
        const std::string mc_tag = "_MC";
        // check base ends with "_MC" and chop those last 3 characters off
        if (base.size() >= mc_tag.size() &&
            base.compare(base.size() - mc_tag.size(), mc_tag.size(), mc_tag) == 0) {
            base.resize(base.size() - mc_tag.size());
        }
        return base + (isData ? "_DATA" : "_MC") + "_L1L2L3Res_" + algo;
    }
}

RNode applyJetEnergyCorrections(const std::unordered_map<std::string, correction::CorrectionSet>& cset_jerc,
                                const std::unordered_map<std::string, std::string>& jec_prefix_map,
                                const std::unordered_map<std::string, std::string>& jec_suffix_map,
                                RNode df, bool isData) {

    auto compute_factor = [cset_jerc, jec_prefix_map, jec_suffix_map, isData](
            std::string year,
            RVec<float> pt, RVec<float> eta, RVec<float> area, RVec<float> rawFactor,
            float rho, unsigned int run) {

        RVec<float> factor(pt.size(), 1.0f);
        if (pt.empty()) return factor;

        auto cs_it = cset_jerc.find(year);
        auto pf_it = jec_prefix_map.find(year);
        auto sf_it = jec_suffix_map.find(year);
        if (cs_it == cset_jerc.end() || pf_it == jec_prefix_map.end() || sf_it == jec_suffix_map.end()) {
            static std::unordered_set<std::string> warned_years;
            if (warned_years.insert(year).second) {
                std::cout << "Warning: JEC inputs missing for year " << year
                          << ". Skipping nominal JEC re-application." << std::endl;
            }
            return factor;
        }

        const std::string compound_name = makeCompoundJECName(pf_it->second, sf_it->second, isData);
        std::shared_ptr<const correction::CompoundCorrection> compound;
        try {
            compound = cs_it->second.compound().at(compound_name);
        } catch (const std::out_of_range&) {
            static std::unordered_set<std::string> warned_keys;
            if (warned_keys.insert(year + "/" + compound_name).second) {
                std::cout << "Warning: compound JEC '" << compound_name
                          << "' not found for year " << year
                          << ". Skipping nominal JEC re-application." << std::endl;
            }
            return factor;
        }

        for (size_t i = 0; i < pt.size(); ++i) {
            float pt_raw = (1.0f - rawFactor[i]) * pt[i];
            double sf = 1.0;
            try {
                if (isData) {
                    sf = compound->evaluate({(double)area[i], (double)eta[i],
                                             (double)pt_raw,  (double)rho,
                                             (double)run});
                } else {
                    sf = compound->evaluate({(double)area[i], (double)eta[i],
                                             (double)pt_raw,  (double)rho});
                }
            } catch (const std::exception& e) {
                sf = 1.0;
            }
            // Multiplicative factor on the *NanoAOD* pt: Jet_pt_new = Jet_pt * factor
            //   = Jet_pt * (pt_raw / Jet_pt) * sf = (1 - rawFactor) * sf
            factor[i] = (1.0f - rawFactor[i]) * static_cast<float>(sf);
        }
        return factor;
    };

    auto met_propagate = [](RVec<float> pt_old, RVec<float> phi,
                            RVec<float> factor, RVec<float> neEmEF, RVec<float> chEmEF,
                            float met_pt, float met_phi) {
        // Type-I MET propagation: replace the contribution of NanoAOD-corrected jets
        // (pt_old) with re-corrected jets (pt_new = pt_old * factor) in the MET sum.
        double px = met_pt * std::cos(met_phi);
        double py = met_pt * std::sin(met_phi);
        for (size_t i = 0; i < pt_old.size(); ++i) {
            const float pt_new = pt_old[i] * factor[i];
            if (pt_new <= 15.0f) continue;                       // Type-I threshold
            if ((neEmEF[i] + chEmEF[i]) > 0.9f) continue;        // EM-fraction veto
            const double dpt = static_cast<double>(pt_old[i] - pt_new);
            px += dpt * std::cos(phi[i]);
            py += dpt * std::sin(phi[i]);
        }
        return std::make_pair(static_cast<float>(std::sqrt(px*px + py*py)),
                              static_cast<float>(std::atan2(py, px)));
    };

    return df
        .Define("Jet_jecFactor", compute_factor,
                {"year", "Jet_pt", "Jet_eta", "Jet_area", "Jet_rawFactor",
                 "Rho_fixedGridRhoFastjetAll", "run"})
        .Define("_jec_metProp", met_propagate,
                {"Jet_pt", "Jet_phi", "Jet_jecFactor",
                 "Jet_neEmEF", "Jet_chEmEF", "met_pt", "met_phi"})
        .Redefine("met_pt",   "_jec_metProp.first")
        .Redefine("met_phi",  "_jec_metProp.second")
        .Redefine("Jet_pt",   "Jet_pt   * Jet_jecFactor")
        .Redefine("Jet_mass", "Jet_mass * Jet_jecFactor")
        // Reset rawFactor so downstream callers see a self-consistent (pt, mass, rawFactor) triple.
        // pt_new * (1 - rawFactor_new) == pt_raw  =>  rawFactor_new = 1 - (1-rawFactor_old)/Jet_jecFactor
        .Redefine("Jet_rawFactor",
                  [](RVec<float> rawFactor, RVec<float> factor) {
                      RVec<float> out(rawFactor.size());
                      for (size_t i = 0; i < rawFactor.size(); ++i) {
                          out[i] = (factor[i] > 0.f)
                                       ? 1.0f - (1.0f - rawFactor[i]) / factor[i]
                                       : rawFactor[i];
                      }
                      return out;
                  },
                  {"Jet_rawFactor", "Jet_jecFactor"});
}

/*
############################################
FAT JET (AK8) ENERGY CORRECTIONS — nominal 
############################################

Same recipe as AK4 (raw recovery + L1L2L3Res compound from fatJet_jerc.json.gz, with the
AK8PFPuppi algo). Operates on FatJet_* branches. Does NOT propagate to MET — PuppiMET is
built from AK4 jets only.
*/

RNode applyFatJetEnergyCorrections(const std::unordered_map<std::string, correction::CorrectionSet>& cset_jerc,
                                   const std::unordered_map<std::string, std::string>& jec_prefix_map,
                                   const std::unordered_map<std::string, std::string>& jec_suffix_map,
                                   RNode df, bool isData) {

    auto compute_factor = [cset_jerc, jec_prefix_map, jec_suffix_map, isData](
            std::string year,
            RVec<float> pt, RVec<float> eta, RVec<float> area, RVec<float> rawFactor,
            float rho, unsigned int run) {

        RVec<float> factor(pt.size(), 1.0f);
        if (pt.empty()) return factor;

        auto cs_it = cset_jerc.find(year);
        auto pf_it = jec_prefix_map.find(year);
        auto sf_it = jec_suffix_map.find(year);
        if (cs_it == cset_jerc.end() || pf_it == jec_prefix_map.end() || sf_it == jec_suffix_map.end()) {
            static std::unordered_set<std::string> warned_years;
            if (warned_years.insert(year).second) {
                std::cout << "Warning: AK8 JEC inputs missing for year " << year
                          << ". Skipping nominal AK8 JEC re-application." << std::endl;
            }
            return factor;
        }

        const std::string compound_name = makeCompoundJECName(pf_it->second, sf_it->second, isData);
        std::shared_ptr<const correction::CompoundCorrection> compound;
        try {
            compound = cs_it->second.compound().at(compound_name);
        } catch (const std::out_of_range&) {
            static std::unordered_set<std::string> warned_keys;
            if (warned_keys.insert(year + "/" + compound_name).second) {
                std::cout << "Warning: AK8 compound JEC '" << compound_name
                          << "' not found for year " << year
                          << ". Skipping nominal AK8 JEC re-application." << std::endl;
            }
            return factor;
        }

        for (size_t i = 0; i < pt.size(); ++i) {
            float pt_raw = (1.0f - rawFactor[i]) * pt[i];
            double sf = 1.0;
            try {
                if (isData) {
                    sf = compound->evaluate({(double)area[i], (double)eta[i],
                                             (double)pt_raw,  (double)rho,
                                             (double)run});
                } else {
                    sf = compound->evaluate({(double)area[i], (double)eta[i],
                                             (double)pt_raw,  (double)rho});
                }
            } catch (const std::exception&) {
                sf = 1.0;
            }
            factor[i] = (1.0f - rawFactor[i]) * static_cast<float>(sf);
        }
        return factor;
    };

    return df
        .Define("FatJet_jecFactor", compute_factor,
                {"year", "FatJet_pt", "FatJet_eta", "FatJet_area", "FatJet_rawFactor",
                 "Rho_fixedGridRhoFastjetAll", "run"})
        .Redefine("FatJet_pt",   "FatJet_pt   * FatJet_jecFactor")
        .Redefine("FatJet_mass", "FatJet_mass * FatJet_jecFactor")
        // Keep FatJet_rawFactor self-consistent against the new pt/mass.
        .Redefine("FatJet_rawFactor",
                  [](RVec<float> rawFactor, RVec<float> factor) {
                      RVec<float> out(rawFactor.size());
                      for (size_t i = 0; i < rawFactor.size(); ++i) {
                          out[i] = (factor[i] > 0.f)
                                       ? 1.0f - (1.0f - rawFactor[i]) / factor[i]
                                       : rawFactor[i];
                      }
                      return out;
                  },
                  {"FatJet_rawFactor", "FatJet_jecFactor"});
}

/*
############################################
JET ENERGY SCALE — per-source uncertainties (up/down) into suffixed columns
############################################

Per-source JES uncertainty applied as a multiplicative shift on top of the nominal JEC.
Run *after* applyJetEnergyCorrections / applyFatJetEnergyCorrections. Each call writes
one set of suffixed columns; the driver applyJESVariations loops over the 11 Regrouped V2
sources × {Up, Dn} = 22 variations.

  AK4: Jet_pt_<sfx>, Jet_mass_<sfx>, met_pt_<sfx>, met_phi_<sfx>   (Type-I MET propagated)
  AK8: FatJet_pt_<sfx>, FatJet_mass_<sfx>                          (no MET propagation;
                                                                    FatJet_msoftdrop has
                                                                    its own JEC recipe)
  <sfx> = "jes" + column_label + direction   e.g. "jesAbsoluteUp", "jesAbsoluteYearDn"

Source list lives in kJESRegroupedSources below. Lookup keys in jet_jerc.json.gz are
  <prefix>_<base_source>[_<year_token>]_<suffix>
where year_token comes from jetEnergyCorrections_yearToken when the source is year-
decorrelated (e.g. Regrouped_Absolute_2018 vs the year-correlated Regrouped_Absolute).
*/

namespace {
// Returns a closure giving per-jet shift factors (1 + sign*u) for use with RDF::Define.
auto makeJESShiftEvaluator(
        const std::unordered_map<std::string, correction::CorrectionSet>& cset_jerc,
        const std::unordered_map<std::string, std::string>& jec_prefix_map,
        const std::unordered_map<std::string, std::string>& jec_suffix_map,
        const std::unordered_map<std::string, std::string>& year_token_map,
        std::string base_source, bool yearDecorrelated, float sign) {
    return [cset_jerc, jec_prefix_map, jec_suffix_map, year_token_map,
            base_source, yearDecorrelated, sign](
            std::string year, RVec<float> pt, RVec<float> eta) {
        RVec<float> out(pt.size(), 1.0f);
        if (pt.empty()) return out;
        auto cs_it = cset_jerc.find(year);
        auto pf_it = jec_prefix_map.find(year);
        auto sf_it = jec_suffix_map.find(year);
        if (cs_it == cset_jerc.end() || pf_it == jec_prefix_map.end() || sf_it == jec_suffix_map.end()) {
            return out;
        }
        std::string src = base_source;
        if (yearDecorrelated) {
            auto yt_it = year_token_map.find(year);
            if (yt_it == year_token_map.end()) return out;
            src += "_" + yt_it->second;
        }
        const std::string unc_name = pf_it->second + "_" + src + "_" + sf_it->second;
        std::shared_ptr<const correction::Correction> unc;
        try {
            unc = cs_it->second.at(unc_name);
        } catch (const std::out_of_range&) {
            static std::unordered_set<std::string> warned;
            if (warned.insert(year + "/" + unc_name).second) {
                std::cout << "Warning: JEC uncertainty source '" << unc_name
                          << "' not found for year " << year << "." << std::endl;
            }
            return out;
        }
        for (size_t i = 0; i < pt.size(); ++i) {
            float u = unc->evaluate({(double)eta[i], (double)pt[i]});
            out[i] = 1.0f + sign * u;
        }
        return out;
    };
}

// Type-I-style MET shift: rebuilds (met_pt, met_phi) from (pt_old - pt_new) per AK4 jet,
// applying the standard >15 GeV / emf<0.9 filter.
auto jesMetPropagator() {
    return [](RVec<float> pt_old, RVec<float> pt_new, RVec<float> phi,
              RVec<float> neEmEF, RVec<float> chEmEF,
              float met_pt, float met_phi) {
        double px = met_pt * std::cos(met_phi);
        double py = met_pt * std::sin(met_phi);
        for (size_t i = 0; i < pt_old.size(); ++i) {
            if (pt_new[i] <= 15.0f) continue;
            if ((neEmEF[i] + chEmEF[i]) > 0.9f) continue;
            const double dpt = static_cast<double>(pt_old[i] - pt_new[i]);
            px += dpt * std::cos(phi[i]);
            py += dpt * std::sin(phi[i]);
        }
        return std::make_pair(static_cast<float>(std::sqrt(px*px + py*py)),
                              static_cast<float>(std::atan2(py, px)));
    };
}
} // anonymous namespace

RNode defineAK4JESVariation(
        const std::unordered_map<std::string, correction::CorrectionSet>& cset_jerc,
        const std::unordered_map<std::string, std::string>& jec_prefix_map,
        const std::unordered_map<std::string, std::string>& jec_suffix_map,
        const std::unordered_map<std::string, std::string>& year_token_map,
        RNode df, const std::string& column_label, const std::string& base_source,
        bool yearDecorrelated, const std::string& direction) {
    if (direction != "Up" && direction != "Dn") return df;
    const float sign = (direction == "Up") ? +1.0f : -1.0f;
    const std::string sfx = "jes" + column_label + direction;
    const std::string shiftCol = "_Jet_shift_" + sfx;
    const std::string metPropCol = "_metProp_" + sfx;

    auto eval_factor = makeJESShiftEvaluator(cset_jerc, jec_prefix_map, jec_suffix_map,
                                             year_token_map, base_source, yearDecorrelated, sign);
    auto met_prop = jesMetPropagator();

    return df
        .Define(shiftCol, eval_factor, {"year", "Jet_pt", "Jet_eta"})
        .Define("Jet_pt_"   + sfx, "Jet_pt   * " + shiftCol)
        .Define("Jet_mass_" + sfx, "Jet_mass * " + shiftCol)
        .Define(metPropCol, met_prop,
                {"Jet_pt", "Jet_pt_" + sfx, "Jet_phi", "Jet_neEmEF", "Jet_chEmEF",
                 "met_pt", "met_phi"})
        .Define("met_pt_"  + sfx, [](std::pair<float,float> p){ return p.first;  }, {metPropCol})
        .Define("met_phi_" + sfx, [](std::pair<float,float> p){ return p.second; }, {metPropCol});
}

RNode defineFatJetJESVariation(
        const std::unordered_map<std::string, correction::CorrectionSet>& cset_jerc,
        const std::unordered_map<std::string, std::string>& jec_prefix_map,
        const std::unordered_map<std::string, std::string>& jec_suffix_map,
        const std::unordered_map<std::string, std::string>& year_token_map,
        RNode df, const std::string& column_label, const std::string& base_source,
        bool yearDecorrelated, const std::string& direction) {
    if (direction != "Up" && direction != "Dn") return df;
    const float sign = (direction == "Up") ? +1.0f : -1.0f;
    const std::string sfx = "jes" + column_label + direction;
    const std::string shiftCol = "_FatJet_shift_" + sfx;

    auto eval_factor = makeJESShiftEvaluator(cset_jerc, jec_prefix_map, jec_suffix_map,
                                             year_token_map, base_source, yearDecorrelated, sign);

    return df
        .Define(shiftCol, eval_factor, {"year", "FatJet_pt", "FatJet_eta"})
        .Define("FatJet_pt_"   + sfx, "FatJet_pt   * " + shiftCol)
        .Define("FatJet_mass_" + sfx, "FatJet_mass * " + shiftCol);
}

namespace {
struct JESSourceSpec {
    const char* label;          // column-name suffix, e.g. "Absolute" or "AbsoluteYear"
    const char* base_source;    // JEC source key root, e.g. "Regrouped_Absolute"
    bool yearDecorrelated;      // append _<year_token> at lookup time?
};

// 11 Regrouped V2 sources — confirmed against jet_jerc.json.gz / fatJet_jerc.json.gz
// for all 9 campaigns on 2026-05-14. JME-POG standard set; combine treats the *_Year
// entries as uncorrelated across years (handled downstream via the `year` hist axis).
static const std::vector<JESSourceSpec> kJESRegroupedSources = {
    {"Absolute",            "Regrouped_Absolute",       false},
    {"BBEC1",               "Regrouped_BBEC1",          false},
    {"EC2",                 "Regrouped_EC2",            false},
    {"HF",                  "Regrouped_HF",             false},
    {"RelativeBal",         "Regrouped_RelativeBal",    false},
    {"FlavorQCD",           "Regrouped_FlavorQCD",      false},
    {"AbsoluteYear",        "Regrouped_Absolute",       true},
    {"BBEC1Year",           "Regrouped_BBEC1",          true},
    {"EC2Year",             "Regrouped_EC2",            true},
    {"HFYear",              "Regrouped_HF",             true},
    {"RelativeSampleYear",  "Regrouped_RelativeSample", true},
};
static const std::array<const char*, 2> kJESDirections = {"Up", "Dn"};
} // anonymous namespace

RNode applyJESVariations(RNode df) {
    for (const auto& src : kJESRegroupedSources) {
        for (const char* dir : kJESDirections) {
            df = defineAK4JESVariation(
                jetEnergyCorrections,
                jetEnergyCorrections_JEC_prefix,
                jetEnergyCorrections_JEC_suffix,
                jetEnergyCorrections_yearToken,
                df, src.label, src.base_source, src.yearDecorrelated, dir);
            df = defineFatJetJESVariation(
                fatJetEnergyCorrections,
                fatJetEnergyCorrections_JEC_prefix,
                fatJetEnergyCorrections_JEC_suffix,
                jetEnergyCorrections_yearToken,
                df, src.label, src.base_source, src.yearDecorrelated, dir);
        }
    }
    return df;
}

static bool g_storeSysts = true;
void setStoreSysts(bool v) { g_storeSysts = v; }

std::vector<std::string> jesVariationSuffixes() {
    if (!g_storeSysts) return {};
    std::vector<std::string> out;
    out.reserve(kJESRegroupedSources.size() * kJESDirections.size());
    for (const auto& src : kJESRegroupedSources) {
        for (const char* dir : kJESDirections) {
            out.push_back(std::string("jes") + src.label + dir);
        }
    }
    return out;
}

/*
############################################
JET ENERGY RESOLUTION — hybrid smearing (MC only)
############################################

Per AK4 jet (using the JEC-corrected pt that this function consumes):
  res = PtResolution(eta, pt, rho)                        // relative resolution
  sf  = ScaleFactor(eta, "nom"|"up"|"down")               // data/MC width ratio
  pt_smear = JERSmear(pt, eta, gen_pt_or_-1, rho, event, res, sf)
               // hybrid: scaling when |pt - pt_gen| < 3*res*pt and gen-matched (genJetIdx>=0),
               // stochastic Gaussian otherwise. Always >= 0.
  pt_new   = pt * pt_smear ; mass_new = mass * pt_smear

Propagates the smearing to met_pt / met_phi via Type-I-style recipe.
*/

RNode applyJetEnergyResolution(const std::unordered_map<std::string, correction::CorrectionSet>& cset_jerc,
                               const std::unordered_map<std::string, correction::CorrectionSet>& cset_jer_smear,
                               const std::unordered_map<std::string, std::string>& jer_res_map,
                               const std::unordered_map<std::string, std::string>& jer_sf_map,
                               RNode df, std::string variation) {

    const std::string vary = (variation == "nominal") ? "nom" : variation;

    auto compute_smear = [cset_jerc, cset_jer_smear, jer_res_map, jer_sf_map, vary](
            std::string year,
            RVec<float> pt, RVec<float> eta,
            RVec<short> genJet_idx, RVec<float> genJet_pt,
            float rho, unsigned long long event) {

        RVec<float> factor(pt.size(), 1.0f);
        if (pt.empty()) return factor;

        auto cs_it = cset_jerc.find(year);
        auto sm_it = cset_jer_smear.find("jer_smear");
        auto rn_it = jer_res_map.find(year);
        auto sn_it = jer_sf_map.find(year);
        if (cs_it == cset_jerc.end() || sm_it == cset_jer_smear.end()
            || rn_it == jer_res_map.end() || sn_it == jer_sf_map.end()) {
            static std::unordered_set<std::string> warned;
            if (warned.insert(year).second) {
                std::cout << "Warning: JER inputs missing for year " << year
                          << ". Skipping JER smearing." << std::endl;
            }
            return factor;
        }

        std::shared_ptr<const correction::Correction> res, sf, smear;
        try {
            res   = cs_it->second.at(rn_it->second);
            sf    = cs_it->second.at(sn_it->second);
            smear = sm_it->second.at("JERSmear");
        } catch (const std::out_of_range&) {
            static std::unordered_set<std::string> warned;
            if (warned.insert(year).second) {
                std::cout << "Warning: JER res/sf/smear key missing for year " << year
                          << " (res=" << rn_it->second << " sf=" << sn_it->second
                          << "). Skipping JER smearing." << std::endl;
            }
            return factor;
        }

        const bool sf_has_pt = (sf->inputs().size() >= 3);
        for (size_t i = 0; i < pt.size(); ++i) {
            const float gen_pt = (genJet_idx[i] >= 0 && genJet_idx[i] < (int)genJet_pt.size())
                                     ? genJet_pt[genJet_idx[i]] : -1.0f;
            double r  = res->evaluate({(double)eta[i], (double)pt[i], (double)rho});
            double s  = sf_has_pt ? sf->evaluate({(double)eta[i], (double)pt[i], vary})
                                  : sf->evaluate({(double)eta[i], vary});
            double sm = smear->evaluate({(double)pt[i], (double)eta[i], (double)gen_pt,
                                         (double)rho, (int)event, r, s});
            factor[i] = static_cast<float>(sm);
        }
        return factor;
    };

    auto met_propagate_jer = [](RVec<float> pt_pre, RVec<float> factor, RVec<float> phi,
                                RVec<float> neEmEF, RVec<float> chEmEF,
                                float met_pt, float met_phi) {
        double px = met_pt * std::cos(met_phi);
        double py = met_pt * std::sin(met_phi);
        for (size_t i = 0; i < pt_pre.size(); ++i) {
            const float pt_new = pt_pre[i] * factor[i];
            if (pt_new <= 15.0f) continue;
            if ((neEmEF[i] + chEmEF[i]) > 0.9f) continue;
            const double dpt = static_cast<double>(pt_pre[i] - pt_new);
            px += dpt * std::cos(phi[i]);
            py += dpt * std::sin(phi[i]);
        }
        return std::make_pair(static_cast<float>(std::sqrt(px*px + py*py)),
                              static_cast<float>(std::atan2(py, px)));
    };

    return df
        .Define("Jet_jerFactor", compute_smear,
                {"year", "Jet_pt", "Jet_eta", "Jet_genJetIdx", "GenJet_pt",
                 "Rho_fixedGridRhoFastjetAll", "event"})
        .Define("_Jet_pt_preJER", "Jet_pt")
        .Define("_jer_metProp", met_propagate_jer,
                {"_Jet_pt_preJER", "Jet_jerFactor", "Jet_phi",
                 "Jet_neEmEF", "Jet_chEmEF", "met_pt", "met_phi"})
        .Redefine("met_pt",   "_jer_metProp.first")
        .Redefine("met_phi",  "_jer_metProp.second")
        .Redefine("Jet_pt",   "Jet_pt   * Jet_jerFactor")
        .Redefine("Jet_mass", "Jet_mass * Jet_jerFactor");
}

/*
############################################
FAT JET (AK8) ENERGY RESOLUTION — hybrid smearing (MC only)
############################################

Same recipe as AK4 but reading FatJet_* + GenJetAK8_* and the AK8PFPuppi resolution/SF
keys. Does NOT propagate to MET.
*/

RNode applyFatJetEnergyResolution(const std::unordered_map<std::string, correction::CorrectionSet>& cset_jerc,
                                  const std::unordered_map<std::string, correction::CorrectionSet>& cset_jer_smear,
                                  const std::unordered_map<std::string, std::string>& jer_res_map,
                                  const std::unordered_map<std::string, std::string>& jer_sf_map,
                                  RNode df, std::string variation) {

    const std::string vary = (variation == "nominal") ? "nom" : variation;

    auto compute_smear = [cset_jerc, cset_jer_smear, jer_res_map, jer_sf_map, vary](
            std::string year,
            RVec<float> pt, RVec<float> eta,
            RVec<short> genJet_idx, RVec<float> genJet_pt,
            float rho, unsigned long long event) {

        RVec<float> factor(pt.size(), 1.0f);
        if (pt.empty()) return factor;

        auto cs_it = cset_jerc.find(year);
        auto sm_it = cset_jer_smear.find("jer_smear");
        auto rn_it = jer_res_map.find(year);
        auto sn_it = jer_sf_map.find(year);
        if (cs_it == cset_jerc.end() || sm_it == cset_jer_smear.end()
            || rn_it == jer_res_map.end() || sn_it == jer_sf_map.end()) {
            static std::unordered_set<std::string> warned;
            if (warned.insert(year).second) {
                std::cout << "Warning: AK8 JER inputs missing for year " << year
                          << ". Skipping AK8 JER smearing." << std::endl;
            }
            return factor;
        }

        std::shared_ptr<const correction::Correction> res, sf, smear;
        try {
            res   = cs_it->second.at(rn_it->second);
            sf    = cs_it->second.at(sn_it->second);
            smear = sm_it->second.at("JERSmear");
        } catch (const std::out_of_range&) {
            static std::unordered_set<std::string> warned;
            if (warned.insert(year).second) {
                std::cout << "Warning: AK8 JER res/sf/smear key missing for year " << year
                          << " (res=" << rn_it->second << " sf=" << sn_it->second
                          << "). Skipping AK8 JER smearing." << std::endl;
            }
            return factor;
        }

        const bool sf_has_pt = (sf->inputs().size() >= 3);
        for (size_t i = 0; i < pt.size(); ++i) {
            const float gen_pt = (genJet_idx[i] >= 0 && genJet_idx[i] < (int)genJet_pt.size())
                                     ? genJet_pt[genJet_idx[i]] : -1.0f;
            double r  = res->evaluate({(double)eta[i], (double)pt[i], (double)rho});
            double s  = sf_has_pt ? sf->evaluate({(double)eta[i], (double)pt[i], vary})
                                  : sf->evaluate({(double)eta[i], vary});
            double sm = smear->evaluate({(double)pt[i], (double)eta[i], (double)gen_pt,
                                         (double)rho, (int)event, r, s});
            factor[i] = static_cast<float>(sm);
        }
        return factor;
    };

    return df
        .Define("FatJet_jerFactor", compute_smear,
                {"year", "FatJet_pt", "FatJet_eta", "FatJet_genJetAK8Idx", "GenJetAK8_pt",
                 "Rho_fixedGridRhoFastjetAll", "event"})
        .Redefine("FatJet_pt",   "FatJet_pt   * FatJet_jerFactor")
        .Redefine("FatJet_mass", "FatJet_mass * FatJet_jerFactor");
}

/*
############################################
JET VETO MAPS
############################################
*/

RNode applyJetVetoMaps(RNode df) {
    auto eval_correction = [] (std::string year, RVec<float> pt, RVec<float> eta, RVec<float> phi, RVec<float> jet_id, RVec<float> jet_nuEmEF, RVec<float> jet_chEmEF) {
        RVec<bool> jet_veto_map;
        
        if (jetVetoMaps.find(year) == jetVetoMaps.end()) {
            static std::unordered_set<std::string> warned_years;
            if (warned_years.find(year) == warned_years.end()) {
                std::cout << "Warning: Jet veto map for year " << year << " not found. Setting all jets to not vetoed." << std::endl;
                warned_years.insert(year);
            }
            for (size_t i = 0; i < eta.size(); i++) {
                jet_veto_map.push_back(false);
            }
            return jet_veto_map;
        }

        for (size_t i = 0; i < eta.size(); i++) {
            float eta_ = eta[i];
            if (std::abs(eta_) > 5.19) {
                eta_ = 5.19 * (eta_ > 0 ? 1 : -1);
            }
            bool is_vetoed = jetVetoMaps.at(year).at(jetVetoMap_names.at(year))->evaluate({"jetvetomap", eta_, phi[i]}) != 0;
            if (is_vetoed && (pt[i] > 15.0 && jet_id[i] == 6 && (jet_nuEmEF[i] + jet_chEmEF[i]) < 0.9)) {
                jet_veto_map.push_back(true);
            } else {
                jet_veto_map.push_back(false);
            }
        }

        return jet_veto_map;
    };
    
    return df.Define("Jet_vetoMap", eval_correction, {"year", "Jet_pt", "Jet_eta", "Jet_phi", "Jet_jetId", "Jet_neEmEF", "Jet_chEmEF"});
}

/*
############################################
ELECTRON SCALE AND SMEARING CORRECTIONS
############################################
*/

RNode applyElectronScaleAndSmearing(RNode df, bool isData) {
    auto eval_data_scale = [](std::string year, RVec<float> pt, RVec<float> eta, RVec<float> deltaEtaSC, RVec<float> r9, RVec<unsigned char> seedGain, unsigned int run) {
        RVec<float> scales(pt.size(), 1.0);
        if (electronSSCorrections.find(year) == electronSSCorrections.end()) {
            static std::unordered_set<std::string> warned_years;
            if (warned_years.find(year) == warned_years.end()) {
                std::cout << "Warning: Electron SS correction set for year " << year << " not found. Skipping scale corrections." << std::endl;
                warned_years.insert(year);
            }
            return scales;
        }

        for (size_t i = 0; i < pt.size(); i++) {
            float scEta = eta[i] + deltaEtaSC[i];
            scales[i] = electronSSCorrections.at(year).compound().at("Scale")->evaluate({"scale", static_cast<double>(run), scEta, r9[i], pt[i], static_cast<double>(seedGain[i])}); 
        }
        return scales;
    };

    auto eval_data_smearWidth = [](std::string year, RVec<float> pt, RVec<float> eta, RVec<float> deltaEtaSC, RVec<float> r9, RVec<float> scales) {
        RVec<float> widths(pt.size(), 0.0);
        if (electronSSCorrections.find(year) == electronSSCorrections.end()) {
            return widths;
        }

        for (size_t i = 0; i < pt.size(); i++) {
            float scEta = eta[i] + deltaEtaSC[i];
            float pt_corr = pt[i] * scales[i];
            widths[i] = electronSSCorrections.at(year).at("SmearAndSyst")->evaluate({"smear", pt_corr, r9[i], scEta});
        }
        return widths;
    };

    auto eval_mc_smear_factor = [](std::string year, RVec<float> pt, RVec<float> eta, RVec<float> deltaEtaSC, RVec<float> r9, unsigned long long event) {
        RVec<float> smear_factors(pt.size(), 1.0);
        if (electronSSCorrections.find(year) == electronSSCorrections.end()) {
            static std::unordered_set<std::string> warned_years;
            if (warned_years.find(year) == warned_years.end()) {
                std::cout << "Warning: Electron SS correction set for year " << year << " not found. Skipping smearing corrections." << std::endl;
                warned_years.insert(year);
            }
            return smear_factors;
        }

        TRandom3 rng(event);
        for (size_t i = 0; i < pt.size(); i++) {
            float scEta = eta[i] + deltaEtaSC[i];
            float smear = electronSSCorrections.at(year).at("SmearAndSyst")->evaluate({"smear", pt[i], r9[i], scEta});
            smear_factors[i] = 1.0 + smear * rng.Gaus(0.0, 1.0);
        }
        return smear_factors;
    };

    auto eval_mc_smearWidth = [](std::string year, RVec<float> pt, RVec<float> eta, RVec<float> deltaEtaSC, RVec<float> r9) {
        RVec<float> widths(pt.size(), 0.0);
        if (electronSSCorrections.find(year) == electronSSCorrections.end()) {
            return widths;
        }

        for (size_t i = 0; i < pt.size(); i++) {
            float scEta = eta[i] + deltaEtaSC[i];
            widths[i] = electronSSCorrections.at(year).at("SmearAndSyst")->evaluate({"smear", pt[i], r9[i], scEta});
        }
        return widths;
    };

    auto calc_energyErr = [](RVec<float> pt, RVec<float> eta, RVec<float> energyErr, RVec<float> widths, RVec<float> factors) {
        RVec<float> new_energyErr(pt.size(), 0.0);
        for (size_t i = 0; i < pt.size(); i++) {
            float energy = pt[i] * TMath::CosH(eta[i]);
            new_energyErr[i] = TMath::Sqrt(energyErr[i] * energyErr[i] + (energy * widths[i]) * (energy * widths[i])) * factors[i];
        }
        return new_energyErr;
    };

    if (isData) {
        return df.Define("_ele_scale", eval_data_scale, {"year", "Electron_pt", "Electron_eta", "Electron_deltaEtaSC", "Electron_r9", "Electron_seedGain", "run"})
                 .Define("_ele_smearWidth", eval_data_smearWidth, {"year", "Electron_pt", "Electron_eta", "Electron_deltaEtaSC", "Electron_r9", "_ele_scale"})
                 .Redefine("Electron_energyErr", calc_energyErr, {"Electron_pt", "Electron_eta", "Electron_energyErr", "_ele_smearWidth", "_ele_scale"})
                 .Redefine("Electron_pt", "Electron_pt * _ele_scale")
                 .Redefine("Electron_mass", "Electron_mass * _ele_scale");
    } else {
        return df.Define("_ele_smearFactor", eval_mc_smear_factor, {"year", "Electron_pt", "Electron_eta", "Electron_deltaEtaSC", "Electron_r9", "event"})
                 .Define("_ele_smearWidth", eval_mc_smearWidth, {"year", "Electron_pt", "Electron_eta", "Electron_deltaEtaSC", "Electron_r9"})
                 .Redefine("Electron_energyErr", calc_energyErr, {"Electron_pt", "Electron_eta", "Electron_energyErr", "_ele_smearWidth", "_ele_smearFactor"})
                 .Redefine("Electron_pt", "Electron_pt * _ele_smearFactor")
                 .Redefine("Electron_mass", "Electron_mass * _ele_smearFactor");
    }
}

/*
############################################
GENERAL CORRECTIONS
############################################
*/

RNode applyDataCorrections(RNode df_) {
    auto df = applyMETPhiCorrections(df_, true);
    // Re-apply the latest JEC stack (raw recovery + L1*L2*L3*Residual compound).
    // AK4 propagates to MET; AK8 does not (PuppiMET only includes AK4 Type-I).
    df = applyJetEnergyCorrections(jetEnergyCorrections,
                                   jetEnergyCorrections_JEC_prefix,
                                   jetEnergyCorrections_JEC_suffix,
                                   df, /*isData=*/true);
    df = applyFatJetEnergyCorrections(fatJetEnergyCorrections,
                                      fatJetEnergyCorrections_JEC_prefix,
                                      fatJetEnergyCorrections_JEC_suffix,
                                      df, /*isData=*/true);
    df = HEMCorrection(df, true);
    df = applyElectronScaleAndSmearing(df, true);
    return df;
}

RNode applyMCCorrections(RNode df_) {
    auto df = applyMETPhiCorrections(df_, false);
    // Nominal JEC for AK4 and AK8.
    df = applyJetEnergyCorrections(jetEnergyCorrections,
                                   jetEnergyCorrections_JEC_prefix,
                                   jetEnergyCorrections_JEC_suffix,
                                   df, /*isData=*/false);
    df = applyFatJetEnergyCorrections(fatJetEnergyCorrections,
                                      fatJetEnergyCorrections_JEC_prefix,
                                      fatJetEnergyCorrections_JEC_suffix,
                                      df, /*isData=*/false);
    // JER smearing on top of JEC-corrected jets.
    df = applyJetEnergyResolution(jetEnergyCorrections,
                                  jetEnergyResolution_smear,
                                  jetEnergyResolution_JER_res_name,
                                  jetEnergyResolution_JER_sf_name,
                                  df, /*variation=*/"nominal");
    df = applyFatJetEnergyResolution(fatJetEnergyCorrections,
                                     jetEnergyResolution_smear,
                                     fatJetEnergyResolution_JER_res_name,
                                     fatJetEnergyResolution_JER_sf_name,
                                     df, /*variation=*/"nominal");
    // JES uncertainties — 11 Regrouped V2 sources × {Up, Dn} as suffixed branches.
    // Must run after the nominal JEC+JER so the shift sits on top of the corrected baseline.
    if (g_storeSysts) df = applyJESVariations(df);
    // JMS / JMR on FatJet_msoftdrop. Currently no-op: centrals are placeholders
    // (shift = 0 GeV, factor = 1.0). Replace jetMassScale_central / jetMassResolution_central
    // in corrections.h with the values derived from the per-analysis calibration fit.
    df = applyJetMassScale(jetMassScale_central, df);
    df = applyJetMassResolution(jetMassResolution_central, jetMassResolution_sigmaRel_central, df);
    df = HEMCorrection(df, false);
    df = applyElectronScaleAndSmearing(df, false);
    return df;
}
