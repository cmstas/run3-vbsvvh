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

// /*
// ############################################
// JET MASS SCALE CORRECTIONS
// ############################################
// */

// RNode applyJMSCorrections(std::unordered_map<std::string, correction::CorrectionSet> cset_jms, RNode df, std::string variation) { 
//     auto eval_correction = [cset_jms, variation] (std::string year, float mass) {
//         if (cset_jms.find(year) == cset_jms.end()) {
//             return mass;
//         }
//         double scaleVal = 1. + 0.05 * cset_jms.at(year).at("JMS")->evaluate({year, variation});
//         // https://docs.google.com/presentation/d/1C7CqO3Wv3-lYd7vw4IQXq69wmULesTsSniFXXM__ReU
//         return mass * scaleVal;
//     };
//     return df.Redefine("Hbbmass", eval_correction, {"year", "GHiggs_mass"})
//              .Redefine("Wjetmass", eval_correction, {"year", "GW_mass"});
// }

// /*
// ############################################
// JET MASS RESOLUTION CORRECTIONS
// ############################################
// */

// RNode applyJMRCorrections(std::unordered_map<std::string, correction::CorrectionSet> cset_jmr, RNode df, std::string variation) {
//     auto eval_correction = [cset_jmr, variation] (std::string year, float mass, unsigned int lumi, unsigned long long event) {
//         if (cset_jmr.find(year) == cset_jmr.end()) {
//             return mass;
//         }
//         TRandom3 rnd((lumi << 10) + event);
//         return rnd.Gaus(1, 0.1 * cset_jmr.at(year).at("JMR")->evaluate({year, variation})) * mass;
//     };
//     return df.Redefine("Hbbmass", eval_correction, {"year", "GHiggs_mass", "luminosityBlock", "event"})
//              .Redefine("Wjetmass", eval_correction, {"year", "GW_mass", "luminosityBlock", "event"});
// }

/*
############################################
HEM Corrections
############################################
*/

RNode HEMCorrection(RNode df, bool isData) {
    auto HEMCorrections = [isData](unsigned int run, unsigned long long event, std::string sample_year, RVec<float> pt, RVec<float> eta, RVec<float> phi, RVec<float> jet_id) {
        RVec<bool> jet_mask;
        if (sample_year.find("2016") != std::string::npos) {
            jet_mask = (jet_id >= 1 && pt > 15.0);
        } else {
            jet_mask = (jet_id >= 2 && pt > 15.0);
        }
        // Need to check if there is dependence on jet ID in v15
        auto eta_ = eta[jet_mask];
        auto phi_ = phi[jet_mask];
        if (sample_year == "2018" && ((isData && run >= 319077) || (!isData && event % 1961 < 1286))) {
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
            return std::make_pair(pt_corr, phi_corr);
        }
        
        std::string pt_corr_name = isData ? "pt_metphicorr_puppimet_data" : "pt_metphicorr_puppimet_mc";
        std::string phi_corr_name = isData ? "phi_metphicorr_puppimet_data" : "phi_metphicorr_puppimet_mc";
        
        pt_corr = metCorrections.at(year).at(pt_corr_name)->evaluate({pt, phi, static_cast<double>(npvs), static_cast<double>(run)});
        phi_corr = metCorrections.at(year).at(phi_corr_name)->evaluate({pt, phi, static_cast<double>(npvs), static_cast<double>(run)});
        
        return std::make_pair(pt_corr, phi_corr);
    };
    return df.Define("_MET_phicorr", eval_correction, {"year", "PuppiMET_pt", "PuppiMET_phi", "PV_npvs", "run"})
            .Redefine("PuppiMET_pt", "_MET_phicorr.first")
            .Redefine("PuppiMET_phi", "_MET_phicorr.second");
}

/*
############################################
MET UNCLUSTERED CORRECTIONS
############################################
*/

RNode applyMETUnclusteredCorrections(RNode df, std::string variation) {
    if (variation == "up") {
        return df.Define("_MET_uncert_dx", "PuppiMET_pt * TMath::Cos(PuppiMET_phi) + MET_MetUnclustEnUpDeltaX")
                .Define("_MET_uncert_dy", "PuppiMET_pt * TMath::Sin(PuppiMET_phi) + MET_MetUnclustEnUpDeltaY")
                .Redefine("PuppiMET_pt", "TMath::Sqrt(_MET_uncert_dx * _MET_uncert_dx + _MET_uncert_dy * _MET_uncert_dy)");
    }
    else if (variation == "down") {
        return df.Define("_MET_uncert_dx", "PuppiMET_pt * TMath::Cos(PuppiMET_phi) - MET_MetUnclustEnUpDeltaX")
                .Define("_MET_uncert_dy", "PuppiMET_pt * TMath::Sin(PuppiMET_phi) - MET_MetUnclustEnUpDeltaY")
                .Redefine("PuppiMET_pt", "TMath::Sqrt(_MET_uncert_dx * _MET_uncert_dx + _MET_uncert_dy * _MET_uncert_dy)");
    }
    return df;
}

/*
############################################
JET ENERGY CORRECTIONS
############################################
*/

RNode applyJetEnergyCorrections(std::unordered_map<std::string, correction::CorrectionSet> cset_jerc, std::unordered_map<std::string, std::string> jec_prefix_map, std::unordered_map<std::string, std::string> jec_suffix_map, RNode df, std::string JEC_type, std::string variation) {
    auto eval_correction = [cset_jerc, jec_prefix_map, jec_suffix_map, JEC_type, variation] (std::string year, RVec<float> pt, RVec<float> eta, RVec<float> var) {
        RVec<float> jec_factors;

        if (var.size() == 0) {
            return var;
        }
        
        if (cset_jerc.find(year) == cset_jerc.end()) {
            static std::unordered_set<std::string> warned_years;
            if (warned_years.find(year) == warned_years.end()) {
                std::cout << "Warning: JEC correction set for year " << year << " not found. Skipping JEC corrections." << std::endl;
                warned_years.insert(year);
            }
            return var;
        }
        
        if (jec_prefix_map.find(year) == jec_prefix_map.end()) {
            static std::unordered_set<std::string> warned_years;
            if (warned_years.find(year) == warned_years.end()) {
                std::cout << "Warning: JEC prefix map for year " << year << " not found. Skipping JEC corrections." << std::endl;
                warned_years.insert(year);
            }
            return var;
        }
        
        if (jec_suffix_map.find(year) == jec_suffix_map.end()) {
            static std::unordered_set<std::string> warned_years;
            if (warned_years.find(year) == warned_years.end()) {
                std::cout << "Warning: JEC suffix map for year " << year << " not found. Skipping JEC corrections." << std::endl;
                warned_years.insert(year);
            }
            return var;
        }

        // Construct JEC correction name based on year
        std::string JEC_name = jec_prefix_map.at(year) + "_" + JEC_type + "_" + jec_suffix_map.at(year);

        for (size_t i = 0; i < var.size(); i++) {
            if (variation == "up") {
                jec_factors.push_back(var[i] * (1 + cset_jerc.at(year).at(JEC_name)->evaluate({eta[i], pt[i]})));
            }
            else if (variation == "down") {
                jec_factors.push_back(var[i] * (1 - cset_jerc.at(year).at(JEC_name)->evaluate({eta[i], pt[i]})));
            }
            else {
                return var;
            }
        }
        return jec_factors;
    };
    
    auto df_jetcorr = df.Redefine("Jet_pt", eval_correction, {"year", "Jet_pt", "Jet_eta", "Jet_pt"})
                        .Redefine("Jet_mass", eval_correction, {"year", "Jet_pt", "Jet_eta", "Jet_mass"})
                        .Redefine("FatJet_pt", eval_correction, {"year", "FatJet_pt", "FatJet_eta", "FatJet_pt"})
                        .Redefine("FatJet_mass", eval_correction, {"year", "FatJet_pt", "FatJet_eta", "FatJet_mass"});

    auto correctmet = [JEC_type](std::string year, RVec<float> Jet_pt, RVec<float> jet_phi, RVec<float> jet_pt, float PuppiMET_pt, float PuppiMET_phi) {
        if (Jet_pt.empty()) {
            return PuppiMET_pt;
        }
        float px = PuppiMET_pt * TMath::Cos(PuppiMET_phi);
        float py = PuppiMET_pt * TMath::Sin(PuppiMET_phi);
        for (size_t i = 0; i < Jet_pt.size(); i++) {
            px += (Jet_pt[i] - jet_pt[i]) * TMath::Cos(jet_phi[i]);
            py += (Jet_pt[i] - jet_pt[i]) * TMath::Sin(jet_phi[i]);
        }
        return (float)TMath::Sqrt(px * px + py * py);
    };

    return df_jetcorr.Redefine("PuppiMET_pt", correctmet, {"year", "Jet_pt", "Jet_phi", "Jet_pt", "PuppiMET_pt", "PuppiMET_phi"});
}

/*
############################################
JET ENERGY RESOLUTION
############################################
*/

RNode applyJetEnergyResolution(std::unordered_map<std::string, correction::CorrectionSet> cset_jerc, std::unordered_map<std::string, correction::CorrectionSet> cset_jer_smear, std::unordered_map<std::string, std::string> jer_res_map, std::unordered_map<std::string, std::string> jer_sf_map, RNode df, std::string variation) {
    auto eval_correction = [cset_jerc, cset_jer_smear, jer_res_map, jer_sf_map, variation] (std::string year, RVec<float> pt, RVec<float> eta, RVec<int> genJet_idx, RVec<float> genJet_pt, float rho, unsigned long long event, RVec<float> var) {
        RVec<float> jer_factors;
        std::string vary = (variation == "nominal") ? "nom" : variation;
        
        if (var.size() == 0) {
            return var;
        }
        
        if (cset_jerc.find(year) == cset_jerc.end()) {
            static std::unordered_set<std::string> warned_years;
            if (warned_years.find(year) == warned_years.end()) {
                std::cout << "Warning: JER correction set for year " << year << " not found. Skipping JER corrections." << std::endl;
                warned_years.insert(year);
            }
            return var;
        }
        
        if (cset_jer_smear.find("jer_smear") == cset_jer_smear.end()) {
            static bool warned = false;
            if (!warned) {
                std::cout << "Warning: JER smear correction set not found. Skipping JER corrections." << std::endl;
                warned = true;
            }
            return var;
        }
        
        if (jer_res_map.find(year) == jer_res_map.end()) {
            static std::unordered_set<std::string> warned_years;
            if (warned_years.find(year) == warned_years.end()) {
                std::cout << "Warning: JER resolution map for year " << year << " not found. Skipping JER corrections." << std::endl;
                warned_years.insert(year);
            }
            return var;
        }
        
        if (jer_sf_map.find(year) == jer_sf_map.end()) {
            static std::unordered_set<std::string> warned_years;
            if (warned_years.find(year) == warned_years.end()) {
                std::cout << "Warning: JER scale factor map for year " << year << " not found. Skipping JER corrections." << std::endl;
                warned_years.insert(year);
            }
            return var;
        }

        std::string jer_res_name = jer_res_map.at(year);
        std::string jer_sf_name = jer_sf_map.at(year);

        for (size_t i = 0; i < var.size(); i++) {
            float genjetpt = genJet_idx[i] >= 0 ? genJet_pt[genJet_idx[i]] : -1;
            float jer = cset_jerc.at(year).at(jer_res_name)->evaluate({eta[i], pt[i], rho});
            float jer_sf = cset_jerc.at(year).at(jer_sf_name)->evaluate({eta[i], vary});
            jer_factors.push_back(var[i] * cset_jer_smear.at("jer_smear").at("JERSmear")->evaluate({pt[i], eta[i], genjetpt, rho, (int)event, jer, jer_sf}));
        }
        
        return jer_factors;
    };
    
    return df.Redefine("Jet_pt", eval_correction, {"year", "Jet_pt", "Jet_eta", "Jet_genJetIdx", "GenJet_pt", "fixedGridRhoFastjetAll", "event", "Jet_pt"})
            .Redefine("Jet_mass", eval_correction, {"year", "Jet_pt", "Jet_eta", "Jet_genJetIdx", "GenJet_pt", "fixedGridRhoFastjetAll", "event", "Jet_mass"});
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
GENERAL CORRECTIONS
############################################
*/

RNode applyDataCorrections(RNode df_) {
    auto df = applyMETPhiCorrections(df_, true);
    //df = HEMCorrection(df, true);
    return df;
}

RNode applyMCCorrections(RNode df_) {
    auto df = applyMETPhiCorrections(df_, false);
    //df = HEMCorrection(df, false);
    return df;
}
