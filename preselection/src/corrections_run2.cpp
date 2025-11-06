#include "corrections_run2.h"

namespace Run2
{

float looseDFBtagWP(std::string year){
    if(year == "2016preVFP")
        return 0.0508;
    if(year == "2016postVFP")
        return 0.0480;
    if(year == "2017")
        return 0.0532;
    if(year == "2018")
        return 0.0490;
    return -1;
}

float mediumDFBtagWP(std::string year){
    if(year == "2016preVFP")
        return 0.2598;
    if(year == "2016postVFP")
        return 0.2489;
    if(year == "2017")
        return 0.3040;
    if(year == "2018")
        return 0.2783;
    return -1;
}

float tightDFBtagWP(std::string year){
    if(year == "2016preVFP")
        return 0.6502;
    if(year == "2016postVFP")
        return 0.6377;
    if(year == "2017")
        return 0.7476;
    if(year == "2018")
        return 0.7100;
    return -1;
}

// CORRECTIONS
RNode defineCorrectedCols(RNode df) {
    return df.Define("CorrJet_pt", "Jet_pt")
            .Define("CorrJet_mass", "Jet_mass")
            .Define("CorrFatJet_pt", "FatJet_pt")
            .Define("CorrFatJet_mass", "FatJet_mass")
            .Define("CorrMET_pt", "MET_pt")
            .Define("CorrMET_phi", "MET_phi");
}

RNode HEMCorrection(RNode df) {
    auto df1 = df.Define("HEMJets", "CorrJet_pt > 15 && "
                "((is2016 && Jet_jetId >= 1) || (!is2016 && Jet_jetId >= 2)) && "
                "(CorrJet_pt >= 50 || (CorrJet_pt < 50 && Jet_puId != 0))")
            .Define("HEMJet_pt", "CorrJet_pt[HEMJets]")
            .Define("HEMJet_eta", "Jet_eta[HEMJets]")
            .Define("HEMJet_phi", "Jet_phi[HEMJets]");

    auto HEMCorrections = [](unsigned int run, unsigned long long event, std::string year, std::string category, RVec<float> pt, RVec<float> eta, RVec<float> phi) {
        bool isData = false;
        if (category == "data") isData = true;
        if (year == "2018" && ((isData && run >= 319077) || (!isData && event % 1961 < 1286))) {
            for (size_t i = 0; i < pt.size(); i++) {
                if (eta[i] > -3.2 && eta[i] < -1.3 && phi[i] > -1.57 && phi[i] < -0.87) {
                    return 0.0;
                }
            }
        }
        return 1.0;
    };

    return df1.Define("HEMweight", HEMCorrections, {"run", "event", "year", "category", "HEMJet_pt", "HEMJet_eta", "HEMJet_phi"});
}

RNode JMS_Corrections(correction::CorrectionSet cset_jet_mass_scale, RNode df, std::string variation) { 
    auto eval_correction = [cset_jet_mass_scale, variation] (std::string year, float mass) {
        double scaleVal = 1. + 0.05 * cset_jet_mass_scale.at("JMS")->evaluate({year, variation});
        // https://docs.google.com/presentation/d/1C7CqO3Wv3-lYd7vw4IQXq69wmULesTsSniFXXM__ReU
        return mass * scaleVal;
    };
    return df.Redefine("Hbbmass", eval_correction, {"year", "GHiggs_mass"})
             .Redefine("Wjetmass", eval_correction, {"year", "GW_mass"});
}

RNode JMR_Corrections(correction::CorrectionSet cset_jet_mass_resolution, RNode df, std::string variation) {
    auto eval_correction = [cset_jet_mass_resolution, variation] (std::string year, float mass, unsigned int lumi, unsigned long long event) {
        TRandom3 rnd((lumi << 10) + event);
        return rnd.Gaus(1, 0.1 * cset_jet_mass_resolution.at("JMR")->evaluate({year, variation})) * mass;
    };
    return df.Redefine("Hbbmass", eval_correction, {"year", "GHiggs_mass", "luminosityBlock", "event"})
             .Redefine("Wjetmass", eval_correction, {"year", "GW_mass", "luminosityBlock", "event"});
}

RNode METPhiCorrections(correction::CorrectionSet cset_met_2016preVFP, correction::CorrectionSet cset_met_2016postVFP, correction::CorrectionSet cset_met_2017, correction::CorrectionSet cset_met_2018, RNode df) {
    // auto eval_correction = [cset_met_2016preVFP, cset_met_2016postVFP, cset_met_2017, cset_met_2018] (std::string category, std::string year, float pt, float phi, int npvs, long long run) {
    //     bool isData = false;
    //     double pt_corr = 0, phi_corr = 0;
    //     if (category == "data") isData = true;
    //     if (year == "2016preVFP") {
    //         if (isData) {
    //             pt_corr = cset_met_2016preVFP.at("pt_metphicorr_pfmet_data")->evaluate({pt, phi, npvs, run});
    //             phi_corr = cset_met_2016preVFP.at("phi_metphicorr_pfmet_data")->evaluate({pt, phi, npvs, run});
    //         }
    //         else {
    //             pt_corr = cset_met_2016preVFP.at("pt_metphicorr_pfmet_mc")->evaluate({pt, phi, npvs, run});
    //             phi_corr = cset_met_2016preVFP.at("phi_metphicorr_pfmet_mc")->evaluate({pt, phi, npvs, run});
    //         }
    //     }
    //     else if (year == "2016postVFP") {
    //         if (isData) {
    //             pt_corr = cset_met_2016postVFP.at("pt_metphicorr_pfmet_data")->evaluate({pt, phi, npvs, run});
    //             phi_corr = cset_met_2016postVFP.at("phi_metphicorr_pfmet_data")->evaluate({pt, phi, npvs, run});
    //         }
    //         else {
    //             pt_corr = cset_met_2016postVFP.at("pt_metphicorr_pfmet_mc")->evaluate({pt, phi, npvs, run});
    //             phi_corr = cset_met_2016postVFP.at("phi_metphicorr_pfmet_mc")->evaluate({pt, phi, npvs, run});
    //         }
    //     }
    //     else if (year == "2017") {
    //         if (isData) {
    //             pt_corr = cset_met_2017.at("pt_metphicorr_pfmet_data")->evaluate({pt, phi, npvs, run});
    //             phi_corr = cset_met_2017.at("phi_metphicorr_pfmet_data")->evaluate({pt, phi, npvs, run});
    //         }
    //         else {
    //             pt_corr = cset_met_2017.at("pt_metphicorr_pfmet_mc")->evaluate({pt, phi, npvs, run});
    //             phi_corr = cset_met_2017.at("phi_metphicorr_pfmet_mc")->evaluate({pt, phi, npvs, run});
    //         }
    //     }
    //     else if (year == "2018") {
    //         if (isData) {
    //             pt_corr = cset_met_2018.at("pt_metphicorr_pfmet_data")->evaluate({pt, phi, npvs, run});
    //             phi_corr = cset_met_2018.at("phi_metphicorr_pfmet_data")->evaluate({pt, phi, npvs, run});
    //         }
    //         else {
    //             pt_corr = cset_met_2018.at("pt_metphicorr_pfmet_mc")->evaluate({pt, phi, npvs, run});
    //             phi_corr = cset_met_2018.at("phi_metphicorr_pfmet_mc")->evaluate({pt, phi, npvs, run});
    //         }
    //     }
    //     return std::make_pair(pt_corr, phi_corr);
    // };
    // return df.Define("MET_phicorr", eval_correction, {"category", "year", "MET_pt", "MET_phi", "PV_npvs", "run"})
    //         .Redefine("CorrMET_pt", "MET_phicorr.first")
    //         .Redefine("CorrMET_phi", "MET_phicorr.second");
    return df;
}

RNode METUnclusteredCorrections(RNode df, std::string variation) {
    if (variation == "up") {
        return df.Define("MET_uncert_dx", "CorrMET_pt * TMath::Cos(CorrMET_phi) + MET_MetUnclustEnUpDeltaX")
                .Define("MET_uncert_dy", "CorrMET_pt * TMath::Sin(CorrMET_phi) + MET_MetUnclustEnUpDeltaY")
                .Redefine("CorrMET_pt", "TMath::Sqrt(MET_uncert_dx * MET_uncert_dx + MET_uncert_dy * MET_uncert_dy)");
    }
    else if (variation == "down") {
        return df.Define("MET_uncert_dx", "CorrMET_pt * TMath::Cos(CorrMET_phi) - MET_MetUnclustEnUpDeltaX")
                .Define("MET_uncert_dy", "CorrMET_pt * TMath::Sin(CorrMET_phi) - MET_MetUnclustEnUpDeltaY")
                .Redefine("CorrMET_pt", "TMath::Sqrt(MET_uncert_dx * MET_uncert_dx + MET_uncert_dy * MET_uncert_dy)");
    }
    return df;
}

RNode JetEnergyCorrection(correction::CorrectionSet cset_jerc_2016preVFP, correction::CorrectionSet cset_jerc_2016postVFP, correction::CorrectionSet cset_jerc_2017, correction::CorrectionSet cset_jerc_2018, RNode df, std::string JEC_type, std::string variation) {
    auto eval_correction = [cset_jerc_2016preVFP, cset_jerc_2016postVFP, cset_jerc_2017, cset_jerc_2018, JEC_type, variation] (std::string year, RVec<float> pt, RVec<float> eta, RVec<float> var) {
        RVec<float> jec_factors;

        std::string JEC;

        if (var.size() == 0) {
            return var;
        }

        for (size_t i = 0; i < var.size(); i++) {
            if (year == "2016preVFP") {
                if (JEC_type.find("2016post") != std::string::npos || JEC_type.find("2017") != std::string::npos || JEC_type.find("2018") != std::string::npos) {
                    return var;
                }
                else {
                    JEC = std::string("Summer19UL16APV_V7_MC_") + JEC_type + std::string("_AK4PFchs");
                    if (variation == "up") {
                        jec_factors.push_back(var[i] * (1 + cset_jerc_2016preVFP.at(JEC)->evaluate({eta[i], pt[i]})));
                    }
                    else if (variation == "down") {
                        jec_factors.push_back(var[i] * (1 - cset_jerc_2016preVFP.at(JEC)->evaluate({eta[i], pt[i]})));
                    }
                    else {
                        return var;
                    }
                }
            }
            else if (year == "2016postVFP") {
                if (JEC_type.find("2016pre") != std::string::npos || JEC_type.find("2017") != std::string::npos || JEC_type.find("2018") != std::string::npos) {
                    return var;
                }
                else {
                    JEC = std::string("Summer19UL16_V7_MC_") + JEC_type + std::string("_AK4PFchs");
                    if (variation == "up") {
                        jec_factors.push_back(var[i] * (1 + cset_jerc_2016postVFP.at(JEC)->evaluate({eta[i], pt[i]})));
                    }
                    else if (variation == "down") {
                        jec_factors.push_back(var[i] * (1 - cset_jerc_2016postVFP.at(JEC)->evaluate({eta[i], pt[i]})));
                    }
                    else {
                        return var;
                    }
                }
            }
            else if (year == "2017") {
                if (JEC_type.find("2016") != std::string::npos || JEC_type.find("2018") != std::string::npos) {
                    return var;
                }
                else {
                    JEC = std::string("Summer19UL17_V5_MC_") + JEC_type + std::string("_AK4PFchs");
                    if (variation == "up") {
                        jec_factors.push_back(var[i] * (1 + cset_jerc_2017.at(JEC)->evaluate({eta[i], pt[i]})));
                    }
                    else if (variation == "down") {
                        jec_factors.push_back(var[i] * (1 - cset_jerc_2017.at(JEC)->evaluate({eta[i], pt[i]})));
                    }
                    else {
                        return var;
                    }
                }
            }
            else if (year == "2018") {
                if (JEC_type.find("2016") != std::string::npos || JEC_type.find("2017") != std::string::npos) {
                    return var;
                }
                else {
                    JEC = std::string("Summer19UL18_V5_MC_") + JEC_type + std::string("_AK4PFchs");
                    if (variation == "up") {
                        jec_factors.push_back(var[i] * (1 + cset_jerc_2018.at(JEC)->evaluate({eta[i], pt[i]})));
                    }
                    else if (variation == "down") {
                        jec_factors.push_back(var[i] * (1 - cset_jerc_2018.at(JEC)->evaluate({eta[i], pt[i]})));
                    }
                    else {
                        return var;
                    }
                }
            }
            else {
                return var;
            }
        }
        return jec_factors;
    };
    auto df_jetcorr = df.Redefine("CorrJet_pt", eval_correction, {"year", "Jet_pt", "Jet_eta", "Jet_pt"})
                        .Redefine("CorrJet_mass", eval_correction, {"year", "Jet_pt", "Jet_eta", "Jet_mass"})
                        .Redefine("CorrFatJet_pt", eval_correction, {"year", "FatJet_pt", "FatJet_eta", "FatJet_mass"})
                        .Redefine("CorrFatJet_mass", eval_correction, {"year", "FatJet_pt", "FatJet_eta", "FatJet_mass"});

    auto correctmet = [JEC_type](std::string year, RVec<float> corrjet_pt, RVec<float> jet_phi, RVec<float> jet_pt, float MET_pt, float MET_phi) {
        if (year == "2016preVFP") {
            if (JEC_type.find("2016post") != std::string::npos || JEC_type.find("2017") != std::string::npos || JEC_type.find("2018") != std::string::npos) {
                return MET_pt;
            }
        }
        else if (year == "2016postVFP") {
            if (JEC_type.find("2016pre") != std::string::npos || JEC_type.find("2017") != std::string::npos || JEC_type.find("2018") != std::string::npos) {
                return MET_pt;
            }
        }
        else if (year == "2017") {
            if (JEC_type.find("2016") != std::string::npos || JEC_type.find("2018") != std::string::npos) {
                return MET_pt;
            }
        }
        else if (year == "2018") {
            if (JEC_type.find("2016") != std::string::npos || JEC_type.find("2017") != std::string::npos) {
                return MET_pt;
            }
        }
        if (corrjet_pt.empty()) {
            return MET_pt;
        }
        float px = MET_pt * TMath::Cos(MET_phi);
        float py = MET_pt * TMath::Sin(MET_phi);
        for (size_t i = 0; i < corrjet_pt.size(); i++) {
            px += (corrjet_pt[i] - jet_pt[i]) * TMath::Cos(jet_phi[i]);
            py += (corrjet_pt[i] - jet_pt[i]) * TMath::Sin(jet_phi[i]);
        }
        return (float)TMath::Sqrt(px * px + py * py);
    };

    return df_jetcorr.Redefine("CorrMET_pt", correctmet, {"year", "CorrJet_pt", "Jet_phi", "Jet_pt", "MET_pt", "MET_phi"});
}

RNode JetEnergyResolution(correction::CorrectionSet cset_jerc_2016preVFP, correction::CorrectionSet cset_jerc_2016postVFP, correction::CorrectionSet cset_jerc_2017, correction::CorrectionSet cset_jerc_2018, correction::CorrectionSet cset_jer_smear, RNode df, std::string variation) {
    auto eval_correction = [cset_jerc_2016preVFP, cset_jerc_2016postVFP, cset_jerc_2017, cset_jerc_2018, cset_jer_smear, variation] (std::string year, RVec<float> pt, RVec<float> eta, RVec<int> genJet_idx, RVec<float> genJet_pt, float rho, unsigned long long event, RVec<float> var) {
        RVec<float> jer_factors;
        std::string vary;
        // for jer json
        if (variation == "nominal") {
            vary = "nom";
        }
        else {
            vary = variation;
        }
        float jer_sf;
        float jer;
        float genjetpt;
        if (var.size() == 0) {
            return var;
        }
        for (size_t i = 0; i < var.size(); i++) {
            genjetpt = genJet_idx[i] >= 0 ? genJet_pt[genJet_idx[i]] : -1;
            if (year == "2016preVFP") {
                jer = cset_jerc_2016preVFP.at("Summer20UL16APV_JRV3_MC_PtResolution_AK4PFchs")->evaluate({eta[i], pt[i], rho});
                jer_sf = cset_jerc_2016preVFP.at("Summer20UL16APV_JRV3_MC_ScaleFactor_AK4PFchs")->evaluate({eta[i], vary});
                jer_factors.push_back(var[i] * cset_jer_smear.at("JERSmear")->evaluate({pt[i], eta[i], genjetpt, rho, (int)event, jer, jer_sf}));
            }
            else if (year == "2016postVFP") {
                jer = cset_jerc_2016postVFP.at("Summer20UL16_JRV3_MC_PtResolution_AK4PFchs")->evaluate({eta[i], pt[i], rho});
                jer_sf = cset_jerc_2016postVFP.at("Summer20UL16_JRV3_MC_ScaleFactor_AK4PFchs")->evaluate({eta[i], vary});
                jer_factors.push_back(var[i] * cset_jer_smear.at("JERSmear")->evaluate({pt[i], eta[i], genjetpt, rho, (int)event, jer, jer_sf}));
            }
            else if (year == "2017") {
                jer = cset_jerc_2017.at("Summer19UL17_JRV2_MC_PtResolution_AK4PFchs")->evaluate({eta[i], pt[i], rho});
                jer_sf = cset_jerc_2017.at("Summer19UL17_JRV2_MC_ScaleFactor_AK4PFchs")->evaluate({eta[i], vary});
                jer_factors.push_back(var[i] * cset_jer_smear.at("JERSmear")->evaluate({pt[i], eta[i], genjetpt, rho, (int)event, jer, jer_sf}));
            }
            else if (year == "2018") {
                jer = cset_jerc_2018.at("Summer19UL18_JRV2_MC_PtResolution_AK4PFchs")->evaluate({eta[i], pt[i], rho});
                jer_sf = cset_jerc_2018.at("Summer19UL18_JRV2_MC_ScaleFactor_AK4PFchs")->evaluate({eta[i], vary});
                jer_factors.push_back(var[i] * cset_jer_smear.at("JERSmear")->evaluate({pt[i], eta[i], genjetpt, rho, (int)event, jer, jer_sf}));
            }
            else {
                return var;
            }
        }
        return jer_factors;
    };
    return df.Redefine("CorrJet_pt", eval_correction, {"year", "Jet_pt", "Jet_eta", "Jet_genJetIdx", "GenJet_pt", "fixedGridRhoFastjetAll", "event", "Jet_pt"})
            .Redefine("CorrJet_mass", eval_correction, {"year", "Jet_pt", "Jet_eta", "Jet_genJetIdx", "GenJet_pt", "fixedGridRhoFastjetAll", "event", "Jet_mass"});
}   

} // end namespace Run2