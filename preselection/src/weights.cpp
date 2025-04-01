#include "weights.h"


/*
############################################
GOLDEN JSON
############################################
*/

RNode applyGoldenJSONWeight(lumiMask golden, RNode df){
    auto goldenjson = [&golden](unsigned int &run, unsigned int &luminosityBlock){return golden.accept(run, luminosityBlock);};
    return df.Define("goldenJSON", goldenjson, {"run", "luminosityBlock"});
}

/*
############################################
PILEUP SFs
############################################
*/

RNode applyPileupScaleFactors(std::unordered_map<std::string, correction::CorrectionSet> cset_pileup, std::unordered_map<std::string, std::string> year_map, RNode df) {
    auto eval_correction = [&cset_pileup, &year_map] (std::string year, float ntrueint) {
        RVec<double> pileup_weights;
        auto correctionset = cset_pileup.at(year).at(year_map.at(year));
        pileup_weights.push_back(correctionset->evaluate({ntrueint, "nominal"}));
        pileup_weights.push_back(correctionset->evaluate({ntrueint, "up"}));
        pileup_weights.push_back(correctionset->evaluate({ntrueint, "down"}));
        return pileup_weights;
    };
    return df.Define("PILEUP_weight", eval_correction, {"sample_year", "Pileup_nTrueInt"});
}

/*
############################################
MUON SFs
############################################
*/

RNode applyMuonIDScaleFactors(std::unordered_map<std::string, correction::CorrectionSet> cset_muon, std::unordered_map<std::string, std::string> year_map, RNode df) {
    auto eval_correction = [&cset_muon, &year_map] (std::string year, const RVec<float> eta, const RVec<float> pt) {
        RVec<double> muon_sf_weights = {1., 1., 1.};
        if (eta.size() == 0) {
            return muon_sf_weights;
        }
        auto correctionset = cset_muon.at(year).at(year_map.at(year));
        for (size_t i = 0; i < eta.size(); i++) {
            muon_sf_weights[0] *= correctionset->evaluate({abs(eta[i]), pt[i], "nominal"});
            muon_sf_weights[1] *= correctionset->evaluate({abs(eta[i]), pt[i], "systup"});
            muon_sf_weights[2] *= correctionset->evaluate({abs(eta[i]), pt[i], "systdown"});
        }
        return muon_sf_weights;
    };
    return df.Define("MUONID_scalefactors", eval_correction, {"sample_year", "muon_eta", "muon_pt"});
}

RNode applyMuonRecoScaleFactors(std::unordered_map<std::string, correction::CorrectionSet> cset_muon, std::unordered_map<std::string, std::string> year_map, RNode df) {
    auto eval_correction = [&cset_muon, &year_map] (std::string year, const RVec<float> eta, const RVec<float> pt) {
        RVec<double> muon_sf_weights = {1., 1., 1.};
        if (eta.size() == 0) {
            return muon_sf_weights;
        }
        auto correctionset = cset_muon.at(year).at(year_map.at(year));
        for (size_t i = 0; i < eta.size(); i++) {
            muon_sf_weights[0] *= correctionset->evaluate({abs(eta[i]), pt[i], "nominal"});
            muon_sf_weights[1] *= correctionset->evaluate({abs(eta[i]), pt[i], "systup"});
            muon_sf_weights[2] *= correctionset->evaluate({abs(eta[i]), pt[i], "systdown"});
        }
        return muon_sf_weights;
    };
    return df.Define("MUONRECO_scalefactors", eval_correction, {"sample_year", "muon_eta", "muon_pt"});
}

RNode applyMuonTriggerScaleFactors(std::unordered_map<std::string, correction::CorrectionSet> cset_muon, std::unordered_map<std::string, std::string> year_map, RNode df) {
    auto eval_correction = [&cset_muon, &year_map] (std::string year, const RVec<float> eta, const RVec<float> pt) {
        RVec<double> muon_sf_weights = {1., 1., 1.};
        if (eta.size() == 0) {
            return muon_sf_weights;
        }
        auto correctionset = cset_muon.at(year).at(year_map.at(year));
        for (size_t i = 0; i < eta.size(); i++) {
            muon_sf_weights[0] *= correctionset->evaluate({abs(eta[i]), pt[i], "nominal"});
            muon_sf_weights[1] *= correctionset->evaluate({abs(eta[i]), pt[i], "systup"});
            muon_sf_weights[2] *= correctionset->evaluate({abs(eta[i]), pt[i], "systdown"});
        }
        return muon_sf_weights;
    };
    return df.Define("MUONTRIGGER_scalefactors", eval_correction, {"sample_year", "muon_eta", "muon_pt"});
}

/*
############################################
ELECTRON SFs
############################################
*/

RNode applyElectronIDScaleFactors(std::unordered_map<std::string, correction::CorrectionSet> cset_electron, std::unordered_map<std::string, std::string> year_map, RNode df) {
    auto eval_correction = [&cset_electron, &year_map] (std::string year, const RVec<float> eta, const RVec<float> pt) {
        RVec<double> electron_sf_weights = {1., 1., 1.};
        if (eta.size() == 0) {
            return electron_sf_weights;
        }
        auto correctionset = cset_electron.at(year).at(year_map.at(year));
        for (size_t i = 0; i < eta.size(); i++) {
            electron_sf_weights[0] *= correctionset->evaluate({year, "sf", "Tight", eta[i], pt[i]});
            electron_sf_weights[1] *= correctionset->evaluate({year, "sfup", "Tight", eta[i], pt[i]});
            electron_sf_weights[2] *= correctionset->evaluate({year, "sfdown", "Tight", eta[i], pt[i]});
        }
        return electron_sf_weights;
    };
    return df.Define("ELECTRONID_scalefactors", eval_correction, {"sample_year", "electron_sceta", "electron_pt"});
}

RNode applyElectronRecoScaleFactors(std::unordered_map<std::string, correction::CorrectionSet> cset_electron, std::unordered_map<std::string, std::string> year_map, RNode df) {
    auto eval_correction = [&cset_electron, &year_map] (std::string year, const RVec<float> eta, const RVec<float> pt) {
        RVec<double> electron_sf_weights = {1., 1., 1.};
        if (eta.size() == 0) {
            return electron_sf_weights;
        }
        auto correctionset = cset_electron.at(year).at(year_map.at(year));
        for (size_t i = 0; i < eta.size(); i++) {
            if (pt[i] >= 20 && pt[i] < 75) {
                electron_sf_weights[0] *= correctionset->evaluate({year, "sf", "Reco20to75", eta[i], pt[i]});
                electron_sf_weights[1] *= correctionset->evaluate({year, "sfup", "Reco20to75", eta[i], pt[i]});
                electron_sf_weights[2] *= correctionset->evaluate({year, "sfdown", "Reco20to75", eta[i], pt[i]});
            }
            else if (pt[i] >= 75) {
                electron_sf_weights[0] *= correctionset->evaluate({year, "sf", "RecoAbove75", eta[i], pt[i]});
                electron_sf_weights[1] *= correctionset->evaluate({year, "sfup", "RecoAbove75", eta[i], pt[i]});
                electron_sf_weights[2] *= correctionset->evaluate({year, "sfdown", "RecoAbove75", eta[i], pt[i]});
            }
        }
        return electron_sf_weights;
    };
    return df.Define("ELECTRONRECO_scalefactors", eval_correction, {"sample_year", "electron_sceta", "electron_pt"});
}

RNode applyElectronTriggerScaleFactors(std::unordered_map<std::string, correction::CorrectionSet> cset_electron, std::unordered_map<std::string, std::string> year_map, RNode df) {
    auto eval_correction = [&cset_electron, &year_map] (std::string year, const RVec<float> eta, const RVec<float> pt) {
        RVec<double> electron_sf_weights = {1., 1., 1.};
        if (eta.size() == 0) {
            return electron_sf_weights;
        }
        auto correctionset = cset_electron.at(year).at(year_map.at(year));
        for (size_t i = 0; i < eta.size(); i++) {
            electron_sf_weights[0] *= correctionset->evaluate({year, "sf", "HLT_SF_Ele30_TightID", eta[i], pt[i]});
            electron_sf_weights[1] *= correctionset->evaluate({year, "sfup", "HLT_SF_Ele30_TightID", eta[i], pt[i]});
            electron_sf_weights[2] *= correctionset->evaluate({year, "sfdown", "HLT_SF_Ele30_TightID", eta[i], pt[i]});
        }
        return electron_sf_weights;
    };
    return df.Define("ELECTRONTRIGGER_scalefactors", eval_correction, {"sample_year", "electron_sceta", "electron_pt"});
}

/*
############################################
OTHER SFs
############################################
*/

RNode applyL1PreFiringReweighting(RNode df){
    auto eval_correction = [] (float L1prefire, float L1prefireup, float L1prefiredown) {
        return RVec<float>{L1prefire, L1prefireup, L1prefiredown};
    };
    return df.Define("L1PREFIRING_weight", eval_correction, {"L1PreFiringWeight_Nom", "L1PreFiringWeight_Up", "L1PreFiringWeight_Dn"});
}

RNode applyPSWeight_FSR(RNode df) {
    auto eval_correction = [] (const RVec<float> PSWeight) {
        return RVec<float>{1., PSWeight[1], PSWeight[3]};
    };
    return df.Define("PSWeight_FSR", eval_correction, {"PSWeight"});
}

RNode applyPSWeight_ISR(RNode df) {
    auto eval_correction = [] (const RVec<float> PSWeight) {
        return RVec<float>{1., PSWeight[0], PSWeight[2]};
    };
    return df.Define("PSWeight_ISR", eval_correction, {"PSWeight"});
}

RNode applyLHEScaleWeight_muF(RNode df) {
    auto eval_correction = [] (const RVec<float> LHEScaleWeight) {
        return RVec<float>{1., LHEScaleWeight[5], LHEScaleWeight[3]};
    };
    return df.Define("LHEScaleWeight_muF", eval_correction, {"LHEScaleWeight"});
}

RNode applyLHEScaleWeight_muR(RNode df) {
    auto eval_correction = [] (const RVec<float> LHEScaleWeight) {
        return RVec<float>{1., LHEScaleWeight[7], LHEScaleWeight[1]};
    };
    return df.Define("LHEScaleWeight_muR", eval_correction, {"LHEScaleWeight"});
}

RNode applyLHEWeights_pdf(RNode df) {
    auto eval_correction = [] (const RVec<float> LHEWeights) {
        RVec<float> PDFWeights = {1., 1., 1.};
        float PDFUncValue = 0.0;
        for (const auto& weight : LHEWeights) {
            PDFUncValue += (weight - 1) * (weight - 1);
        }
        PDFUncValue = sqrt(PDFUncValue);

        PDFWeights[1] = (1 + PDFUncValue);
        PDFWeights[2] = (1 - PDFUncValue);
        return PDFWeights;
    };  
    return df.Define("LHEWeights_pdf", eval_correction, {"LHEPdfWeight"});
}

RNode applyDataWeights(RNode df_) {
    auto df = applyGoldenJSONWeight(LumiMask, df_);
    return df.Define("weight", "goldenJSON");
}

RNode applyMCWeights(RNode df_) {
    auto df = applyPileupScaleFactors(pileupScaleFactors, pileupScaleFactors_yearmap, df_);
    df = applyMuonIDScaleFactors(muonScaleFactors, muonIDScaleFactors_yearmap, df);
    df = applyMuonRecoScaleFactors(muonScaleFactors, muonRecoScaleFactors_yearmap, df);
    df = applyMuonTriggerScaleFactors(muonScaleFactors, muonTriggerScaleFactors_yearmap, df);

    df = applyElectronIDScaleFactors(electronScaleFactors, electronScaleFactors_yearmap, df);
    df = applyElectronRecoScaleFactors(electronScaleFactors, electronScaleFactors_yearmap, df);
    df = applyElectronTriggerScaleFactors(electronTriggerScaleFactors, electronTriggerScaleFactors_yearmap, df);

    return df.Define("weight", 
        "PILEUP_weight[0] * "
        "MUONID_scalefactors[0] * "
        "MUONRECO_scalefactors[0] * "
        "MUONTRIGGER_scalefactors[0] * "
        "ELECTRONID_scalefactors[0] * "
        "ELECTRONRECO_scalefactors[0] * "
        "ELECTRONTRIGGER_scalefactors[0] * "
        "genWeight * "
        "xsec_weight");
}