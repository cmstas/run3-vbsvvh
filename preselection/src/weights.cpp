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
PARTICLE NET SFs
############################################
*/

// RNode PNET_W_ScaleFactors_2016preVFP(correction::CorrectionSet cset_pnet_w, RNode df) {
//     auto eval_correction = [cset_pnet_w] (std::string year, float pt) {
//         RVec<double> pnet_w = {1., 1., 1.};
//         if (year != "2016preVFP") {
//             return pnet_w;
//         }
//         pnet_w[0] *= cset_pnet_w.at("PNET_W")->evaluate({pt, year, "nominal"});
//         pnet_w[1] *= cset_pnet_w.at("PNET_W")->evaluate({pt, year, "up"});
//         pnet_w[2] *= cset_pnet_w.at("PNET_W")->evaluate({pt, year, "down"});
//         return pnet_w;
//     };
//     return df.Define("particlenet_w_weight_2016preVFP", eval_correction, {"sample_year", "GW_pt"});
// }

// RNode PNET_W_ScaleFactors_2016postVFP(correction::CorrectionSet cset_pnet_w, RNode df) {
//     auto eval_correction = [cset_pnet_w] (std::string year, float pt) {
//         RVec<double> pnet_w = {1., 1., 1.};
//         if (year != "2016postVFP") {
//             return pnet_w;
//         }
//         pnet_w[0] *= cset_pnet_w.at("PNET_W")->evaluate({pt, year, "nominal"});
//         pnet_w[1] *= cset_pnet_w.at("PNET_W")->evaluate({pt, year, "up"});
//         pnet_w[2] *= cset_pnet_w.at("PNET_W")->evaluate({pt, year, "down"});
//         return pnet_w;
//     };
//     return df.Define("particlenet_w_weight_2016postVFP", eval_correction, {"sample_year", "GW_pt"});
// }

// RNode PNET_W_ScaleFactors_2017(correction::CorrectionSet cset_pnet_w, RNode df) {
//     auto eval_correction = [cset_pnet_w] (std::string year, float pt) {
//         RVec<double> pnet_w = {1., 1., 1.};
//         if (year != "2017") {
//             return pnet_w;
//         }
//         pnet_w[0] *= cset_pnet_w.at("PNET_W")->evaluate({pt, year, "nominal"});
//         pnet_w[1] *= cset_pnet_w.at("PNET_W")->evaluate({pt, year, "up"});
//         pnet_w[2] *= cset_pnet_w.at("PNET_W")->evaluate({pt, year, "down"});
//         return pnet_w;
//     };
//     return df.Define("particlenet_w_weight_2017", eval_correction, {"sample_year", "GW_pt"});
// }

// RNode PNET_W_ScaleFactors_2018(correction::CorrectionSet cset_pnet_w, RNode df) {
//     auto eval_correction = [cset_pnet_w] (std::string year, float pt) {
//         RVec<double> pnet_w = {1., 1., 1.};
//         if (year != "2018") {
//             return pnet_w;
//         }
//         pnet_w[0] *= cset_pnet_w.at("PNET_W")->evaluate({pt, year, "nominal"});
//         pnet_w[1] *= cset_pnet_w.at("PNET_W")->evaluate({pt, year, "up"});
//         pnet_w[2] *= cset_pnet_w.at("PNET_W")->evaluate({pt, year, "down"});
//         return pnet_w;
//     };
//     return df.Define("particlenet_w_weight_2018", eval_correction, {"sample_year", "GW_pt"});
// }

// RNode PNET_H_ScaleFactors_2016preVFP(correction::CorrectionSet cset_pnet_h, RNode df) {
//     auto eval_correction = [cset_pnet_h] (std::string year, float pt) {
//         RVec<double> pnet_h = {1., 1., 1.};
//         if (year != "2016preVFP") {
//             return pnet_h;
//         }
//         pnet_h[0] *= cset_pnet_h.at("PNET_H")->evaluate({pt, year, "nominal"});
//         pnet_h[1] *= cset_pnet_h.at("PNET_H")->evaluate({pt, year, "up"});
//         pnet_h[2] *= cset_pnet_h.at("PNET_H")->evaluate({pt, year, "down"});
//         return pnet_h;
//     };
//     return df.Define("particlenet_h_weight_2016preVFP", eval_correction, {"sample_year", "GHiggs_pt"});
// }

// RNode PNET_H_ScaleFactors_2016postVFP(correction::CorrectionSet cset_pnet_h, RNode df) {
//     auto eval_correction = [cset_pnet_h] (std::string year, float pt) {
//         RVec<double> pnet_h = {1., 1., 1.};
//         if (year != "2016postVFP") {
//             return pnet_h;
//         }
//         pnet_h[0] *= cset_pnet_h.at("PNET_H")->evaluate({pt, year, "nominal"});
//         pnet_h[1] *= cset_pnet_h.at("PNET_H")->evaluate({pt, year, "up"});
//         pnet_h[2] *= cset_pnet_h.at("PNET_H")->evaluate({pt, year, "down"});
//         return pnet_h;
//     };
//     return df.Define("particlenet_h_weight_2016postVFP", eval_correction, {"sample_year", "GHiggs_pt"});
// }

// RNode PNET_H_ScaleFactors_2017(correction::CorrectionSet cset_pnet_h, RNode df) {
//     auto eval_correction = [cset_pnet_h] (std::string year, float pt) {
//         RVec<double> pnet_h = {1., 1., 1.};
//         if (year != "2017") {
//             return pnet_h;
//         }
//         pnet_h[0] *= cset_pnet_h.at("PNET_H")->evaluate({pt, year, "nominal"});
//         pnet_h[1] *= cset_pnet_h.at("PNET_H")->evaluate({pt, year, "up"});
//         pnet_h[2] *= cset_pnet_h.at("PNET_H")->evaluate({pt, year, "down"});
//         return pnet_h;
//     };
//     return df.Define("particlenet_h_weight_2017", eval_correction, {"sample_year", "GHiggs_pt"});
// }

// RNode PNET_H_ScaleFactors_2018(correction::CorrectionSet cset_pnet_h, RNode df) {
//     auto eval_correction = [cset_pnet_h] (std::string year, float pt) {
//         RVec<double> pnet_h = {1., 1., 1.};
//         if (year != "2018") {
//             return pnet_h;
//         }
//         pnet_h[0] *= cset_pnet_h.at("PNET_H")->evaluate({pt, year, "nominal"});
//         pnet_h[1] *= cset_pnet_h.at("PNET_H")->evaluate({pt, year, "up"});
//         pnet_h[2] *= cset_pnet_h.at("PNET_H")->evaluate({pt, year, "down"});
//         return pnet_h;
//     };
//     return df.Define("particlenet_h_weight_2018", eval_correction, {"sample_year", "GHiggs_pt"});
// }

/*
############################################
B-TAGGING SFs
############################################
*/

// RNode bTaggingScaleFactors_HF(correction::CorrectionSet cset_btag_2016preVFP, correction::CorrectionSet cset_btag_2016postVFP, correction::CorrectionSet cset_btag_2017, correction::CorrectionSet cset_btag_2018, correction::CorrectionSet cset_btag_eff, RNode df) {
//     auto eval_correction = [cset_btag_2016preVFP, cset_btag_2016postVFP, cset_btag_2017, cset_btag_2018, cset_btag_eff] (std::string year, const RVec<float> eta, const RVec<float> pt, const RVec<int> jetflavor) {
//         RVec<double> btag_sf_weights = {1., 1., 1.};
//         if (eta.size() == 0) {
//             return btag_sf_weights;
//         }
//         float num = 1., num_up = 1., num_down = 1., den = 1.;
//         for (size_t i = 0; i < eta.size(); i++) {
//             float btag_sf_tight = 0, btag_sf_loose = 0, btag_sf_tight_up = 0, btag_sf_loose_up = 0, btag_sf_tight_down = 0, btag_sf_loose_down = 0, btag_eff_tight = 0, btag_eff_loose = 0;
//             if (year == "2016preVFP") {
//                 if (jetflavor[i] == 5) {
//                     btag_sf_tight = cset_btag_2016preVFP.at("deepCSV_comb")->evaluate({"central", "T", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_loose = cset_btag_2016preVFP.at("deepCSV_comb")->evaluate({"central", "L", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_tight_up = cset_btag_2016preVFP.at("deepCSV_comb")->evaluate({"up_uncorrelated", "T", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_loose_up = cset_btag_2016preVFP.at("deepCSV_comb")->evaluate({"up_uncorrelated", "L", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_tight_down = cset_btag_2016preVFP.at("deepCSV_comb")->evaluate({"down_uncorrelated", "T", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_loose_down = cset_btag_2016preVFP.at("deepCSV_comb")->evaluate({"down_uncorrelated", "L", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_eff_tight = cset_btag_eff.at("btag_2016preVFP")->evaluate({"B", "T", pt[i], eta[i]});
//                     btag_eff_loose = cset_btag_eff.at("btag_2016preVFP")->evaluate({"B", "L", pt[i], eta[i]});
//                 }
//                 else if (jetflavor[i] == 4) {
//                     btag_sf_tight = cset_btag_2016preVFP.at("deepCSV_comb")->evaluate({"central", "T", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_loose = cset_btag_2016preVFP.at("deepCSV_comb")->evaluate({"central", "L", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_tight_up = cset_btag_2016preVFP.at("deepCSV_comb")->evaluate({"up_uncorrelated", "T", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_loose_up = cset_btag_2016preVFP.at("deepCSV_comb")->evaluate({"up_uncorrelated", "L", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_tight_down = cset_btag_2016preVFP.at("deepCSV_comb")->evaluate({"down_uncorrelated", "T", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_loose_down = cset_btag_2016preVFP.at("deepCSV_comb")->evaluate({"down_uncorrelated", "L", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_eff_tight = cset_btag_eff.at("btag_2016preVFP")->evaluate({"C", "T", pt[i], eta[i]});
//                     btag_eff_loose = cset_btag_eff.at("btag_2016preVFP")->evaluate({"C", "L", pt[i], eta[i]});
//                 }
//                 num *= (btag_sf_tight * btag_eff_tight) * (btag_sf_loose * btag_eff_loose - btag_sf_tight * btag_eff_tight) * (1. - btag_sf_loose * btag_eff_loose);
//                 num_up *= (btag_sf_tight_up * btag_eff_tight) * (btag_sf_loose_up * btag_eff_loose - btag_sf_tight_up * btag_eff_tight) * (1. - btag_sf_loose_up * btag_eff_loose);
//                 num_down *= (btag_sf_tight_down * btag_eff_tight) * (btag_sf_loose_down * btag_eff_loose - btag_sf_tight_down * btag_eff_tight) * (1. - btag_sf_loose_down * btag_eff_loose);
//                 den *= (btag_eff_tight) * (btag_eff_loose - btag_eff_tight) * (1. - btag_eff_loose);
//             }
//             if (year == "2016postVFP") {
//                 if (jetflavor[i] == 5) {
//                     btag_sf_tight = cset_btag_2016postVFP.at("deepCSV_comb")->evaluate({"central", "T", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_loose = cset_btag_2016postVFP.at("deepCSV_comb")->evaluate({"central", "L", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_tight_up = cset_btag_2016postVFP.at("deepCSV_comb")->evaluate({"up_uncorrelated", "T", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_loose_up = cset_btag_2016postVFP.at("deepCSV_comb")->evaluate({"up_uncorrelated", "L", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_tight_down = cset_btag_2016postVFP.at("deepCSV_comb")->evaluate({"down_uncorrelated", "T", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_loose_down = cset_btag_2016postVFP.at("deepCSV_comb")->evaluate({"down_uncorrelated", "L", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_eff_tight = cset_btag_eff.at("btag_2016postVFP")->evaluate({"B", "T", pt[i], eta[i]});
//                     btag_eff_loose = cset_btag_eff.at("btag_2016postVFP")->evaluate({"B", "L", pt[i], eta[i]});
//                 }
//                 else if (jetflavor[i] == 4) {
//                     btag_sf_tight = cset_btag_2016postVFP.at("deepCSV_comb")->evaluate({"central", "T", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_loose = cset_btag_2016postVFP.at("deepCSV_comb")->evaluate({"central", "L", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_tight_up = cset_btag_2016postVFP.at("deepCSV_comb")->evaluate({"up_uncorrelated", "T", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_loose_up = cset_btag_2016postVFP.at("deepCSV_comb")->evaluate({"up_uncorrelated", "L", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_tight_down = cset_btag_2016postVFP.at("deepCSV_comb")->evaluate({"down_uncorrelated", "T", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_loose_down = cset_btag_2016postVFP.at("deepCSV_comb")->evaluate({"down_uncorrelated", "L", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_eff_tight = cset_btag_eff.at("btag_2016postVFP")->evaluate({"C", "T", pt[i], eta[i]});
//                     btag_eff_loose = cset_btag_eff.at("btag_2016postVFP")->evaluate({"C", "L", pt[i], eta[i]});
//                 }
//                 num *= (btag_sf_tight * btag_eff_tight) * (btag_sf_loose * btag_eff_loose - btag_sf_tight * btag_eff_tight) * (1. - btag_sf_loose * btag_eff_loose);
//                 num_up *= (btag_sf_tight_up * btag_eff_tight) * (btag_sf_loose_up * btag_eff_loose - btag_sf_tight_up * btag_eff_tight) * (1. - btag_sf_loose_up * btag_eff_loose);
//                 num_down *= (btag_sf_tight_down * btag_eff_tight) * (btag_sf_loose_down * btag_eff_loose - btag_sf_tight_down * btag_eff_tight) * (1. - btag_sf_loose_down * btag_eff_loose);
//                 den *= (btag_eff_tight) * (btag_eff_loose - btag_eff_tight) * (1. - btag_eff_loose);
//             }
//             if (year == "2017") {
//                 if (jetflavor[i] == 5) {
//                     btag_sf_tight = cset_btag_2017.at("deepCSV_comb")->evaluate({"central", "T", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_loose = cset_btag_2017.at("deepCSV_comb")->evaluate({"central", "L", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_tight_up = cset_btag_2017.at("deepCSV_comb")->evaluate({"up_uncorrelated", "T", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_loose_up = cset_btag_2017.at("deepCSV_comb")->evaluate({"up_uncorrelated", "L", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_tight_down = cset_btag_2017.at("deepCSV_comb")->evaluate({"down_uncorrelated", "T", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_loose_down = cset_btag_2017.at("deepCSV_comb")->evaluate({"down_uncorrelated", "L", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_eff_tight = cset_btag_eff.at("btag_2017")->evaluate({"B", "T", pt[i], eta[i]});
//                     btag_eff_loose = cset_btag_eff.at("btag_2017")->evaluate({"B", "L", pt[i], eta[i]});
//                 }
//                 else if (jetflavor[i] == 4) {
//                     btag_sf_tight = cset_btag_2017.at("deepCSV_comb")->evaluate({"central", "T", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_loose = cset_btag_2017.at("deepCSV_comb")->evaluate({"central", "L", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_tight_up = cset_btag_2017.at("deepCSV_comb")->evaluate({"up_uncorrelated", "T", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_loose_up = cset_btag_2017.at("deepCSV_comb")->evaluate({"up_uncorrelated", "L", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_tight_down = cset_btag_2017.at("deepCSV_comb")->evaluate({"down_uncorrelated", "T", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_loose_down = cset_btag_2017.at("deepCSV_comb")->evaluate({"down_uncorrelated", "L", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_eff_tight = cset_btag_eff.at("btag_2017")->evaluate({"C", "T", pt[i], eta[i]});
//                     btag_eff_loose = cset_btag_eff.at("btag_2017")->evaluate({"C", "L", pt[i], eta[i]});
//                 }
//                 num *= (btag_sf_tight * btag_eff_tight) * (btag_sf_loose * btag_eff_loose - btag_sf_tight * btag_eff_tight) * (1. - btag_sf_loose * btag_eff_loose);
//                 num_up *= (btag_sf_tight_up * btag_eff_tight) * (btag_sf_loose_up * btag_eff_loose - btag_sf_tight_up * btag_eff_tight) * (1. - btag_sf_loose_up * btag_eff_loose);
//                 num_down *= (btag_sf_tight_down * btag_eff_tight) * (btag_sf_loose_down * btag_eff_loose - btag_sf_tight_down * btag_eff_tight) * (1. - btag_sf_loose_down * btag_eff_loose);
//                 den *= (btag_eff_tight) * (btag_eff_loose - btag_eff_tight) * (1. - btag_eff_loose);
//             }
//             if (year == "2018") {
//                 if (jetflavor[i] == 5) {
//                     btag_sf_tight = cset_btag_2018.at("deepCSV_comb")->evaluate({"central", "T", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_loose = cset_btag_2018.at("deepCSV_comb")->evaluate({"central", "L", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_tight_up = cset_btag_2018.at("deepCSV_comb")->evaluate({"up_uncorrelated", "T", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_loose_up = cset_btag_2018.at("deepCSV_comb")->evaluate({"up_uncorrelated", "L", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_tight_down = cset_btag_2018.at("deepCSV_comb")->evaluate({"down_uncorrelated", "T", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_loose_down = cset_btag_2018.at("deepCSV_comb")->evaluate({"down_uncorrelated", "L", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_eff_tight = cset_btag_eff.at("btag_2018")->evaluate({"B", "T", pt[i], eta[i]});
//                     btag_eff_loose = cset_btag_eff.at("btag_2018")->evaluate({"B", "L", pt[i], eta[i]});
//                 }
//                 else if (jetflavor[i] == 4) {
//                     btag_sf_tight = cset_btag_2018.at("deepCSV_comb")->evaluate({"central", "T", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_loose = cset_btag_2018.at("deepCSV_comb")->evaluate({"central", "L", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_tight_up = cset_btag_2018.at("deepCSV_comb")->evaluate({"up_uncorrelated", "T", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_loose_up = cset_btag_2018.at("deepCSV_comb")->evaluate({"up_uncorrelated", "L", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_tight_down = cset_btag_2018.at("deepCSV_comb")->evaluate({"down_uncorrelated", "T", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_loose_down = cset_btag_2018.at("deepCSV_comb")->evaluate({"down_uncorrelated", "L", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_eff_tight = cset_btag_eff.at("btag_2018")->evaluate({"C", "T", pt[i], eta[i]});
//                     btag_eff_loose = cset_btag_eff.at("btag_2018")->evaluate({"C", "L", pt[i], eta[i]});
//                 }
//                 num *= (btag_sf_tight * btag_eff_tight) * (btag_sf_loose * btag_eff_loose - btag_sf_tight * btag_eff_tight) * (1. - btag_sf_loose * btag_eff_loose);
//                 num_up *= (btag_sf_tight_up * btag_eff_tight) * (btag_sf_loose_up * btag_eff_loose - btag_sf_tight_up * btag_eff_tight) * (1. - btag_sf_loose_up * btag_eff_loose);
//                 num_down *= (btag_sf_tight_down * btag_eff_tight) * (btag_sf_loose_down * btag_eff_loose - btag_sf_tight_down * btag_eff_tight) * (1. - btag_sf_loose_down * btag_eff_loose);
//                 den *= (btag_eff_tight) * (btag_eff_loose - btag_eff_tight) * (1. - btag_eff_loose);
//             }
//         }
//         if (den == 0.) {
//             btag_sf_weights[0] = 1.;
//             btag_sf_weights[1] = 1.;
//             btag_sf_weights[2] = 1.;
//         }
//         else {
//             btag_sf_weights[0] = num / den;
//             btag_sf_weights[1] = num_up / den;
//             btag_sf_weights[2] = num_down / den;
//         }
//         return btag_sf_weights;
//     };
//     return df.Define("btagging_scale_factors_HF", eval_correction, {"sample_year", "LnTBJet_eta", "LnTBJet_pt", "LnTBJet_hadronFlavour"});
// }

// RNode bTaggingScaleFactors_LF(correction::CorrectionSet cset_btag_2016preVFP, correction::CorrectionSet cset_btag_2016postVFP, correction::CorrectionSet cset_btag_2017, correction::CorrectionSet cset_btag_2018, correction::CorrectionSet cset_btag_eff, RNode df) {
//     auto eval_correction = [cset_btag_2016preVFP, cset_btag_2016postVFP, cset_btag_2017, cset_btag_2018, cset_btag_eff] (std::string year, const RVec<float> eta, const RVec<float> pt, const RVec<int> jetflavor) {
//         RVec<double> btag_sf_weights = {1., 1., 1.};
//         if (eta.size() == 0) {
//             return btag_sf_weights;
//         }
//         float num = 1., num_up = 1., num_down = 1., den = 1.;
//         for (size_t i = 0; i < eta.size(); i++) {
//             float btag_sf_tight = 0, btag_sf_loose = 0, btag_sf_tight_up = 0, btag_sf_loose_up = 0, btag_sf_tight_down = 0, btag_sf_loose_down = 0, btag_eff_tight = 0, btag_eff_loose = 0;
//             if (year == "2016preVFP") {
//                 if (jetflavor[i] == 0) {
//                     btag_sf_tight = cset_btag_2016preVFP.at("deepCSV_incl")->evaluate({"central", "T", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_loose = cset_btag_2016preVFP.at("deepCSV_incl")->evaluate({"central", "L", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_tight_up = cset_btag_2016preVFP.at("deepCSV_incl")->evaluate({"up_uncorrelated", "T", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_loose_up = cset_btag_2016preVFP.at("deepCSV_incl")->evaluate({"up_uncorrelated", "L", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_tight_down = cset_btag_2016preVFP.at("deepCSV_incl")->evaluate({"down_uncorrelated", "T", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_loose_down = cset_btag_2016preVFP.at("deepCSV_incl")->evaluate({"down_uncorrelated", "L", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_eff_tight = cset_btag_eff.at("btag_2016preVFP")->evaluate({"L", "T", pt[i], eta[i]});
//                     btag_eff_loose = cset_btag_eff.at("btag_2016preVFP")->evaluate({"L", "L", pt[i], eta[i]});
//                 }
//                 num *= (btag_sf_tight * btag_eff_tight) * (btag_sf_loose * btag_eff_loose - btag_sf_tight * btag_eff_tight) * (1. - btag_sf_loose * btag_eff_loose);
//                 num_up *= (btag_sf_tight_up * btag_eff_tight) * (btag_sf_loose_up * btag_eff_loose - btag_sf_tight_up * btag_eff_tight) * (1. - btag_sf_loose_up * btag_eff_loose);
//                 num_down *= (btag_sf_tight_down * btag_eff_tight) * (btag_sf_loose_down * btag_eff_loose - btag_sf_tight_down * btag_eff_tight) * (1. - btag_sf_loose_down * btag_eff_loose);
//                 den *= (btag_eff_tight) * (btag_eff_loose - btag_eff_tight) * (1. - btag_eff_loose);
//             }
//             if (year == "2016postVFP") {
//                 if (jetflavor[i] == 0) {
//                     btag_sf_tight = cset_btag_2016postVFP.at("deepCSV_incl")->evaluate({"central", "T", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_loose = cset_btag_2016postVFP.at("deepCSV_incl")->evaluate({"central", "L", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_tight_up = cset_btag_2016postVFP.at("deepCSV_incl")->evaluate({"up_uncorrelated", "T", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_loose_up = cset_btag_2016postVFP.at("deepCSV_incl")->evaluate({"up_uncorrelated", "L", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_tight_down = cset_btag_2016postVFP.at("deepCSV_incl")->evaluate({"down_uncorrelated", "T", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_loose_down = cset_btag_2016postVFP.at("deepCSV_incl")->evaluate({"down_uncorrelated", "L", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_eff_tight = cset_btag_eff.at("btag_2016postVFP")->evaluate({"L", "T", pt[i], eta[i]});
//                     btag_eff_loose = cset_btag_eff.at("btag_2016postVFP")->evaluate({"L", "L", pt[i], eta[i]});
//                 }
//                 num *= (btag_sf_tight * btag_eff_tight) * (btag_sf_loose * btag_eff_loose - btag_sf_tight * btag_eff_tight) * (1. - btag_sf_loose * btag_eff_loose);
//                 num_up *= (btag_sf_tight_up * btag_eff_tight) * (btag_sf_loose_up * btag_eff_loose - btag_sf_tight_up * btag_eff_tight) * (1. - btag_sf_loose_up * btag_eff_loose);
//                 num_down *= (btag_sf_tight_down * btag_eff_tight) * (btag_sf_loose_down * btag_eff_loose - btag_sf_tight_down * btag_eff_tight) * (1. - btag_sf_loose_down * btag_eff_loose);
//                 den *= (btag_eff_tight) * (btag_eff_loose - btag_eff_tight) * (1. - btag_eff_loose);
//             }
//             if (year == "2017") {
//                 if (jetflavor[i] == 0) {
//                     btag_sf_tight = cset_btag_2017.at("deepCSV_incl")->evaluate({"central", "T", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_loose = cset_btag_2017.at("deepCSV_incl")->evaluate({"central", "L", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_tight_up = cset_btag_2017.at("deepCSV_incl")->evaluate({"up_uncorrelated", "T", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_loose_up = cset_btag_2017.at("deepCSV_incl")->evaluate({"up_uncorrelated", "L", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_tight_down = cset_btag_2017.at("deepCSV_incl")->evaluate({"down_uncorrelated", "T", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_loose_down = cset_btag_2017.at("deepCSV_incl")->evaluate({"down_uncorrelated", "L", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_eff_tight = cset_btag_eff.at("btag_2017")->evaluate({"L", "T", pt[i], eta[i]});
//                     btag_eff_loose = cset_btag_eff.at("btag_2017")->evaluate({"L", "L", pt[i], eta[i]});
//                 }
//                 num *= (btag_sf_tight * btag_eff_tight) * (btag_sf_loose * btag_eff_loose - btag_sf_tight * btag_eff_tight) * (1. - btag_sf_loose * btag_eff_loose);
//                 num_up *= (btag_sf_tight_up * btag_eff_tight) * (btag_sf_loose_up * btag_eff_loose - btag_sf_tight_up * btag_eff_tight) * (1. - btag_sf_loose_up * btag_eff_loose);
//                 num_down *= (btag_sf_tight_down * btag_eff_tight) * (btag_sf_loose_down * btag_eff_loose - btag_sf_tight_down * btag_eff_tight) * (1. - btag_sf_loose_down * btag_eff_loose);
//                 den *= (btag_eff_tight) * (btag_eff_loose - btag_eff_tight) * (1. - btag_eff_loose);
//             }
//             if (year == "2018") {
//                 if (jetflavor[i] == 0) {
//                     btag_sf_tight = cset_btag_2018.at("deepCSV_incl")->evaluate({"central", "T", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_loose = cset_btag_2018.at("deepCSV_incl")->evaluate({"central", "L", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_tight_up = cset_btag_2018.at("deepCSV_incl")->evaluate({"up_uncorrelated", "T", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_loose_up = cset_btag_2018.at("deepCSV_incl")->evaluate({"up_uncorrelated", "L", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_tight_down = cset_btag_2018.at("deepCSV_incl")->evaluate({"down_uncorrelated", "T", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_sf_loose_down = cset_btag_2018.at("deepCSV_incl")->evaluate({"down_uncorrelated", "L", jetflavor[i], abs(eta[i]), pt[i]});
//                     btag_eff_tight = cset_btag_eff.at("btag_2018")->evaluate({"L", "T", pt[i], eta[i]});
//                     btag_eff_loose = cset_btag_eff.at("btag_2018")->evaluate({"L", "L", pt[i], eta[i]});
//                 }
//                 num *= (btag_sf_tight * btag_eff_tight) * (btag_sf_loose * btag_eff_loose - btag_sf_tight * btag_eff_tight) * (1. - btag_sf_loose * btag_eff_loose);
//                 num_up *= (btag_sf_tight_up * btag_eff_tight) * (btag_sf_loose_up * btag_eff_loose - btag_sf_tight_up * btag_eff_tight) * (1. - btag_sf_loose_up * btag_eff_loose);
//                 num_down *= (btag_sf_tight_down * btag_eff_tight) * (btag_sf_loose_down * btag_eff_loose - btag_sf_tight_down * btag_eff_tight) * (1. - btag_sf_loose_down * btag_eff_loose);
//                 den *= (btag_eff_tight) * (btag_eff_loose - btag_eff_tight) * (1. - btag_eff_loose);
//             }
//         }
//         if (den == 0.) {
//             btag_sf_weights[0] = 1.;
//             btag_sf_weights[1] = 1.;
//             btag_sf_weights[2] = 1.;
//         }
//         else {
//             btag_sf_weights[0] = num / den;
//             btag_sf_weights[1] = num_up / den;
//             btag_sf_weights[2] = num_down / den;
//         }
//         return btag_sf_weights;
//     };
//     return df.Define("btagging_scale_factors_LF", eval_correction, {"sample_year", "LnTBJet_eta", "LnTBJet_pt", "LnTBJet_hadronFlavour"});
// }

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