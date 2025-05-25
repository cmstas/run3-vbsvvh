#include "commonSelections.h"
#include "onnx-inferences.hpp"

RNode CommonSelections(RNode df_, bool runSPANet) 
{
    auto df = EventFilters(df_);
    df = AK8JetsSelection(df);
    df = AK4JetsSelection(df);
    df = VBSJetsSelection(df);
    if (runSPANet) {
        df_out = RunSPANetInference(df_out);
    }

    return df;
}

/*
 *   Event filters
 */
RNode EventFilters(RNode df_) 
{
    return df_.Define("passesEventFilters", "Flag_goodVertices && "
            "Flag_globalSuperTightHalo2016Filter && "
            "Flag_EcalDeadCellTriggerPrimitiveFilter && "
            "Flag_BadPFMuonFilter && "
            "Flag_BadPFMuonDzFilter && "
            "Flag_hfNoisyHitsFilter &&"
            "Flag_eeBadScFilter && "
            "Flag_ecalBadCalibFilter");
}


/*
 *   AK8 Jets selection
 */
RNode AK8JetsSelection(RNode df_) 
{
    auto df = df_.Define("goodAK8Jets",
                          "CorrFatJet_pt > 300 && "
                          "abs(FatJet_eta) <= 2.5 && "
                          "FatJet_msoftdrop > 40 && "
                          "FatJet_jetId > 0")
                    //.Define("FatJet_HbbScore", "FatJet_particleNetMD_Xbb / (FatJet_particleNetMD_Xbb + FatJet_particleNetMD_QCD)")
                    //.Define("FatJet_WqqScore", "(FatJet_particleNetMD_Xcc + FatJet_particleNetMD_Xqq) / (FatJet_particleNetMD_Xcc + FatJet_particleNetMD_Xqq + FatJet_particleNetMD_QCD)")
                    .Define("goodAK8Jets_pt", "CorrFatJet_pt[goodAK8Jets]")
                    .Define("goodAK8Jets_eta", "FatJet_eta[goodAK8Jets]")
                    .Define("goodAK8Jets_phi", "FatJet_phi[goodAK8Jets]")
                    .Define("goodAK8Jets_mass", "CorrFatJet_mass[goodAK8Jets]")
                    .Define("goodAK8Jets_msoftdrop", "FatJet_msoftdrop[goodAK8Jets]")
                    .Define("goodAK8Jets_nConstituents", "FatJet_nConstituents[goodAK8Jets]")
                    .Define("ht_goodAK8Jets", "Sum(goodAK8Jets_pt)")
                    .Define("n_goodAK8Jets", "Sum(goodAK8Jets)")
                    .Define("ptSortedGoodAK8Jets", "Argsort(-goodAK8Jets_pt)"); 

    return df;
}

/*
 *   AK4 Jets selection
 */
RNode AK4JetsSelection(RNode df_)
{
    auto df = df_.Define("goodAK4Jets", "CorrJet_pt >= 20"
                                           //" && ((is2016 && Jet_jetId >= 1) || (!is2016 && Jet_jetId >= 2)) && "
                                           //"(CorrJet_pt >= 50 || (CorrJet_pt < 50 && Jet_puId != 0))"
                        )
                    //.Define("ak4tightBjetScore", tightDFBtagWP, {"sample_year"})
                    //.Define("ak4mediumBjetScore", mediumDFBtagWP, {"sample_year"})
                    //.Define("ak4looseBjetScore", looseDFBtagWP, {"sample_year"})
                    //.Define("Jet_isTightBTag", "Jet_btagDeepFlavB > ak4tightBjetScore")
                    //.Define("Jet_isMediumBTag", "Jet_btagDeepFlavB > ak4mediumBjetScore")
                    //.Define("Jet_isLooseBTag", "Jet_btagDeepFlavB > ak4looseBjetScore")
                      .Define("goodAK4Jets_pt", "CorrJet_pt[goodAK4Jets]")
                      .Define("goodAK4Jets_eta", "Jet_eta[goodAK4Jets]")
                      .Define("goodAK4Jets_phi", "Jet_phi[goodAK4Jets]")
                      .Define("goodAK4Jets_mass", "CorrJet_mass[goodAK4Jets]")
                      .Define("ht_goodAK4Jets", "Sum(CorrJet_pt[goodAK4Jets])")
                      .Define("n_goodAK4Jets", "Sum(goodAK4Jets)")
                      .Define("ptSortedGoodAK4Jets", "Argsort(-CorrJet_pt)") 
                      .Define("goodAK4Jets_minDrFromAnyGoodAK8Jet", VfdRfromClosestJet, {"goodAK4Jets_eta", "goodAK4Jets_phi", "goodAK8Jets_eta", "goodAK8Jets_phi"})
                      .Define("goodAK4Jets_passAK8OR", "goodAK4Jets_minDrFromAnyGoodAK8Jet>0.8")
                      .Define("n_goodAK4JetsWithAK8OR", "Sum(goodAK4Jets_passAK8OR)"); 

    return df;
}

/*
 *   VBS Jets selection
 */
RNode VBSJetsSelection(RNode df_)
{
    auto df = df_.Define("goodVBSJets", "CorrJet_pt >= 20 && "
                                        "abs(Jet_eta) <= 4.7"
                                           //" && ((is2016 && Jet_jetId >= 1) || (!is2016 && Jet_jetId >= 2)) && "
                                           //"(CorrJet_pt >= 50 || (CorrJet_pt < 50 && Jet_puId != 0))"
                        )
                      .Define("goodVBSJets_pt", "CorrJet_pt[goodVBSJets]")
                      .Define("goodVBSJets_eta", "Jet_eta[goodVBSJets]")
                      .Define("goodVBSJets_phi", "Jet_phi[goodVBSJets]")
                      .Define("goodVBSJets_mass", "CorrJet_mass[goodVBSJets]");
    return df;
}



RNode RunSPANetInference(RNode df)
{
    const char *model_path = "/home/users/mmazza/public/forAashay/spanet_v31_model/spanet.onnx";
    Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "ONNXInference");

    Ort::SessionOptions session_options;
    session_options.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL); // Max optimization
    session_options.DisableProfiling(); // Explicitly disable profiling, which adds computational overhead
    static Ort::Session session(env, model_path, session_options);

    auto runInference = [&session](const RVecF &ak4_pt, const RVecF &ak4_eta, const RVecF &ak4_phi, const RVecF &ak4_mass,
                                   const RVecI &ak4_isTightBTag, const RVecI &ak4_isMediumBTag, const RVecI &ak4_isLooseBTag,
                                   const RVecF &ak8_pt, const RVecF &ak8_eta, const RVecF &ak8_phi, const RVecF &ak8_msoftdrop,
                                   const RVecF &ak8_HbbScore, const RVecF &ak8_WqqScore, const RVecUCar &ak8_nConstituents,
                                   float met_pt, float ht_ak4, float ht_ak8)
    {
        // Prepare AK4 jet inputs (up to 10 jets, 8 features each)
        std::vector<std::vector<float>> ak4_jets(10, std::vector<float>(8, 0.0f));
        std::vector<bool> ak4_mask(10, false);
        for (size_t i = 0; i < std::min<size_t>(10, ak4_pt.size()); ++i)
        {
            ak4_mask[i] = true;
            ak4_jets[i] = {
                std::log(ak4_mass[i] + 1.0f),            // Feature 5: log(mass+1)
                std::log(ak4_pt[i] + 1.0f),              // Feature 1: log(pT+1)
                ak4_eta[i],                              // Feature 2: eta
                std::sin(ak4_phi[i]),                    // Feature 4: sin(phi)
                std::cos(ak4_phi[i]),                    // Feature 3: cos(phi)
                static_cast<float>(ak4_isTightBTag[i]),  // Feature 6: tight b-tag
                static_cast<float>(ak4_isMediumBTag[i]), // Feature 7: medium b-tag
                static_cast<float>(ak4_isLooseBTag[i])   // Feature 8: loose b-tag
            };
        }

        // Prepare AK8 jet inputs (exactly 3 jets, 8 features each)
        std::vector<std::vector<float>> ak8_jets(3, std::vector<float>(8, 0.0f));
        std::vector<bool> ak8_mask(3, false);
        for (size_t i = 0; i < std::min<size_t>(3, ak8_pt.size()); ++i)
        {
            ak8_mask[i] = true;
            ak8_jets[i] = {
                std::log(ak8_msoftdrop[i] +  1.0f),       // Feature 5: log(msoftdrop+1)
                std::log(ak8_pt[i] +  1.0f),              // Feature 1: log(pT+1)
                ak8_eta[i],                               // Feature 2: eta
                std::sin(ak8_phi[i]),                     // Feature 4: sin(phi)
                std::cos(ak8_phi[i]),                     // Feature 3: cos(phi)
                static_cast<float>(ak8_nConstituents[i]), // Feature 8: nConstituents
                ak8_HbbScore[i],                          // Feature 6: Hbb score
                ak8_WqqScore[i]                           // Feature 7: Wqq score
            };
        }

        // Event-level inputs (adjust as needed)
        float event_input_1 = std::log(met_pt + 1.0f); // log(MET)
        float event_input_2 = std::log(ht_ak4 + 1.0f); // HT from AK4 jets
        float event_input_3 = std::log(ht_ak8 + 1.0f); // HT from AK8 jets
        bool event_mask = true;                        // Event mask is always True

        // Run SPANet inference
        auto outputs = runSPANetInference( session, 
            ak4_jets[0], ak4_jets[1], ak4_jets[2], ak4_jets[3], ak4_jets[4],
            ak4_jets[5], ak4_jets[6], ak4_jets[7], ak4_jets[8], ak4_jets[9],
            ak4_mask,
            ak8_jets[0], ak8_jets[1], ak8_jets[2], ak8_mask,
            event_input_1, event_input_2, event_input_3, event_mask);


        // Flatten all outputs into a single RVec<float>
        RVecF flat_output;
        for (const auto &output : outputs) {
            flat_output.insert(flat_output.end(), output.begin(), output.end());
        }

        
        return flat_output;
    };


 
    const std::vector<std::string> input_columns = {
        "goodAK4Jets_pt", "goodAK4Jets_eta", "goodAK4Jets_phi", "goodAK4Jets_mass",
        "goodAK4Jets_isTightBTag", "goodAK4Jets_isMediumBTag", "goodAK4Jets_isLooseBTag",
        "goodAK8Jets_pt", "goodAK8Jets_eta", "goodAK8Jets_phi", "goodAK8Jets_msoftdrop",
        "goodAK8Jets_HbbScore", "goodAK8Jets_WqqScore", "goodAK8Jets_nConstituents",
        "MET_pt", "ht_goodAK4Jets", "ht_goodAK8Jets"};

    // FIXME: infer from NMAX_AK8 and NMAX_AK4
    const std::vector<std::pair<size_t, size_t>> offsets = {
        {0, 169},   // h_assignment_probability
        {169, 338}, // v1_assignment_probability
        {338, 507}, // v2_assignment_probability
        {507, 520}, // bh_assignment_probability
        {520, 533}, // bv1_assignment_probability
        {533, 546}, // bv2_assignment_probability
        {546, 715}, // vbs_assignment_probability
        {715, 716}, // h_detection_probability
        {716, 717}, // v1_detection_probability
        {717, 718}, // v2_detection_probability
        {718, 719}, // bh_detection_probability
        {719, 720}, // bv1_detection_probability
        {720, 721}, // bv2_detection_probability
        {721, 722}, // vbs_detection_probability
        {722, 724}  // EVENT/isSignal
    };

    const std::vector<std::string> output_names = {
        "h_assignment_probability",
        "v1_assignment_probability",
        "v2_assignment_probability",
        "bh_assignment_probability",
        "bv1_assignment_probability",
        "bv2_assignment_probability",
        "vbs_assignment_probability",
        "h_detection_probability",
        "v1_detection_probability",
        "v2_detection_probability",
        "bh_detection_probability",
        "bv1_detection_probability",
        "bv2_detection_probability",
        "vbs_detection_probability",
        "event_isSignal"
    };

    // Define temporary column for flattened outputs
    auto df_out = df.Define("spanet_outputs", runInference, input_columns);

    // Define final columns by slicing the flattened output
    for (size_t i = 0; i < output_names.size(); ++i)
    {

        df_out = df_out.Define(output_names[i],
                               [start = offsets[i].first, end = offsets[i].second, name = output_names[i], i = i](const RVecF &flat_output)
                               {            
                                   return RVecF(flat_output.begin() + start, flat_output.begin() + end);
                               },
                               {"spanet_outputs"});
    }

    return df_out;
}
