#include "commonSelections.h"

RNode CommonSelections(RNode df_) 
{
    auto df = EventFilters(df_);
    df = AK8JetsSelection(df);
    df = AK4JetsSelection(df);
    df = VBSJetsSelection(df);

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


