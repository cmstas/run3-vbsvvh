#include "selections_run2.h"

namespace Run2
{
    RNode runPreselection(RNode df_, std::string channel, bool noCut)
    {
        auto df = flagSelections(df_);
        df = defineCorrectedCols(df); // this is not doing anything right now, but kept it to avoid having to change branch names everywhere
        //df = redefineGenColumns(df);
        df = AK8JetsPreselection(df);
        df = ElectronPreselection(df);
        df = MuonPreselection(df);
        df = AK4JetsPreselection(df);
        df = VBSJetsPreselection(df);

        df = MuonSelections(df);
        df = ElectronSelections(df);

        if (noCut)
            return df; // for spanet training data

        // df = TriggerSelections(df, channel, TriggerMap);
        //  channel-specific selections
        //  if (channel == "1Lep2FJ") {
        //      df = df.Filter("((nMuon_Loose == 1 && nMuon_Tight == 1 && nElectron_Loose == 0 && nElectron_Tight == 0) || "
        //          "(nMuon_Loose == 0 && nMuon_Tight == 0 && nElectron_Loose == 1 && nElectron_Tight == 1)) && "
        //          "(Lepton_pt[0] > 40)", "C2: 1-lepton selection");
        //  }
        //  else if (channel == "0Lep3FJ") {
        //      df = df.Filter("nMuon_Loose == 0 && nElectron_Loose == 0", "C2: 0-lepton selection");
        //  }

        return df;
    }

    /*
     *  Hack
    */
    // RNode redefineGenColumns(RNode df)
    // {
    //     return df
    //         .Alias("gen_vbs1_idx", "truthVBSq1_idx")
    //         .Alias("gen_vbs2_idx", "truthVBSq2_idx")
    //         .Alias("gen_h_idx", "truthH_idx")
    //         .Alias("gen_b1_idx", "truthH_daught1_idx")
    //         .Alias("gen_b2_idx", "truthH_daught2_idx")
    //         .Alias("gen_v1_idx", "truthV1_idx")
    //         .Alias("gen_v2_idx", "truthV2_idx")
    //         .Alias("gen_v1q1_idx", "truthV1_daught1_idx")
    //         .Alias("gen_v1q2_idx", "truthV1_daught2_idx")
    //         .Alias("gen_v2q1_idx", "truthV2_daught1_idx")
    //         .Alias("gen_v2q2_idx", "truthV2_daught2_idx");
    // }
    /*
     * Corrections 
     */
    RNode defineCorrectedCols(RNode df)
    {
        return df.Define("CorrJet_pt", "Jet_pt")
            .Define("CorrJet_mass", "Jet_mass")
            .Define("CorrFatJet_pt", "FatJet_pt")
            .Define("CorrFatJet_mass", "FatJet_mass")
            .Define("CorrMET_pt", "MET_pt");
    }
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

    /*
     * Triggers
     */
    RNode TriggerSelections(RNode df_, std::string channel, const std::unordered_map<std::string, std::string> &trigger_map)
    {
        if (trigger_map.empty())
        {
            std::cerr << "Warning: No trigger map provided. Skipping trigger selection." << std::endl;
            return df_;
        }
        if (trigger_map.find(channel) == trigger_map.end())
        {
            std::cerr << "Warning: Channel '" << channel << "' not found in trigger map. Skipping trigger selection." << std::endl;
            return df_;
        }

        std::string trigger_condition = trigger_map.at(channel);
        return df_.Filter(trigger_condition, "C1: Trigger Selection");
    }

    /*
     *   Event filters
     */
    RNode flagSelections(RNode df_)
    {
        std::cout << " -> Run commonSelections::flagSelections()" << std::endl;
        auto df = df_.Define("Pass_EventFilters",
                                "Flag_goodVertices &&"
                                "Flag_HBHENoiseFilter &&"
                                "Flag_HBHENoiseIsoFilter &&"
                                "Flag_EcalDeadCellTriggerPrimitiveFilter &&"
                                "Flag_BadPFMuonFilter &&"
                                "Flag_BadPFMuonDzFilter &&"
                                "Flag_hfNoisyHitsFilter &&"
                                "Flag_eeBadScFilter &&"
                                "( (is2016) || Flag_ecalBadCalibFilter) &&" // apply only to 2017, 2018
                                //"( (!isData) || Flag_globalSuperTightHalo2016Filter)" // apply only to data
                                "Flag_globalSuperTightHalo2016Filter");
        return df;
    }

    /*
     *   AK8 Jets selection
     */
    RNode AK8JetsPreselection(RNode df_)
    {
        std::cout << " -> Run commonSelections::AK8JetsPreselection()" << std::endl;
        // Select good AK8 jets
        auto df = df_.Define("goodAK8Jets",
                             "CorrFatJet_pt > 300 && "
                             "abs(FatJet_eta) <= 2.5 && "
                             "FatJet_msoftdrop > 40 && "
                             "FatJet_jetId > 0")
                      .Define("FatJet_HbbScore", "FatJet_particleNetMD_Xbb / (FatJet_particleNetMD_Xbb + FatJet_particleNetMD_QCD)")
                      .Define("FatJet_WqqScore", "(FatJet_particleNetMD_Xcc + FatJet_particleNetMD_Xqq) / (FatJet_particleNetMD_Xcc + FatJet_particleNetMD_Xqq + FatJet_particleNetMD_QCD)")
                      .Define("goodAK8Jets_pt", "CorrFatJet_pt[goodAK8Jets]")
                      .Define("goodAK8Jets_eta", "FatJet_eta[goodAK8Jets]")
                      .Define("goodAK8Jets_phi", "FatJet_phi[goodAK8Jets]")
                      .Define("goodAK8Jets_mass", "CorrFatJet_mass[goodAK8Jets]")
                      .Define("goodAK8Jets_msoftdrop", "FatJet_msoftdrop[goodAK8Jets]")
                      .Define("goodAK8Jets_particleNet_mass", "FatJet_particleNet_mass[goodAK8Jets]")
                      .Define("goodAK8Jets_HbbScore", "FatJet_HbbScore[goodAK8Jets]")
                      .Define("goodAK8Jets_WqqScore", "FatJet_WqqScore[goodAK8Jets]")
                      .Define("goodAK8Jets_nConstituents", "FatJet_nConstituents[goodAK8Jets]")
                      .Define("ht_goodAK8Jets", "Sum(goodAK8Jets_pt)")
                      .Define("n_goodAK8Jets", "Sum(goodAK8Jets)")
                      .Define("ptSortedGoodAK8Jets", "Argsort(-goodAK8Jets_pt)");

        return df;
    }

    /*
     *   Lepton selection
     */
    RNode ElectronPreselection(RNode df_)
    {
        std::cout << " -> Run commonSelections::ElectronPreselection()" << std::endl;
        // Select good AK8 jets
        auto df = df_.Define("goodEle",
                                "Electron_pt > 10 && "
                                "abs(Electron_eta) <= 2.5")
                          .Define("goodElectron_pt", "Electron_pt[goodEle]")
                          .Define("goodElectron_eta", "Electron_eta[goodEle]");

        return df;
    }
    RNode MuonPreselection(RNode df_)
    {
        std::cout << " -> Run commonSelections::MuonPreselection()" << std::endl;
        // Select good AK8 jets
        auto df = df_.Define("goodMuo",
                                "Muon_pt > 10 && "
                                "abs(Muon_eta) <= 2.5")
                          .Define("goodMuon_pt", "Muon_pt[goodMuo]")
                          .Define("goodMuon_eta", "Muon_eta[goodMuo]");

        return df;
    }

    /*
     *   AK4 and VBS Jets selection
     */
    RNode AK4JetsPreselection(RNode df_)
    {
        std::cout << " -> Run commonSelections::AK4JetsPreselection()" << std::endl;
        auto df = df_.Define("ak4tightBjetScore", tightDFBtagWP, {"year"})
                          .Define("ak4mediumBjetScore", mediumDFBtagWP, {"year"})
                          .Define("ak4looseBjetScore", looseDFBtagWP, {"year"})
                          .Define("Jet_isTightBTag", "Jet_btagDeepFlavB > ak4tightBjetScore")
                          .Define("Jet_isMediumBTag", "Jet_btagDeepFlavB > ak4mediumBjetScore")
                          .Define("Jet_isLooseBTag", "Jet_btagDeepFlavB > ak4looseBjetScore")
                          .Define("Jet_minDrFromAnyGoodAK8Jet", dRfromClosestJet, {"Jet_eta", "Jet_phi", "goodAK8Jets_eta", "goodAK8Jets_phi"})
                          .Define("goodAK4Jets", "CorrJet_pt >= 20 && "
                                                 "((is2016 && Jet_jetId >= 1) || (!is2016 && Jet_jetId >= 2)) && "
                                                 "(CorrJet_pt >= 50 || (CorrJet_pt < 50 && Jet_puId != 0))")
                          .Define("goodAK4Jets_pt", "CorrJet_pt[goodAK4Jets]")
                          .Define("goodAK4Jets_eta", "Jet_eta[goodAK4Jets]")
                          .Define("goodAK4Jets_phi", "Jet_phi[goodAK4Jets]")
                          .Define("goodAK4Jets_mass", "CorrJet_mass[goodAK4Jets]")
                          .Define("goodAK4Jets_isTightBTag", "Jet_isTightBTag[goodAK4Jets]")
                          .Define("goodAK4Jets_isMediumBTag", "Jet_isMediumBTag[goodAK4Jets]")
                          .Define("goodAK4Jets_isLooseBTag", "Jet_isLooseBTag[goodAK4Jets]")
                          .Define("ht_goodAK4Jets", "Sum(CorrJet_pt[goodAK4Jets])")
                          .Define("n_goodAK4Jets", "Sum(goodAK4Jets)")
                          .Define("ptSortedGoodAK4Jets", "Argsort(-CorrJet_pt)") // checkme
                          .Define("goodAK4Jets_minDrFromAnyGoodAK8Jet", dRfromClosestJet, {"goodAK4Jets_eta", "goodAK4Jets_phi", "goodAK8Jets_eta", "goodAK8Jets_phi"})
                          .Define("goodAK4Jets_passAK8OverlapRemoval", "goodAK4Jets_minDrFromAnyGoodAK8Jet>0.8")
                          .Define("n_goodAK4JetsWithAK8OverlapRemoval", "Sum(goodAK4Jets_passAK8OverlapRemoval)");

        return df;
    }

    RNode VBSJetsPreselection(RNode df_)
    {
        std::cout << " -> Run commonSelections::VBSJetsPreselection()" << std::endl;
        auto df = df_.Define("goodVBSJets", "CorrJet_pt >= 20 && "
                                               "abs(Jet_eta) <= 4.7 && "
                                               "((is2016 && Jet_jetId >= 1) || (!is2016 && Jet_jetId >= 2)) && "
                                               "(CorrJet_pt >= 50 || (CorrJet_pt < 50 && Jet_puId != 0))")
                          .Define("goodVBSJets_pt", "CorrJet_pt[goodVBSJets]")
                          .Define("goodVBSJets_eta", "Jet_eta[goodVBSJets]")
                          .Define("goodVBSJets_phi", "Jet_phi[goodVBSJets]")
                          .Define("goodVBSJets_mass", "CorrJet_mass[goodVBSJets]");
        return df;
    }

    RNode ElectronSelections(RNode df_)
    {
        return df_.Define("_Electron_SC_eta", "Electron_eta + Electron_deltaEtaSC")
            .Define("_looseElectrons",
                    "Electron_pt > 7 &&"
                    "abs(_Electron_SC_eta) < 2.5 && "
                    "((abs(_Electron_SC_eta) <= 1.479 && abs(Electron_dxy) <= 0.05 && abs(Electron_dz) < 0.1) || (abs(Electron_dxy) <= 0.1 && abs(Electron_dz) < 0.2)) && "
                    "abs(Electron_sip3d) < 8 && "
                    "Electron_cutBased >= 2 && "
                    "Electron_pfRelIso03_all < 0.4 && "
                    "Electron_lostHits <= 1")
            .Define("electron_nloose", "nElectron == 0 ? 0 : Sum(_looseElectrons)")
            .Define("_tightElectrons", "_looseElectrons &&"
                                       "Electron_pt > 30 && "
                                       "Electron_cutBased >= 4 && "
                                       "Electron_pfRelIso03_all < 0.15 && "
                                       "Electron_hoe < 0.1 && "
                                       "Electron_eInvMinusPInv > -0.04 && "
                                       "((abs(_Electron_SC_eta) <= 1.479 && Electron_sieie < 0.011) || Electron_sieie <= 0.030) && "
                                       "Electron_convVeto == true && "
                                       "Electron_tightCharge == 2 && "
                                       "Electron_lostHits == 0")
            .Define("electron_ntight", "nElectron == 0 ? 0 : Sum(_tightElectrons)")
            .Define("vvhTightLepMaskElectron", "_tightElectrons")
            .Define("electron_pt", "Electron_pt[_tightElectrons]")
            .Define("electron_sceta", "_Electron_SC_eta[_tightElectrons]")
            .Define("electron_phi", "Electron_phi[_tightElectrons]")
            .Define("electron_mass", "Electron_mass[_tightElectrons]");
    }

    RNode MuonSelections(RNode df_)
    {
        return df_.Define("_looseMuons",
                          "Muon_pt > 5 && "
                          "Muon_pfIsoId >= 2 && "
                          "abs(Muon_eta) < 2.4 && "
                          "abs(Muon_dxy) < 0.2 && "
                          "abs(Muon_dz) < 0.5 && "
                          "abs(Muon_sip3d) < 8 && "
                          "Muon_looseId == 1")
            .Define("muon_nloose", "nMuon == 0 ? 0 : Sum(_looseMuons)")
            .Define("_tightMuons", "_looseMuons && "
                                   "Muon_pt > 30 && "
                                   "Muon_pfIsoId > 4 && "
                                   "Muon_tightCharge == 2 && "
                                   "Muon_highPurity && "
                                   "Muon_tightId")
            .Define("muon_ntight", "nMuon == 0 ? 0 : Sum(_tightMuons)")
            .Define("vvhTightLepMaskMuon", "_tightMuons")
            .Define("muon_pt", "Muon_pt[_tightMuons]")
            .Define("muon_eta", "Muon_eta[_tightMuons]")
            .Define("muon_phi", "Muon_phi[_tightMuons]")
            .Define("muon_mass", "Muon_mass[_tightMuons]");
    }

} // end namespace