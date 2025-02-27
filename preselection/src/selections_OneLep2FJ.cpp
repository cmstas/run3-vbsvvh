#include "commonSelections.h"
#include "selections_OneLep2FJ.h"

namespace OneLep2FJ {

    RNode TriggerSelections(RNode df_) {
        return df_.Define("passesTriggers", "HLT_Ele30_WPTight_Gsf || HLT_IsoMu24");
    }

    RNode ElectronSelections(RNode df_) {
        return df_.Define("Electron_SC_eta", "Electron_eta + Electron_deltaEtaSC")
                .Define("looseElectrons", 
                    "Electron_pt > 7 &&"
                    "abs(Electron_SC_eta) < 2.5 && "
                    "((abs(Electron_SC_eta) <= 1.479 && abs(Electron_dxy) <= 0.05 && abs(Electron_dz) < 0.1) || (abs(Electron_dxy) <= 0.1 && abs(Electron_dz) < 0.2)) && "
                    "abs(Electron_sip3d) < 8 && "
                    "Electron_cutBased >= 2 && "
                    "Electron_lostHits <= 1")
                .Define("nLooseElectrons", "nElectron == 0 ? 0 : Sum(looseElectrons)")
                .Define("tightElectrons", "looseElectrons &&" 
                    "Electron_cutBased >= 4 && "
                    "Electron_pt > 35 && "
                    "Electron_hoe < 0.1 && "
                    "Electron_eInvMinusPInv > -0.04 && "
                    "((abs(Electron_SC_eta) <= 1.479 && Electron_sieie < 0.011) || Electron_sieie <= 0.030) && "
                    "Electron_convVeto == true && "
                    "Electron_tightCharge == 2 && "
                    "Electron_lostHits == 0")
                .Define("nTightElectrons", "nElectron == 0 ? 0 : Sum(tightElectrons)")
                .Define("electron_pt", "Electron_pt[tightElectrons]")
                .Define("electron_eta", "Electron_eta[tightElectrons]")
                .Define("electron_sceta", "Electron_SC_eta[tightElectrons]")
                .Define("electron_phi", "Electron_phi[tightElectrons]")
                .Define("electron_mass", "Electron_mass[tightElectrons]")
                .Define("electron_charge", "Electron_charge[tightElectrons]")
                .Define("isElectron", "nLooseElectrons == 1 && nTightElectrons == 1");
    }

    RNode MuonSelections(RNode df_) {
        return df_.Define("looseMuons", 
                    "Muon_pt > 5 && "
                    "abs(Muon_eta) < 2.4 && "
                    "abs(Muon_dxy) < 0.2 && "
                    "abs(Muon_dz) < 0.5 && "
                    "abs(Muon_sip3d) < 8 && "
                    "Muon_looseId == 1")
                .Define("nLooseMuons", "nMuon == 0 ? 0 : Sum(looseMuons)")
                .Define("tightMuons", "looseMuons && "
                    "Muon_pt > 35 && "
                    "Muon_pfIsoId > 4 && "
                    "Muon_tightCharge == 2 && "
                    "Muon_highPurity && "
                    "Muon_tightId")
                .Define("nTightMuons", "nMuon == 0 ? 0 : Sum(tightMuons)")
                .Define("muon_pt", "Muon_pt[tightMuons]")
                .Define("muon_eta", "Muon_eta[tightMuons]")
                .Define("muon_phi", "Muon_phi[tightMuons]")
                .Define("muon_mass", "Muon_mass[tightMuons]")
                .Define("muon_charge", "Muon_charge[tightMuons]")
                .Define("isMuon", "nLooseMuons == 1 && nTightMuons == 1");
    }

    RNode LeptonSelections(RNode df_) {
        auto df = ElectronSelections(df_);
        df = MuonSelections(df);
        return df.Define("lepton_pt", "isElectron ? electron_pt[0] : muon_pt[0]")
                .Define("lepton_eta", "isElectron ? electron_eta[0] : muon_eta[0]")
                .Define("lepton_phi", "isElectron ? electron_phi[0] : muon_phi[0]")
                .Define("lepton_mass", "isElectron ? electron_mass[0] : muon_mass[0]");
    }

    RNode HbbSelections(RNode df_) {
        return df_.Define("XbbLepDeltaR", VfDeltaR, {"FatJet_eta", "FatJet_phi", "lepton_eta", "lepton_phi"})
                .Define("Xbbscore", "FatJet_particleNetMD_Xbb / (FatJet_particleNetMD_Xbb + FatJet_particleNetMD_QCD)")
                .Define("XbbJets", 
                    "CorrFatJet_pt > 250 && "
                    "XbbLepDeltaR >= 0.8 && "
                    "CorrFatJet_mass > 50 && "
                    "abs(FatJet_eta) <= 2.5 && "
                    "FatJet_msoftdrop > 40 && "
                    "FatJet_jetId > 0")
                .Define("xbbscore_idx", "Xbbscore.size() != 0 ? ArgMax(Xbbscore[XbbJets]) : 999.0")
                .Define("xbb_score", "xbbscore_idx != 999.0 ? Xbbscore[XbbJets][xbbscore_idx] : -999.0")
                .Define("xbb_pt", "CorrFatJet_pt[XbbJets][xbbscore_idx]")
                .Define("xbb_eta", "FatJet_eta[XbbJets][xbbscore_idx]")
                .Define("xbb_phi", "FatJet_phi[XbbJets][xbbscore_idx]")
                .Define("xbb_mass", "FatJet_particleNet_mass[XbbJets][xbbscore_idx]")
                .Define("xbb_msoftdrop", "FatJet_msoftdrop[XbbJets][xbbscore_idx]")
                .Define("xbb_tau2", "FatJet_tau2[XbbJets][xbbscore_idx]")
                .Define("xbb_tau1", "FatJet_tau1[XbbJets][xbbscore_idx]")
                .Define("xbb_tau21", "xbb_tau2 / xbb_tau1");
    }

    RNode VqqSelections(RNode df_) {
        return df_.Define("XWqqLepDeltaR", VfDeltaR, {"FatJet_eta", "FatJet_phi", "lepton_eta", "lepton_phi"})
                .Define("XWqqXbbDeltaR", VfDeltaR, {"FatJet_eta", "FatJet_phi", "xbb_eta", "xbb_phi"})
                .Define("XWqqJets", 
                    "CorrFatJet_pt > 250 && "
                    "CorrFatJet_mass > 50 && "
                    "XWqqLepDeltaR >= 0.8 && "
                    "XWqqXbbDeltaR >= 0.8 && "
                    "abs(FatJet_eta) <= 2.5 && "
                    "FatJet_msoftdrop > 40 && "
                    "FatJet_jetId > 0")
                .Define("XWqqscore", "(FatJet_particleNetMD_Xqq + FatJet_particleNetMD_Xcc) / (FatJet_particleNetMD_Xqq + FatJet_particleNetMD_Xcc + FatJet_particleNetMD_QCD)")
                .Define("xwqqscore_idx", "XWqqscore.size() != 0 ? ArgMax(XWqqscore[XWqqJets]) : -1")
                .Define("xwqq_score", "xwqqscore_idx != -1 ? XWqqscore[XWqqJets][xwqqscore_idx] : -1")
                .Define("xwqq_pt", "CorrFatJet_pt[XWqqJets][xwqqscore_idx]")
                .Define("xwqq_eta", "FatJet_eta[XWqqJets][xwqqscore_idx]")
                .Define("xwqq_phi", "FatJet_phi[XWqqJets][xwqqscore_idx]")
                .Define("xwqq_mass", "FatJet_particleNet_mass[XWqqJets][xwqqscore_idx]")
                .Define("xwqq_msoftdrop", "FatJet_msoftdrop[XWqqJets][xwqqscore_idx]")
                .Define("xwqq_tau2", "FatJet_tau2[XWqqJets][xwqqscore_idx]")
                .Define("xwqq_tau1", "FatJet_tau1[XWqqJets][xwqqscore_idx]")
                .Define("xwqq_tau21", "xwqq_tau2 / xwqq_tau1");
    }

    // RNode AK4Selections(RNode df_) {
    //     return df_.Define("AK4LepDeltaR", VfDeltaR, {"Jet_eta", "Jet_phi", "lepton_eta", "lepton_phi"})
    //             .Define("AK4XbbDeltaR", VfDeltaR, { "Jet_eta", "Jet_phi", "xbb_eta", "xbb_phi"})
    //             .Define("AK4XWqqDeltaR", VfDeltaR, {"Jet_eta", "Jet_phi", "xwqq_eta", "xwqq_phi"})
    //             .Define("ak4tightBjetScore", tightDFBtagWP, {"sample_year"})
    //             .Define("ak4looseBjetScore", looseDFBtagWP, {"sample_year"})
    //             .Define("goodJets", "CorrJet_pt >= 20 && "
    //                 "abs(Jet_eta) < 2.5 && "
    //                 "AK4LepDeltaR >= 0.4 && "
    //                 "AK4XbbDeltaR >= 0.8 && "
    //                 "AK4XWqqDeltaR >= 0.8 && "
    //                 "((is2016 && Jet_jetId >= 1) || (!is2016 && Jet_jetId >= 2)) && "
    //                 "(CorrJet_pt >= 50 || (CorrJet_pt < 50 && Jet_puId != 0))")
    //             .Define("ak4FromBJet", "goodJets && Jet_btagDeepFlavB > ak4tightBjetScore")
    //             .Define("goodLooseBJets", "goodJets && Jet_btagDeepFlavB > ak4looseBjetScore")
    //             .Define("sortedBJets", "Argsort(-Jet_btagDeepFlavB[goodLooseBJets])")
    //             .Define("GBJet_pt", "Take(CorrJet_pt[goodLooseBJets], sortedBJets)")
    //             .Define("GBJet_eta", "Take(Jet_eta[goodLooseBJets], sortedBJets)")
    //             .Define("GBJet_phi", "Take(Jet_phi[goodLooseBJets], sortedBJets)")
    //             .Define("GBJet_mass", "Take(CorrJet_mass[goodLooseBJets], sortedBJets)")
    //             .Define("GBJet_score", "Take(Jet_btagDeepFlavB[goodLooseBJets], sortedBJets)")
    //             .Define("bjet1pt", "GBJet_pt.size() > 0 ? GBJet_pt[0] : -999")
    //             .Define("bjet1eta", "GBJet_pt.size() > 0 ? GBJet_eta[0] : -999")
    //             .Define("bjet1phi", "GBJet_pt.size() > 0 ? GBJet_phi[0] : -999")
    //             .Define("bjet1score", "GBJet_pt.size() > 0 ? GBJet_score[0] : -999")
    //             .Define("bjet2pt", "GBJet_pt.size() > 1 ? GBJet_pt[1] : -999")
    //             .Define("bjet2eta", "GBJet_pt.size() > 1 ? GBJet_eta[1] : -999")
    //             .Define("bjet2phi", "GBJet_pt.size() > 1 ? GBJet_phi[1] : -999")
    //             .Define("bjet2score", "GBJet_pt.size() > 1 ? GBJet_score[1] : -999")
    //             .Define("Mlb", VfInvariantMass, {"GBJet_pt", "GBJet_eta", "GBJet_phi", "GBJet_mass", "lepton_pt", "lepton_eta", "lepton_phi", "lepton_mass"})
    //             .Define("MinMlbJetIdx", "Mlb.size() != 0 ? ArgMin(Mlb) : -1")
    //             .Define("Mlbminloose", "Mlb.size() != 0 ? Mlb[MinMlbJetIdx] : 1000");
    // }

    // RNode VBSJetsSelections(RNode df_) {
    //     return df_.Define("goodVBSJets", "CorrJet_pt >= 30 && "
    //                 "abs(Jet_eta) <= 4.7 && "
    //                 "AK4LepDeltaR >= 0.4 && "
    //                 "AK4XbbDeltaR >= 0.8 && "
    //                 "AK4XWqqDeltaR >= 0.8 && "
    //                 "((is2016 && Jet_jetId >= 1) || (!is2016 && Jet_jetId >= 2)) && "
    //                 "(CorrJet_pt >= 50 || (CorrJet_pt < 50 && Jet_puId != 0))")
    //             .Define("VBSJets_pt", "CorrJet_pt[goodVBSJets]")
    //             .Define("VBSJets_eta", "Jet_eta[goodVBSJets]")
    //             .Define("VBSJets_phi", "Jet_phi[goodVBSJets]")
    //             .Define("VBSJets_mass", "CorrJet_mass[goodVBSJets]")
    //             .Define("VBSjetidxs", VBS_MaxE, {"VBSJets_pt", "VBSJets_eta", "VBSJets_phi", "VBSJets_mass"})
    //             .Define("vbs1_pt", "VBSJets_pt[VBSjetidxs[0]]")
    //             .Define("vbs1_eta", "VBSJets_eta[VBSjetidxs[0]]")
    //             .Define("vbs1_phi", "VBSJets_phi[VBSjetidxs[0]]")
    //             .Define("vbs1_mass", "VBSJets_mass[VBSjetidxs[0]]")
    //             .Define("vbs2_pt", "VBSJets_pt[VBSjetidxs[1]]")
    //             .Define("vbs2_eta", "VBSJets_eta[VBSjetidxs[1]]")
    //             .Define("vbs2_phi", "VBSJets_phi[VBSjetidxs[1]]")
    //             .Define("vbs2_mass", "VBSJets_mass[VBSjetidxs[1]]")
    //             .Define("vbs_ptjj", "VBSjet1pt + VBSjet2pt")
    //             .Define("vbs_detajj", "abs(VBSjet1eta - VBSjet2eta)")
    //             .Define("vbs_mjj", fInvariantMass, {"VBSjet1pt", "VBSjet1eta", "VBSjet1phi", "VBSjet1mass", "VBSjet2pt", "VBSjet2eta", "VBSjet2phi", "VBSjet2mass"})
    //             .Define("met", "CorrMET_pt")
    //             .Define("ST", "lepton_pt + met + xbb_pt + xwqq_pt");
    // }

    RNode ObjectSelections(RNode df_) {
        auto df = TriggerSelections(df_);
        df = LeptonSelections(df);
        //df = HbbSelections(df);
        //df = VqqSelections(df);
        //df = AK4Selections(df);
        //df = VBSJetsSelections(df);
        return df; 
    }

    RNode EventSelections(RNode df_) {
        return df_.Define("passes_filters_triggers", "passesEventFilters && passesTriggers")
            .Define("passes_lepton_selection", "passes_filters_triggers && "
                "((nLooseMuons == 1 && nTightMuons == 1 && nLooseElectrons == 0 && nTightElectrons == 0) || (nLooseMuons == 0 && nTightMuons == 0 && nLooseElectrons == 1 && nTightElectrons == 1)) && "
                "(lepton_pt > 40)")
            //.Define("passes_hqq_selection", "passes_lepton_selection && HighestHScore > 0")
            //.Define("passes_xwqq_selection", "passes_hqq_selection && HighestWjetScore > 0")
            ;
    }

    RNode runAnalysis(RNode df_) {
        auto df = CommonSelections(df_);
        df = ObjectSelections(df);
        df = EventSelections(df);
        return df;
    }

} // OneLep2FJ