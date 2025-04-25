#include "selections_OneLep2FJ.h"

namespace OneLep2FJ {
    TMVA::Experimental::RBDT vbstagger("VBSTagger", "mva/VBSTagger/OneLep2FJ.root");

    RNode TriggerSelections(RNode df_) {
        return df_.Define("_cut_trigger", "HLT_Ele30_WPTight_Gsf || HLT_IsoMu24");
    }

    RNode runPreselection(RNode df_) {
        // lepton
        auto df = EventFilters(df_);

        // df = df.Define("GenPart_motherPdgId", "Take(GenPart_pdgId, GenPart_genPartIdxMother)")
        //         .Define("gauge_bosons", "Sum((abs(GenPart_pdgId) == 24 || abs(GenPart_pdgId) == 23) && abs(GenPart_motherPdgId) < 6)")
        //         .Filter("gauge_bosons == 2");

        df = TriggerSelections(df);

        df = ElectronSelections(df);
        df = MuonSelections(df);
            
        df = df.Define("_isElectron", "electron_nloose == 1 && electron_ntight == 1")
            .Define("lepton_pt", "_isElectron ? electron_pt[0] : muon_pt[0]")
            .Define("lepton_eta", "_isElectron ? electron_sceta[0] : muon_eta[0]")
            .Define("lepton_phi", "_isElectron ? electron_phi[0] : muon_phi[0]")
            .Define("lepton_mass", "_isElectron ? electron_mass[0] : muon_mass[0]")
            .Define("_cut_lepton", "((muon_nloose == 1 && muon_ntight == 1 && electron_nloose == 0 && electron_ntight == 0) || "
                "(muon_nloose == 0 && muon_ntight == 0 && electron_nloose == 1 && electron_ntight == 1)) && "
                "(lepton_pt > 40)");

        df = AK4JetsSelection(df);
        df = df.Define("_ak4_lep_dR", VdR, {"ak4jet_eta", "ak4jet_phi", "lepton_eta", "lepton_phi"})
            .Define("_better_ak4jets", "Sum(_ak4_lep_dR > 0.4)");

        df = VBSJetSelections(df, vbstagger);
        df = df.Define("_cut_vbs", "_better_ak4jets >= 2 && vbs_score > 0.4");
                
        df = AK8JetsSelection(df);
        df = df.Define("_ak8_lep_dR", VdR, {"ak8jet_eta", "ak8jet_phi", "lepton_eta", "lepton_phi"})
            .Define("_ak8_vbs1_dR", VdR, {"ak8jet_eta", "ak8jet_phi", "vbs1_eta", "vbs1_phi"})
            .Define("_ak8_vbs2_dR", VdR, {"ak8jet_eta", "ak8jet_phi", "vbs2_eta", "vbs2_phi"})
            .Define("_xbb_candidates", "_ak8_lep_dR > 0.8 && _ak8_vbs1_dR > 0.8 && _ak8_vbs2_dR > 0.8")
            .Define("_xbbvsqcd_score", "ak8jet_xbbvsqcd[_xbb_candidates]")
            .Define("_max_hbb_idx", "_xbbvsqcd_score.size() != 0 ? ArgMax(_xbbvsqcd_score) : -999")
            .Define("hbb_score", "_max_hbb_idx != -999 ? _xbbvsqcd_score[_max_hbb_idx] : -999.0f")
            .Define("hbb_pt", "_max_hbb_idx != -999 ? ak8jet_pt[_xbb_candidates][_max_hbb_idx] : -999.0f")
            .Define("hbb_eta", "_max_hbb_idx != -999 ? ak8jet_eta[_xbb_candidates][_max_hbb_idx] : -999.0f")
            .Define("hbb_phi", "_max_hbb_idx != -999 ? ak8jet_phi[_xbb_candidates][_max_hbb_idx] : -999.0f")
            .Define("hbb_msoftdrop", "_max_hbb_idx != -999 ? ak8jet_msoftdrop[_xbb_candidates][_max_hbb_idx] : -999.0f");

        df = df.Define("_ak8_hbb_dR", VdR, {"ak8jet_eta", "ak8jet_phi", "hbb_eta", "hbb_phi"})
            .Define("_xqq_candidates", "_ak8_lep_dR > 0.8 && _ak8_vbs1_dR > 0.8 && _ak8_vbs2_dR > 0.8 && _ak8_hbb_dR > 0.8")
            .Define("_xqqvsqcd_score", "ak8jet_xqqvsqcd[_xqq_candidates]")
            .Define("_max_wqq_idx", "_xqqvsqcd_score.size() != 0 ? ArgMax(_xqqvsqcd_score) : -999")
            .Define("wqq_score", "_max_wqq_idx != -999 ? _xqqvsqcd_score[_max_wqq_idx] : -999.0f")
            .Define("wqq_pt", "_max_wqq_idx != -999 ? ak8jet_pt[_xqq_candidates][_max_wqq_idx] : -999.0f")
            .Define("wqq_eta", "_max_wqq_idx != -999 ? ak8jet_eta[_xqq_candidates][_max_wqq_idx] : -999.0f")
            .Define("wqq_phi", "_max_wqq_idx != -999 ? ak8jet_phi[_xqq_candidates][_max_wqq_idx] : -999.0f")
            .Define("wqq_msoftdrop", "_max_wqq_idx != -999 ? ak8jet_msoftdrop[_xqq_candidates][_max_wqq_idx] : -999.0f");

        df = df.Define("_cut_hbb", "hbb_score > 0.2")
            .Define("_cut_wqq", "wqq_score > 0.2");

        return df;
    }

} // OneLep2FJ