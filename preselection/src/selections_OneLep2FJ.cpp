#include "selections_OneLep2FJ.h"

namespace OneLep2FJ {
    TMVA::Experimental::RBDT vbstagger("VBSTagger", "mva/VBSTagger/OneLep2FJ.root");

    RNode TriggerSelections(RNode df_) {
        return df_.Define("_cut_trigger", "HLT_Ele30_WPTight_Gsf || HLT_IsoMu24");
    }

    RNode runPreselection(RNode df_) {
        // lepton
        auto df = EventFilters(df_);
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

        // df = AK4JetsSelection(df);
        // df = df.Define("_ak4_lep_dR", VdR, {"ak4jet_eta", "ak4jet_phi", "lepton_eta", "lepton_phi"})
        //     .Define("_better_ak4jets", "Sum(_ak4_lep_dR > 0.4)");

        // df = VBSJetSelections(df, vbstagger);
        // df = df.Define("_cut_vbs", "_better_ak4jets >= 2 && vbs_score > 0.4");
                
        df = AK8JetsSelection(df);
        df = df.Define("_ak8_lep_dR", VdR, {"ak8jet_eta", "ak8jet_phi", "lepton_eta", "lepton_phi"})
            // .Define("_ak8_vbs1_dR", VdR, {"ak8jet_eta", "ak8jet_phi", "vbs1_eta", "vbs1_phi"})
            // .Define("_ak8_vbs2_dR", VdR, {"ak8jet_eta", "ak8jet_phi", "vbs2_eta", "vbs2_phi"});
            // .Define("_hbb_candidates", "_ak8_lep_dR > 0.8 && _ak8_vbs1_dR > 0.8 && _ak8_vbs2_dR > 0.8")
            .Define("_hbb_candidates", "_ak8_lep_dR > 0.8")
            .Define("_hbbvsqcd_score", "ak8jet_hbbvsqcd[_hbb_candidates]")
            .Define("_max_hbb_idx", "_hbbvsqcd_score.size() != 0 ? ArgMax(_hbbvsqcd_score) : -999.0")
            .Define("hbb_score", "_max_hbb_idx != -999.0 ? _hbbvsqcd_score[_max_hbb_idx] : -999.0")
            .Define("hbb_pt", "ak8jet_pt[_hbb_candidates][_max_hbb_idx]")
            .Define("hbb_eta", "ak8jet_eta[_hbb_candidates][_max_hbb_idx]")
            .Define("hbb_phi", "ak8jet_phi[_hbb_candidates][_max_hbb_idx]")
            .Define("hbb_msoftdrop", "ak8jet_msoftdrop[_hbb_candidates][_max_hbb_idx]");

        df = df.Define("_ak8_hbb_dR", VdR, {"ak8jet_eta", "ak8jet_phi", "hbb_eta", "hbb_phi"})
            // .Define("_wqq_candidates", "_ak8_lep_dR > 0.8 && _ak8_vbs1_dR > 0.8 && _ak8_vbs2_dR > 0.8 && _ak8_hbb_dR > 0.8")
            .Define("_wqq_candidates", "_ak8_lep_dR > 0.8 && _ak8_hbb_dR > 0.8")
            .Define("_wqqvsqcd_score", "ak8jet_wvsqcd[_wqq_candidates]")
            .Define("_max_wqq_idx", "_wqqvsqcd_score.size() != 0 ? ArgMax(_wqqvsqcd_score) : -999.0")
            .Define("wqq_score", "_max_wqq_idx != -999.0 ? _wqqvsqcd_score[_max_wqq_idx] : -999.0")
            .Define("wqq_pt", "ak8jet_pt[_max_wqq_idx]")
            .Define("wqq_eta", "ak8jet_eta[_max_wqq_idx]")
            .Define("wqq_phi", "ak8jet_phi[_max_wqq_idx]")
            .Define("wqq_msoftdrop", "ak8jet_msoftdrop[_max_wqq_idx]");

        df = df.Define("_cut_hbb", "hbb_score > 0.5")
            .Define("_cut_wqq", "wqq_score > 0.5");
        
        df = AK4JetsSelection(df);
        df = df.Define("_ak4_lep_dR", VdR, {"ak4jet_eta", "ak4jet_phi", "lepton_eta", "lepton_phi"})
            .Define("ak4_hbb_dR", VdR, {"ak4jet_eta", "ak4jet_phi", "hbb_eta", "hbb_phi"})
            .Define("ak4_wqq_dR", VdR, {"ak4jet_eta", "ak4jet_phi", "wqq_eta", "wqq_phi"})
            // .Define("_better_ak4jets", "Sum(_ak4_lep_dR > 0.4)");
            .Define("_better_ak4jets", "Sum(_ak4_lep_dR > 0.4 && ak4_hbb_dR > 0.8 && ak4_wqq_dR > 0.8)");

        df = VBSJetSelections(df, vbstagger);
        df = df.Define("_cut_vbs", "_better_ak4jets >= 2 && vbs_score > 0.4");

        return df;
    }

} // OneLep2FJ