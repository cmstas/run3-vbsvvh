#include "selections_OneLep2FJ.h"

namespace OneLep2FJ {
    TMVA::Experimental::RBDT vbstagger("VBSTagger", "mva/VBSTagger/OneLep2FJ.root");
    SPANet::SPANetInference spanet_inference("/home/users/aaarora/phys/run3/SPANet/spanet.onnx");

    RNode TriggerSelections(RNode df_) {
        return df_.Define("_cut_trigger", "HLT_Ele30_WPTight_Gsf || HLT_IsoMu24");
    }

    RNode runPreselection(RNode df_) {
        // lepton
        auto df = EventFilters(df_);
        
        //cast GenPart_genPartIdxMother to int
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
        
        df = df.Define("_cut_hbb", "hbb_score > 0")
            .Define("_cut_wqq", "wqq_score > 0");

        df = df.Define("st", "lepton_pt + vbs1_pt + vbs2_pt + hbb_pt + wqq_pt")
            .Define("_cut_st", "st > 1000")
            .Define("_cut_hbb_score", "hbb_score > 0.2")
            .Define("_cut_wqq_score", "wqq_score > 0.2");

        df = df.Define("hbb_idx", "_max_hbb_idx")
            .Define("wqq_idx", "_max_wqq_idx");

        df = spanet_inference.RunSPANetInference(df);

        df = df.Define("spanet_vbs1_idx", "spanet_outputs[4][1]")
            .Define("spanet_vbs2_idx", "spanet_outputs[4][2]")
            .Define("spanet_h1_idx", "spanet_outputs[0][1]")
            .Define("spanet_h2_idx", "spanet_outputs[0][2]")
            .Define("spanet_v1_idx", "spanet_outputs[1][1]")
            .Define("spanet_v2_idx", "spanet_outputs[1][2]")
            .Define("spanet_bh_idx", "spanet_outputs[2][1]-10")
            .Define("spanet_bv_idx", "spanet_outputs[2][1]-10");

        // VBS jets matching to generator-level partons
        df = df.Define("truth_vbs1_eta", "GenPart_eta[5]")
            .Define("truth_vbs1_phi", "GenPart_phi[5]")
            .Define("truth_vbs2_eta", "GenPart_eta[6]")
            .Define("truth_vbs2_phi", "GenPart_phi[6]")
            .Define("vbs1_dR", VdR, {"ak4jet_eta", "ak4jet_phi", "truth_vbs1_eta", "truth_vbs1_phi"})
            .Define("vbs2_dR", VdR, {"ak4jet_eta", "ak4jet_phi", "truth_vbs2_eta", "truth_vbs2_phi"})
            .Define("vbs1_closest_jet", "ROOT::VecOps::Min(vbs1_dR)")
            .Define("vbs2_closest_jet", "ROOT::VecOps::Min(vbs2_dR)")
            .Define("truth_vbs1_idx", "(vbs1_closest_jet < 0.4) ? (int)ROOT::VecOps::ArgMin(vbs1_dR) : -1")
            .Define("truth_vbs2_idx", "(vbs2_closest_jet < 0.4) ? (int)ROOT::VecOps::ArgMin(vbs2_dR) : -1");

        // Boosted Higgs matching
        df = df.Define("truth_higgs_eta", "GenPart_eta[4]")
            .Define("truth_higgs_phi", "GenPart_phi[4]")
            .Define("higgs_dR", VdR, {"ak8jet_eta", "ak8jet_phi", "truth_higgs_eta", "truth_higgs_phi"})
            .Define("higgs_closest_fatjet", "ROOT::VecOps::Min(higgs_dR)")
            .Define("truth_bh_idx", "(higgs_closest_fatjet < 0.8) ? (int)ROOT::VecOps::ArgMin(higgs_dR) : -1");
            
        // Resolved Higgs matching (b-quarks from Higgs decay)
        df = df.Define("truth_b1_eta", "GenPart_eta[abs(GenPart_pdgId) == 5 && GenPart_genPartIdxMother == 4][0]")
            .Define("truth_b1_phi", "GenPart_phi[abs(GenPart_pdgId) == 5 && GenPart_genPartIdxMother == 4][0]")
            .Define("truth_b2_eta", "GenPart_eta[abs(GenPart_pdgId) == 5 && GenPart_genPartIdxMother == 4][1]")
            .Define("truth_b2_phi", "GenPart_phi[abs(GenPart_pdgId) == 5 && GenPart_genPartIdxMother == 4][1]")
            .Define("b1_dR", VdR, {"ak4jet_eta", "ak4jet_phi", "truth_b1_eta", "truth_b1_phi"})
            .Define("b2_dR", VdR, {"ak4jet_eta", "ak4jet_phi", "truth_b2_eta", "truth_b2_phi"})
            .Define("b1_closest_jet", "ROOT::VecOps::Min(b1_dR)")
            .Define("b2_closest_jet", "ROOT::VecOps::Min(b2_dR)")
            .Define("truth_h1_idx", "(b1_closest_jet < 0.4) ? (int)ROOT::VecOps::ArgMin(b1_dR) : -1")
            .Define("truth_h2_idx", "(b2_closest_jet < 0.4) ? (int)ROOT::VecOps::ArgMin(b2_dR) : -1");

        // Get the hadronic vector boson index (W or Z)
        df = df.Define("V_idx", get_hadronic_gauge_boson_idx, {"GenPart_pdgId", "GenPart_genPartIdxMother"})
            .Filter("V_idx == 2 || V_idx == 3");

        // Boosted vector boson matching
        df = df.Define("V_eta", "GenPart_eta[V_idx]")
            .Define("V_phi", "GenPart_phi[V_idx]")
            .Define("V_dR", VdR, {"ak8jet_eta", "ak8jet_phi", "V_eta", "V_phi"})
            .Define("V_closest_fatjet", "V_dR.size() > 0 ? ROOT::VecOps::Min(V_dR) : 999.0")
            .Define("truth_bv_idx", "(V_closest_fatjet < 0.8) ? (int)ROOT::VecOps::ArgMin(V_dR) : -1");
        
        // Resolved vector boson matching (quarks from W/Z decay)
        df = df.Define("truth_q1_eta", "GenPart_eta[abs(GenPart_pdgId) <= 5 && GenPart_genPartIdxMother == V_idx][0]")
            .Define("truth_q1_phi", "GenPart_phi[abs(GenPart_pdgId) <= 5 && GenPart_genPartIdxMother == V_idx][0]")
            .Define("truth_q2_eta", "GenPart_eta[abs(GenPart_pdgId) <= 5 && GenPart_genPartIdxMother == V_idx][1]")
            .Define("truth_q2_phi", "GenPart_phi[abs(GenPart_pdgId) <= 5 && GenPart_genPartIdxMother == V_idx][1]")
            .Define("q1_dR", VdR, {"ak4jet_eta", "ak4jet_phi", "truth_q1_eta", "truth_q1_phi"})
            .Define("q2_dR", VdR, {"ak4jet_eta", "ak4jet_phi", "truth_q2_eta", "truth_q2_phi"})
            .Define("q1_closest_jet", "ROOT::VecOps::Min(q1_dR)")
            .Define("q2_closest_jet", "ROOT::VecOps::Min(q2_dR)")
            .Define("truth_v1_idx", "(q1_closest_jet < 0.4) ? (int)ROOT::VecOps::ArgMin(q1_dR) : -1")
            .Define("truth_v2_idx", "(q2_closest_jet < 0.4) ? (int)ROOT::VecOps::ArgMin(q2_dR) : -1");

        return df;
    }

} // OneLep2FJ