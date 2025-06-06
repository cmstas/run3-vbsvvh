#include "selections_OneLep2FJ.h"

namespace OneLep2FJ {
    TMVA::Experimental::RBDT vbstagger("VBSTagger", "mva/VBSTagger/OneLep2FJ.root");
    SPANet::SPANetInference spanet_inference("/home/users/aaarora/phys/run3/SPANet/spanet.onnx", 8);

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
        // df = AK8JetsSelection(df);

        // df = spanet_inference.RunSPANetInference(df);
        // // spanet output order: h, v, bh, bv, vbs

        // df = df.Define("h_score", "spanet_h_detection")
        //     .Define("h1_idx", "spanet_h_assignment[1]")
        //     .Define("h2_idx", "spanet_h_assignment[2]")
        //     .Define("v_score", "spanet_v_detection")
        //     .Define("v1_idx", "spanet_v_assignment[1]")
        //     .Define("v2_idx", "spanet_v_assignment[2]")
        //     .Define("bh_score", "spanet_bh_detection")
        //     .Define("bh_idx", "spanet_bh_assignment[1]")
        //     .Define("bv_score", "spanet_bv_detection")
        //     .Define("bv_idx", "spanet_bv_assignment[1]")
        //     .Define("vbs_score", "spanet_vbs_detection")
        //     .Define("vbs1_idx", "spanet_vbs_assignment[1]")
        //     .Define("vbs2_idx", "spanet_vbs_assignment[2]");

        // // define vbs variables
        // df = df.Define("vbs1_pt", "ak4jet_pt[vbs1_idx]")
        //     .Define("vbs1_eta", "ak4jet_eta[vbs1_idx]")
        //     .Define("vbs1_phi", "ak4jet_phi[vbs1_idx]")
        //     .Define("vbs1_mass", "ak4jet_mass[vbs1_idx]")
        //     .Define("vbs2_pt", "ak4jet_pt[vbs2_idx]")
        //     .Define("vbs2_eta", "ak4jet_eta[vbs2_idx]")
        //     .Define("vbs2_phi", "ak4jet_phi[vbs2_idx]")
        //     .Define("vbs2_mass", "ak4jet_mass[vbs2_idx]")
        //     .Define("vbs_detajj", "abs(vbs1_eta - vbs2_eta)")
        //     .Define("vbs_mjj", "(ROOT::Math::PtEtaPhiMVector(vbs1_pt, vbs1_eta, vbs1_phi, vbs1_mass) + "
        //         "ROOT::Math::PtEtaPhiMVector(vbs2_pt, vbs2_eta, vbs2_phi, vbs2_mass)).M()");

        // df = df.Define("hbb_pt", "ak8jet_pt[bh_idx]")
        //     .Define("hbb_eta", "ak8jet_eta[bh_idx]")
        //     .Define("hbb_phi", "ak8jet_phi[bh_idx]")
        //     .Define("hbb_mass", "ak8jet_mass[bh_idx]")
        //     .Define("hbb_xbbscore", "ak8jet_xbbvsqcd[bh_idx]");

        // df = df.Define("wqq_pt", "ak8jet_pt[bv_idx]")
        //     .Define("wqq_eta", "ak8jet_eta[bv_idx]")
        //     .Define("wqq_phi", "ak8jet_phi[bv_idx]")
        //     .Define("wqq_mass", "ak8jet_mass[bv_idx]")
        //     .Define("wqq_xwqqscore", "ak8jet_xqqvsqcd[bv_idx]");

        // df = df.Define("st", "lepton_pt + vbs1_pt + vbs2_pt + hbb_pt + wqq_pt")
        //     .Define("_cut_st", "st > 1000");

        df = AK4JetsSelection(df);
        df = df.Define("_ak4_lep_dR", VdR, {"ak4jet_eta", "ak4jet_phi", "lepton_eta", "lepton_phi"})
            .Define("_better_ak4jets", "Sum(_ak4_lep_dR > 0.4)");
        
        df = VBSJetSelections(df, vbstagger);
                
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
        
        df = df.Define("hbb_idx", "_max_hbb_idx")
            .Define("wqq_idx", "_max_wqq_idx");

        df = spanet_inference.RunSPANetInference(df);

        df = df.Define("spanet_vbs1_idx", "spanet_vbs_assignment.size() > 1 ? static_cast<int>(spanet_vbs_assignment[1]) : -1")
            .Define("spanet_vbs2_idx", "spanet_vbs_assignment.size() > 2 ? static_cast<int>(spanet_vbs_assignment[2]) : -1")
            .Define("spanet_h1_idx", "spanet_h_assignment.size() > 1 ? static_cast<int>(spanet_h_assignment[1]) : -1")
            .Define("spanet_h2_idx", "spanet_h_assignment.size() > 2 ? static_cast<int>(spanet_h_assignment[2]) : -1")
            .Define("spanet_v1_idx", "spanet_v_assignment.size() > 1 ? static_cast<int>(spanet_v_assignment[1]) : -1")
            .Define("spanet_v2_idx", "spanet_v_assignment.size() > 2 ? static_cast<int>(spanet_v_assignment[2]) : -1")
            .Define("spanet_bh_idx", "spanet_bh_assignment.size() > 1 ? static_cast<int>(spanet_bh_assignment[1])-10 : -1")
            .Define("spanet_bv_idx", "spanet_bv_assignment.size() > 1 ? static_cast<int>(spanet_bv_assignment[1])-10 : -1");

        df = df.Define("spanet_vbs1_pt", "spanet_vbs1_idx >= 0 && spanet_vbs1_idx < ak4jet_pt.size() ? ak4jet_pt[spanet_vbs1_idx] : -999.0f")
            .Define("spanet_vbs1_eta", "spanet_vbs1_idx >= 0 && spanet_vbs1_idx < ak4jet_eta.size() ? ak4jet_eta[spanet_vbs1_idx] : -999.0f")
            .Define("spanet_vbs1_phi", "spanet_vbs1_idx >= 0 && spanet_vbs1_idx < ak4jet_phi.size() ? ak4jet_phi[spanet_vbs1_idx] : -999.0f")
            .Define("spanet_vbs1_mass", "spanet_vbs1_idx >= 0 && spanet_vbs1_idx < ak4jet_mass.size() ? ak4jet_mass[spanet_vbs1_idx] : -999.0f")
            .Define("spanet_vbs2_pt", "spanet_vbs2_idx >= 0 && spanet_vbs2_idx < ak4jet_pt.size() ? ak4jet_pt[spanet_vbs2_idx] : -999.0f")
            .Define("spanet_vbs2_eta", "spanet_vbs2_idx >= 0 && spanet_vbs2_idx < ak4jet_eta.size() ? ak4jet_eta[spanet_vbs2_idx] : -999.0f")
            .Define("spanet_vbs2_phi", "spanet_vbs2_idx >= 0 && spanet_vbs2_idx < ak4jet_phi.size() ? ak4jet_phi[spanet_vbs2_idx] : -999.0f")
            .Define("spanet_vbs2_mass", "spanet_vbs2_idx >= 0 && spanet_vbs2_idx < ak4jet_mass.size() ? ak4jet_mass[spanet_vbs2_idx] : -999.0f")
            .Define("spanet_vbs_detajj", "spanet_vbs1_idx >= 0 && spanet_vbs2_idx >= 0 ? abs(spanet_vbs1_eta - spanet_vbs2_eta) : -999.0f")
            .Define("spanet_vbs_mjj", "spanet_vbs1_idx >= 0 && spanet_vbs2_idx >= 0 ? (ROOT::Math::PtEtaPhiMVector(spanet_vbs1_pt, spanet_vbs1_eta, spanet_vbs1_phi, spanet_vbs1_mass) + "
                "ROOT::Math::PtEtaPhiMVector(spanet_vbs2_pt, spanet_vbs2_eta, spanet_vbs2_phi, spanet_vbs2_mass)).M() : -999.0f");

        df = df.Define("spanet_hbb_pt", "spanet_bh_idx >= 0 && spanet_bh_idx < ak8jet_pt.size() ? ak8jet_pt[spanet_bh_idx] : -999.0f")
            .Define("spanet_hbb_eta", "spanet_bh_idx >= 0 && spanet_bh_idx < ak8jet_eta.size() ? ak8jet_eta[spanet_bh_idx] : -999.0f")
            .Define("spanet_hbb_phi", "spanet_bh_idx >= 0 && spanet_bh_idx < ak8jet_phi.size() ? ak8jet_phi[spanet_bh_idx] : -999.0f")
            .Define("spanet_hbb_mass", "spanet_bh_idx >= 0 && spanet_bh_idx < ak8jet_mass.size() ? ak8jet_mass[spanet_bh_idx] : -999.0f")
            .Define("spanet_hbb_score", "spanet_bh_idx >= 0 && spanet_bh_idx < ak8jet_xbbvsqcd.size() ? ak8jet_xbbvsqcd[spanet_bh_idx] : -999.0f");

        df = df.Define("spanet_wqq_pt", "spanet_bv_idx >= 0 && spanet_bv_idx < ak8jet_pt.size() ? ak8jet_pt[spanet_bv_idx] : -999.0f")
            .Define("spanet_wqq_eta", "spanet_bv_idx >= 0 && spanet_bv_idx < ak8jet_eta.size() ? ak8jet_eta[spanet_bv_idx] : -999.0f")
            .Define("spanet_wqq_phi", "spanet_bv_idx >= 0 && spanet_bv_idx < ak8jet_phi.size() ? ak8jet_phi[spanet_bv_idx] : -999.0f")
            .Define("spanet_wqq_mass", "spanet_bv_idx >= 0 && spanet_bv_idx < ak8jet_mass.size() ? ak8jet_mass[spanet_bv_idx] : -999.0f")
            .Define("spanet_wqq_score", "spanet_bv_idx >= 0 && spanet_bv_idx < ak8jet_xqqvsqcd.size() ? ak8jet_xqqvsqcd[spanet_bv_idx] : -999.0f");


        // // GEN LEVEL MATCHING FOR VALIDATION ONLY
        // // VBS jets matching to generator-level partons
        // df = df.Define("truth_vbs1_eta", "GenPart_eta[5]")
        //     .Define("truth_vbs1_phi", "GenPart_phi[5]")
        //     .Define("truth_vbs2_eta", "GenPart_eta[6]")
        //     .Define("truth_vbs2_phi", "GenPart_phi[6]")
        //     .Define("vbs1_dR", VdR, {"ak4jet_eta", "ak4jet_phi", "truth_vbs1_eta", "truth_vbs1_phi"})
        //     .Define("vbs2_dR", VdR, {"ak4jet_eta", "ak4jet_phi", "truth_vbs2_eta", "truth_vbs2_phi"})
        //     .Define("vbs1_closest_jet", "ROOT::VecOps::Min(vbs1_dR)")
        //     .Define("vbs2_closest_jet", "ROOT::VecOps::Min(vbs2_dR)")
        //     .Define("truth_vbs1_idx", "(vbs1_closest_jet < 0.4) ? (int)ROOT::VecOps::ArgMin(vbs1_dR) : -1")
        //     .Define("truth_vbs2_idx", "(vbs2_closest_jet < 0.4) ? (int)ROOT::VecOps::ArgMin(vbs2_dR) : -1");

        // // Boosted Higgs matching
        // df = df.Define("truth_higgs_eta", "GenPart_eta[4]")
        //     .Define("truth_higgs_phi", "GenPart_phi[4]")
        //     .Define("higgs_dR", VdR, {"ak8jet_eta", "ak8jet_phi", "truth_higgs_eta", "truth_higgs_phi"})
        //     .Define("higgs_closest_fatjet", "ROOT::VecOps::Min(higgs_dR)")
        //     .Define("truth_bh_idx", "(higgs_closest_fatjet < 0.8) ? (int)ROOT::VecOps::ArgMin(higgs_dR) : -1");
            
        // // Resolved Higgs matching (b-quarks from Higgs decay)
        // df = df.Define("truth_b1_eta", "GenPart_eta[abs(GenPart_pdgId) == 5 && GenPart_genPartIdxMother == 4][0]")
        //     .Define("truth_b1_phi", "GenPart_phi[abs(GenPart_pdgId) == 5 && GenPart_genPartIdxMother == 4][0]")
        //     .Define("truth_b2_eta", "GenPart_eta[abs(GenPart_pdgId) == 5 && GenPart_genPartIdxMother == 4][1]")
        //     .Define("truth_b2_phi", "GenPart_phi[abs(GenPart_pdgId) == 5 && GenPart_genPartIdxMother == 4][1]")
        //     .Define("b1_dR", VdR, {"ak4jet_eta", "ak4jet_phi", "truth_b1_eta", "truth_b1_phi"})
        //     .Define("b2_dR", VdR, {"ak4jet_eta", "ak4jet_phi", "truth_b2_eta", "truth_b2_phi"})
        //     .Define("b1_closest_jet", "ROOT::VecOps::Min(b1_dR)")
        //     .Define("b2_closest_jet", "ROOT::VecOps::Min(b2_dR)")
        //     .Define("truth_h1_idx", "(b1_closest_jet < 0.4) ? (int)ROOT::VecOps::ArgMin(b1_dR) : -1")
        //     .Define("truth_h2_idx", "(b2_closest_jet < 0.4) ? (int)ROOT::VecOps::ArgMin(b2_dR) : -1");

        // // Get the hadronic vector boson index (W or Z)
        // df = df.Define("V_idx", get_hadronic_gauge_boson_idx, {"GenPart_pdgId", "GenPart_genPartIdxMother"})
        //     .Filter("V_idx == 2 || V_idx == 3");

        // // Boosted vector boson matching
        // df = df.Define("V_eta", "GenPart_eta[V_idx]")
        //     .Define("V_phi", "GenPart_phi[V_idx]")
        //     .Define("V_dR", VdR, {"ak8jet_eta", "ak8jet_phi", "V_eta", "V_phi"})
        //     .Define("V_closest_fatjet", "V_dR.size() > 0 ? ROOT::VecOps::Min(V_dR) : 999.0")
        //     .Define("truth_bv_idx", "(V_closest_fatjet < 0.8) ? (int)ROOT::VecOps::ArgMin(V_dR) : -1");
        
        // // Resolved vector boson matching (quarks from W/Z decay)
        // df = df.Define("truth_q1_eta", "GenPart_eta[abs(GenPart_pdgId) <= 5 && GenPart_genPartIdxMother == V_idx][0]")
        //     .Define("truth_q1_phi", "GenPart_phi[abs(GenPart_pdgId) <= 5 && GenPart_genPartIdxMother == V_idx][0]")
        //     .Define("truth_q2_eta", "GenPart_eta[abs(GenPart_pdgId) <= 5 && GenPart_genPartIdxMother == V_idx][1]")
        //     .Define("truth_q2_phi", "GenPart_phi[abs(GenPart_pdgId) <= 5 && GenPart_genPartIdxMother == V_idx][1]")
        //     .Define("q1_dR", VdR, {"ak4jet_eta", "ak4jet_phi", "truth_q1_eta", "truth_q1_phi"})
        //     .Define("q2_dR", VdR, {"ak4jet_eta", "ak4jet_phi", "truth_q2_eta", "truth_q2_phi"})
        //     .Define("q1_closest_jet", "ROOT::VecOps::Min(q1_dR)")
        //     .Define("q2_closest_jet", "ROOT::VecOps::Min(q2_dR)")
        //     .Define("truth_v1_idx", "(q1_closest_jet < 0.4) ? (int)ROOT::VecOps::ArgMin(q1_dR) : -1")
        //     .Define("truth_v2_idx", "(q2_closest_jet < 0.4) ? (int)ROOT::VecOps::ArgMin(q2_dR) : -1");

        return df;
    }

} // OneLep2FJ