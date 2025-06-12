#include "selections_OneLep2FJ.h"

namespace OneLep2FJ {
    TMVA::Experimental::RBDT vbstagger("VBSTagger", "mva/VBSTagger/OneLep2FJ.root");
    SPANet::SPANetInference spanet_inference("/home/users/aaarora/phys/run3/SPANet/spanet.onnx", 512);

    RNode TriggerSelections(RNode df_) {
        return df_.Define("_cut_trigger", "HLT_Ele30_WPTight_Gsf || HLT_IsoMu24");
    }

    RNode LeptonSelections(RNode df_) {
        return df_.Define("_isElectron", "nElectron_Loose == 1 && nElectron_Tight == 1")
            .Define("lepton_pt", "_isElectron ? Electron_pt[0] : Muon_pt[0]")
            .Define("lepton_eta", "_isElectron ? Electron_eta[0] : Muon_eta[0]")
            .Define("lepton_phi", "_isElectron ? Electron_phi[0] : Muon_phi[0]")
            .Define("lepton_mass", "_isElectron ? Electron_mass[0] : Muon_mass[0]");
    }

    RNode AK4JetsSelection(RNode df_) {
        auto df = CommonSelections::AK4JetsSelection(df_);
        return df.Define("_ak4_lep_dR", VdR, {"Jet_eta", "Jet_phi", "lepton_eta", "lepton_phi"})
            .Redefine("Jet_pt", "Jet_pt[_ak4_lep_dR > 0.4]")
            .Redefine("Jet_eta", "Jet_eta[_ak4_lep_dR > 0.4]")
            .Redefine("Jet_phi", "Jet_phi[_ak4_lep_dR > 0.4]")
            .Redefine("Jet_mass", "Jet_mass[_ak4_lep_dR > 0.4]")
            .Redefine("Jet_jetId", "Jet_jetId[_ak4_lep_dR > 0.4]")
            .Redefine("Jet_btagDeepFlavB", "Jet_btagDeepFlavB[_ak4_lep_dR > 0.4]")
            .Redefine("nJet", "Sum(_ak4_lep_dR > 0.4)")
            .Redefine("Jet_ht", "Sum(Jet_pt)")
            .Define("Jet_isTightBTag", "Jet_btagDeepFlavB > 0.7183")
            .Define("Jet_isMediumBTag", "Jet_btagDeepFlavB > 0.3086")
            .Define("Jet_isLooseBTag", "Jet_btagDeepFlavB > 0.0583");
    }
    
    RNode AK8JetsSelection(RNode df_) {
        auto df = CommonSelections::AK8JetsSelection(df_);
        return df.Define("_ak8_lep_dR", VdR, {"FatJet_eta", "FatJet_phi", "lepton_eta", "lepton_phi"})
            .Redefine("FatJet_pt", "FatJet_pt[_ak8_lep_dR > 0.8]")
            .Redefine("FatJet_eta", "FatJet_eta[_ak8_lep_dR > 0.8]")
            .Redefine("FatJet_phi", "FatJet_phi[_ak8_lep_dR > 0.8]")
            .Redefine("FatJet_mass", "FatJet_mass[_ak8_lep_dR > 0.8]")
            .Redefine("FatJet_msoftdrop", "FatJet_msoftdrop[_ak8_lep_dR > 0.8]")
            .Redefine("FatJet_nConstituents", "FatJet_nConstituents[_ak8_lep_dR > 0.8]")
            .Redefine("FatJet_particleNet_XbbVsQCD", "FatJet_particleNet_XbbVsQCD[_ak8_lep_dR > 0.8]")
            .Redefine("FatJet_particleNet_XqqVsQCD", "FatJet_particleNet_XqqVsQCD[_ak8_lep_dR > 0.8]")
            .Redefine("nFatJet", "Sum(_ak8_lep_dR > 0.8)")
            .Redefine("FatJet_ht", "Sum(FatJet_pt)");
    }
        
    RNode runPreselection(RNode df_) {
        // lepton
        auto df = CommonSelections::EventFilters(df_);
        df = TriggerSelections(df);

        df = CommonSelections::ElectronSelections(df);
        df = CommonSelections::MuonSelections(df);
        df = LeptonSelections(df);

        df = df.Define("_cut_lepton", "((nMuon_Loose == 1 && nMuon_Tight == 1 && nElectron_Loose == 0 && nElectron_Tight == 0) || "
            "(nMuon_Loose == 0 && nMuon_Tight == 0 && nElectron_Loose == 1 && nElectron_Tight == 1)) && "
            "(lepton_pt > 40)");

        // SPANET SELECTIONS
        df = AK4JetsSelection(df);
        df = CommonSelections::VBSJetSelections(df, vbstagger);

        df = AK8JetsSelection(df);

        df = spanet_inference.RunSPANetInference(df);
        // spanet output order: h, v, bh, bv, vbs
        
        // pick boosted h and v candidate first (depending on whichever has higher detection score)
        // then check if other boosted candidate is present and doesn't have the same idx as the first one
        // if second candidate is not present, then pick the resolved candidates 

        df = df.Define("vbs1_idx", "spanet_vbs_assignment[0][1]")
            .Define("vbs2_idx", "spanet_vbs_assignment[0][2]")
            .Define("spanet_vbs1_eta", "Jet_eta[vbs1_idx]")
            .Define("spanet_vbs1_phi", "Jet_phi[vbs1_idx]")
            .Define("spanet_vbs2_eta", "Jet_eta[vbs2_idx]")
            .Define("spanet_vbs2_phi", "Jet_phi[vbs2_idx]");

        df = df.Define("bh_idx", "ROOT::VecOps::DeltaR(FatJet_eta[spanet_bh_assignment[0][1]], spanet_vbs1_eta, FatJet_phi[spanet_bh_assignment[0][1]], spanet_vbs1_phi) > 0.8 && "
                "ROOT::VecOps::DeltaR(FatJet_eta[spanet_bh_assignment[0][1]], spanet_vbs2_eta, FatJet_phi[spanet_bh_assignment[0][1]], spanet_vbs2_phi) > 0.8 ? spanet_bh_assignment[0][1] : "
                "spanet_bh_assignment[1][1] == 0 ? -1 : "
                "ROOT::VecOps::DeltaR(FatJet_eta[spanet_bh_assignment[1][1]], spanet_vbs1_eta, FatJet_phi[spanet_bh_assignment[1][1]], spanet_vbs1_phi) > 0.8 && "
                "ROOT::VecOps::DeltaR(FatJet_eta[spanet_bh_assignment[1][1]], spanet_vbs2_eta, FatJet_phi[spanet_bh_assignment[1][1]], spanet_vbs2_phi) > 0.8 ? spanet_bh_assignment[1][1] : "
                "spanet_bh_assignment[2][1] == 0 ? -1 : "
                "ROOT::VecOps::DeltaR(FatJet_eta[spanet_bh_assignment[2][1]], spanet_vbs1_eta, FatJet_phi[spanet_bh_assignment[2][1]], spanet_vbs1_phi) > 0.8 && "
                "ROOT::VecOps::DeltaR(FatJet_eta[spanet_bh_assignment[2][1]], spanet_vbs2_eta, FatJet_phi[spanet_bh_assignment[2][1]], spanet_vbs2_phi) > 0.8 ? spanet_bh_assignment[2][1] : -1")
            .Define("spanet_bh_eta", "bh_idx != -1 ? FatJet_eta[bh_idx] : -999.0f")
            .Define("spanet_bh_phi", "bh_idx != -1 ? FatJet_phi[bh_idx] : -999.0f")
            .Define("bv_idx", "spanet_bv_assignment[0][1] != bh_idx && "
                "ROOT::VecOps::DeltaR(FatJet_eta[spanet_bv_assignment[0][1]], spanet_vbs1_eta, FatJet_phi[spanet_bv_assignment[0][1]], spanet_vbs1_phi) > 0.8 && "
                "ROOT::VecOps::DeltaR(FatJet_eta[spanet_bv_assignment[0][1]], spanet_vbs2_eta, FatJet_phi[spanet_bv_assignment[0][1]], spanet_vbs2_phi) > 0.8 ? spanet_bv_assignment[0][1] : "
                "spanet_bv_assignment[1][1] == 0 ? -1 : "
                "spanet_bv_assignment[1][1] != bh_idx && "
                "ROOT::VecOps::DeltaR(FatJet_eta[spanet_bv_assignment[1][1]], spanet_vbs1_eta, FatJet_phi[spanet_bv_assignment[1][1]], spanet_vbs1_phi) > 0.8 && "
                "ROOT::VecOps::DeltaR(FatJet_eta[spanet_bv_assignment[1][1]], spanet_vbs2_eta, FatJet_phi[spanet_bv_assignment[1][1]], spanet_vbs2_phi) > 0.8 ? spanet_bv_assignment[1][1] : "
                "spanet_bv_assignment[2][1] == 0 ? -1 : "
                "spanet_bv_assignment[2][1] != bh_idx && "
                "ROOT::VecOps::DeltaR(FatJet_eta[spanet_bv_assignment[2][1]], spanet_vbs1_eta, FatJet_phi[spanet_bv_assignment[2][1]], spanet_vbs1_phi) > 0.8 && "
                "ROOT::VecOps::DeltaR(FatJet_eta[spanet_bv_assignment[2][1]], spanet_vbs2_eta, FatJet_phi[spanet_bv_assignment[2][1]], spanet_vbs2_phi) > 0.8 ? spanet_bv_assignment[2][1] : -1")
            .Define("spanet_bv_eta", "bv_idx != -1 ? FatJet_eta[bv_idx] : -999.0f")
            .Define("spanet_bv_phi", "bv_idx != -1 ? FatJet_phi[bv_idx] : -999.0f")
            .Define("h1_idx", "spanet_h_assignment[0][1] != vbs1_idx && spanet_h_assignment[0][1] != vbs2_idx && spanet_h_assignment[0][2] != vbs1_idx && spanet_h_assignment[0][2] != vbs2_idx ? spanet_h_assignment[0][1] : spanet_h_assignment[1][1]")
            .Define("h2_idx", "spanet_h_assignment[0][2] != vbs1_idx && spanet_h_assignment[0][2] != vbs2_idx && spanet_h_assignment[0][1] != vbs1_idx && spanet_h_assignment[0][1] != vbs2_idx ? spanet_h_assignment[0][2] : spanet_h_assignment[1][2]")
            .Define("v1_idx", "spanet_v_assignment[0][1] != vbs1_idx && spanet_v_assignment[0][1] != vbs2_idx && spanet_v_assignment[0][2] != vbs1_idx && spanet_v_assignment[0][2] != vbs2_idx ? spanet_v_assignment[0][1] : spanet_v_assignment[1][1]")
            .Define("v2_idx", "spanet_v_assignment[0][2] != vbs1_idx && spanet_v_assignment[0][2] != vbs2_idx && spanet_v_assignment[0][1] != vbs1_idx && spanet_v_assignment[0][1] != vbs2_idx ? spanet_v_assignment[0][2] : spanet_v_assignment[1][2]");
        
        df = df.Define("isBoostedVH", "nFatJet >= 2 && bh_idx != -1 && bv_idx != -1")
            .Define("isBoostedHResolvedV", "nFatJet == 1 && bh_idx != -1 && bv_idx == -1")
            .Define("isBoostedVResolvedH", "nFatJet == 1 && bv_idx != -1 && bh_idx == -1");

        // OLD SELECTIONS
        df = df.Define("_detajj_vbs_idxs", VBS_MaxEtaJJ, {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass"})
            .Define("old_vbs1_idx", "_detajj_vbs_idxs[0]")
            .Define("old_vbs2_idx", "_detajj_vbs_idxs[1]");
                
        df = df.Define("__ak8_lep_dR", VdR, {"FatJet_eta", "FatJet_phi", "lepton_eta", "lepton_phi"})
            .Define("_ak8_vbs1_dR", VdR, {"FatJet_eta", "FatJet_phi", "vbs1_eta", "vbs1_phi"})
            .Define("_ak8_vbs2_dR", VdR, {"FatJet_eta", "FatJet_phi", "vbs2_eta", "vbs2_phi"})
            .Define("_xbb_candidates", "__ak8_lep_dR > 0.8 && _ak8_vbs1_dR > 0.8 && _ak8_vbs2_dR > 0.8")
            .Define("_xbbvsqcd_score", "FatJet_particleNet_XbbVsQCD[_xbb_candidates]")
            .Define("old_hbb_idx", "_xbbvsqcd_score.size() != 0 ? ArgMax(_xbbvsqcd_score) : -999")
            .Define("hbb_score", "old_hbb_idx != -999 ? _xbbvsqcd_score[old_hbb_idx] : -999.0f")
            .Define("hbb_pt", "old_hbb_idx != -999 ? FatJet_pt[_xbb_candidates][old_hbb_idx] : -999.0f")
            .Define("hbb_eta", "old_hbb_idx != -999 ? FatJet_eta[_xbb_candidates][old_hbb_idx] : -999.0f")
            .Define("hbb_phi", "old_hbb_idx != -999 ? FatJet_phi[_xbb_candidates][old_hbb_idx] : -999.0f")
            .Define("hbb_msoftdrop", "old_hbb_idx != -999 ? FatJet_msoftdrop[_xbb_candidates][old_hbb_idx] : -999.0f");

        df = df.Define("_ak8_hbb_dR", VdR, {"FatJet_eta", "FatJet_phi", "hbb_eta", "hbb_phi"})
            .Define("_xqq_candidates", "__ak8_lep_dR > 0.8 && _ak8_vbs1_dR > 0.8 && _ak8_vbs2_dR > 0.8 && _ak8_hbb_dR > 0.8")
            .Define("_xqqvsqcd_score", "FatJet_particleNet_XqqVsQCD[_xqq_candidates]")
            .Define("old_wqq_idx", "_xqqvsqcd_score.size() != 0 ? ArgMax(_xqqvsqcd_score) : -999")
            .Define("wqq_score", "old_wqq_idx != -999 ? _xqqvsqcd_score[old_wqq_idx] : -999.0f")
            .Define("wqq_pt", "old_wqq_idx != -999 ? FatJet_pt[_xqq_candidates][old_wqq_idx] : -999.0f")
            .Define("wqq_eta", "old_wqq_idx != -999 ? FatJet_eta[_xqq_candidates][old_wqq_idx] : -999.0f")
            .Define("wqq_phi", "old_wqq_idx != -999 ? FatJet_phi[_xqq_candidates][old_wqq_idx] : -999.0f")
            .Define("wqq_msoftdrop", "old_wqq_idx != -999 ? FatJet_msoftdrop[_xqq_candidates][old_wqq_idx] : -999.0f");
        
        // GEN LEVEL MATCHING FOR VALIDATION ONLY
        // VBS jets matching to generator-level partons
        df = df.Define("truth_vbs1_eta", "GenPart_eta[5]")
            .Define("truth_vbs1_phi", "GenPart_phi[5]")
            .Define("truth_vbs2_eta", "GenPart_eta[6]")
            .Define("truth_vbs2_phi", "GenPart_phi[6]")
            .Define("vbs1_dR", VdR, {"Jet_eta", "Jet_phi", "truth_vbs1_eta", "truth_vbs1_phi"})
            .Define("vbs2_dR", VdR, {"Jet_eta", "Jet_phi", "truth_vbs2_eta", "truth_vbs2_phi"})
            .Define("vbs1_closest_jet", "ROOT::VecOps::Min(vbs1_dR)")
            .Define("vbs2_closest_jet", "ROOT::VecOps::Min(vbs2_dR)")
            .Define("truth_vbs1_idx", "(vbs1_closest_jet < 0.4) ? (int)ROOT::VecOps::ArgMin(vbs1_dR) : -1")
            .Define("truth_vbs2_idx", "(vbs2_closest_jet < 0.4) ? (int)ROOT::VecOps::ArgMin(vbs2_dR) : -1");

        // Boosted Higgs matching
        df = df.Define("truth_higgs_eta", "GenPart_eta[4]")
            .Define("truth_higgs_phi", "GenPart_phi[4]")
            .Define("higgs_dR", VdR, {"FatJet_eta", "FatJet_phi", "truth_higgs_eta", "truth_higgs_phi"})
            .Define("higgs_closest_fatjet", "ROOT::VecOps::Min(higgs_dR)")
            .Define("truth_bh_idx", "(higgs_closest_fatjet < 0.8) ? (int)ROOT::VecOps::ArgMin(higgs_dR) : -1");
            
        // Resolved Higgs matching (b-quarks from Higgs decay)
        df = df.Define("truth_b1_eta", "GenPart_eta[abs(GenPart_pdgId) == 5 && GenPart_genPartIdxMother == 4][0]")
            .Define("truth_b1_phi", "GenPart_phi[abs(GenPart_pdgId) == 5 && GenPart_genPartIdxMother == 4][0]")
            .Define("truth_b2_eta", "GenPart_eta[abs(GenPart_pdgId) == 5 && GenPart_genPartIdxMother == 4][1]")
            .Define("truth_b2_phi", "GenPart_phi[abs(GenPart_pdgId) == 5 && GenPart_genPartIdxMother == 4][1]")
            .Define("b1_dR", VdR, {"Jet_eta", "Jet_phi", "truth_b1_eta", "truth_b1_phi"})
            .Define("b2_dR", VdR, {"Jet_eta", "Jet_phi", "truth_b2_eta", "truth_b2_phi"})
            .Define("b1_closest_jet", "ROOT::VecOps::Min(b1_dR)")
            .Define("b2_closest_jet", "ROOT::VecOps::Min(b2_dR)")
            .Define("truth_h1_idx", "(b1_closest_jet < 0.4) ? (int)ROOT::VecOps::ArgMin(b1_dR) : -1")
            .Define("truth_h2_idx", "(b2_closest_jet < 0.4) ? (int)ROOT::VecOps::ArgMin(b2_dR) : -1");

        // Get the hadronic vector boson index (W or Z)
        df = df.Define("V_idx", get_hadronic_gauge_boson_idx, {"GenPart_pdgId", "GenPart_genPartIdxMother"});

        // Boosted vector boson matching
        df = df.Define("V_eta", "GenPart_eta[V_idx]")
            .Define("V_phi", "GenPart_phi[V_idx]")
            .Define("V_dR", VdR, {"FatJet_eta", "FatJet_phi", "V_eta", "V_phi"})
            .Define("V_closest_fatjet", "V_dR.size() > 0 ? ROOT::VecOps::Min(V_dR) : 999.0")
            .Define("truth_bv_idx", "(V_closest_fatjet < 0.8) ? (int)ROOT::VecOps::ArgMin(V_dR) : -1");
        
        // Resolved vector boson matching (quarks from W/Z decay)
        df = df.Define("truth_q1_eta", "GenPart_eta[abs(GenPart_pdgId) <= 5 && GenPart_genPartIdxMother == V_idx][0]")
            .Define("truth_q1_phi", "GenPart_phi[abs(GenPart_pdgId) <= 5 && GenPart_genPartIdxMother == V_idx][0]")
            .Define("truth_q2_eta", "GenPart_eta[abs(GenPart_pdgId) <= 5 && GenPart_genPartIdxMother == V_idx][1]")
            .Define("truth_q2_phi", "GenPart_phi[abs(GenPart_pdgId) <= 5 && GenPart_genPartIdxMother == V_idx][1]")
            .Define("q1_dR", VdR, {"Jet_eta", "Jet_phi", "truth_q1_eta", "truth_q1_phi"})
            .Define("q2_dR", VdR, {"Jet_eta", "Jet_phi", "truth_q2_eta", "truth_q2_phi"})
            .Define("q1_closest_jet", "ROOT::VecOps::Min(q1_dR)")
            .Define("q2_closest_jet", "ROOT::VecOps::Min(q2_dR)")
            .Define("truth_v1_idx", "(q1_closest_jet < 0.4) ? (int)ROOT::VecOps::ArgMin(q1_dR) : -1")
            .Define("truth_v2_idx", "(q2_closest_jet < 0.4) ? (int)ROOT::VecOps::ArgMin(q2_dR) : -1");

        return df;
    }

} // OneLep2FJ