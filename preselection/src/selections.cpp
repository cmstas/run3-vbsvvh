#include "selections.h"

RNode EventFilters(RNode df_) {
    return df_.Define("_cut_filters", "Flag_goodVertices && "
            "Flag_globalSuperTightHalo2016Filter && "
            "Flag_EcalDeadCellTriggerPrimitiveFilter && "
            "Flag_BadPFMuonFilter && "
            "Flag_BadPFMuonDzFilter && "
            "Flag_hfNoisyHitsFilter &&"
            "Flag_eeBadScFilter && "
            "Flag_ecalBadCalibFilter");
}

RNode TriggerSelections(RNode df_, std::string channel, const std::unordered_map<std::string, std::string>& trigger_map) {
    if (trigger_map.empty()) {
        std::cerr << "Warning: No trigger map provided. Skipping trigger selection." << std::endl;
        return df_;
    }
    if (trigger_map.find(channel) == trigger_map.end()) {
        std::cerr << "Warning: Channel '" << channel << "' not found in trigger map. Skipping trigger selection." << std::endl;
        return df_;
    }

    std::string trigger_condition = trigger_map.at(channel);
    return df_.Define("_cut_trigger", trigger_condition);
}

RNode ElectronSelections(RNode df_) {
    auto df = df_.Define("Electron_SC_eta", "Electron_eta + Electron_deltaEtaSC")
        .Define("_looseElectrons", 
            "Electron_pt > 7 &&"
            "abs(Electron_SC_eta) < 2.5 && "
            "((abs(Electron_SC_eta) <= 1.479 && abs(Electron_dxy) <= 0.05 && abs(Electron_dz) < 0.1) || (abs(Electron_dxy) <= 0.1 && abs(Electron_dz) < 0.2)) && "
            "abs(Electron_sip3d) < 8 && "
            "Electron_cutBased >= 2 && "
            "Electron_pfRelIso03_all < 0.4 && "
            "Electron_lostHits <= 1")
        .Define("_tightElectrons", "_looseElectrons &&" 
            "Electron_pt > 30 && "
            "Electron_cutBased >= 4 && "
            "Electron_pfRelIso03_all < 0.15 && "
            "Electron_hoe < 0.1 && "
            "Electron_eInvMinusPInv > -0.04 && "
            "((abs(Electron_SC_eta) <= 1.479 && Electron_sieie < 0.011) || Electron_sieie <= 0.030) && "
            "Electron_convVeto == true && "
            "Electron_tightCharge == 2 && "
            "Electron_lostHits == 0")
        .Define("nElectron_Loose", "nElectron == 0 ? 0 : Sum(_looseElectrons)")
        .Define("nElectron_Tight", "nElectron_Loose == 0 ? 0 : Sum(_tightElectrons)")
        .Define("vvhTightLepMaskElectron", "_tightElectrons");
    return applyObjectMask(df, "_tightElectrons", "Electron");
}

RNode MuonSelections(RNode df_) {
    auto df = df_.Define("_looseMuons", 
            "Muon_pt > 5 && "
            "Muon_pfIsoId >= 2 && "
            "abs(Muon_eta) < 2.4 && "
            "abs(Muon_dxy) < 0.2 && "
            "abs(Muon_dz) < 0.5 && "
            "abs(Muon_sip3d) < 8 && "
            "Muon_looseId == 1")
        .Define("_tightMuons", "_looseMuons && "
            "Muon_pt > 30 && "
            "Muon_pfIsoId > 4 && "
            "Muon_tightCharge == 2 && "
            "Muon_highPurity && "
            "Muon_tightId")
        .Define("nMuon_Loose", "nMuon == 0 ? 0 : Sum(_looseMuons)")
        .Define("nMuon_Tight", "nMuon_Loose == 0 ? 0 : Sum(_tightMuons)")
        .Define("vvhTightLepMaskMuon", "_tightMuons");
    return applyObjectMask(df, "_tightMuons", "Muon");
}

RNode LeptonSelections(RNode df_) {
    auto df = ElectronSelections(df_);
    df = MuonSelections(df);
    return df.Define("Lepton_pt", "Concatenate(Electron_pt, Muon_pt)")
        .Define("Lepton_eta", "Concatenate(Electron_eta, Muon_eta)")
        .Define("Lepton_phi", "Concatenate(Electron_phi, Muon_phi)")
        .Define("Lepton_mass", "Concatenate(Electron_mass, Muon_mass)");
}

RNode AK4JetsSelection(RNode df_) {
    auto df = df_.Define("_dR_ak4_lep", VVdR, {"Jet_eta", "Jet_phi", "Lepton_eta", "Lepton_phi"})
        .Define("_good_ak4jets", " _dR_ak4_lep > 0.4 && "
            "Jet_pt > 20 && "
            "abs(Jet_eta) < 4.7 && "
            "Jet_jetId >= 2")
        .Define("Jet_isTightBTag", "Jet_btagDeepFlavB > 0.6708")
        .Define("Jet_isMediumBTag", "Jet_btagDeepFlavB > 0.2480")
        .Define("Jet_isLooseBTag", "Jet_btagDeepFlavB > 0.0485");
        // .Define("Jet_isTightBTag", "Jet_btagDeepFlavB > 0.7183")
        // .Define("Jet_isMediumBTag", "Jet_btagDeepFlavB > 0.3086")
        // .Define("Jet_isLooseBTag", "Jet_btagDeepFlavB > 0.0583");
    df = applyObjectMask(df, "_good_ak4jets", "Jet");
    return df;
}

RNode AK8JetsSelection(RNode df_) {
    auto df = df_.Define("_dR_ak8_lep", VVdR, {"FatJet_eta", "FatJet_phi", "Lepton_eta", "Lepton_phi"})
        .Define("_good_ak8jets", "_dR_ak8_lep > 0.4 && "
            "FatJet_pt > 250 && "
            "abs(FatJet_eta) <= 2.5 && "
            "FatJet_mass > 50 && "
            "FatJet_msoftdrop > 40 && "
            "FatJet_jetId > 0");
    df = applyObjectMask(df, "_good_ak8jets", "FatJet");
    return df;
}

RNode runPreselection(RNode df_, std::string channel, SPANet::SPANetInference &spanet_session) {
    auto df = EventFilters(df_);
    df = TriggerSelections(df, channel, TriggerMap);
    df = LeptonSelections(df);
    df = AK4JetsSelection(df);
    df = AK8JetsSelection(df);

    // channel-specific selections
    if (channel == "1Lep2FJ") {
        df = df.Define("_cut_lepton", "((nMuon_Loose == 1 && nMuon_Tight == 1 && nElectron_Loose == 0 && nElectron_Tight == 0) || "
            "(nMuon_Loose == 0 && nMuon_Tight == 0 && nElectron_Loose == 1 && nElectron_Tight == 1)) && "
            "(Lepton_pt[0] > 40)");
    }

    df = spanet_session.RunSPANetInference(df);
    df = ParseSpanet(df);

    df = GenLevelSelections(df);

    return df;
}

RNode ParseSpanet(RNode df_){
    auto df = df_.Define("_all_assignments", assign_all_objects, {
        "spanet_vbs_assignment", "spanet_h_assignment", "spanet_bh_assignment", 
        "spanet_v1_assignment", "spanet_v2_assignment", "spanet_bv1_assignment", "spanet_bv2_assignment",
        "spanet_vbs_detection", "spanet_h_detection", "spanet_bh_detection", 
        "spanet_v1_detection", "spanet_v2_detection", "spanet_bv1_detection", "spanet_bv2_detection",
        "Jet_eta", "Jet_phi", "FatJet_eta", "FatJet_phi"
    })
        .Define("vbs1_idx", "_all_assignments[0]")
        .Define("vbs2_idx", "_all_assignments[1]")
        .Define("h1_idx", "_all_assignments[2]")
        .Define("h2_idx", "_all_assignments[3]")
        .Define("bh_idx", "_all_assignments[4]")
        .Define("v1_j1_idx", "_all_assignments[5]")
        .Define("v1_j2_idx", "_all_assignments[6]")
        .Define("v2_j1_idx", "_all_assignments[7]")
        .Define("v2_j2_idx", "_all_assignments[8]")
        .Define("bv1_idx", "_all_assignments[9]")
        .Define("bv2_idx", "_all_assignments[10]");

    // VBS jet variables
    df = df.Define("spanet_vbs1_pt", "vbs1_idx >= 0 ? Jet_pt[vbs1_idx] : -999.0f")
            .Define("spanet_vbs1_eta", "vbs1_idx >= 0 ? Jet_eta[vbs1_idx] : -999.0f")
            .Define("spanet_vbs1_phi", "vbs1_idx >= 0 ? Jet_phi[vbs1_idx] : -999.0f")
            .Define("spanet_vbs1_mass", "vbs1_idx >= 0 ? Jet_mass[vbs1_idx] : -999.0f")
            .Define("spanet_vbs2_pt", "vbs2_idx >= 0 ? Jet_pt[vbs2_idx] : -999.0f")
            .Define("spanet_vbs2_eta", "vbs2_idx >= 0 ? Jet_eta[vbs2_idx] : -999.0f")
            .Define("spanet_vbs2_phi", "vbs2_idx >= 0 ? Jet_phi[vbs2_idx] : -999.0f")
            .Define("spanet_vbs2_mass", "vbs2_idx >= 0 ? Jet_mass[vbs2_idx] : -999.0f")
            .Define("spanet_vbs_detajj", "vbs1_idx >= 0 && vbs2_idx >= 0 ? abs(spanet_vbs1_eta - spanet_vbs2_eta) : -999.0f")
            .Define("spanet_vbs_mjj", "vbs1_idx >= 0 && vbs2_idx >= 0 ? (ROOT::Math::PtEtaPhiMVector(spanet_vbs1_pt, spanet_vbs1_eta, spanet_vbs1_phi, spanet_vbs1_mass) + "
                "ROOT::Math::PtEtaPhiMVector(spanet_vbs2_pt, spanet_vbs2_eta, spanet_vbs2_phi, spanet_vbs2_mass)).M() : -999.0f");

    // Resolved Higgs variables
    df = df.Define("spanet_h1_pt", "h1_idx >= 0 ? Jet_pt[h1_idx] : -999.0f")
            .Define("spanet_h1_eta", "h1_idx >= 0 ? Jet_eta[h1_idx] : -999.0f")
            .Define("spanet_h1_phi", "h1_idx >= 0 ? Jet_phi[h1_idx] : -999.0f")
            .Define("spanet_h1_mass", "h1_idx >= 0 ? Jet_mass[h1_idx] : -999.0f")
            .Define("spanet_h2_pt", "h2_idx >= 0 ? Jet_pt[h2_idx] : -999.0f")
            .Define("spanet_h2_eta", "h2_idx >= 0 ? Jet_eta[h2_idx] : -999.0f")
            .Define("spanet_h2_phi", "h2_idx >= 0 ? Jet_phi[h2_idx] : -999.0f")
            .Define("spanet_h2_mass", "h2_idx >= 0 ? Jet_mass[h2_idx] : -999.0f")
            .Define("spanet_h_mjj", "h1_idx >= 0 && h2_idx >= 0 ? (ROOT::Math::PtEtaPhiMVector(spanet_h1_pt, spanet_h1_eta, spanet_h1_phi, spanet_h1_mass) + "
                "ROOT::Math::PtEtaPhiMVector(spanet_h2_pt, spanet_h2_eta, spanet_h2_phi, spanet_h2_mass)).M() : -999.0f");

    // Boosted Higgs variables
    df = df.Define("spanet_bh_eta", "bh_idx >= 0 ? FatJet_eta[bh_idx] : -999.0f")
            .Define("spanet_bh_phi", "bh_idx >= 0 ? FatJet_phi[bh_idx] : -999.0f")
            .Define("spanet_bh_msoftdrop", "bh_idx >= 0 ? FatJet_msoftdrop[bh_idx] : -999.0f")
            .Define("spanet_bh_pt", "bh_idx >= 0 ? FatJet_pt[bh_idx] : -999.0f")
            .Define("spanet_bh_score", "bh_idx >= 0 ? FatJet_globalParT3_Xbb[bh_idx] / (FatJet_globalParT3_Xbb[bh_idx] + FatJet_globalParT3_QCD[bh_idx]) : -999.0f");

    // Resolved V1 variables
    df = df.Define("spanet_v1_j1_pt", "v1_j1_idx >= 0 ? Jet_pt[v1_j1_idx] : -999.0f")
            .Define("spanet_v1_j1_eta", "v1_j1_idx >= 0 ? Jet_eta[v1_j1_idx] : -999.0f")
            .Define("spanet_v1_j1_phi", "v1_j1_idx >= 0 ? Jet_phi[v1_j1_idx] : -999.0f")
            .Define("spanet_v1_j1_mass", "v1_j1_idx >= 0 ? Jet_mass[v1_j1_idx] : -999.0f")
            .Define("spanet_v1_j2_pt", "v1_j2_idx >= 0 ? Jet_pt[v1_j2_idx] : -999.0f")
            .Define("spanet_v1_j2_eta", "v1_j2_idx >= 0 ? Jet_eta[v1_j2_idx] : -999.0f")
            .Define("spanet_v1_j2_phi", "v1_j2_idx >= 0 ? Jet_phi[v1_j2_idx] : -999.0f")
            .Define("spanet_v1_j2_mass", "v1_j2_idx >= 0 ? Jet_mass[v1_j2_idx] : -999.0f")
            .Define("spanet_v1_mjj", "v1_j1_idx >= 0 && v1_j2_idx >= 0 ? (ROOT::Math::PtEtaPhiMVector(spanet_v1_j1_pt, spanet_v1_j1_eta, spanet_v1_j1_phi, spanet_v1_j1_mass) + "
                "ROOT::Math::PtEtaPhiMVector(spanet_v1_j2_pt, spanet_v1_j2_eta, spanet_v1_j2_phi, spanet_v1_j2_mass)).M() : -999.0f");

    // Resolved V2 variables
    df = df.Define("spanet_v2_j1_pt", "v2_j1_idx >= 0 ? Jet_pt[v2_j1_idx] : -999.0f")
            .Define("spanet_v2_j1_eta", "v2_j1_idx >= 0 ? Jet_eta[v2_j1_idx] : -999.0f")
            .Define("spanet_v2_j1_phi", "v2_j1_idx >= 0 ? Jet_phi[v2_j1_idx] : -999.0f")
            .Define("spanet_v2_j1_mass", "v2_j1_idx >= 0 ? Jet_mass[v2_j1_idx] : -999.0f")
            .Define("spanet_v2_j2_pt", "v2_j2_idx >= 0 ? Jet_pt[v2_j2_idx] : -999.0f")
            .Define("spanet_v2_j2_eta", "v2_j2_idx >= 0 ? Jet_eta[v2_j2_idx] : -999.0f")
            .Define("spanet_v2_j2_phi", "v2_j2_idx >= 0 ? Jet_phi[v2_j2_idx] : -999.0f")
            .Define("spanet_v2_j2_mass", "v2_j2_idx >= 0 ? Jet_mass[v2_j2_idx] : -999.0f")
            .Define("spanet_v2_mjj", "v2_j1_idx >= 0 && v2_j2_idx >= 0 ? (ROOT::Math::PtEtaPhiMVector(spanet_v2_j1_pt, spanet_v2_j1_eta, spanet_v2_j1_phi, spanet_v2_j1_mass) + "
                "ROOT::Math::PtEtaPhiMVector(spanet_v2_j2_pt, spanet_v2_j2_eta, spanet_v2_j2_phi, spanet_v2_j2_mass)).M() : -999.0f");

    // Boosted V1 variables
    df = df.Define("spanet_bv1_eta", "bv1_idx >= 0 ? FatJet_eta[bv1_idx] : -999.0f")
            .Define("spanet_bv1_phi", "bv1_idx >= 0 ? FatJet_phi[bv1_idx] : -999.0f")
            .Define("spanet_bv1_msoftdrop", "bv1_idx >= 0 ? FatJet_msoftdrop[bv1_idx] : -999.0f")
            .Define("spanet_bv1_pt", "bv1_idx >= 0 ? FatJet_pt[bv1_idx] : -999.0f")
            .Define("spanet_bv1_score", "bv1_idx >= 0 ? FatJet_particleNet_XqqVsQCD[bv1_idx] : -999.0f")
            .Define("spanet_bv1_w_score", "bv1_idx >= 0 ? FatJet_globalParT3_Xqq[bv1_idx] / (FatJet_globalParT3_Xqq[bv1_idx] + FatJet_globalParT3_Xcs[bv1_idx] + FatJet_globalParT3_QCD[bv1_idx]) : -999.0f");

    // Boosted V2 final_variables
    df = df.Define("spanet_bv2_eta", "bv2_idx >= 0 ? FatJet_eta[bv2_idx] : -999.0f")
            .Define("spanet_bv2_phi", "bv2_idx >= 0 ? FatJet_phi[bv2_idx] : -999.0f")
            .Define("spanet_bv2_msoftdrop", "bv2_idx >= 0 ? FatJet_msoftdrop[bv2_idx] : -999.0f")
            .Define("spanet_bv2_pt", "bv2_idx >= 0 ? FatJet_pt[bv2_idx] : -999.0f")
            .Define("spanet_bv2_score", "bv2_idx >= 0 ? FatJet_particleNet_XqqVsQCD[bv2_idx] : -999.0f")
            .Define("spanet_bv2_w_score", "bv2_idx >= 0 ? FatJet_globalParT3_Xqq[bv2_idx] / (FatJet_globalParT3_Xqq[bv2_idx] + FatJet_globalParT3_Xcs[bv2_idx] + FatJet_globalParT3_QCD[bv2_idx]) : -999.0f");

    return df;
}

RNode GenLevelSelections(RNode df_) {
    auto df = df_.Define("higgs_info", get_higgs_boson_idx, {"GenPart_pdgId", "GenPart_status", "GenPart_genPartIdxMother"})
            .Define("Higgs_idx", "higgs_info[0]")
            .Define("Higgs_eta", "Higgs_idx >= 0 ? GenPart_eta[Higgs_idx] : -999.0f")
            .Define("Higgs_phi", "Higgs_idx >= 0 ? GenPart_phi[Higgs_idx] : -999.0f")
            .Define("b1_eta", "higgs_info[1] >= 0 ? GenPart_eta[higgs_info[1]] : -999.0f")
            .Define("b1_phi", "higgs_info[1] >= 0 ? GenPart_phi[higgs_info[1]] : -999.0f")
            .Define("b2_eta", "higgs_info[2] >= 0 ? GenPart_eta[higgs_info[2]] : -999.0f")
            .Define("b2_phi", "higgs_info[2] >= 0 ? GenPart_phi[higgs_info[2]] : -999.0f");

    df = df.Define("vboson_info", get_v_boson_idx, {"GenPart_pdgId", "GenPart_status", "GenPart_genPartIdxMother"})
            .Define("V1_idx", "vboson_info[0]")
            .Define("V1_eta", "V1_idx >= 0 ? GenPart_eta[V1_idx] : -999.0f")
            .Define("V1_phi", "V1_idx >= 0 ? GenPart_phi[V1_idx] : -999.0f")
            .Define("V1_q1_eta", "vboson_info[1] >= 0 ? GenPart_eta[vboson_info[1]] : -999.0f")
            .Define("V1_q1_phi", "vboson_info[1] >= 0 ? GenPart_phi[vboson_info[1]] : -999.0f")
            .Define("V1_q2_eta", "vboson_info[2] >= 0 ? GenPart_eta[vboson_info[2]] : -999.0f")
            .Define("V1_q2_phi", "vboson_info[2] >= 0 ? GenPart_phi[vboson_info[2]] : -999.0f")
            .Define("V2_idx", "vboson_info[3]")
            .Define("V2_eta", "V2_idx >= 0 ? GenPart_eta[V2_idx] : -999.0f")
            .Define("V2_phi", "V2_idx >= 0 ? GenPart_phi[V2_idx] : -999.0f")
            .Define("V2_q1_eta", "vboson_info[4] >= 0 ? GenPart_eta[vboson_info[4]] : -999.0f")
            .Define("V2_q1_phi", "vboson_info[4] >= 0 ? GenPart_phi[vboson_info[4]] : -999.0f")
            .Define("V2_q2_eta", "vboson_info[5] >= 0 ? GenPart_eta[vboson_info[5]] : -999.0f")
            .Define("V2_q2_phi", "vboson_info[5] >= 0 ? GenPart_phi[vboson_info[5]] : -999.0f");

    df = df.Define("vbs_info", get_vbs_quarks_idxs, {"GenPart_pdgId", "GenPart_status", "GenPart_genPartIdxMother"})
            .Define("VBSJet1_eta", "vbs_info[0] >= 0 ? GenPart_eta[vbs_info[0]] : -999.0f")
            .Define("VBSJet1_phi", "vbs_info[0] >= 0 ? GenPart_phi[vbs_info[0]] : -999.0f")
            .Define("VBSJet2_eta", "vbs_info[1] >= 0 ? GenPart_eta[vbs_info[1]] : -999.0f")
            .Define("VBSJet2_phi", "vbs_info[1] >= 0 ? GenPart_phi[vbs_info[1]] : -999.0f");

    df = df.Define("vbs1_dR", get_dR, {"VBSJet1_eta", "VBSJet1_phi", "Jet_eta", "Jet_phi"})
            .Define("vbs2_dR", get_dR, {"VBSJet2_eta", "VBSJet2_phi", "Jet_eta", "Jet_phi"})
            .Define("empty_exclusions", "ROOT::RVec<int>{}")
            .Define("vbs1_idx_temp", find_matching_jet, {"vbs1_dR", "empty_exclusions"})
            .Define("vbs1_exclusions", "ROOT::RVec<int>{vbs1_idx_temp}")
            .Define("vbs2_idx_temp", find_matching_jet, {"vbs2_dR", "vbs1_exclusions"})
            .Define("gen_vbs1_idx", "vbs1_idx_temp >= 0 && vbs1_idx_temp < 10 ? vbs1_idx_temp : -1")
            .Define("gen_vbs2_idx", "vbs2_idx_temp >= 0 && vbs2_idx_temp < 10 ? vbs2_idx_temp : -1")
            .Filter("(gen_vbs1_idx != -1 && gen_vbs2_idx != -1)");

    df = df.Define("hbb_fatjet_dR", get_dR, {"Higgs_eta", "Higgs_phi", "FatJet_eta", "FatJet_phi"})
            .Define("hbb_dR", "ROOT::VecOps::DeltaR(b1_eta, b1_phi, b2_eta, b2_phi)")
            .Define("empty_fatjet_exclusions", "ROOT::RVec<int>{}")
            .Define("hbb_fatjet_idx_temp", find_matching_fatjet_conditional, {"Higgs_idx", "hbb_fatjet_dR", "empty_fatjet_exclusions"})
            .Define("hbb_fatjet_candidate_b1_dR", "hbb_fatjet_idx_temp >= 0 && hbb_fatjet_idx_temp < FatJet_eta.size() ? ROOT::VecOps::DeltaR(FatJet_eta[hbb_fatjet_idx_temp], FatJet_phi[hbb_fatjet_idx_temp], b1_eta, b1_phi) : 999.0")
            .Define("hbb_fatjet_candidate_b2_dR", "hbb_fatjet_idx_temp >= 0 && hbb_fatjet_idx_temp < FatJet_eta.size() ? ROOT::VecOps::DeltaR(FatJet_eta[hbb_fatjet_idx_temp], FatJet_phi[hbb_fatjet_idx_temp], b2_eta, b2_phi) : 999.0")
            .Define("hbb_isBoosted", "hbb_fatjet_idx_temp != -1 && hbb_fatjet_idx_temp < 3 && hbb_dR < 0.8 && hbb_fatjet_candidate_b1_dR < 0.8 && hbb_fatjet_candidate_b2_dR < 0.8")
            .Define("gen_bh_idx", "hbb_isBoosted ? hbb_fatjet_idx_temp : -1");

    df = df.Define("b1_jet_dR", get_dR, {"b1_eta", "b1_phi", "Jet_eta", "Jet_phi"})
            .Define("b2_jet_dR", get_dR, {"b2_eta", "b2_phi", "Jet_eta", "Jet_phi"})
            .Define("excluded_jets_for_hbb", "ROOT::RVec<int>{gen_vbs1_idx, gen_vbs2_idx}")
            .Define("h1_idx_temp", find_matching_jet_conditional, {"Higgs_idx", "b1_jet_dR", "excluded_jets_for_hbb"})
            .Define("excluded_jets_for_hbb2", "ROOT::RVec<int>{gen_vbs1_idx, gen_vbs2_idx, h1_idx_temp}")
            .Define("tempJet_hbb2_idx", find_matching_jet_conditional, {"Higgs_idx", "b2_jet_dR", "excluded_jets_for_hbb2"})
            .Define("gen_h1_idx", "h1_idx_temp >= 0 && h1_idx_temp < 10 ? h1_idx_temp : -1")
            .Define("gen_h2_idx", "tempJet_hbb2_idx >= 0 && tempJet_hbb2_idx < 10 ? tempJet_hbb2_idx : -1")
            .Define("hbb_isResolved", "!hbb_isBoosted && gen_h1_idx != -1 && gen_h2_idx != -1");

    df = df.Define("v1qq_fatjet_dR", get_dR_conditional, {"V1_idx", "V1_eta", "V1_phi", "FatJet_eta", "FatJet_phi"})
            .Define("v1qq_dR", "V1_idx != -1 ? ROOT::VecOps::DeltaR(V1_q1_eta, V1_q1_phi, V1_q2_eta, V1_q2_phi) : 999.0")
            .Define("excluded_fatjets_for_v1", "ROOT::RVec<int>{gen_bh_idx}")
            .Define("v1qq_fatjet_idx_temp", find_matching_fatjet_conditional, {"V1_idx", "v1qq_fatjet_dR", "excluded_fatjets_for_v1"})
            .Define("v1qq_fatjet_candidate_q1_dR", "v1qq_fatjet_idx_temp >= 0 && v1qq_fatjet_idx_temp < FatJet_eta.size() ? ROOT::VecOps::DeltaR(FatJet_eta[v1qq_fatjet_idx_temp], FatJet_phi[v1qq_fatjet_idx_temp], V1_q1_eta, V1_q1_phi) : 999.0")
            .Define("v1qq_fatjet_candidate_q2_dR", "v1qq_fatjet_idx_temp >= 0 && v1qq_fatjet_idx_temp < FatJet_eta.size() ? ROOT::VecOps::DeltaR(FatJet_eta[v1qq_fatjet_idx_temp], FatJet_phi[v1qq_fatjet_idx_temp], V1_q2_eta, V1_q2_phi) : 999.0")
            .Define("v1qq_isBoosted", "V1_idx != -1 && v1qq_fatjet_idx_temp != -1 && v1qq_fatjet_idx_temp < 3 && v1qq_dR < 0.8 && v1qq_fatjet_candidate_q1_dR < 0.8 && v1qq_fatjet_candidate_q2_dR < 0.8")
            .Define("gen_bv1_idx", "v1qq_isBoosted ? v1qq_fatjet_idx_temp : -1");

    df = df.Define("v2qq_fatjet_dR", get_dR_conditional, {"V2_idx", "V2_eta", "V2_phi", "FatJet_eta", "FatJet_phi"})
            .Define("v2qq_dR", "V2_idx != -1 ? ROOT::VecOps::DeltaR(V2_q1_eta, V2_q1_phi, V2_q2_eta, V2_q2_phi) : 999.0")
            .Define("excluded_fatjets_for_v2", "ROOT::RVec<int>{gen_bh_idx, gen_bv1_idx}")
            .Define("v2qq_fatjet_idx_temp", find_matching_fatjet_conditional, {"V2_idx", "v2qq_fatjet_dR", "excluded_fatjets_for_v2"})
            .Define("v2qq_fatjet_candidate_q1_dR", "v2qq_fatjet_idx_temp >= 0 && v2qq_fatjet_idx_temp < FatJet_eta.size() ? ROOT::VecOps::DeltaR(FatJet_eta[v2qq_fatjet_idx_temp], FatJet_phi[v2qq_fatjet_idx_temp], V2_q1_eta, V2_q1_phi) : 999.0")
            .Define("v2qq_fatjet_candidate_q2_dR", "v2qq_fatjet_idx_temp >= 0 && v2qq_fatjet_idx_temp < FatJet_eta.size() ? ROOT::VecOps::DeltaR(FatJet_eta[v2qq_fatjet_idx_temp], FatJet_phi[v2qq_fatjet_idx_temp], V2_q2_eta, V2_q2_phi) : 999.0")
            .Define("v2qq_isBoosted", "V2_idx != -1 && v2qq_fatjet_idx_temp != -1 && v2qq_fatjet_idx_temp < 3 && v2qq_dR < 0.8 && v2qq_fatjet_candidate_q1_dR < 0.8 && v2qq_fatjet_candidate_q2_dR < 0.8")
            .Define("gen_bv2_idx", "v2qq_isBoosted ? v2qq_fatjet_idx_temp : -1");

    df = df.Define("v1q1_jet_dR", get_dR_conditional, {"V1_idx", "V1_q1_eta", "V1_q1_phi", "Jet_eta", "Jet_phi"})
            .Define("v1q2_jet_dR", get_dR_conditional, {"V1_idx", "V1_q2_eta", "V1_q2_phi", "Jet_eta", "Jet_phi"})
            .Define("excluded_jets_for_v1q1", "ROOT::RVec<int>{gen_vbs1_idx, gen_vbs2_idx, gen_h1_idx, gen_h2_idx}")
            .Define("rV1_idx_temp", find_matching_jet_conditional, {"V1_idx", "v1q1_jet_dR", "excluded_jets_for_v1q1"})
            .Define("excluded_jets_for_v1q2", "ROOT::RVec<int>{gen_vbs1_idx, gen_vbs2_idx, gen_h1_idx, gen_h2_idx, rV1_idx_temp}")
            .Define("rV2_idx_temp", find_matching_jet_conditional, {"V1_idx", "v1q2_jet_dR", "excluded_jets_for_v1q2"})
            .Define("gen_v1_j1_idx", "rV1_idx_temp >= 0 && rV1_idx_temp < 10 ? rV1_idx_temp : -1")
            .Define("gen_v1_j2_idx", "rV2_idx_temp >= 0 && rV2_idx_temp < 10 ? rV2_idx_temp : -1")
            .Define("v1qq_isResolved", "V1_idx != -1 && !v1qq_isBoosted && gen_v1_j1_idx != -1 && gen_v1_j2_idx != -1");

    df = df.Define("v2q1_jet_dR", get_dR_conditional, {"V2_idx", "V2_q1_eta", "V2_q1_phi", "Jet_eta", "Jet_phi"})
            .Define("v2q2_jet_dR", get_dR_conditional, {"V2_idx", "V2_q2_eta", "V2_q2_phi", "Jet_eta", "Jet_phi"})
            .Define("excluded_jets_for_v2q1", "ROOT::RVec<int>{gen_vbs1_idx, gen_vbs2_idx, gen_h1_idx, gen_h2_idx, gen_v1_j1_idx, gen_v1_j2_idx}")
            .Define("tempJet_v2q1_idx", find_matching_jet_conditional, {"V2_idx", "v2q1_jet_dR", "excluded_jets_for_v2q1"})
            .Define("excluded_jets_for_v2q2", "ROOT::RVec<int>{gen_vbs1_idx, gen_vbs2_idx, gen_h1_idx, gen_h2_idx, gen_v1_j1_idx, gen_v1_j2_idx, tempJet_v2q1_idx}")
            .Define("tJet_v2q2_idx_temp", find_matching_jet_conditional, {"V2_idx", "v2q2_jet_dR", "excluded_jets_for_v2q2"})
            .Define("gen_v2_j1_idx", "tempJet_v2q1_idx >= 0 && tempJet_v2q1_idx < 10 ? tempJet_v2q1_idx : -1")
            .Define("gen_v2_j2_idx", "tJet_v2q2_idx_temp >= 0 && tJet_v2q2_idx_temp < 10 ? tJet_v2q2_idx_temp : -1")
            .Define("v2qq_isResolved", "V2_idx != -1 && !v2qq_isBoosted && gen_v2_j1_idx != -1 && gen_v2_j2_idx != -1");

    return df;
}
