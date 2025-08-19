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
        .Define("nElectron_Tight", "nElectron_Loose == 0 ? 0 : Sum(_tightElectrons)");
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
        .Define("nMuon_Tight", "nMuon_Loose == 0 ? 0 : Sum(_tightMuons)");
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

    return df;
}