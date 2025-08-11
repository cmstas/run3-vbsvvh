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

    df = GenLevelSelections(df);

    return df;
}















RNode GenLevelSelections(RNode df_) {
    auto df = df_.Define("Higgs_idx", get_higgs_boson_idx, {"GenPart_pdgId", "GenPart_genPartIdxMother"})
            .Define("Higgs_eta", "GenPart_eta[Higgs_idx]")
            .Define("Higgs_phi", "GenPart_phi[Higgs_idx]")
            .Define("VBSJet1_eta", "GenPart_eta[5]")
            .Define("VBSJet1_phi", "GenPart_phi[5]")
            .Define("VBSJet2_eta", "GenPart_eta[6]")
            .Define("VBSJet2_phi", "GenPart_phi[6]")
            .Define("num_v", num_hadronic_gauge_bosons, {"GenPart_pdgId", "GenPart_genPartIdxMother"})
            .Define("V1_idx", "num_v == 1 ? 2 : num_v == 2 ? 3 : num_v == 3 ? 2 : -1")
            .Define("V2_idx", "num_v == 3 ? 3 : -1")
            .Define("V1_eta", "V1_idx != -1 ? GenPart_eta[V1_idx] : -999")
            .Define("V1_phi", "V1_idx != -1 ? GenPart_phi[V1_idx] : -999")
            .Define("V2_eta", "V2_idx != -1 ? GenPart_eta[V2_idx] : -999")
            .Define("V2_phi", "V2_idx != -1 ? GenPart_phi[V2_idx] : -999")
            .Define("qs_from_v1", "abs(GenPart_pdgId) <= 5 && GenPart_genPartIdxMother == V1_idx")
            .Define("qs_from_v2", "abs(GenPart_pdgId) <= 5 && GenPart_genPartIdxMother == V2_idx")
            .Define("V1_q1_eta", "GenPart_eta[qs_from_v1][0]")
            .Define("V1_q1_phi", "GenPart_phi[qs_from_v1][0]")
            .Define("V1_q2_eta", "GenPart_eta[qs_from_v1][1]")
            .Define("V1_q2_phi", "GenPart_phi[qs_from_v1][1]")
            .Define("V2_q1_eta", "GenPart_eta[qs_from_v2][0]")
            .Define("V2_q1_phi", "GenPart_phi[qs_from_v2][0]")
            .Define("V2_q2_eta", "GenPart_eta[qs_from_v2][1]")
            .Define("V2_q2_phi", "GenPart_phi[qs_from_v2][1]")
            .Define("bs_from_higgs", "abs(GenPart_pdgId) == 5 && GenPart_genPartIdxMother == Higgs_idx")
            .Define("b1_eta", "GenPart_eta[bs_from_higgs][0]")
            .Define("b1_phi", "GenPart_phi[bs_from_higgs][0]")
            .Define("b2_eta", "GenPart_eta[bs_from_higgs][1]")
            .Define("b2_phi", "GenPart_phi[bs_from_higgs][1]");

    // VBS matching
    df = df.Filter("Jet_pt.size() > 0")
            .Define("vbs1_dR", get_dR, {"VBSJet1_eta", "VBSJet1_phi", "Jet_eta", "Jet_phi"})
            .Define("vbs2_dR", get_dR, {"VBSJet2_eta", "VBSJet2_phi", "Jet_eta", "Jet_phi"})
            .Define("empty_exclusions", "ROOT::RVec<int>{}")
            .Define("vbs1_idx_temp", find_matching_jet, {"vbs1_dR", "empty_exclusions"})
            .Define("vbs1_exclusions", "ROOT::RVec<int>{vbs1_idx_temp}")
            .Define("vbs2_idx_temp", find_matching_jet, {"vbs2_dR", "vbs1_exclusions"})
            .Define("vbs1_idx", "vbs1_idx_temp >= 0 && vbs1_idx_temp < 10 ? vbs1_idx_temp : -1")
            .Define("vbs2_idx", "vbs2_idx_temp >= 0 && vbs2_idx_temp < 10 ? vbs2_idx_temp : -1")
            .Filter("(vbs1_idx != -1 && vbs2_idx != -1)");

    // Boosted Hbb matching
    df = df.Filter("FatJet_pt.size() > 0")
            .Define("hbb_fatjet_dR", get_dR, {"Higgs_eta", "Higgs_phi", "FatJet_eta", "FatJet_phi"})
            .Define("hbb_dR", "ROOT::VecOps::DeltaR(b1_eta, b1_phi, b2_eta, b2_phi)")
            .Define("empty_fatjet_exclusions", "ROOT::RVec<int>{}")
            .Define("hbb_fatjet_idx_temp", find_matching_fatjet_conditional, {"Higgs_idx", "hbb_fatjet_dR", "empty_fatjet_exclusions"})
            .Define("hbb_fatjet_candidate_b1_dR", "hbb_fatjet_idx_temp >= 0 && hbb_fatjet_idx_temp < FatJet_eta.size() ? ROOT::VecOps::DeltaR(FatJet_eta[hbb_fatjet_idx_temp], FatJet_phi[hbb_fatjet_idx_temp], b1_eta, b1_phi) : 999.0")
            .Define("hbb_fatjet_candidate_b2_dR", "hbb_fatjet_idx_temp >= 0 && hbb_fatjet_idx_temp < FatJet_eta.size() ? ROOT::VecOps::DeltaR(FatJet_eta[hbb_fatjet_idx_temp], FatJet_phi[hbb_fatjet_idx_temp], b2_eta, b2_phi) : 999.0")
            .Define("hbb_isBoosted", "hbb_fatjet_idx_temp != -1 && hbb_fatjet_idx_temp < 3 && hbb_dR < 0.8 && hbb_fatjet_candidate_b1_dR < 0.8 && hbb_fatjet_candidate_b2_dR < 0.8")
            .Define("bH_idx", "hbb_isBoosted ? hbb_fatjet_idx_temp : -1");

    // Resolved Hbb matching
    df = df.Define("b1_jet_dR", get_dR, {"b1_eta", "b1_phi", "Jet_eta", "Jet_phi"})
          .Define("b2_jet_dR", get_dR, {"b2_eta", "b2_phi", "Jet_eta", "Jet_phi"})
          .Define("excluded_jets_for_hbb", "ROOT::RVec<int>{vbs1_idx, vbs2_idx}")
          .Define("rH1_idx_temp", find_matching_jet_conditional, {"Higgs_idx", "b1_jet_dR", "excluded_jets_for_hbb"})
          .Define("excluded_jets_for_hbb2", "ROOT::RVec<int>{vbs1_idx, vbs2_idx, rH1_idx_temp}")
          .Define("Jet_hbb2_idx_temp", find_matching_jet_conditional, {"Higgs_idx", "b2_jet_dR", "excluded_jets_for_hbb2"})
          .Define("rH1_idx", "rH1_idx_temp >= 0 && rH1_idx_temp < 10 ? rH1_idx_temp : -1")
          .Define("rH2_idx", "Jet_hbb2_idx_temp >= 0 && Jet_hbb2_idx_temp < 10 ? Jet_hbb2_idx_temp : -1")
          .Define("hbb_isResolved", "!hbb_isBoosted && rH1_idx != -1 && rH2_idx != -1");

    // Boosted V1qq matching
    df = df.Define("v1qq_fatjet_dR", get_dR_conditional, {"V1_idx", "V1_eta", "V1_phi", "FatJet_eta", "FatJet_phi"})
              .Define("v1qq_dR", "V1_idx != -1 ? ROOT::VecOps::DeltaR(V1_q1_eta, V1_q1_phi, V1_q2_eta, V1_q2_phi) : 999.0")
              .Define("excluded_fatjets_for_v1", "ROOT::RVec<int>{bH_idx}")
              .Define("v1qq_fatjet_idx_temp", find_matching_fatjet_conditional, {"V1_idx", "v1qq_fatjet_dR", "excluded_fatjets_for_v1"})
              .Define("v1qq_fatjet_candidate_q1_dR", "v1qq_fatjet_idx_temp >= 0 && v1qq_fatjet_idx_temp < FatJet_eta.size() ? ROOT::VecOps::DeltaR(FatJet_eta[v1qq_fatjet_idx_temp], FatJet_phi[v1qq_fatjet_idx_temp], V1_q1_eta, V1_q1_phi) : 999.0")
              .Define("v1qq_fatjet_candidate_q2_dR", "v1qq_fatjet_idx_temp >= 0 && v1qq_fatjet_idx_temp < FatJet_eta.size() ? ROOT::VecOps::DeltaR(FatJet_eta[v1qq_fatjet_idx_temp], FatJet_phi[v1qq_fatjet_idx_temp], V1_q2_eta, V1_q2_phi) : 999.0")
              .Define("v1qq_isBoosted", "V1_idx != -1 && v1qq_fatjet_idx_temp != -1 && v1qq_fatjet_idx_temp < 3 && v1qq_dR < 0.8 && v1qq_fatjet_candidate_q1_dR < 0.8 && v1qq_fatjet_candidate_q2_dR < 0.8")
              .Define("bV1_idx", "v1qq_isBoosted ? v1qq_fatjet_idx_temp : -1");

    // Boosted V2qq matching
    df = df.Define("v2qq_fatjet_dR", get_dR_conditional, {"V2_idx", "V2_eta", "V2_phi", "FatJet_eta", "FatJet_phi"})
            .Define("v2qq_dR", "V2_idx != -1 ? ROOT::VecOps::DeltaR(V2_q1_eta, V2_q1_phi, V2_q2_eta, V2_q2_phi) : 999.0")
            .Define("excluded_fatjets_for_v2", "ROOT::RVec<int>{bH_idx, bV1_idx}")
            .Define("v2qq_fatjet_idx_temp", find_matching_fatjet_conditional, {"V2_idx", "v2qq_fatjet_dR", "excluded_fatjets_for_v2"})
            .Define("v2qq_fatjet_candidate_q1_dR", "v2qq_fatjet_idx_temp >= 0 && v2qq_fatjet_idx_temp < FatJet_eta.size() ? ROOT::VecOps::DeltaR(FatJet_eta[v2qq_fatjet_idx_temp], FatJet_phi[v2qq_fatjet_idx_temp], V2_q1_eta, V2_q1_phi) : 999.0")
            .Define("v2qq_fatjet_candidate_q2_dR", "v2qq_fatjet_idx_temp >= 0 && v2qq_fatjet_idx_temp < FatJet_eta.size() ? ROOT::VecOps::DeltaR(FatJet_eta[v2qq_fatjet_idx_temp], FatJet_phi[v2qq_fatjet_idx_temp], V2_q2_eta, V2_q2_phi) : 999.0")
            .Define("v2qq_isBoosted", "V2_idx != -1 && v2qq_fatjet_idx_temp != -1 && v2qq_fatjet_idx_temp < 3 && v2qq_dR < 0.8 && v2qq_fatjet_candidate_q1_dR < 0.8 && v2qq_fatjet_candidate_q2_dR < 0.8")
            .Define("bV2_idx", "v2qq_isBoosted ? v2qq_fatjet_idx_temp : -1");

    // Resolved V1qq matching
    df = df.Define("v1q1_jet_dR", get_dR_conditional, {"V1_idx", "V1_q1_eta", "V1_q1_phi", "Jet_eta", "Jet_phi"})
            .Define("v1q2_jet_dR", get_dR_conditional, {"V1_idx", "V1_q2_eta", "V1_q2_phi", "Jet_eta", "Jet_phi"})
            .Define("excluded_jets_for_v1q1", "ROOT::RVec<int>{vbs1_idx, vbs2_idx, rH1_idx, rH2_idx}")
            .Define("rV1_idx_temp", find_matching_jet_conditional, {"V1_idx", "v1q1_jet_dR", "excluded_jets_for_v1q1"})
            .Define("excluded_jets_for_v1q2", "ROOT::RVec<int>{vbs1_idx, vbs2_idx, rH1_idx, rH2_idx, rV1_idx_temp}")
            .Define("rV2_idx_temp", find_matching_jet_conditional, {"V1_idx", "v1q2_jet_dR", "excluded_jets_for_v1q2"})
            .Define("rV11_idx", "rV1_idx_temp >= 0 && rV1_idx_temp < 10 ? rV1_idx_temp : -1")
            .Define("rV12_idx", "rV2_idx_temp >= 0 && rV2_idx_temp < 10 ? rV2_idx_temp : -1")
            .Define("v1qq_isResolved", "V1_idx != -1 && !v1qq_isBoosted && rV11_idx != -1 && rV12_idx != -1");

    // Resolved V2qq matching
    df = df.Define("v2q1_jet_dR", get_dR_conditional, {"V2_idx", "V2_q1_eta", "V2_q1_phi", "Jet_eta", "Jet_phi"})
            .Define("v2q2_jet_dR", get_dR_conditional, {"V2_idx", "V2_q2_eta", "V2_q2_phi", "Jet_eta", "Jet_phi"})
            .Define("excluded_jets_for_v2q1", "ROOT::RVec<int>{vbs1_idx, vbs2_idx, rH1_idx, rH2_idx, rV11_idx, rV12_idx}")
            .Define("Jet_v2q1_idx_temp", find_matching_jet_conditional, {"V2_idx", "v2q1_jet_dR", "excluded_jets_for_v2q1"})
            .Define("excluded_jets_for_v2q2", "ROOT::RVec<int>{vbs1_idx, vbs2_idx, rH1_idx, rH2_idx, rV11_idx, rV12_idx, Jet_v2q1_idx_temp}")
            .Define("Jet_v2q2_idx_temp", find_matching_jet_conditional, {"V2_idx", "v2q2_jet_dR", "excluded_jets_for_v2q2"})
            .Define("rV21_idx", "Jet_v2q1_idx_temp >= 0 && Jet_v2q1_idx_temp < 10 ? Jet_v2q1_idx_temp : -1")
            .Define("rV22_idx", "Jet_v2q2_idx_temp >= 0 && Jet_v2q2_idx_temp < 10 ? Jet_v2q2_idx_temp : -1")
            .Define("v2qq_isResolved", "V2_idx != -1 && !v2qq_isBoosted && rV21_idx != -1 && rV22_idx != -1");

    return df;
}
