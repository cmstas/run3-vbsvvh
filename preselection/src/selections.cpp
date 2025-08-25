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
            "Muon_looseId")
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
        .Define("_good_ak8jets", "_dR_ak8_lep > 0.8 && "
            "FatJet_pt > 250 && "
            "abs(FatJet_eta) <= 2.5 && "
            "FatJet_msoftdrop > 40 && "
            "FatJet_jetId > 0");
    df = applyObjectMask(df, "_good_ak8jets", "FatJet");
    return df;
}

RNode runPreselection(RNode df_, std::string channel) {
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
    else if (channel == "0Lep3FJ") {
        df = df.Define("_cut_lepton", "nMuon_Loose == 0 && nElectron_Loose == 0");
    }
    return df;
}