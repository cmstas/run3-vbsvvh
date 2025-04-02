#include "commonSelections.h"

RNode EventFilters(RNode df_) {
    return df_.Filter("Flag_goodVertices && "
            "Flag_globalSuperTightHalo2016Filter && "
            "Flag_EcalDeadCellTriggerPrimitiveFilter && "
            "Flag_BadPFMuonFilter && "
            "Flag_BadPFMuonDzFilter && "
            "Flag_hfNoisyHitsFilter &&"
            "Flag_eeBadScFilter && "
            "Flag_ecalBadCalibFilter");
}

RNode ElectronSelections(RNode df_) {
    return df_.Define("Electron_SC_eta", "Electron_eta + Electron_deltaEtaSC")
        .Define("looseElectrons", 
            "Electron_pt > 7 &&"
            "abs(Electron_SC_eta) < 2.5 && "
            "((abs(Electron_SC_eta) <= 1.479 && abs(Electron_dxy) <= 0.05 && abs(Electron_dz) < 0.1) || (abs(Electron_dxy) <= 0.1 && abs(Electron_dz) < 0.2)) && "
            "abs(Electron_sip3d) < 8 && "
            "Electron_cutBased >= 2 && "
            "Electron_pfRelIso03_all < 0.4 && "
            "Electron_lostHits <= 1")
        .Define("electron_nloose", "nElectron == 0 ? 0 : Sum(looseElectrons)")
        .Define("tightElectrons", "looseElectrons &&" 
            "Electron_pt > 30 && "
            "Electron_cutBased >= 4 && "
            "Electron_pfRelIso03_all < 0.15 && "
            "Electron_hoe < 0.1 && "
            "Electron_eInvMinusPInv > -0.04 && "
            "((abs(Electron_SC_eta) <= 1.479 && Electron_sieie < 0.011) || Electron_sieie <= 0.030) && "
            "Electron_convVeto == true && "
            "Electron_tightCharge == 2 && "
            "Electron_lostHits == 0")
        .Define("electron_ntight", "nElectron == 0 ? 0 : Sum(tightElectrons)")
        .Define("electron_pt", "Electron_pt[tightElectrons]")
        .Define("electron_eta", "Electron_eta[tightElectrons]")
        .Define("electron_sceta", "Electron_SC_eta[tightElectrons]")
        .Define("electron_phi", "Electron_phi[tightElectrons]")
        .Define("electron_mass", "Electron_mass[tightElectrons]")
        .Define("electron_charge", "Electron_charge[tightElectrons]");
}

RNode MuonSelections(RNode df_) {
    return df_.Define("looseMuons", 
            "Muon_pt > 5 && "
            "Muon_pfIsoId >= 2 && "
            "abs(Muon_eta) < 2.4 && "
            "abs(Muon_dxy) < 0.2 && "
            "abs(Muon_dz) < 0.5 && "
            "abs(Muon_sip3d) < 8 && "
            "Muon_looseId == 1")
        .Define("muon_nloose", "nMuon == 0 ? 0 : Sum(looseMuons)")
        .Define("tightMuons", "looseMuons && "
            "Muon_pt > 30 && "
            "Muon_pfIsoId > 4 && "
            "Muon_tightCharge == 2 && "
            "Muon_highPurity && "
            "Muon_tightId")
        .Define("muon_ntight", "nMuon == 0 ? 0 : Sum(tightMuons)")
        .Define("muon_pt", "Muon_pt[tightMuons]")
        .Define("muon_eta", "Muon_eta[tightMuons]")
        .Define("muon_phi", "Muon_phi[tightMuons]")
        .Define("muon_mass", "Muon_mass[tightMuons]")
        .Define("muon_charge", "Muon_charge[tightMuons]");
}

RNode AK8JetsSelection(RNode df_) {
    return df_.Define("goodAK8Jets",
            "CorrFatJet_pt > 250 && "
            "abs(FatJet_eta) <= 2.5 && "
            "FatJet_mass > 50 && "
            "FatJet_msoftdrop > 40 && "
            "FatJet_jetId > 0")
        .Define("ak8jet_pt", "CorrFatJet_pt[goodAK8Jets]")
        .Define("ak8jet_eta", "FatJet_eta[goodAK8Jets]")
        .Define("ak8jet_phi", "FatJet_phi[goodAK8Jets]")
        .Define("ak8jet_mass", "CorrFatJet_mass[goodAK8Jets]")
        .Define("ak8jet_msoftdrop", "FatJet_msoftdrop[goodAK8Jets]")
        .Define("ak8jet_nConstituents", "FatJet_nConstituents[goodAK8Jets]")
        .Define("ak8jet_hbbvsqcd", "FatJet_particleNetWithMass_HbbvsQCD[goodAK8Jets]")
        .Define("ak8jet_wvsqcd", "FatJet_particleNetWithMass_WvsQCD[goodAK8Jets]")
        .Define("ak8jet_zvsqcd", "FatJet_particleNetWithMass_ZvsQCD[goodAK8Jets]")
        .Define("ak8jet_ht", "Sum(ak8jet_pt)")
        .Define("ak8jet_n", "Sum(goodAK8Jets)");
}

RNode AK4JetsSelection(RNode df_) {
    return df_.Define("goodAK4Jets", "Jet_pt > 30 && "
            "abs(Jet_eta) > 1.47 && "
            "abs(Jet_eta) < 4.7 && "
            "Jet_jetId >= 2")
        .Define("ak4jet_pt", "CorrJet_pt[goodAK4Jets]")
        .Define("ak4jet_eta", "Jet_eta[goodAK4Jets]")
        .Define("ak4jet_phi", "Jet_phi[goodAK4Jets]")
        .Define("ak4jet_mass", "CorrJet_mass[goodAK4Jets]")
        .Define("ak4jet_ht", "Sum(ak4jet_pt)")
        .Define("ak4jet_n", "Sum(goodAK4Jets)");
}

RNode VBSJetSelections(RNode df_, TMVA::Experimental::RBDT &vbstagger) {
    return df_.Define("ak4jet_pair_idx", "ROOT::VecOps::Combinations(ak4jet_pt, 2)")
        .Define("jet1_pt", "ROOT::VecOps::Take(ak4jet_pt, ak4jet_pair_idx[0])")
        .Define("jet1_eta", "ROOT::VecOps::Take(ak4jet_eta, ak4jet_pair_idx[0])")
        .Define("jet1_phi", "ROOT::VecOps::Take(ak4jet_phi, ak4jet_pair_idx[0])")
        .Define("jet1_mass", "ROOT::VecOps::Take(ak4jet_mass, ak4jet_pair_idx[0])")
        .Define("jet2_pt", "ROOT::VecOps::Take(ak4jet_pt, ak4jet_pair_idx[1])")
        .Define("jet2_eta", "ROOT::VecOps::Take(ak4jet_eta, ak4jet_pair_idx[1])")
        .Define("jet2_phi", "ROOT::VecOps::Take(ak4jet_phi, ak4jet_pair_idx[1])")
        .Define("jet2_mass", "ROOT::VecOps::Take(ak4jet_mass, ak4jet_pair_idx[1])")
        .Define("pt_jj", "jet1_pt + jet2_pt")
        .Define("deta_jj", "abs(jet1_eta - jet2_eta)")
        .Define("dphi_jj", "ROOT::VecOps::DeltaPhi(jet1_phi, jet2_phi)")
        .Define("m_jj", "ROOT::VecOps::InvariantMasses(jet1_pt, jet1_eta, jet1_phi, jet1_mass, jet2_pt, jet2_eta, jet2_phi, jet2_mass)")
        .Define("VBSTaggerScores", Compute<12, float>(vbstagger), {"jet1_pt", "jet2_pt", "jet1_eta", "jet2_eta", "jet1_phi", "jet2_phi", "jet1_mass", "jet2_mass", "pt_jj", "deta_jj", "dphi_jj", "m_jj"})
        .Define("VBSTagPairIdx", "ArgMax(VBSTaggerScores)")
        .Define("vbs1_idx", "ak4jet_pair_idx[0][VBSTagPairIdx]")
        .Define("vbs2_idx", "ak4jet_pair_idx[1][VBSTagPairIdx]")
        .Define("vbs_score", "VBSTaggerScores[VBSTagPairIdx]")
        .Define("vbs1_pt", "ak4jet_pt[vbs1_idx]")
        .Define("vbs2_pt", "ak4jet_pt[vbs2_idx]")
        .Define("vbs1_eta", "ak4jet_eta[vbs1_idx]")
        .Define("vbs2_eta", "ak4jet_eta[vbs2_idx]")
        .Define("vbs1_phi", "ak4jet_phi[vbs1_idx]")
        .Define("vbs2_phi", "ak4jet_phi[vbs2_idx]")
        .Define("vbs1_mass", "ak4jet_mass[vbs1_idx]")
        .Define("vbs2_mass", "ak4jet_mass[vbs2_idx]");
}