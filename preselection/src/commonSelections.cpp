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
    return df_.Define("_Electron_SC_eta", "Electron_eta + Electron_deltaEtaSC")
        .Define("_looseElectrons", 
            "Electron_pt > 7 &&"
            "abs(_Electron_SC_eta) < 2.5 && "
            "((abs(_Electron_SC_eta) <= 1.479 && abs(Electron_dxy) <= 0.05 && abs(Electron_dz) < 0.1) || (abs(Electron_dxy) <= 0.1 && abs(Electron_dz) < 0.2)) && "
            "abs(Electron_sip3d) < 8 && "
            "Electron_cutBased >= 2 && "
            "Electron_pfRelIso03_all < 0.4 && "
            "Electron_lostHits <= 1")
        .Define("electron_nloose", "nElectron == 0 ? 0 : Sum(_looseElectrons)")
        .Define("_tightElectrons", "_looseElectrons &&" 
            "Electron_pt > 30 && "
            "Electron_cutBased >= 4 && "
            "Electron_pfRelIso03_all < 0.15 && "
            "Electron_hoe < 0.1 && "
            "Electron_eInvMinusPInv > -0.04 && "
            "((abs(_Electron_SC_eta) <= 1.479 && Electron_sieie < 0.011) || Electron_sieie <= 0.030) && "
            "Electron_convVeto == true && "
            "Electron_tightCharge == 2 && "
            "Electron_lostHits == 0")
        .Define("electron_ntight", "nElectron == 0 ? 0 : Sum(_tightElectrons)")
        .Define("electron_pt", "Electron_pt[_tightElectrons]")
        .Define("electron_eta", "Electron_eta[_tightElectrons]")
        .Define("electron_sceta", "_Electron_SC_eta[_tightElectrons]")
        .Define("electron_phi", "Electron_phi[_tightElectrons]")
        .Define("electron_mass", "Electron_mass[_tightElectrons]")
        .Define("electron_charge", "Electron_charge[_tightElectrons]");
}

RNode MuonSelections(RNode df_) {
    return df_.Define("_looseMuons", 
            "Muon_pt > 5 && "
            "Muon_pfIsoId >= 2 && "
            "abs(Muon_eta) < 2.4 && "
            "abs(Muon_dxy) < 0.2 && "
            "abs(Muon_dz) < 0.5 && "
            "abs(Muon_sip3d) < 8 && "
            "Muon_looseId == 1")
        .Define("muon_nloose", "nMuon == 0 ? 0 : Sum(_looseMuons)")
        .Define("_tightMuons", "_looseMuons && "
            "Muon_pt > 30 && "
            "Muon_pfIsoId > 4 && "
            "Muon_tightCharge == 2 && "
            "Muon_highPurity && "
            "Muon_tightId")
        .Define("muon_ntight", "nMuon == 0 ? 0 : Sum(_tightMuons)")
        .Define("muon_pt", "Muon_pt[_tightMuons]")
        .Define("muon_eta", "Muon_eta[_tightMuons]")
        .Define("muon_phi", "Muon_phi[_tightMuons]")
        .Define("muon_mass", "Muon_mass[_tightMuons]")
        .Define("muon_charge", "Muon_charge[_tightMuons]");
}

RNode AK8JetsSelection(RNode df_) {
    return df_.Define("_goodAK8Jets",
            "FatJet_pt > 250 && "
            "abs(FatJet_eta) <= 2.5 && "
            "FatJet_mass > 50 && "
            "FatJet_msoftdrop > 40 && "
            "FatJet_jetId > 0")
        .Define("ak8jet_pt", "FatJet_pt[_goodAK8Jets]")
        .Define("ak8jet_eta", "FatJet_eta[_goodAK8Jets]")
        .Define("ak8jet_phi", "FatJet_phi[_goodAK8Jets]")
        .Define("ak8jet_mass", "FatJet_mass[_goodAK8Jets]")
        .Define("ak8jet_msoftdrop", "FatJet_msoftdrop[_goodAK8Jets]")
        .Define("ak8jet_nConstituents", "FatJet_nConstituents[_goodAK8Jets]")
        .Define("ak8jet_hbbvsqcd", "FatJet_particleNetWithMass_HbbvsQCD[_goodAK8Jets]")
        .Define("ak8jet_wvsqcd", "FatJet_particleNetWithMass_WvsQCD[_goodAK8Jets]")
        .Define("ak8jet_zvsqcd", "FatJet_particleNetWithMass_ZvsQCD[_goodAK8Jets]")
        .Define("ak8jet_ht", "Sum(ak8jet_pt)")
        .Define("ak8jet_n", "Sum(_goodAK8Jets)");
}

RNode AK4JetsSelection(RNode df_) {
    return df_.Define("_goodAK4Jets", "Jet_pt > 30 && "
            "abs(Jet_eta) < 4.7 && "
            "Jet_jetId >= 2")
        .Define("ak4jet_pt", "Jet_pt[_goodAK4Jets]")
        .Define("ak4jet_eta", "Jet_eta[_goodAK4Jets]")
        .Define("ak4jet_phi", "Jet_phi[_goodAK4Jets]")
        .Define("ak4jet_mass", "Jet_mass[_goodAK4Jets]")
        .Define("ak4jet_ht", "Sum(ak4jet_pt)")
        .Define("ak4jet_n", "Sum(_goodAK4Jets)");
}

RNode VBSJetSelections(RNode df_, TMVA::Experimental::RBDT &vbstagger) {
    return df_.Define("_vbsjets_pair_idx", getVBSPairs, {"_goodAK4Jets", "ak4jet_pt"})
        .Define("_comb_1_pt", "ROOT::VecOps::Take(ak4jet_pt, _vbsjets_pair_idx[0])")
        .Define("_comb_1_eta", "ROOT::VecOps::Take(ak4jet_eta, _vbsjets_pair_idx[0])")
        .Define("_comb_1_phi", "ROOT::VecOps::Take(ak4jet_phi, _vbsjets_pair_idx[0])")
        .Define("_comb_1_mass", "ROOT::VecOps::Take(ak4jet_mass, _vbsjets_pair_idx[0])")
        .Define("_comb_2_pt", "ROOT::VecOps::Take(ak4jet_pt, _vbsjets_pair_idx[1])")
        .Define("_comb_2_eta", "ROOT::VecOps::Take(ak4jet_eta, _vbsjets_pair_idx[1])")
        .Define("_comb_2_phi", "ROOT::VecOps::Take(ak4jet_phi, _vbsjets_pair_idx[1])")
        .Define("_comb_2_mass", "ROOT::VecOps::Take(ak4jet_mass, _vbsjets_pair_idx[1])")
        .Define("_pt_jj", "_comb_1_pt + _comb_2_pt")
        .Define("_deta_jj", "abs(_comb_1_eta - _comb_2_eta)")
        .Define("_dphi_jj", "ROOT::VecOps::DeltaPhi(_comb_1_phi, _comb_2_phi)")
        .Define("_m_jj", "ROOT::VecOps::InvariantMasses(_comb_1_pt, _comb_1_eta, _comb_1_phi, _comb_1_mass, _comb_2_pt, _comb_2_eta, _comb_2_phi, _comb_2_mass)")
        .Define("vbs_tagger_score", Compute<12, float>(vbstagger), {"_comb_1_pt", "_comb_2_pt", "_comb_1_eta", "_comb_2_eta", "_comb_1_phi", "_comb_2_phi", "_comb_1_mass", "_comb_2_mass", "_pt_jj", "_deta_jj", "_dphi_jj", "_m_jj"})
        .Define("_vbs_tag_idx", "ROOT::VecOps::ArgMax(vbs_tagger_score)")
        .Define("vbs1_idx", "_vbsjets_pair_idx[0][_vbs_tag_idx]")
        .Define("vbs2_idx", "_vbsjets_pair_idx[1][_vbs_tag_idx]")
        .Define("vbs_score", "vbs_tagger_score[_vbs_tag_idx]")
        .Define("vbs1_pt", "ak4jet_pt[vbs1_idx]")
        .Define("vbs2_pt", "ak4jet_pt[vbs2_idx]")
        .Define("vbs1_eta", "ak4jet_eta[vbs1_idx]")
        .Define("vbs2_eta", "ak4jet_eta[vbs2_idx]")
        .Define("vbs1_phi", "ak4jet_phi[vbs1_idx]")
        .Define("vbs2_phi", "ak4jet_phi[vbs2_idx]")
        .Define("vbs1_mass", "ak4jet_mass[vbs1_idx]")
        .Define("vbs2_mass", "ak4jet_mass[vbs2_idx]");
}

