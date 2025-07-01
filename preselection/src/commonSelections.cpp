#include "commonSelections.h"

namespace CommonSelections {
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

    RNode ElectronSelections(RNode df_) {
        return df_.Define("Electron_SC_eta", "Electron_eta + Electron_deltaEtaSC")
            .Define("_looseElectrons", 
                "Electron_pt > 7 &&"
                "abs(Electron_SC_eta) < 2.5 && "
                "((abs(Electron_SC_eta) <= 1.479 && abs(Electron_dxy) <= 0.05 && abs(Electron_dz) < 0.1) || (abs(Electron_dxy) <= 0.1 && abs(Electron_dz) < 0.2)) && "
                "abs(Electron_sip3d) < 8 && "
                "Electron_cutBased >= 2 && "
                "Electron_pfRelIso03_all < 0.4 && "
                "Electron_lostHits <= 1")
            .Define("nElectron_Loose", "nElectron == 0 ? 0 : Sum(_looseElectrons)")
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
            .Define("nElectron_Tight", "nElectron_Loose == 0 ? 0 : Sum(_tightElectrons)")
            .Redefine("Electron_pt", "Electron_pt[_tightElectrons]")
            .Redefine("Electron_eta", "Electron_eta[_tightElectrons]")
            .Redefine("Electron_SC_eta", "Electron_SC_eta[_tightElectrons]")
            .Redefine("Electron_phi", "Electron_phi[_tightElectrons]")
            .Redefine("Electron_mass", "Electron_mass[_tightElectrons]");
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
            .Define("nMuon_Loose", "nMuon == 0 ? 0 : Sum(_looseMuons)")
            .Define("_tightMuons", "_looseMuons && "
                "Muon_pt > 30 && "
                "Muon_pfIsoId > 4 && "
                "Muon_tightCharge == 2 && "
                "Muon_highPurity && "
                "Muon_tightId")
            .Define("nMuon_Tight", "nMuon_Loose == 0 ? 0 : Sum(_tightMuons)")
            .Redefine("Muon_pt", "Muon_pt[_tightMuons]")
            .Redefine("Muon_eta", "Muon_eta[_tightMuons]")
            .Redefine("Muon_phi", "Muon_phi[_tightMuons]")
            .Redefine("Muon_mass", "Muon_mass[_tightMuons]");
    }

    RNode AK4JetsSelection(RNode df_) {
        return df_.Define("_good_ak4jets", "Jet_pt > 30 && "
                "abs(Jet_eta) < 4.7 && "
                "Jet_jetId >= 2")
            .Redefine("nJet", "Sum(_good_ak4jets)")
            .Redefine("Jet_pt", "Jet_pt[_good_ak4jets]")
            .Redefine("Jet_eta", "Jet_eta[_good_ak4jets]")
            .Redefine("Jet_phi", "Jet_phi[_good_ak4jets]")
            .Redefine("Jet_mass", "Jet_mass[_good_ak4jets]")
            .Redefine("Jet_jetId", "Jet_jetId[_good_ak4jets]")
            .Redefine("Jet_btagDeepFlavB", "Jet_btagDeepFlavB[_good_ak4jets]")
            .Define("Jet_ht", "Sum(Jet_pt)");
    }

    RNode AK8JetsSelection(RNode df_) {
        return df_.Define("_good_ak8jets",
                "FatJet_pt > 250 && "
                "abs(FatJet_eta) <= 2.5 && "
                "FatJet_mass > 50 && "
                "FatJet_msoftdrop > 40 && "
                "FatJet_jetId > 0")
            .Redefine("nFatJet", "Sum(_good_ak8jets)")
            .Redefine("FatJet_pt", "FatJet_pt[_good_ak8jets]")
            .Redefine("FatJet_eta", "FatJet_eta[_good_ak8jets]")
            .Redefine("FatJet_phi", "FatJet_phi[_good_ak8jets]")
            .Redefine("FatJet_mass", "FatJet_mass[_good_ak8jets]")
            .Redefine("FatJet_msoftdrop", "FatJet_msoftdrop[_good_ak8jets]")
            .Redefine("FatJet_nConstituents", "FatJet_nConstituents[_good_ak8jets]")
            .Redefine("FatJet_particleNet_XbbVsQCD", "FatJet_particleNet_XbbVsQCD[_good_ak8jets]")
            .Redefine("FatJet_particleNet_XqqVsQCD", "FatJet_particleNet_XqqVsQCD[_good_ak8jets]")
            .Redefine("FatJet_particleNetWithMass_WvsQCD", "FatJet_particleNetWithMass_WvsQCD[_good_ak8jets]")
            .Redefine("FatJet_particleNetWithMass_ZvsQCD", "FatJet_particleNetWithMass_ZvsQCD[_good_ak8jets]")
            .Redefine("FatJet_particleNetWithMass_HbbvsQCD", "FatJet_particleNetWithMass_HbbvsQCD[_good_ak8jets]")
            .Define("FatJet_ht", "Sum(FatJet_pt[_good_ak8jets])");
    }

    RNode VBSJetSelections(RNode df_, TMVA::Experimental::RBDT &vbstagger) {
        return df_.Define("_good_vbsjets", "Jet_pt > 30 && "
                "abs(Jet_eta) < 4.7 && "
                "Jet_jetId >= 2")
            .Define("_vbsjet_pt", "Jet_pt[_good_vbsjets]")
            .Define("_vbsjet_eta", "Jet_eta[_good_vbsjets]")
            .Define("_vbsjet_phi", "Jet_phi[_good_vbsjets]")
            .Define("_vbsjet_mass", "Jet_mass[_good_vbsjets]")
            .Define("_vbsjets_pair_idx", getVBSPairs, {"_good_vbsjets", "_vbsjet_pt"})
            .Define("_comb_1_pt", "ROOT::VecOps::Take(_vbsjet_pt, _vbsjets_pair_idx[0])")
            .Define("_comb_1_eta", "ROOT::VecOps::Take(_vbsjet_eta, _vbsjets_pair_idx[0])")
            .Define("_comb_1_phi", "ROOT::VecOps::Take(_vbsjet_phi, _vbsjets_pair_idx[0])")
            .Define("_comb_1_mass", "ROOT::VecOps::Take(_vbsjet_mass, _vbsjets_pair_idx[0])")
            .Define("_comb_2_pt", "ROOT::VecOps::Take(_vbsjet_pt, _vbsjets_pair_idx[1])")
            .Define("_comb_2_eta", "ROOT::VecOps::Take(_vbsjet_eta, _vbsjets_pair_idx[1])")
            .Define("_comb_2_phi", "ROOT::VecOps::Take(_vbsjet_phi, _vbsjets_pair_idx[1])")
            .Define("_comb_2_mass", "ROOT::VecOps::Take(_vbsjet_mass, _vbsjets_pair_idx[1])")
            .Define("_pt_jj", "_comb_1_pt + _comb_2_pt")
            .Define("_deta_jj", "abs(_comb_1_eta - _comb_2_eta)")
            .Define("_dphi_jj", "ROOT::VecOps::DeltaPhi(_comb_1_phi, _comb_2_phi)")
            .Define("_m_jj", "ROOT::VecOps::InvariantMasses(_comb_1_pt, _comb_1_eta, _comb_1_phi, _comb_1_mass, _comb_2_pt, _comb_2_eta, _comb_2_phi, _comb_2_mass)")
            .Define("_vbs_tagger_score", Compute<12, float>(vbstagger), {"_comb_1_pt", "_comb_2_pt", "_comb_1_eta", "_comb_2_eta", "_comb_1_phi", "_comb_2_phi", "_comb_1_mass", "_comb_2_mass", "_pt_jj", "_deta_jj", "_dphi_jj", "_m_jj"})
            .Define("_vbs_tag_idx", "ROOT::VecOps::ArgMax(_vbs_tagger_score)")
            .Define("bdt_vbs1_idx", "_vbsjets_pair_idx[0][_vbs_tag_idx]")
            .Define("bdt_vbs2_idx", "_vbsjets_pair_idx[1][_vbs_tag_idx]")
            .Define("vbs_score", "_vbs_tagger_score[_vbs_tag_idx]")
            .Define("vbs1_pt", "_comb_1_pt[_vbs_tag_idx]")
            .Define("vbs1_eta", "_comb_1_eta[_vbs_tag_idx]")
            .Define("vbs1_phi", "_comb_1_phi[_vbs_tag_idx]")
            .Define("vbs1_mass", "_comb_1_mass[_vbs_tag_idx]")
            .Define("vbs2_pt", "_comb_2_pt[_vbs_tag_idx]")
            .Define("vbs2_eta", "_comb_2_eta[_vbs_tag_idx]")
            .Define("vbs2_phi", "_comb_2_phi[_vbs_tag_idx]")
            .Define("vbs2_mass", "_comb_2_mass[_vbs_tag_idx]")
            .Define("vbs_ptjj", "_pt_jj[_vbs_tag_idx]")
            .Define("vbs_detajj", "_deta_jj[_vbs_tag_idx]")
            .Define("vbs_dphijj", "_dphi_jj[_vbs_tag_idx]")
            .Define("vbs_mjj", "_m_jj[_vbs_tag_idx]");
    }
}