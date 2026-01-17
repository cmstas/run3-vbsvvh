#include "selections.h"

RNode TriggerSelections(RNode df_, std::string channel, const std::unordered_map<std::string, std::string> &trigger_map)
{
    if (trigger_map.empty())
    {
        std::cerr << "    Warning: No trigger map provided. Skipping trigger selection." << std::endl;
        return df_;
    }
    if (trigger_map.find(channel) == trigger_map.end())
    {
        std::cerr << "    Warning: Channel '" << channel << "' not found in trigger map. Skipping trigger selection." << std::endl;
        return df_;
    }

    std::string trigger_condition = trigger_map.at(channel);
    return df_.Filter(trigger_condition, "C1: Trigger Selection");
}

RNode ElectronSelections(RNode df_)
{
    auto df = df_.Define("Electron_SC_eta", "Electron_eta + Electron_deltaEtaSC")
                  .Define("_looseElectrons", "Electron_pt > 7 &&"
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
    return applyObjectMaskNewAffix(df, "_tightElectrons", "Electron", "electron");
}

RNode MuonSelections(RNode df_)
{
    auto df = df_.Define("_looseMuons", "Muon_pt > 5 && "
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
                  .Define("nMuon_Tight", "nMuon_Loose == 0 ? 0 : Sum(_tightMuons)")
                  .Define("vvhTightLepMaskMuon", "_tightMuons");
    return applyObjectMaskNewAffix(df, "_tightMuons", "Muon", "muon");
}

RNode LeptonSelections(RNode df_)
{
    auto df = ElectronSelections(df_);
    df = MuonSelections(df);
    return df.Define("lepton_pt", "Concatenate(electron_pt, muon_pt)")
            .Define("_leptonSorted", "Argsort(-lepton_pt)")
            .Redefine("lepton_pt", "Take(lepton_pt, _leptonSorted)")
            .Define("lepton_eta", "Take(Concatenate(electron_eta, muon_eta), _leptonSorted)")
            .Define("lepton_phi", "Take(Concatenate(electron_phi, muon_phi), _leptonSorted)")
            .Define("lepton_mass", "Take(Concatenate(electron_mass, muon_mass), _leptonSorted)")
            .Define("lepton_charge", "Take(Concatenate(electron_charge, muon_charge), _leptonSorted)");
}

RNode AK4JetsSelection(RNode df_)
{
    auto df = df_.Define("_dR_ak4_lep", VVdR, {"Jet_eta", "Jet_phi", "lepton_eta", "lepton_phi"})
                .Define("_good_ak4jets", " _dR_ak4_lep > 0.4 && "
                                        "Jet_pt > 20 && "
                                        "((abs(Jet_eta) <= 2.5) || (abs(Jet_eta) >= 3.0 && abs(Jet_eta) < 5.0) || "
                                        "((Jet_pt > 30 && (Jet_neHEF < ((Jet_pt - 30) / 20) && ((Jet_puIdDisc > 0.3) || " 
                                        "(Jet_puIdDisc > 2 * Jet_chHEF)) && (Jet_puIdDisc > (-1 / 35 * (Jet_pt - 30) + 0.4)) && "
                                        "(Jet_puIdDisc > -0.6) && "
                                        "(Jet_muEF < 0.4))) || "
                                        "(Jet_pt < 30 && Jet_pt > 20 && (Jet_neHEF < ((Jet_pt - 20) / 15) && ((Jet_puIdDisc > 0.3) || " 
                                        "(Jet_puIdDisc > 2 * Jet_chHEF)) && (Jet_puIdDisc > (-1 / 25 * (Jet_pt - 20) + 0.4)) && "
                                        "(Jet_puIdDisc > -0.6) && "
                                        "(Jet_muEF < 0.4))))) && "
                                        "((is2016 && Jet_jetId >= 1) || (!is2016 && Jet_jetId >= 2))");
    
    df = applyJetVetoMaps(df);
    df = df.Redefine("_good_ak4jets", "_good_ak4jets && !Jet_vetoMap");

    df = df.Define("Jet_isTightBTag", "Jet_btagUParTAK4B > 0.4648")
            .Define("Jet_isMediumBTag", "Jet_btagUParTAK4B > 0.1272")
            .Define("Jet_isLooseBTag", "Jet_btagUParTAK4B > 0.0246");

    df = applyObjectMaskNewAffix(df, "_good_ak4jets", "Jet", "jet");
    df = df.Define("ht_jets", "Sum(jet_pt)");
    return df;
}

RNode AK8JetsSelection(RNode df_)
{
    auto df = df_.Define("_dR_ak8_lep", VVdR, {"FatJet_eta", "FatJet_phi", "lepton_eta", "lepton_phi"})
                  .Define("_good_ak8jets", "_dR_ak8_lep > 0.8 && "
                                           "FatJet_pt > 250 && "
                                           "abs(FatJet_eta) <= 2.5 && "
                                           "FatJet_msoftdrop > 40 && "
                                           "FatJet_jetId > 0");

    df = applyObjectMaskNewAffix(df, "_good_ak8jets", "FatJet", "fatjet");
    df = df.Define("ht_fatjets", "Sum(fatjet_pt)");

    return df;
}

RNode runPreselection(RNode df_, std::string channel, bool noCut)
{
    auto df = LeptonSelections(df_);
    df = AK4JetsSelection(df);
    df = AK8JetsSelection(df);
    df = df.Define("jet_minDrFromAnyGoodFatJet", dRfromClosestJet, {"jet_eta", "jet_phi", "fatjet_eta", "fatjet_phi"})
            .Define("jet_passFatJetOverlapRemoval", "jet_minDrFromAnyGoodFatJet>0.8");
    
    if (noCut)
        return df; // for spanet training data

    df = TriggerSelections(df, channel, TriggerMap);
    if (channel == "1Lep2FJ")
    {
        df = df.Filter("((nMuon_Loose == 1 && nMuon_Tight == 1 && nElectron_Loose == 0 && nElectron_Tight == 0) || "
                       "(nMuon_Loose == 0 && nMuon_Tight == 0 && nElectron_Loose == 1 && nElectron_Tight == 1)) && "
                       "(lepton_pt[0] > 40)",
                       "C2: 1-lepton selection");
    }
    else if (channel == "0Lep3FJ")
    {
        df = df.Filter("nMuon_Loose == 0 && nElectron_Loose == 0", "C2: 0-lepton selection");
    }
    return df;
}
