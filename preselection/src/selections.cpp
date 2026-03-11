#include "selections.h"
#include "cutflow.h"

RNode TriggerSelections(RNode df_, std::string channel, const std::unordered_map<std::string, std::string> &trigger_map)
{
    if (trigger_map.find(channel) == trigger_map.end())
    {
        std::cerr << "    Warning: Channel '" << channel << "' not found in trigger map. Skipping trigger selection." << std::endl;
        return df_;
    }
    std::string trigger_condition = trigger_map.at(channel);
    std::regex hlt_regex("HLT_[a-zA-Z0-9_]+");
    auto hlt_begin = std::sregex_iterator(trigger_condition.begin(), trigger_condition.end(), hlt_regex);
    auto hlt_end = std::sregex_iterator();
    for (auto it = hlt_begin; it != hlt_end; ++it)
    {
        std::string hlt_name = it->str(0);
        auto column_names = df_.GetColumnNames();
        if (std::find(column_names.begin(), column_names.end(), hlt_name) == column_names.end()){
            df_ = df_.DefaultValueFor(hlt_name, (int)0);
        }
    }
    
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
                                             "Electron_pt > 10 && "
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
                                         "Muon_pt > 10 && "
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
                .Define("_good_ak4jets", "_dR_ak4_lep > 0.4 &&"
                                        "((isRun3 && ((Jet_pt > 20 && (abs(Jet_eta) <= 2.5 || abs(Jet_eta) >= 3.0) && abs(Jet_eta) < 5.0) || "
                                        "(Jet_pt > 50 && abs(Jet_eta) > 2.5 && abs(Jet_eta) < 3.0))) ||" // horn removal JME recommendation for Run3
                                        "(isRun2 && Jet_pt > 20 && abs(Jet_eta) < 5.0)) &&  "
                                        "((is2016 && Jet_jetId >= 1) || (!is2016 && Jet_jetId >= 2))");
    
    df = applyJetVetoMaps(df);
    df = df.Filter("(isRun2) || (isRun3 && !Any(Jet_vetoMap))"); // for Run3, events with any jet in veto region are removed
    df = df.Redefine("_good_ak4jets", "_good_ak4jets && !Jet_vetoMap");

    df = df.Define("Jet_isTightBTag", isbTagTight, {"year", "Jet_btagUParTAK4B"})
            .Define("Jet_isMediumBTag", isbTagMedium, {"year", "Jet_btagUParTAK4B"})
            .Define("Jet_isLooseBTag", isbTagLoose, {"year", "Jet_btagUParTAK4B"});

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
                                           "FatJet_jetId > 0")
                  .Define("nFatJets", "Sum(_good_ak8jets)");

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
    
    Cutflow::Add(df_, "All events");
    df = TriggerSelections(df, channel, TriggerMap);
    Cutflow::Add(df, "C1: Trigger selection");

    if (channel == "all_events"){
        df = df.Filter(
            "nMuon_Loose > -1", // Probably there is a better way to write a pass through
            "C2: all_events"
        );
    }
    else if (channel == "0lep_0FJ"){
        df = df.Filter(
            "((nMuon_Loose == 0) && (nElectron_Loose == 0)) &&"
            "(nFatJets == 0)",
            "C2: 0lep_0FJ"
        );
    }
    else if (channel == "0lep_1FJ"){
        df = df.Filter(
            "((nMuon_Loose == 0) && (nElectron_Loose == 0)) &&"
            "(nFatJets == 1)",
            "C2: 0lep_1FJ"
        );
    }
    else if (channel == "0lep_2FJ"){
        df = df.Filter(
            "((nMuon_Loose == 0) && (nElectron_Loose == 0)) &&"
            "(nFatJets == 2)",
            "C2: 0lep_2FJ"
        );
    }
    else if (channel == "0lep_3FJ"){
        df = df.Filter(
            "((nMuon_Loose == 0) && (nElectron_Loose == 0)) &&"
            "(nFatJets == 3)",
            "C2: 0lep_3FJ"
        );
    }
    else if (channel == "1lep_1FJ"){
        df = df.Filter("((nMuon_Loose == 1 && nMuon_Tight == 1 && nElectron_Loose == 0 && nElectron_Tight == 0) || "
                       "(nMuon_Loose == 0 && nMuon_Tight == 0 && nElectron_Loose == 1 && nElectron_Tight == 1)) && "
                       "(lepton_pt[0] > 40)");
        Cutflow::Add(df, "C2: 1-lepton selection");
        df = df.Filter("nfatjet == 1");
        Cutflow::Add(df, "C3: exactly 1 fat jet");
        df = df.Filter("njet >= 4");
        Cutflow::Add(df, "C4: at-least 4 jets");
    }
    else if (channel == "1lep_2FJ"){
    	df = df.Filter("((nMuon_Loose == 1 && nMuon_Tight == 1 && nElectron_Loose == 0 && nElectron_Tight == 0) || "
                       "(nMuon_Loose == 0 && nMuon_Tight == 0 && nElectron_Loose == 1 && nElectron_Tight == 1)) && "
                       "(lepton_pt[0] > 40)");
        Cutflow::Add(df, "C2: 1-lepton selection");
        df = df.Filter("nfatjet >= 2");
        Cutflow::Add(df, "C3: at-least 2 fat jets");
        df = df.Filter("njet >= 2");
        Cutflow::Add(df, "C4: at-least 2 jets");
    }
    else if (channel == "2lep_1FJ"){
        df = df.Filter(
            "((nMuon_Loose + nElectron_Loose) == 2) &&"
            "(nFatJets == 1)",
            "C2: 2lep_1FJ"
        );
    }
    else if (channel == "2lep_2FJ"){
        df = df.Filter(
            "((nMuon_Loose + nElectron_Loose) == 2) &&"
            "(nFatJets == 2)",
            "C2: 2lep_2FJ"
        );
    }
    //else if (channel == "2lepSS"){
    //    df = df.Filter(
    //        "((nMuon_Loose + nElectron_Loose) == 2) &&"
    //        //same sign requirement
    //        "C2: 2lepSS"
    //    );
    //}
    else if (channel == "3lep"){
        df = df.Filter(
            "((nMuon_Loose + nElectron_Loose) == 3)",
            "C2: 3lep"
        );
    }
    else if (channel == "4lep"){
        df = df.Filter(
            "((nMuon_Loose + nElectron_Loose) == 4)",
            "C2: 4lep"
        );
    }
    return df;
}
