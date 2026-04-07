#include "selections.h"
#include "cutflow.h"

RNode TriggerSelections(RNode df_, std::string trigger_logic_string) {

    // Add default values to the df for all triggers that show up in the trigger_logic_string
    std::string trigger_condition = trigger_logic_string;
    std::regex hlt_regex("HLT_[a-zA-Z0-9_]+");
    auto hlt_begin = std::sregex_iterator(trigger_condition.begin(), trigger_condition.end(), hlt_regex);
    auto hlt_end = std::sregex_iterator();
    std::vector<std::string> seen_triggers;

    for (auto it = hlt_begin; it != hlt_end; ++it)
    {
        std::string hlt_name = it->str(0);
        if (std::find(seen_triggers.begin(), seen_triggers.end(), hlt_name) != seen_triggers.end())
            continue;
        seen_triggers.push_back(hlt_name);
        df_ = df_.DefaultValueFor(hlt_name, (bool)false);
    }

    // Filter based on this the trigger logic
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

    if (noCut) return df;

    // Passthrough
    if (channel == "all_events"){
        df = df.Filter(
            "nMuon_Loose > -1", // Probably there is a better way to write a pass through
            "C2: all_events"
        );
    }

    // 0lep_0FJ
    else if (channel == "0lep_0FJ"){

        df = TriggerSelections(df,trigger_logic_string_ht);
        Cutflow::Add(df, "C1: Trigger selection");

        df = df.Filter(
            "((nMuon_Loose == 0) && (nElectron_Loose == 0)) &&"
            "(nFatJets == 0)",
            "C2: 0lep_0FJ"
        );
        df = df.Define("met_significance", "PuppiMET_significance")
                .Define("met_uncorrPt", "PuppiMET_pt")
                .Define("met_uncorrPhi", "PuppiMET_phi");
    }

    // 0lep_1FJ
    else if (channel == "0lep_1FJ"){

        df = TriggerSelections(df,trigger_logic_string_ht);
        Cutflow::Add(df, "C1: Trigger selection");

        df = df.Filter(
            "((nMuon_Loose == 0) && (nElectron_Loose == 0)) &&"
            "(nFatJets == 1)",
            "C2: 0lep_1FJ"
        );
        df = df.Define("met_significance", "PuppiMET_significance")
                .Define("met_uncorrPt", "PuppiMET_pt")
                .Define("met_uncorrPhi", "PuppiMET_phi");

    }

    // 0lep_1FJ_met
    else if (channel == "0lep_1FJ_met"){

        df = TriggerSelections(df,trigger_logic_string_met);
        Cutflow::Add(df, "C1: Trigger selection");

        df = df.Filter(
            "((nMuon_Loose == 0) && (nElectron_Loose == 0)) &&"
            "(nFatJets == 1)",
            "C2: 0lep_1FJ"
        );
    }

    // 0lep_2FJ
    else if (channel == "0lep_2FJ"){

        df = TriggerSelections(df,trigger_logic_string_ht);
        Cutflow::Add(df, "C1: Trigger selection");

        df = df.Filter(
            "((nMuon_Loose == 0) && (nElectron_Loose == 0)) &&"
            "(nFatJets == 2)",
            "C2: 0lep_2FJ"
        );
    }

    // 0lep_2FJ_met
    else if (channel == "0lep_2FJ_met"){

        df = TriggerSelections(df,trigger_logic_string_met);
        Cutflow::Add(df, "C1: Trigger selection");

        df = df.Filter(
            "((nMuon_Loose == 0) && (nElectron_Loose == 0)) &&"
            "(nFatJets == 2)",
            "C2: 0lep_2FJ"
        );
    }

    // 0lep_3FJ
    else if (channel == "0lep_3FJ"){

        df = TriggerSelections(df,trigger_logic_string_ht);
        Cutflow::Add(df, "C1: Trigger selection");

        df = df.Filter(
            "((nMuon_Loose == 0) && (nElectron_Loose == 0)) &&"
            "(nFatJets == 3)",
            "C2: 0lep_3FJ"
        );
    }

    // 1lep_1FJ
    else if (channel == "1lep_1FJ"){

        df = TriggerSelections(df,trigger_logic_string_singlelep);
        Cutflow::Add(df, "C1: Trigger selection");

        df = df.Filter("((nMuon_Loose == 1 && nMuon_Tight == 1 && nElectron_Loose == 0 && nElectron_Tight == 0) || "
                       "(nMuon_Loose == 0 && nMuon_Tight == 0 && nElectron_Loose == 1 && nElectron_Tight == 1)) && "
                       "(lepton_pt[0] > 40)");
        Cutflow::Add(df, "C2: 1-lepton selection");

        df = df.Filter("nfatjet == 1");
        Cutflow::Add(df, "C3: exactly 1 fat jet");
        
        df = df.Filter("njet >= 4");
        Cutflow::Add(df, "C4: at-least 4 jets");

        df = df.Define("_vbs_candidate_jet_pairs", VBSBDTInfer, {"jet_pt", "jet_eta", "jet_phi", "jet_mass", "isRun2"})
            .Define("vbs_jet1_idx", "_vbs_candidate_jet_pairs[0]")
            .Define("vbs_jet2_idx", "_vbs_candidate_jet_pairs[1]")
            .Define("vbs_jet1_pt", "vbs_jet1_idx != -1 ? jet_pt[vbs_jet1_idx] : -999.0f")
            .Define("vbs_jet1_eta", "vbs_jet1_idx != -1 ? jet_eta[vbs_jet1_idx] : -999.0f")
            .Define("vbs_jet1_phi", "vbs_jet1_idx != -1 ? jet_phi[vbs_jet1_idx] : -999.0f")
            .Define("vbs_jet1_mass", "vbs_jet1_idx != -1 ? jet_mass[vbs_jet1_idx] : -999.0f")
            .Define("vbs_jet2_pt", "vbs_jet2_idx != -1 ? jet_pt[vbs_jet2_idx] : -999.0f")
            .Define("vbs_jet2_eta", "vbs_jet2_idx != -1 ? jet_eta[vbs_jet2_idx] : -999.0f")
            .Define("vbs_jet2_phi", "vbs_jet2_idx != -1 ? jet_phi[vbs_jet2_idx] : -999.0f")
            .Define("vbs_jet2_mass", "vbs_jet2_idx != -1 ? jet_mass[vbs_jet2_idx] : -999.0f")
            .Define("vbs_mjj", "(vbs_jet1_idx != -1 && vbs_jet2_idx != -1) ? (ROOT::Math::PtEtaPhiMVector(vbs_jet1_pt, vbs_jet1_eta, vbs_jet1_phi, vbs_jet1_mass) + "
                "ROOT::Math::PtEtaPhiMVector(vbs_jet2_pt, vbs_jet2_eta, vbs_jet2_phi, vbs_jet2_mass)).M() : -999.0f")
            .Define("vbs_detajj", "(vbs_jet1_idx != -1 && vbs_jet2_idx != -1) ? abs(vbs_jet1_eta - vbs_jet2_eta) : -999.0f")
            .Define("vbs_candidate_found", "vbs_jet1_idx != -1 && vbs_jet2_idx != -1")
            .Filter("vbs_candidate_found");

        Cutflow::Add(df, "C5: VBS candidate selection");
        
        df = df.Define("_fatjet_vbs1_dR", VdR, {"fatjet_eta", "fatjet_phi", "vbs_jet1_eta", "vbs_jet1_phi"})
            .Define("_fatjet_vbs2_dR", VdR, {"fatjet_eta", "fatjet_phi", "vbs_jet2_eta", "vbs_jet2_phi"})
            .Define("_boosted_candidate_jets", 
                "_fatjet_vbs1_dR >= 0.8 && "
                "_fatjet_vbs2_dR >= 0.8")
            .Define("fatjet_HvsQCD", "fatjet_globalParT3_Xbb / (fatjet_globalParT3_Xbb + fatjet_globalParT3_QCD)")
            .Define("fatjet_VvsQCD", "fatjet_globalParT3_Xqq / (fatjet_globalParT3_Xqq + fatjet_globalParT3_Xcs + fatjet_globalParT3_QCD)")
            .Define("fatjet_is_h", "Sum(_boosted_candidate_jets) > 0 && (fatjet_HvsQCD[_boosted_candidate_jets][0] > fatjet_VvsQCD[_boosted_candidate_jets][0])")
            .Define("fatjet_is_v", "!fatjet_is_h")
            .Define("boosted_h_candidate_eta", "fatjet_is_h ? fatjet_eta[0] : -999.0f")
            .Define("boosted_h_candidate_phi", "fatjet_is_h ? fatjet_phi[0] : -999.0f")
            .Define("boosted_h_candidate_mass", "fatjet_is_h ? fatjet_mass[0] : -999.0f")
            .Define("boosted_h_candidate_pt", "fatjet_is_h ? fatjet_pt[0] : -999.0f")
            .Define("boosted_h_candidate_tau21", "fatjet_is_h ? fatjet_tau2[0] / fatjet_tau1[0] : -999.0f")
            .Define("boosted_h_candidate_score", "fatjet_is_h ? fatjet_HvsQCD[0] : -999.0f")
            .Define("boosted_v_candidate_eta", "fatjet_is_v ? fatjet_eta[0] : -999.0f")
            .Define("boosted_v_candidate_phi", "fatjet_is_v ? fatjet_phi[0] : -999.0f")
            .Define("boosted_v_candidate_mass", "fatjet_is_v ? fatjet_mass[0] : -999.0f")
            .Define("boosted_v_candidate_pt", "fatjet_is_v ? fatjet_pt[0] : -999.0f")
            .Define("boosted_v_candidate_tau21", "fatjet_is_v ? fatjet_tau2[0] / fatjet_tau1[0] : -999.0f")
            .Define("boosted_v_candidate_score", "fatjet_is_v ? fatjet_VvsQCD[0] : -999.0f")
            .Filter("Sum(_boosted_candidate_jets) > 0");

        Cutflow::Add(df, "C6: Boosted candidate selection");

        df = df.Define("jet_v_dR", VdR, {"jet_eta", "jet_phi", "boosted_v_candidate_eta", "boosted_v_candidate_phi"})
            .Define("jet_h_dR", VdR, {"jet_eta", "jet_phi", "boosted_h_candidate_eta", "boosted_h_candidate_phi"})
            .Define("jet_vbs1_dR", VdR, {"jet_eta", "jet_phi", "vbs_jet1_eta", "vbs_jet1_phi"})
            .Define("jet_vbs2_dR", VdR, {"jet_eta", "jet_phi", "vbs_jet2_eta", "vbs_jet2_phi"})
            .Define("_resolved_candidate_jets",
                "abs(jet_eta) <= 2.5 && "
                "jet_v_dR >= 0.8 && "
                "jet_h_dR >= 0.8 && "
                "jet_vbs1_dR >= 0.4 && "
                "jet_vbs2_dR >= 0.4 ")
            .Define("_resolved_candidate_pt", "jet_pt[_resolved_candidate_jets]")
            .Define("_resolved_candidate_eta", "jet_eta[_resolved_candidate_jets]")
            .Define("_resolved_candidate_phi", "jet_phi[_resolved_candidate_jets]")
            .Define("_resolved_candidate_mass", "jet_mass[_resolved_candidate_jets]")
            .Define("_resolved_candidate_btag", "jet_btagUParTAK4B[_resolved_candidate_jets]")
            .Define("_resolved_candidate_pairs", getJetPairs, {"_resolved_candidate_pt"})
            .Define("_resolved_candidate_pairs1_pt", "Take(_resolved_candidate_pt, _resolved_candidate_pairs[0], -999.0f)")
            .Define("_resolved_candidate_pairs2_pt", "Take(_resolved_candidate_pt, _resolved_candidate_pairs[1], -999.0f)")
            .Define("_resolved_candidate_pairs1_eta", "Take(_resolved_candidate_eta, _resolved_candidate_pairs[0], -999.0f)")
            .Define("_resolved_candidate_pairs2_eta", "Take(_resolved_candidate_eta, _resolved_candidate_pairs[1], -999.0f)")
            .Define("_resolved_candidate_pairs1_phi", "Take(_resolved_candidate_phi, _resolved_candidate_pairs[0], -999.0f)")
            .Define("_resolved_candidate_pairs2_phi", "Take(_resolved_candidate_phi, _resolved_candidate_pairs[1], -999.0f)")
            .Define("_resolved_candidate_pairs1_mass", "Take(_resolved_candidate_mass, _resolved_candidate_pairs[0], -999.0f)")
            .Define("_resolved_candidate_pairs2_mass", "Take(_resolved_candidate_mass, _resolved_candidate_pairs[1], -999.0f)")
            .Define("_resolved_candidate_pairs1_btag", "Take(_resolved_candidate_btag, _resolved_candidate_pairs[0], -999.0f)")
            .Define("_resolved_candidate_pairs2_btag", "Take(_resolved_candidate_btag, _resolved_candidate_pairs[1], -999.0f)")
            .Define("_resolved_candidate_ptjj", VVInvariantPt, {"_resolved_candidate_pairs1_pt", "_resolved_candidate_pairs1_eta", "_resolved_candidate_pairs1_phi", "_resolved_candidate_pairs1_mass", "_resolved_candidate_pairs2_pt", "_resolved_candidate_pairs2_eta", "_resolved_candidate_pairs2_phi", "_resolved_candidate_pairs2_mass"})
            .Define("resolved_candidate_mjj", VVInvariantMass, {"_resolved_candidate_pairs1_pt", "_resolved_candidate_pairs1_eta", "_resolved_candidate_pairs1_phi", "_resolved_candidate_pairs1_mass", "_resolved_candidate_pairs2_pt", "_resolved_candidate_pairs2_eta", "_resolved_candidate_pairs2_phi", "_resolved_candidate_pairs2_mass"})
            .Define("resolved_candidate_dR", VVDeltaR, {"_resolved_candidate_pairs1_eta", "_resolved_candidate_pairs1_phi", "_resolved_candidate_pairs2_eta", "_resolved_candidate_pairs2_phi"})
            .Define("_sorted_resolved_dR", "Argsort(resolved_candidate_dR)")
            .Define("resolved_mjj_1", "resolved_candidate_mjj.size() > 0 ? resolved_candidate_mjj[_sorted_resolved_dR[0]] : -999.0f")
            .Define("resolved_mjj_2", "resolved_candidate_mjj.size() > 1 ? resolved_candidate_mjj[_sorted_resolved_dR[1]] : -999.0f")
            .Define("resolved_mjj_3", "resolved_candidate_mjj.size() > 2 ? resolved_candidate_mjj[_sorted_resolved_dR[2]] : -999.0f")
            .Define("resolved_mjj_4", "resolved_candidate_mjj.size() > 3 ? resolved_candidate_mjj[_sorted_resolved_dR[3]] : -999.0f")
            .Define("resolved_mjj_5", "resolved_candidate_mjj.size() > 4 ? resolved_candidate_mjj[_sorted_resolved_dR[4]] : -999.0f")
            .Define("resolved_dR_1", "resolved_candidate_dR.size() > 0 ? resolved_candidate_dR[_sorted_resolved_dR[0]] : -999.0f")
            .Define("resolved_dR_2", "resolved_candidate_dR.size() > 1 ? resolved_candidate_dR[_sorted_resolved_dR[1]] : -999.0f")
            .Define("resolved_dR_3", "resolved_candidate_dR.size() > 2 ? resolved_candidate_dR[_sorted_resolved_dR[2]] : -999.0f")
            .Define("resolved_dR_4", "resolved_candidate_dR.size() > 3 ? resolved_candidate_dR[_sorted_resolved_dR[3]] : -999.0f")
            .Define("resolved_dR_5", "resolved_candidate_dR.size() > 4 ? resolved_candidate_dR[_sorted_resolved_dR[4]] : -999.0f")
            .Define("resolved_ptjj_1", "_resolved_candidate_ptjj.size() > 0 ? _resolved_candidate_ptjj[_sorted_resolved_dR[0]] : -999.0f")
            .Define("resolved_ptjj_2", "_resolved_candidate_ptjj.size() > 1 ? _resolved_candidate_ptjj[_sorted_resolved_dR[1]] : -999.0f")
            .Define("resolved_ptjj_3", "_resolved_candidate_ptjj.size() > 2 ? _resolved_candidate_ptjj[_sorted_resolved_dR[2]] : -999.0f")
            .Define("resolved_ptjj_4", "_resolved_candidate_ptjj.size() > 3 ? _resolved_candidate_ptjj[_sorted_resolved_dR[3]] : -999.0f")
            .Define("resolved_ptjj_5", "_resolved_candidate_ptjj.size() > 4 ? _resolved_candidate_ptjj[_sorted_resolved_dR[4]] : -999.0f")
            .Filter("Sum(_resolved_candidate_jets) >= 2");

            Cutflow::Add(df, "C7: Resolved candidate selection");
    }

    else if (channel == "1lep_2FJ"){

        df = TriggerSelections(df,trigger_logic_string_singlelep);
        Cutflow::Add(df, "C1: Trigger selection");

    	df = df.Filter("((nMuon_Loose == 1 && nMuon_Tight == 1 && nElectron_Loose == 0 && nElectron_Tight == 0) || "
                       "(nMuon_Loose == 0 && nMuon_Tight == 0 && nElectron_Loose == 1 && nElectron_Tight == 1)) && "
                       "(lepton_pt[0] > 40)");
        Cutflow::Add(df, "C2: 1-lepton selection");

        df = df.Filter("nfatjet >= 2");
        Cutflow::Add(df, "C3: at-least 2 fat jets");

        df = df.Filter("njet >= 2");
        Cutflow::Add(df, "C4: at-least 2 jets");

        df = df.Define("_vbs_candidate_jet_pairs", VBSBDTInfer, {"jet_pt", "jet_eta", "jet_phi", "jet_mass", "isRun2"})
            .Define("vbs_jet1_idx", "_vbs_candidate_jet_pairs[0]")
            .Define("vbs_jet2_idx", "_vbs_candidate_jet_pairs[1]")
            .Define("vbs_jet1_pt", "vbs_jet1_idx != -1 ? jet_pt[vbs_jet1_idx] : -999.0f")
            .Define("vbs_jet1_eta", "vbs_jet1_idx != -1 ? jet_eta[vbs_jet1_idx] : -999.0f")
            .Define("vbs_jet1_phi", "vbs_jet1_idx != -1 ? jet_phi[vbs_jet1_idx] : -999.0f")
            .Define("vbs_jet1_mass", "vbs_jet1_idx != -1 ? jet_mass[vbs_jet1_idx] : -999.0f")
            .Define("vbs_jet2_pt", "vbs_jet2_idx != -1 ? jet_pt[vbs_jet2_idx] : -999.0f")
            .Define("vbs_jet2_eta", "vbs_jet2_idx != -1 ? jet_eta[vbs_jet2_idx] : -999.0f")
            .Define("vbs_jet2_phi", "vbs_jet2_idx != -1 ? jet_phi[vbs_jet2_idx] : -999.0f")
            .Define("vbs_jet2_mass", "vbs_jet2_idx != -1 ? jet_mass[vbs_jet2_idx] : -999.0f")
            .Define("vbs_mjj", "(vbs_jet1_idx != -1 && vbs_jet2_idx != -1) ? (ROOT::Math::PtEtaPhiMVector(vbs_jet1_pt, vbs_jet1_eta, vbs_jet1_phi, vbs_jet1_mass) + "
                "ROOT::Math::PtEtaPhiMVector(vbs_jet2_pt, vbs_jet2_eta, vbs_jet2_phi, vbs_jet2_mass)).M() : -999.0f")
            .Define("vbs_detajj", "(vbs_jet1_idx != -1 && vbs_jet2_idx != -1) ? abs(vbs_jet1_eta - vbs_jet2_eta) : -999.0f")
            .Define("vbs_candidate_found", "vbs_jet1_idx != -1 && vbs_jet2_idx != -1")
            .Filter("vbs_candidate_found");

        Cutflow::Add(df, "C5: VBS candidate selection");

        df = df.Define("_fatjet_vbs1_dR", VdR, {"fatjet_eta", "fatjet_phi", "vbs_jet1_eta", "vbs_jet1_phi"})
            .Define("_fatjet_vbs2_dR", VdR, {"fatjet_eta", "fatjet_phi", "vbs_jet2_eta", "vbs_jet2_phi"})
            .Define("_boosted_h_candidate_jets", 
                "_fatjet_vbs1_dR >= 0.8 && "
                "_fatjet_vbs2_dR >= 0.8")
            .Define("fatjet_HvsQCD", "fatjet_globalParT3_Xbb / (fatjet_globalParT3_Xbb + fatjet_globalParT3_QCD)")
            .Define("fatjet_VvsQCD", "fatjet_globalParT3_Xqq / (fatjet_globalParT3_Xqq + fatjet_globalParT3_Xcs + fatjet_globalParT3_QCD)")
            .Define("_best_h_idx", "fatjet_HvsQCD.size() != 0 ? ArgMax(fatjet_HvsQCD[_boosted_h_candidate_jets]) : 999.0")
            .Define("boosted_h_candidate_score", "_best_h_idx != 999.0 ? fatjet_HvsQCD[_boosted_h_candidate_jets][_best_h_idx] : -999.0f")
            .Define("boosted_h_candidate_found", "boosted_h_candidate_score > 0")
            .Define("boosted_h_candidate_eta", "boosted_h_candidate_found ? fatjet_eta[_boosted_h_candidate_jets][_best_h_idx] : -999.0f")
            .Define("boosted_h_candidate_phi", "boosted_h_candidate_found ? fatjet_phi[_boosted_h_candidate_jets][_best_h_idx] : -999.0f")
            .Define("boosted_h_candidate_mass", "boosted_h_candidate_found ? fatjet_mass[_boosted_h_candidate_jets][_best_h_idx] : -999.0f")
            .Define("boosted_h_candidate_pt", "boosted_h_candidate_found ? fatjet_pt[_boosted_h_candidate_jets][_best_h_idx] : -999.0f")
            .Define("boosted_h_candidate_tau21", "boosted_h_candidate_found ? fatjet_tau2[_boosted_h_candidate_jets][_best_h_idx] / fatjet_tau1[_boosted_h_candidate_jets][_best_h_idx] : -999.0f")
            .Filter("boosted_h_candidate_found");

        Cutflow::Add(df, "C6: Boosted Higgs candidate selection");

        df = df.Define("_fatjet_h_dR", VdR, {"fatjet_eta", "fatjet_phi", "boosted_h_candidate_eta", "boosted_h_candidate_phi"})
            .Define("_boosted_v_candidate_jets", 
                "_fatjet_h_dR >= 0.8 && "
                "_fatjet_vbs1_dR >= 0.8 && "
                "_fatjet_vbs2_dR >= 0.8")
            .Define("_best_w_idx", "fatjet_VvsQCD.size() != 0 ? ArgMax(fatjet_VvsQCD[_boosted_v_candidate_jets]) : -1")
            .Define("boosted_v_candidate_score", "_best_w_idx != -1 ? fatjet_VvsQCD[_boosted_v_candidate_jets][_best_w_idx] : -999.0f")
            .Define("boosted_v_candidate_found", "boosted_v_candidate_score > 0")
            .Define("boosted_v_candidate_eta", "boosted_v_candidate_found ? fatjet_eta[_boosted_v_candidate_jets][_best_w_idx] : -999.0f")
            .Define("boosted_v_candidate_phi", "boosted_v_candidate_found ? fatjet_phi[_boosted_v_candidate_jets][_best_w_idx] : -999.0f")
            .Define("boosted_v_candidate_mass", "boosted_v_candidate_found ? fatjet_mass[_boosted_v_candidate_jets][_best_w_idx] : -999.0f")
            .Define("boosted_v_candidate_pt", "boosted_v_candidate_found ? fatjet_pt[_boosted_v_candidate_jets][_best_w_idx] : -999.0f")
            .Define("boosted_v_candidate_tau21", "boosted_v_candidate_found ? fatjet_tau2[_boosted_v_candidate_jets][_best_w_idx] / fatjet_tau1[_boosted_v_candidate_jets][_best_w_idx] : -999.0f")
            .Filter("boosted_v_candidate_found");

        Cutflow::Add(df, "C7: Boosted Vector boson candidate selection");
    }

    // 2lepSS
    else if (channel == "2lepSS"){

        df = TriggerSelections(df,trigger_logic_string_multilep);
        Cutflow::Add(df, "C1: Trigger selection");

        df = df.Filter(
            "((nMuon_Loose + nElectron_Loose) == 2)",
            //TODO implement a same sign requirement
            "C2: 2lepSS"
        );
    }

    // 2lep_1FJ (currently shared between OF and SF)
    else if (channel == "2lep_1FJ"){

        df = TriggerSelections(df,trigger_logic_string_multilep);
        Cutflow::Add(df, "C1: Trigger selection");

        df = df.Filter(
            "((nMuon_Loose + nElectron_Loose) == 2) &&"
            "(nFatJets == 1)",
            "C2: 2lep_1FJ"
        );
    }

    // 2lep_2FJ
    else if (channel == "2lep_2FJ"){

        df = TriggerSelections(df,trigger_logic_string_multilep);
        Cutflow::Add(df, "C1: Trigger selection");

        df = df.Filter(
            "((nMuon_Loose + nElectron_Loose) == 2) &&"
            "(nFatJets == 2)",
            "C2: 2lep_2FJ"
        );
    }


    // 3lep
    else if (channel == "3lep"){

        df = TriggerSelections(df,trigger_logic_string_multilep);
        Cutflow::Add(df, "C1: Trigger selection");

        df = df.Filter(
            "((nMuon_Loose + nElectron_Loose) == 3)",
            "C2: 3lep"
        );
    }

    // 4lep
    else if (channel == "4lep"){

        df = TriggerSelections(df,trigger_logic_string_multilep);
        Cutflow::Add(df, "C1: Trigger selection");

        df = df.Filter(
            "((nMuon_Loose + nElectron_Loose) == 4)",
            "C2: 4lep"
        );
    }

    return df;
}
