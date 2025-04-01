#include "commonSelections.h"
#include "selections_OneLep2FJ.h"

namespace OneLep2FJ {
    RNode TriggerSelections(RNode df_) {
        return df_.Define("passes_triggers", "HLT_Ele30_WPTight_Gsf || HLT_IsoMu24");
    }

    RNode ObjectSelections(RNode df_) {
        // lepton
        auto df = df_.Define("isElectron", "nLooseElectrons == 1 && nTightElectrons == 1")
            .Define("isMuon", "nLooseMuons == 1 && nTightMuons == 1")
            .Define("lepton_pt", "isElectron ? electron_pt[0] : muon_pt[0]")
            .Define("lepton_eta", "isElectron ? electron_eta[0] : muon_eta[0]")
            .Define("lepton_phi", "isElectron ? electron_phi[0] : muon_phi[0]")
            .Define("lepton_mass", "isElectron ? electron_mass[0] : muon_mass[0]");

        df = df.Define("xbbscore_idx", "ak8jet_xbbvsqcd.size() != 0 ? ArgMax(ak8jet_xbbvsqcd) : -999.0")
            .Define("xbb_score", "xbbscore_idx != -999.0 ? ak8jet_xbbvsqcd[xbbscore_idx] : -999.0")
            .Define("xbb_pt", "ak8jet_pt[xbbscore_idx]")
            .Define("xbb_eta", "ak8jet_eta[xbbscore_idx]")
            .Define("xbb_phi", "ak8jet_phi[xbbscore_idx]")
            .Define("xbb_msoftdrop", "ak8jet_msoftdrop[xbbscore_idx]");

        // jets
        return df;
    }

    RNode EventSelections(RNode df_) {
        auto df = TriggerSelections(df_);
        df = df.Define("passes_lepton_selection",
            "((nLooseMuons == 1 && nTightMuons == 1 && nLooseElectrons == 0 && nTightElectrons == 0) || (nLooseMuons == 0 && nTightMuons == 0 && nLooseElectrons == 1 && nTightElectrons == 1)) && "
            "(lepton_pt > 40)");
        return df;
    }

    RNode runPreselection(RNode df_) {
        auto df = CommonSelections(df_);
        df = ObjectSelections(df);
        df = EventSelections(df);
        return df;
    }
} // OneLep2FJ