#include "commonSelections.h"
#include "selections_OneLep2FJ.h"

namespace OneLep2FJ {
    TMVA::Experimental::RBDT vbstagger("VBSTagger", "mva/VBSTagger/OneLep2FJ.root");

    RNode TriggerSelections(RNode df_) {
        return df_.Filter("HLT_Ele30_WPTight_Gsf || HLT_IsoMu24");
    }

    RNode runPreselection(RNode df_) {
        // lepton
        auto df = EventFilters(df_);
        df = TriggerSelections(df);

        df = ElectronSelections(df);
        df = MuonSelections(df);

        df = df.Define("_isElectron", "electron_nloose == 1 && electron_ntight == 1")
            .Define("lepton_pt", "_isElectron ? electron_pt[0] : muon_pt[0]")
            .Define("lepton_eta", "_isElectron ? electron_eta[0] : muon_eta[0]")
            .Define("lepton_phi", "_isElectron ? electron_phi[0] : muon_phi[0]")
            .Define("lepton_mass", "_isElectron ? electron_mass[0] : muon_mass[0]")
            .Filter("((muon_nloose == 1 && muon_ntight == 1 && electron_nloose == 0 && electron_ntight == 0) || "
                "(muon_nloose == 0 && muon_ntight == 0 && electron_nloose == 1 && electron_ntight == 1)) && "
                "(lepton_pt > 40)");

        df = AK8JetsSelection(df);
        df = AK4JetsSelection(df);
        df = VBSJetSelections(df, vbstagger);

        df = df.Filter("ak4jet_n >= 2");

        df = df.Define("max_hbb_idx", "ak8jet_hbbvsqcd.size() != 0 ? ArgMax(ak8jet_hbbvsqcd) : -999.0")
            .Define("hbb_score", "max_hbb_idx != -999.0 ? ak8jet_hbbvsqcd[max_hbb_idx] : -999.0")
            .Define("hbb_pt", "ak8jet_pt[max_hbb_idx]")
            .Define("hbb_eta", "ak8jet_eta[max_hbb_idx]")
            .Define("hbb_phi", "ak8jet_phi[max_hbb_idx]")
            .Define("hbb_msoftdrop", "ak8jet_msoftdrop[max_hbb_idx]");
        // jets
        return df;
    }

} // OneLep2FJ