#include "selections.h"
#include "cutflow.h"
#include "weights.h"

namespace {
// AK4 good-jet selection as a string template, parametrised by the pt-column name.
// Used for both the nominal mask (Jet_pt) and the per-source JES variations
// (Jet_pt_jesAbsoluteUp, ...). Must be combined with `&& !Jet_vetoMap` after the veto
// map has been applied (the nominal path does this via Redefine; the variation loop
// appends it explicitly).
inline std::string ak4GoodJetSelectionExpr(const std::string& ptCol) {
    return "_dR_ak4_lep > 0.4 &&"
           "((isRun3 && ((" + ptCol + " > 20 && (abs(Jet_eta) <= 2.5 || abs(Jet_eta) >= 3.0) && abs(Jet_eta) < 5.0) || "
                        "(" + ptCol + " > 50 && abs(Jet_eta) > 2.5 && abs(Jet_eta) < 3.0))) ||"
           "(isRun2 && " + ptCol + " > 20 && abs(Jet_eta) < 5.0)) &&  "
           "(Jet_jetId >= 2)";
}

inline std::string ak8GoodJetSelectionExpr(const std::string& ptCol) {
    return "_dR_ak8_lep > 0.8 && "
           + ptCol + " > 250 && "
           "abs(FatJet_eta) <= 2.5 && "
           "FatJet_msoftdrop > 40 && "
           "FatJet_jetId > 0";
}

// Kinematic-variation suffixes (JES + JER) whose varied jet columns actually exist in
// this graph. Presence-checked rather than assumed: variations are MC-only (and only
// with --systs given), so on data or in the default nominal-only mode this returns empty and the
// selection stays nominal-only. Deriving the set from the columns the correction path
// actually defined means there is no flag to keep in sync with that path.
inline std::vector<std::string> activeKinVariations(RNode df) {
    std::vector<std::string> out;
    for (const auto& sfx : kinematicVariationSuffixes()) {
        if (df.HasColumn("Jet_pt_" + sfx) && df.HasColumn("FatJet_pt_" + sfx))
            out.push_back(sfx);
    }
    return out;
}

// Per-variation channel-pass-flag generation. Each channel filter that depends on a
// jet-multiplicity cut emits a boolean per event for each variation v ∈ {nom} ∪
// activeKinVariations() - passes_<channel>_<v> - then the channel filter is OR over
// those flags. This keeps events that pass the channel under ANY variation, which is
// the only way a single output file can carry per-variation histograms downstream.
//
// `expr_for(sfx)` returns the full pass condition for variation `sfx`. sfx == "" means
// nominal; the lambda should return the nominal condition string in that case.
template <typename F>
RNode definePerVariationPassFlags(RNode df, const std::string& channel, F&& expr_for) {
    df = df.Define("passes_" + channel + "_nom", expr_for(std::string{}));
    for (const auto& sfx : activeKinVariations(df)) {
        df = df.Define("passes_" + channel + "_" + sfx, expr_for(sfx));
    }
    return df;
}

inline std::string orPassExpr(RNode df, const std::string& channel) {
    std::string out = "passes_" + channel + "_nom";
    for (const auto& sfx : activeKinVariations(df)) {
        out += " || passes_" + channel + "_" + sfx;
    }
    return out;
}

// ---------------------------------------------------------------------------
// 1lep_WSF probe gen-truth. Lives here, not in genSelections.cpp, because it reads
// the flattened probe (fj0_eta / fj0_phi) that only the 1lep_WSF channel body
// defines -- the precondition and the thing that establishes it stay together.
// Internal linkage: nothing outside this channel has any use for it.
//
// Strategy: JET-FIRST. Ask what truth object, if any, is contained in the leading
// fat jet -- never start from the sample's expected topology and check containment
// afterwards. The second direction forces an a priori choice exactly where it is
// least safe: in tW->4Q, TTW-WtoQQ, TTWW and TTWZ a top-W and a prompt W are both
// hadronic in the same event, and any fixed priority between them is applied blind.
//
// IDENTITY BEFORE CONTAINMENT. Decide what the jet *is* first; the containment test
// then follows from the identity, rather than one dR threshold standing in for both
// questions. The answer is one code:
//
//   fj0_gen_match   0 nothing   1 prompt W   2 W from top   3 top
//
// 0 covers both "no hadronic W in the record" and "one exists but is not in the
// cone" -- in both cases the jet corresponds to nothing in the hard process, and
// dr_Wq_max (999 vs finite-but->=0.8) separates them.
//
// NanoAOD v15 types, verified against the v30 skims:
//   GenPart_pdgId            Int_t    -> int
//   GenPart_genPartIdxMother Short_t  -> short
//   GenPart_statusFlags      UShort_t -> unsigned short
//   GenPart_status           Int_t    -> int
// isLastCopy is bit 13 of statusFlags.
//
// Downstream (jmsr.selection.categorize):
//   wmerged   = gen_match in {1,2}   (or {2} alone, if prompt W is excluded)
//   tmerged   = gen_match == 3
//   unmatched = otherwise
namespace wsf {

// The AK8 radius. "Inside this jet" has exactly one meaning, so this is a
// definition rather than a tunable working point -- but every continuous dR is
// emitted anyway, so a consumer wanting a different radius never needs a re-run.
constexpr float kMatchDR = 0.8f;

constexpr float kNone = 999.f;

enum : int { kMatchNone = 0, kMatchPromptW = 1, kMatchTopW = 2, kMatchTop = 3 };

} // namespace wsf

ROOT::RVec<float> wsf_truth_vars(
    float probe_eta, float probe_phi,
    const ROOT::RVec<int>& pdg,
    const ROOT::RVec<short>& mom,
    const ROOT::RVec<unsigned short>& flags,
    const ROOT::RVec<int>& status,
    const ROOT::RVec<float>& geta,
    const ROOT::RVec<float>& gphi)
{
    const int n = static_cast<int>(pdg.size());
    auto isLastCopy = [&](int i) { return (flags[i] >> 13) & 1; };
    auto dRprobe    = [&](int i) { return ROOT::VecOps::DeltaR(probe_eta, geta[i], probe_phi, gphi[i]); };

    // Walk the mother chain from i up to a |pdg|==6 top; return its index or -1.
    auto topAncestor = [&](int i) -> int {
        int cur = static_cast<int>(mom[i]);
        for (int guard = 0; cur >= 0 && cur < n && guard < 50; ++guard) {
            if (std::abs(pdg[cur]) == 6) return cur;
            cur = static_cast<int>(mom[cur]);
        }
        return -1;
    };

    // The quark daughters of the boson at index w. A last-copy W decays either to
    // qq' or to lnu, never both -- downstream radiation attaches to the quarks, not
    // back to the W -- so no charged-lepton guard is needed here. (Measured on 240k
    // last-copy W in tW: the multiplicity is exactly 0 or 2, never 1 or 3, and zero
    // W carry both a quark and a charged lepton.) A count != 2 is a record we do not
    // understand, so it is rejected rather than silently truncated to the first two.
    auto quarkDaughters = [&](int w, int& q1, int& q2) -> bool {
        q1 = -1; q2 = -1; int nq = 0;
        for (int j = 0; j < n; ++j) {
            if (static_cast<int>(mom[j]) != w) continue;
            const int ad = std::abs(pdg[j]);
            if (ad < 1 || ad > 5) continue;
            ++nq;
            if      (q1 < 0) q1 = j;
            else if (q2 < 0) q2 = j;
        }
        return nq == 2;
    };

    // ---- choose the candidate, across provenances, by best match to THIS jet ----
    //
    // Only |pdg|==24. A merged hadronic Z is a 91 GeV object and this is a W
    // calibration, so it is not a candidate; it surfaces as a bump in the unmatched
    // template, where it is visible, rather than inside the template whose peak
    // position is the quantity being measured.
    int   best_w = -1, best_q1 = -1, best_q2 = -1;
    float best_dr = 1e9f;
    for (int i = 0; i < n; ++i) {
        if (std::abs(pdg[i]) != 24 || !isLastCopy(i)) continue;
        int q1, q2;
        if (!quarkDaughters(i, q1, q2)) continue;   // leptonic, or a record we reject
        const float dr = std::max(dRprobe(q1), dRprobe(q2));
        if (dr < best_dr) { best_dr = dr; best_w = i; best_q1 = q1; best_q2 = q2; }
    }

    // ---- contamination probe: nearest hard-process b, whatever its origin ----
    //
    // Status 23 is an outgoing particle of the hard process, which catches both a
    // top's b and the initial-state b of 5-flavour-scheme tW. This is what makes
    // "merged W plus a stray b" measurable downstream instead of invisible; it never
    // enters the categorisation.
    float dr_b_near = wsf::kNone;
    for (int j = 0; j < n; ++j) {
        if (std::abs(pdg[j]) != 5 || status[j] != 23) continue;
        dr_b_near = std::min(dr_b_near, dRprobe(j));
    }

    if (best_w < 0) {
        return {static_cast<float>(wsf::kMatchNone), wsf::kNone, wsf::kNone,
                wsf::kNone, dr_b_near, 0.f, 0.f};
    }

    const float dr_Wq_max   = best_dr;
    const float dr_W_daus   = ROOT::VecOps::DeltaR(geta[best_q1], geta[best_q2],
                                                   gphi[best_q1], gphi[best_q2]);
    const int   Wq_max_pdg  = (dRprobe(best_q1) >= dRprobe(best_q2)) ? pdg[best_q1] : pdg[best_q2];

    // ---- the top's partner quark ----
    //
    // t -> q W is one vertex, so the partner is always emitted with the W -- but it
    // is not always a b. Accepting only |pdg|==5 loses CKM-suppressed t->sW / t->dW
    // (measured in tW: 0.1485% s, 0.0071% d, total 0.156% against
    // |Vts|^2 + |Vtd|^2 = 0.171%), and the old code emitted 999 for those, which
    // reads downstream as "the partner is outside the jet" and makes the event a
    // merged W. tt is generated with t->bW forced, which is why this is invisible
    // there. Hence any DOWN-TYPE daughter.
    const int top = topAncestor(best_w);
    float dr_top_q  = wsf::kNone;
    int   top_q_pdg = 0;
    if (top >= 0) {
        for (int j = 0; j < n; ++j) {
            if (static_cast<int>(mom[j]) != top) continue;
            const int ad = std::abs(pdg[j]);
            if (ad != 1 && ad != 3 && ad != 5) continue;
            dr_top_q  = dRprobe(j);
            top_q_pdg = pdg[j];
            break;
        }
    }

    // ---- identity ----
    int match;
    if      (dr_Wq_max >= wsf::kMatchDR) match = wsf::kMatchNone;    // found, not contained
    else if (top < 0)                    match = wsf::kMatchPromptW;
    else if (dr_top_q <= wsf::kMatchDR)  match = wsf::kMatchTop;
    else                                 match = wsf::kMatchTopW;

    return {static_cast<float>(match), dr_Wq_max, dr_top_q, dr_W_daus, dr_b_near,
            static_cast<float>(Wq_max_pdg), static_cast<float>(top_q_pdg)};
}

} // anonymous namespace

// MET filters
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
RNode METFilters(RNode df_) {
    return df_.Filter("(isRun3 && (Flag_goodVertices && "
            "Flag_globalSuperTightHalo2016Filter && "
            "Flag_EcalDeadCellTriggerPrimitiveFilter && "
            "Flag_BadPFMuonFilter && "
            "Flag_BadPFMuonDzFilter && "
            "Flag_hfNoisyHitsFilter &&"
            "Flag_eeBadScFilter && "
            "Flag_ecalBadCalibFilter)) || "
        "(isRun2 && (Flag_goodVertices && "
            "Flag_globalSuperTightHalo2016Filter && "
            "Flag_HBHENoiseFilter && "
            "Flag_HBHENoiseIsoFilter && "
            "Flag_EcalDeadCellTriggerPrimitiveFilter && "
            "Flag_BadPFMuonFilter &&"
            "Flag_BadPFMuonDzFilter && "
            "Flag_eeBadScFilter && "
            "(!(is2017 || is2018) || Flag_ecalBadCalibFilter)))"
        );
}

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


// Ele selection
RNode ElectronSelections(RNode df_)
{
    auto df = df_.Define("Electron_SC_eta", "Electron_eta + Electron_deltaEtaSC");

    // Veto selection
    df = df.Define(
        "_vetoElectrons",
        "Electron_pt > 10 &&"
        "abs(Electron_SC_eta) < 2.5 && "
        "((abs(Electron_SC_eta) <= 1.479 && abs(Electron_dxy) <= 0.05 && abs(Electron_dz) < 0.1) || ((abs(Electron_SC_eta) > 1.479) && abs(Electron_dxy) <= 0.1 && abs(Electron_dz) < 0.2)) && "
        "Electron_cutBased >= 1"
    );

    // Loose selection
    df = df.Define(
        "_looseElectrons",
        "_vetoElectrons &&"
        "Electron_cutBased >= 2"
    );

    // We will write out the electron object
    df = applyObjectMaskNewAffix(df, "_vetoElectrons", "Electron", "electron");

    // Append the masks to electron object
    df = df.Define("electron_isLoose",  "electron_cutBased >=2");
    df = df.Define("electron_isMedium", "electron_cutBased >=3");
    df = df.Define("electron_isTight",  "electron_cutBased >=4");

    // Define the counts of each
    df = df.Define("nElectron_Veto", "nElectron == 0 ? 0 : Sum(_vetoElectrons)");
    df = df.Define("nElectron_Loose", "Sum(electron_isLoose)");
    df = df.Define("nElectron_Tight", "Sum(electron_isTight)");

    return df;
}

// Muon selections
RNode MuonSelections(RNode df_)
{

    // Loose selection
    auto df = df_.Define(
        "_looseMuons",
        "Muon_pt > 10 && "
        "Muon_pfIsoId >= 2 && "
        "abs(Muon_eta) < 2.4 && "
        "abs(Muon_dxy) < 0.2 && "
        "abs(Muon_dz) < 0.5 && "
        "abs(Muon_sip3d) < 8 && "
        "Muon_looseId"
    );

    // Medium selection
    df = df.Define(
        "_mediumMuons",
        "_looseMuons && "
        "Muon_pfIsoId >= 3 && "
        "Muon_mediumId"
    );

    // We will write out the muon object
    df = applyObjectMaskNewAffix(df, "_looseMuons", "Muon", "muon");

    // Append the masks to electron object
    df = df.Define("muon_isMedium",    "muon_mediumId && (muon_pfIsoId >=3)");
    df = df.Define("muon_isTight",     "muon_mediumId && (muon_pfIsoId >=4)");
    df = df.Define("muon_isVeryTight", "muon_mediumId && (muon_pfIsoId >=5)");

    // Define the counts of each
    df = df.Define("nMuon_Loose", "nMuon == 0 ? 0 : Sum(_looseMuons)");
    df = df.Define("nMuon_Medium", "Sum(muon_isMedium)");
    df = df.Define("nMuon_Tight", "Sum(muon_isTight)");

    return df;
}

// Put ele and mu together into lepton collection
RNode LeptonSelections(RNode df_)
{
    auto df = ElectronSelections(df_);
    df = MuonSelections(df);

    // Counting lepton collections: used for SF calculation and channel counting
    //     - These are the leptons that define the channel (loose WP for both flavors)
    //     - Prefixed with _ so they are not written to the output snapshot
    //     - Muons: loose muons (muon_* collection is already loose)
    //     - Electrons: loose electrons (i.e., loose subset of the veto electron collection)
    df = df.Define("_muonSel_pt",  "muon_pt")
           .Define("_muonSel_eta", "muon_eta");
    df = df.Define("_electronSel_pt",     "electron_pt[electron_isLoose]")
           .Define("_electronSel_eta",    "electron_eta[electron_isLoose]")
           .Define("_electronSel_SC_eta", "electron_SC_eta[electron_isLoose]");
    // Counter for these selection level leptons
    df = df.Define("nLep_Sel", "nMuon_Loose + nElectron_Loose");

    return df.Define("lepton_pt", "Concatenate(electron_pt, muon_pt)")
            .Define("_leptonSorted", "Argsort(-lepton_pt)")
            .Redefine("lepton_pt", "Take(lepton_pt, _leptonSorted)")
            .Define("lepton_eta", "Take(Concatenate(electron_eta, muon_eta), _leptonSorted)")
            .Define("lepton_phi", "Take(Concatenate(electron_phi, muon_phi), _leptonSorted)")
            .Define("lepton_mass", "Take(Concatenate(electron_mass, muon_mass), _leptonSorted)")
            .Define("lepton_charge", "Take(Concatenate(electron_charge, muon_charge), _leptonSorted)");
}


// Small-radius AK4 jet selection
RNode AK4JetsSelection(RNode df_, bool cleanAgainstFJ, std::string affix)
{
    auto df = DefineOrRedefine(df_, "_dR_ak4_lep", VVdR, {"Jet_eta", "Jet_phi", "lepton_eta", "lepton_phi"});

    // Optionally define the dR between AK4 jets and fat jets, and include it in the AK4 selection if cleanAgainstFJ is set
    std::string fj_cut = "";
    if (cleanAgainstFJ) {
        df = df.Define("_dR_ak4_fatjet", VVdR, {"Jet_eta", "Jet_phi", "fatjet_eta", "fatjet_phi"});
        fj_cut = "_dR_ak4_fatjet > 0.8 && ";
    }

    df = DefineOrRedefine(df, "_good_ak4jets", "_dR_ak4_lep > 0.4 && " + fj_cut +
                                    "((isRun3 && ((Jet_pt > 20 && (abs(Jet_eta) <= 2.5 || abs(Jet_eta) >= 3.0) && abs(Jet_eta) < 5.0) || "
                                    "(Jet_pt > 50 && abs(Jet_eta) > 2.5 && abs(Jet_eta) < 3.0))) ||"
                                    "(isRun2 && Jet_pt > 20 && abs(Jet_eta) < 5.0)) &&  "
                                    "(Jet_jetId >= 2)"); // NanoAOD jetID convention https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD#Jets
                                                        // should still work for skims < v28, which set jetId==3 : "pass tight ID, fail tightLepVeto", jetId==7 : "pass tight and tightLepVeto ID"

    // Jet veto maps. Run 3 drops the whole event (the veto exists to remove the spurious
    // MET a bad-region jet induces, which masking the jet alone would not undo); Run 2
    // drops only the jets. Both are skipped under --no_jetveto, which leaves Jet_vetoMap
    // defined but unapplied -- DEBUGGING ONLY, see corrections.h.
    const std::string vetoTerm = applyJetVetoMapsEnabled() ? " && !Jet_vetoMap" : "";
    if (applyJetVetoMapsEnabled()) {
        df = df.Filter("(isRun2) || (isRun3 && !Any(Jet_vetoMap))");
        df = df.Redefine("_good_ak4jets", "_good_ak4jets && !Jet_vetoMap");
    }

    df = applyObjectMaskNewAffix(df, "_good_ak4jets", "Jet", affix);
    df = df.Define("ht_" + affix, "Sum(" + affix + "_pt)");

    // Nominal + per-variation good-jet masks stored at the NanoAOD (capital-collection) level.
    // Defined once, alongside the primary FJ-cleaned "jet" collection, and AFTER
    // applyObjectMaskNewAffix so they are not aliased into lowercase jet_isGood columns.
    // Downstream coffea uses Jet_isGood (nominal) and Jet_isGood_<sfx> (variations)
    // uniformly to construct jet collections from the shared Jet_* arrays.
    //
    // Each per-variation mask is fat-jet-cleaned against the good fat jets FOR THAT
    // variation (FatJet_isGood_<sfx>), so Jet_isGood_<sfx> parallels the FJ-cleaned nominal
    // Jet_isGood (= _good_ak4jets for the "jet" affix). Only the good-fatjet membership
    // changes with the variation; the ak4-ak8 dR geometry itself is JEC-invariant (eta/phi
    // are unchanged by JEC). VVdR returns 999 for every AK4 jet when a variation has no good
    // fat jets, so the >0.8 cut then passes (no cleaning) -- matching the nominal fj_cut.
    //
    // njet_<sfx> = Sum(Jet_isGood_<sfx>) then parallels the nominal `njet` (defined inside
    // applyObjectMaskNewAffix, also FJ-cleaned) and feeds the per-variation channel-pass flags.
    if (affix == "jet") {
        df = df.Define("Jet_isGood", "_good_ak4jets");
        for (const auto& sfx : activeKinVariations(df)) {
            df = df.Define("_fatjet_eta_" + sfx, "FatJet_eta[FatJet_isGood_" + sfx + "]");
            df = df.Define("_fatjet_phi_" + sfx, "FatJet_phi[FatJet_isGood_" + sfx + "]");
            df = df.Define("_dR_ak4_fatjet_" + sfx, VVdR,
                           {"Jet_eta", "Jet_phi", "_fatjet_eta_" + sfx, "_fatjet_phi_" + sfx});
            df = df.Define("Jet_isGood_" + sfx,
                           ak4GoodJetSelectionExpr("Jet_pt_" + sfx)
                               + vetoTerm + " && _dR_ak4_fatjet_" + sfx + " > 0.8");
            df = df.Define("njet_" + sfx, "Sum(Jet_isGood_" + sfx + ")");
        }
    }
    // Nominal + per-variation masks for the non-FJ-cleaned collection, used by
    // VBSTagging in channels that tag on "jetNoFJClean" (1lep_2FJ). No FJ cleaning to
    // re-evaluate, so the mask is just the pt-varied good-jet expression + veto map.
    // The Jet_ prefix means the masks are written to the output (snapshot keeps all
    // Jet_* columns), mirroring Jet_isGood / Jet_isGood_<sfx> for the cleaned set.
    else if (affix == "jetNoFJClean") {
        df = df.Define("Jet_isGoodNoFJClean", "_good_ak4jets");
        for (const auto& sfx : activeKinVariations(df)) {
            df = df.Define("Jet_isGoodNoFJClean_" + sfx,
                           ak4GoodJetSelectionExpr("Jet_pt_" + sfx) + vetoTerm);
            df = df.Define("njetNoFJClean_" + sfx, "Sum(Jet_isGoodNoFJClean_" + sfx + ")");
        }
    }
    return df;
}

// Define per-jet properties (b-tagging, veto maps) independent of jet selection
RNode AK4JetProperties(RNode df_)
{
    df_ = applyJetVetoMaps(df_);
    return df_.Define("Jet_isTightBTag", isbTagTight, {"year", "Jet_btagUParTAK4B"})
              .Define("Jet_isMediumBTag", isbTagMedium, {"year", "Jet_btagUParTAK4B"})
              .Define("Jet_isLooseBTag", isbTagLoose, {"year", "Jet_btagUParTAK4B"})
              .Define("Jet_isExtraTightBTag", isbTagExtraTight, {"year", "Jet_btagUParTAK4B"})
              .Define("Jet_isExtraExtraTightBTag", isbTagExtraExtraTight, {"year", "Jet_btagUParTAK4B"});
}


// Fat jet AK8 selection
RNode AK8JetsSelection(RNode df_)
{
    auto df = df_.Define("_dR_ak8_lep", VVdR, {"FatJet_eta", "FatJet_phi", "lepton_eta", "lepton_phi"})
                  .Define("_good_ak8jets", ak8GoodJetSelectionExpr("FatJet_pt"));

    // GloParT regressed mass: the calibrated AK8 mass used by the analysis. massCorrX2p is a
    // ratio predicting the particle-level mass from the RAW jet, so it multiplies the raw mass
    // snapshotted before JEC (_FatJet_rawmass, from snapshotRawJetKinematics in corrections.cpp)
    // rather than the JEC/JER-corrected FatJet_mass. It replaces the JEC/JER mass calibration
    // for AK8 and is not stacked on top of it. Defined before applyObjectMaskNewAffix so the
    // masked fatjet_massGloParT3 alias is created with the rest of the collection.
    df = df.Define("FatJet_massGloParT3", "FatJet_globalParT3_massCorrX2p * _FatJet_rawmass");

    // GloParT tagging discriminants, built from the raw score outputs as signal / (signal + QCD):
    //   Hbb = Xbb / (Xbb + QCD)
    //   Wqq = (Xqq/3 + Xcs) / (Xqq/3 + Xcs + QCD)
    // Unlike the mass these are kinematics-independent, so they need no raw-quantity snapshot.
    // Jets with a vanishing denominator (no score at all) get -1 rather than a NaN.
    auto vsQCD = [](const RVec<float>& sig, const RVec<float>& qcd) {
        RVec<float> out(sig.size(), -1.0f);
        for (size_t i = 0; i < sig.size(); ++i) {
            const float den = sig[i] + qcd[i];
            if (den > 0.0f) out[i] = sig[i] / den;
        }
        return out;
    };
    df = df.Define("_FatJet_Vqq_score", "FatJet_globalParT3_Xqq / 3.0f + FatJet_globalParT3_Xcs")
           .Define("FatJet_Hbb", vsQCD, {"FatJet_globalParT3_Xbb", "FatJet_globalParT3_QCD"})
           .Define("FatJet_Wqq", vsQCD, {"_FatJet_Vqq_score",     "FatJet_globalParT3_QCD"});

    df = applyObjectMaskNewAffix(df, "_good_ak8jets", "FatJet", "fatjet");
    df = df.Define("ht_fatjets", "Sum(fatjet_pt)");

    // Nominal good-fat-jet mask - defined AFTER applyObjectMaskNewAffix so it is not
    // aliased into fatjet_isGood. FatJet_isGood (nominal) mirrors the per-variation
    // FatJet_isGood_<sfx> branches; coffea uses them uniformly. The nominal count is
    // already available as "nfatjet" (from applyObjectMaskNewAffix above), matching the
    // AK4 convention (njet / njet_<sfx>), so no separate nominal count is defined here.
    df = df.Define("FatJet_isGood", "_good_ak8jets");

    // Per-variation good-fatjet masks + counts. JES/JER affect FatJet_pt only;
    // FatJet_msoftdrop is unchanged. Defined AFTER applyObjectMaskNewAffix for the same
    // reason as the AK4 case. Named nfatjet_<sfx> to parallel the nominal "nfatjet" and
    // AK4's njet / njet_<sfx>.
    for (const auto& sfx : activeKinVariations(df)) {
        df = df.Define("FatJet_isGood_" + sfx, ak8GoodJetSelectionExpr("FatJet_pt_" + sfx));
        df = df.Define("nfatjet_" + sfx, "Sum(FatJet_isGood_" + sfx + ")");
    }
    return df;
}

// Perform VBS jet tagging via VBSBDTInfer and compute VBS pair kinematics
// The nominal or a variation with < 2 taggable jets gets sentinel values instead (indices -1, score/kinematics -999).
RNode VBSTagging(RNode df_, std::string jetCollectionName = "jet")
{
    // All-jet-level mask stem for the tagging collection. The per-variation masks are
    // defined in AK4JetsSelection for both stems.
    const bool fjCleaned = (jetCollectionName == "jet");
    const std::string maskStem = fjCleaned ? "Jet_isGood" : "Jet_isGoodNoFJClean";

    // Continue to store the full four vector of the two VBS candidates when running
    // on nominal for backwards compatibility.
    auto df = df_
        .Define("_vbs_candidate_jet_pairs", VBSBDTInfer, {jetCollectionName + "_pt", jetCollectionName + "_eta", jetCollectionName + "_phi", jetCollectionName + "_mass", "isRun2"})
        .Define("vbs_jet1_idx", "static_cast<int>(_vbs_candidate_jet_pairs[0])") // indices in the "good jet" collection
        .Define("vbs_jet2_idx", "static_cast<int>(_vbs_candidate_jet_pairs[1])")
        .Define("vbs_score", "vbs_jet1_idx >= 0 ? _vbs_candidate_jet_pairs[2] : -999.f")
        .Define("vbs_jet1_pt",   "vbs_jet1_idx >= 0 ? " + jetCollectionName + "_pt[vbs_jet1_idx]   : -999.f")
        .Define("vbs_jet1_eta",  "vbs_jet1_idx >= 0 ? " + jetCollectionName + "_eta[vbs_jet1_idx]  : -999.f")
        .Define("vbs_jet1_phi",  "vbs_jet1_idx >= 0 ? " + jetCollectionName + "_phi[vbs_jet1_idx]  : -999.f")
        .Define("vbs_jet1_mass", "vbs_jet1_idx >= 0 ? " + jetCollectionName + "_mass[vbs_jet1_idx] : -999.f")
        .Define("vbs_jet2_pt",   "vbs_jet2_idx >= 0 ? " + jetCollectionName + "_pt[vbs_jet2_idx]   : -999.f")
        .Define("vbs_jet2_eta",  "vbs_jet2_idx >= 0 ? " + jetCollectionName + "_eta[vbs_jet2_idx]  : -999.f")
        .Define("vbs_jet2_phi",  "vbs_jet2_idx >= 0 ? " + jetCollectionName + "_phi[vbs_jet2_idx]  : -999.f")
        .Define("vbs_jet2_mass", "vbs_jet2_idx >= 0 ? " + jetCollectionName + "_mass[vbs_jet2_idx] : -999.f")
        .Define("vbs_mjj",    "vbs_jet1_idx >= 0 ? (ROOT::Math::PtEtaPhiMVector(vbs_jet1_pt, vbs_jet1_eta, vbs_jet1_phi, vbs_jet1_mass) + "
                              "ROOT::Math::PtEtaPhiMVector(vbs_jet2_pt, vbs_jet2_eta, vbs_jet2_phi, vbs_jet2_mass)).M() : -999.")
        .Define("vbs_detajj", "vbs_jet1_idx >= 0 ? abs(vbs_jet1_eta - vbs_jet2_eta) : -999.");

    // Nominal VBS-jet indices in the all-jet (NanoAOD Jet_*) collection, corresponding
    // to the same physical jets as vbs_jet1/2_idx; stored so downstream
    // handles nominal and variations uniformly via Jet_*[vbs_jet*_Jetidx*].
    df = df.Define("vbs_jet1_Jetidx", "vbs_jet1_idx >= 0 ? static_cast<int>(ROOT::VecOps::Nonzero(" + maskStem + ")[vbs_jet1_idx]) : -1")
           .Define("vbs_jet2_Jetidx", "vbs_jet2_idx >= 0 ? static_cast<int>(ROOT::VecOps::Nonzero(" + maskStem + ")[vbs_jet2_idx]) : -1");

    // Per-variation VBS tagging: re-run the BDT on the variation's good-jet set with its
    // varied pt/mass (eta/phi are JEC-invariant). Only the two all-jet-frame indices and
    // the score are stored per variation - downstream can rebuild the kinematics from
    // Jet_*[vbs_jet*_Jetidx*], keeping the branch count small.
    for (const auto& sfx : activeKinVariations(df)) {
        df = df.Define("_vbs_pair_" + sfx, VBSBDTInferMasked,
                       {maskStem + "_" + sfx, "Jet_pt_" + sfx, "Jet_eta", "Jet_phi", "Jet_mass_" + sfx, "isRun2"});
        df = df.Define("vbs_jet1_Jetidx_" + sfx, "static_cast<int>(_vbs_pair_" + sfx + "[0])");
        df = df.Define("vbs_jet2_Jetidx_" + sfx, "static_cast<int>(_vbs_pair_" + sfx + "[1])");
        df = df.Define("vbs_score_" + sfx, "_vbs_pair_" + sfx + "[2]");
    }
    return df;
}



///////////////// Main channel selection block /////////////////
RNode runPreselection(RNode df_, std::string channel, bool noCut, bool isData)
{

    Cutflow::Add(df_, "All events");

    // Standard MET filters
    auto df = METFilters(df_);
    Cutflow::Add(df, "C0: MET filters");

    // Perform lepton selection first (these have highest priority)
    df = LeptonSelections(df);

    // Fat jet selection (cleaned against leptons)
    df = AK8JetsSelection(df);

    // Do the ak4 jet selection
    df = AK4JetProperties(df);
    df = AK4JetsSelection(df, /*cleanAgainstFJ=*/true, "jet");
    df = AK4JetsSelection(df, /*cleanAgainstFJ=*/false, "jetNoFJClean");

    // Get ahold of the relevant lepton SFs (will choose which ones to actually apply in each channel)
    df = applyMuonWorkingPointSFs(df, isData, {
        "_weight_muon_looseid_looseiso",
        "_weight_muon_mediumid_tightiso",
        "_weight_muon_tightid_tightiso"
    });
    df = applyElectronWorkingPointSFs(df, isData, {
        "_weight_electron_reco_looseid",
        "_weight_electron_reco_tightid"
    });


    // Passthrough
    if (noCut) return df;
    if (channel == "all_events"){
        df = df.Filter(
            "nMuon_Loose > -1", // Probably there is a better way to write a pass through
            "C2: all_events"
        );
    }

    ////////// Each channel selections //////////

    // 0lep_0FJ
    else if (channel == "0lep_0FJ"){

        df = VBSTagging(df);
        Cutflow::Add(df, "VBS pair candidate found");

        df = TriggerSelections(df,trigger_logic_string_0lep0FJ);
        Cutflow::Add(df, "C1: Trigger selection");

        // Channel orthogonality selection
        df = definePerVariationPassFlags(df, "0lep_0FJ", [](const std::string& sfx){
            const std::string fjCountCol = sfx.empty() ? "nfatjet" : "nfatjet_" + sfx;
            return "(nLep_Sel == 0) && (" + fjCountCol + " == 0)";
        });
        df = df.Filter(orPassExpr(df, "0lep_0FJ"), "C2: 0lep_0FJ");

        df = df.Define("met_significance", "PuppiMET_significance")
                .Define("met_uncorrPt", "PuppiMET_pt")
                .Define("met_uncorrPhi", "PuppiMET_phi");
    }

    // 0lep_1FJ
    else if (channel == "0lep_1FJ"){

        df = VBSTagging(df);
        Cutflow::Add(df, "VBS pair candidate found");

        df = TriggerSelections(df,trigger_logic_string_0lep1FJ);
        Cutflow::Add(df, "C1: Trigger selection");

        // Channel orthogonality selection
        df = definePerVariationPassFlags(df, "0lep_1FJ", [](const std::string& sfx){
            const std::string fjCountCol = sfx.empty() ? "nfatjet" : "nfatjet_" + sfx;
            return "(nLep_Sel == 0) && (" + fjCountCol + " == 1)";
        });
        df = df.Filter(orPassExpr(df, "0lep_1FJ"), "C2: 0lep_1FJ");
    }

    // 0lep_1FJ_met
    else if (channel == "0lep_1FJ_met"){

        df = VBSTagging(df);
        Cutflow::Add(df, "VBS pair candidate found");

        // HLT_PFMETNoMu120_PFMHTNoMu120_IDTight is required here, not deferred to
        // the offline analysis. Pass --store_hlt to also write out the remaining
        // MET path decisions for trigger studies; they cannot widen this channel.
        df = TriggerSelections(df,trigger_logic_string_met);
        Cutflow::Add(df, "C1: Trigger selection");

        // Channel orthogonality selection
        df = definePerVariationPassFlags(df, "0lep_1FJ_met", [](const std::string& sfx){
            const std::string fjCountCol = sfx.empty() ? "nfatjet" : "nfatjet_" + sfx;
            return "(nLep_Sel == 0) && (" + fjCountCol + " == 1)";
        });
        df = df.Filter(orPassExpr(df, "0lep_1FJ_met"), "C2: 0lep_1FJ");

        df = df.DefaultValueFor("Pileup_nTrueInt", (float)-1.f)
                .Define("met_significance", "PuppiMET_significance")
                .Define("met_uncorrPt", "PuppiMET_pt")
                .Define("met_uncorrPhi", "PuppiMET_phi")
                .Define("pileup_nTrueInt", "Pileup_nTrueInt")
                .Define("pv_npvsGood", "PV_npvsGood");
    }

    // 0lep_2FJ
    else if (channel == "0lep_2FJ"){

        df = VBSTagging(df);
        Cutflow::Add(df, "VBS pair candidate found");

        df = TriggerSelections(df,trigger_logic_string_ht);
        Cutflow::Add(df, "C1: Trigger selection");

        // Channel orthogonality selection
        df = definePerVariationPassFlags(df, "0lep_2FJ", [](const std::string& sfx){
            const std::string fjCountCol = sfx.empty() ? "nfatjet" : "nfatjet_" + sfx;
            return "(nLep_Sel == 0) && (" + fjCountCol + " == 2)";
        });
        df = df.Filter(orPassExpr(df, "0lep_2FJ"), "C2: 0lep_2FJ");
    }

    // 0lep_2FJ_met
    else if (channel == "0lep_2FJ_met"){

        df = VBSTagging(df);
        Cutflow::Add(df, "VBS pair candidate found");

        // HLT_PFMETNoMu120_PFMHTNoMu120_IDTight is required here, not deferred to
        // the offline analysis. Pass --store_hlt to also write out the remaining
        // MET path decisions for trigger studies; they cannot widen this channel.
        df = TriggerSelections(df,trigger_logic_string_met);
        Cutflow::Add(df, "C1: Trigger selection");

        // Channel orthogonality selection
        df = definePerVariationPassFlags(df, "0lep_2FJ_met", [](const std::string& sfx){
            const std::string fjCountCol = sfx.empty() ? "nfatjet" : "nfatjet_" + sfx;
            return "(nLep_Sel == 0) && (" + fjCountCol + " == 2)";
        });
        df = df.Filter(orPassExpr(df, "0lep_2FJ_met"), "C2: 0lep_2FJ");

        df = df.DefaultValueFor("Pileup_nTrueInt", (float)-1.f)
                .Define("met_significance", "PuppiMET_significance")
                .Define("met_uncorrPt", "PuppiMET_pt")
                .Define("met_uncorrPhi", "PuppiMET_phi")
                .Define("pileup_nTrueInt", "Pileup_nTrueInt")
                .Define("pv_npvsGood", "PV_npvsGood");
    }

    // 0lep_3FJ
    else if (channel == "0lep_3FJ"){

        df = VBSTagging(df);
        Cutflow::Add(df, "VBS pair candidate found");

        df = TriggerSelections(df,trigger_logic_string_ht);
        Cutflow::Add(df, "C1: Trigger selection");

        // Channel orthogonality selection
        df = definePerVariationPassFlags(df, "0lep_3FJ", [](const std::string& sfx){
            const std::string fjCountCol = sfx.empty() ? "nfatjet" : "nfatjet_" + sfx;
            return "(nLep_Sel == 0) && (" + fjCountCol + " == 3)";
        });
        df = df.Filter(orPassExpr(df, "0lep_3FJ"), "C2: 0lep_3FJ");
    }

    // 1lep_1FJ - fatjet + njet cuts must combine inside a single per-variation pass flag,
    // because OR-ing over two consecutive jet-count filters would mix variations across
    // different cuts. Cutflow C3 + C4 (separate fatjet/njet steps) collapse into a single
    // "any-variation jet selection" entry.
    else if (channel == "1lep_1FJ"){

        df = VBSTagging(df);
        Cutflow::Add(df, "VBS pair candidate found");

        df = TriggerSelections(df,trigger_logic_string_singlelep);
        Cutflow::Add(df, "C1: Trigger selection");

        df = lepSFWrapper(df,isData, /*ele_sf_name=*/ "_weight_electron_reco_tightid", /*muo_sf_name=*/ "_weight_muon_tightid_tightiso", /*include_trigger_sf=*/ true);

        // Channel orthogonality selection
        df = df.Filter("(nLep_Sel == 1)", "C2: 1lep_1FJ");

        df = df.Filter("((nMuon_Loose == 1 && nMuon_Tight == 1 && nElectron_Loose == 0 && nElectron_Tight == 0) || "
                       "(nMuon_Loose == 0 && nMuon_Tight == 0 && nElectron_Loose == 1 && nElectron_Tight == 1)) && "
                       "(lepton_pt[0] > 40)");
        Cutflow::Add(df, "C2: 1-lepton selection");

        df = definePerVariationPassFlags(df, "1lep_1FJ", [](const std::string& sfx){
            const std::string fjCountCol = sfx.empty() ? "nfatjet" : "nfatjet_" + sfx;
            const std::string jCountCol  = sfx.empty() ? "njet"     : "njet_"     + sfx;
            return "(" + fjCountCol + " == 1) && (" + jCountCol + " >= 4)";
        });
        df = df.Filter(orPassExpr(df, "1lep_1FJ"), "C3: jet selection (any variation)");
    }

    // 1lep_2FJ - same caveat as 1lep_1FJ (C3 + C4 collapsed).
    else if (channel == "1lep_2FJ"){

        df = VBSTagging(df, "jetNoFJClean");
        Cutflow::Add(df, "VBS pair candidate found");

        df = TriggerSelections(df,trigger_logic_string_singlelep);
        Cutflow::Add(df, "C1: Trigger selection");

        df = lepSFWrapper(df,isData, /*ele_sf_name=*/ "_weight_electron_reco_tightid", /*muo_sf_name=*/ "_weight_muon_tightid_tightiso", /*include_trigger_sf=*/ true);

        // Channel orthogonality selection
        df = df.Filter("(nLep_Sel == 1)", "C2: 1lep_2FJ");

        df = df.Filter("((nMuon_Loose == 1 && nMuon_Tight == 1 && nElectron_Veto == 0 && nElectron_Loose == 0 && nElectron_Tight == 0) || "
                       "(nMuon_Loose == 0 && nMuon_Tight == 0 && nElectron_Veto == 1 && nElectron_Loose == 1 && nElectron_Tight == 1)) && "
                       "(lepton_pt[0] > 40)");

        Cutflow::Add(df, "C2: 1-lepton selection");

        df = definePerVariationPassFlags(df, "1lep_2FJ", [](const std::string& sfx){
            const std::string fjCountCol = sfx.empty() ? "nfatjet" : "nfatjet_" + sfx;
            const std::string jCountCol  = sfx.empty() ? "njet"     : "njet_"     + sfx;
            return "(" + fjCountCol + " >= 2) && (" + jCountCol + " >= 2)";
        });
        df = df.Filter(orPassExpr(df, "1lep_2FJ"), "C3: jet selection (any variation)");
    }

    // 1lep_WSF - GloParT W-tagging and JMS/JMR calibration region.
    //
    // Standalone tag-and-probe measurement region whose output feeds the jmsr fit and the Wqq
    // tagging-SF fit. Because it is not part of the analysis, run_rdf.py keeps it
    // out of the "--channels all" expansion -- see NON_ANALYSIS_CHANNELS there.
    //
    else if (channel == "1lep_WSF"){

        df = TriggerSelections(df, trigger_logic_string_1lepWSF);
        Cutflow::Add(df, "C1: single-muon trigger");

        // nMuon_Tight is mediumId && pfIsoId>=4 (MuonSelections above), so the matching
        // SF payload is mediumid_tightiso, NOT the tightid_tightiso the 1lep_* channels use.
        df = lepSFWrapper(df, isData,
                          /*ele_sf_name=*/ "_weight_electron_reco_tightid",
                          /*muo_sf_name=*/ "_weight_muon_mediumid_tightiso",
                          /*include_trigger_sf=*/ true);

        // Tag: exactly one muon in the event, and it is tight. nMuon_Tight alone would
        // not veto a second, softer loose muon, and muon_pt[0] indexes the loose
        // collection -- so nMuon_Loose == 1 is what makes [0] the tag muon.
        df = df.Filter("nMuon_Loose == 1 && nMuon_Tight == 1 && "
                       "nElectron_Loose == 0 && muon_pt[0] > 30", "C2: single muon tag");
        Cutflow::Add(df, "C2: single muon tag (pT>30)");

        // Probe. Hard nominal cut, NOT definePerVariationPassFlags/orPassExpr: the
        // probe payload below indexes fatjet_*[0], which is only defined once this has
        // passed nominally. An OR-over-variations filter would admit events with
        // nfatjet == 0 nominally and read out of bounds. The cost is that under --systs
        // the flat fj0_* columns stay nominal; the JES-varied FatJet_* arrays and
        // FatJet_isGood_<sfx> masks are still written, so a varied probe is
        // reconstructible downstream.
        df = df.Filter("nfatjet >= 1", "C3: >=1 fat jet (probe)");
        Cutflow::Add(df, "C3: >=1 fat jet");

        // Type-1 rebuilt MET. The W-SF fit re-cuts at 40 downstream.
        df = df.Filter("met_pt > 30", "C4: MET > 30");
        Cutflow::Add(df, "C4: MET > 30");

        // Suppresses W+jets, whose untagged low-mass tail is what rails the jmsr fit.
        df = df.Filter("Sum(jet_isMediumBTag) >= 1", "C5: >=1 medium b-tag");
        Cutflow::Add(df, "C5: >=1 medium b-tag");

        // ---- probe payload: Defines only, no further cuts (produce loose, cut downstream) ----
        // Probe = leading good AK8, already dR>0.8 isolated from the lepton by AK8JetsSelection.
        df = df.Define("fj0_pt",     "fatjet_pt[0]")
               .Define("fj0_eta",    "fatjet_eta[0]")
               .Define("fj0_phi",    "fatjet_phi[0]")
               .Define("fj0_sdmass", "fatjet_msoftdrop[0]")
               .Define("fj0_massGloParT3", "fatjet_massGloParT3[0]")
               .Define("fj0_Wqq", "fatjet_Wqq[0]")
               .Define("fj0_Hbb", "fatjet_Hbb[0]")
               // Tag-side b multiplicity. jet_* is the FJ-CLEANED AK4 collection
               // (dR>0.8 from every good fat jet, AK4JetsSelection above), so this counts
               // b-tags OUTSIDE the probe -- it is not the b content of the fat jet, which
               // is what fj0_dr_T_b gives in MC. Named without the fj0_ prefix for that
               // reason. Cast to int because Sum(RVec<bool>) returns std::size_t, whose
               // branch type uproot has no interpretation rule for.
               .Define("nMediumBTag_outsideFJ", "(int)Sum(jet_isMediumBTag)");

        // Leptonic-W pT: the tag-side boost. In semileptonic tt the two tops recoil
        // against each other, so requiring a boosted leptonic W (the fit cuts at 150
        // downstream) is what makes the probe AK8 contain a genuinely merged hadronic W.
        df = df.Define("_lepW_px", "muon_pt[0]*cos(muon_phi[0]) + met_pt*cos(met_phi)")
               .Define("_lepW_py", "muon_pt[0]*sin(muon_phi[0]) + met_pt*sin(met_phi)")
               .Define("leptonicW_pt", "sqrt(_lepW_px*_lepW_px + _lepW_py*_lepW_py)");

        // Probe gen-truth
        if (!isData) {
            df = df.Define("_wsf_truth", wsf_truth_vars,
                           {"fj0_eta", "fj0_phi", "GenPart_pdgId", "GenPart_genPartIdxMother",
                            "GenPart_statusFlags", "GenPart_status", "GenPart_eta", "GenPart_phi"})
                   .Define("fj0_gen_match",     "(int)_wsf_truth[0]")
                   .Define("fj0_dr_Wq_max",     "_wsf_truth[1]")
                   .Define("fj0_dr_top_q",      "_wsf_truth[2]")
                   .Define("fj0_dr_W_daus",     "_wsf_truth[3]")
                   .Define("fj0_dr_b_near",     "_wsf_truth[4]")
                   .Define("fj0_Wq_max_pdgId",  "(int)_wsf_truth[5]")
                   .Define("fj0_top_q_pdgId",   "(int)_wsf_truth[6]");
        }
    }

    // 2lepSS
    else if (channel == "2lepSS"){

        df = VBSTagging(df);
        Cutflow::Add(df, "VBS pair candidate found");

        df = TriggerSelections(df,trigger_logic_string_multilep);
        Cutflow::Add(df, "C1: Trigger selection");

        df = lepSFWrapper(df, isData, /*ele_sf_name=*/ "_weight_electron_reco_looseid", /*muo_sf_name=*/ "_weight_muon_looseid_looseiso");

        // Channel orthogonality selection
        //TODO implement a same sign requirement
        df = df.Filter("(nLep_Sel == 2)", "C2: 2lepSS");

    }

    // 2lep_1FJ (currently shared between OF and SF)
    else if (channel == "2lep_1FJ"){

        df = VBSTagging(df);
        Cutflow::Add(df, "VBS pair candidate found");

        df = TriggerSelections(df,trigger_logic_string_multilep);
        Cutflow::Add(df, "C1: Trigger selection");

        df = lepSFWrapper(df, isData, /*ele_sf_name=*/ "_weight_electron_reco_looseid", /*muo_sf_name=*/ "_weight_muon_mediumid_tightiso");

        // Channel orthogonality selection
        df = definePerVariationPassFlags(df, "2lep_1FJ", [](const std::string& sfx){
            const std::string fjCountCol = sfx.empty() ? "nfatjet" : "nfatjet_" + sfx;
            return "(nLep_Sel == 2) && (" + fjCountCol + " == 1)";
        });
        df = df.Filter(orPassExpr(df, "2lep_1FJ"), "C2: 2lep_1FJ");
    }

    // 2lep_2FJ
    else if (channel == "2lep_2FJ"){

        df = VBSTagging(df);
        Cutflow::Add(df, "VBS pair candidate found");

        df = TriggerSelections(df,trigger_logic_string_multilep);
        Cutflow::Add(df, "C1: Trigger selection");

        df = lepSFWrapper(df, isData, /*ele_sf_name=*/ "_weight_electron_reco_looseid", /*muo_sf_name=*/ "_weight_muon_looseid_looseiso");

        // Channel orthogonality selection
        df = definePerVariationPassFlags(df, "2lep_2FJ", [](const std::string& sfx){
            const std::string fjCountCol = sfx.empty() ? "nfatjet" : "nfatjet_" + sfx;
            return "(nLep_Sel == 2) && (" + fjCountCol + " == 2)";
        });
        df = df.Filter(orPassExpr(df, "2lep_2FJ"), "C2: 2lep_2FJ");
    }


    // 3lep
    else if (channel == "3lep"){

        df = VBSTagging(df);
        Cutflow::Add(df, "VBS pair candidate found");

        df = TriggerSelections(df,trigger_logic_string_multilep);
        Cutflow::Add(df, "C1: Trigger selection");

        df = lepSFWrapper(df, isData, /*ele_sf_name=*/ "_weight_electron_reco_looseid", /*muo_sf_name=*/ "_weight_muon_mediumid_tightiso");

        // Channel orthogonality selection
        df = df.Filter("(nLep_Sel == 3)", "C2: 3lep");

    }

    // 4lep
    else if (channel == "4lep"){

        df = VBSTagging(df);
        Cutflow::Add(df, "VBS pair candidate found");

        df = TriggerSelections(df,trigger_logic_string_multilep);
        Cutflow::Add(df, "C1: Trigger selection");

        df = lepSFWrapper(df, isData, /*ele_sf_name=*/ "_weight_electron_reco_looseid", /*muo_sf_name=*/ "_weight_muon_looseid_looseiso");

        // Channel orthogonality selection
        df = df.Filter("(nLep_Sel == 4)", "C2: 4lep");

    }

    return df;
}
