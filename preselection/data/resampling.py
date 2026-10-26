import json
import os
import re
import subprocess
import sys
from array import array

import ROOT as r

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

r.EnableImplicitMT(32)

# ---------------------------------------------------------------------------
# CLI: run period + test mode
# ---------------------------------------------------------------------------
_args = [a.lower() for a in sys.argv[1:]]
RUN = "run2" if "run2" in _args else "run3"
TEST = "test" in _args

# ---------------------------------------------------------------------------
# Binning
# ---------------------------------------------------------------------------
# Marginal score axis: full [0, 1], finely binned.
N_SCORE_BINS = 100
SCORE_EDGES = array("d", list(np.linspace(0.0, 1.0, N_SCORE_BINS + 1)))
# 2D joint score axis: coarser so the H-V plane stays populated in thin bins.
N_JOINT_BINS = 50
# pT bins: [250, 300, 350, 500, 750, 1000, +inf)  (last bin is "1000+").
PT_EDGES = array("d", [250.0, 300.0, 350.0, 500.0, 750.0, 1000.0, 1.0e5])
# |eta| bins by detector region: barrel (< 1.479, ECAL EB/EE transition) and
# endcap (1.479 - 2.5, the good-fat-jet acceptance edge).  Filled with |eta|.
ETA_EDGES = array("d", [0.0, 1.479, 2.5])

SCORES = ["HvsQCD", "VvsQCD"]

_HERE = os.path.dirname(os.path.abspath(__file__))      # preselection/data
_PRESEL = os.path.dirname(_HERE)                        # preselection
_SUFFIX = "" if RUN == "run3" else f"_{RUN}"
OUT_ROOT = os.path.join(_HERE, f"resampling_pdfs{_SUFFIX}.root")
OUT_PNG = os.path.join(_HERE, f"resampling_pdfs{_SUFFIX}.png")
OUT_JOINT_PNG = os.path.join(_HERE, f"resampling_joint{_SUFFIX}.png")

n_pt = len(PT_EDGES) - 1
n_eta = len(ETA_EDGES) - 1

# ---------------------------------------------------------------------------
# MET-filter lists (match preselection/src/selections.cpp::METFilters)
# ---------------------------------------------------------------------------
RUN3_MET = (
    "Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && "
    "Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter && "
    "Flag_BadPFMuonDzFilter && Flag_hfNoisyHitsFilter && "
    "Flag_eeBadScFilter && Flag_ecalBadCalibFilter"
)
RUN2_MET_COMMON = (
    "Flag_goodVertices && Flag_globalSuperTightHalo2016Filter && "
    "Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && "
    "Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_BadPFMuonFilter && "
    "Flag_BadPFMuonDzFilter && Flag_eeBadScFilter"
)  # 2017/2018 additionally require Flag_ecalBadCalibFilter

# HT trigger, verbatim from preselection/src/selections.h (trigger_logic_string_ht).
# Applied via trigger_selection(), which -- like the C++ TriggerSelections -- defines
# any HLT branch absent from a file as false (so e.g. a 2016 file that only has
# HLT_PFHT900 reduces to HLT_PFHT900). `is2016` is supplied per era group.
TRIGGER_LOGIC_HT = "(is2016 && HLT_PFHT800||HLT_PFHT900 ) || (!is2016 && HLT_PFHT1050)"
HT_TRIGGER_BRANCHES = ["HLT_PFHT800", "HLT_PFHT900", "HLT_PFHT1050"]

# Hadronic (HT-triggered) primary dataset per run; the other PDs in the spec (MET, Muon,
# EGamma) are on different trigger paths and have no place in this CR.
HADRONIC_PD = {"run3": "JetMET", "run2": "JetHT"}


def spec_path(run):
    """The production spec JSON the preselection itself is run on for this CR's skim."""
    return os.path.join(_PRESEL, "etc", "old_config", "0Lep3FJ", f"0Lep3FJ_{run}-data.json")


def era_met_filters(era):
    """MET-filter list for a JERC era key (matches selections.cpp::METFilters)."""
    if era.startswith("2016"):
        return RUN2_MET_COMMON
    if era in ("2017", "2018"):
        return RUN2_MET_COMMON + " && Flag_ecalBadCalibFilter"
    return RUN3_MET


def groups_for_run(run):
    """Era groups (files, era key, MET filter, is2016) to sum.

    Taken from the production spec JSON rather than a path glob so that the input file
    list is exactly the preselection's, and so that every file carries the `year` key
    that selects its JEC payload (2022Re-recoBCD, 2023PromptD, ...) -- a glob cannot tell
    Run2022D from Run2022E, which sit in different JME era directories."""
    if run not in HADRONIC_PD:
        raise ValueError(f"unknown run '{run}'")
    pd = HADRONIC_PD[run]
    with open(spec_path(run)) as fp:
        spec = json.load(fp)
    by_era = {}
    for name, entry in spec["samples"].items():
        if not name.startswith(pd):
            continue
        by_era.setdefault(entry["metadata"]["year"], []).extend(entry["files"])
    if not by_era:
        raise RuntimeError(f"no {pd} samples found in {spec_path(run)}")
    return [dict(name=f"{era}-{pd}", files=sorted(files), era=era,
                 met=era_met_filters(era), is2016=era.startswith("2016"),
                 # 2016 switched HT menu mid-year, so those files need the
                 # trigger-branch partitioning below; Run 3 JetMET is uniform.
                 split_by_trigger=(run == "run2"))
            for era, files in sorted(by_era.items())]


# ---------------------------------------------------------------------------
# Nominal AK8 jet energy corrections
# ---------------------------------------------------------------------------
# Port of corrections.cpp::applyFatJetEnergyCorrections for data: recover the raw pT via
# FatJet_rawFactor, then re-apply the era's DATA L1L2L3Res compound from the pinned
# fatJet_jerc.json.gz.  Same by-name argument resolution as the C++ (the compound's input
# list is era-dependent -- 2023BPix/2024/2025 also take JetPhi -- so positional args are
# not safe).
JME_BASE = "/cvmfs/cms-griddata.cern.ch/cat/metadata/JME/"
CORRECTIONS_CPP = os.path.join(_PRESEL, "src", "corrections.cpp")


def parse_jerc_era_table():
    """era -> (JME era dir, pinned snapshot, JEC tag), read out of corrections.cpp.

    Parsed instead of copied: correctionlib has no "give me the newest tag" API, so the
    snapshot directory and the tag string must be bumped together in that one table.  A
    hand-kept duplicate here would eventually pin a different JEC than the preselection
    applies, which is exactly the mismatch these templates exist to avoid."""
    with open(CORRECTIONS_CPP) as fp:
        src = fp.read()
    try:
        block = src.split("static const std::map<std::string, EraJERC> table = {", 1)[1]
        block = block.split("\n    };", 1)[0]
    except IndexError:
        raise RuntimeError(f"could not locate eraJERCTable() in {CORRECTIONS_CPP}")
    rows = re.findall(
        r'\{"([^"]+)",\s*\{"([^"]+)",\s*"([^"]+)",\s*"([^"]+)",\s*"([^"]+)",\s*"([^"]+)"\}\}',
        block)
    if not rows:
        raise RuntimeError(f"eraJERCTable() in {CORRECTIONS_CPP} did not parse -- layout changed?")
    return {era: (jme_dir, snapshot, jec_tag)
            for era, jme_dir, snapshot, jec_tag, _jer_tag, _year_token in rows}


ERA_JERC = parse_jerc_era_table()


def load_correctionlib():
    """Make correctionlib's C++ API callable from cling.

    The shared library has to be loaded *through ROOT*: importing the correctionlib python
    module instead dlopens it outside cling's JIT search path, and the JIT then fails to
    materialise every symbol in the declaration below."""
    incdir = subprocess.check_output(["correction", "config", "--incdir"], text=True).strip()
    lib = os.path.join(os.path.dirname(incdir), "lib", "libcorrectionlib.so")
    if r.gSystem.Load(lib) < 0:
        raise RuntimeError(f"failed to load {lib} (is the CMSSW environment set up?)")
    r.gInterpreter.AddIncludePath(incdir)
    r.gInterpreter.Declare('#include "correction.h"')


load_correctionlib()

# Function-local statics, not namespace-scope globals: cling fails to run the static
# initialisers of namespace-scope std::map<..., unique_ptr<...>> in a JIT'd module.
r.gInterpreter.Declare("""
namespace jecfj {

using ROOT::VecOps::RVec;

inline std::vector<std::unique_ptr<correction::CorrectionSet>>& ownedSets() {
    static std::vector<std::unique_ptr<correction::CorrectionSet>> v;
    return v;
}
inline std::map<std::string, correction::CompoundCorrection::Ref>& compounds() {
    static std::map<std::string, correction::CompoundCorrection::Ref> m;
    return m;
}

void registerEra(const std::string& era, const std::string& file, const std::string& compound) {
    auto cs = correction::CorrectionSet::from_file(file);
    compounds()[era] = cs->compound().at(compound);
    ownedSets().push_back(std::move(cs));
}

// factor such that FatJet_pt * factor == the preselection's corrected FatJet_pt:
//   (1 - rawFactor) undoes the JEC baked into NanoAOD, then the compound is applied.
RVec<float> jecFactor(const std::string& era,
                      const RVec<float>& pt, const RVec<float>& eta, const RVec<float>& phi,
                      const RVec<float>& area, const RVec<float>& rawFactor,
                      float rho, unsigned int run) {
    RVec<float> factor(pt.size(), 1.0f);
    if (pt.empty()) return factor;
    auto it = compounds().find(era);
    if (it == compounds().end())
        throw std::runtime_error("resampling: no AK8 JEC compound registered for era " + era);
    const auto& comp = *it->second;
    std::vector<correction::Variable::Type> args;
    args.reserve(comp.inputs().size());
    for (size_t i = 0; i < pt.size(); ++i) {
        const float pt_raw = (1.0f - rawFactor[i]) * pt[i];
        args.clear();
        for (const auto& v : comp.inputs()) {
            const std::string n = v.name();
            if      (n == "JetA")   args.push_back((double)area[i]);
            else if (n == "JetEta") args.push_back((double)eta[i]);
            else if (n == "JetPt")  args.push_back((double)pt_raw);
            else if (n == "Rho")    args.push_back((double)rho);
            else if (n == "JetPhi") args.push_back((double)phi[i]);
            else if (n == "run")    args.push_back((double)run);
            else throw std::runtime_error("resampling: unexpected JEC compound input " + n);
        }
        factor[i] = (1.0f - rawFactor[i]) * static_cast<float>(comp.evaluate(args));
    }
    return factor;
}

}
""")

_registered_eras = set()


def register_jec_era(era):
    if era in _registered_eras:
        return
    if era not in ERA_JERC:
        raise RuntimeError(f"era '{era}' from the spec JSON is not in "
                           f"corrections.cpp::eraJERCTable() -- the two are out of sync")
    jme_dir, snapshot, jec_tag = ERA_JERC[era]
    payload = f"{JME_BASE}{jme_dir}/{snapshot}/fatJet_jerc.json.gz"
    compound = f"{jec_tag}_DATA_L1L2L3Res_AK8PFPuppi"
    r.jecfj.registerEra(era, payload, compound)
    _registered_eras.add(era)
    print(f"[resampling] JEC {era}: {compound}  <-  {jme_dir}/{snapshot}")


def partition_by_trigger_branches(files):
    """Group files by which of HT_TRIGGER_BRANCHES they contain, so that within
    each partition the trigger branches are homogeneous (required for
    define-missing-as-false to be valid in a single RDataFrame).

    All files in one dataset directory share a trigger menu, so we probe just one
    file per directory instead of opening every file (much faster over ceph)."""
    import uproot
    from collections import defaultdict

    by_dir = defaultdict(list)
    for f in files:
        by_dir[os.path.dirname(f)].append(f)

    partitions = {}
    dropped = 0
    for d, dfiles in by_dir.items():
        present = None
        for probe in dfiles:  # first readable file in the dir
            try:
                keys = set(uproot.open(probe)["Events"].keys())
            except Exception:
                continue
            present = tuple(b for b in HT_TRIGGER_BRANCHES if b in keys)
            break
        if not present:
            dropped += len(dfiles)
            continue
        partitions.setdefault(present, []).extend(dfiles)
    if dropped:
        print(f"[resampling] WARNING: {dropped} files had no HT trigger branch and were dropped")
    return partitions


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
r.gInterpreter.Declare("""
using namespace ROOT::VecOps;

RVec<float> VVdR(const RVec<float>& vec_eta1, const RVec<float>& vec_phi1, const RVec<float>& vec_eta2, const RVec<float>& vec_phi2) {
    if (vec_eta1.empty()) {
        return RVec<float>();
    }
    if (vec_eta2.empty()) {
        return RVec<float>(vec_eta1.size(), 999.0f);
    }
    RVec<float> out(vec_eta1.size());
    for (size_t i = 0; i < vec_eta1.size(); i++) {
        float mindR = 999.;
        for (size_t j = 0; j < vec_eta2.size(); j++) {
            float dR = ROOT::VecOps::DeltaR(vec_eta1[i], vec_eta2[j], vec_phi1[i], vec_phi2[j]);
            if (dR < mindR) {
                mindR = dR;
            }
        }
        out[i] = mindR;
    }
    return out;
}
""")


def pt_label(lo, hi):
    hi_s = "Inf" if hi >= 1.0e5 else str(int(hi))
    return f"{int(lo)}to{hi_s}"


def eta_label(lo, hi):
    return f"{lo:g}to{hi:g}".replace(".", "p")


def trigger_selection(df, is2016):
    """Clone of preselection/src/selections.cpp::TriggerSelections for the HT path:
    define any HLT branch in the logic string that is absent as false, supply the
    is2016 flag, then filter on the verbatim trigger string. Requires the trigger
    branches to be homogeneous across the RDataFrame's files (see partition_*)."""
    cols = set(str(c) for c in df.GetColumnNames())
    for hlt in HT_TRIGGER_BRANCHES:
        if hlt not in cols:
            df = df.Define(hlt, "false")
    if "is2016" in cols:
        df = df.Redefine("is2016", "true" if is2016 else "false")
    else:
        df = df.Define("is2016", "true" if is2016 else "false")
    return df.Filter(TRIGGER_LOGIC_HT, "HT trigger")


def apply_jec(df, era):
    """Nominal AK8 JEC, redefined in place on FatJet_pt exactly as the preselection does.

    FatJet_mass is not touched (the preselection scales it too, but it never enters this
    script), and neither is FatJet_msoftdrop -- the preselection leaves that uncorrected
    as well, so the msoftdrop > 40 object cut below already matches."""
    register_jec_era(era)
    return (df.Define("_jecEra", f'std::string("{era}")')
              .Define("FatJet_jecFactor",
                      "jecfj::jecFactor(_jecEra, FatJet_pt, FatJet_eta, FatJet_phi, "
                      "FatJet_area, FatJet_rawFactor, Rho_fixedGridRhoFastjetAll, run)")
              .Redefine("FatJet_pt", "FatJet_pt * FatJet_jecFactor"))


def apply_selection(df):
    """0-lepton, exactly-2-good-fatjet CR selection (identical for Run 2 / Run 3)."""
    # Electron selections
    df = (df.Define("Electron_SC_eta", "Electron_eta + Electron_deltaEtaSC")
        .Define("_vetoElectrons",
            "Electron_pt > 10 && "
            "abs(Electron_SC_eta) < 2.5 && "
            "((abs(Electron_SC_eta) <= 1.479 && abs(Electron_dxy) <= 0.05 && abs(Electron_dz) < 0.1) || "
            "((abs(Electron_SC_eta) > 1.479) && abs(Electron_dxy) <= 0.1 && abs(Electron_dz) < 0.2)) && "
            "Electron_cutBased >= 1")
        .Define("_looseElectrons", "_vetoElectrons && Electron_cutBased >= 2")
        .Define("nElectron_Veto", "nElectron == 0 ? 0 : Sum(_vetoElectrons)")
        .Define("nElectron_Loose", "nElectron_Veto == 0 ? 0 : Sum(_looseElectrons)")
        .Define("electron_pt", "Electron_pt[_vetoElectrons]")
        .Define("electron_eta", "Electron_eta[_vetoElectrons]")
        .Define("electron_phi", "Electron_phi[_vetoElectrons]")
        .Define("electron_mass", "Electron_mass[_vetoElectrons]"))

    # Muon selections
    df = (df.Define("_looseMuons",
            "Muon_pt > 10 && "
            "Muon_pfIsoId >= 2 && "
            "abs(Muon_eta) < 2.4 && "
            "abs(Muon_dxy) < 0.2 && "
            "abs(Muon_dz) < 0.5 && "
            "abs(Muon_sip3d) < 8 && "
            "Muon_looseId")
        .Define("nMuon_Loose", "nMuon == 0 ? 0 : Sum(_looseMuons)")
        .Define("muon_pt", "Muon_pt[_looseMuons]")
        .Define("muon_eta", "Muon_eta[_looseMuons]")
        .Define("muon_phi", "Muon_phi[_looseMuons]")
        .Define("muon_mass", "Muon_mass[_looseMuons]"))

    # Combined lepton collection (pt-sorted), used to clean fat jets against leptons
    df = (df.Define("lepton_pt", "Concatenate(electron_pt, muon_pt)")
        .Define("_leptonSorted", "Argsort(-lepton_pt)")
        .Redefine("lepton_pt", "Take(lepton_pt, _leptonSorted)")
        .Define("lepton_eta", "Take(Concatenate(electron_eta, muon_eta), _leptonSorted)")
        .Define("lepton_phi", "Take(Concatenate(electron_phi, muon_phi), _leptonSorted)"))

    # AK8 fat jet selection (matches the analysis "good" fat jet definition)
    df = (df.Define("_dR_ak8_lep", "VVdR(FatJet_eta, FatJet_phi, lepton_eta, lepton_phi)")
        .Define("_good_ak8jets",
            "_dR_ak8_lep > 0.8 && "
            "FatJet_pt > 250 && "
            "abs(FatJet_eta) <= 2.5 && "
            "FatJet_msoftdrop > 40 && "
            "FatJet_jetId > 0")
        .Define("nfatjet", "Sum(_good_ak8jets)")
        .Define("FatJet_HvsQCD", "FatJet_globalParT3_Xbb / (FatJet_globalParT3_Xbb + FatJet_globalParT3_QCD)")
        .Define("FatJet_VvsQCD", "(FatJet_globalParT3_Xqq/3 + FatJet_globalParT3_Xcs) / (FatJet_globalParT3_Xqq/3 + FatJet_globalParT3_Xcs + FatJet_globalParT3_QCD)")
        .Define("fatjet_pt", "FatJet_pt[_good_ak8jets]")
        .Define("fatjet_eta", "FatJet_eta[_good_ak8jets]")
        .Define("fatjet_abseta", "abs(fatjet_eta)")
        .Define("fatjet_HvsQCD", "FatJet_HvsQCD[_good_ak8jets]")
        .Define("fatjet_VvsQCD", "FatJet_VvsQCD[_good_ak8jets]"))

    # 0-lepton, exactly-2-good-fatjet control region (orthogonal to the >=3 FJ channel)
    df = df.Filter("nMuon_Loose == 0 && nElectron_Loose == 0", "0-lepton veto")
    df = df.Filter("nfatjet == 2", "exactly 2 good fat jets")
    return df


def book_group(files, met_expr, is2016, era):
    """Apply MET filters + HT trigger + JEC + selection and book the marginal/joint histos."""
    df = r.RDataFrame("Events", files)
    r.RDF.Experimental.AddProgressBar(df)
    df = df.Filter(met_expr, "MET filters")
    df = trigger_selection(df, is2016)
    # Between the event filters and the jet selection -- the preselection corrects before
    # both, but this ordering matters here: a few skim events carry a FatJet_jetId vector
    # shorter than nFatJet (e.g. run 391572 / event 140368253 in 2025B has 44 fat jets up
    # to 27 TeV and an empty jetId), which makes the RVec && in _good_ak8jets throw on a
    # size mismatch.  Those are noise events that the MET filters reject, so building the
    # jet columns only downstream of the filters keeps them out of the event loop.
    df = apply_jec(df, era)
    df = apply_selection(df)

    # Marginal TH3s (score, pT, |eta|); RVec columns -> one fill per fat jet.
    h3_ptrs = {}
    for score in SCORES:
        tmpl = r.TH3D(
            f"h3_{score}",
            f"{score};GloParT {score};fat jet p_{{T}} [GeV];fat jet |#eta|",
            N_SCORE_BINS, SCORE_EDGES, n_pt, PT_EDGES, n_eta, ETA_EDGES,
        )
        h3_ptrs[score] = df.Histo3D(r.RDF.TH3DModel(tmpl), f"fatjet_{score}", "fatjet_pt", "fatjet_abseta")

    # Eta-integrated joint (HvsQCD x VvsQCD) per pT bin -> H-V correlation.
    joint_ptrs = {}
    for ipt in range(1, n_pt + 1):
        lo_pt, hi_pt = PT_EDGES[ipt - 1], PT_EDGES[ipt]
        tag = f"p{ipt}"
        df = (df.Define(f"_mask_{tag}", f"fatjet_pt >= {lo_pt} && fatjet_pt < {hi_pt} && fatjet_abseta < {ETA_EDGES[n_eta]}")
                .Define(f"_jh_{tag}", f"fatjet_HvsQCD[_mask_{tag}]")
                .Define(f"_jv_{tag}", f"fatjet_VvsQCD[_mask_{tag}]"))
        lab = pt_label(lo_pt, hi_pt)
        model = r.RDF.TH2DModel(
            f"joint_HV_pt{lab}", f"pT {lab};GloParT HvsQCD;GloParT VvsQCD",
            N_JOINT_BINS, 0.0, 1.0, N_JOINT_BINS, 0.0, 1.0,
        )
        joint_ptrs[ipt] = df.Histo2D(model, f"_jh_{tag}", f"_jv_{tag}")

    return h3_ptrs, joint_ptrs, df.Report()


# ---------------------------------------------------------------------------
# Run over all groups and sum the histograms
# ---------------------------------------------------------------------------
groups = groups_for_run(RUN)
if TEST:
    for g in groups:
        g["files"] = g["files"][:8]

# Expand each era group into trigger-branch-homogeneous tasks so the canonical
# HT trigger string can be applied uniformly (missing triggers -> false).
tasks = []
for g in groups:
    if g["split_by_trigger"]:
        for present, fl in sorted(partition_by_trigger_branches(g["files"]).items()):
            tasks.append(dict(name=f"{g['name']} [{'|'.join(present)}]", files=fl,
                              era=g["era"], met=g["met"], is2016=g["is2016"]))
    else:
        tasks.append(dict(name=g["name"], files=g["files"],
                          era=g["era"], met=g["met"], is2016=g["is2016"]))

print(f"[resampling] period={RUN}  tasks:")
for t in tasks:
    print(f"  - {t['name']}: {len(t['files'])} files (era={t['era']}, is2016={t['is2016']})")

h3 = {score: None for score in SCORES}
joint = {ipt: None for ipt in range(1, n_pt + 1)}

for g in tasks:
    if not g["files"]:
        print(f"[resampling] WARNING: task {g['name']} has no files, skipping")
        continue
    h3_ptrs, joint_ptrs, report = book_group(g["files"], g["met"], g["is2016"], g["era"])
    # Trigger this group's event loop and accumulate.
    for score in SCORES:
        h = h3_ptrs[score].GetValue()
        if h3[score] is None:
            h3[score] = h.Clone(f"h3_{score}")
            h3[score].SetDirectory(0)
        else:
            h3[score].Add(h)
    for ipt in range(1, n_pt + 1):
        hj = joint_ptrs[ipt].GetValue()
        if joint[ipt] is None:
            joint[ipt] = hj.Clone(hj.GetName())
            joint[ipt].SetDirectory(0)
        else:
            joint[ipt].Add(hj)
    print(f"\n[resampling] cutflow for group {g['name']}:")
    report.Print()

# ---------------------------------------------------------------------------
# Write everything
# ---------------------------------------------------------------------------
fout = r.TFile(OUT_ROOT, "RECREATE")
pdfs = {score: {} for score in SCORES}  # (ipt, ieta) -> normalised TH1D (for plotting)

for score in SCORES:
    h = h3[score]
    h.Write()
    for ipt in range(1, n_pt + 1):
        for ieta in range(1, n_eta + 1):
            lab = f"{score}_pt{pt_label(PT_EDGES[ipt-1], PT_EDGES[ipt])}_eta{eta_label(ETA_EDGES[ieta-1], ETA_EDGES[ieta])}"
            proj = h.ProjectionX(f"pdf_{lab}", ipt, ipt, ieta, ieta)
            integral = proj.Integral()
            if integral > 0:
                proj.Scale(1.0 / integral)
            proj.SetTitle(lab)
            proj.Write()
            pdfs[score][(ipt, ieta)] = proj

for ipt in range(1, n_pt + 1):
    joint[ipt].Write()  # raw counts; GetRandom uses the integral

# OUT_ROOT already lives at preselection/data/resampling_pdfs*.root, which is the
# relative path the C++ QCD score resampling (utils.cpp applyQCDScoreResampling) loads
# and the path condor/submit.py packages into the job tarball -- nothing to copy.
print(f"[resampling] wrote {OUT_ROOT}")

# ---------------------------------------------------------------------------
# Summaries: statistics per (pT, |eta|) bin, and H-V correlation per pT bin
# ---------------------------------------------------------------------------
print("\nfat jets per (pT, |eta|) bin [HvsQCD TH3]:")
href = h3["HvsQCD"]
print("  pT \\ |eta|   " + "  ".join(
    f"{eta_label(ETA_EDGES[j-1], ETA_EDGES[j]):>12}" for j in range(1, n_eta + 1)))
for ipt in range(1, n_pt + 1):
    row = f"  {pt_label(PT_EDGES[ipt-1], PT_EDGES[ipt]):>12}"
    for ieta in range(1, n_eta + 1):
        row += f"  {href.ProjectionX('_tmp', ipt, ipt, ieta, ieta).Integral():12.0f}"
    print(row)

print("\nH-V correlation (eta-integrated) per pT bin:")
for ipt in range(1, n_pt + 1):
    hj = joint[ipt]
    print(f"  pT {pt_label(PT_EDGES[ipt-1], PT_EDGES[ipt]):>10}: "
          f"corr={hj.GetCorrelationFactor():+.3f}  N={hj.Integral():.0f}")

# ---------------------------------------------------------------------------
# Plot 1: marginal score PDFs overlaid across pT, one panel per (score, |eta|)
# ---------------------------------------------------------------------------
def th1_edges_contents(h):
    nb = h.GetNbinsX()
    edges = np.array([h.GetBinLowEdge(i) for i in range(1, nb + 2)])
    vals = np.array([h.GetBinContent(i) for i in range(1, nb + 1)])
    return edges, vals

fig, axes = plt.subplots(len(SCORES), n_eta, figsize=(7 * n_eta, 5 * len(SCORES)), squeeze=False)
cmap = plt.cm.viridis(np.linspace(0, 0.9, n_pt))
for i, score in enumerate(SCORES):
    for j in range(n_eta):
        ax = axes[i][j]
        for ipt in range(1, n_pt + 1):
            proj = pdfs[score][(ipt, j + 1)]
            if proj.Integral() <= 0:
                continue
            edges, vals = th1_edges_contents(proj)
            ax.stairs(vals, edges, color=cmap[ipt - 1],
                      label=f"pT {pt_label(PT_EDGES[ipt-1], PT_EDGES[ipt])}")
        ax.set_yscale("log")
        ax.set_xlabel(f"GloParT {score}")
        ax.set_ylabel("a.u. (unit norm)")
        ax.set_title(f"[{RUN}] {score}, |eta| {eta_label(ETA_EDGES[j], ETA_EDGES[j+1])}", fontsize=12)
        ax.legend(fontsize=8, ncol=2)
plt.tight_layout()
plt.savefig(OUT_PNG, dpi=120)
print(f"[resampling] wrote {OUT_PNG}")

# ---------------------------------------------------------------------------
# Plot 2: eta-integrated joint (HvsQCD x VvsQCD) per pT bin
# ---------------------------------------------------------------------------
def th2_array(h):
    nx, ny = h.GetNbinsX(), h.GetNbinsY()
    return np.array([[h.GetBinContent(ix, iy) for ix in range(1, nx + 1)]
                     for iy in range(1, ny + 1)])

ncol = 3
nrow = int(np.ceil(n_pt / ncol))
fig2, axes2 = plt.subplots(nrow, ncol, figsize=(5 * ncol, 4.2 * nrow), squeeze=False)
for ipt in range(1, n_pt + 1):
    ax = axes2[(ipt - 1) // ncol][(ipt - 1) % ncol]
    a = th2_array(joint[ipt])
    a = np.ma.masked_where(a <= 0, a)
    im = ax.imshow(a, origin="lower", extent=[0, 1, 0, 1], aspect="auto",
                   norm=LogNorm(), cmap="viridis")
    ax.set_xlabel("GloParT HvsQCD")
    ax.set_ylabel("GloParT VvsQCD")
    ax.set_title(f"[{RUN}] pT {pt_label(PT_EDGES[ipt-1], PT_EDGES[ipt])} "
                 f"(corr {joint[ipt].GetCorrelationFactor():+.3f})", fontsize=11)
    fig2.colorbar(im, ax=ax, fraction=0.046)
for k in range(n_pt, nrow * ncol):
    axes2[k // ncol][k % ncol].axis("off")
plt.tight_layout()
plt.savefig(OUT_JOINT_PNG, dpi=120)
print(f"[resampling] wrote {OUT_JOINT_PNG}")
