"""
Extract data-driven GloParT score templates for QCD resampling (Run 2 or Run 3).

Motivation
----------
The GloParT tagger scores (HvsQCD, VvsQCD) are badly mismodelled for QCD MC.
To fix this we resample the QCD-MC scores from *data* templates measured in a
QCD-enriched control region, in bins of fat-jet pT and |eta|.

The two scores are *correlated* per jet (they share the QCD term in the
denominator: a strongly b-like jet has suppressed Xqq/Xcs), so they must not be
drawn independently.  Because the analysis picks the H candidate first
(ArgMax HvsQCD) and the V candidates from the remainder, we keep a *factorised*
template that gets P(H | pT,|eta|) exactly right and draws V conditionally on H:

    P(H, V | pT, |eta|)  ~=  P(H | pT, |eta|)  x  P(V | H, pT)

The P(V | H, pT) piece is eta-integrated (the H-V correlation itself is nearly
eta-independent, and integrating keeps the 2D template well populated).

Control region
--------------
- 0 leptons, exactly 2 *good* (analysis-selected) fat jets  -> orthogonal to the
  signal channel (>= 3 good fat jets).
- Hadronic data with the 0-lepton channel HT trigger:
    * Run 3: JetMET, 2022-2025, HLT_PFHT1050.
    * Run 2: JetHT, 2016 (HLT_PFHT800 || HLT_PFHT900) + 2017-2018 (HLT_PFHT1050),
      summed.  The MET-filter list also differs between the two runs.

Note: the input skim (`..._0Lep3FJ`) has a hard *raw* nFatJet >= 3 cut, so the
2-good-fatjet events selected here always have a 3rd fat jet that fails the good
selection.  The residual bias on per-jet score shapes is expected to be small.

Usage
-----
    python resampling.py [run2|run3] [test]     # default run3; "test" -> few files

Output
------
  resampling_pdfs.root       (run3) / resampling_pdfs_run2.root       (run2)
    * h3_HvsQCD, h3_VvsQCD          : TH3D (score, pT, |eta|) marginal counts
    * pdf_<score>_pt..._eta...      : unit-normalised TH1D per (pT,|eta|) bin
    * joint_HV_pt...                : TH2D (HvsQCD x VvsQCD), eta-integrated, one
                                      per pT bin  -> the H-V correlation template

Sampling recipe (per QCD-MC fat jet with a given pT, |eta|)
-----------------------------------------------------------
  1. locate its (pT, |eta|) bin.
  2. H* = pdf_HvsQCD_<that bin>.GetRandom()             # eta-dependent
  3. in joint_HV_pt<that pT>, take the column at H* (ProjectionY of that x-bin);
     V* = that_projection.GetRandom()                   # V conditional on H
     -> if the H* column is empty (thin tail), fall back to the full ProjectionY
        of joint_HV_pt<that pT>, or to pdf_VvsQCD_<that bin>.
"""

import os
import sys
from glob import glob
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

_HERE = os.path.dirname(os.path.abspath(__file__))
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

SKIM_BASE = "/ceph/cms/store/user/aaarora/VBSVVH_skim_v30"


def groups_for_run(run):
    """Return the list of era groups (files, MET filter, is2016) to sum."""
    if run == "run3":
        base = f"{SKIM_BASE}/Run3_Data_v15_v30_0Lep3FJ"
        files = []
        for era in ("2022", "2023", "2024", "2025"):
            files += glob(f"{base}/JetMET*Run{era}*/*.root")
        return [dict(name="Run3-JetMET", files=sorted(files),
                     met=RUN3_MET, is2016=False, split_by_trigger=False)]
    if run == "run2":
        base = f"{SKIM_BASE}/Run2_Data_v15_v30_0Lep3FJ"
        f2016 = sorted(glob(f"{base}/JetHT_Run2016*/*.root"))
        f1718 = sorted(glob(f"{base}/JetHT_Run2017*/*.root")
                       + glob(f"{base}/JetHT_Run2018*/*.root"))
        return [
            dict(name="2016-JetHT", files=f2016,
                 met=RUN2_MET_COMMON, is2016=True, split_by_trigger=True),
            dict(name="2017_2018-JetHT", files=f1718,
                 met=RUN2_MET_COMMON + " && Flag_ecalBadCalibFilter",
                 is2016=False, split_by_trigger=True),
        ]
    raise ValueError(f"unknown run '{run}'")


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


def book_group(files, met_expr, is2016):
    """Apply MET filters + HT trigger + selection and book the marginal/joint histos."""
    df = r.RDataFrame("Events", files)
    r.RDF.Experimental.AddProgressBar(df)
    df = df.Filter(met_expr, "MET filters")
    df = trigger_selection(df, is2016)
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
            tasks.append(dict(name=f"{g['name']} [{'|'.join(present)}]",
                              files=fl, met=g["met"], is2016=g["is2016"]))
    else:
        tasks.append(dict(name=g["name"], files=g["files"],
                          met=g["met"], is2016=g["is2016"]))

print(f"[resampling] period={RUN}  tasks:")
for t in tasks:
    print(f"  - {t['name']}: {len(t['files'])} files (is2016={t['is2016']})")

h3 = {score: None for score in SCORES}
joint = {ipt: None for ipt in range(1, n_pt + 1)}

for g in tasks:
    if not g["files"]:
        print(f"[resampling] WARNING: task {g['name']} has no files, skipping")
        continue
    h3_ptrs, joint_ptrs, report = book_group(g["files"], g["met"], g["is2016"])
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
