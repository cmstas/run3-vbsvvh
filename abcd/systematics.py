"""The per-event systematic weight variations carried through to the datacards.

The preselection stores each of these as a length-3 branch
``[nominal, up, down]`` and folds *only the nominal element* into the event
weight (see ``preselection/src/weights.cpp`` ``applyMCWeights``):

    weight = xsecweight * genWeight * prod(weight_<syst>[0])

So a varied event weight is ``weight * weight_<syst>[k] / weight_<syst>[0]``.
We store those two ratios rather than the raw elements: they are what the
datacard multiplies by, they are independent of the nominal weight, and they
sit near 1.0 so they compress well.

``weight_btagging_sf_HF``/``_LF`` and ``weight_l1prefiring`` are commented out
of the preselection (both the Define and the weight product), so they are not
in the ntuples and are absent here. ``weight_ewk`` is a plain scalar with no
up/down and is likewise not a systematic.
"""

# Branch name -> (nuisance stem, tag with the era?)
#
# Names are scoped per channel and per scan to match the existing
# CMS_{proc}_{scan}_signal_Region* stat nuisances, so nothing is shared between
# cards: every channel x scan floats its own copy. The era tag on the
# experimental entries is descriptive only, since the proc name already pins
# the era.
SYST_WEIGHTS = {
    "weight_muF":             ("scale_muF",      False),
    "weight_muR":             ("scale_muR",      False),
    "weight_PSISR":           ("ps_isr",         False),
    "weight_PSFSR":           ("ps_fsr",         False),
    "weight_pileup":          ("pileup",         True),
    "weight_muonid":          ("eff_m_id",       True),
    "weight_muonreco":        ("eff_m_reco",     True),
    "weight_muontrigger":     ("eff_m_trigger",  True),
    "weight_electronid":      ("eff_e_id",       True),
    "weight_electronreco":    ("eff_e_reco",     True),
    "weight_electrontrigger": ("eff_e_trigger",  True),
}

# Variations treated as acceptance-only: the inclusive yield change is divided
# out so the nuisance carries only migration into the ABCD regions. The signal
# cross section is the parameter of interest, so its normalization must not be
# constrained here as well.
ACCEPTANCE_ONLY = {"weight_muF", "weight_muR", "weight_PSISR", "weight_PSFSR"}

UP_SUFFIX = "_syst_up"
DN_SUFFIX = "_syst_dn"


def ratio_columns(branch):
    """The (up, down) ratio column names written for one systematic branch."""
    return branch + UP_SUFFIX, branch + DN_SUFFIX


def all_ratio_columns():
    cols = []
    for branch in SYST_WEIGHTS:
        cols.extend(ratio_columns(branch))
    return cols


def era_suffix(proc_name):
    """Map a process name like '0lep_3fj_r3' to its centre-of-mass energy tag."""
    name = str(proc_name).lower()
    if "r2" in name.split("_"):
        return "13TeV"
    if "r3" in name.split("_"):
        return "13p6TeV"
    return None


def nuisance_name(branch, proc_name, scan_name):
    """Combine-facing nuisance name, e.g. CMS_1lep_2fj_r3_Scan3_eff_m_id_13p6TeV.

    Scoped by channel and scan to match the signal stat nuisances alongside it.
    """
    stem, tag_era = SYST_WEIGHTS[branch]
    name = f"CMS_{proc_name}_{scan_name}_{stem}"
    if tag_era:
        suffix = era_suffix(proc_name)
        if suffix:
            name = f"{name}_{suffix}"
    return name
