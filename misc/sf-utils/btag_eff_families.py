"""Conservative MC sample families used for b-tag efficiency production."""

# Keep this list ordered: specific samples must precede broad process prefixes.
FAMILY_RULES = (
    ("VBS_VVH", ("VBSWWH_", "VBSWZH_", "VBSZZH_")),
    ("VBS_SSWW", ("VBS-SSWW",)),
    ("WWJJ", ("WWJJ",)),
    ("ZZJJ", ("ZZJJ",)),
    ("ggZZ_continuum", ("GluGluToContinto2Z", "GluGlutoContinto2Z")),
    ("ggH", ("GluGluH-",)),
    ("VH", ("GluGluZH", "WplusH", "WminusH", "ZH-")),
    ("ttH", ("TTH-",)),
    ("ttV_rare", ("TTWW", "TTWZ", "TTW-", "TZQB")),
    ("ttbar", ("TTto", "TTLL", "TTLNu")),
    ("tW", ("TWminus", "TbarWplus")),
    ("single_top", ("TBbarQ", "TbarBQ", "TBbarto", "TbarBto")),
    ("triboson", ("WWW-", "WWZ-", "WZZ-", "ZZZ-")),
    ("WZ", ("WZ_", "WZto")),
    ("WW", ("WW_", "WWto")),
    ("ZZ", ("ZZ_", "ZZto")),
    ("DY", ("DYto",)),
    ("W_leptonic", ("WtoLNu",)),
    ("W_hadronic", ("Wto2Q",)),
    ("Z_hadronic", ("Zto2Q",)),
    ("QCD_HT", ("QCD-4Jets",)),
    ("QCD_PT", ("QCD_Bin",)),
)


def sample_family(sample):
    """Return the first matching family; rule order resolves intentional overlaps."""
    for family, needles in FAMILY_RULES:
        if any(needle in sample for needle in needles):
            return family
    raise ValueError(f"No b-tag efficiency family is configured for sample {sample!r}")
