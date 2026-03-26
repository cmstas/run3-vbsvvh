# Script to build a giant string of logic for RDF to evaluate

DS_DICT = {
    "2016" : {
        "ds_prio_lst" : ["DoubleMuon", "MuonEG", "DoubleEG", "SingleMuon", "SingleElectron"],
        "ds_trg_dict" : {
            "DoubleMuon" : [
                "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
                "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",
                "Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL",
                "Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",
                "TripleMu_12_10_5",
            ],
            "MuonEG" : [
                "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
                "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                "Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                "Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
                "Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                "Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL",
                "Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ",
                "Mu8_DiEle12_CaloIdL_TrackIdL",
                "DiMu9_Ele9_CaloIdL_TrackIdL",
            ],
            "DoubleEG" : [
                "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
                "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                "Ele16_Ele12_Ele8_CaloIdL_TrackIdL",
            ],
            "SingleMuon" : [
                "IsoMu24",
                "IsoTkMu24",
                "IsoMu22_eta2p1",
                "IsoTkMu22_eta2p1",
                "IsoMu22",
                "IsoTkMu22",
                "IsoMu27",
            ],
            "SingleElectron" : [
                'Ele27_WPTight_Gsf',
                "Ele25_eta2p1_WPTight_Gsf",
                "Ele27_eta2p1_WPLoose_Gsf",
            ],
        },
    },

    "2017" : {
        "ds_prio_lst" : ["DoubleMuon", "MuonEG", "DoubleEG", "SingleMuon", "SingleElectron"],
        "ds_trg_dict" : {
            "DoubleMuon" : [
                "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
                "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
                "TripleMu_12_10_5",
            ],
            "MuonEG" : [
                "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
                "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                "Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                "Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                "Mu8_DiEle12_CaloIdL_TrackIdL",
                "Mu8_DiEle12_CaloIdL_TrackIdL_DZ", # Note: Listed in Andrew's thesis, but not TOP-19-001 AN
                "DiMu9_Ele9_CaloIdL_TrackIdL_DZ",
            ],
            "DoubleEG" : [
                "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
                "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                "Ele16_Ele12_Ele8_CaloIdL_TrackIdL",
            ],
            "SingleMuon" : [
                "IsoMu24",
                "IsoMu27",
            ],
            "SingleElectron" : [
                "Ele32_WPTight_Gsf",
                "Ele35_WPTight_Gsf",
            ],
        },
    },

    "2018" : {
        "ds_prio_lst" : ["DoubleMuon", "MuonEG", "SingleMuon", "EGamma"],
        "ds_trg_dict" : {
            "DoubleMuon" : [
                "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
                "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
                "TripleMu_12_10_5",
            ],
            "MuonEG" : [
                "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
                "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                "Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                "Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                "Mu8_DiEle12_CaloIdL_TrackIdL",
                "Mu8_DiEle12_CaloIdL_TrackIdL_DZ",
                "DiMu9_Ele9_CaloIdL_TrackIdL_DZ",
            ],
            "SingleMuon" : [
                "IsoMu24",
                "IsoMu27",
            ],
            "EGamma" : [
                "Ele32_WPTight_Gsf",
                "Ele35_WPTight_Gsf",
                "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
                "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                "Ele16_Ele12_Ele8_CaloIdL_TrackIdL",
            ],
        },
    },

}

# Takes a lits of triggers, makes a string of "(trg1 == true) || (trg2 == true)..." 
def get_or_of_trgs(trg_lst):
    out_str = ""
    for i,trg in enumerate(trg_lst):
        if i < (len(trg_lst) - 1):
            out_str = out_str + f"({trg} == true) || "
        else:
            out_str = out_str + f"({trg} == true)"
    out_str = f"({out_str})"
    return out_str


# Given a ds name and a priority ordered list of ds names, return all ds with higher priority
#     - Used for finding the list of triggers to check for overlap removal
#     - We will extract the priority ordered list of ds names from DS_DICT given a year
#     - Then e.g., if my_ds is dsC, and ds_prio_lst is [dsA, dsB, dsC, dsD], will return trgs for dsA and dsB
#     - Assumes my_ds only shows up once in the list (should always be true)
def get_higher_priority_ds_trgs(my_ds, year):

    # Get the ds_prio_lst from the DS_DICT
    ds_prio_lst = DS_DICT[year]["ds_prio_lst"]

    # Get the triggers in ds that are higher priority than my_ds
    out_lst = []
    for ds_name in ds_prio_lst:
        if ds_name == my_ds: break
        else: out_lst = out_lst + DS_DICT[year]["ds_trg_dict"][ds_name]

    return out_lst


# Takes a year and a list of ds names, finds the triggers for these from DS_DICT
def get_trg_lst(lst_of_ds, year):
    out_lst = []
    for ds_name in lst_of_ds:
        out_lst.append(DS_DICT[year]["ds_trg_dict"])
    return out_lst



############# Main function #############

def main():

    # The final output string
    out_str = ""

    # Loop over the years
    for year in DS_DICT:

        # Check for self-consistency in the DS_DICT
        ds_names_prio_lst = DS_DICT[year]["ds_prio_lst"]
        ds_names_from_keys = DS_DICT[year]["ds_trg_dict"].keys()
        if len(ds_names_prio_lst) != len(ds_names_from_keys): raise Exception(f"Mismatch length between ds_prio_lst and ds_trg_dict keys in year: {year}")
        if set(ds_names_prio_lst) != set(ds_names_from_keys): raise Exception(f"Mismatch names between ds_prio_lst and ds_trg_dict keys in year: {year}")

        # The string we will build up for this year
        passes_no_overlap = ""

        # Loop over datasets in this year
        for ds_name in DS_DICT[year]["ds_trg_dict"]:

            # Grab the list of triggers for this ds
            trgs_for_this_ds = DS_DICT[year]["ds_trg_dict"][ds_name]

            # Find the list of ds names that are higher priority than this one
            trgs_for_higher_priority_ds = get_higher_priority_ds_trgs(ds_name,year)

            # Build a string of ORs between all of the triggers
            trg_passes   = get_or_of_trgs(trgs_for_this_ds)
            trg_overlaps = get_or_of_trgs(trgs_for_higher_priority_ds)

            # Append this to the string for this year (note short dataset name e.g. "MuonEG" is called shortname in the RDF)
            passs_no_overlap = passes_no_overlap + f"( (((shortname=={ds_name}) || !isData) & {trg_passes}) & !({trg_overlaps} & isData) )"

        out_str = out_str + f"(is{year} && {passs_no_overlap})"

    print(out_str)


main()





