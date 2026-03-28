# Script to build a giant string of logic for RDF to evaluate

DS_DICT = {
    "2016" : {
        "ds_prio_lst" : ["DoubleMuon", "MuonEG", "DoubleEG", "SingleMuon", "SingleElectron"],
        "ds_trg_dict" : {
            "DoubleMuon" : [
                "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
                "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",
                "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL",
                "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",
                "HLT_TripleMu_12_10_5",
            ],
            "MuonEG" : [
                "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
                "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
                "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL",
                "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ",
                "HLT_Mu8_DiEle12_CaloIdL_TrackIdL",
                "HLT_DiMu9_Ele9_CaloIdL_TrackIdL",
            ],
            "DoubleEG" : [
                "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
                "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                "HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL",
            ],
            "SingleMuon" : [
                "HLT_IsoMu24",
                "HLT_IsoTkMu24",
                "HLT_IsoMu22_eta2p1",
                "HLT_IsoTkMu22_eta2p1",
                "HLT_IsoMu22",
                "HLT_IsoTkMu22",
                "HLT_IsoMu27",
            ],
            "SingleElectron" : [
                'HLT_Ele27_WPTight_Gsf',
                "HLT_Ele25_eta2p1_WPTight_Gsf",
                "HLT_Ele27_eta2p1_WPLoose_Gsf",
            ],
        },
    },

    "2017" : {
        "ds_prio_lst" : ["DoubleMuon", "MuonEG", "DoubleEG", "SingleMuon", "SingleElectron"],
        "ds_trg_dict" : {
            "DoubleMuon" : [
                "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
                "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
                "HLT_TripleMu_12_10_5",
            ],
            "MuonEG" : [
                "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
                "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                "HLT_Mu8_DiEle12_CaloIdL_TrackIdL",
                "HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ", # Note: Listed in Andrew's thesis, but not TOP-19-001 AN
                "HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ",
            ],
            "DoubleEG" : [
                "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
                "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                "HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL",
            ],
            "SingleMuon" : [
                "HLT_IsoMu24",
                "HLT_IsoMu27",
            ],
            "SingleElectron" : [
                "HLT_Ele32_WPTight_Gsf",
                "HLT_Ele35_WPTight_Gsf",
            ],
        },
    },

    "2018" : {
        "ds_prio_lst" : ["DoubleMuon", "MuonEG", "SingleMuon", "EGamma"],
        "ds_trg_dict" : {
            "DoubleMuon" : [
                "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
                "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
                "HLT_TripleMu_12_10_5",
            ],
            "MuonEG" : [
                "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
                "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                "HLT_Mu8_DiEle12_CaloIdL_TrackIdL",
                "HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ",
                "HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ",
            ],
            "SingleMuon" : [
                "HLT_IsoMu24",
                "HLT_IsoMu27",
            ],
            "EGamma" : [
                "HLT_Ele32_WPTight_Gsf",
                "HLT_Ele35_WPTight_Gsf",
                "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
                "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                "HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL",
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

    if out_str == "":
        return None
    else:
        return f"({out_str})"


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



############# Main function #############

def main():

    # The final output string
    out_str = ""

    # Loop over the years
    for i,year in enumerate(DS_DICT):

        # Check for self-consistency in the DS_DICT
        ds_names_prio_lst = DS_DICT[year]["ds_prio_lst"]
        ds_names_from_keys = DS_DICT[year]["ds_trg_dict"].keys()
        if len(ds_names_prio_lst) != len(ds_names_from_keys): raise Exception(f"Mismatch length between ds_prio_lst and ds_trg_dict keys in year: {year}")
        if set(ds_names_prio_lst) != set(ds_names_from_keys): raise Exception(f"Mismatch names between ds_prio_lst and ds_trg_dict keys in year: {year}")

        # The string we will build up for this year
        passes_no_overlap = ""

        # Loop over datasets in this year
        for j, ds_name in enumerate(DS_DICT[year]["ds_trg_dict"]):

            # Grab the list of triggers for this ds
            trgs_for_this_ds = DS_DICT[year]["ds_trg_dict"][ds_name]

            # Find the list of ds names that are higher priority than this one
            trgs_for_higher_priority_ds = get_higher_priority_ds_trgs(ds_name,year)

            # Build a string of ORs between all of the triggers
            trg_passes   = get_or_of_trgs(trgs_for_this_ds)
            trg_overlaps = get_or_of_trgs(trgs_for_higher_priority_ds)

            # Append this to the string for this year (note short dataset name e.g. "MuonEG" is called shortname in the RDF)
            if trg_overlaps is not None:
                passes_no_overlap = passes_no_overlap + f"( (((shortname==\\\"{ds_name}\\\") || !isData) && {trg_passes}) && !({trg_overlaps} && isData) )"
            else:
                passes_no_overlap = passes_no_overlap + f"( ((shortname==\\\"{ds_name}\\\") || !isData)  && {trg_passes} )"

            # Append and OR if this is not the last one
            if j < (len(DS_DICT[year]["ds_trg_dict"]) - 1):
                passes_no_overlap = passes_no_overlap + " || "
            # Otherwise if we're done, wrap the whole thing in parentheses
            else:
                passes_no_overlap = f"({passes_no_overlap})"

        # Append to the final out string
        out_str = out_str + f"(is{year} && {passes_no_overlap})"

        # Append the OR if this is not the last one
        if i < len(DS_DICT)-1:
            out_str = out_str + " || "



    print(out_str)


main()





