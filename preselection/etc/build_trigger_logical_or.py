# Script to build a giant string of logic for RDF to evaluate


DS_DICT_HT = {
    "2016" : {
        "ds_prio_lst" : None,
        "ds_trg_dict" : {
            "ds_name_1" : [
                "trg_name_1",
            ],
        },
    },
    "2017" : {
        "ds_prio_lst" : None,
        "ds_trg_dict" : {
            "ds_name_1" : [
                "trg_name_1",
            ],
        },
    },
    "2018" : {
        "ds_prio_lst" : None,
        "ds_trg_dict" : {
            "ds_name_1" : [
                "trg_name_1",
            ],
        },
    },
}


DS_DICT_MET = {
    "2016" : {
        "ds_prio_lst" : None,
        "ds_trg_dict" : {
            "ds_name_1" : [
                "trg_name_1",
            ],
        },
    },
    "2017" : {
        "ds_prio_lst" : None,
        "ds_trg_dict" : {
            "ds_name_1" : [
                "trg_name_1",
            ],
        },
    },
    "2018" : {
        "ds_prio_lst" : None,
        "ds_trg_dict" : {
            "ds_name_1" : [
                "trg_name_1",
            ],
        },
    },
}

DS_DICT_SINGLELEP = {
    "2016" : {
        "ds_prio_lst" : ["SingleMuon", "SingleElectron"],
        "ds_trg_dict" : {
            "SingleMuon" : [
                "HLT_IsoMu24",
                "HLT_IsoTkMu24",
            ],
            "SingleElectron" : [
                "HLT_Ele27_eta2p1_WPTight_Gsf",
            ],
        },
    },
    "2017" : {
        "ds_prio_lst" : ["SingleMuon", "SingleElectron"],
        "ds_trg_dict" : {
            "SingleMuon" : [
                "HLT_IsoMu27",
            ],
            "SingleElectron" : [
                "HLT_Ele32_WPTight_Gsf_L1DoubleEG"
            ],
        },
    },
    "2018" : {
        "ds_prio_lst" : ["SingleMuon", "EGamma"],
        "ds_trg_dict" : {
            "SingleMuon" : [
                "HLT_IsoMu24",
            ],
            "EGamma" : [
                "HLT_Ele32_WPTight_Gsf",
            ],
        },
    },
    "2024" : {
        "ds_prio_lst" : ["Muon", "EGamma"],
        "ds_trg_dict" : {
            "Muon" : [
                "HLT_IsoMu24",
            ],
            "EGamma" : [
                "HLT_Ele30_WPTight_Gsf",
            ],
        },
    },

}


DS_DICT_MULTILEP = {
    "2016" : {
        "ds_prio_lst" : ["DoubleMuon", "MuonEG", "DoubleEG", "SingleMuon", "SingleElectron"],
        "ds_trg_dict" : {
            "DoubleMuon" : [
                "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
                "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",
                "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL",
                "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",
            ],
            "MuonEG" : [
                "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
                "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL",
                "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ",
            ],
            "DoubleEG" : [
                "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
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
                "HLT_Ele27_WPTight_Gsf",
                "HLT_Ele25_eta2p1_WPTight_Gsf",
                "HLT_Ele27_eta2p1_WPLoose_Gsf",
            ],
        },
    },

    "2017" : {
        "ds_prio_lst" : ["DoubleMuon", "MuonEG", "DoubleEG", "SingleMuon", "SingleElectron"],
        "ds_trg_dict" : {
            "DoubleMuon" : [
                "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",
            ],
            "MuonEG" : [
                "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
            ],
            "DoubleEG" : [
                "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
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
                "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
            ],
            "MuonEG" : [
                "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
            ],
            "SingleMuon" : [
                "HLT_IsoMu24",
                "HLT_IsoMu27",
            ],
            "EGamma" : [
                "HLT_Ele32_WPTight_Gsf",
                "HLT_Ele35_WPTight_Gsf",
                "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"
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
#     - We will extract the priority ordered list of ds names from ds_dict given a year
#     - Then e.g., if my_ds is dsC, and ds_prio_lst is [dsA, dsB, dsC, dsD], will return trgs for dsA and dsB
#     - Assumes my_ds only shows up once in the list (should always be true)
def get_higher_priority_ds_trgs(ds_dict, my_ds, year):

    # Get the ds_prio_lst from the ds_dict
    ds_prio_lst = ds_dict[year]["ds_prio_lst"]

    # Get the triggers in ds that are higher priority than my_ds
    out_lst = []
    for ds_name in ds_prio_lst:
        if ds_name == my_ds: break
        else: out_lst = out_lst + ds_dict[year]["ds_trg_dict"][ds_name]

    return out_lst



# Main wrapper function for making the logical OR string for a given dataset dictionary
def dump_logical_or_string(ds_dict,do_overlap_remoal):

    # The final output string
    out_str = ""

    # Loop over the years
    for i,year in enumerate(ds_dict):

        # Check for self-consistency in the ds_dict for dataset priority
        if do_overlap_remoal:
            ds_names_prio_lst = ds_dict[year]["ds_prio_lst"]
            ds_names_from_keys = ds_dict[year]["ds_trg_dict"].keys()
            if len(ds_names_prio_lst) != len(ds_names_from_keys): raise Exception(f"Mismatch length between ds_prio_lst and ds_trg_dict keys in year: {year}")
            if set(ds_names_prio_lst) != set(ds_names_from_keys): raise Exception(f"Mismatch names between ds_prio_lst and ds_trg_dict keys in year: {year}")

        # The string we will build up for this year
        passes_no_overlap = ""

        # Loop over datasets in this year
        for j, ds_name in enumerate(ds_dict[year]["ds_trg_dict"]):

            # Grab the list of triggers for this ds
            trgs_for_this_ds = ds_dict[year]["ds_trg_dict"][ds_name]

            # Build a string of ORs between all of the triggers that pass
            trg_passes   = get_or_of_trgs(trgs_for_this_ds)

            # Build a string of ORs between all of the triggers that overlap
            if do_overlap_remoal:
                trgs_for_higher_priority_ds = get_higher_priority_ds_trgs(ds_dict, ds_name,year)
                trg_overlaps = get_or_of_trgs(trgs_for_higher_priority_ds)
            else:
                trg_overlaps = None


            # Append this to the string for this year (note short dataset name e.g. "MuonEG" is called shortname in the RDF)
            if trg_overlaps is None:
                # No overlap to remove (either this is the highest priority dataset, or we are not doing trigger overlap removal)
                passes_no_overlap = passes_no_overlap + f"( ((shortname==\\\"{ds_name}\\\") || !isData)  && {trg_passes} )"
            else:
                # We are removing overlap, so build a "passes trigger and !overlap" type of string
                passes_no_overlap = passes_no_overlap + f"( (((shortname==\\\"{ds_name}\\\") || !isData) && {trg_passes}) && !({trg_overlaps} && isData) )"


            # Append and OR if this is not the last one
            if j < (len(ds_dict[year]["ds_trg_dict"]) - 1):
                passes_no_overlap = passes_no_overlap + " || "
            # Otherwise if we're done, wrap the whole thing in parentheses
            else:
                passes_no_overlap = f"({passes_no_overlap})"

        # Append to the final out string
        out_str = out_str + f"(is{year} && {passes_no_overlap})"

        # Append the OR if this is not the last one
        if i < len(ds_dict)-1:
            out_str = out_str + " || "


    # Dump the final output string to the screen
    print("")
    print(out_str)



############# Main function #############

def main():

    dump_logical_or_string(DS_DICT_HT,do_overlap_remoal=False)
    dump_logical_or_string(DS_DICT_MET,do_overlap_remoal=False)
    dump_logical_or_string(DS_DICT_SINGLELEP,do_overlap_remoal=True)
    dump_logical_or_string(DS_DICT_MULTILEP,do_overlap_remoal=True)

main()





