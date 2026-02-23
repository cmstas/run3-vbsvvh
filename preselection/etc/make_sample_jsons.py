import os
import uproot
import json

import xsec_ref
import dataset_names_ref

# Which skims serve which analysis channels 
#    0lep_0FJ : 0lep_0FJ_6j,
#    0lep_1FJ : 0lep_1FJ_met, 0lep_1FJ_4j,
#    0lep_2FJ : 0lep_2FJ_met, 0lep_2FJ_2j,
#    0lep_3FJ : 0lep_3FJ,
#    1lep_1FJ : 1lep_1FJ,1lep_2FJ,
#    2lep_1FJ : 2lepOSOF_1FJ, 2lepOSSF_1FJ,
#    2lep_2FJ : 2lepOSSF_2FJ,
#    3lepSS   : 2lepSS
#    3lep     : 3lep,
#    4lep     : 4lep,


# Dictionary of the paths to the skim sets
SKIM_PATH_DICT = {
    "all_events" : {
        ("run2", "sig_sm")  : "/ceph/cms/store/user/mmazza/SignalGeneration/VBSVVH_VBSCuts_13TeV_4f_LO_MG_2_9_18_c2v_1p0_c3_10p0_c2Vc3scan_slc7_amd64_gcc10_CMSSW_12_4_8",
        #("run3", "sig_sm")  : "",
    },
    "0lep_0FJ" : {
        ("run2", "bkg")  : "/ceph/cms/store/user/mdittric/skim/nanoaodv15_bkg_0Lep0FJ_11Feb2026_v4",
        #("run2", "data") : "",
        #("run3", "bkg")  : "",
        #("run3", "data") : "",
    },
    "0lep_1FJ" : {
        ("run2", "bkg")  : "/ceph/cms/store/user/mdittric/skim/nanoaodv15_bkg_0Lep1FJ_11Feb2026_v4",
        #("run2", "data") : "",
        #("run3", "bkg")  : "",
        #("run3", "data") : "",
    },
    "0lep_2FJ" : {
        ("run2", "bkg")  : "/ceph/cms/store/user/mdittric/skim/nanoaodv15_bkg_0Lep2FJ_11Feb2026_v4",
        #("run2", "data") : "",
        #("run3", "bkg")  : "",
        #("run3", "data") : "",
    },
    "0lep_3FJ" : {
        ("run2", "bkg")  : "/ceph/cms/store/user/mdittric/skim/nanoaodv15_bkg_0Lep3FJ_11Feb2026_v4",
        #("run2", "data") : "",
        #("run3", "bkg")  : "",
        #("run3", "data") : "",
    },
    "1lep_1FJ" : {
        ("run2", "bkg")  : "/ceph/cms/store/user/mdittric/skim/nanoaodv15_bkg_1Lep1FJ_11Feb2026_v4",
        #("run2", "data") : "",
        #("run3", "bkg")  : "",
        #("run3", "data") : "",
    },
    "2lep_1FJ" : {
        ("run2", "bkg")  : "/ceph/cms/store/user/mdittric/skim/nanoaodv15_bkg_2Lep1FJ_11Feb2026_v4",
        #("run2", "data") : "",
        #("run3", "bkg")  : "",
        #("run3", "data") : "",
    },
    "2lep_2FJ" : {
        ("run2", "bkg")  : "/ceph/cms/store/user/mdittric/skim/nanoaodv15_bkg_2Lep2FJ_11Feb2026_v4",
        #("run2", "data") : "",
        #("run3", "bkg")  : "",
        #("run3", "data") : "",
    },
    "2lep_SS" : {
        #("run2", "bkg")  : "",
        #("run2", "data") : "",
        #("run3", "bkg")  : "",
        #("run3", "data") : "",
    },
    "3lep"    : {
        ("run2", "bkg")  : "/ceph/cms/store/user/mdittric/skim/nanoaodv15_bkg_3Lep_11Feb2026_v4",
        #("run2", "data") : "",
        #("run3", "bkg")  : "",
        #("run3", "data") : "",
    },
    "4lep"    : {
        ("run2", "bkg")  : "/ceph/cms/store/user/mdittric/skim/nanoaodv15_bkg_4Lep_11Feb2026_v4",
        #("run2", "data") : "",
        #("run3", "bkg")  : "",
        #("run3", "data") : "",
    },
}

LUMI_DICT = {
    "UL16APV": 19.52,
    "UL16": 16.81,
    "UL17": 41.48,
    "UL18": 59.83,
    "2022": 7.9804,
    "2022EE": 26.6717,
    "2023": 17.794,
    "2023BPix": 9.451,
    "2024": 109.08,
}


################# Helper functions #################

# Get the list of root files in a given dir
#     - Input should be a full path to the given dir
#     - Return the list of root files joined with the full path
def get_root_file_lst(path_to_dir):
    out_lst = []
    file_name_lst = os.listdir(path_to_dir)
    for fname in file_name_lst:
        if not fname.endswith(".root"): continue
        fname_fullpath = os.path.join(path_to_dir,fname)
        out_lst.append(fname_fullpath)
    out_lst.sort()
    return out_lst


# Match a sample name to an xsec in the given list
# Return the name of the xsec and the value
def match_xsec(dataset_name,xsec_dict):
    n_matches = 0
    for xsec_name in xsec_dict:
        if dataset_name.startswith(xsec_name):
            dataset_name_short = xsec_name
            dataset_xsec = xsec_dict[xsec_name]
            n_matches = n_matches + 1

    # Throw and error if we did not find a match
    if n_matches == 0:
        raise Exception(f"Failed to find xsec match for sample: \"{dataset_name}\"")
    if n_matches > 1:
        raise Exception(f"Found more than one matching xsec for sample: \"{dataset_name}\"")

    return [dataset_name_short,dataset_xsec]


# Get sum of weights for files in a list
# File list should be full path
def get_sow(list_of_files):
    sumw_tot = 0
    for filepath in list_of_files:
        with uproot.open(filepath) as f:
            sumw = sum(f["Runs"]["genEventSumw"].array())
            sumw_tot = sumw_tot + sumw
    return(sumw_tot)


# Strip everything before the a given string in all paths in a list
# Assumes all paths are the same up to the given string
def strip_prefixes(fullpaths_lst,split_on="store"):
    out_lst =[]
    for fullpath in fullpaths_lst:
        before,after = fullpath.split(split_on)
        after = split_on+after
        out_lst.append(after)
    return [before, out_lst]



################# Main wrapper for producing json #################

# Create the dict to dump to the json for a given sample
#     - dataset_info should be from the dataset ref file
#     - path should be absolute
#     - kind should be e.g. "bkg"
#     - xsec_dict should be from the xsec ref file for this kind
def make_json_dict_for_dataset(dataset_info, path, kind, xsec_dict):
        out_dict = {}

        # Get year for this dataset
        year = dataset_info["year"]

        # Get the lumi for this year
        lumi = LUMI_DICT[year]

        # Get name of this dataset
        dataset_name = dataset_info["dataset_name"]

        ##################################################################
        ### NOTE Remove this when skim formatting is updated next time ###
        tag_tmp = path.split("/")[-1] 
        tag_tmp = tag_tmp[:-1] + "2"
        if kind == "bkg": dataset_name = f"{dataset_name}_{tag_tmp}"
        ##################################################################

        # Get the xsec name and value for this dataset
        dataset_name_short, xsec_val = match_xsec(dataset_name,xsec_dict)

        # Get full paths to all root files for this dataset
        file_fullpath_lst = get_root_file_lst(os.path.join(path,dataset_name))
        #print(file_fullpath_lst)

        # Get the sum of weights for all of the files in this dataset
        #sumw = get_sow(file_fullpath_lst)
        sumw = 0

        # Get rid of the local prefix
        local_prefix, file_fullpath_lst = strip_prefixes(file_fullpath_lst)

        # Fill the out dict
        out_dict["trees"] = "Events"
        out_dict["metadata"] = {
            "category" : kind,
            "year" : year,
            "xsec" : xsec_val,
            "lumi" : lumi,
            "sumw" : sumw,
            "local_prefix" : local_prefix,
        }
        out_dict["files"] = file_fullpath_lst

        # Return what will be the key and value of the full json
        return [dataset_name_short, out_dict]



################# Main function #################

def main():

    # Loop over the skim sets (e.g., 3lep)
    for skim_set_name in SKIM_PATH_DICT:
        print(f"\nSkim set: {skim_set_name}")

        # Loop over the kinds of samples for each skim (e.g., bkg)
        for run_tag,kind in SKIM_PATH_DICT[skim_set_name]:

            # Get the set of datasets and xsecs for this kind of sample
            datasets_lst = dataset_names_ref.datasets[(run_tag,kind)]
            xsec_dict = xsec_ref.xsec_dict[kind]

            # Get the path to this set of skims
            path = SKIM_PATH_DICT[skim_set_name][(run_tag,kind)]

            # Loop over all of the datasets and build up a dict we will dump to json
            dict_for_json = {}
            print(f"\n{run_tag} {kind}: {len(datasets_lst)} total datasets.")
            for i,dataset_info in enumerate(datasets_lst):
                print(f"{i+1}/{len(datasets_lst)}: {dataset_info['dataset_name']}")

                # Get the key and value for the output json
                short_name, info_dict = make_json_dict_for_dataset(dataset_info, path, kind, xsec_dict)

                # Add this into the full dict
                dict_for_json[short_name] = info_dict

            # Dump the dict to an output json
            with open(f"input_{run_tag}_{kind}_{skim_set_name}.json", "w") as fp:
                json.dump(dict_for_json, fp, indent=4)



main()
