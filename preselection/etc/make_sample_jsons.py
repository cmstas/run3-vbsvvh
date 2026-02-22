import os
import uproot
import json

import xsec_ref
import dataset_names


# Dictionary of the skim sets and the analysis channels they serve
# NOTE: 2lepSS has no skim?
SKIM_INFO_DICT = {
    "0lep_0FJ" : {
        "ana_chans": ["0lep_0FJ_6j"],
        "path"     : "/ceph/cms/store/user/mdittric/skim/nanoaodv15_bkg_0Lep0FJ_11Feb2026_v4",
    },
    "0lep_1FJ" : {
        "ana_chans": ["0lep_1FJ_met", "0lep_1FJ_4j"],
        "path"     : "/ceph/cms/store/user/mdittric/skim/nanoaodv15_bkg_0Lep1FJ_11Feb2026_v4",
    },
    "0lep_2FJ" : {
        "ana_chans": ["0lep_2FJ_met", "0lep_2FJ_2j"],
        "path"     : "/ceph/cms/store/user/mdittric/skim/nanoaodv15_bkg_0Lep2FJ_11Feb2026_v4",
    },
    "0lep_3FJ" : {
        "ana_chans": ["0lep_3FJ"],
        "path"     : "/ceph/cms/store/user/mdittric/skim/nanoaodv15_bkg_0Lep3FJ_11Feb2026_v4",
    },
    "1lep_1FJ" : {
        "ana_chans": ["1lep_1FJ","1lep_2FJ"],
        "path"     : "/ceph/cms/store/user/mdittric/skim/nanoaodv15_bkg_1Lep1FJ_11Feb2026_v4",
    },
    "2lep_1FJ" : {
        "ana_chans": ["2lepOSOF_1FJ", "2lepOSSF_1FJ"],
        "path"     : "/ceph/cms/store/user/mdittric/skim/nanoaodv15_bkg_2Lep1FJ_11Feb2026_v4",
    },
    "2lep_2FJ" : {
        "ana_chans": ["2lepOSSF_2FJ"],
        "path"     : "/ceph/cms/store/user/mdittric/skim/nanoaodv15_bkg_2Lep2FJ_11Feb2026_v4",
    },
    "3lep"    : {
        "ana_chans": ["3lep"],
        "path"     : "/ceph/cms/store/user/mdittric/skim/nanoaodv15_bkg_3Lep_11Feb2026_v4",
    },
    "4lep"    : {
        "ana_chans": ["4lep"],
        "path"     : "/ceph/cms/store/user/mdittric/skim/nanoaodv15_bkg_4Lep_11Feb2026_v4",
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



################# Main wrapper for producing json #################

# Create a json for a given sample
#     - Input is the full path to the sample in question
def make_json_for_dataset(path_to_dataset,year,xsec_dict):

        out_dict = {}

        # Get dataset_name from /full/path/to/dataset_name
        dataset_name = path_to_dataset.split("/")[-1]

        # Get full paths to all root files for this dataset
        file_fullpath_lst = get_root_file_lst(path_to_dataset)
        #print(file_fullpath_lst)

        # Get the xsec for this dataset
        dataset_name_short, xsec_val = match_xsec(dataset_name,xsec_dict)

        # Get the lumi
        lumi = LUMI_DICT[year]

        # Get the sum of weights
        sumw = get_sow(file_fullpath_lst)

        # Fill the out dict
        out_dict["trees"] = "Events"
        out_dict["files"] = file_fullpath_lst
        out_dict["metadata"] = {
            "category" : "bkg", # TODO fix
            "year" : year,
            "xsec" : xsec_val,
            "lumi" : lumi,
            "sumw" : sumw,
        }

        return [dataset_name_short, out_dict]



################# Main function #################

def main():

    xsec_dict = xsec_ref.xsec_dict["bkg_run2"]

    cat = "bkg"

    # Loop over the skim sets
    for skim_set_name in SKIM_INFO_DICT:
        print(skim_set_name)
        if skim_set_name == "0lep_0FJ": continue
        path = SKIM_INFO_DICT[skim_set_name]["path"]
        skim_fulpath_dict = {}
        n_datasets = len(dataset_names.dataset_dict_run2_bkg)

        # Get the modified tag at the end of the skim name, remove this when skim formatting is updated
        # NOTE: Remove this when skim formatting is updated in later versions
        tag_tmp = path.split("/")[-1] 
        tag_tmp = tag_tmp[:-1] + "2"

        # Loop over the set of datasets
        for i,dataset_info in enumerate(dataset_names.dataset_dict_run2_bkg):

            dataset_name = dataset_info["dataset_name"]

            # Modify the dataset name with the second tag
            # NOTE Remove this when skim formatting is updated in later versions
            dataset_name_tmp = f"{dataset_name}_{tag_tmp}"

            year = dataset_info["year"]

            # Get the full path to this dataset, and create the json
            path_full = os.path.join(path,dataset_name_tmp)

            short_name, info_dict = make_json_for_dataset(path_full,year,xsec_dict)

            # Add the dict for this dataset into the full dict
            skim_fulpath_dict[short_name] = info_dict

            print(f"{i+1}/{n_datasets}: {dataset_name_tmp}")

        with open(f"{skim_set_name}_run2.json", "w") as fp:
            json.dump(skim_fulpath_dict, fp, indent=4)



main()
