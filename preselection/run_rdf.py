import argparse
import json
import os
from datetime import datetime
import subprocess


# Merge the input jsons into one dictionary
def merge_jsons(input_paths_lst):
    out_dict = {"samples": {}}

    def _append_info_to_dict(the_dict,path_to_json):
        # Modify the dict in place with the info from a given json
        with open(path_to_json, 'r') as file:
            data = json.load(file)["samples"]
            for k,v in data.items():
                if k in the_dict: raise Exception(f"ERROR: key \"{k}\" already exsists in the output dict")
                the_dict["samples"][k] = v 

    # Loop over paths
    for path in input_paths_lst:

        # This isn't a path! It's a json file
        if path.endswith(".json"):
            _append_info_to_dict(out_dict,path)

        # This is a dir
        elif os.path.isdir(path):

            # Loop over files in this path
            for filename in os.listdir(path):

                # Skip any non json files
                if not filename.endswith('.json'): continue

                # Append the info from this json into the out dict
                _append_info_to_dict(out_dict,os.path.join(path,filename))

        # It's not a json or a dir, what are you trying to pass?
        else:
            raise Exception("Unknown input type, please pass json files or a directory with json files in it.")

    return out_dict


# Apply prefixes to samples in the input jsons
def apply_prefix(in_dict,prefix):
    for dataset_name in in_dict["samples"]:
        flst_with_prefixes = []
        for fname in in_dict["samples"][dataset_name]["files"]:
            flst_with_prefixes.append(os.path.join(prefix+fname))
            flst_with_prefixes.sort()
        in_dict["samples"][dataset_name]["files"] = flst_with_prefixes



# Main function
def main():

    # Set up the command line parser
    parser = argparse.ArgumentParser()
    parser.add_argument('jsons', nargs='+',    help = 'Input json file(s) containing files and metadata')
    parser.add_argument('-m', "--mode",        help = "Which mode to run in (local or condor)", choices=["local","condor"])
    parser.add_argument('-o', "--output-path", help = "Output directory")
    parser.add_argument('-n', "--output-name", help = "Output name")
    parser.add_argument('-r', "--run",         help = "Which run (2 or 3)", choices=["2","3"])
    parser.add_argument('-p', "--prefix",      help = "Prefix to append to the file paths", default=None)
    parser.add_argument('-j', "--n-cores",     help = "Number of cores to use for local execution")
    args = parser.parse_args()

    # Merge the input jsons
    merged_json_dict = merge_jsons(args.jsons)

    # Prepend the appropriate prefix to all files in the input json
    if args.prefix is not None:
        apply_prefix(merged_json_dict,args.prefix)

    # Dump the merged content to a json file in merged_jsons/
    if not os.path.exists("merged_jsons"): os.mkdir("merged_jsons")
    timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
    merged_json_name = os.path.join("merged_jsons",f"merged_{args.output_name}_{timestamp}.json")
    print(f"\nWriting merged file \"{merged_json_name}\"")
    with open(merged_json_name, 'w') as outfile:
        json.dump(merged_json_dict, outfile, indent=4)

    # Construct the bash run command
    if args.mode == "local":
        command = f"bin/runAnalysis -i {merged_json_name} -o {args.output_path} -a 1Lep2FJ -j 8 --run_number {args.run}"
        print(f"Running command \"{command}\"...\n")
        os.system(command)
    elif args.mode == "condor":
        command = f"python condor/submit.py -i {merged_json_name} -o {args.output_path} -a {args.output_name}"
        print(f"Running command \"{command}\"...\n")
        #os.system(command)


main()

