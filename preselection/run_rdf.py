import argparse
import json
import os
from datetime import datetime
import subprocess


# Merge the input jsons into one dictionary
def merge_jsons(input_paths_lst):
    out_dict = {"samples": {}}

    # Loop over paths
    for path in input_paths_lst:

        # Loop over files in this path
        for filename in os.listdir(path):

            # Skip non json files
            if not filename.endswith('.json'): continue

            # Append the info from this json into the out dict
            with open(os.path.join(path, filename), 'r') as file:
                data = json.load(file)["samples"]
                for k,v in data.items():
                    if k in out_dict: raise Exception(f"ERROR: key \"{k}\" already exsists in the output dict")
                    out_dict["samples"][k] = v 

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
    parser.add_argument('-p', "--prefix",      help = "Prefix to append to the file paths")
    parser.add_argument('-j', "--n-cores",     help = "Number of cores to use for local execution")
    args = parser.parse_args()

    # Merge the input jsons
    merged_json_dict = merge_jsons(args.jsons)

    # Prepend the appropriate prefix to all files in the input json
    apply_prefix(merged_json_dict,args.prefix)

    # Dump the merged content to a json file
    now = datetime.now()
    timestamp = now.strftime("%Y%m%d%H%M%S")
    merged_json_name = f"merged_{args.output_name}_{timestamp}.json"
    print(f"\nWriting merged file \"{merged_json_name}\"")
    with open(merged_json_name, 'w') as outfile:
        json.dump({"samples": merged_json_dict}, outfile, indent=4)

    # Construct the bash run command
    if args.mode == "local":
        command = f"bin/runAnalysis -i {merged_json_name} -o {args.output_path} -t {args.output_name} -j 8"
        print(f"Running command \"{command}\"")
        #os.system(command)
    elif args.mode == "condor":
        command = f"python condor/submit.py -i {merged_json_name} -o {args.output_path} -t {args.output_name}"
        print(f"Running command \"{command}\"")
        #os.system(command)


main()

