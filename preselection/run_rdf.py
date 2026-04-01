import argparse
import json
import os
from datetime import datetime
import subprocess

# The known analysis channels (key) and which skim they use (value)
ANA_CHANNELS = {
        "0lep_0FJ"     : "0lep_0FJ",
        "0lep_1FJ"     : "0lep_1FJ",
        "0lep_1FJ_met" : "0lep_1FJ",
        "0lep_2FJ"     : "0lep_2FJ",
        "0lep_2FJ_met" : "0lep_2FJ",
        "0lep_3FJ"     : "0lep_3FJ",
        "1lep_1FJ"     : "1lep_1FJ",
        "1lep_2FJ"     : "1lep_1FJ",
        #"2lepSS"       : "2lepSS", # DNE yet
        "2lep_1FJ"     : "2lep_1FJ", # Analysis channel shared between SF and OF
        "2lep_2FJ"     : "2lep_2FJ",
        "3lep"         : "3lep",
        "4lep"         : "4lep",
}

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


# Do a check of the merged json to make sure it does not have more than one kind (if we're in local mode)
# This is because the runAnalysis assumes all inputs are of the same kind
def check_inputs(merged_json_dict,mode):
    if mode == "local":
        kind_lst = []
        for ds in merged_json_dict["samples"].keys():
            kind_lst.append(merged_json_dict["samples"][ds]["metadata"]["kind"])
        if len(kind_lst) != len(set(kind_lst)):
            raise Exception("Error, more than one kind of input is specified, not able to run runAnalysis over multiple kinds")



################### Main function ###################
def main():

    # Set up the command line parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--jsons', nargs='+',    help = 'Input json file(s) containing files and metadata')
    parser.add_argument('-c', '--channels', nargs='+', help = 'Which analysis selection channels to run')
    parser.add_argument('-m', '--mode',                help = 'Which mode to run in (local or condor)', choices=['local','condor'])
    parser.add_argument('-o', '--outpath',             help = 'Output directory', default=".")
    parser.add_argument('-n', '--outname',             help = 'Output name', default="rdf_output")
    parser.add_argument('-r', '--run',                 help = 'Which run (2 or 3)', choices=['2','3'])
    parser.add_argument('-p', '--prefix',              help = 'Prefix to append to the file paths', default="/ceph/cms/")
    parser.add_argument('-j', '--n-cores',             help = 'Number of cores to use for local execution', default=64)
    parser.add_argument('-d', '--dry-run',             help = 'Do not actually execute the run command', action='store_true')
    args = parser.parse_args()

    # Get the list of channels to run over (if we ask for "all", use known analysis channels)
    if args.channels == ["all"]: channels_to_run = ANA_CHANNELS.keys()
    else: channels_to_run = args.channels

    # Run RDF once for each specified channel
    for chan_name in channels_to_run:
        print(f"\nRunning for channel {chan_name}\n")

        # Make sure this is a channel RDF knows about
        if (chan_name not in ANA_CHANNELS.keys()) and (chan_name != "all_events"):
            raise Exception(f"Error unknown channel: {chan_name}")

        # Merge the input jsons
        if args.jsons is not None:
            merged_json_dict = merge_jsons(args.jsons)
        else:
            # If no input jsons specified, look them up based on analysis channel
            # Assume we want signal, data, and bkg
            jsons = [
                "etc/input_sample_jsons/sig_c2v1p0_c3_1p0/all_events/",
                "etc/input_sample_jsons/sig_c2v1p5_c3_1p0/all_events/",
                "etc/input_sample_jsons/sig_c2v1p0_c3_10p0/all_events/",
                f"etc/input_sample_jsons/bkg/{ANA_CHANNELS[chan_name]}",
                f"etc/input_sample_jsons/data/{ANA_CHANNELS[chan_name]}"
            ]
            merged_json_dict = merge_jsons(jsons)

        # Prepend the appropriate prefix to all files in the input json
        if args.prefix is not None:
            apply_prefix(merged_json_dict,args.prefix)

        # Construct the output name, given the channel and output tag
        outname = f"{chan_name}_{args.outname}"

        # Dump the merged content to a json file in merged_jsons/
        if not os.path.exists("merged_jsons"): os.mkdir("merged_jsons")
        timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
        merged_json_name = os.path.join("merged_jsons",f"merged_{outname}_{timestamp}.json")
        print(f"  -> Writing merged file \"{merged_json_name}\"")
        with open(merged_json_name, 'w') as outfile:
            json.dump(merged_json_dict, outfile, indent=4)

        # Make sure we are not passing more than one kind to RDF for one runAnalysis
        check_inputs(merged_json_dict,args.mode)

        # Make an output directory out of outpath/outname
        outdir = os.path.join(args.outpath,outname)
        if not os.path.isdir(outdir): os.mkdir(outdir)
        print(f"  -> RDF output will be located in: {outdir}")

        # Construct the bash run command
        if args.mode == "local":
            command = f"bin/runAnalysis -i {merged_json_name} -o {outdir} -n {args.outname} -a {chan_name} -j {args.n_cores} --run_number {args.run} --progress"
            print(f"  -> Now running command \"{command}\"...\n")
            if not args.dry_run: os.system(command)
        elif args.mode == "condor":
            # TODO how to pass output dir to the condor run script
            command = f"python3 condor/submit.py -c {merged_json_name} -a {chan_name} --run_number {args.run}"
            print(f"  -> Running command \"{command}\"...\n")
            if not args.dry_run: os.system(command)

        print("Done!")

main()

