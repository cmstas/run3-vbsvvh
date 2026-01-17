from glob import glob
import os
import json
import re
import argparse
import uproot
import concurrent.futures

LUMI_MAP = {
    "2016preVFP": 19.52,
    "2016postVFP": 16.81,
    "2017": 41.48,
    "2018": 59.83,
    "2022Re-recoBCD": 7.9804,
    "2022Re-recoE+PromptFG": 26.6717,
    "2023PromptC": 18.063,
    "2023PromptD": 9.693,
    "2024Prompt": 109.08
}

class Config:
    def __init__(self, sample_category : str, channel : str, samples: str, xsecs: dict, nthreads: int = 8):
        self.sample_category = sample_category
        self.samples = sorted(glob(samples))
        self.config = {"samples": {}}
        self.process_samples(xsecs, nthreads)
        self.write_config(f"{channel}-{sample_category}.json")

    def write_config(self, output_file):
        with open(output_file, "w") as f:
            json.dump(self.config, f, indent=4)

    def extract_sample_year(self, sample):
        if self.sample_category == "data":
            if any(run in sample for run in ["Run2022C", "Run2022D"]):
                return "2022Re-recoBCD"
            elif any(run in sample for run in ["Run2022E", "Run2022F", "Run2022G"]):
                return "2022Re-recoE+PromptFG"
            elif any(run in sample for run in ["Run2023C"]):
                return "2023PromptC"
            elif any(run in sample for run in ["Run2023D"]):
                return "2023PromptD"
            elif any(run in sample for run in ["Run2024B", "Run2024C", "Run2024D", "Run2024E", "Run2024F", "Run2024G", "Run2024H", "Run2024I"]):
                return "2024Prompt"
        else:
            # 2016
            if "UL16" in sample and "APV" in sample:
                return "2016preVFP"
            elif "UL16" in sample:
                return "2016postVFP"
            elif ("HIPM" in sample and 
                any(run in sample for run in ["Run2016B", "Run2016C", "Run2016D", "Run2016E", "Run2016F"])):
                return "2016preVFP"
            elif (not "HIPM" in sample and 
                any(run in sample for run in ["Run2016F", "Run2016G", "Run2016H"])):
                return "2016postVFP"
            # 2017
            elif any(tag in sample for tag in ["UL17", "UL2017", "Run2017"]):
                return "2017"
            # 2018
            elif any(tag in sample for tag in ["UL18", "UL2018", "Run2018"]):
                return "2018"
            # 2022
            elif "22EE" in sample:
                return "2022Re-recoE+PromptFG"
            elif "22" in sample:
                return "2022Re-recoBCD"
            elif "23BPix" in sample:
                return "2023PromptD"
            elif "23" in sample:
                return "2023PromptC"
            elif "24" in sample:
                return "2024Prompt"
            else:
                raise ValueError(f"Error: year not found for {sample}")

    def get_sample_name_and_xsec(self, dataset_name, xsec_dict):
        if self.sample_category == "data":
            return ("data", 1)
        dataset_name_short = ""
        # Loop through the xsec dict and get the name that matches this dataset
        # Find the xsec name that matches this dataset
        matching_keys = [key for key in xsec_dict if dataset_name.startswith(key)]

        if len(matching_keys) != 1:
            raise ValueError(f"Error: could not find unique matching xsec name for dataset {dataset_name}. Found matches: {matching_keys}")
        
        dataset_name_short = matching_keys[0]
        xsec = xsec_dict[dataset_name_short]

        return (dataset_name_short, xsec)

    @staticmethod
    def extract_mc_sample_type(sample_name):
        sample_type_mapping = {
            "VBSWWH_OS": "VBSWWH_OS",
            "VBSWWH_SS": "VBSWWH_SS",
            "VBSZZH": "VBSZZH",
            "VBSWZH": "VBSWZH", 
            "DY": "DY",
            "TTto": "TTbar",
            "TTTo": "TTbar", # Run2 naming
            "TT": "ttX", # if made it past other two, likely ttX
            "ST_": "SingleTop",
            "Wto": "WJets",
            "WJets": "WJets", # Run2 naming
            "EWK": "EWK",
            "Zto": "ZJets",
            "QCD": "QCD",
            "WW": "Boson",
            "WZ": "Boson",
            "ZZ": "Boson",
            "ZH": "Boson",
            "Wminus": "Boson",
            "Wplus": "Boson",
        }
        for key, value in sample_type_mapping.items():
            if key in sample_name:
                return value
        return "Other"

    @staticmethod
    def extract_data_sample_type(sample_name):
        if "Muon" in sample_name or "SingleMuon" in sample_name:
            return "Muon"
        elif "EGamma" in sample_name or "SingleElectron" in sample_name:
            return "Electron"
        elif "MuonEG" in sample_name or "DoubleMuon" in sample_name or "DoubleEG" in sample_name:
            return "DiLepton"
        elif "JetMET" in sample_name or "JetHT" in sample_name or "MET" in sample_name:
            return "JetMET"
        else:
            return "Other"

    @staticmethod
    def get_sample_name(sample):
        if "data" in sample:
            match = re.search(r"Run(\d{4}[A-Z])/(.+)/([A-Za-z]+)", sample)
            return f"Run{match.group(1)}-{match.group(3)}"
        else:
            return re.search(r"/([^/]+)_TuneCP5", sample).group(1)

    def process_samples(self, xsecs, nthreads):
        for sample in self.samples:
            dataset_name = os.path.basename(sample)
            sample_name = self.get_sample_name(sample)
            try:
                _, xsec = self.get_sample_name_and_xsec(dataset_name, xsecs)
            except Exception as e:
                print(f"    -> Skipping {dataset_name} as {e}.")
                continue
            sample_year = self.extract_sample_year(sample)
            num_events = 0
            files_path = f"{sample}/*.root"
            if self.sample_category != "data":
                sample_type = self.extract_mc_sample_type(sample_name)
                files = glob(files_path)
                with concurrent.futures.ProcessPoolExecutor(max_workers=nthreads) as executor:
                    results = list(executor.map(self._process_file, files))
                num_events = sum(results)
                try:
                    lumi = LUMI_MAP[sample_year]
                except KeyError as e:
                    print(f"    -> Luminosity for year {sample_year} not found. Skipping {dataset_name}.")
                    continue
            else:
                sample_type = self.extract_data_sample_type(sample_name)
                num_events = 1.0
                lumi = 1.0
            self.config["samples"].update(
                {
                    f"{sample_name}_{sample_year}": {
                        "trees": ["Events"],
                        "files": [files_path],
                        "metadata": {
                            "category": self.sample_category,
                            "year": sample_year,
                            "type": sample_type,
                            "xsec": float(xsec),
                            "nevents": num_events,
                            "lumi": lumi,
                        }
                    }
                }
            )
    
    @staticmethod
    def _process_file(file):
        try:
            with uproot.open(file) as upf:
                # check if file is corrupted
                assert upf["Events"]
                return sum(upf["Runs"]["genEventSumw"].array())
        except Exception as e:
            with open("corrupt_files", "a") as f:
                f.write(file + "\n")
                return 0
            
if __name__ == "__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--channel", type=str, help="channel", required=True)
    argparser.add_argument("--category", help="categories: bkg, sig or data", required=True)
    argparser.add_argument("--config", help="paths file", default="paths.json")
    argparser.add_argument("--xsecs", help="xsec file", default="xsecs.json")
    argparser.add_argument("-n", "--nthreads", type=int, default=8, help="number of threads to use for processing files")
    args = argparser.parse_args()

    if not args.category or not any(cat in args.category for cat in ["bkg", "sig", "data"]):
        raise ValueError("Please provide a valid category")

    with open(args.xsecs, "r") as f_xsecs:
        xsecs = json.load(f_xsecs)

    with open(args.config, "r") as f:
        skim_paths = json.load(f)

    if args.channel not in skim_paths:
        raise ValueError(f"No skim paths found for channel {args.channel}")
    if args.category not in skim_paths[args.channel]:
        raise ValueError(f"No skim paths found for category {args.category}")
    
    config = Config(args.category, args.channel, skim_paths[args.channel][args.category], xsecs, args.nthreads)
