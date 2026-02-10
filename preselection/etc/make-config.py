from glob import glob
import os
import json
import re
import argparse
import uproot
import numpy as np
import concurrent.futures
from tqdm import tqdm

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

ERA_TO_YEAR = {
    # Data eras
    "Run2022C": "2022Re-recoBCD", "Run2022D": "2022Re-recoBCD",
    "Run2022E": "2022Re-recoE+PromptFG", "Run2022F": "2022Re-recoE+PromptFG", "Run2022G": "2022Re-recoE+PromptFG",
    "Run2023C": "2023PromptC",
    "Run2023D": "2023PromptD",
    "Run2024B": "2024Prompt", "Run2024C": "2024Prompt", "Run2024D": "2024Prompt", "Run2024E": "2024Prompt",
    "Run2024F": "2024Prompt", "Run2024G": "2024Prompt", "Run2024H": "2024Prompt", "Run2024I": "2024Prompt",
}

CHANNEL_DATA_MAP = {
    "2Lep": ["MuonEG", "DoubleMuon", "DoubleEG"],
    "1Lep": ["EGamma", "SingleElectron", "Muon", "SingleMuon"],
    "0Lep": ["JetMET", "JetHT"],
}

MC_ERA_PATTERNS = [
    # (pattern, condition, year) - condition is optional extra check
    (r"UL16.*APV", None, "2016preVFP"),
    (r"UL16", None, "2016postVFP"),
    (r"HIPM.*(Run2016[B-F])", None, "2016preVFP"),
    (r"Run2016[F-H]", lambda s: "HIPM" not in s, "2016postVFP"),
    (r"UL17|UL2017|Run2017", None, "2017"),
    (r"UL18|UL2018|Run2018", None, "2018"),
    (r"22EE", None, "2022Re-recoE+PromptFG"),
    (r"22", None, "2022Re-recoBCD"),
    (r"23BPix", None, "2023PromptD"),
    (r"23", None, "2023PromptC"),
    (r"24", None, "2024Prompt"),
]

class Config:
    def __init__(self, sample_category : str, channel : str, samples: str, xsecs: dict, nthreads: int = 8, n_files: int = -1):
        self.sample_category = sample_category
        self.channel = channel
        self.channel_type = self.extract_channel_type(channel)
        self.samples = sorted(glob(samples))
        self.n_files = n_files
        self.configs = {}  # Dictionary to store configs grouped by era/type
        self.process_samples(xsecs, nthreads, n_files)
        self.write_configs()

    def write_configs(self):
        for group_key, config_data in self.configs.items():
            output_file = f"{self.channel}-{self.sample_category}-{group_key}.json"
            with open(output_file, "w") as f:
                json.dump({"samples": config_data}, f, indent=4)
            print(f"Written config to {output_file}")

    def extract_sample_era(self, sample):
        """Extract era from sample. For data, this is the run (e.g., Run2022C). For MC, same as year."""
        if self.sample_category == "data":
            # Extract run from data sample path
            match = re.search(r"(Run\d{4}[A-Z])", sample)
            if match:
                return match.group(1)
            else:
                raise ValueError(f"Error: could not extract era from data sample {sample}")
        else:
            return self.extract_sample_year(sample)

    def extract_sample_year(self, sample):
        if self.sample_category == "data":
            era = self.extract_sample_era(sample)
            if era in ERA_TO_YEAR:
                return ERA_TO_YEAR[era]
            else:
                raise ValueError(f"Error: year not found for era {era}")
        else:
            for pattern, condition, year in MC_ERA_PATTERNS:
                if re.search(pattern, sample):
                    if condition is None or condition(sample):
                        return year
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
    # we need to come up with a better way to do this
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
            "WW": "Boson",
            "WZ": "Boson",
            "ZZ": "Boson",
            "ZH": "Boson",
            "Wminus": "Boson",
            "Wplus": "Boson",
            "Wto": "WJets",
            "WJets": "WJets", # Run2 naming
            "EWK": "EWK",
            "Zto": "ZJets",
            "QCD": "QCD",
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
            match = re.search(r"Run(\d{4}[A-Z])/(.+)", sample)
            return f"Run{match.group(1)}-{match.group(2)}"
        else:
            return re.search(r"/([^/]+)_TuneCP5", sample).group(1)

    @staticmethod
    def extract_channel_type(channel):
        """Extract channel type (0Lep, 1Lep, 2Lep) from channel name like 1Lep2FJ, 0Lep3FJ, etc."""
        match = re.match(r"(\dLep)", channel)
        if match:
            return match.group(1)
        raise ValueError(f"Could not extract channel type from {channel}")

    def is_valid_data_sample(self, sample_name):
        """Check if data sample is valid for the current channel based on CHANNEL_DATA_MAP."""
        if self.sample_category != "data":
            return True
        allowed_datasets = CHANNEL_DATA_MAP.get(self.channel_type, [])
        return any(dataset in sample_name for dataset in allowed_datasets)

    def process_samples(self, xsecs, nthreads, n_files):
        for sample in tqdm(self.samples, desc=f"Processing {self.sample_category} samples"):
            dataset_name = os.path.basename(sample)
            sample_name = self.get_sample_name(sample)
            
            # Filter data samples based on channel
            if self.sample_category == "data" and not self.is_valid_data_sample(sample):
                tqdm.write(f"    -> Skipping {dataset_name} as it's not valid for channel {self.channel}.")
                continue
            
            try:
                _, xsec = self.get_sample_name_and_xsec(dataset_name, xsecs)
            except Exception as e:
                tqdm.write(f"    -> Skipping {dataset_name} as {e}.")
                continue
            sample_year = self.extract_sample_year(sample)
            sample_era = self.extract_sample_era(sample)
            num_events = 0
            files_path = f"{sample}/*.root"
            if self.sample_category != "data":
                sample_type = self.extract_mc_sample_type(sample_name)
                files = np.random.choice(glob(files_path), size=min(n_files if n_files > 0 else 100, len(glob(files_path))), replace=False)
                with concurrent.futures.ProcessPoolExecutor(max_workers=nthreads) as executor:
                    results = list(executor.map(self._process_file, files))
                num_events = sum(results)
                try:
                    lumi = LUMI_MAP[sample_year]
                except KeyError as e:
                    tqdm.write(f"    -> Luminosity for year {sample_year} not found. Skipping {dataset_name}.")
                    continue
                # Don't group signal samples ?
                if self.sample_category == "sig":
                    group_key = "all"
                else:
                    group_key = sample_type  # Group by sample type for MC background
            else:
                sample_type = self.extract_data_sample_type(sample_name)
                num_events = 1.0
                lumi = 1.0
                group_key = sample_era  # Group by era for data
            
            if group_key not in self.configs:
                self.configs[group_key] = {}
            
            self.configs[group_key].update(
                {
                    f"{sample_name}_{sample_year}": {
                        "trees": ["Events"],
                        "files": list(files) if self.n_files != -1 else [files_path],
                        "metadata": {
                            "category": self.sample_category,
                            "year": sample_year,
                            "era": sample_era,
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
    argparser.add_argument("--n-files", type=int, default=-1, help="number of files to process per sample")
    argparser.add_argument("-n", "--nthreads", type=int, default=8, help="number of threads to use for processing files")
    args = argparser.parse_args()

    if not args.category or not any(cat in args.category for cat in ["bkg", "sig", "data"]):
        raise ValueError("Please provide a valid category")
    
    if args.category == "data" and args.n_files != -1:
        print("Warning: --n-files should not be used for data category. Ignoring the provided value.")
        args.n_files = -1

    with open(args.xsecs, "r") as f_xsecs:
        xsecs = json.load(f_xsecs)

    with open(args.config, "r") as f:
        skim_paths = json.load(f)

    if args.channel not in skim_paths:
        raise ValueError(f"No skim paths found for channel {args.channel}")
    if args.category not in skim_paths[args.channel]:
        raise ValueError(f"No skim paths found for category {args.category}")
    
    config = Config(args.category, args.channel, skim_paths[args.channel][args.category], xsecs, args.nthreads, args.n_files)