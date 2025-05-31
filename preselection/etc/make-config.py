from glob import glob
import json
import re
import argparse
import uproot
import concurrent.futures

class Config:
    def __init__(self, sample_category : str, channel : str, samples: str, nthreads: int = 8):
        self.sample_category = sample_category
        self.samples = sorted(glob(samples))
        self.config = {"samples": {}}

        self.process_samples(xsecs, nthreads)
        self.write_config(f"{channel}-{sample_category}.json")

    def write_config(self, output_file):
        with open(output_file, "w") as f:
            json.dump(self.config, f, indent=4)

    @staticmethod
    def extract_sample_year(sample):
        # Data
        if any(run in sample for run in ["Run2022C", "Run2022D"]):
            return "2022Re-recoBCD"
        elif any(run in sample for run in ["Run2022E", "Run2022F", "Run2022G"]):
            return "2022Re-recoE+PromptFG"
        # MC
        elif "22EE" in sample:
            return "2022Re-recoE+PromptFG"
        elif "22" in sample:
            return "2022Re-recoBCD"
        else:
            raise ValueError(f"Error: year not found for {sample}")

    @staticmethod
    def get_xsec_weight(xsecs, sample):
        if sample in xsecs:
            return xsecs[sample]
        else:
            raise ValueError(f"xsec not found for {sample}")

    @staticmethod
    def get_lumi(year):
        lumi = {
            "2022Re-recoBCD": 7.9804,
            "2022Re-recoE+PromptFG": 26.6717
        }
        if year in lumi:
            return lumi[year]
        else:
            raise ValueError(f"lumi not found for {year}")

    @staticmethod
    def extract_mc_sample_type(sample_name):
        sample_type_mapping = {
            "DY": "DY",
            "TTto": "TTbar",
            "Wto": "WJets",
            "VBS": "VBS", 
            "QCD": "QCD",
            "WW": "Boson",
            "WZ": "Boson"
        }
        for key, value in sample_type_mapping.items():
            if key in sample_name:
                return value
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
            try:
                sample_name = self.get_sample_name(sample)
                sample_year = self.extract_sample_year(sample)
                xsec = self.get_xsec_weight(xsecs, sample_name) if self.sample_category != "data" else 1.0
                num_events = 0
                files_path = f"{sample}/*.root"
                if self.sample_category != "data":
                    files = glob(files_path)
                    with concurrent.futures.ProcessPoolExecutor(max_workers=nthreads) as executor:
                        results = list(executor.map(self._process_file, files))
                    num_events = sum(results)
                else:
                    num_events = 1.0
                self.config["samples"].update(
                    {
                        f"{sample_name}_{sample_year}": {
                            "trees": ["Events"],
                            "files": [files_path],
                            "metadata": {
                                "sample_category": self.sample_category,
                                "sample_year": sample_year,
                                "sample_type": self.extract_mc_sample_type(sample_name) if self.sample_category != "data" else "Muon" if "Muon" in sample_name else "Electron",
                                "xsec": xsec,
                                "lumi": self.get_lumi(sample_year) if self.sample_category != "data" else 1.0,
                                "nevents": num_events
                            }
                        }
                    }
                )
            except Exception as e:
                print(f"Error in {sample}: {e}")
    
    @staticmethod
    def _process_file(file):
        with uproot.open(file) as upf:
            return sum(upf["Runs"]["genEventSumw"].array())

if __name__ == "__main__":
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--channel", type=str, help="channel", required=True)
    argparser.add_argument("--category", help="categories: bkg, sig or data", required=True)
    argparser.add_argument("--config", help="paths file", default="paths.json")
    argparser.add_argument("-n, --nthreads", type=int, default=8, help="number of threads to use for processing files")
    args = argparser.parse_args()

    if not args.category or not any(cat in args.category for cat in ["bkg", "sig", "data"]):
        raise ValueError("Please provide a valid category")

    with open("xsecs.json", "r") as f_xsecs:
        xsecs = json.load(f_xsecs)

    with open(args.config, "r") as f:
        skim_paths = json.load(f)

    if args.channel not in skim_paths:
        raise ValueError(f"No skim paths found for channel {args.channel}")
    if args.category not in skim_paths[args.channel]:
        raise ValueError(f"No skim paths found for category {args.category}")
    
    config = Config(args.category, args.channel, skim_paths[args.channel][args.category])
