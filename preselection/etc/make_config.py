"""
Generate JSON configuration files in RDataFrame FromSpec format for the run3-vbsvvh preselection framework.

Usage:
    python3 make_config.py --channel 0Lep2FJ_run3 --category sig
    python3 make_config.py --channel 0Lep2FJ_run2 --category bkg -o mybkg.json

"""

from glob import glob
import os
import json
import re
import argparse
from pathlib import Path
from typing import Dict, List, Tuple

import concurrent.futures
import uproot


# ============================================================================
# Constants
# ============================================================================

SCRIPT_DIR = Path(__file__).parent.resolve()

LUMI = {
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

SAMPLE_TYPE_MAPPING = {
    "VBSWWH_OS": "VBSWWH_OS",
    "VBSWWH_SS": "VBSWWH_SS",
    "VBSZZH": "VBSZZH",
    "VBSWZH": "VBSWZH",
    "DY": "DY",
    "TTto": "TTbar",
    "TTTo": "TTbar",  # Run2 naming
    "TT": "ttX",  # if made it past other two, likely ttX
    "ST_": "SingleTop",
    "Wto": "WJets",
    "WJets": "WJets",  # Run2 naming
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


# ============================================================================
# Config Generation
# ============================================================================

class ConfigGenerator:
    """Generate JSON config for the analysis framework."""

    def __init__(
        self,
        category: str,
        channel: str,
        sample_paths: str,
        xsecs: Dict[str, float],
        nthreads: int = 8
    ):
        self.category = category
        self.channel = channel
        self.xsecs = xsecs
        self.nthreads = nthreads
        self.is_data = (category == "data")

        # Discover samples from glob pattern
        self.sample_dirs = self._expand_glob(sample_paths)

        # Build config
        self.config = {"samples": {}}
        self.corrupt_files = []

    # ========================================================================
    # Instance methods that use self.is_data
    # ========================================================================

    def extract_year(self, sample_path: str) -> str:
        """Extract year/era from sample path."""
        # Data patterns
        if self.is_data:
            if any(run in sample_path for run in ["Run2022C", "Run2022D"]):
                return "2022Re-recoBCD"
            elif any(run in sample_path for run in ["Run2022E", "Run2022F", "Run2022G"]):
                return "2022Re-recoE+PromptFG"
            elif any(run in sample_path for run in ["Run2023C"]):
                return "2023PromptC"
            elif any(run in sample_path for run in ["Run2023D"]):
                return "2023PromptD"
            elif any(run in sample_path for run in ["Run2024B", "Run2024C", "Run2024D", "Run2024E", "Run2024F", "Run2024G", "Run2024H", "Run2024I"]):
                return "2024Prompt"

        # MC patterns - Run 2
        if "UL16" in sample_path and "APV" in sample_path:
            return "2016preVFP"
        elif "UL16" in sample_path and "APV" not in sample_path:
            return "2016postVFP"
        elif ("HIPM" in sample_path and
              any(run in sample_path for run in ["Run2016B", "Run2016C", "Run2016D", "Run2016E", "Run2016F"])):
            return "2016preVFP"
        elif ("HIPM" not in sample_path and
              any(run in sample_path for run in ["Run2016F", "Run2016G", "Run2016H"])):
            return "2016postVFP"
        elif any(x in sample_path for x in ["UL17", "UL2017", "Run2017"]):
            return "2017"
        elif any(x in sample_path for x in ["UL18", "UL2018", "Run2018"]):
            return "2018"

        # MC patterns - Run 3
        elif "22EE" in sample_path:
            return "2022Re-recoE+PromptFG"
        elif "22" in sample_path and "Run2" not in sample_path:
            return "2022Re-recoBCD"
        elif "23BPix" in sample_path:
            return "2023PromptD"
        elif "23" in sample_path and "Run2" not in sample_path:
            return "2023PromptC"
        elif "24" in sample_path and "Run2" not in sample_path:
            return "2024Prompt"
        else:
            raise ValueError(f"Could not extract year from: {sample_path}")

    # ========================================================================
    # Static utility methods
    # ========================================================================

    @staticmethod
    def extract_mc_sample_type(sample_name: str) -> str:
        """Extract sample type (DY, TTbar, QCD, etc.) from sample name."""
        for key, value in SAMPLE_TYPE_MAPPING.items():
            if key in sample_name:
                return value
        return "Other"

    @staticmethod
    def extract_data_sample_type(sample_name: str) -> str:
        """Extract data sample type from sample name."""
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

    def get_sample_name(self, sample_path: str) -> str:
        """Extract clean sample name from path."""
        if self.is_data:
            match = re.search(r"Run(\d{4}[A-Z])/(.+)/([A-Za-z]+)", sample_path)
            if not match:
                raise ValueError(f"Could not extract data sample name from: {sample_path}")
            return f"Run{match.group(1)}-{match.group(3)}"
        else:
            match = re.search(r"/([^/]+)_TuneCP5", sample_path)
            if not match:
                raise ValueError(f"Could not extract MC sample name from: {sample_path}")
            return match.group(1)

    def match_xsec(self, dataset_name: str) -> Tuple[str, float]:
        """Match dataset name to cross-section."""
        if self.is_data:
            return ("data", 1.0)

        matches = []
        for xsec_name, xsec_value in self.xsecs.items():
            if dataset_name.startswith(xsec_name):
                matches.append((xsec_name, xsec_value))

        if len(matches) != 1:
            raise ValueError(f"No xsec match for: {dataset_name}")

        return matches[0]

    @staticmethod
    def process_file_nevents(filepath: str) -> float:
        """Get nevents (genEventSumw) from a single ROOT file."""
        try:
            with uproot.open(filepath) as f:
                if "Events" not in f:
                    print(f"    WARNING: No Events tree in {filepath}")
                    return 0.0
                if "Runs" not in f:
                    print(f"    WARNING: No Runs tree in {filepath}")
                    return 0.0
                if "genEventSumw" not in f["Runs"].keys():
                    print(f"    WARNING: No genEventSumw in {filepath}")
                    return 0.0
                return float(sum(f["Runs"]["genEventSumw"].array()))
        except Exception as e:
            print(f"    WARNING: Error reading {filepath}: {e}")
            return 0.0

    @staticmethod
    def compute_nevents(files: List[str], nthreads: int = 8) -> Tuple[float, List[str]]:
        """
        Compute total nevents for a list of files.

        Returns (total_nevents, list_of_corrupt_files)
        """
        total = 0.0
        corrupt_files = []

        # Process files in parallel
        with concurrent.futures.ProcessPoolExecutor(max_workers=nthreads) as executor:
            results = list(executor.map(ConfigGenerator.process_file_nevents, files))

        for filepath, nevents in zip(files, results):
            if nevents == 0.0:
                corrupt_files.append(filepath)
            else:
                total += nevents

        return total, corrupt_files

    # ========================================================================
    # Instance methods
    # ========================================================================

    def _expand_glob(self, pattern: str) -> List[str]:
        """Expand glob pattern, handling brace expansion."""
        patterns = self._brace_expand(pattern)
        results = []
        for p in patterns:
            results.extend(sorted(glob(p)))
        return sorted(set(results))

    def _brace_expand(self, pattern: str) -> List[str]:
        """Expand bash-like brace expressions."""
        match = re.search(r"\{([^}]+)\}", pattern)
        if not match:
            return [pattern]

        options = match.group(1).split(",")
        results = []
        for opt in options:
            expanded = pattern[:match.start()] + opt + pattern[match.end():]
            results.extend(self._brace_expand(expanded))
        return results

    def generate(self) -> Dict:
        """Generate the config dictionary."""
        for sample_dir in self.sample_dirs:
            try:
                self._process_sample(sample_dir)
            except Exception as e:
                print(f"  Skipping {sample_dir}: {e}")

        return self.config

    def _process_sample(self, sample_dir: str):
        """Process a single sample directory."""
        dataset_name = os.path.basename(sample_dir.rstrip('/'))

        # Skip PT-binned QCD samples (use HT-binned instead)
        if dataset_name.startswith("QCD_Bin-PT"):
            print(f"Skipping PT-binned QCD sample: {dataset_name}")
            return

        print(f"Processing: {dataset_name}")

        # Extract metadata
        sample_name = self.get_sample_name(sample_dir)

        try:
            xsec_name, xsec = self.match_xsec(dataset_name)
        except ValueError as e:
            print(f"  -> Skipping {dataset_name}: {e}")
            return

        year = self.extract_year(sample_dir)

        if self.is_data:
            sample_type = self.extract_data_sample_type(sample_name)
        else:
            sample_type = self.extract_mc_sample_type(sample_name)

        # Get files - strip trailing slash to avoid double slashes in pattern
        sample_dir_clean = sample_dir.rstrip('/')
        files_pattern = f"{sample_dir_clean}/*.root"
        files = sorted(glob(files_pattern))

        if not files:
            print(f"  -> No files found")
            return

        # Compute nevents
        if self.is_data:
            nevents = 1.0
        else:
            nevents, corrupt = self.compute_nevents(files, self.nthreads)
            self.corrupt_files.extend(corrupt)
            if corrupt:
                print(f"  -> {len(corrupt)} corrupt files")

        # Get luminosity
        lumi = LUMI.get(year, 1.0) if not self.is_data else 1.0

        # Build sample entry
        # Use glob pattern for files (as expected by RDataFrame FromSpec)
        sample_key = f"{sample_name}_{year}"
        self.config["samples"][sample_key] = {
            "trees": ["Events"],
            "files": [f"{sample_dir_clean}/*.root"],
            "metadata": {
                "category": self.category,
                "year": year,
                "type": sample_type,
                "xsec": float(xsec),
                "lumi": lumi,
                "nevents": nevents
            }
        }
        print(f"  -> {sample_key}: {len(files)} files, {nevents:.0f} events")

    def write(self, output_file: str):
        """Write config to file."""
        with open(output_file, 'w') as f:
            json.dump(self.config, f, indent=4)
        print(f"Wrote config to: {output_file}")

        # Write corrupt files list if any
        if self.corrupt_files:
            corrupt_file = output_file.replace('.json', '_corrupt.txt')
            with open(corrupt_file, 'w') as f:
                f.write('\n'.join(self.corrupt_files))
            print(f"Wrote corrupt files list to: {corrupt_file}")


# ============================================================================
# Main
# ============================================================================

def infer_run_from_channel(channel: str) -> int:
    """Infer run number from channel name (expects _run2 or _run3 suffix)."""
    if "_run2" in channel.lower():
        return 2
    elif "_run3" in channel.lower():
        return 3
    else:
        raise ValueError(f"Cannot infer run from channel '{channel}'. Expected '_run2' or '_run3' suffix.")


def main():
    parser = argparse.ArgumentParser(
        description="Generate config files for run3-vbsvvh preselection",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--channel", required=True,
                        help="Analysis channel (e.g., 0Lep2FJ_run3, 1Lep2FJ_run2)")
    parser.add_argument("--category", required=True,
                        help="Sample category (sig, bkg, data)")
    parser.add_argument("--paths", default="paths.json",
                        help="Paths config file (default: paths.json)")
    parser.add_argument("--xsecs", default=None,
                        help="Cross-sections file (default: auto-select based on run)")
    parser.add_argument("-n", "--nthreads", type=int, default=8,
                        help="Number of threads for nevents computation")
    parser.add_argument("-o", "--output", default=None,
                        help="Output file (default: {channel}-{category}.json)")
    args = parser.parse_args()

    # Validate category
    if not any(cat in args.category for cat in ["bkg", "sig", "data"]):
        print(f"Error: Invalid category '{args.category}'. Must contain 'bkg', 'sig', or 'data'.")
        return 1

    # Infer run from channel name
    try:
        run_number = infer_run_from_channel(args.channel)
    except ValueError as e:
        print(f"Error: {e}")
        return 1

    # Determine xsecs file
    if args.xsecs:
        xsecs_file = args.xsecs
    else:
        xsecs_file = "xsecs_13TeV.json" if run_number == 2 else "xsecs_13p6TeV.json"

    xsecs_path = SCRIPT_DIR / xsecs_file
    if not xsecs_path.exists():
        print(f"Error: xsecs file not found: {xsecs_path}")
        return 1

    with open(xsecs_path) as f:
        xsecs = json.load(f)

    # Load paths config
    paths_path = SCRIPT_DIR / args.paths
    if not paths_path.exists():
        print(f"Error: paths file not found: {paths_path}")
        return 1

    with open(paths_path) as f:
        paths = json.load(f)

    # Find sample paths for this channel/category
    if args.channel not in paths:
        print(f"Error: channel '{args.channel}' not found in paths config")
        print(f"Available channels: {list(paths.keys())}")
        return 1

    if args.category not in paths[args.channel]:
        print(f"Error: category '{args.category}' not found for channel '{args.channel}'")
        print(f"Available categories: {list(paths[args.channel].keys())}")
        return 1

    sample_paths = paths[args.channel][args.category]
    if not sample_paths:
        print(f"Error: empty path for {args.channel}/{args.category}")
        return 1

    # Generate config
    print(f"Channel: {args.channel}")
    print(f"Category: {args.category}")
    print(f"Run: {run_number}")
    print(f"Xsecs: {xsecs_file}")
    print(f"Sample paths: {sample_paths}")
    print()

    generator = ConfigGenerator(
        category=args.category,
        channel=args.channel,
        sample_paths=sample_paths,
        xsecs=xsecs,
        nthreads=args.nthreads,
    )

    generator.generate()

    # Write output (always to etc/ directory)
    output_filename = args.output or f"{args.channel}-{args.category}.json"
    output_file = SCRIPT_DIR / output_filename
    generator.write(str(output_file))

    return 0


if __name__ == "__main__":
    exit(main())
