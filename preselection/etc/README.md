# Config Generation 

This directory contains tools for generating JSON configuration files used by the preselection framework.

## Directory Structure

```
etc/
├── make_config.py      # Config generator script
├── paths.json          # Sample path definitions
├── xsecs_13TeV.json    # Cross-sections for Run 2 (13 TeV)
├── xsecs_13p6TeV.json  # Cross-sections for Run 3 (13.6 TeV)
└── README.md           # This file
```

## Quick Start

1. Create a `paths.json` file defining where your channel's samples are located:

```json
{
    "0Lep2FJ_run3": {
        "bkg": "/ceph/cms/store/user/aaarora/skims_v3/bkg/*/*",
        "data": "/ceph/cms/store/user/aaarora/skims_v3/data/*/*",
        "sig": "/ceph/cms/store/user/aaarora/skims_v3/sig/Run3Summer24/*"
    }
}
```

2. Generate a config:

```bash
python3 make_config.py --channel 0Lep2FJ_run3 --category sig
```

This creates `0Lep2FJ_run3-sig.json` ready for use with the preselection framework.

## Usage

```
python3 make_config.py --channel CHANNEL --category CATEGORY [options]
```

### Required Arguments

| Argument | Description |
|----------|-------------|
| `--channel` | Analysis channel with run suffix (e.g., `0Lep2FJ_run3`, `1Lep2FJ_run2`) |
| `--category` | Sample category (e.g., `sig`, `bkg`, `bkg_QCD`, `data`) |

The run number is inferred from the channel name (`_run2` or `_run3` suffix) and used to select the appropriate cross-section file.

### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--paths` | `paths.json` | Path to the paths configuration file |
| `--xsecs` | auto | Cross-sections file (auto-selects based on run) |
| `-n`, `--nthreads` | `8` | Number of threads for nevents computation |
| `-o`, `--output` | auto | Output filename (default: `{channel}-{category}.json`) |



## paths.json Format

The paths file maps channels (with run suffix) and categories to glob patterns:

```json
{
    "CHANNEL_runN": {
        "CATEGORY": "GLOB_PATTERN"
    }
}
```

- **Channel naming**: All channels must include `_run2` or `_run3` suffix
- **Glob patterns**: Standard glob syntax with bash-style brace expansion
  - `*` matches anything
  - `{A,B,C}` expands to match A, B, or C



## Output Format

The generated JSON follows the RDataFrame FromSpec format:

```json
{
    "samples": {
        "SAMPLE_NAME_YEAR": {
            "trees": ["Events"],
            "files": ["/path/to/sample/*.root"],
            "metadata": {
                "category": "sig|bkg|data",
                "year": "2022Re-recoBCD|2022Re-recoE+PromptFG|...",
                "type": "TTbar|QCD|...",
                "xsec": 0.00212,
                "lumi": 7.9804,
                "nevents": 123456.0
            }
        }
    }
}
```

## Cross-Section Files

Cross-sections are defined in JSON files:

- `xsecs_13TeV.json` - Run 2 samples (13 TeV)
- `xsecs_13p6TeV.json` - Run 3 samples (13.6 TeV)

The script matches sample names to cross-sections by prefix. If it does not find an exact match, it throws an error.

Example for editing a xsecs file:

```json
{
    "SAMPLE_PREFIX": xsec_value_in_pb,
    "TTto4Q": 431.5,
    "VBSWWH_OSWW_C2V1p0_13p6TeV_5f_LO": 0.00212
}
```

## Troubleshooting

**No files found for a sample:**
- Check that the glob pattern in `paths.json` is correct
- Verify the files exist at the specified location

**No xsec match:**
- Add the sample's cross-section to the appropriate xsecs file
- The sample name must start with a key in the xsecs file

**Could not extract year:**
- The sample path must contain year identifiers (e.g., `22EE`, `UL17`, `Run2022C`)
- Check the `extract_year()` function for supported patterns

**WARNING: No genEventSumw:**
- The ROOT file may be corrupted or have a non-standard structure
- Corrupt files are listed in a `*_corrupt.txt` file alongside the config
