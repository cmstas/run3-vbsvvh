# Vector Boson Scattering Analysis Framework

## Preselection Framework

### Directory Structure

- `src/`: Contains the source code for the preselection framework.
- `include/`: Contains the header files miscellaneous imports.
- `etc/`: Contains configuration scripts and JSON files for defining inputs.
- `bin/`: Directory where the compiled binary will be placed.
- `build/`: Directory where the object files will be placed.

### Prerequisites
Included in CMSSW 15_0_4:
- Python 3.x
- ROOT
- Uproot
- Boost library (for GoldenJSON)
- Correctionlib

A script to set up the environment is included in `misc/env.sh`

### Compilation

To set up and compile the preselection framework, navigate to the top level directory and run:

```bash
cd misc/CMSSW_15_0_4/src/
cmsenv
cd ../../../preselection/
make
```

This will compile the source files and place the binary in the bin/ directory.

---
### Configuration
The inputs for the preselection framework are defined using JSON files. The make-config.py script in the etc/ directory is used to generate these JSON configuration files, specifics will have to be adapted by the analyzer.

#### Generating Configuration Files
To generate a configuration file, run the make-config.py script with the appropriate category (bkg, sig, or data). The script expects to be run from within `etc/`.

```bash
python3 make-config.py --category <category> --channel <channel>
```

For example, to generate a configuration file for background samples:

```bash
python3 make-config.py --category bkg --channel 1Lep2FJ
```

This will create a JSON file in the etc/ directory with the necessary configuration for the specified category.

---
### Running the Preselection
To run the preselection, use the compiled binary and provide the input specification and output file:

```bash
bin/runAnalysis -i <input_spec.json> -a <channel> -n <n_threads>
```

For example,

```bash
bin/runAnalysis -i etc/1Lep2FJ-bkg.json -a 1Lep2FJ -n 64
```

More command line options can be found by running 

```bash
bin/runAnalysis --help
```
