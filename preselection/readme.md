# Vector Boson Scattering Analysis Framework

## Preselection Framework

### Directory Structure

- `src/`: Contains the source code for the preselection framework.
- `include/`: Contains the header files miscellaneous imports.
- `etc/`: Contains configuration scripts and JSON files for defining inputs.
- `bin/`: Directory where the compiled binary will be placed.
- `build/`: Directory where the object files will be placed.

### Prerequisites

- Python 3.x
- ROOT
- Uproot
- Boost library (for GoldenJSON)
- Correctionlib

A script to set up the environment is included in `misc/env.sh`

### Compilation

To compile the preselection framework, navigate to the `preselection` directory and run:

```bash
make
```

This will compile the source files and place the binary in the bin/ directory.

---
### Configuration
The inputs for the preselection framework are defined using JSON files. The make-config.py script in the etc/ directory is used to generate these JSON configuration files, specifics will have to be adapted by the analyzer.

#### Generating Configuration Files
To generate a configuration file, run the make-config.py script with the appropriate category (bkg, sig, or data):

```bash
python3 etc/make-config.py --category <category> --channel <channel>
```

For example, to generate a configuration file for background samples:

```bash
python3 etc/make-config.py --category bkg --channel OneLep2FJ
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
bin/runAnalysis -i etc/OneLep2FJ-bkg.json -a OneLep2FJ -n 64
```

More command line options can be found by running 

```bash
bin/runAnalysis --help
```