# Vector Boson Scattering Analysis Framework

## Preselection Framework

### Directory Structure

- `src/`: Contains the source code for the preselection framework.
- `include/`: Contains the header files miscellaneous imports.
- `etc/`: Contains configuration scripts and JSON files for defining inputs.
- `bin/`: Directory where the compiled binary will be placed.
- `build/`: Directory where the object files will be placed.

### Prerequisites

- ROOT
- GCC compiler
- Python 3.x
- Uproot
- Boost library (for GoldenJSON)

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
python3 etc/make-config.py --category <category>
```

For example, to generate a configuration file for background samples:

```bash
python3 etc/make-config.py --category bkg
```

This will create a JSON file in the etc/ directory with the necessary configuration for the specified category.

---
### Running the Preselection
To run the preselection, use the compiled binary and provide the input specification and output file:

```bash
bin/runAnalysis --spec <input_spec> --output <output_file>
```

#### Command Line Arguments:
```
--spec: Path to the JSON configuration file defining the input specification.
--output: Path to the output file where the results will be saved.
--nthreads: Number of threads to use for parallel processing (optional).
```
