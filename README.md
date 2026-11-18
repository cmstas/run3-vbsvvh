# Vector Boson Scattering Analysis Framework

Repository based on RDF for VBS analysis. 
See [preselection](https://github.com/cmstas/run3-vbsvvh/tree/main/preselection#preselection-framework) for preprocessing, and [plotting](https://github.com/cmstas/run3-vbsvvh/tree/main/plotter#plotting-script-documentation) (or [ewkcoffea](https://github.com/cmstas/ewkcoffea)) for histogramming.

## 1. Preselection
Below outlines the preselection chain. User only needs to run single command `run_rdf.sh` once, and everything is made.
1. User runs `run_rdf.sh`.
2. `run_rdf.sh` runs `run_rdf.py` 6 times (sig, bkg, data)$\times$(run2, run3)
    - Each of the 6 function call runs for both 1FJ and 2FJ channels.
    - Input files are found in [`preselection/etc/input_sample_jsons`](preselection/etc/input_sample_jsons).
    - Output are saved in `/groups/cjessop/users/jkil/HVV_2L_RDF`.
3. `run_rdf.py` runs [`condor/submit.py`](preselection/condor/submit.py).
    - `submit.py` script compiles and runs `main.cpp` through `executable.sh`.
    - Details are found [here](preselection/README.md#notes-on-running).
4. `submit.py` runs `src/main.cpp`.

## 2. makeDataset.py
- Three outputs:
    1) `abcd/dataset/c2v*_c3_*/run*_*FJ_train.pt`
    2) `abcd/dataset/c2v*_c3_*/run*_*FJ_valid.pt`
    3) `abcd/dataset/c2v*_c3_*/run*_*FJ.parq`

## 3. 