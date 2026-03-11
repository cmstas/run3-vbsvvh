The explanation for the older config making script is `README_for_make_config.md`

The script used for making the current set of jsons is `make_sample_jsons.py` (run via `python make_sample_jsons.py`), and is described here along with the other relevant files:
* `make_sample_jsons.py`: This is the main script. It specifies the paths to the skim sets. 
* `dataset_names_ref.py` This includes the names of the datasets for signal, background, and data. Here we can also store a reference of exactly which sample names to which we should apply the EWK corrections. 
* `dataset_sumw_ref.json`: This file stores the sum of the weights for the datasets. Eventually this will not be needed (since the skimmer will output this information into the skim directory, and we can simply read the values from there). 
* `xsec_ref.py`: This lists the cross sections for the samples. The names specified here are the "short names" we use to refer to each dataset uniquely (and should exactly match the start of the full dataset names), but note it does not specify a year. 

The output of the `make_sample_jsons.py` are json files for each dataset, located in the `input_sample_jsons` directory. 
