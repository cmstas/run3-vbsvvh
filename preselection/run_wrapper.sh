# Wrapper script for running the full analysis

# If running locally the /ceph prefix, otherwise use the redirector
PREFIX="/ceph/cms/"
#PREFIX="root://redirector.t2.ucsd.edu:1095//"

# Run over a single file locally
# TODO add example

# Run at scale over signal
python run_rdf.py etc/input_sample_jsons/sig_sm/all_events --prefix $PREFIX -n r2_sig_sm -m condor -r 2

# Run at scale over the 10 sets of data/bkg skims
python run_rdf.py etc/input_sample_jsons/bkg/0lep_0FJ etc/input_sample_jsons/data/0lep_0FJ --prefix $PREFIX -n r2_0lep_0FJ -m condor -r 2
python run_rdf.py etc/input_sample_jsons/bkg/0lep_1FJ etc/input_sample_jsons/data/0lep_1FJ --prefix $PREFIX -n r2_0lep_1FJ -m condor -r 2
python run_rdf.py etc/input_sample_jsons/bkg/0lep_2FJ etc/input_sample_jsons/data/0lep_2FJ --prefix $PREFIX -n r2_0lep_2FJ -m condor -r 2
python run_rdf.py etc/input_sample_jsons/bkg/0lep_3FJ etc/input_sample_jsons/data/0lep_3FJ --prefix $PREFIX -n r2_0lep_3FJ -m condor -r 2
python run_rdf.py etc/input_sample_jsons/bkg/1lep_1FJ etc/input_sample_jsons/data/1lep_1FJ --prefix $PREFIX -n r2_1lep_1FJ -m condor -r 2
python run_rdf.py etc/input_sample_jsons/bkg/2lep_1FJ etc/input_sample_jsons/data/2lep_1FJ --prefix $PREFIX -n r2_2lep_1FJ -m condor -r 2
python run_rdf.py etc/input_sample_jsons/bkg/2lep_2FJ etc/input_sample_jsons/data/2lep_2FJ --prefix $PREFIX -n r2_2lep_2FJ -m condor -r 2
python run_rdf.py etc/input_sample_jsons/bkg/2lepSS   etc/input_sample_jsons/data/2lepSS   --prefix $PREFIX -n r2_2lepSS   -m condor -r 2 # DNE yet
python run_rdf.py etc/input_sample_jsons/bkg/3lep     etc/input_sample_jsons/data/3lep     --prefix $PREFIX -n r2_3lep     -m condor -r 2
python run_rdf.py etc/input_sample_jsons/bkg/4lep     etc/input_sample_jsons/data/4lep     --prefix $PREFIX -n r2_4lep     -m condor -r 2
