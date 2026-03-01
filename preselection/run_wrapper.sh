# Wrapper script for running the full analysis

# If running locally the /ceph prefix, otherwise use the redirector
PREFIX="/ceph/cms/"
#PREFIX="root://redirector.t2.ucsd.edu:1095//"

# Define an output dir, e.g. on UAF could be:
OUT_DIR="/data/userdata/$USER/vbs_vvh_rdf/"
mkdir $OUT_DIR

# Run over a single file locally (for testing sig, bkg, data)
#python run_rdf.py etc/input_sample_jsons/sig_c2v1p0_c3_1p0/all_events/2017_VBSWZH_c2v1p0_c3_1p0.json --prefix $PREFIX -o $OUT_DIR -n test_small -a all_events -m local -r 2 -j 1
#python run_rdf.py etc/input_sample_jsons/bkg/1lep_1FJ/2018_WZTo1L1Nu2Q_4f_TuneCP5_13TeV.json --prefix $PREFIX -o $OUT_DIR -n test_small -a 1lep_1FJ -m local -r 2 -j 1
#python run_rdf.py etc/input_sample_jsons/data/0lep_2FJ/2016postVFP_JetHT_Run2016G-UL2016_NanoAODv15-v1_NANOAOD.json --prefix $PREFIX -o $OUT_DIR -n test_small -a 0lep_2FJ -m local -r 2 -j 1

# Run at scale over signal
python run_rdf.py etc/input_sample_jsons/sig_c2v1p0_c3_1p0/all_events/ -p $PREFIX -o $OUT_DIR -n sigsm -a all_events -m condor -r 2

# Run at scale over the 10 sets of data/bkg skims
python run_rdf.py etc/input_sample_jsons/bkg/0lep_0FJ etc/input_sample_jsons/data/0lep_0FJ -p $PREFIX -o $OUT_DIR -n r2_0lep_0FJ -a 0lep_0FJ -m condor -r 2
python run_rdf.py etc/input_sample_jsons/bkg/0lep_1FJ etc/input_sample_jsons/data/0lep_1FJ -p $PREFIX -o $OUT_DIR -n r2_0lep_1FJ -a 0lep_1FJ -m condor -r 2
python run_rdf.py etc/input_sample_jsons/bkg/0lep_2FJ etc/input_sample_jsons/data/0lep_2FJ -p $PREFIX -o $OUT_DIR -n r2_0lep_2FJ -a 0lep_2FJ -m condor -r 2
python run_rdf.py etc/input_sample_jsons/bkg/0lep_3FJ etc/input_sample_jsons/data/0lep_3FJ -p $PREFIX -o $OUT_DIR -n r2_0lep_3FJ -a 0lep_3FJ -m condor -r 2
python run_rdf.py etc/input_sample_jsons/bkg/1lep_1FJ etc/input_sample_jsons/data/1lep_1FJ -p $PREFIX -o $OUT_DIR -n r2_1lep_1FJ -a 1lep_1FJ -m condor -r 2
python run_rdf.py etc/input_sample_jsons/bkg/2lep_1FJ etc/input_sample_jsons/data/2lep_1FJ -p $PREFIX -o $OUT_DIR -n r2_2lep_1FJ -a 2lep_1FJ -m condor -r 2
python run_rdf.py etc/input_sample_jsons/bkg/2lep_2FJ etc/input_sample_jsons/data/2lep_2FJ -p $PREFIX -o $OUT_DIR -n r2_2lep_2FJ -a 2lep_2FJ -m condor -r 2
python run_rdf.py etc/input_sample_jsons/bkg/2lepSS   etc/input_sample_jsons/data/2lepSS   -p $PREFIX -o $OUT_DIR -n r2_2lepSS   -a 2lepSS   -m condor -r 2 # DNE yet
python run_rdf.py etc/input_sample_jsons/bkg/3lep     etc/input_sample_jsons/data/3lep     -p $PREFIX -o $OUT_DIR -n r2_3lep     -a 3lep     -m condor -r 2
python run_rdf.py etc/input_sample_jsons/bkg/4lep     etc/input_sample_jsons/data/4lep     -p $PREFIX -o $OUT_DIR -n r2_4lep     -a 4lep     -m condor -r 2
