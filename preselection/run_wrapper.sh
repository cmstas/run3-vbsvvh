# Wrapper script for running the full analysis

# Specify a local prefix or xrd prefix to reach the files
#PREFIX="/cmsuf/data/" # Local HPG
#PREFIX="/ceph/cms/" # Local UAF
PREFIX="root://cmsio2.rc.ufl.edu//" # Redirector HPG
#PREFIX="root://redirector.t2.ucsd.edu:1095//" # Redirector UAF

# Define an output dir
#OUT_DIR="/cmsuf/data/store/user/t2/users/$USER/vbs_vvh_rdf/test" # Example for HPG
OUT_DIR="/data/userdata/$USER/vbs_vvh_rdf/test" # Example for UAF
mkdir $OUT_DIR

# Examples of running single datase locally for testing
#python3 run_rdf.py -i etc/input_sample_jsons/sig_c2v1p0_c3_1p0/all_events/2017_VBSWZH_c2v1p0_c3_1p0.json --prefix $PREFIX -o $OUT_DIR -n test_small -c all_events -m local -r 2 -j 16
#python3 run_rdf.py -i etc/input_sample_jsons/bkg/3lep/2018_WZTo1L1Nu2Q_4f_TuneCP5_13TeV.json --prefix $PREFIX -o $OUT_DIR -n test_small -c 3lep -m local -r 2 -j 16
#python3 run_rdf.py -i etc/input_sample_jsons/data/2lep_1FJ/2016postVFP_MuonEG_Run2016H-UL2016_NanoAODv15-v1_NANOAOD.json --prefix $PREFIX -o $OUT_DIR -n test_small -c 2lep_1FJ -m local -r 2 -j 16

# Examples of running at scale locally (appropriate for signal, and appropriate for data/bkg in channels that are small enough)
#python3 run_rdf.py -i etc/input_sample_jsons/sig_c2v1p0_c3_1p0/all_events/ --prefix $PREFIX -o $OUT_DIR -n r2sigSM -c 3lep 2lep_1FJ -m local -r 2 -j 16
#python3 run_rdf.py -i etc/input_sample_jsons/data/2lep_1FJ/                --prefix $PREFIX -o $OUT_DIR -n r2data  -c 2lep_1FJ      -m local -r 2 -j 16


# Run at scale over all channels with condor
python3 run_rdf.py --channels all -p $PREFIX -o $OUT_DIR -m condor -n r2_test -r 2



