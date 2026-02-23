# Wrapper script for running the full analysis

# If running locally the /ceph prefix, otherwise use the redirector
PREFIX="/ceph/cms/"
#PREFIX="root://redirector.t2.ucsd.edu:1095"

# Run at scale over signal
python condor/submit.py etc/input_sig_sm_all_events.json --prefix PREFIX --tag sig_sm -r 2

# Run at scale over the 10 sets of data/bkg skims
python condor/submit.py etc/input_r2_bkg_0lep_0FJ.json etc/input_r2_data_0lep_0FJ.json --prefix PREFIX --tag r2_0lep_0FJ -r 2
python condor/submit.py etc/input_r2_bkg_0lep_1FJ.json etc/input_r2_data_0lep_1FJ.json --prefix PREFIX --tag r2_0lep_1FJ -r 2
python condor/submit.py etc/input_r2_bkg_0lep_2FJ.json etc/input_r2_data_0lep_2FJ.json --prefix PREFIX --tag r2_0lep_2FJ -r 2
python condor/submit.py etc/input_r2_bkg_0lep_3FJ.json etc/input_r2_data_0lep_3FJ.json --prefix PREFIX --tag r2_0lep_3FJ -r 2
python condor/submit.py etc/input_r2_bkg_1lep_1FJ.json etc/input_r2_data_1lep_1FJ.json --prefix PREFIX --tag r2_1lep_1FJ -r 2
python condor/submit.py etc/input_r2_bkg_2lep_1FJ.json etc/input_r2_data_2lep_1FJ.json --prefix PREFIX --tag r2_2lep_1FJ -r 2
python condor/submit.py etc/input_r2_bkg_2lep_2FJ.json etc/input_r2_data_2lep_2FJ.json --prefix PREFIX --tag r2_2lep_2FJ -r 2
python condor/submit.py etc/input_r2_bkg_2lepSS.json   etc/input_r2_data_2lepSS.json   --prefix PREFIX --tag r2_2lepSS   -r 2 # DNE yet
python condor/submit.py etc/input_r2_bkg_3lep.json     etc/input_r2_data_3lep.json     --prefix PREFIX --tag r2_3lep     -r 2
python condor/submit.py etc/input_r2_bkg_4lep.json     etc/input_r2_data_4lep.json     --prefix PREFIX --tag r2_4lep     -r 2
