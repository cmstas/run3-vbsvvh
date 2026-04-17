# Wrapper script for running the full analysis via SLURM on HPG

PREFIX="/cmsuf/data/"
OUT_DIR="/blue/avery/$USER/samples/run3-vbsvvh/"
mkdir -p $OUT_DIR

# Signal (all_events channel)
python3 run_rdf.py etc/input_sample_jsons/sig_c2v1p0_c3_1p0/all_events/  -p $PREFIX -o $OUT_DIR -n sig_sm  -a all_events -m slurm -r 2 -f 1
python3 run_rdf.py etc/input_sample_jsons/sig_c2v1p5_c3_1p0/all_events/  -p $PREFIX -o $OUT_DIR -n sig_c2v -a all_events -m slurm -r 2 -f 1
python3 run_rdf.py etc/input_sample_jsons/sig_c2v1p0_c3_10p0/all_events/ -p $PREFIX -o $OUT_DIR -n sig_kl  -a all_events -m slurm -r 2 -f 1

# Bkg + data (per-channel skims)

python3 run_rdf.py etc/input_sample_jsons/bkg/0lep_1FJ etc/input_sample_jsons/data/0lep_1FJ -p $PREFIX -o $OUT_DIR -n r2_0lep_1FJ -a 0lep_1FJ -m slurm -r 2 -f 1
python3 run_rdf.py etc/input_sample_jsons/bkg/0lep_2FJ etc/input_sample_jsons/data/0lep_2FJ -p $PREFIX -o $OUT_DIR -n r2_0lep_2FJ -a 0lep_2FJ -m slurm -r 2 -f 1
