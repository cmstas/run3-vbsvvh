# Wrapper script for running the full analysis.
# Auto-picks slurm (if sbatch is available) or condor for batch submission.

set -e

# ---- Site paths: uncomment the set that matches your machine ----

# HPG local filesystem
PREFIX="/cmsuf/data/"
OUT_DIR="/blue/avery/$USER/samples/run3-vbsvvh"

# HPG via xrd
#PREFIX="root://cmsio2.rc.ufl.edu//"
#OUT_DIR="/blue/avery/$USER/samples/run3-vbsvvh"

# UAF local filesystem
#PREFIX="/ceph/cms/"
#OUT_DIR="/data/userdata/$USER/vbs_vvh_rdf"

# UAF via xrd
#PREFIX="root://redirector.t2.ucsd.edu:1095//"
#OUT_DIR="/data/userdata/$USER/vbs_vvh_rdf"

mkdir -p "$OUT_DIR"

# ---- Pick batch mode based on what's available on the system ----
if command -v sbatch &>/dev/null; then
    MODE="slurm"
else
    MODE="condor"
fi
echo "Batch mode: $MODE"

# ---- Examples of local runs for testing (uncomment one; skip the batch block below) ----
#python3 run_rdf.py -i etc/input_sample_jsons/run2/sig_c2v1p0_c3_1p0/all_events/2017_VBSWZH_c2v1p0_c3_1p0.json -p $PREFIX -o $OUT_DIR -n test -c 3lep      -m local -r 2 -j 16
#python3 run_rdf.py -i etc/input_sample_jsons/run2/sig_c2v1p0_c3_1p0/all_events/                              -p $PREFIX -o $OUT_DIR -n r2sigSM -c 2lep_1FJ -m local -r 2 -j 16

# ---- Batch submission: signal + per-channel bkg/data, for Run 2 and Run 3 ----
for RUN in 2 3; do
    RUN_BASE="etc/input_sample_jsons/run${RUN}"

    # Signal (all_events channel)
    python3 run_rdf.py -i ${RUN_BASE}/sig_c2v1p0_c3_1p0/all_events/  -p $PREFIX -o $OUT_DIR -n r${RUN}_sig_sm  -c all_events -m $MODE -r $RUN -f 1
    python3 run_rdf.py -i ${RUN_BASE}/sig_c2v1p5_c3_1p0/all_events/  -p $PREFIX -o $OUT_DIR -n r${RUN}_sig_c2v -c all_events -m $MODE -r $RUN -f 1
    python3 run_rdf.py -i ${RUN_BASE}/sig_c2v1p0_c3_10p0/all_events/ -p $PREFIX -o $OUT_DIR -n r${RUN}_sig_kl  -c all_events -m $MODE -r $RUN -f 1

    # Bkg + data (per-channel skims)
    for CHAN in 0lep_1FJ 0lep_2FJ; do
        python3 run_rdf.py -i ${RUN_BASE}/bkg/${CHAN} ${RUN_BASE}/data/${CHAN} \
            -p $PREFIX -o $OUT_DIR -n r${RUN}_${CHAN} -c $CHAN -m $MODE -r $RUN -f 1
    done
done
