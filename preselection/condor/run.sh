#!/bin/bash

# epoch time
DATE=$(date +%s)

# ./submit.py -c /home/users/aaarora/phys/run3/cmstas-run3-vbsvvh/preselection/etc/old_config/1Lep2FJ/1Lep2FJ_run2-sig.json -a 1lep_2FJ -r 2 -t r2_1lep_2fj_sig_$DATE -j 4 -m 8 --files-per-job 5 --cutflow
# ./submit.py -c /home/users/aaarora/phys/run3/cmstas-run3-vbsvvh/preselection/etc/old_config/1Lep2FJ/1Lep2FJ_run2-sig.json -a 1lep_1FJ -r 2 -t r2_1lep_1fj_sig_$DATE -j 4 -m 8 --files-per-job 5 --cutflow
# ./submit.py -c /home/users/aaarora/phys/run3/cmstas-run3-vbsvvh/preselection/etc/old_config/0Lep3FJ/0Lep3FJ_run2-sig.json -a 0lep_3FJ -r 2 -t r2_0Lep_3fj_sig_$DATE -j 4 -m 8 --files-per-job 5 --cutflow

./submit.py -c /home/users/aaarora/phys/run3/cmstas-run3-vbsvvh/preselection/etc/old_config/1Lep2FJ/1Lep2FJ_run3-sig.json -a 1lep_2FJ -r 3 -t r3_1lep_2fj_sig_$DATE -j 4 -m 8 --files-per-job 5 --cutflow
./submit.py -c /home/users/aaarora/phys/run3/cmstas-run3-vbsvvh/preselection/etc/old_config/1Lep2FJ/1Lep2FJ_run3-sig.json -a 1lep_1FJ -r 3 -t r3_1lep_1fj_sig_$DATE -j 4 -m 8 --files-per-job 5 --cutflow
./submit.py -c /home/users/aaarora/phys/run3/cmstas-run3-vbsvvh/preselection/etc/old_config/0Lep3FJ/0Lep3FJ_run3-sig.json -a 0lep_3FJ -r 3 -t r3_0lep_3fj_sig_$DATE -j 4 -m 8 --files-per-job 5 --cutflow

# ./submit.py -c /home/users/aaarora/phys/run3/cmstas-run3-vbsvvh/preselection/etc/old_config/1Lep2FJ/1Lep2FJ_run2-data.json -a 1lep_2FJ -r 2 -t r2_1lep_2fj_data_$DATE -j 4 -m 8 --files-per-job 5 --cutflow
# ./submit.py -c /home/users/aaarora/phys/run3/cmstas-run3-vbsvvh/preselection/etc/old_config/1Lep2FJ/1Lep2FJ_run2-data.json -a 1lep_1FJ -r 2 -t r2_1lep_1fj_data_$DATE -j 4 -m 8 --files-per-job 5 --cutflow
./submit.py -c /home/users/aaarora/phys/run3/cmstas-run3-vbsvvh/preselection/etc/old_config/0Lep3FJ/0Lep3FJ_run2-data.json -a 0lep_3FJ -r 2 -t r2_0lep_3fj_data_$DATE -j 4 -m 8 --files-per-job 5 --cutflow

# ./submit.py -c /home/users/aaarora/phys/run3/cmstas-run3-vbsvvh/preselection/etc/old_config/1Lep2FJ/1Lep2FJ_run3-data.json -a 1lep_2FJ -r 3 -t r3_1lep_2fj_data_$DATE -j 4 -m 8 --files-per-job 5 --cutflow
# ./submit.py -c /home/users/aaarora/phys/run3/cmstas-run3-vbsvvh/preselection/etc/old_config/1Lep2FJ/1Lep2FJ_run3-data.json -a 1lep_1FJ -r 3 -t r3_1lep_1fj_data_$DATE -j 4 -m 8 --files-per-job 5 --cutflow
./submit.py -c /home/users/aaarora/phys/run3/cmstas-run3-vbsvvh/preselection/etc/old_config/0Lep3FJ/0Lep3FJ_run3-data.json -a 0lep_3FJ -r 3 -t r3_0lep_3fj_data_$DATE -j 4 -m 8 --files-per-job 5 --cutflow

# ./submit.py -c /home/users/aaarora/phys/run3/cmstas-run3-vbsvvh/preselection/etc/old_config/1Lep2FJ/1Lep2FJ_run2-bkg.json -a 1lep_2FJ -r 2 -t r2_1lep_2fj_bkg_$DATE -j 4 -m 8 --files-per-job 4 --cutflow
# ./submit.py -c /home/users/aaarora/phys/run3/cmstas-run3-vbsvvh/preselection/etc/old_config/1Lep2FJ/1Lep2FJ_run2-bkg.json -a 1lep_1FJ -r 2 -t r2_1lep_1fj_bkg_$DATE -j 4 -m 8 --files-per-job 4 --cutflow
# ./submit.py -c /home/users/aaarora/phys/run3/cmstas-run3-vbsvvh/preselection/etc/old_config/0Lep3FJ/0Lep3FJ_run2-bkg.json -a 0lep_3FJ -r 2 -t r2_0lep_3fj_bkg_$DATE -j 4 -m 8 --files-per-job 4 --cutflow

# ./submit.py -c /home/users/aaarora/phys/run3/cmstas-run3-vbsvvh/preselection/etc/old_config/1Lep2FJ/1Lep2FJ_run3-bkg.json -a 1lep_2FJ -r 3 -t r3_1lep_2fj_bkg_$DATE -j 4 -m 8 --files-per-job 4 --cutflow
# ./submit.py -c /home/users/aaarora/phys/run3/cmstas-run3-vbsvvh/preselection/etc/old_config/1Lep2FJ/1Lep2FJ_run3-bkg.json -a 1lep_1FJ -r 3 -t r3_1lep_1fj_bkg_$DATE -j 4 -m 8 --files-per-job 4 --cutflow
# ./submit.py -c /home/users/aaarora/phys/run3/cmstas-run3-vbsvvh/preselection/etc/old_config/0Lep3FJ/0Lep3FJ_run3-bkg.json -a 0lep_3FJ -r 3 -t r3_0lep_3fj_bkg_$DATE -j 4 -m 8 --files-per-job 4 --cutflow
