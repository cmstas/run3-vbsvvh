#!/bin/bash

# python3 2fj_analysis.py -i output/1LEP_2FJ_RUN2/single/version_0/predictions_single.csv -d output/1LEP_2FJ_RUN2/single/version_0/predictions_single_data.csv -n 5 &
# python3 2fj_analysis.py -i output/1LEP_2FJ_RUN3/single/version_4/predictions_single.csv -d output/1LEP_2FJ_RUN3/single/version_4/predictions_single_data.csv -n 5 &
python3 1fj_analysis.py -i output/1LEP_1FJ_RUN2/single/version_0/predictions_single.csv -d output/1LEP_1FJ_RUN2/single/version_0/predictions_single_data.csv -n 5 &
python3 1fj_analysis.py -i output/1LEP_1FJ_RUN3/single/version_0/predictions_single.csv -d output/1LEP_1FJ_RUN3/single/version_0/predictions_single_data.csv -n 5 &
# python3 3fj_analysis.py -i output/0LEP_3FJ_RUN2/single/version_0/predictions_single.csv -d output/0LEP_3FJ_RUN2/single/version_0/predictions_single_data.csv -n 5 &
# python3 3fj_analysis.py -i output/0LEP_3FJ_RUN3/single/version_0/predictions_single.csv -d output/0LEP_3FJ_RUN3/single/version_0/predictions_single_data.csv -n 5 &