#!/bin/bash

# Return the path of the latest version_N training dir (highest N) that has predictions.
latest() {
    ls -d "$1"/version_* 2>/dev/null \
        | sort -t_ -k2 -n \
        | while read -r d; do [ -f "$d/predictions_single.csv" ] && echo "$d"; done \
        | tail -2 | head -1
}

BASE=$(latest output/1LEP_2FJ_RUN2/single); echo "Using $(basename $BASE) for 1LEP_2FJ_RUN2"; python3 2fj_analysis.py -i $BASE/predictions_single.csv -d $BASE/predictions_single_data.csv -n 5 --data-driven > config_1lep_boosted_run2.log &
BASE=$(latest output/1LEP_2FJ_RUN3/single); echo "Using $(basename $BASE) for 1LEP_2FJ_RUN3"; python3 2fj_analysis.py -i $BASE/predictions_single.csv -d $BASE/predictions_single_data.csv -n 5 --data-driven > config_1lep_boosted_run3.log &
BASE=$(latest output/1LEP_1FJ_RUN2/single); echo "Using $(basename $BASE) for 1LEP_1FJ_RUN2"; python3 1fj_analysis.py -i $BASE/predictions_single.csv -d $BASE/predictions_single_data.csv -n 5 --data-driven > config_1lep_semiresolved_run2.log &
BASE=$(latest output/1LEP_1FJ_RUN3/single); echo "Using $(basename $BASE) for 1LEP_1FJ_RUN3"; python3 1fj_analysis.py -i $BASE/predictions_single.csv -d $BASE/predictions_single_data.csv -n 5 --data-driven > config_1lep_semiresolved_run3.log &
BASE=$(latest output/0LEP_3FJ_RUN2/single); echo "Using $(basename $BASE) for 0LEP_3FJ_RUN2"; python3 3fj_analysis.py -i $BASE/predictions_single.csv -d $BASE/predictions_single_data.csv -n 5 --data-driven > config_0lep_boosted_run2.log &
BASE=$(latest output/0LEP_3FJ_RUN3/single); echo "Using $(basename $BASE) for 0LEP_3FJ_RUN3"; python3 3fj_analysis.py -i $BASE/predictions_single.csv -d $BASE/predictions_single_data.csv -n 5 --data-driven > config_0lep_boosted_run3.log &
