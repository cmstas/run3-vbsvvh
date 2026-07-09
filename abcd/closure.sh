#!/bin/bash

REGION=${REGION:-all}
NSCAN=${NSCAN:-50}

# Latest version_N dir (highest N) that has BOTH the data predictions and a
# regions.yaml (i.e. inference-on-data and signal-region optimization have run).
latest() {
    ls -d "$1"/version_* 2>/dev/null \
        | sort -t_ -k2 -n \
        | while read -r d; do [ -f "$d/predictions_single_data.csv" ] && [ -f "$d/regions.yaml" ] && echo "$d"; done \
        | tail -1
}

# run <output_tag> <channel> <logname>
run() {
    BASE=$(latest "output/$1/single")
    if [ -z "$BASE" ]; then
        echo "[skip] $1: no version_* with predictions_single_data.csv + regions.yaml" >&2
        return
    fi
    echo "[run]  $1 ($2) -> $BASE"
    python3 closure.py -i "$BASE/predictions_single_data.csv" \
        --channel "$2" \
        --regions "$BASE/regions.yaml" \
        --region "$REGION" \
        -n "$NSCAN" \
        -o "$BASE/closure" > "closure_$3.log" 2>&1 &
}

run 1LEP_2FJ_RUN2 2fj 1lep_boosted_run2
run 1LEP_2FJ_RUN3 2fj 1lep_boosted_run3
run 1LEP_1FJ_RUN2 1fj 1lep_semiresolved_run2
run 1LEP_1FJ_RUN3 1fj 1lep_semiresolved_run3
run 0LEP_3FJ_RUN2 3fj 0lep_boosted_run2
run 0LEP_3FJ_RUN3 3fj 0lep_boosted_run3

wait
echo "All closure jobs finished. Summaries: output/*/single/version_*/closure/closure_summary_*.yaml"
