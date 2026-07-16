#!/bin/bash

for cfg in single/*0lep*; do nohup python3 main.py --config $cfg --flavor single --infer >> $(basename $cfg .yaml).log & done
for cfg in single/*0lep*; do nohup python3 main.py --config $cfg --flavor single --infer --data >> $(basename $cfg .yaml).log & done
