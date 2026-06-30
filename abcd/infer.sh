#!/bin/bash

for cfg in single/*; do nohup python3 main.py --config $cfg --flavor single --infer > $(basename $cfg .yaml) & done
for cfg in single/*; do nohup python3 main.py --config $cfg --flavor single --infer --data > $(basename $cfg .yaml) & done
