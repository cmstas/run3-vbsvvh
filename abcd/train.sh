#!/bin/bash

for cfg in single/*0lep*; do nohup python3 main.py --config $cfg --flavor single > $(basename $cfg .yaml).log & done