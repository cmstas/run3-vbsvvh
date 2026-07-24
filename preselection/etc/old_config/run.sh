#!/bin/bash

python3 make_config.py --channel 1Lep2FJ_run3 -n 64 --category sig --xsecs xsecs_13p6TeV.json
python3 make_config.py --channel 1Lep2FJ_run3 -n 64 --category bkg --xsecs xsecs_13p6TeV.json
python3 make_config.py --channel 1Lep2FJ_run3 -n 64 --category data --xsecs xsecs_13p6TeV.json

python3 make_config.py --channel 0Lep3FJ_run3 -n 64 --category sig --xsecs xsecs_13p6TeV.json
python3 make_config.py --channel 0Lep3FJ_run3 -n 64 --category bkg --xsecs xsecs_13p6TeV.json
python3 make_config.py --channel 0Lep3FJ_run3 -n 64 --category data --xsecs xsecs_13p6TeV.json

python3 make_config.py --channel 1Lep2FJ_run2 -n 64 --category sig --xsecs xsecs_13TeV.json
python3 make_config.py --channel 1Lep2FJ_run2 -n 64 --category bkg --xsecs xsecs_13TeV.json
python3 make_config.py --channel 1Lep2FJ_run2 -n 64 --category data --xsecs xsecs_13TeV.json

python3 make_config.py --channel 0Lep3FJ_run2 -n 64 --category sig --xsecs xsecs_13TeV.json
python3 make_config.py --channel 0Lep3FJ_run2 -n 64 --category bkg --xsecs xsecs_13TeV.json
python3 make_config.py --channel 0Lep3FJ_run2 -n 64 --category data --xsecs xsecs_13TeV.json
