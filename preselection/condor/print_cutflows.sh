#!/bin/bash

echo "Run3 SIG"
python3 sum_cutflows.py -p jobs/*sig*r3_0lep_3fj*/*1p5*/

echo "Run3 BKG"
python3 sum_cutflows.py -p jobs/*bkg*r3_0lep_3fj*/*/

echo "Run3 DATA"
python3 sum_cutflows.py -p jobs/*data*r3_0lep_3fj*/*/

echo "Run2 SIG"
python3 sum_cutflows.py -p jobs/*sig*r2_0Lep_3fj*/*1p5*/

echo "Run2 BKG"
python3 sum_cutflows.py -p jobs/*bkg*r2_0lep_3fj*/*/

echo "Run2 DATA"
python3 sum_cutflows.py -p jobs/*data*r2_0lep_3fj*/*/

