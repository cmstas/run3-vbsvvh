#!/bin/bash

# 0Lep3FJ

echo "0Lep3FJ Run3 SIG k2V=1.0, kL=1.0"
python3 sum_cutflows.py -p jobs/*sig*r3_0lep_3fj*/*C2V_1p0_C3_1p0*/
echo "0Lep3FJ Run3 SIG k2V=1.5, kL=1.0"
python3 sum_cutflows.py -p jobs/*sig*r3_0lep_3fj*/*1p5*/
echo "0Lep3FJ Run3 SIG k2V=1.0, kL=10"
python3 sum_cutflows.py -p jobs/*sig*r3_0lep_3fj*/*10p0*/
echo "0Lep3FJ Run3 SIG k2V=0.25, kL=1.0"
python3 sum_cutflows.py -p jobs/*sig*r3_0lep_3fj*/*0p25*/
echo "0Lep3FJ Run3 SIG k2V=0.75, kL=1.0"
python3 sum_cutflows.py -p jobs/*sig*r3_0lep_3fj*/*0p75*/
echo "0Lep3FJ Run3 BKG"
python3 sum_cutflows.py -p jobs/*bkg*r3_0lep_3fj*/*/
echo "0Lep3FJ Run3 DATA"
python3 sum_cutflows.py -p jobs/*data*r3_0lep_3fj*/*/


echo "0Lep3FJ Run2 SIG k2V=1.0, kL=1.0"
python3 sum_cutflows.py -p jobs/*sig*r2_0Lep_3fj*/*C2V_1p0_C3_1p0*/
echo "0Lep3FJ Run2 SIG k2V=1.5, kL=1.0"
python3 sum_cutflows.py -p jobs/*sig*r2_0Lep_3fj*/*1p5*/
echo "0Lep3FJ Run2 SIG k2V=1.0, kL=10"
python3 sum_cutflows.py -p jobs/*sig*r2_0Lep_3fj*/*10p0*/
echo "0Lep3FJ Run2 BKG"
python3 sum_cutflows.py -p jobs/*bkg*r2_0lep_3fj*/*/
echo "0Lep3FJ Run2 DATA"
python3 sum_cutflows.py -p jobs/*data*r2_0lep_3fj*/*/

# 1Lep1FJ

echo "1Lep1FJ Run3 SIG k2V=1.0, kL=1.0"
python3 sum_cutflows.py -p jobs/*sig*r3_1lep_1fj*/*C2V_1p0_C3_1p0*/
echo "1Lep1FJ Run3 SIG k2V=1.5, kL=1.0"
python3 sum_cutflows.py -p jobs/*sig*r3_1lep_1fj*/*1p5*/
echo "1Lep1FJ Run3 SIG k2V=1.0, kL=10"
python3 sum_cutflows.py -p jobs/*sig*r3_1lep_1fj*/*10p0*/
echo "1Lep1FJ Run3 SIG k2V=0.25, kL=1.0"
python3 sum_cutflows.py -p jobs/*sig*r3_1lep_1fj*/*0p25*/
echo "1Lep1FJ Run3 SIG k2V=0.75, kL=1.0"
python3 sum_cutflows.py -p jobs/*sig*r3_1lep_1fj*/*0p75*/
echo "1Lep1FJ Run3 BKG"
python3 sum_cutflows.py -p jobs/*bkg*r3_1lep_1fj*/*/
echo "1Lep1FJ Run3 DATA"
python3 sum_cutflows.py -p jobs/*data*r3_1lep_1fj*/*/


echo "1Lep1FJ Run2 SIG k2V=1.0, kL=1.0"
python3 sum_cutflows.py -p jobs/*sig*r2_1lep_1fj*/*C2V_1p0_C3_1p0*/
echo "1Lep1FJ Run2 SIG k2V=1.5, kL=1.0"
python3 sum_cutflows.py -p jobs/*sig*r2_1lep_1fj*/*1p5*/
echo "1Lep1FJ Run2 SIG k2V=1.0, kL=10"
python3 sum_cutflows.py -p jobs/*sig*r2_1lep_1fj*/*10p0*/
echo "1Lep1FJ Run2 BKG"
python3 sum_cutflows.py -p jobs/*bkg*r2_1lep_1fj*/*/
echo "1Lep1FJ Run2 DATA"
python3 sum_cutflows.py -p jobs/*data*r2_1lep_1fj*/*/ 

# 1Lep2FJ
echo "1Lep2FJ Run3 SIG k2V=1.0, kL=1.0"
python3 sum_cutflows.py -p jobs/*sig*r3_1lep_2fj*/*C2V_1p0_C3_1p0*/
echo "1Lep2FJ Run3 SIG k2V=1.5, kL=1.0"
python3 sum_cutflows.py -p jobs/*sig*r3_1lep_2fj*/*1p5*/
echo "1Lep2FJ Run3 SIG k2V=1.0, kL=10"
python3 sum_cutflows.py -p jobs/*sig*r3_1lep_2fj*/*10p0*/
echo "1Lep2FJ Run3 SIG k2V=0.25, kL=1.0"
python3 sum_cutflows.py -p jobs/*sig*r3_1lep_2fj*/*0p25*/
echo "1Lep2FJ Run3 SIG k2V=0.75, kL=1.0"
python3 sum_cutflows.py -p jobs/*sig*r3_1lep_2fj*/*0p75*/
echo "1Lep2FJ Run3 BKG"
python3 sum_cutflows.py -p jobs/*bkg*r3_1lep_2fj*/*/
echo "1Lep2FJ Run3 DATA"
python3 sum_cutflows.py -p jobs/*data*r3_1lep_2fj*/*/

echo "1Lep2FJ Run2 SIG k2V=1.0, kL=1.0"
python3 sum_cutflows.py -p jobs/*sig*r2_1lep_2fj*/*C2V_1p0_C3_1p0*/
echo "1Lep2FJ Run2 SIG k2V=1.5, kL=1.0"
python3 sum_cutflows.py -p jobs/*sig*r2_1lep_2fj*/*1p5*/
echo "1Lep2FJ Run2 SIG k2V=1.0, kL=10"
python3 sum_cutflows.py -p jobs/*sig*r2_1lep_2fj*/*10p0*/
echo "1Lep2FJ Run2 BKG"
python3 sum_cutflows.py -p jobs/*bkg*r2_1lep_2fj*/*/
echo "1Lep2FJ Run2 DATA"
python3 sum_cutflows.py -p jobs/*data*r2_1lep_2fj*/*/