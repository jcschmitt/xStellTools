#!/bin/tcsh
#This file is called submit-script.sh
#SBATCH --partition=general
#SBATCH --time=24:00:00   # run time in days-hh:mm:ss
#SBATCH --ntasks=8          # default 16 if this line not specified
#SBATCH --mem-per-cpu=2gb
#SBATCH -o out.slurm-%j.%N #
#SBATCH -e err.slurm-%j.%N #
#Make sure to change the above two lines to reflect your appropriate
# file locations for standard error and output

source /p/stellopt/ANALYSIS/jschmitt/src/module_fotd
rm log.vmec
mpirun -n 8 --mca mtl ^psm --mca btl ^openib  /p/stellopt/ANALYSIS/jschmitt/src/stellopt/VMEC2000/Release/xvmec2000 input.aten_a3b25 >& log.vmec
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz0_nboz0 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz0_nboz1 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz0_nboz2 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz0_nboz3
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz0_nboz4 & 
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz0_nboz5 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz0_nboz6 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz0_nboz7
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz0_nboz8 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz1_nboz0 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz1_nboz1 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz1_nboz2
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz1_nboz3 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz1_nboz4 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz1_nboz5 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz1_nboz6
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz1_nboz7 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz1_nboz8 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz2_nboz0 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz2_nboz1
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz2_nboz2 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz2_nboz3 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz2_nboz4 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz2_nboz5
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz2_nboz6 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz2_nboz7 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz2_nboz8 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz3_nboz0
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz3_nboz1 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz3_nboz2 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz3_nboz3 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz3_nboz4
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz3_nboz5 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz3_nboz6 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz3_nboz7 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz3_nboz8
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz4_nboz0 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz4_nboz1 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz4_nboz2 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz4_nboz3
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz4_nboz4 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz4_nboz5 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz4_nboz6 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz4_nboz7
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz4_nboz8 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz5_nboz0 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz5_nboz1 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz5_nboz2
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz5_nboz3 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz5_nboz4 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz5_nboz5 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz5_nboz6
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz5_nboz7 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz5_nboz8 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz6_nboz0 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz6_nboz1
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz6_nboz2 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz6_nboz3 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz6_nboz4 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz6_nboz5
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz6_nboz6 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz6_nboz7 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz6_nboz8 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz7_nboz0
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz7_nboz1 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz7_nboz2 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz7_nboz3 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz7_nboz4
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz7_nboz5 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz7_nboz6 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz7_nboz7 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz7_nboz8
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz8_nboz0 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz8_nboz1 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz8_nboz2 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz8_nboz3
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz8_nboz4 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz8_nboz5 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz8_nboz6 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz8_nboz7 &
/p/stellopt/ANALYSIS/jschmitt/src/stellopt/BOOZ_XFORM/Release/xbooz_xform in_booz.aten_a3b25_mboz8_nboz8
 
