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
mpirun -n 8 --mca mtl ^psm --mca btl ^openib  /p/stellopt/ANALYSIS/jschmitt/src/STELLOPT/STELLOPTV2/Release/xstelloptv2 input.stell0 >& log.stellopt
stellopt_tidyup
cp input.stell0.00000 input.stell0.rerun
mpirun -n 8 --mca mtl ^psm --mca btl ^openib  /p/stellopt/ANALYSIS/jschmitt/src/STELLOPT/VMEC2000/Release/xvmec2000 input.stell0.rerun >& log.vmec
/p/stellopt/ANALYSIS/jschmitt/src/STELLOPT/COBRAVMEC/Release/xcobravmec in_cobra.stell0.rerun >& log.cobra

