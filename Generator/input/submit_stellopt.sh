#!/bin/sh

mpiexec -n 8 --allow-run-as-root /bin/xstelloptv2 input.stell0
rm wout_stell0_opt* regcoil_nescout.stell0_opt*

