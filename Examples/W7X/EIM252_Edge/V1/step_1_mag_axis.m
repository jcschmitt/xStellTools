coilset = 'coilset_w7x_v3';
coil_current_array = [12989 12989 12989 12989 12989 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
rGuess = 6.235

 [rAxis, current_2p5_tesla] = find_magnetic_axis(coilset, coil_current_array, rGuess)
 
%<----Finished axis search
%<----Axis location: r = 6.2371
%<----Coil current array for 2.5 T on-axis at phi=0 (Amps): 14457.8053      14457.8053      14457.8053      14457.8053      14457.8053               0               0               0               0               0               0               0               0               0               0               0               0               0               0               0               0               0               0               0               0               0               0               0               0
%rAxis =     6.237141085932668e+00