coilset = 'coilset_w7x_v3';
coil_current_array = [13546 13546 13546 13546 13546 -4191 -4191 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
rGuess = 6.235

 [rAxis, current_2p5_tesla] = find_magnetic_axis(coilset, coil_current_array, rGuess)
 

<----Finished axis search
<----Axis location: r = 6.1669
<----Coil current array for 2.5 T on-axis at phi=0 (Amps): 14570.2583      14570.2583      14570.2583      14570.2583      14570.2583     -4507.89551     -4507.89551               0               0               0               0               0               0               0               0               0               0               0               0               0               0               0               0               0               0               0               0               0               0
rAxis =
     6.166929644955605e+00

