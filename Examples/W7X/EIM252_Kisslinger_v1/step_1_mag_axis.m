coilset = 'coilset_w7x_Kisslinger_ideal_v1';
coil_current_array = [12989 12989 12989 12989 12989 0. 0. 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
rGuess = 5.90;

 [rAxis, current_2p5_tesla] = find_magnetic_axis(coilset, coil_current_array, rGuess)
 
