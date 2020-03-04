coilset = 'coilset_w7x_v3';
coil_current_array = [13546 13546 13546 13546 13546 -4191 -4191 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
rGuess = 5.8

 [rAxis, current_2p5_tesla] = find_magnetic_axis(coilset, coil_current_array, rGuess)
 

%<----Finished axis search
%<----Axis location: r = 5.7663
%<----Coil current array for 2.5 T on-axis at phi=0 (Amps): 12148.8792      12148.8792      12148.8792      12148.8792      12148.8792     -3758.74448     -3758.74448               0               0               0               0               0               0               0               0               0               0               0               0               0               0               0               0               0               0               0               0               0               0
%rAxis =
%     5.766325108150850e+00

