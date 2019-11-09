rGuess = 2.4;
coilset = 'coilset_wista_15';
coil_current_array = 632955.5;
[rAxis, current_2p5_tesla] = find_magnetic_axis(coilset, coil_current_array, rGuess)
% 
% <----Finished axis search
% <----Axis location: r = 2.4415
% <----Coil current array for 2.5 T on-axis at phi=0 (Amps): 632955.5269
% rAxis =
%      2.441491543155680e+00
% current_2p5_tesla =
%      6.329555268532927e+05