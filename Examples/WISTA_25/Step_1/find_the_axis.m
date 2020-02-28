rGuess = 2.4;
coilset = 'coilset_wista_25';
coil_current_array = 639791.2;
[rAxis, current_2p5_tesla] = find_magnetic_axis(coilset, coil_current_array, rGuess)


% <----Finished axis search
% <----Axis location: r = 2.4449
% <----Coil current array for 2.5 T on-axis at phi=0 (Amps): 639791.2046
% rAxis =
%      2.444949148094872e+00
% current_2p5_tesla =
%      6.397912045987240e+05
