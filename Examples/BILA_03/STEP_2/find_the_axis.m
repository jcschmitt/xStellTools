rGuess = 2.4;
coilset = 'coilset_bila_03_stelsym';
coil_current_array = 497292.;

[rAxis, current_2p5_tesla] = find_magnetic_axis(coilset, coil_current_array, rGuess)

% 
% <----Finished axis search
% <----Axis location: r = 2.3889
% <----Coil current array for 2.5 T on-axis at phi=0 (Amps): 521815.4147
% rAxis =
%      2.388919066924667e+00
% current_2p5_tesla =
%      5.218154146907073e+05
