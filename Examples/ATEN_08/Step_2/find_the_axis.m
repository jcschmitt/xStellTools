rGuess = 2.4;
coilset = 'coilset_aten25t_08_stelsym';
coil_current_array = 635981.;

[rAxis, current_2p5_tesla] = find_magnetic_axis(coilset, coil_current_array, rGuess)

% 
% <----Finished axis search
% <----Axis location: r = 2.4416
% <----Coil current array for 2.5 T on-axis at phi=0 (Amps): 632451.6968
% rAxis =
%      2.441615140689149e+00
% current_2p5_tesla =
%      6.324516967742629e+05
