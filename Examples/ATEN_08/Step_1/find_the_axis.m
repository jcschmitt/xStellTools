rGuess = 2.4;
coilset = 'coilset_aten25t_08';
coil_current_array = 635981.;

[rAxis, current_2p5_tesla] = find_magnetic_axis(coilset, coil_current_array, rGuess)

% 
% <----Finished axis search
% <----Axis location: r = 2.4416
% <----Coil current array for 2.5 T on-axis at phi=0 (Amps): 632468.6738
% rAxis =
%      2.441634981219627e+00
% current_2p5_tesla =
%      6.324686738183743e+055
