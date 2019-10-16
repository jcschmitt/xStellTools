rGuess = 2.41;
coilset = 'coilset_bila_03_stelsym';
coil_current_array = 497292.;

[rAxis, current_2p5_tesla] = find_magnetic_axis(coilset, coil_current_array, rGuess)
% 
% <----Finished axis search
% <----Axis location: r = 2.3881
% <----Coil current array for 2.5 T on-axis at phi=0 (Amps): 519966.7677
% rAxis =
%      2.388121565178873e+00
% current_2p5_tesla =
%      5.199667677158224e+05
%      