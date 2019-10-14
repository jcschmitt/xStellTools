rGuess = 2.4;
coilset = 'coilset_bila_03';
coil_current_array = 497292.;

[rAxis, current_2p5_tesla] = find_magnetic_axis(coilset, coil_current_array, rGuess)


% <----Finished axis search
% <----Axis location: r = 2.3885
% <----Coil current array for 2.5 T on-axis at phi=0 (Amps): 520145.3395
% rAxis =
%      2.388490331785801e+00
% current_2p5_tesla =
%      5.201453395455460e+05
