rGuess = 1.44;
coilset = 'coilset_hsx_complete';
coil_current_array = [10000. 0. 0. 0. 0. 0. 0.];
[rAxis, current_2p5_tesla] = find_magnetic_axis(coilset, coil_current_array, rGuess)

% <----Finished axis search
% <----Axis location: r = 1.4454
% <----Coil current array for 2.5 T on-axis at phi=0 (Amps): 26804.3132               0               0               0               0               0               0
% rAxis =
%      1.445362672888437e+00
% current_2p5_tesla =
%   Columns 1 through 3
%      2.680431324216060e+04                         0                         0
%   Columns 4 through 6
%                          0                         0                         0
%   Column 7

0