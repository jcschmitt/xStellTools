coilset = 'coilset_hsx_complete';

% Hill 5%
coil_current_array = [10116. 7081. 7081. 7081. 7081. 7081. 7081.];


rGuess = 1.44;

 [rAxis, current_2p5_tesla] = find_magnetic_axis(coilset, coil_current_array, rGuess)
 
% rAxis =   1.439699403602296