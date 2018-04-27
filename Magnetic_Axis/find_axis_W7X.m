function [rAxis, current_2p5_tesla] = find_axis_W7X(coil_current_array, rGuess)
% function [rAxis, current_2p5_tesla] = find_axis_W7X(current, coil_current_array, rGuess)


disp(['<----coil_current_array: ' num2str(coil_current_array)]);
disp(['<----Initial guess r: ' num2str(rGuess)]);

OPTIONS = optimset('MaxFunEvals', 1e4, 'TolX', 1e-9, 'TolFun', 1e-9, 'MaxIter', 1e4);
lb = 5.7;
ub = 6.1;

disp('<----Starting axis search');
rAxis = lsqnonlin(@axis_finder_function, rGuess, lb, ub, OPTIONS, coil_current_array);
disp('<----Finished axis search');
disp(['<----Axis location: r = ' num2str(rAxis)]);
    
[bx, by, bz] = calc_b_W7X_RPhiZ(rAxis, 0, 0, coil_current_array);
modB = sqrt(bx^2 + by^2 + bz^2);
current_2p5_tesla = (2.5 / modB) * coil_current_array;
disp(['<----Coil current array for 2.5 T on-axis at phi=0:' ...
    num2str(current_2p5_tesla)]);



function [return_value] = axis_finder_function(rStart, coil_current_array)
% function [delta_after_one_pass] = axis_finder_function(rStart, current, taper, USE_GRID_INTERP);
relTol = 1e-7; % Pretty tight-used for single machine line following
absTol = 1e-9;
options = odeset('RelTol', relTol, 'AbsTol', absTol, 'Stats', 'on');

phiSpan = [0:(pi/2):(2*pi)];
coords0 = [rStart 0];

disp(['<----Performing Line Following.  r = ' num2str(rStart)]);

[phi, coords] = ...
    ode113(@LineFollowDerivs_axis, phiSpan, coords0, options, coil_current_array);
rEnd = coords(end, 1); zEnd = coords(end, 2);
% r = coords(:, 1); z = coords(:, 2);
disp(['<----Done with single pass.']);
disp(['<----r = ' num2str(rEnd)]);
disp(['<----z = ' num2str(zEnd)]);
% disp(['r = ' num2str(r')]);
% disp(['z = ' num2str(z')]);

% area_after_one_pass = polyarea(r, z)

delta_after_one_pass = sqrt( (rStart - rEnd)^2 + (zEnd)^2 )

return_value = delta_after_one_pass;
% return_value = area_after_one_pass;

%==========================================================================
%==========================================================================

function [dcoords_dphi] = LineFollowDerivs_axis(phi, coords, coil_current_array)
% need to return drdphi and dzdphi (which are cylindrical coordinates...)
r = coords(1);
z = coords(2);

% Biot-Savart code
[bx, by, bz] = calc_b_W7X_RPhiZ(r, phi, z, coil_current_array);

br = bx*cos(phi) + by*sin(phi);
bphi_over_r = ( -bx*sin(phi) + by*cos(phi) ) / r;
dcoords_dphi = [br bz]' / bphi_over_r;


