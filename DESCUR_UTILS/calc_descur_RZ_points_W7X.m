function [surf_x, surf_y, surf_z, phiout] = ...
    calc_descur_RZ_points_W7X(rboundary, nphi, ntheta, coil_current_array)
% function [surf_x, surf_y, surf_z, phiout] = ...
%     calc_descur_RZ_points_W7X(rboundary, nphi, ntheta, coil_current_array)
% inputs:
%   rboundary- R coord of boundary at the boxport
%   nphi- number of toridal cuts
%   ntheta- number of points on each cut
%   current- main coil current in amps (Check the following claim: + for clockwise, - for counter-clockwise)
%
% outputs:
%   surf_x, surf_y, surf_z: x,y,z coords of the points on the surface
%   phiout - the phi coordinate of each point
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For time-keeping puposes
tic;

% Starting position of field line to follow [r0 z0]
x0 = [rboundary 0];

% Values of phi at which B is to be calculated. (At the symmetry plane)
phi_points= 0 : (2*pi/nphi) : (ntheta*2*pi);

% Integrate the ode from field_line_derivs.
% NOTE: the tolerances here make a big difference, so be careful! These values are
% hand-picked, but they should be check/verified for other cases. 
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);

% The output 'xout' is [r z]
% Perform the field line following.
[phiout, xout] = ode45(@field_line_derivs, phi_points, x0, options, coil_current_array); 

% Extract the data into usful arrays
r_surf = xout(:,1);
z_surf = xout(:,2);
surf_x = [r_surf.*cos(phiout)]';
surf_y = [r_surf.*sin(phiout)]';
surf_z = [z_surf]';

% Report the time that was required (wall time)
toc


function [dxdz] = field_line_derivs(phi, x, coil_current_array)
% function [dxdz] = field_line_derivs(phi, x, coil_current_array)
%
% This returns the ode for a field line, which can be solved with one 
% of matlab's ode solvers.
%
% input:
%   phi: the independent variable to be integrate over
%   x: column vector containing r and z (x = [r z])
%   current: main field coil current
%   taper: taper array for auxilliary coil currents
%
% output:
%   dxdz: row vector containg the derivatives of r and z w/ respect to 
%         phi. (dxdz = [drdphi dzdphi]')
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = x(1);
z = x(2);

[bx, by, bz] = calc_b_W7X_RPhiZ(r, phi, z, coil_current_array);

br = bx * cos(phi) + by * sin(phi);
bphi = -bx * sin(phi) + by * cos(phi);

% Return derivs
drdphi = (r .* br) ./ bphi;
dzdphi = (r .* bz) ./ bphi;
dxdz(1) = drdphi;
dxdz(2) = dzdphi;
dxdz = dxdz';