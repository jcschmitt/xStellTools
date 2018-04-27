function [phi, coords] = ... 
    LineFollow_QHS46_34742_v6a(current, rStart, zStart, phiStart, phiInc, phiEnd, relTol, absTol)
% function [phi, coords] = ... 
%     LineFollow_QHS46_34742_v6a(current, rStart, zStart, phiStart,
%     phiInc, phiEnd, relTol, absTol)
% 
% 1) Generates data to help determine the values of Iota (rotational
%       transform) of the magnetic field  of HSX.
%
% --------------
%   Input
% --------------
% current       main coil current, in amps
% rStart        r coordinate of starting point of field line
% zStart        z coordinate of starting point of field line
% phiStart      phi of starting point (usually zero)
% phiInc        increment of phi for each point returned in phi and coords
%               If phiInc > phiEnd, an error occurs
%               If phiInc > 0.5 * phiEnd, the returned values of
%               phi and coords will be the ones that were chosen by the ODE
%               routine, not by phiInc.  The last point will be phiInc
%               If phiInc < 0.5 * phiEnd, the returned values of
%               phi and coords will be the ones for
%               phiStart:phiInc:phiEnd
%               If phiInc = 0, an error occurs
% phiEnd        the largest angle to integrate to. 2*pi = 1 transit.
% relTol        relative tolerance for ODE routine
% absTol        absolute tolerance for oDE routine
%
% -------------
%  Output
% -------------
% phi           phi will contain either the values of phi that are chosen
%               by the user (phiStart:phiInc:phiEnd), or it will
%               contain all of the points chosen by the ODE routine
% coords        coords will contain the r and z coordinates along with 
%               Int(dl/B) along the field line.
%               to the values of phi
%
% -------------
% Requires
% -------------
%
%%%%%%%%%%%%%%%%%%%%%%%
% Last updated on:
%%%%%%%%%%%%%%%%%%%%%%%

%---------------------------
% Options for the ODE solver
%---------------------------
% relTol = 1e-8; % Real tight-used for mulitple machine line following
% absTol = 1e-10;
% relTol = 1e-7; % Pretty tight-used for single machine line following
% absTol = 1e-9;
% relTol = 1e-3;  % Real loose-used for testing
% absTol = 1e-5;

if (nargin < 8) % if you don't specify the tolerances
    disp('Warning.  Tolerances not specified.');
	relTol = 1e-3;  % Real loose-used for testing
	absTol = 1e-5;
	disp('Choosing lowest tolerances.');
else
    disp('Using specified tolerances.');
end

options = odeset('RelTol', relTol, 'AbsTol', absTol, 'Stats', 'on');

phiSpan = [phiStart:phiInc:phiEnd];
coords0 = [rStart, zStart, 0];

[phi, coords] = ...
    ode113(@LineFollowDerivs, phiSpan, coords0, options, current);

%==========================================================================
%==========================================================================

function [dcoords_dphi] = LineFollowDerivs(phi, coords, current)
% need to return drdphi and dzdphi (which are cylindrical coordinates...)
r = coords(1);
z = coords(2);
% The following line calls a modified version of bs_dervs_aux_mex that
% includes the Earth's magnetic field.  Not yet implemented correctly in
% this revision
% [bx, by, bz] = bs_derivs_aux_mex(r, phi, z, earthField, current, taper);
% If taper == [0 0 0 0 0 0], don't input the taper array into the
% solver--The solver runs faster.
% Biot-Savart code
[bx, by, bz] = calc_b_QHS46_RPhiZ(r, phi, z, current);
% Grid-Interpolation code
% [bx, by, bz, delb] = hsxfield_grid_interp(r, phi, z, current, taper);

br = bx*cos(phi) + by*sin(phi);
bphi_over_r = ( -bx*sin(phi) + by*cos(phi) ) / r;
% r_over_bphi = r / bphi;  % This is dl/B
dcoords_dphi = [br bz 1]' / bphi_over_r;
