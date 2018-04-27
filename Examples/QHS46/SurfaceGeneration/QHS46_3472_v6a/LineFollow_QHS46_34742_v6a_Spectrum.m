function [chi, coords] = ... 
    LineFollow_QHS46_34742_v6a_Spectrum(spectrumType, current, rStart, ...
    phiStart, zStart, chiStart, chiInc, chiEnd, relTol, absTol)
% function [chi, coords] = ... 
%     LineFollow_QHS46_34742_v6a_Spectrum(spectrumType, current, rStart, phiStart, zStart, chiStart, chiInc, chiEnd, relTol, absTol);
% 
% 1) Generates data to help determine the Hamada or Boozer spectrum of the magnetic
%       field .
%
% --------------
%   Input
% --------------
% spectrumType  'Hamada' or 'Boozer'
% current       main coil current, in amps
% rStart        r coordinate of starting point of field line
% zStart        z coordinate of starting point of field line
% chiStart      chi of starting point (usually zero)
% chiInc        increment of chi for each point returned in phi and coords
%               If phiInc > phiEnd, an error occurs
%               If phiInc > 0.5 * phiEnd, the returned values of
%               phi and coords will be the ones that were chosen by the ODE
%               routine, not by phiInc.  The last point will be phiInc
%               If phiInc < 0.5 * phiEnd, the returned values of
%               phi and coords will be the ones for
%               phiStart:phiInc:phiEnd
%               If phiInc = 0, an error occurs
% chiEnd        the largest angle to integrate to. 2*pi = 1 transit.
% relTol        relative tolerance for ODE routine
% absTol        absolute tolerance for oDE routine
%
% -------------
%  Output
% -------------
% chi           phi will contain either the values of phi that are chosen
%               by the user (phiStart:phiInc:phiEnd), or it will
%               contain all of the points chosen by the ODE routine
% coords        coords will contain the r, phi and z coordinates corresponding
%               to the values of chi
%
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

if (nargin < 9) % if you don't specify the tolerances
    disp('Warning.  Tolerances not specified.');
	relTol = 1e-3;  % Real loose-used for testing
	absTol = 1e-5;
	disp('Choosing lowest tolerances.');
else
    disp('Using specified tolerances.');
end


options = odeset('RelTol', relTol, 'AbsTol', absTol, 'Stats', 'on');

chiSpan = [chiStart:chiInc:chiEnd];
coords0 = [rStart, phiStart, zStart];

if strcmp(lower(spectrumType), 'hamada')
    [chi, coords] = ...
        ode113(@LineFollowDerivs_Hamada, chiSpan, coords0, options, current);
elseif strcmp(lower(spectrumType), 'boozer')
    [chi, coords] = ...
        ode113(@LineFollowDerivs_Boozer, chiSpan, coords0, options, current);
else
    error('Unknown spectrum request');
end
    

%==========================================================================
%==========================================================================

function [dcoords_dchi] = LineFollowDerivs_Hamada(chi, coords, current)
% need to return drdphi and dzdphi (which are cylindrical coordinates...)
r = coords(1);
phi = coords(2);
z = coords(3);


% Biot-Savart code
[bx, by, bz] = calc_b_QHS46_RPhiZ(r, phi, z, current);

br = bx*cos(phi) + by*sin(phi);
bphi_over_r = (-bx*sin(phi) + by*cos(phi)) / r;
dcoords_dchi = [br bphi_over_r bz]';

function [dcoords_dchi] = LineFollowDerivs_Boozer(chi, coords, current)
% need to return drdphi and dzdphi (which are cylindrical coordinates...)
r = coords(1);
phi = coords(2);
z = coords(3);

% Biot-Savart code
[bx, by, bz] = calc_b_QHS46_RPhiZ(r, phi, z, current);

B_squared = bx.^2 + by.^2 + bz.^2;
br = bx*cos(phi) + by*sin(phi);
bphi_over_r = (-bx*sin(phi) + by*cos(phi)) / r;
dcoords_dchi = [br bphi_over_r bz]' / B_squared;
