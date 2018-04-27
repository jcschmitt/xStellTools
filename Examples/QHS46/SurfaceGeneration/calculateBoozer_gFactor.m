function g_Boozer = ... 
    calculateBoozer_gFactor(earthField, current, taper, rStart, relTol, absTol);

if (nargin < 5) % if you don't specify the tolerances
    disp('Warning.  Tolerances not specified.');
	relTol = 1e-3;  % Real loose-used for testing
	absTol = 1e-5;
	disp('Choosing lowest tolerances.');
else
    disp('Using specified tolerances.');
end


options = odeset('RelTol', relTol, 'AbsTol', absTol, 'Stats', 'on');

phiSpan = [0 2*pi];
coords0 = [rStart 0 0];  % [r z g_Boozer]

[phi, coords] = ...
    ode113(@Boozer_gFactorDeriv, phiSpan, coords0, options, earthField, current, taper);

g_Boozer = coords(end,3) / (2*pi);
disp(['Boozer g-factor = ', num2str(g_Boozer) ]);
    

%==========================================================================
%==========================================================================

function [dcoords_dphi] = Boozer_gFactorDeriv(phi, coords, earthField, current, taper);
% need to return drdphi and dzdphi (which are cylindrical coordinates...)
r = coords(1);
z = coords(2);
% The following line calls a modified version of bs_dervs_aux_mex that
% includes the Earth's magnetic field.  Not yet implemented correctly in
% this revision
% [bx, by, bz] = bs_derivs_aux_mex(r, phi, z, earthField, current, taper);
% If taper == [0 0 0 0 0 0], don't input the taper array into the
% solver--The solver runs faster.
if (taper == [0 0 0 0 0 0])
    [bx, by, bz] = calc_b_bs_JL(r, phi, z, current);
else
    [bx, by, bz] = calc_b_bs_JL(r, phi, z, current, taper);
end

B_squared = bx.^2 + by.^2 + bz.^2;
br = bx*cos(phi) + by*sin(phi);
bphi_over_r = (-bx*sin(phi) + by*cos(phi)) / r;
dr_dphi = br / bphi_over_r;
dz_dphi = bz / bphi_over_r;
dg_dphi = B_squared / bphi_over_r;
dcoords_dphi = [dr_dphi dz_dphi dg_dphi]';
