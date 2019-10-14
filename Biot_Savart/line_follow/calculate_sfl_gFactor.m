function g_Boozer = ... 
    calculate_sfl_gFactor(coilsetID, coilCurrents, rStart, relTol, absTol)

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
    ode113(@Boozer_gFactorDeriv, phiSpan, coords0, options, coilsetID, coilCurrents);

g_Boozer = coords(end,3) / (2*pi);
disp(['Boozer g-factor = ', num2str(g_Boozer) ]);
    

%==========================================================================
%==========================================================================

function [dcoords_dphi] = Boozer_gFactorDeriv(phi, coords, coilsetID, coilCurrents)
% need to return drdphi and dzdphi (which are cylindrical coordinates...)
r = coords(1);
z = coords(2);

% Biot-Savart code
[bx, by, bz] = calc_b_RPhiZ(coilsetID, r, phi, z, coilCurrents);

B_squared = bx.^2 + by.^2 + bz.^2;
br = bx*cos(phi) + by*sin(phi);
bphi_over_r = (-bx*sin(phi) + by*cos(phi)) / r;
dr_dphi = br / bphi_over_r;
dz_dphi = bz / bphi_over_r;
dg_dphi = B_squared / bphi_over_r;
dcoords_dphi = [dr_dphi dz_dphi dg_dphi]';
