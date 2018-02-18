function [Bx, By, Bz, BR, BPhi] = calc_b_HSX(P_x, P_y, P_z, main_current, taper)
% function [Bx, By, Bz, Br, Bphi] = calc_b_W7X(P_R, P_Phi, P_Z,
% main_current, taper)
%
% Input: P_R, P_Phi, P_Z:  Cylindrical coordinates of the observation point
%        main_current:  The current in the main field coils (per turn)
%        taper: 'Taper Array' to define the ratio of the coil courrent in
%        the auxilliarry coilset: 
%               [Coil_1 Coil_2 Coil_3 Coil_4 Coil_5 Coil_6]
%  
%
% Output:  Bx, By, Bz, BR, BPhi:  Magnetic flux density (T/m) at
%          observation point 
%
% Calculates cartesian B field from given cartesian lab frame
% coordinates and coil current array.  If the observation point lies
% on a current filament (or extremely close to it) then the contribution
% from that filament is ignored.
%
% The HSX coil fileset (coils.mat aux_coils.mat) must be available to this
% function
%
%

% The details of the HSX coils are stored in this persistent array
persistent HSX_coils

if isempty(HSX_coils)  % load the HSX coil information if necessary
    HSX_coils = load_HSX_coils;
end

if nargin < 4
    % If function is called without a taper array, assume all zeros (QHS)
    taper = [0 0 0 0 0 0];
end

% Factor of '14' comes from # of windings in the aux coils that are NOT
% accounted for in the coil description.
aux_current_array = 14 * main_current * taper;

% Calculate the field components due to individual coil sets
[Bx_Coil_main, By_Coil_main, Bz_Coil_main] = ...
    calc_B_BiotSavart(P_x, P_y, P_z, HSX_coils.main, main_current);

[Bx_AuxCoil_1, By_AuxCoil_1, Bz_AuxCoil_1] = ...
    calc_B_BiotSavart(P_x, P_y, P_z, HSX_coils.aux1, aux_current_array(1));

[Bx_AuxCoil_2, By_AuxCoil_2, Bz_AuxCoil_2] = ...
    calc_B_BiotSavart(P_x, P_y, P_z, HSX_coils.aux2, aux_current_array(2));

[Bx_AuxCoil_3, By_AuxCoil_3, Bz_AuxCoil_3] = ...
    calc_B_BiotSavart(P_x, P_y, P_z, HSX_coils.aux3, aux_current_array(3));

[Bx_AuxCoil_4, By_AuxCoil_4, Bz_AuxCoil_4] = ...
    calc_B_BiotSavart(P_x, P_y, P_z, HSX_coils.aux4, aux_current_array(4));

[Bx_AuxCoil_5, By_AuxCoil_5, Bz_AuxCoil_5] = ...
    calc_B_BiotSavart(P_x, P_y, P_z, HSX_coils.aux5, aux_current_array(5));

[Bx_AuxCoil_6, By_AuxCoil_6, Bz_AuxCoil_6] = ...
    calc_B_BiotSavart(P_x, P_y, P_z, HSX_coils.aux6, aux_current_array(6));

% Add them together.
Bx = Bx_AuxCoil_1 +  Bx_AuxCoil_2 +  Bx_AuxCoil_3 + Bx_AuxCoil_4 + ...
    Bx_AuxCoil_5 + Bx_AuxCoil_6 + Bx_Coil_main;

By = By_AuxCoil_1 +  By_AuxCoil_2 + By_AuxCoil_3 + By_AuxCoil_4 + ...
    By_AuxCoil_5 + By_AuxCoil_6 + By_Coil_main;

Bz = Bz_AuxCoil_1 + Bz_AuxCoil_2 + Bz_AuxCoil_3 +  Bz_AuxCoil_4 + ...
    Bz_AuxCoil_5 + Bz_AuxCoil_6 + Bz_Coil_main;

% The other useful values.
P_phi = atan2(P_y, P_x);
BR = Bx .* cos(P_phi) + By .* sin(P_phi);
BPhi = By .* cos(P_phi) - Bx .* sin(P_phi);


