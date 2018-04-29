function [Bx,By,Bz,BR,BPhi] = calc_b_QHS46(P_x, P_y, P_z, coil_current_array)
%function [Bx,By,Bz,Br,Bphi] = calc_b_QHS46(P_x, P_y, P_z, 

% Input: Px, Py, Pz:  Cartesian coordinates of the observation point
%        coil_current_array:  A number: Coil_1
%
% Output:  Bx, By, Bz:  Magnetic flux density (T/m) at observation point
%          Br, Bphi: Bx and By are used to find radial and toroidal
%          components.
%
% Calculates cartesian B field from given cartesian lab frame
% coordinates and coil current array.  If the observation point lies
% on a current filament (or extremely close) the contribution from that
% filament is ignored.
%
% The file QHS46_xxx_coils.mat must be available.
%

% The details of the QHS46 coils are stored in this persistent array
persistent QHS46_coils

if isempty(QHS46_coils)  % load the QHS46 coil information if necessary
    QHS46_coils = load_QHS46_coils;
end

if length(coil_current_array) ~= 1
    error('The coir_current_arrray must have 1 element');
end

% Calculate the B field due to each individual coil set
%

% WF_ : winding factors
WF_Coil_1 = 1;

[Bx_Coil_1, By_Coil_1, Bz_Coil_1] = ...
    calc_B_BiotSavart(P_x, P_y, P_z, QHS46_coils.Coil_1, coil_current_array(1));


% Add them together.
Bx = WF_Coil_1 * Bx_Coil_1;

By = WF_Coil_1 * By_Coil_1;

Bz = WF_Coil_1 * Bz_Coil_1;


% The other useful values.
P_phi = atan2(P_y, P_x);
BR = Bx .* cos(P_phi) + By .* sin(P_phi);
BPhi = By .* cos(P_phi) - Bx .* sin(P_phi);
