function [Bx,By,Bz,BR,BPhi] = calc_b_LTX(P_x, P_y, P_z, coil_current_array)
%function [Bx,By,Bz,Br,Bphi] = calc_b_LTX(P_x, P_y, P_z, [TF OH RedUpper RedLower ...
% OrangeUpper etc. ])
% Input: Px, Py, Pz:  Cartesian coordinates of the observation point
%        coil_current_array:  An array of 14 elements that specify the coil
%        current in the following order:
%                [ToroidalField OhmicHeating Red_Upper Red_Lower Orange_U
%                Orange_L Yellow_U Yellow_L Green_U Green_L Blue_U Blue_L
%                Internal_U Internal_L]
%
% Output:  Bx, By, Bz:  Magnetic flux density (T/m) at observation point
%          Br, Bphi: Bx and By are used to find radial and toroidal
%          components.
%
% Note: Positive current gives counter-clockwise current in poloidal field
% coils and OH.  A negative current in the toroidal field coils produces
% a toroidal field that points in the clockwise direction when viewed from
% above.  This is the standard way that the LTX fields are operated.  A
% NEGATIVE CURRENT IN THE TOROIDAL FIELD COILS IS THE STANDARD OPERATION OF
% LTX.
%
% Calculates cartesian B field from given cartesian lab frame
% coordinates and coil current array.  If the observation point lies
% on a current filament (or extremely close to one) the contribution from
% that filament is ignored.
%
% The file LTX_coils.mat must be available.
%

% The details of the LTX coils are stored in this persistent array
persistent LTX_coils

if isempty(LTX_coils)  % load the LTX coil information if necessary
    LTX_coils = load_LTX_coils;
end

if length(coil_current_array) ~= 14
    error('The coir_current_arrray must have 14 elements');
end

% Calculate the B field due to each individual coil set
[Bx_TF, By_TF, Bz_TF] = ...
    calc_B_BiotSavart(P_x, P_y, P_z, LTX_coils.TF, coil_current_array(1));

[Bx_OH, By_OH, Bz_OH] = ...
    calc_B_BiotSavart(P_x, P_y, P_z, LTX_coils.OH_coil, coil_current_array(2));

[Bx_red_upper, By_red_upper, Bz_red_upper] = ...
    calc_B_BiotSavart(P_x, P_y, P_z, LTX_coils.red_coil_upper, coil_current_array(3));

[Bx_red_lower, By_red_lower, Bz_red_lower] = ...
    calc_B_BiotSavart(P_x, P_y, P_z, LTX_coils.red_coil_lower, coil_current_array(4));

[Bx_orange_upper, By_orange_upper, Bz_orange_upper] = ...
    calc_B_BiotSavart(P_x, P_y, P_z, LTX_coils.orange_coil_upper, coil_current_array(5));

[Bx_orange_lower, By_orange_lower, Bz_orange_lower] = ...
    calc_B_BiotSavart(P_x, P_y, P_z, LTX_coils.orange_coil_lower, coil_current_array(6));

[Bx_yellow_upper, By_yellow_upper, Bz_yellow_upper] = ...
    calc_B_BiotSavart(P_x, P_y, P_z, LTX_coils.yellow_coil_upper, coil_current_array(7));

[Bx_yellow_lower, By_yellow_lower, Bz_yellow_lower] = ...
    calc_B_BiotSavart(P_x, P_y, P_z, LTX_coils.yellow_coil_lower, coil_current_array(8));

[Bx_green_upper, By_green_upper, Bz_green_upper] = ...
    calc_B_BiotSavart(P_x, P_y, P_z, LTX_coils.green_coil_upper, coil_current_array(9));

[Bx_green_lower, By_green_lower, Bz_green_lower] = ...
    calc_B_BiotSavart(P_x, P_y, P_z, LTX_coils.green_coil_lower, coil_current_array(10));

[Bx_blue_upper, By_blue_upper, Bz_blue_upper] = ...
    calc_B_BiotSavart(P_x, P_y, P_z, LTX_coils.blue_coil_upper, coil_current_array(11));

[Bx_blue_lower, By_blue_lower, Bz_blue_lower] = ...
    calc_B_BiotSavart(P_x, P_y, P_z, LTX_coils.blue_coil_lower, coil_current_array(12));

[Bx_internal_upper, By_internal_upper, Bz_internal_upper] = ...
    calc_B_BiotSavart(P_x, P_y, P_z, LTX_coils.internal_coil_upper, coil_current_array(13));

[Bx_internal_lower, By_internal_lower, Bz_internal_lower] = ...
    calc_B_BiotSavart(P_x, P_y, P_z, LTX_coils.internal_coil_lower, coil_current_array(14));

% Add them together.
Bx = Bx_TF + Bx_OH + Bx_red_upper + Bx_red_lower + Bx_orange_upper + ...
    Bx_orange_lower + Bx_yellow_upper + Bx_yellow_lower + Bx_green_upper + ...
    Bx_green_lower + Bx_blue_upper  + Bx_blue_lower  + Bx_internal_upper + ...
    Bx_internal_lower;
By = By_TF + By_OH + By_red_upper + By_red_lower + By_orange_upper + ...
    By_orange_lower + By_yellow_upper + By_yellow_lower + By_green_upper + ...
    By_green_lower + By_blue_upper  + By_blue_lower  + By_internal_upper + ...
    By_internal_lower;
Bz = Bz_TF + Bz_OH + Bz_red_upper + Bz_red_lower + Bz_orange_upper + ...
    Bz_orange_lower + Bz_yellow_upper + Bz_yellow_lower + Bz_green_upper + ...
    Bz_green_lower + Bz_blue_upper  + Bz_blue_lower  + Bz_internal_upper + ...
    Bz_internal_lower;

P_phi = atan2(P_y, P_x);
BR = Bx .* cos(P_phi) + By .* sin(P_phi);
BPhi = By .* cos(P_phi) - Bx .* sin(P_phi);
