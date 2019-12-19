function [Bx,By,Bz,BR,BPhi] = calc_b(Coilset, P_x, P_y, P_z, coil_current_array)
%function [Bx,By,Bz,Br,Bphi] = calc_b(Coilset, P_x, P_y, P_z,
%coil_current_array)

% Inputs: Coilset: Name of coilset to use for Biot-Savart calcuations
%         Px, Py, Pz:  Cartesian coordinates of the observation point
%         coil_current_array:  An array of 1 or more numbers (coilset
%         specific)
%
% Output:  Bx, By, Bz:  Magnetic flux density (T/m) at observation point
%          Br, Bphi: Radial and toroidal components.
%
% Calculates cartesian B field from given cartesian lab frame
% coordinates and coil current array.  If the observation point lies
% on a current filament (or extremely close) the contribution from that
% filament is ignored.
%
% The file coilset = 'coilset_xyz.mat' must be available in the path
% somewhere
%
% Example:
% [bx, by, bz, br, bphi] = calc_b('coilset_wistell_a_004', xx, 0 ,0, 4.7e5 * 5984/4444);




% The details of the Coilset coils are stored in this persistent array
persistent Field_coils coilnames Coilset_Loaded

if isempty(Field_coils)  % load the coil information if necessary
    % No coilset loaded -> Reload coilset
    Reload_coilset = 1;
elseif strcmpi(Coilset_Loaded, Coilset)
    % Requested coilset is laoded -> do nothing
    Reload_coilset = 0;
else
    % Else -> Requestion coilset is not loaded -> Reload coilset
    Reload_coilset = 1;
end

if Reload_coilset
    [Field_coils, coilnames] = load_field_coils(Coilset);
    Coilset_Loaded = Coilset;
end

if length(coil_current_array) < 1
    error('The coil_current_arrray must have at least 1 element');
end

% Calculate the B field due to each individual coil in the coilset
%

num_coils = length(Field_coils.winding_factors);
% Easy case: One coil
if (num_coils == 1)
    coilname = coilnames{1};
    [Bx, By, Bz] = ...
        eval(['calc_B_BiotSavart(P_x, P_y, P_z, Field_coils.', coilname ', coil_current_array);']);
    Bx = Field_coils.winding_factors * Bx;
    By = Field_coils.winding_factors * By;
    Bz = Field_coils.winding_factors * Bz;
else
    % Initialize the variables to 0
    Bx_Coil = 0 * Field_coils.winding_factors;
    By_Coil = 0 * Field_coils.winding_factors;
    Bz_Coil = 0 * Field_coils.winding_factors;
    for ii = 1:num_coils
        coilname = Field_coils.coil_order{ii};
        [Bx_Coil(ii), By_Coil(ii), Bz_Coil(ii)] = ...
            eval(['calc_B_BiotSavart(P_x, P_y, P_z, Field_coils.', coilname ', coil_current_array(ii));']);
    end
    Bx = sum(Field_coils.winding_factors .* Bx_Coil);
    By = sum(Field_coils.winding_factors .* By_Coil);
    Bz = sum(Field_coils.winding_factors .* Bz_Coil);
end

% The other useful values.
P_phi = atan2(P_y, P_x);
BR = Bx .* cos(P_phi) + By .* sin(P_phi);
BPhi = By .* cos(P_phi) - Bx .* sin(P_phi);
