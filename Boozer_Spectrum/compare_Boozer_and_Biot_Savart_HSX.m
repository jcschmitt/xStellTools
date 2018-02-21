function [figure_handle] = ...
    compare_Boozer_and_Biot_Savart_HSX(Boozer_output_file, ...
    surface_to_plot, modes_to_include, current, taper)
%
% Usage:  
% compare_Boozer_and_Biot_Savart_HSX(...
%     'boozmn_QHS_Rstart_1_513_6x8.nc', 6, [1:50], 10722, [0 0 0 0 0 0])
%
% Compare Boozer/VMEC (xbooz_xform) output with field line following code.
% This integrates the field line along increments of |B|*dl.
% A conversion from chi to phi is done by dividing by the boozer g-factor.
%
% Must have read_boozer.m (part of stellopt matlab package - source should
% go here) 
%
% Modified from previous versions by J. Schmitt and J. Talmadge.
%
% Please report bugs to jcschmitt@wisc.edu

relTol = 1e-10;
absTol = 1e-12;

%**************************************************************************
% Load up xbooz_xform output.  
boozer_data = read_Boozer_output(Boozer_output_file);
ns_b = boozer_data.ns_b; % number of surfs
disp(['<----Found ' num2str(ns_b) ' surfaces in ' Boozer_output_file]);

mnboz_b = boozer_data.mnboz_b; % " " of total modes = (nboz_b * 2 + 1) * (mboz_b) - nboz_b
if strcmp(modes_to_include, 'all')
    modes_to_include = 1:mnboz_b;
end

ixm_b = boozer_data.xm_b; % poloidal mode numbers
ixn_b = boozer_data.xn_b; % toroidal mode numbers
bmnc_b = boozer_data.bmnc_b; % mode magnitudes (signed)
iota_b = boozer_data.iota_b; % boozer iota
rmnc_b = boozer_data.rmnc_b; % cos terms of r-expansion
bvco_b = boozer_data.bvco_b;  % The Boozer g factor - depends on toroidal flux

% Total number of divisions per field period
num_divisions_fp = 1001;

% Set the number of tooridal field periods to trace the field line
tor_periods = 1;

phi_extent = 2.*pi*tor_periods;
num_divisions = num_divisions_fp * tor_periods;

boozer_g = abs(bvco_b(surface_to_plot));
chi_extent = phi_extent * boozer_g;

% [ixn_b' ixm_b' bmn_b(:,1)] would give list of tor, pol mode number and
% magnitude for surface 1
% [ixn_b' ixm_b' bmn_b(:,2)] "  " for surface 2.  etc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Re-sort the amplitude data into nice arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_of_mode = zeros(1,mnboz_b);
mode_label = cell(1,mnboz_b);

for ii = 1:mnboz_b
    % Sorted amplitude of modes according to highest value over all
    % surface.
    %    max_of_mode(ii) = max( abs( bmn_b(:, ii) ));
    % Sort modes only over the surface you want to look at.
    % max_of_mode(ii) = max( abs( bmn_b(surf_to_plot, ii) )); %old way
    max_of_mode(ii) = max( abs( bmnc_b(surface_to_plot, ii) ));
    mode_label{ii} = [ '(' num2str(ixn_b(ii)) ',' num2str(ixm_b(ii)) ')' ];
end

[~, sorted_indices] = sort(max_of_mode, 'descend');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section is for making plots of |B| along a field line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%First calculate |B| along field line from reconstruction of spectrum
tor_angle_modB = linspace(0, phi_extent, num_divisions);
pol_angle_modB = iota_b(surface_to_plot) * (tor_angle_modB);
modB_recon_modB = zeros(size(pol_angle_modB));

disp('Using the following modes');
disp('Index |  m  |  n | bmn ')
for ii = (modes_to_include)
    disp([num2str(ii) ' | ' num2str(ixm_b(sorted_indices(ii))) ' | ' ...
        num2str(ixn_b(sorted_indices(ii))) ' | ' ...
        num2str(bmnc_b(surface_to_plot, sorted_indices(ii)))]);
end

for ii = (modes_to_include)
    angle_value_modB = (ixm_b(sorted_indices(ii)) * pol_angle_modB - ...
        ixn_b(sorted_indices(ii)) * tor_angle_modB);
    modB_recon_modB = modB_recon_modB + ...
        bmnc_b(surface_to_plot, sorted_indices(ii)) * cos(angle_value_modB);
end

disp('Completed calculating |B| along field line from Boozer spectrum');


%*****************************************************************
% First, follow B to determine coordinates, with boozer spacing***
%*****************************************************************
% Find a starting location: 
% A convenient location is at the Boxport (phi=0), midplane (z=0), where
% R = sum of cos-components at phi = 0;
rStart = sum(rmnc_b(surface_to_plot, :));
phiStart = 0;
zStart = 0;

% Set up the initial condition vector for the ODE solver.
initial = [rStart phiStart zStart];   
disp(['Starting coordinates (R, Phi, Z) = ' num2str(initial)]);

% Set up integration tolerances. Adjusting these will affect the
% accuracy/speed of the calculation
int_options = odeset('AbsTol', absTol, 'RelTol', relTol, 'NormControl', 'off');

% Set up chispace vector
chispace = linspace(0, chi_extent, num_divisions);

% Follow the vacuum field lines using the Boozer derivaties
% This creates a function that is used by the ode solver. I need to pass in
% the current and taper variables, and this is the easiest way in MATLAB.
myode_fcn = @(my_chi, my_X) Boozer_chi_derivatives_HSX(my_chi, my_X, ...
    current, taper);

[chi_out, X_out] = ode113(myode_fcn, chispace, initial, ...
    int_options);
%chif = 0;X = [rStart phiStart zStart];
% Extract components from the solution matrix Xpace
R_out = X_out(:,1);
phi_out = X_out(:,2);
Z_out = X_out(:,3);

numPointsAlongFieldLine = length(chi_out);

% preallocate memory for B components
bx = zeros(1, numPointsAlongFieldLine);
by = zeros(1, numPointsAlongFieldLine);
bz = zeros(1, numPointsAlongFieldLine);

% Now go back and find the magnitude of the field for each point on the
% field line
disp('Field line following complete. Calculating |B| along line')
for kk = 1:numPointsAlongFieldLine
    [bx(kk), by(kk), bz(kk)] = ...
        calc_b_HSX_RPhiZ(R_out(kk), phi_out(kk), Z_out(kk), current, taper);
end
modB_grid = sqrt(bx.^2 + by.^2 + bz.^2);

% Make some plots
figure_handle = figure;
box on; hold on;
plot(tor_angle_modB, modB_recon_modB, 'b', 'LineWidth', 2)
plot(chi_out / boozer_g, modB_grid, 'r', 'Linewidth',2);

xlabel('Phi')
ylabel('|B|')

legend(['Boozer, with ' num2str(length(modes_to_include)) ' modes'], ...
    'Field line following using chi/boozerg')
title([strrep(Boozer_output_file, '_', '\_') '  Surface #' ...
    num2str(surface_to_plot)]);


