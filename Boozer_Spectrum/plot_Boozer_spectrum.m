function [] = compare_Boozer_and_Biot_Savart(fname_ext, surf_to_plot, ...
    modes_to_display)
%
% Usage:  
% compare_Boozer_and_Biot_Savart('boozmn_QHS_Rstart_1_513_6x8.nc', 6, [1:50])
%
% Compare Boozer/VMEC (xbooz_xform) output with field line following code
% It integrates the field line along increments of B * dl and then converts
% from Chi back to phi afterwards by dividing by the  boozer g factor.
%
% Must have read_boozer.m (part of stellopt matlab package on course site)
%
% Modified from previous versions by J. Schmitt and J. Talmadge.
% Intended for use during ECE-908
%
% Please report bugs to jcschmitt@wisc.edu

global current
global taper
current=10722; % QHS 1-Tesla on-axis
taper=0.*[1 1 1 -1 -1 -1]; % Standard QHS configuration is all zero.

relTol = 1e-10;
absTol = 1e-12;

%**************************************************************************
% Load up vmec xbooz_xform output.  
% Must have booz_xform (part of stellopt package) for details
boozer_data = read_boozer(fname_ext);
ns_b = boozer_data.ns; % number of surfs
mboz_b = boozer_data.mpol; % number of poloidal modes
nboz_b = boozer_data.ntor; % " " of toroidal modes
mnboz_b = boozer_data.mn_modes; % " " of total modes = (nboz_b * 2 + 1) * (mboz_b) - nboz_b
ixm_b = boozer_data.xm; % poloidal mode numbers
ixn_b = boozer_data.xn; % toroidal mode numbers
bmn_b = boozer_data.bmnc; % mode magnitudes (signed)
phi_b = boozer_data.phi; % enclosed toroidal flux (as a function of radius)
gmn_b = boozer_data.gmnc; % boozer jacobian
iota_b = boozer_data.iota; % boozer iota
rmnc_b = boozer_data.rmnc; % cos terms of r-expansion
zmns_b = boozer_data.zmns; % sin terms of z-expansion
bvco_b = boozer_data.bvco;  % The Boozer g factor 
buco_b = boozer_data.buco;  % The Boozer I factor

% Modify the settings below to change integration path lengths and path
% length spacing

% set up spacing in chi space by specifying roughly
% how far in toroidal phi you want to go. Also specify the number of
% data points you want on the field line

% Total number of divisions per field period
num_divisions_fp = 1001;

% Set the number of tooridal field periods to trace the field line
tor_periods = 1;

phi_extent = 2.*pi*tor_periods;
num_divisions = num_divisions_fp * tor_periods;

boozer_g = abs(bvco_b(surf_to_plot));
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
    max_of_mode(ii) = max( abs( bmn_b(ii, surf_to_plot) ));
    mode_label{ii} = [ '(' num2str(ixn_b(ii)) ',' num2str(ixm_b(ii)) ')' ];
end

[~, sorted_indices] = sort(max_of_mode, 'descend');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section is for making plots of |B| along a field line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%First calculate |B| along field line from reconstruction of spectrum
tor_angle_modB = linspace(0,phi_extent,num_divisions);
pol_angle_modB = iota_b(surf_to_plot) * (tor_angle_modB);
modB_recon_modB = zeros(size(pol_angle_modB));

disp('Using the following modes');
disp('Index |  m  |  n | bmn ')
for ii = (modes_to_display)
    disp([num2str(ii) ' | ' num2str(ixm_b(sorted_indices(ii))) ' | ' ...
        num2str(ixn_b(sorted_indices(ii))) ' | ' ...
        num2str(bmn_b(sorted_indices(ii), surf_to_plot))]);
end

for ii = (modes_to_display)
    angle_value_modB = (ixm_b(sorted_indices(ii)) * pol_angle_modB - ...
        ixn_b(sorted_indices(ii)) * tor_angle_modB);
    modB_recon_modB = modB_recon_modB + ...
        bmn_b(sorted_indices(ii), surf_to_plot) * cos(angle_value_modB);
end

disp('Completed calculating |B| along field line from Boozer spectrum');


%*****************************************************************
% First, follow B to determine coordinates, with boozer spacing***
%*****************************************************************
% Find a starting location: 
% A convenient location is at the Boxport (phi=0), midplane (z=0), where
% R = sum of cos-components at phi = 0;
rStart = sum(rmnc_b(:, surf_to_plot));
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
[Chif, Xspace] = ode113('chi_Boozer_derivatives', chispace, initial, ...
    int_options);

% Extract components from the solution matrix Xpace
rspace = Xspace(:,1);
phispace = Xspace(:,2);
zspace = Xspace(:,3);

numPointsAlongFieldLine = length(Chif);

% preallocate memory for B components
bx = zeros(1,numPointsAlongFieldLine);
by = bx; bz = bx;
modB_grid = bx;

% Now go back and find the magnitude of the field for each point on the
% field line
disp('Field line following complete. Calculating |B| along line')
for kk = 1:numPointsAlongFieldLine
    [bx(kk), by(kk), bz(kk), modB_grid(kk)] = ...
        calc_b_bs(rspace(kk), phispace(kk), zspace(kk), current, taper);
end

% Make some plots
figure;
box on; hold on
plot(tor_angle_modB, modB_recon_modB, 'b', 'LineWidth', 2)
plot(Chif/boozer_g, modB_grid,'r', 'Linewidth',2);

xlabel('Phi')
ylabel('|B|')

legend(['Boozer, with ' num2str(length(modes_to_display)) ' modes'], ...
    'Field line following using chi/boozerg')


