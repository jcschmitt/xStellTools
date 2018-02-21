function [fig_handle_1, fig_handle_2] = ...
    plot_Boozer_spectrum(Boozer_output_file, surface_to_plot, ...
    modes_to_display, PLOT_SPECTRUM, PLOT_SURFACE_2D, ...
    color_list, symbol_cells, line_cells, colormap_in)
% function [x_3D, y_3D, z_3D, modB_recon_3D] = ...
%     plot_boozer_spectrum(Boozer_output_file, surface_to_plot, ...
%     modes_to_display, PLOT_SPECTRUM, PLOT_SURFACE_2D, ...
%     color_list, symbol_cells, line_cells, colormap_in)
%
% Usage: plot_Boozer_spectrum('boozmn_qhs.nc', 99, 1:20, 1, 1);
%            or
%        plot_Boozer_spectrum('boozmn_qhs.nc', 75, 3:8, 1, 1);
%
% Inputs
%     Boozer_output_file: The full name of the Boozer output file.
%     'boozmn___.nc'
%     surf_to_plot: The index of the surface you wist to plot
%     modes_to_display: A 'mask' array of the modes you wish to display.
%     Modes are sorted from largest to smallest, and the index '1' is the
%     largest, index '2' is the next largest, etc.
%     color_list: A list of colors for the different spectrum. If more
%     modes are selected than colors, the colors are cylced.
%     symbol_cells: A cell of symbols for the different spectrum. If more
%     modes are selected than COLORS (yes, the check is on colors), the
%     symbols are cylced.
%     line_cells: A cell of line styles.
%     colormap_in: A preferred colormap. Default = 'jet'
%======================================================================


if nargin < 8
    % color_list =  'bkgrcmybkgrcmy';
    % symbol_list = '...xosdxxx....';
    disp('<----Using default color/marker/line scheme.')
    color_list =  'bbbbbbbbbbbbb';
    symbol_cells = {'none', 'none', 'none', 'none', 'none', 'none', ...
        'none', 'none', 'none', 'none', 'none', 'none', 'none'};
    line_cells = {'-', '-', '-', '-', '-', '-', '-', '-', '-', '-', ...
        '-', '-', '-'};
    colormap_in = 'jet';
end

% These variable control the generation of figures
if nargin < 5
    PLOT_SURFACE_2D = 1;
end
if nargin < 4
    PLOT_SPECTRUM = 1;
end

% This controls the visibility (on/off) of mode #'s on the figures.
SHOW_LABELS = 1; % 1=on  0=off

% This controls whether the x-axis is 'rho=sqrt(s)' or 's=normalized flux'
% 1 -> x-axis is 'rho'; 0 -> x-axis is 's'
PLOT_VS_R = 1;

if nargin < 3
    modes_to_display = 1:8;
end

if nargin < 2
    error('K----Don''''t be lazy.  Enter a surface index #');
end

boozer_data = read_Boozer_output(Boozer_output_file);
ns_b = boozer_data.ns_b; % number of surfs
disp(['<----Found ' num2str(ns_b) ' surfaces in ' Boozer_output_file]);

mnboz_b = boozer_data.mnboz_b; % " " of total modes = (nboz_b * 2 + 1) * (mboz_b) - nboz_b
ixm_b = boozer_data.xm_b; % poloidal mode numbers
ixn_b = boozer_data.xn_b; % toroidal mode numbers
bmnc_b = boozer_data.bmnc_b; % mode magnitudes (signed)
phi_b = boozer_data.phi_b; % enclosed toroidal flux (as a function of radius)
iota_b = boozer_data.iota_b; % boozer iota

% [ixn_b ixm_b bmnc_b(1,:)'] would give list of tor, pol mode number and
% magnitude for surface 1
% [ixn_b ixm_b bmnc_b(2,:)'] "  " for surface 2.  etc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Re-sort the amplitude data into nice arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_of_mode = zeros(1,mnboz_b);
mode_label = cell(1,mnboz_b);

for ii = 1:mnboz_b
    % to use the 'global' max components
    %     max_of_mode(ii) = max( abs( bmnc_b(:, ii) ));
    % to use the 'surface' maximum components
    max_of_mode(ii) = max( abs( bmnc_b(surface_to_plot, ii) ));
    
    mode_label{ii} = [ '(' num2str(ixn_b(ii)) ',' num2str(ixm_b(ii)) ')' ];
end

[~, sorted_indices] = sort(max_of_mode, 'descend');

if PLOT_SPECTRUM
    % Make plots of the radial profile of the largest modes (sorted
    % according to the values of the mode on the 'surface_to_plot'
    if (modes_to_display ~= 0)
        disp(['Displaying modes with indices: [' num2str(sorted_indices(modes_to_display)) '] ']);
        disp(['Displaying modes with m nums:  [' num2str(ixm_b(sorted_indices(modes_to_display))') '] ']);
        disp(['Displaying modes with n nums:  [' num2str(ixn_b(sorted_indices(modes_to_display))') '] ']);
        fig_handle_1 = figure;box on;
        legend_text = [];
        
        x_values = phi_b(1:end)/phi_b(end);
        if PLOT_VS_R
            % rho = sqrt(s)
            x_values = sqrt(x_values);
            XLABEL = '\rho';
        else
            XLABEL = 's = \Psi/\Psi_{LCFS}';
        end
        %     XAXIS = ns_b;
        
        for ii = modes_to_display
            jj = mod(ii-1, length(color_list)) + 1;
            plot(x_values(2:end), bmnc_b(:, sorted_indices(ii)), ...
                'Color', color_list(jj), ...
                'Marker', symbol_cells{jj}, 'LineStyle', line_cells{jj});
            hold on
            if SHOW_LABELS
                text( x_values(end)+.01*(x_values(end) - x_values(1)), bmnc_b(end, sorted_indices(ii)), ...
                    mode_label(sorted_indices(ii)), 'Color', color_list(jj));
            end
            legend_text = [legend_text ; mode_label(sorted_indices(ii))];
        end
        
        xlabel(XLABEL);
        ylabel('Spectrum amplitude');
        cur_axis = axis;
        axis([cur_axis(1) cur_axis(2)*1.015 cur_axis(3) cur_axis(4)]);
        title(strrep([Boozer_output_file], '_', '\_'));
        legend(legend_text);
        
        
        for ii = modes_to_display
            disp(['Mode with index: [' num2str(sorted_indices(ii)) '] ']);
            disp(['Mode with (n,m) =  (' num2str(ixn_b(sorted_indices(ii))) ', ' ...
                num2str(ixm_b(sorted_indices(ii))') ') ']);
            for jj = surface_to_plot
                disp([' xaxis val = ' num2str(x_values(jj))]);
                disp([' b_nm = ' num2str(bmnc_b(jj, sorted_indices(ii))) ]);
            end
        end
    end
end


if PLOT_SURFACE_2D
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  this section is for making the |B| plot on the flux surface
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    num_divisions = 400;
    num_contour_divisions = 10;
    tor_angle = linspace(0,2*pi, num_divisions);
    pol_angle = linspace(0,2*pi, num_divisions);
    
    [tor_matrix, pol_matrix] = meshgrid(tor_angle, pol_angle);
    
    modB_recon_3D = zeros(size(pol_matrix));
    
    for ii = (modes_to_display)
        angle_value = (ixm_b(sorted_indices(ii)) * pol_matrix - ixn_b(sorted_indices(ii)) * tor_matrix);
        modB_recon_3D = modB_recon_3D + bmnc_b(surface_to_plot, sorted_indices(ii)) * cos(angle_value);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % this section sets up the B- lines for the contour plot
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tor_angle = linspace(0,2*pi, num_divisions);
    pol_angle = linspace(0,2*pi, num_divisions);
    
    num_field_lines = 9;
    B_0 = -2*pi; B_end = 2*pi;
    BLine_0 = linspace(B_0, B_end, num_field_lines);
    
    % BLine_0 = 0;
    BLines_tor = tor_angle;
    BLines_pol = zeros(num_field_lines, num_divisions);
    
    for ii = 1:length(BLine_0)
        BLines_pol(ii,:) = iota_b(surface_to_plot) * (BLines_tor - BLine_0(ii));
    end
    
    fig_handle_2 = figure;
    contourf(tor_angle, pol_angle, modB_recon_3D, num_contour_divisions);
    colorbar;
    colormap(colormap_in);
    hold on;box on;
    for ii = 1:length(BLine_0)
        plot3(BLines_tor, BLines_pol(ii, :),2*ones(size(BLines_tor)), 'k:', 'LineWidth', 2);
    end
    xlabel('\zeta_{Boozer}','FontSize', 14);
    ylabel('\theta_{Boozer}','FontSize', 14);
    
    view(2)
    set(gca,'FontSize', 14);
    title(['|B| on surf \rho = ' num2str(sqrt(phi_b(surface_to_plot)/ phi_b(end)))])
    hold on;
    axis([0 2*pi 0 2*pi]);
end


