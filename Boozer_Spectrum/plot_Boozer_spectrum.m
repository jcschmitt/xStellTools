function [fig_handle_1, fig_handle_2] = ...
    plot_Boozer_spectrum(Boozer_output_file, surface_to_plot, ...
    modes_to_display, sorting_surface, PLOT_SPECTRUM, PLOT_SURFACE_2D, ...
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


if nargin < 10
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
if nargin < 6
    PLOT_SURFACE_2D = 1;
end
if nargin < 5
    PLOT_SPECTRUM = 1;
end
if nargin < 4
    sorting_surface = 1;
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
nfp_b = boozer_data.nfp_b; % number of surfs
disp(['<----Found ' num2str(ns_b) ' surfaces in ' Boozer_output_file]);

lasym_b = boozer_data.lasym; % " " of total modes = (nboz_b * 2 + 1) * (mboz_b) - nboz_b

mnboz_b = boozer_data.mnboz_b; % " " of total modes = (nboz_b * 2 + 1) * (mboz_b) - nboz_b
ixm_b = boozer_data.xm_b; % poloidal mode numbers
ixn_b = boozer_data.xn_b; % toroidal mode numbers
bmnc_b = boozer_data.bmnc_b; % mode magnitudes (signed) stellarator symmetric terms
phi_b = boozer_data.phi_b; % enclosed toroidal flux (as a function of radius)
iota_b = boozer_data.iota_b; % boozer iota

if lasym_b
    bmns_b = boozer_data.bmns_b; % mode magnitudes (signed) stellarator symmetric terms
end

% [ixn_b ixm_b bmnc_b(1,:)'] would give list of tor, pol mode number and
% magnitude for surface 1
% [ixn_b ixm_b bmnc_b(2,:)'] "  " for surface 2.  etc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Re-sort the amplitude data into nice arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_max_stelsym = -inf;
max_max_nonstelsym = -inf;

max_of_mode_stelsym = zeros(1,mnboz_b);
mode_label_stelsym = cell(1,mnboz_b);
mode_label_nonstelsym = cell(1,mnboz_b);

% stellarator symmetric modes
for ii = 1:mnboz_b
    % to use the 'global' max components
    max_of_mode_stelsym(ii) = max( abs( bmnc_b(:, ii) ));
    if lasym_b
        max_of_mode_nonstelsym(ii) = max( abs( bmns_b(:, ii) ));
    end
    % to use the 'surface' maximum components instead of the global max
    % max_of_mode(ii) = max( abs( bmnc_b(surface_to_plot, ii) ));
    % max_of_mode_stelsym(ii) = max( abs( bmnc_b(surface_to_plot, ii) ));
    % max_of_mode_nonstelsym(ii) = max( abs( bmns_b(surface_to_plot, ii) ));
    
    % to use inner surface (I think)
    % max_max = max(max_max, abs(bmnc_b(1,ii)));
    max_max_stelsym = max(max_max_stelsym, abs(bmnc_b(1,ii)));
    if lasym_b
        max_max_nonstelsym = max(max_max_nonstelsym, abs(bmns_b(1,ii)));
    end  
    
    mode_label_stelsym{ii} = [ '(' num2str(ixn_b(ii)) ',' num2str(ixm_b(ii)) ')' ];
    mode_label_nonstelsym{ii} = [ '(' num2str(ixn_b(ii)) ',' num2str(ixm_b(ii)) ')' ];
    %if (ixn_b(ii) == 48 && ixm_b(ii) ==0)
    %    keyboard
    %end
end

if (max_max_nonstelsym > max_max_stelsym)
    error('<----How did this happen?');
end

[~, sorted_indices_stelsym] = sort(max_of_mode_stelsym, 'descend');
    if lasym_b

[~, sorted_indices_nonstelsym] = sort(max_of_mode_nonstelsym, 'descend');
    end
    
if PLOT_SPECTRUM
    % Make plots of the radial profile of the largest modes (sorted
    % according to the values of the mode on the 'surface_to_plot'
    if (modes_to_display ~= 0)
        disp(['Displaying modes with indices: [' num2str(sorted_indices_stelsym(modes_to_display)) '] ']);
        disp(['Displaying modes with m nums:  [' num2str(ixm_b(sorted_indices_stelsym(modes_to_display))') '] ']);
        disp(['Displaying modes with n nums:  [' num2str(ixn_b(sorted_indices_stelsym(modes_to_display))') '] ']);
        
        x_values = phi_b(1:end)/phi_b(end);
        if PLOT_VS_R
            % rho = sqrt(s)
            x_values = sqrt(x_values);
            XLABEL = '\rho';
        else
            XLABEL = 's = \Psi/\Psi_{LCFS}';
        end
        
        % stellarator symmetric components
        fig_handle_1 = figure;box on;
        legend_text = [];
        
        %     XAXIS = ns_b;
        for ii = modes_to_display
            jj = mod(ii-1, length(color_list)) + 1;
            plot(x_values(2:end), bmnc_b(:, sorted_indices_stelsym(ii)) / max_max_stelsym, ...
                'Color', color_list(jj), ...
                'Marker', symbol_cells{jj}, 'LineStyle', line_cells{jj});
            hold on;grid on
            
            if SHOW_LABELS
                text( x_values(end)+.01*(x_values(end) - x_values(1)), bmnc_b(end, sorted_indices_stelsym(ii)) / max_max_stelsym, ...
                    mode_label_stelsym(sorted_indices_stelsym(ii)), 'Color', color_list(jj));
            end
            legend_text = [legend_text ; mode_label_stelsym(sorted_indices_stelsym(ii))];
        end
        
        xlabel(XLABEL);
        ylabel('Symmetric Spectrum Amplitude');
        
        cur_axis = axis;
        axis([cur_axis(1) cur_axis(2)*1.015 cur_axis(3) cur_axis(4)]);
        title(strrep([Boozer_output_file], '_', '\_'));
        legend(legend_text);
        
        
        % non-stellarator symmetric components
        if lasym_b
            fig_handle_1b = figure;box on;
            legend_text_b = [];
            
            for ii = modes_to_display
                jj = mod(ii-1, length(color_list)) + 1;
                plot(x_values(2:end), bmns_b(:, sorted_indices_nonstelsym(ii)) / max_max_stelsym, ...
                    'Color', color_list(jj), ...
                    'Marker', symbol_cells{jj}, 'LineStyle', line_cells{jj});
                hold on;grid on
                
                if SHOW_LABELS
                    text( x_values(end)+.01*(x_values(end) - x_values(1)), bmns_b(end, sorted_indices_nonstelsym(ii)) / max_max_stelsym, ...
                        mode_label_nonstelsym(sorted_indices_nonstelsym(ii)), 'Color', color_list(jj));
                end
                legend_text_b = [legend_text_b ; mode_label_nonstelsym(sorted_indices_nonstelsym(ii))];
            end
            
            xlabel(XLABEL);
            ylabel('Non-Symm Spectrum Amplitude');
            cur_axis = axis;
            axis([cur_axis(1) cur_axis(2)*1.015 cur_axis(3) cur_axis(4)]);
            title(strrep([Boozer_output_file], '_', '\_'));
            legend(legend_text_b);
        end
        
        for ii = modes_to_display
            disp(['ii = ' num2str(ii)])
            disp(['StelSym Mode with index: [' num2str(sorted_indices_stelsym(ii)) '] ']);
            disp(['StelSym Mode with (n,m) =  (' num2str(ixn_b(sorted_indices_stelsym(ii))) ', ' ...
                num2str(ixm_b(sorted_indices_stelsym(ii))') ') ']);
            for jj = surface_to_plot
                disp([' xaxis val = ' num2str(x_values(jj))]);
                disp([' b_nmc = ' num2str(bmnc_b(jj, sorted_indices_stelsym(ii))) ]);
            end
        end
        if lasym_b
                for ii = modes_to_display
                disp(['ii = ' num2str(ii)])
                disp(['NonStelSym Mode with index: [' num2str(sorted_indices_nonstelsym(ii)) '] ']);
                disp(['NonStelSym Mode with (n,m) =  (' num2str(ixn_b(sorted_indices_nonstelsym(ii))) ', ' ...
                    num2str(ixm_b(sorted_indices_nonstelsym(ii))') ') ']);
                for jj = surface_to_plot
                    disp([' xaxis val = ' num2str(x_values(jj))]);
                    disp([' b_nms = ' num2str(bmns_b(jj, sorted_indices_nonstelsym(ii))) ]);
                end
            end
        end
        
    end
end
try
%    make_my_plot_pretty3;
catch
end

if PLOT_SURFACE_2D
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  this section is for making the |B| plot on the flux surface
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    num_divisions = 400;
    num_contour_divisions = 10;
    %tor_angle = linspace(0,2*pi, num_divisions);
    %pol_angle = linspace(0,2*pi, num_divisions);
    tor_angle = linspace(0,2*pi/nfp_b, num_divisions);
    pol_angle = linspace(0,2*pi, num_divisions);
    
    [tor_matrix, pol_matrix] = meshgrid(tor_angle, pol_angle);
    
    modB_recon_3D = zeros(size(pol_matrix));
    
    if ~lasym_b
        for ii = (modes_to_display)
            angle_value = (ixm_b(sorted_indices_stelsym(ii)) * pol_matrix - ixn_b(sorted_indices_stelsym(ii)) * tor_matrix);
            modB_recon_3D = modB_recon_3D + bmnc_b(surface_to_plot, sorted_indices_stelsym(ii)) * cos(angle_value);
        end
    else
        for ii = (modes_to_display)
            angle_value = (ixm_b(sorted_indices_stelsym(ii)) * pol_matrix - ixn_b(sorted_indices_stelsym(ii)) * tor_matrix);
            modB_recon_3D = modB_recon_3D + ...
                bmnc_b(surface_to_plot, sorted_indices_stelsym(ii)) * ...
                cos(angle_value) + ...
                bmns_b(surface_to_plot, sorted_indices_nonstelsym(ii)) * ...
                sin(angle_value);
        end
        
    end
    
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % this section sets up the B- lines for the contour plot
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tor_angle = linspace(0,2*pi/double(nfp_b), num_divisions);
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
        axis([0 2*pi/nfp_b 0 2*pi]);
    end
    
    
