function do_the_descur_stuff(coilsetID, coilCurrents, configuration_name, title_str, ...
    Rlaunch, num_field_periods, nphi_step1, ntheta_step1, num_launch_pts, ...
    nphi_step2, ntheta_step2, items_to_do)

if nargin < 12
    items_to_do = 1:num_launch_pts
end

%num_launch_pts = length(Rlaunch);

% make filenames for LCFS and more
fluxout_filename = ['Flux_' configuration_name];

for ii = items_to_do
    LCFS_filename{ii} = ['LCFS_' configuration_name '_Rlaunch_' num2str(round(1000*Rlaunch(ii))) 'mm.mat'];
    LCFS_pathname{ii} = './';
    outcurve_filename{ii} = ['outcurve_Rlaunch_' num2str(round(1000*Rlaunch(ii))) 'mm'];
    if (ii == 1)
        legend_str{ii} = ['Rlaunch = ' num2str(round(1000*Rlaunch(ii))) 'mm'];
    else
        legend_str{ii} = ['          ' num2str(round(1000*Rlaunch(ii))) 'mm'];
    end
end

% Make poincare plots near edge.
if 1
    
    parfor ii = items_to_do
        [allsurf_x{ii}, allsurf_y{ii}, allsurf_z{ii}, allsurf_phi{ii}] = ...
            calc_descur_RZ_points(coilsetID, Rlaunch(ii), nphi_step1, ntheta_step1, coilCurrents);
    end
    for ii = items_to_do
        surf_x = allsurf_x{ii};
        surf_y = allsurf_y{ii};
        surf_z = allsurf_z{ii};
        surf_phi = allsurf_phi{ii};
        save(LCFS_filename{ii}, 'Rlaunch', 'coilCurrents', ...
            'coilsetID', 'configuration_name', 'nphi_step1', 'ntheta_step1', ...
            'surf_x', 'surf_y', 'surf_z', 'surf_phi');
    end
end

surf_x = 0; surf_y = 0; surf_z = 0;
for ii = items_to_do
    load(LCFS_filename{ii})
    surf_x_flux{ii} = surf_x;
    surf_y_flux{ii} = surf_y;
    surf_z_flux{ii} = surf_z;
end

% Calculate toroidal flux (this should be pretty quick)
if 1
    tor_flux = 0 * Rlaunch;
    for ii = items_to_do
       
    %parfor ii = items_to_do
        disp(['<----Starting index #' num2str(ii)]);
        R = sqrt(surf_x_flux{ii}(1:nphi_step2:end).^2 + surf_y_flux{ii}(1:nphi_step2:end).^2);
        Z = surf_z_flux{ii}(1:nphi_step2:end);
        tor_flux(ii) = calculate_toroidal_flux(coilsetID, coilCurrents, R, Z, 1);
        save(fluxout_filename, 'Rlaunch', 'tor_flux');
    end
    
    save(fluxout_filename, 'Rlaunch', 'tor_flux');
end

if 1
    % Convert poincare plots to rzdata input for descur
    for ii = items_to_do
        
        % returns the rz and descur command names
        disp(['<----Starting index #' num2str(ii)]);
        [outfilename{ii}, descurcmdfilename{ii}] = ...
            make_descur_file(LCFS_filename{ii}, LCFS_pathname{ii}, ...
            num_field_periods, nphi_step2, ntheta_step2);
    end
    figure(401);
    legend(legend_str);
    grid on; axis equal
    axis equal
    axis tight
    title(title_str)
    make_my_plot_pretty3
    set(gcf, 'Position', [200 200 750 550])
    
end

if 1
    % Run descur and move files
    for ii = items_to_do
        descur_path = '~/src/stellinstall/trunk/build/bin/xdescur'
        %descur_path = '~/src/stellinstall/trunk/build/bin/xdescur'
        cmd1 = [descur_path ' < ' descurcmdfilename{ii}];
        [status, cmdout] = system(cmd1, '-echo');
        cmd2 = ['mv outcurve ' outcurve_filename{ii}];
        [status, cmdout] = system(cmd2, '-echo');
    end
end

if 1
    % Check descur
    for ii = items_to_do
        check_descur_output(LCFS_filename{ii}, outcurve_filename{ii}, ...
            num_field_periods)
    end
end
