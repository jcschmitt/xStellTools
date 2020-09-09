function do_the_descur_stuff(coilsetID, coilCurrents, configuration_name, title_str, ...
    Rlaunch, Zlaunch, num_field_periods, nphi_step1, ntheta_step1, num_launch_pts, ...
    nphi_step2, ntheta_step2, items_to_do)

if nargin < 13
    disp('<---Specific indices were not given. All launch point included in scan.');
    items_to_do = 1:num_launch_pts
else
    disp(['<----Adjusting num_launch_pts from ', num2str(num_launch_pts),  ...
        ' to ',  num2str(length(items_to_do))]);
    num_launch_pts = length(items_to_do);
end

%num_launch_pts = length(Rlaunch);

% make filenames for LCFS and more
fluxout_filename = ['Flux_' configuration_name];
iotaout_filename = ['Iota_' configuration_name];

for ii = items_to_do
    LCFS_filename{ii} = ['LCFS_' configuration_name '_Rlaunch_' num2str(round(1000*Rlaunch(ii))) 'mm.mat'];
    LCFS_pathname{ii} = './';
    outcurve_filename{ii} = ['outcurve_Rlaunch_' num2str(round(1000*Rlaunch(ii))) 'mm'];
    sortedataout_filename{ii} =['sorted_data_out_Rlaunch_' num2str(round(1000*Rlaunch(ii))) 'mm'];
    if (ii == 1)
        legend_str{ii} = ['Rlaunch = ' num2str(round(1000*Rlaunch(ii))) 'mm'];
    else
        legend_str{ii} = ['          ' num2str(round(1000*Rlaunch(ii))) 'mm'];
    end
end

% Make poincare plots near edge.
if 1
    
    parfor ii = items_to_do
         poincare_help(coilsetID, Rlaunch(ii), Zlaunch(ii), nphi_step1, ntheta_step1, coilCurrents, LCFS_filename{ii}, configuration_name)
    end
end

surf_x = 0; surf_y = 0; surf_z = 0;surf_phi = 0;
for ii = items_to_do
    load(LCFS_filename{ii})
    surf_x_flux{ii} = surf_x;
    surf_y_flux{ii} = surf_y;
    surf_z_flux{ii} = surf_z;
    surf_phi_flux{ii} = surf_phi;
end


% Calculate iota, save it,
if 0
    
    iota = calculate_iota_nonstelsym(surf_x_flux, surf_y_flux, surf_z_flux, surf_phi, nphi_step1, num_launch_pts);
    save(iotaout_filename, 'Rlaunch', 'Zlaunch', 'iota');
end

if 0
    % plot iota
    iota_data = load(iotaout_filename, 'Rlaunch', 'Zlaunch', 'iota');
    figure; box on; hold on;
    plot(iota_data.Rlaunch(items_to_do), iota_data.iota(items_to_do), 'o--');
end

% Calculate toroidal flux (this should be pretty quick)
if 0
    tor_flux = 0 * Rlaunch;
    for ii = items_to_do
        
        %parfor ii = items_to_do
        disp(['<----Starting index #' num2str(ii)]);
        R = sqrt(surf_x_flux{ii}(1:nphi_step2:end).^2 + surf_y_flux{ii}(1:nphi_step2:end).^2);
        Z = surf_z_flux{ii}(1:nphi_step2:end);
        tor_flux(ii) = calculate_toroidal_flux(coilsetID, coilCurrents, R, Z, 1);
        save(fluxout_filename, 'Rlaunch', 'Zlaunch', 'tor_flux');
    end
    
    save(fluxout_filename, 'Rlaunch', 'tor_flux');
end

if 0
    % Convert poincare plots to rzdata input for descur
    for ii = items_to_do
        
        % returns the rz and descur command names
        disp(['<----Starting index #' num2str(ii)]);
        surf_RZ_data{ii} = get_descur_data(LCFS_filename{ii}, LCFS_pathname{ii}, ...
            num_field_periods, nphi_step2, ntheta_step2);
        
        [outfilename{ii}, descurcmdfilename{ii}] = ...
            make_descur_file(LCFS_filename{ii}, LCFS_pathname{ii}, ...
            num_field_periods, nphi_step2, ntheta_step2);
        %close all
        pause(1)
    end
    %figure(401);
    %legend(legend_str);
    %%grid on; axis equal
    %axis equal
    %axis tight
    %title(title_str)
    %make_my_plot_pretty3
    %set(gcf, 'Position', [200 200 750 550])
    
    % end
    %
    if 0
        % Run descur and move files
        for ii = items_to_do
            descur_path = '~/src/stellinstall/trunk/build/bin/xdescur'
            %descur_path = '~/src/stellinstall/trunk/build/bin/xdescur'
            cmd1 = [descur_path ' < ' descurcmdfilename{ii}];
            [status, cmdout] = system(cmd1, '-echo');
            cmd2 = ['mv outcurve ' outcurve_filename{ii}];
            [status, cmdout] = system(cmd2, '-echo');
            cmd2 = ['mv sorted_data_out ' sortedataout_filename{ii}];
            [status, cmdout] = system(cmd2, '-echo');
        end
    end
    
    if 0
        % do my own 2-d fourier fitting
        for ii = items_to_do
            
            surf_RZ_data{ii};
            descur_matlab(surf_RZ_data{ii});
            keyboard;
        end
    end
    
end


if 0
    % Check descur
    for ii = items_to_do
        check_descur_output(LCFS_filename{ii}, outcurve_filename{ii}, ...
            num_field_periods)
    end
end

function output = poincare_help(coilsetID, Rlaunch, Zlaunch, nphi_step1, ntheta_step1, coilCurrents, LCFS_filename, configuration_name)
 
        [surf_x, surf_y, surf_z, surf_phi] = ...
            calc_descur_RZ_points_RZ(coilsetID, Rlaunch, Zlaunch, nphi_step1, ntheta_step1, coilCurrents);
        save(LCFS_filename, 'Rlaunch', 'Zlaunch', 'coilCurrents', ...
            'coilsetID', 'configuration_name', 'nphi_step1', 'ntheta_step1', ...
            'surf_x', 'surf_y', 'surf_z', 'surf_phi');

