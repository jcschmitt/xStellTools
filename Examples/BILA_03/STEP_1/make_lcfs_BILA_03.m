function make_lcfs_BILA_03

% definitions
coilsetID = 'coilset_bila_03';
configuration_name = 'bila_03';
%coilCurrents = 497292.;
coilCurrents = 1.;
title_str='BILA 03';

Rlaunch = 2.529 + 0.002*[0:13];
num_launch_pts = length(Rlaunch);

num_field_periods = 5;
nphi_step1 = 300;
ntheta_step1 = 200;
nphi_step2 = 60;
ntheta_step2 = 1000;


%=====

% make filenames for LCFS and more
fluxout_filename = ['Flux_' configuration_name];

for ii = 1:num_launch_pts
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
if 0
    for ii = 1:num_launch_pts
        this_R = Rlaunch(ii);
        [surf_x, surf_y, surf_z, surf_phi] = ...
            calc_descur_RZ_points(coilsetID, this_R, nphi_step1, ntheta_step1, coilCurrents);
        save(LCFS_filename{ii}, 'Rlaunch', 'coilCurrents', ...
            'coilsetID', 'configuration_name', 'nphi_step1', 'ntheta_step1', ...
            'surf_x', 'surf_y', 'surf_z', 'surf_phi', 'this_R');
    end
end

surf_x = 0; surf_y = 0; surf_z = 0;
for ii = 1:num_launch_pts
    load(LCFS_filename{ii})
    surf_x_flux{ii} = surf_x;
    surf_y_flux{ii} = surf_y;
    surf_z_flux{ii} = surf_z;
end

% Calculate toroidal flux
if 0
    tor_flux = 0 * Rlaunch;
    for ii = 1:num_launch_pts
        R = sqrt(surf_x_flux{ii}(1:nphi_step2:end).^2 + surf_y_flux{ii}(1:nphi_step2:end).^2);
        Z = surf_z_flux{ii}(1:nphi_step2:end);
        tor_flux(ii) = calc_toroidal_flux(coilsetID, coilCurrents, R, Z, 1);
    end
    
    save(fluxout_filename, 'Rlaunch', 'tor_flux');
end

if 0
    % Convert poincare plots to rzdata input for descur
    for ii = 1:num_launch_pts
        % returns the rz and descur command names
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

if 0
    % Run descur and move files
    for ii = 1:num_launch_pts
        cmd1 = ['~/src/stellinstall/trunk/build/bin/xdescur < ' ...
            descurcmdfilename{ii}];
        [status, cmdout] = system(cmd1, '-echo');
        cmd2 = ['mv outcurve ' outcurve_filename{ii}];
        [status, cmdout] = system(cmd2, '-echo');
    end
end

if 1
    % Check descur
    for ii = 1:num_launch_pts
        check_descur_output(LCFS_filename{ii}, outcurve_filename{ii}, ...
            num_field_periods)
    end
end


