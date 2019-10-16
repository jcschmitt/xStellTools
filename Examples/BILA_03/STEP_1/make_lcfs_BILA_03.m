% definitions
coilsetID = 'coilset_bila_03';
configuration_name = 'bila_03';
coilCurrents = 497292.;

Rlaunch = 2.529 + 0.002*[0:13];
num_launch_pts = length(Rlaunch);

num_field_periods = 5;
nphi_step1 = 300;
ntheta_step1 = 200;
nphi_step2 = 60;
ntheta_step2 = 1000;


% make filenames for LCFS and more
fluxout_filename = ['Flux_' configuration_name];

for ii = 1:num_launch_pts
    LCFS_filename{ii} = ['LCFS_' configuration_name '_Rlaunch_' num2str(round(1000*Rlaunch(ii))) 'mm.mat'];
    LCFS_pathname{ii} = './';
    outcurve_filename{ii} = ['outcurve_Rlaunch_' num2str(round(1000*Rlaunch(ii))) 'mm'];
end

% Make poincare plots near edge.
parfor ii = 1:num_launch_pts
    this_R = Rlaunch(ii);
    [surf_x, surf_y, surf_z, surf_phi] = ...
        calc_descur_RZ_points(coilsetID, this_R, nphi_step1, ntheta_step1, coilCurrents);
    save(LCFS_filename{ii}, 'Rlaunch', 'coilCurrents', ...
        'coilsetID', 'configuration_name', 'nphi_step1', 'ntheta_step1', ...
        'surf_x', 'surf_y', 'surf_z', 'surf_phi', 'this_R');
end

% Calculate toroidal flux
tor_flux = 0 * Rlaunch;
parfor ii = 1:num_launch_pts
    load(LCFS_filename{ii})
    R = sqrt(surf_x(1:nphi_step2:end).^2 + surf_y(1:nphi_step2:end).^2);
    Z = surf_z(1:nphi_step2:end)
    tor_flux(ii) = calc_toroidal_flux(coilsetID, coilCurrents, R, Z, MAKE_PLOT)
end

save(fluxout_filename, 'Rlaunch', 'tor_flux');

if 0
    % Convert poincare plots to rzdata input for descur
    for ii = 1:num_launch_pts
        % returns the rz and descur command names
        [outfilename{ii}, descurcmdfilename{ii}] = make_descur_file(LCFS_filename, LCFS_pathname, num_field_periods, nphi_step2, ntheta_step2);
    end
    
    % Run descur and copy files
    for ii = 1:num_launch_pts
        cmd1 = ['~/src/stellinstall/trunk/build/bin/xdescur < ' ...
            descurcmdfilename{ii}];
        [status, cmdout] = system(cmd1, '-echo');
        cmd2 = ['mv outcurve ' outcurve_filename{ii}];
        [status, cmdout] = system(cmd2, '-echo');        
    end
end

% Check descur
for ii = 1:num_launch_pts
    check_descur_output(LCFS_pathname{ii}, outcurve_file, num_field_periods)
end



