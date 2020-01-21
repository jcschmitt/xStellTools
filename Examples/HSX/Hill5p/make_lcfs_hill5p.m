function make_lcfs_hill5p

% Hill 5%

% definitions
coilsetID = 'coilset_hsx_complete';
configuration_name = 'hill_5p';
coilCurrents = [10116. 708. 708. 708. 708. 708. 708.];
%coilCurrents = 1.;
title_str='HILL 5%';

Rlaunch = 2.58 + 0.001 * (0:23);

num_field_periods = 4;
nphi_step1 = 240;
ntheta_step1 = 250;

num_launch_pts = length(Rlaunch);
nphi_step2 = nphi_step1 / num_field_periods;
ntheta_step2 = ntheta_step1 * num_field_periods;

do_the_descur_stuff(coilsetID, coilCurrents, configuration_name, title_str, ...
    Rlaunch, num_field_periods, nphi_step1, ntheta_step1, num_launch_pts, ...
    nphi_step2, ntheta_step2);

%=====

