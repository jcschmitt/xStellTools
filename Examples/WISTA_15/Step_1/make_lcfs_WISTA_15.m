function make_lcfs_WISTA_15

% definitions
coilsetID = 'coilset_wista_15';
configuration_name = 'wista_15';
coilCurrents = 632955.5;
%coilCurrents = 1.;
title_str='WISTA 15';

Rlaunch = 2.58 + 0.002 * (0:13);

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

