function make_lcfs_eim252_edge

% definitions
coilsetID = 'coilset_w7x_v3';
configuration_name = 'eim252_edge';
coilCurrents  = [12989 12989 12989 12989 12989 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
%coilCurrents = 1.;
title_str='W7X EDGE';

Rlaunch = [6.238:.001:6.27];
num_field_periods = 1;
nphi_step1 = 360;
ntheta_step1 = 250;

num_launch_pts = length(Rlaunch);
nphi_step2 = nphi_step1 / num_field_periods;
ntheta_step2 = ntheta_step1 * num_field_periods;

do_the_descur_stuff(coilsetID, coilCurrents, configuration_name, title_str, ...
    Rlaunch, num_field_periods, nphi_step1, ntheta_step1, num_launch_pts, ...
    nphi_step2, ntheta_step2, 13);

%=====

