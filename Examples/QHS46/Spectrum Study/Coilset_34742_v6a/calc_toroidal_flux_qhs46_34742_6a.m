function [] = calc_toroidal_flux_mirfl14_4p()

file_list{1} = 'LCFS_MirFL14_4p_1_503.mat';
file_list{2} = 'LCFS_MirFL14_4p_1_504.mat';
file_list{3} = 'LCFS_MirFL14_4p_1_505.mat';
file_list{4} = 'LCFS_MirFL14_4p_1_506.mat';
file_list{5} = 'LCFS_MirFL14_4p_1_507.mat';
file_list{6} = 'LCFS_MirFL14_4p_1_508.mat';
file_list{7} = 'LCFS_MirFL14_4p_1_509.mat';
file_list{8} = 'LCFS_MirFL14_4p_1_510.mat';
file_list{9} = 'LCFS_MirFL14_4p_1_511.mat';
file_list{10} = 'LCFS_MirFL14_4p_1_512.mat';
file_list{11} = 'LCFS_MirFL14_4p_1_513.mat';
file_list{12} = 'LCFS_MirFL14_4p_1_514.mat';
file_list{13} = 'LCFS_MirFL14_4p_1_515.mat';
file_list{14} = 'LCFS_MirFL14_4p_1_516.mat';
file_list{15} = 'LCFS_MirFL14_4p_1_517.mat';
file_list{16} = 'LCFS_MirFL14_4p_1_518.mat';

R_LauchPoint = 1.503:.001:1.518;
save_name = 'tor_flux_mirFL14_4p';

flux = zeros(size(file_list));
for ii = 1:length(file_list)
    file_data = load(file_list{ii});
    RR = sqrt(file_data.surf_x.^2 + file_data.surf_y.^2);
    ZZ = file_data.surf_z;
    flux(ii) = calc_flux_HSX(RR(1:50:end), ...
        ZZ(1:50:end), 10000, .01*[-1 1 1 1 -1 -1]);

    save(save_name, 'flux', 'R_LauchPoint');
end

