function [zetamat, thetamat, theta_star_vmec, zeta_boozer, theta_boozer] = ...
    generate_vmec_and_boozer_coords(file_vmec_outputa, surfs_to_make, numtheta, numzeta)
% function [zetamat, thetamat, theta_star_vmec, zeta_boozer, theta_boozer] = ...
%     generate_vmec_and_boozer_coords(file_vmec_outputa, surfs_to_make, numtheta, numzeta)
%
% INPUTS
% file_vmec_outputa: the vmec output extension, i.e. 'QHS' for wout_QHS.nc
% surfs_to_make:  the surface #s you wish to calculate, i.e., 1:201
% numtheta: number of poloidal angle divisions (64 works well)
% numzeta:  number of toroidal angle divisions (400 works well)
%
% OUTPUTS:
% zetamat: a matrix of values of the laboratory poloidal angle.
% thetamat: a matrix of values of the laboratory toroidal angle.
% theta_star_vmec: cells of matrics.  one cell/surface. each matrix
% contains the vmec poloidal angle corresponding to the laboratory
% poloidal, toroidal angle matrics
% zeta_boozer: cells of matrices.  each matrix contains the boozer toroidal
% angle corresponding to the laboratory poloidal&toroidal angle matrices.
% theta_boozer: cells of matices.  contains the boozer poloidal angle
%
%  Adapted from BOOZ_XFORM
%
%=============================================
% Updates
%=============================================
%
% 2009/12/02 - Added native netcdf support for Matlab 7.6+ - JW Radder
% 2009/11/18 - version 1.0a - JC Schmitt
%=============================================

% get matlab version, verison >= 7.6 has native netcdf support
matlabVersion = version;

theta=linspace(0, 2*pi, numtheta);
zeta=linspace(0, pi/2, numzeta);
% [thetamat, zetamat] = meshgrid(theta, zeta);
[zetamat, thetamat] = meshgrid(zeta, theta);

% load the variables
file_vmec_output = ['wout_' file_vmec_outputa '.nc'];
% input_extension_vmec = getnc(file_vmec_output,'input_extension');
save_data_filename = ['vmec_and_boozer_coords_' file_vmec_outputa];

if (7.2 <= str2num(matlabVersion(1:3))) & (str2num(matlabVersion(1:3)) < 7.61)
    %USED
    disp('Using 3rd party netcdf support');
    ns_vmec = getnc(file_vmec_output,'ns');
    % mpol_vmec = getnc(file_vmec_output,'mpol');
    % ntor_vmec = getnc(file_vmec_output,'ntor');
    mnmax_vmec = getnc(file_vmec_output,'mnmax');
    mnmax_nyq_vmec = getnc(file_vmec_output,'mnmax_nyq');
    
    %(mn_mode)  poloidal mode numbers
    xm_vmec = getnc(file_vmec_output,'xm');
    %(mn_mode)  toroidal mode numbers
    xn_vmec = getnc(file_vmec_output,'xn');
    
    xm_nyq_vmec = getnc(file_vmec_output,'xm_nyq');      %(mn_mode)
    xn_nyq_vmec = getnc(file_vmec_output,'xn_nyq');      %(mn_mode)
    
    %(radius mn_mode ) cosmn component of cylindrical R, full mesh in [m]
    rmnc_vmec = getnc(file_vmec_output,'rmnc');
    %(radius mn_mode ) sinmn component of cylindrical Z, full mesh in [m]
    zmns_vmec = getnc(file_vmec_output,'zmns');
    %(radius mn_mode ) sinmn component of lambda, half mesh
    lmns_vmec = getnc(file_vmec_output,'lmns');
    %(radius mn_mode_nyq ) cosmn component of jacobian, half mesh
    gmnc_vmec = getnc(file_vmec_output,'gmnc');
    %(radius mn_mode_nyq ) cosmn component of mod-B, half mesh
    bmnc_vmec = getnc(file_vmec_output,'bmnc');
    %(radius mn_mode_nyq ) cosmn component of bsubu (poloidal covariant),
    % half mesh
    bsubumnc_vmec = getnc(file_vmec_output, 'bsubumnc');
    %(radius mn_mode_nyq ) cosmn component of bsubv (toroidal covariant),
    % half mesh
    bsubvmnc_vmec = getnc(file_vmec_output, 'bsubvmnc');
    iota_vmec = getnc(file_vmec_output,'iotaf');  % full mesh iota
    % extcur = getnc(file_vmec_output,'extcur')
    
elseif (str2num(matlabVersion(1:3)) >= 7.6 )
    disp('Using matlab built-in netcdf support');
    ncID = netcdf.open(file_vmec_output,'nc_nowrite');
    
    ns_vmecID = netcdf.inqVarID(ncID, 'ns');
    ns_vmec = netcdf.getVar(ncID, ns_vmecID);
    clear ns_vmecID;
    
    mnmax_vmecID = netcdf.inqVarID(ncID, 'mnmax');
    mnmax_vmec = netcdf.getVar(ncID, mnmax_vmecID);
    clear mnmax_vmecID;
    
    mnmax_nyq_vmecID = netcdf.inqVarID(ncID, 'mnmax_nyq');
    mnmax_nyq_vmec = netcdf.getVar(ncID, mnmax_nyq_vmecID);
    clear mnmax_nyq_vmecID;
    
    %(mn_mode)  poloidal mode numbers
    xm_vmecID = netcdf.inqVarID(ncID, 'xm');
    xm_vmec = netcdf.getVar(ncID, xm_vmecID);
    clear xm_vmecID;
    
    %(mn_mode)  toroidal mode numbers
    xn_vmecID = netcdf.inqVarID(ncID, 'xn');
    xn_vmec = netcdf.getVar(ncID, xn_vmecID);
    clear xn_vmecID;
    
    %(mn_mode)
    xm_nyq_vmecID = netcdf.inqVarID(ncID, 'xm_nyq');
    xm_nyq_vmec = netcdf.getVar(ncID, xm_nyq_vmecID);
    clear xm_nyq_vmecID;
    
    %(mn_mode)
    xn_nyq_vmecID = netcdf.inqVarID(ncID, 'xn_nyq');
    xn_nyq_vmec = netcdf.getVar(ncID, xn_nyq_vmecID);
    clear xn_nyq_vmecID;
    
    %(radius mn_mode ) cosmn component of cylindrical R, full mesh in [m]
    rmnc_vmecID = netcdf.inqVarID(ncID, 'rmnc');
    rmnc_vmec = netcdf.getVar(ncID, rmnc_vmecID)';
    clear rmnc_vmecID;
    
    %(radius mn_mode ) sinmn component of cylindrical Z, full mesh in [m]
    zmns_vmecID = netcdf.inqVarID(ncID, 'zmns');
    zmns_vmec = netcdf.getVar(ncID, zmns_vmecID)';
    clear zmns_vmecID;
    
    %(radius mn_mode ) sinmn component of lambda, half mesh
    lmns_vmecID = netcdf.inqVarID(ncID, 'lmns');
    lmns_vmec = netcdf.getVar(ncID, lmns_vmecID)';
    clear lmns_vmecID;
    
    %(radius mn_mode_nyq ) cosmn component of jacobian, half mesh
    gmnc_vmecID = netcdf.inqVarID(ncID, 'gmnc');
    gmnc_vmec = netcdf.getVar(ncID, gmnc_vmecID)';
    clear gmnc_vmecID;
    
    %(radius mn_mode_nyq ) cosmn component of mod-B, half mesh
    bmnc_vmecID = netcdf.inqVarID(ncID, 'bmnc');
    bmnc_vmec = netcdf.getVar(ncID, bmnc_vmecID)';
    clear bmnc_vmecID;
    
    %(radius mn_mode_nyq ) cosmn component of bsubu (poloidal covariant),
    % half mesh
    bsubumnc_vmecID = netcdf.inqVarID(ncID, 'bsubumnc');
    bsubumnc_vmec = netcdf.getVar(ncID, bsubumnc_vmecID)';
    clear bsubumnc_vmecID;
    
    %(radius mn_mode_nyq ) cosmn component of bsubv (toroidal covariant),
    % half mesh
    bsubvmnc_vmecID = netcdf.inqVarID(ncID, 'bsubvmnc');
    bsubvmnc_vmec = netcdf.getVar(ncID, bsubvmnc_vmecID)';
    clear bsubvmnc_vmecID;
    
    % full mesh iota
    iota_vmecID = netcdf.inqVarID(ncID, 'iotaf');
    iota_vmec = netcdf.getVar(ncID, iota_vmecID);
    clear iota_vmecID;
    
    % extcur = getnc(file_vmec_output,'extcur')
    extcurID = netcdf.inqVarID(ncID, 'extcur');
    extcur = netcdf.getVar(ncID, extcurID);
    clear extcurID;
    
    netcdf.close(ncID);
    
end;

% initialize variable arrays (cells)
R_vmec = cell(1,ns_vmec);
Z_vmec = cell(1,ns_vmec);
B_vmec = cell(1,ns_vmec);
g_vmec = cell(1,ns_vmec);
lambda_vmec = cell(1,ns_vmec);
theta_star_vmec = cell(1,ns_vmec);
gpsi = cell(1,ns_vmec);
Ipsi = cell(1,ns_vmec);
p_boozer = cell(1,ns_vmec);
u_boozer = cell(1,ns_vmec);
v_boozer = cell(1,ns_vmec);
zeta_boozer = cell(1,ns_vmec);
theta_boozer = cell(1,ns_vmec);
w = cell(1,ns_vmec);
Pmn = cell(1,ns_vmec);

bsubtmn = bsubumnc_vmec; % boozer labels: u (poloidal) -> t (theta)
bsubzmn = bsubvmnc_vmec; %                v (toroidal) -> z (zeta)

for surf_ind=surfs_to_make
    disp(['Working on surface# ' num2str(surf_ind)]);
    pause(0.01)
    % initializing matrices (pre-allocate memory)
    B_vmec{surf_ind} = 0*thetamat;
    g_vmec{surf_ind} = 0*thetamat;
    R_vmec{surf_ind} = 0*thetamat;
    Z_vmec{surf_ind} = 0*thetamat;
    lambda_vmec{surf_ind} = 0*thetamat;
    theta_star_vmec{surf_ind} = 0*thetamat;
    p_boozer{surf_ind} = 0*thetamat;
    w{surf_ind} = 0*thetamat;
    Pmn{surf_ind} = zeros(1,mnmax_nyq_vmec);
    u_boozer{surf_ind} = 0*thetamat;
    v_boozer{surf_ind} = 0*thetamat;
    zeta_boozer{surf_ind} = 0*thetamat;
    theta_boozer{surf_ind} = 0*thetamat;
    
    for mn_ind=1:mnmax_nyq_vmec
        cosmn=cos(xm_nyq_vmec(mn_ind)*thetamat - xn_nyq_vmec(mn_ind)*zetamat);
        sinmn=sin(xm_nyq_vmec(mn_ind)*thetamat - xn_nyq_vmec(mn_ind)*zetamat);
        
        B_vmec{surf_ind} = B_vmec{surf_ind} + bmnc_vmec(surf_ind,mn_ind) * cosmn;
        g_vmec{surf_ind} = g_vmec{surf_ind} + gmnc_vmec(surf_ind,mn_ind) * cosmn;
        
        if xn_nyq_vmec(mn_ind) ~= 0
            Pmn{surf_ind}(mn_ind) = - bsubzmn(surf_ind, mn_ind) / xn_nyq_vmec(mn_ind);
        elseif xm_nyq_vmec(mn_ind) ~= 0
            Pmn{surf_ind}(mn_ind) = bsubtmn(surf_ind, mn_ind) / xm_nyq_vmec(mn_ind);
            %             Pmn{surf_ind}(mn_ind) = -bsubtmn(surf_ind, mn_ind) / xm_nyq_vmec(mn_ind);
        else
            Pmn{surf_ind}(mn_ind) = 0;
            gpsi{surf_ind} = bsubzmn(surf_ind, mn_ind);
            Ipsi{surf_ind} = bsubtmn(surf_ind, mn_ind);
        end
        
        % jc schmitt
%         w{surf_ind} = w{surf_ind} + Pmn{surf_ind}(mn_ind) * sinmn;
        w{surf_ind} = w{surf_ind} + Pmn{surf_ind}(mn_ind) * sinmn;
        
    end
    
    for mn_ind=1:mnmax_vmec
        cosmn=cos(xm_vmec(mn_ind)*thetamat - xn_vmec(mn_ind)*zetamat);
        sinmn=sin(xm_vmec(mn_ind)*thetamat - xn_vmec(mn_ind)*zetamat);
        
        R_vmec{surf_ind} = R_vmec{surf_ind} + rmnc_vmec(surf_ind,mn_ind) * cosmn;
        Z_vmec{surf_ind} = Z_vmec{surf_ind} + zmns_vmec(surf_ind,mn_ind) * sinmn;
        lambda_vmec{surf_ind} = lambda_vmec{surf_ind} + lmns_vmec(surf_ind,mn_ind) * sinmn;
    end
    
    x_vmec{surf_ind} = R_vmec{surf_ind} .* cos(zetamat);
    y_vmec{surf_ind} = R_vmec{surf_ind} .* sin(zetamat);
    
    jacfac = gpsi{surf_ind} + iota_vmec(surf_ind) * Ipsi{surf_ind};
    dem = 1/jacfac;
        gpsi_1 = gpsi{surf_ind} * dem;
    hiota_1 = iota_vmec(surf_ind) * dem;
        u_boozer{surf_ind} = gpsi_1 * lambda_vmec{surf_ind} + ...
            hiota_1 * w{surf_ind};
%     u_boozer{surf_ind} = dem * (gpsi{surf_ind} * lambda_vmec{surf_ind} + ...
%         hiota_1 * w{surf_ind});
    v_boozer{surf_ind} = dem * ( w{surf_ind} - ...
        Ipsi{surf_ind} * lambda_vmec{surf_ind} );
    
    theta_star_vmec{surf_ind} = thetamat + lambda_vmec{surf_ind};
    zeta_boozer{surf_ind} = zetamat + v_boozer{surf_ind};
    theta_boozer{surf_ind} = thetamat + lambda_vmec{surf_ind} + ...
        u_boozer{surf_ind};
end


% save(save_data_filename)
% disp(['Finished.  Saved files to ' save_data_filename '.mat'])
disp(['Finished.  Did not save files'])

% make some plots
for surf_ind=surfs_to_make(end)
    
    figure
    surf(x_vmec{surf_ind}, y_vmec{surf_ind}, Z_vmec{surf_ind}, B_vmec{surf_ind}, 'FaceColor', 'flat','edgecolor', 'none');
    axis equal
    xlabel('x'); ylabel('y'); zlabel('z');
    view(3);colorbar;view([125 25])

        figure
    surf(x_vmec{surf_ind}, y_vmec{surf_ind}, Z_vmec{surf_ind}, lambda_vmec{surf_ind}, 'FaceColor', 'flat','edgecolor', 'none');
    axis equal
    xlabel('x'); ylabel('y'); zlabel('z');
    view(3);colorbar;view([125 25])

            figure
    surf(x_vmec{surf_ind}, y_vmec{surf_ind}, Z_vmec{surf_ind}, theta_star_vmec{surf_ind}, 'FaceColor', 'flat','edgecolor', 'none');
    axis equal
    xlabel('x'); ylabel('y'); zlabel('z');
    view(3);colorbar;view([125 25])

    
    B_line_zeta_0 = 0;  B_line_zeta_end = .5*pi;
    B_line_theta_0 = 0; B_line_theta_end = iota_vmec(surf_ind) * ...
        B_line_zeta_end + B_line_theta_0;
    
    
    figure
%     mesh(zetamat, thetamat, B_vmec{surf_ind}, 'FaceColor', 'flat')
    meshc(zetamat, thetamat, B_vmec{end})
    hold on;
    plot3([B_line_zeta_0 B_line_zeta_end], [B_line_theta_0 B_line_theta_end], [0 0], 'k--', 'Linewidth', 2);
    view(2)
    axis([0 pi/2 0 2*pi])
    xlabel('Toroidal VMEC Angle  (rad)')
    ylabel('Poloidal VMEC Angle, \theta  (rad)')
    title('|B| on LCFS from VMEC')
    colorbar
    
    figure
%     mesh(zetamat, theta_star_vmec{surf_ind}, B_vmec{surf_ind}, 'FaceColor', 'flat')
    meshc(zetamat, theta_star_vmec{end}, B_vmec{end})
    hold on;
    plot3([B_line_zeta_0 B_line_zeta_end], [B_line_theta_0 B_line_theta_end], [0 0], 'k--', 'Linewidth', 2);
    view(2)
    axis([0 pi/2 0 2*pi])
    xlabel('Toroidal VMEC Angle  (rad)')
    ylabel('Poloidal VMEC Angle, \theta^*  (rad)')
    title('|B| on LCFS from VMEC')
    colorbar
    
    figure
    mesh(zeta_boozer{surf_ind}, theta_boozer{surf_ind}, B_vmec{surf_ind}, 'FaceColor', 'flat')
    hold on;
    plot3([B_line_zeta_0 B_line_zeta_end], [B_line_theta_0 B_line_theta_end], [0 0], 'k--', 'Linewidth', 2);
    view(2)
    axis([0 pi/2 0 2*pi])
    xlabel('Toroidal Boozer Angle  (rad)')
    ylabel('Poloidal Boozer Angle  (rad)')
    title('|B| on LCFS from VMEC')
    colorbar
    
    figure
    % contour(zeta_boozer{end}, theta_boozer{end}, B_vmec{end})
    contourf(zetamat, thetamat, B_vmec{surf_ind}, 10)
    hold on;
    plot3([B_line_zeta_0 B_line_zeta_end], [B_line_theta_0 B_line_theta_end], [0 0], 'k--', 'Linewidth', 2);
    view(2)
    axis([0 pi/2 0 2*pi])
    xlabel('Toroidal VMEC Angle  (rad)')
    ylabel('Poloidal VMEC Angle, \theta  (rad)')
    title('|B| on LCFS from VMEC')
    colorbar
    
    figure
    % mesh(zetamat, theta_star_vmec{end}, B_vmec{end}, 'FaceColor', 'flat')
    contourf(zetamat, theta_star_vmec{surf_ind}, B_vmec{surf_ind}, 10)
    hold on;
    plot3([B_line_zeta_0 B_line_zeta_end], [B_line_theta_0 B_line_theta_end], [0 0], 'k--', 'Linewidth', 2);
    view(2)
    axis([0 pi/2 0 2*pi])
    xlabel('Toroidal VMEC Angle  (rad)')
    ylabel('Poloidal VMEC Angle, \theta^*  (rad)')
    title('|B| on LCFS from VMEC')
    colorbar
    
    figure
    % mesh(zeta_boozer{end}, theta_boozer{end}, B_vmec{end}, 'FaceColor', 'flat')
    contourf(zeta_boozer{surf_ind}, theta_boozer{surf_ind}, B_vmec{surf_ind}, 10)
    %     contour(zeta_boozer{surf_ind}, theta_boozer{surf_ind}, B_vmec{surf_ind}, 10)
    hold on;
    plot3([B_line_zeta_0 B_line_zeta_end], [B_line_theta_0 B_line_theta_end], [0 0], 'k--', 'Linewidth', 2);
    view(2)
    axis([0 pi/2 0 2*pi])
    xlabel('Toroidal Boozer Angle  (rad)')
    ylabel('Poloidal Boozer Angle  (rad)')
    title('|B| on LCFS from VMEC')
    colorbar
end
