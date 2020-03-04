function [] = VMEC_COORDS_DEMO_2
% the demo
file_vmec_output = 'wout_fnm_edge_6256mm.nc';
TOL_INV = 1e-6;

% close all;

grid_size_1 = 21;
grid_size_2 = 21;
grid_size_3 = 1;

% choose your demo grid. (demo_grid_x, demo_grid_y, demo_grid_z

%demo_grid_x = linspace(6.1, 6.3, grid_size_1);
%demo_grid_y = 0 * demo_grid_x;
%demo_grid_z{1} = 0 * demo_grid_x;

[demo_grid_x, demo_grid_z{1}] = meshgrid(linspace(6.1, 6.3, grid_size_1), ...
    linspace(-.4, .4, grid_size_2));
% [demo_grid_x, demo_grid_z{1}] = meshgrid(linspace(1.3, 1.5, grid_size_1), ...
%     linspace(-.4, .4, grid_size_2));
demo_grid_y = 0*demo_grid_z{1};

% [demo_grid_y, demo_grid_z] = meshgrid(linspace(1.43, 1.47, grid_size_1), ...
%     linspace(-.25, .25, grid_size_2));
% demo_grid_x = 0*demo_grid_z;


% reshape the input grids into an array or vector
x_in = reshape(demo_grid_x,1,grid_size_1 * grid_size_2);
y_in = reshape(demo_grid_y,1,grid_size_1 * grid_size_2);


for jj = 1:grid_size_3
    z_in = reshape(demo_grid_z{jj},1,grid_size_1 * grid_size_2);
    
    % go get the quantities from VMEC output
    [r, z, rho, theta, r_rho, ...
        r_theta, r_phi, z_rho, ...
        z_theta, z_phi, ...
        lambda_phi, lambda_theta, ...
        b_r, b_z, b_phi, ...
        sqrt_g, bsq_out, b1_out, b_x, b_y, ...
        Psi_Phi, Psi_R, Psi_Z, Psi_x, Psi_y] = ...
        get_VMEC_derived_values(file_vmec_output, x_in, y_in, z_in, TOL_INV);
    
    r_bs = sqrt(x_in.^2 + y_in.^2);
    phi_bs = atan2(y_in, x_in);
    z_bs = z_in;
    
    tic
    % check against grid_interp for reconstruction of B
    for ii = 1:length(x_in)
        if isnan(rho(ii)) || rho(ii) > 1
            bx_bs(ii) = NaN;
            by_bs(ii) = NaN;
            bz_bs(ii) = NaN;
        else
            try
                [bx_bs(ii),by_bs(ii),bz_bs(ii)] =...
                    hsxfield_grid_interp(r_bs(ii), phi_bs(ii), z_bs(ii), ...
                    10722, [0 0 0 0 0 0]);
            catch
                bx_bs(ii) = NaN;
                by_bs(ii) = NaN;
                bz_bs(ii) = NaN;
                %             [bx_bs(ii),by_bs(ii),bz_bs(ii)] =...
                %                 bs_derivs_aux_mex(r_bs(ii), phi_bs(ii), z_bs(ii), ...
                %                 10722, [0 0 0 0 0 0]);
            end
        end
    end
    br_bs = bx_bs.*cos(phi_bs) + by_bs.*sin(phi_bs);
    bphi_bs = -bx_bs .* sin(phi_bs) + by_bs.*cos(phi_bs);
    disp(['Time in hsxfield_grid_interp and/or bs_derivs_aux_mex = ' num2str(toc) ' seconds.']);
    
    % reshape the output grids back to match the original grid
    r_plot = reshape(r, grid_size_1, grid_size_2);
    z_plot = reshape(z, grid_size_1, grid_size_2);
    rho_plot = reshape(rho, grid_size_1, grid_size_2);
    theta_plot = reshape(theta, grid_size_1, grid_size_2);
    dr_drho_plot = reshape(r_rho, grid_size_1, grid_size_2);
    dr_dtheta_plot = reshape(r_theta, grid_size_1, grid_size_2);
    dr_dphi_plot = reshape(r_phi, grid_size_1, grid_size_2);
    dz_drho_plot = reshape(z_rho, grid_size_1, grid_size_2);
    dz_dtheta_plot = reshape(z_theta, grid_size_1, grid_size_2);
    dz_dphi_plot = reshape(z_phi, grid_size_1, grid_size_2);
    dlambda_dphi_plot = reshape(lambda_phi, grid_size_1, grid_size_2);
    dlambda_dtheta_plot = reshape(lambda_theta, grid_size_1, grid_size_2);
    b_r_plot = reshape(b_r, grid_size_1, grid_size_2);
    b_z_plot = reshape(b_z, grid_size_1, grid_size_2);
    b_phi_plot = reshape(b_phi, grid_size_1, grid_size_2);
    sqrtg_plot = reshape(sqrt_g, grid_size_1, grid_size_2);
    b_x_plot = reshape(b_x, grid_size_1, grid_size_2);
    b_y_plot = reshape(b_y, grid_size_1, grid_size_2);
    Psi_Phi_plot = reshape(Psi_Phi, grid_size_1, grid_size_2);
    Psi_R_plot = reshape(Psi_R, grid_size_1, grid_size_2);
    Psi_Z_plot = reshape(Psi_Z, grid_size_1, grid_size_2);
    Psi_x_plot = reshape(Psi_x, grid_size_1, grid_size_2);
    Psi_y_plot = reshape(Psi_y, grid_size_1, grid_size_2);
    
    bx_bs_plot = reshape(bx_bs, grid_size_1, grid_size_2);
    by_bs_plot = reshape(by_bs, grid_size_1, grid_size_2);
    bz_bs_plot = reshape(bz_bs, grid_size_1, grid_size_2);
    br_bs_plot = reshape(br_bs, grid_size_1, grid_size_2);
    bphi_bs_plot = reshape(bphi_bs, grid_size_1, grid_size_2);
    
     fig_ind = 1001;figure(fig_ind);surf(demo_grid_x, demo_grid_y, demo_grid_z{jj}, r_plot, 'FaceColor', 'Flat', 'EdgeColor', 'None');title([strrep('r_plot', '_', '\_')]);colorbar;hold on;axis equal;
%     fig_ind = fig_ind + 1;figure(fig_ind);surf(demo_grid_x, demo_grid_y, demo_grid_z{jj}, z_plot, 'FaceColor', 'Flat', 'EdgeColor', 'None');title([strrep('z_plot', '_', '\_')]);colorbar;hold on;axis equal;
    fig_ind = fig_ind + 1;figure(fig_ind);surf(demo_grid_x, demo_grid_y, demo_grid_z{jj}, rho_plot, 'FaceColor', 'Flat', 'EdgeColor', 'None');title([strrep('rho_plot', '_', '\_')]);colorbar;hold on;axis equal;
%     fig_ind = fig_ind + 1;figure(fig_ind);surf(demo_grid_x, demo_grid_y, demo_grid_z{jj}, theta_plot, 'FaceColor', 'Flat', 'EdgeColor', 'None');title([strrep('theta_plot', '_', '\_')]);colorbar;hold on;axis equal;
%     fig_ind = fig_ind + 1;figure(fig_ind);surf(demo_grid_x, demo_grid_y, demo_grid_z{jj}, dr_drho_plot, 'FaceColor', 'Flat', 'EdgeColor', 'None');title([strrep('dr_drho_plot', '_', '\_')]);colorbar;hold on;axis equal;
%     fig_ind = fig_ind + 1;figure(fig_ind);surf(demo_grid_x, demo_grid_y, demo_grid_z{jj}, dr_dtheta_plot, 'FaceColor', 'Flat', 'EdgeColor', 'None');title([strrep('dr_dtheta_plot', '_', '\_')]);colorbar;hold on;axis equal;
%     fig_ind = fig_ind + 1;figure(fig_ind);surf(demo_grid_x, demo_grid_y, demo_grid_z{jj}, dr_dphi_plot, 'FaceColor', 'Flat', 'EdgeColor', 'None');title([strrep('dr_dphi_plot', '_', '\_')]);colorbar;hold on;axis equal;
%     fig_ind = fig_ind + 1;figure(fig_ind);surf(demo_grid_x, demo_grid_y, demo_grid_z{jj}, dz_drho_plot, 'FaceColor', 'Flat', 'EdgeColor', 'None');title([strrep('dz_drho_plot', '_', '\_')]);colorbar;hold on;axis equal;
%     fig_ind = fig_ind + 1;figure(fig_ind);surf(demo_grid_x, demo_grid_y, demo_grid_z{jj}, dz_dtheta_plot, 'FaceColor', 'Flat', 'EdgeColor', 'None');title([strrep('dz_dtheta_plot', '_', '\_')]);colorbar;hold on;axis equal;
%     fig_ind = fig_ind + 1;figure(fig_ind);surf(demo_grid_x, demo_grid_y, demo_grid_z{jj}, dz_dphi_plot, 'FaceColor', 'Flat', 'EdgeColor', 'None');title([strrep('dz_dphi_plot', '_', '\_')]);colorbar;hold on;axis equal;
%     fig_ind = fig_ind + 1;figure(fig_ind);surf(demo_grid_x, demo_grid_y, demo_grid_z{jj}, dlambda_dphi_plot, 'FaceColor', 'Flat', 'EdgeColor', 'None');title([strrep('dlambda_dphi_plot', '_', '\_')]);colorbar;hold on;axis equal;
%     fig_ind = fig_ind + 1;figure(fig_ind);surf(demo_grid_x, demo_grid_y, demo_grid_z{jj}, dlambda_dtheta_plot, 'FaceColor', 'Flat', 'EdgeColor', 'None');title([strrep('dlambda_dtheta_plot', '_', '\_')]);colorbar;hold on;axis equal;
%     fig_ind = fig_ind + 1;figure(fig_ind);subplot(1,2,1);surf(demo_grid_x, demo_grid_y, demo_grid_z{jj}, b_r_plot, 'FaceColor', 'Flat', 'EdgeColor', 'None');title([strrep('b_r_plot', '_', '\_')]);colorbar;hold on;axis equal;
%     subplot(1,2,2);surf(demo_grid_x, demo_grid_y, demo_grid_z{jj}, br_bs_plot, 'FaceColor', 'Flat', 'EdgeColor', 'None');title([strrep('br_bs_or_gi_plot', '_', '\_')]);colorbar;hold on;axis equal;
%     fig_ind = fig_ind + 1;figure(fig_ind);subplot(1,2,1);surf(demo_grid_x, demo_grid_y, demo_grid_z{jj}, b_phi_plot, 'FaceColor', 'Flat', 'EdgeColor', 'None');title([strrep('b_phi_plot', '_', '\_')]);colorbar;hold on;axis equal;
%     subplot(1,2,2);surf(demo_grid_x, demo_grid_y, demo_grid_z{jj}, bphi_bs_plot, 'FaceColor', 'Flat', 'EdgeColor', 'None');title([strrep('bphi_bs_or_gi_plot', '_', '\_')]);colorbar;hold on;axis equal;
%     fig_ind = fig_ind + 1;figure(fig_ind);subplot(1,2,1);surf(demo_grid_x, demo_grid_y, demo_grid_z{jj}, b_z_plot, 'FaceColor', 'Flat', 'EdgeColor', 'None');title([strrep('b_z_plot', '_', '\_')]);colorbar;axis equal;hold on
%     subplot(1,2,2);surf(demo_grid_x, demo_grid_y, demo_grid_z{jj}, bz_bs_plot, 'FaceColor', 'Flat', 'EdgeColor', 'None');title([strrep('bz_bs_or_gi_plot', '_', '\_')]);colorbar;axis equal;hold on
%     fig_ind = fig_ind + 1;figure(fig_ind);surf(demo_grid_x, demo_grid_y, demo_grid_z{jj}, sqrtg_plot, 'FaceColor', 'Flat', 'EdgeColor', 'None');title([strrep('sqrtg_plot', '_', '\_')]);colorbar;axis equal;hold on
%     fig_ind = fig_ind + 1;figure(fig_ind);subplot(1,2,1);surf(demo_grid_x, demo_grid_y, demo_grid_z{jj}, b_x_plot, 'FaceColor', 'Flat', 'EdgeColor', 'None');title([strrep('b_x_plot', '_', '\_')]);colorbar;axis equal;hold on
%     subplot(1,2,2);surf(demo_grid_x, demo_grid_y, demo_grid_z{jj}, bx_bs_plot, 'FaceColor', 'Flat', 'EdgeColor', 'None');title([strrep('bx_bs_or_gi_plot', '_', '\_')]);colorbar;axis equal;hold on
%     fig_ind = fig_ind + 1;figure(fig_ind);subplot(1,2,1);surf(demo_grid_x, demo_grid_y, demo_grid_z{jj}, b_y_plot, 'FaceColor', 'Flat', 'EdgeColor', 'None');title([strrep('b_y_plot', '_', '\_')]);colorbar;axis equal;hold on
%     subplot(1,2,2);surf(demo_grid_x, demo_grid_y, demo_grid_z{jj}, by_bs_plot, 'FaceColor', 'Flat', 'EdgeColor', 'None');title([strrep('by_bs_or_gi_plot', '_', '\_')]);colorbar;axis equal;hold on
%     fig_ind = fig_ind + 1;figure(fig_ind);surf(demo_grid_x, demo_grid_y, demo_grid_z{jj}, Psi_Phi_plot, 'FaceColor', 'Flat', 'EdgeColor', 'None');title([strrep('Psi_Phi_plot', '_', '\_')]);colorbar;axis equal;axis equal;hold on
%     fig_ind = fig_ind + 1;figure(fig_ind);surf(demo_grid_x, demo_grid_y, demo_grid_z{jj}, Psi_R_plot, 'FaceColor', 'Flat', 'EdgeColor', 'None');title([strrep('Psi_R_plot', '_', '\_')]);colorbar;axis equal;hold on
%     fig_ind = fig_ind + 1;figure(fig_ind);surf(demo_grid_x, demo_grid_y, demo_grid_z{jj}, Psi_Z_plot, 'FaceColor', 'Flat', 'EdgeColor', 'None');title([strrep('Psi_Z_plot', '_', '\_')]);colorbar;axis equal;hold on
%     fig_ind = fig_ind + 1;figure(fig_ind);surf(demo_grid_x, demo_grid_y, demo_grid_z{jj}, Psi_x_plot, 'FaceColor', 'Flat', 'EdgeColor', 'None');title([strrep('Psi_x_plot', '_', '\_')]);colorbar;axis equal;hold on
%     fig_ind = fig_ind + 1;figure(fig_ind);surf(demo_grid_x, demo_grid_y, demo_grid_z{jj}, Psi_y_plot, 'FaceColor', 'Flat', 'EdgeColor', 'None');title([strrep('Psi_y_plot', '_', '\_')]);colorbar;axis equal;hold on
%    figure;
%    plot(demo_grid_x, Psi_R_plot, 'o:');
    
end
    
    
