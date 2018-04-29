function [modB_recon_2D, modB_recon_3D, zeta, alpha] = ...
    reconstructModB(chi, modB, spectrumType, dV_dPsi, g_Boozer, n_values, m_values, ...
    pk_n_sorted, pk_m_sorted, num_components_reconstruct, nm_amp, iota, PLOT_2D, PLOT_3D)
% reconstitute the B-curve from its spactral components

if (nargin < 13)
    PLOT_2D = 1;
end

if (nargin < 14)
    PLOT_3D = 1;
end

if strcmp(lower(spectrumType), 'hamada')
    %         n_minus_miota = omega * dV_dPsi / (2 * pi);
    omega_factor = (2*pi)/dV_dPsi;
elseif strcmp(lower(spectrumType), 'boozer')
    %         n_minus_miota = omega * g_Boozer;
    omega_factor = 1/g_Boozer;
end

modB_recon_2D = zeros(size(chi));

for ii = 1:length(num_components_reconstruct)
    omega_value = (n_values(pk_n_sorted(ii)) - iota * m_values(pk_m_sorted(ii))) * omega_factor;
    modB_recon_2D = modB_recon_2D + nm_amp( pk_n_sorted(ii), pk_m_sorted(ii) ) * cos(chi*omega_value);
end

zeta_0 = 0; alpha_0=0;
zeta_end = 2*pi; alpha_end = 2*pi;
zeta_div = 5e-3;
alpha_div = 5e-3;
zeta = [zeta_0:zeta_div:zeta_end];
alpha = alpha_0:alpha_div:alpha_end;
zeta_matrix = zeta' * ones(1,length(alpha));
alpha_matrix = ones(length(zeta),1) * alpha;
modB_recon_3D = zeros(length(zeta), length(alpha));
modB_recon_3D_Tok = modB_recon_3D;

for ii = 1:length(num_components_reconstruct)
    angle_value = (m_values(pk_m_sorted(ii)) * zeta_matrix - n_values(pk_n_sorted(ii)) * alpha_matrix);
    modB_recon_3D = modB_recon_3D + nm_amp( pk_n_sorted(ii), pk_m_sorted(ii) ) * cos(angle_value);
end


angle_value_Tok = zeta_matrix ;
modB_recon_3D_Tok = nm_amp( pk_n_sorted(1), pk_m_sorted(1) ) + ...
    nm_amp( pk_n_sorted(2), pk_m_sorted(2) )* cos(angle_value_Tok);

% a bunch
% B_0 = -2*pi; B_end = 2*pi;B_div = pi/4;
% BLine_0 = B_0:B_div:B_end;
% just one
BLine_0 = 0;
BLines_alpha = alpha;

for ii = 1:length(BLine_0)
    BLines_zeta(ii,:) = iota * (BLines_alpha - BLine_0(ii));
end

cyl_x = .1*cos(zeta_matrix); cyl_y = .1*sin(zeta_matrix); cyl_z =4*alpha_matrix/(2*pi);
cyl_surf_color = modB_recon_3D;
BLines_cyl_x = .1*cos(BLines_zeta);
BLines_cyl_y = .1*sin(BLines_zeta);

% figure
% surf(alpha, zeta, alpha_matrix);
%     xlabel('\alpha');
%     ylabel('\zeta');

if PLOT_2D
    figure();
    plot(chi, modB_recon_2D, 'b');
    hold on
    len_chi = length(chi);
    plot(chi', modB(len_chi-1:end), 'k:');

end

if PLOT_3D
    figure
    subplot(2,1,1);
    surf(alpha, zeta, modB_recon_3D, 'linestyle', 'none');
    xlabel('\alpha');
    ylabel('\zeta');
    colormap('winter');
    subplot(2,1,2);
    contour(alpha, zeta, modB_recon_3D);
    xlabel('\alpha');
    ylabel('\zeta');
    view(2)
    hold on;
    
    for ii = 1:length(BLine_0)
        plot(BLines_alpha, BLines_zeta(ii, :), 'k:');
    end
    % % %     figure
    % % %     surface(cyl_z, cyl_x, cyl_y, cyl_surf_color, 'linestyle', 'none');
    % % %     hold on;
    % % %     for ii = 1:length(BLine_0)
    % % %         plot3(cyl_z, BLines_cyl_x(ii, :), BLines_cyl_y(ii, :), 'k--');
    % % %     end
    % % %     view(3)
    % % %     axis equal
    % % %     colormap('copper');
    % % %     xlabel('Toroidal direction');
    % % %     ylabel('x-dir');
    % % %     zlabel('y-dir');
    % % %     figure;
    % % %     surface(cyl_z, cyl_x, cyl_y, cyl_surf_color, 'linestyle', 'none');
    % % %     hold on;
    % % %     view(3)
    % % %     axis equal
    % % %     colormap('hsv');
    % % %     xlabel('Toroidal direction');
    % % %     for ii = 1:length(BLine_0)
    % % %         plot3(cyl_z, BLines_cyl_x(ii, :), BLines_cyl_y(ii, :), 'k--');
    % % %     end

    
    figure
    contour(alpha, zeta, modB_recon_3D, 6);
    xlabel('\zeta');ylabel('\theta')
    hold on;
    for ii = 1:length(BLine_0)
        plot3(BLines_alpha, BLines_zeta(ii, :),2*ones(size(BLines_alpha)), 'k:', 'LineWidth', 2);
    end
    colormap('Winter')
   
    figure
 
    surf(alpha, zeta, modB_recon_3D, 'linestyle', 'none');
%     contour(alpha, zeta, modB_recon_3D);
    xlabel('\alpha');
    ylabel('\zeta');
    view(2)
    hold on;
    for ii = 1:length(BLine_0)
        plot3(BLines_alpha, BLines_zeta(ii, :),2*ones(size(BLines_alpha)), 'k:', 'LineWidth', 2);
    end

    axis([0 2*pi 0 2*pi]);
    figure
 
    surf(alpha, zeta, modB_recon_3D_Tok, 'linestyle', 'none');
    xlabel('\alpha');
    ylabel('\zeta');
    view(2)
    hold on;
    for ii = 1:length(BLine_0)
        plot3(BLines_alpha, BLines_zeta(ii, :),2*ones(size(BLines_alpha)), 'k:', 'LineWidth', 2);
    end

    axis([0 2*pi 0 2*pi]);
end

