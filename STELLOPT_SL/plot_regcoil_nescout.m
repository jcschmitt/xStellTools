function  plot_regcoil_nescout( data, m_max, n_max )
%PLOT_REGCOIL_NESCOUT
% For plotting the output from stellopt/regcoil/nescout runs

if nargin < 3
    n_max = inf;
end
if nargin < 2
    m_max = inf;
end

% For the 3-D surface
numzeta = 64;
numtheta = 64;
numzetal = numzeta * data.nfp;
zeta = linspace(0,2*pi/data.nfp, numzeta);
zetal = linspace(0,2*pi, numzetal);
theta = linspace(0,2*pi,numtheta);

%[theta_map1fp, zeta_map1fp] = meshgrid(theta, zeta);
[theta_mapallfp, zeta_mapallfp] = meshgrid(theta, zetal);

%surf_x_1fp = 0 * theta_map1fp;
%surf_y_1fp = 0 * theta_map1fp;
%surf_z_1fp = 0 * theta_map1fp;
surf_x_allfp = 0 * theta_mapallfp;
surf_y_allfp = 0 * theta_mapallfp;
surf_z_allfp = 0 * theta_mapallfp;

% quantities used later
m_all = [];
n_all = [];
rmnc_all = [];
zmns_all = [];

for ii = 1:numzetal
    for jj = 1:numtheta
        for kk = 1:data.mn_surf
            m = data.xm_surf(kk);
            n = data.xn_surf(kk);
            % stellarator symmetry only for now...
            if ( (abs(m) < m_max) && (abs(n) < n_max) )
                rmnc = data.rmnc_surf(kk);
                zmns = data.zmns_surf(kk);
            else
                rmnc = 0; zmns = 0;
            end
            angle = m*theta(jj) + n*zetal(ii);
            % angle = m*theta(jj) + data.nfp*zetal(ii);
            cosangle = cos(angle);
            sinangle = sin(angle);
            coszeta = cos(zetal(ii));
            sinzeta = sin(zetal(ii));
            surf_x_allfp(ii, jj) = surf_x_allfp(ii,jj) + ...
                rmnc * cosangle * coszeta;
            surf_y_allfp(ii, jj) = surf_y_allfp(ii,jj) + ...
                rmnc * cosangle * sinzeta;
            surf_z_allfp(ii, jj) = surf_z_allfp(ii,jj) + ...
                zmns * sinangle;
            if ( (ii == 1) && (jj == 1))
                m_all(kk) = m;
                n_all(kk) = n;
                rmnc_all(kk) = rmnc;
                zmns_all(kk) = zmns;
            end
        end
    end
end

figure
%subplot(2,1,1);
%plot3(surf_x_allfp, surf_y_allfp, surf_z_allfp, '.');
%axis equal;axis tight
%subplot(2,1,2);
surf(surf_x_allfp, surf_y_allfp, surf_z_allfp, zeros(size(surf_z_allfp)));
axis equal;axis tight
%surf(surf_x_allfp, surf_y_allfp, surf_z_allfp);

% For the 2-d spectral plot
grid_factor = 15;
exp_factor = 0.2;
contour_m_array_limit = max(abs(m));
contour_n_array_limit = max(abs(n));
m_array = linspace(-contour_m_array_limit, contour_m_array_limit, ...
    contour_m_array_limit*grid_factor);
n_array = linspace(-contour_n_array_limit, contour_n_array_limit, ...
    contour_n_array_limit*grid_factor);
[m_plot, n_plot] = meshgrid(m_array, n_array);

magnitude_plot_rmnc = 0* m_plot;
magnitude_plot_zmns = 0* m_plot;

for ii = 1:length(m_all)
    if 1 %~(( m_all(ii) == 0) && (n_all(ii) == 0))
        mask = 0*m_plot;
        mask(find( ((m_plot - m_all(ii)).^2 + ...
                    (n_plot - n_all(ii)).^2 )<0.45)) = 1; 
        exp_add_rmnc = rmnc_all(ii) * mask .* exp(...
            -( ( (m_plot - m_all(ii)).^2 + (n_plot - n_all(ii)).^2 ) / ...
            exp_factor^2 ) );
        exp_add_zmns = zmns_all(ii) * mask .* exp(...
            -( ( (m_plot - m_all(ii)).^2 + (n_plot - n_all(ii)).^2 ) / ...
            exp_factor^2 ) );
        
        magnitude_plot_rmnc = magnitude_plot_rmnc + exp_add_rmnc;
        magnitude_plot_zmns = magnitude_plot_zmns + exp_add_zmns;
    end
    
end

figure;
subplot(1,2,1);
box on;
surf(m_array, n_array,log(abs( magnitude_plot_rmnc)), 'EdgeColor', 'Interp', ...
    'FaceColor', 'Interp');
xlabel('mpol');
ylabel('ntor');
title('RMNC');
colormap('jet');
view(2)
colorbar
caxis([-20 -2])

subplot(1,2,2);
box on;
surf(m_array, n_array, log(abs(magnitude_plot_zmns)), 'EdgeColor', 'Interp', ...
    'FaceColor', 'Interp');
xlabel('mpol');
ylabel('ntor');
title('ZMNS');
colormap('jet');
view(2)
colorbar
caxis([-20 -2])







