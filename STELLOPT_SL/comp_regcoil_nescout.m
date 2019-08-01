function figure_handle = comp_regcoil_nescout( data1, data2, m_max, n_max  )
%PLOT_REGCOIL_NESCOUT
% For plot the diagnostics output from stellopt/regcoil/nescout runs

if nargin < 4
    n_max = inf;
end
if nargin < 3
    m_max = inf;
end
AU_Blue = [12 35 64]/255;
AU_Orange = [242 101 34]/255;

%numzeta = 32;
%numtheta = 32;

numzeta = 64;
numtheta = 64;

%numzeta = 128;
%numtheta = 256;

%numzeta = 128;
%numtheta = 256;
numzetal = numzeta * data1.nfp;
zeta = linspace(0,2*pi/data1.nfp, numzeta);
zetal = linspace(0,2*pi, numzetal);
theta = linspace(0,2*pi,numtheta);

[surf_x_allfp_1, surf_y_allfp_1, surf_z_allfp_1] = make_surf(numzetal, ...
    numtheta, zetal, theta, data1, m_max, n_max);
[surf_x_allfp_2, surf_y_allfp_2, surf_z_allfp_2] = make_surf(numzetal, ...
    numtheta, zetal, theta, data2, m_max, n_max);

% distance between mesh/grid point locations
diff_1 = sqrt((surf_x_allfp_1 - surf_x_allfp_2).^2 + ...
    (surf_y_allfp_1 - surf_y_allfp_2).^2 + ...
    (surf_z_allfp_1 - surf_z_allfp_2).^2) .* ...
    sign((surf_z_allfp_2 - surf_z_allfp_1));
mymask = ones(size(diff_1));
mymask(find(surf_z_allfp_1 < 0)) = NaN;

%if 1
figure_handle = figure;

subplot(2,2,1); % superimpose both surfaces.
surf(surf_x_allfp_1, surf_y_allfp_1, surf_z_allfp_1, ...
    ones(size(surf_z_allfp_1)), 'EdgeColor', 'None', 'FaceColor', AU_Blue);
hold on;box on;
surf(surf_x_allfp_2, surf_y_allfp_2, surf_z_allfp_2, ...
    ones(size(surf_z_allfp_2)), 'EdgeColor', 'None', 'FaceColor', AU_Orange);
axis equal;axis tight;view(3);
title('Blue = Init. Surface, Orange = New Surface')
box on ; hold on
grid on

%figure
subplot(2,2,2);

for ii = 0:1.5:3
    zeta_ind = round(numzetal*ii/ (6 * data1.nfp) );
    if (zeta_ind == 0)
        zeta_ind = 1; % otherwise, you will try to access an illegal array index
    end
    %subplot(1,2,1); % superimpose both surfaces.
    R_plot = sqrt(surf_x_allfp_1(zeta_ind,:).^2 +...
        surf_y_allfp_1(zeta_ind,:).^2);
    Z_plot = surf_z_allfp_1(zeta_ind,:);
    plot(R_plot, Z_plot, 'Color', AU_Blue);
    hold on;box on;
    R_plot = sqrt(surf_x_allfp_2(zeta_ind,:).^2 +...
        surf_y_allfp_2(zeta_ind,:).^2);
    Z_plot = surf_z_allfp_2(zeta_ind,:);
    plot(R_plot, Z_plot, 'Color', AU_Orange);
    axis equal;axis tight;
    title('Blue = Init. Surface, Orange = New Surface')
    grid on;
    set(get(gca, 'Children'), 'Linewidth', 2);
    grid on
    
end


set(gcf, 'Position', [200 100 800 600]);

scale_factor = 100;

%figure
subplot(2,2,3);
hold on;box on;
surf(surf_x_allfp_2, surf_y_allfp_2, surf_z_allfp_2, ...
    scale_factor*diff_1.*mymask, 'EdgeColor', 'interp', 'FaceColor', 'interp');
colormap('jet');
%colorbar
hold on;box on;
axis equal;axis tight;view(3)
title('\Delta  [cm]')
grid on
%end


%figure
subplot(2,2,4);
hold on;box on;
surf(surf_x_allfp_2, surf_y_allfp_2, surf_z_allfp_2, ...
    scale_factor*diff_1.*mymask, 'EdgeColor', 'interp', 'FaceColor', 'interp');
colormap('jet');colorbar
title('\Delta  [cm]')
hold on;box on;
axis equal;axis tight;view(3)
%title('Difference in mesh points')
grid on
view(2)
%end


function [surf_x, surf_y, surf_z] = make_surf(numzetal, numtheta, ...
            zetal, theta, data, m_max, n_max)

[theta_mapallfp, zeta_mapallfp] = meshgrid(theta, zetal);

surf_x = 0 * theta_mapallfp;
surf_y = 0 * theta_mapallfp;
surf_z = 0 * theta_mapallfp;


for ii = 1:numzetal
    for jj = 1:numtheta
        for kk = 1:data.mn_surf
            m = data.xm_surf(kk);
            n = data.xn_surf(kk);
            % stellarator symmetry only for now...
            if ( (abs(m) <= m_max) && (abs(n) <= n_max) )
                rmnc = data.rmnc_surf(kk);
                zmns = data.zmns_surf(kk);
                angle = m*theta(jj) + n*zetal(ii);
                % angle = m*theta(jj) + data.nfp*zetal(ii);
                cosangle = cos(angle);
                sinangle = sin(angle);
                coszeta = cos(zetal(ii));
                sinzeta = sin(zetal(ii));
                surf_x(ii, jj) = surf_x(ii,jj) + ...
                    rmnc * cosangle * coszeta;
                surf_y(ii, jj) = surf_y(ii,jj) + ...
                    rmnc * cosangle * sinzeta;
                surf_z(ii, jj) = surf_z(ii,jj) + ...
                    zmns * sinangle;
            end
        end
    end
end

