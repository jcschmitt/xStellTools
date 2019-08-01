function [surf_R_allfp_1, surf_Z_allfp_1] = reFFT_nescin(data1, mmax, nmax)


numzeta = nmax*2;
numtheta = mmax*2;
numzetal = numzeta * data1.nfp;
zeta = linspace(0,2*pi/data1.nfp, numzeta);
zetal = linspace(0,2*pi, numzetal);
theta = linspace(0,2*pi,numtheta);

%[theta_map1fp, zeta_map1fp] = meshgrid(theta, zeta);
[theta_mapallfp, zeta_mapallfp] = meshgrid(theta, zetal);

%surf_x_1fp = 0 * theta_map1fp;
%surf_y_1fp = 0 * theta_map1fp;
%surf_z_1fp = 0 * theta_map1fp;
surf_R_allfp_1 = 0 * theta_mapallfp;
surf_Z_allfp_1 = 0 * theta_mapallfp;

% reconstruct the surface
for ii = 1:numzetal
    for jj = 1:numtheta
        for kk = 1:data1.mn_surf
            m = data1.xm_surf(kk);
            n = data1.xn_surf(kk);
            % stellarator symmetry only for now...
                rmnc = data1.rmnc_surf(kk);
                zmns = data1.zmns_surf(kk);
                angle = m*theta(jj) + n*zetal(ii);
                % angle = m*theta(jj) + data.nfp*zetal(ii);
                cosangle = cos(angle);
                sinangle = sin(angle);
                coszeta = cos(zetal(ii));
                sinzeta = sin(zetal(ii));
                surf_R_allfp_1(ii, jj) = surf_R_allfp_1(ii,jj) + ...
                    rmnc * cosangle;
                surf_Z_allfp_1(ii, jj) = surf_Z_allfp_1(ii,jj) + ...
                    zmns * sinangle;
        end
    end
end



