function [r_eff, Vprime, B, Bsquared, S11, S12, S21, S22] = ...
    generate_susceptance_matrix_v3_RHS(vmec_filename, num_surfs_in)
% function [r_eff, Vprime, B, Bsquared, S11, S12, S21, S22] = ...
%     generate_susceptance_matrix(vmec_filename, num_surfs_in)

% under construction.  working on speeding this up and making it work for
% the RHS output from VMEC.

DEBUG_MODE = 0;
PLOTSTUFF = 1;  % plot the output values
PLOTMORESTUFF = 1;  %  plot more details of the calculation

if nargin < 2
    stop
end

% load the variables

if strfind(vmec_filename, '.nc') % is it an .nc file?
    %     disp('Checking file');
    %     input_extension_vmec = getnc(file_vmec_output,'input_extension');
    disp('Opening .nc file');
    ncID = netcdf.open(vmec_filename,'nc_nowrite');
    %USED
    disp('Loading data from VMEC wout file');
    vmecVarID = netcdf.inqVarID(ncID,'ns');
    ns_vmec = double(netcdf.getVar(ncID,vmecVarID))
    
    vmecVarID = netcdf.inqVarID(ncID,'nfp');
    nfp_vmec = double(netcdf.getVar(ncID,vmecVarID));
    
    vmecVarID = netcdf.inqVarID(ncID,'phi');
    phitot_vmec = double(netcdf.getVar(ncID,vmecVarID))';
    
    vmecVarID = netcdf.inqVarID(ncID,'Aminor_p');
    a0 = double(netcdf.getVar(ncID,vmecVarID));
    
    vmecVarID = netcdf.inqVarID(ncID,'Rmajor_p');
    R0 = double(netcdf.getVar(ncID,vmecVarID));
    
    vmecVarID = netcdf.inqVarID(ncID,'mpol');
    mpol_vmec = double(netcdf.getVar(ncID,vmecVarID));
    
    vmecVarID = netcdf.inqVarID(ncID,'ntor');
    ntor_vmec = double(netcdf.getVar(ncID,vmecVarID));
    
    vmecVarID = netcdf.inqVarID(ncID,'mnmax');
    mnmax_vmec = double(netcdf.getVar(ncID,vmecVarID));
    
    vmecVarID = netcdf.inqVarID(ncID,'xn');
    xn_vmec_in = double(netcdf.getVar(ncID,vmecVarID))';
    %disp(strcat('<----Scaling xn_vmec by nfp=', num2str(nfp_vmec)));
    %xn_vmec = xn_vmec_in / nfp_vmec;
    xn_vmec = xn_vmec_in;
    
    vmecVarID = netcdf.inqVarID(ncID,'xm');
    xm_vmec = double(netcdf.getVar(ncID,vmecVarID))';
    
    vmecVarID = netcdf.inqVarID(ncID,'rmnc');
    rmnc_vmec = double(netcdf.getVar(ncID,vmecVarID))';
    
    vmecVarID = netcdf.inqVarID(ncID,'zmns');
    zmns_vmec = double(netcdf.getVar(ncID,vmecVarID))';
    
    vmecVarID = netcdf.inqVarID(ncID,'lmns');
    lmns_vmec = double(netcdf.getVar(ncID,vmecVarID))';
    
    vmecVarID = netcdf.inqVarID(ncID,'iotaf');
    iota_vmec = double(netcdf.getVar(ncID,vmecVarID))';


else
    warning('Not gonna deal w/ old ASCII VMEC files.');
end

disp('Done loading data from VMEC wout file');
% a0 = 0.121499015773253
% R0 = 1.21136371077895
tor_nodes_per_period = ntor_vmec * 8 + 1
pol_nodes = mpol_vmec * 8 + 1
radial_surfs = num_surfs_in
% radial_surfs = ns_vmec

% set up radial array
s_vmec = ((1:ns_vmec)-1)/(ns_vmec-1);
rho_vmec = sqrt(s_vmec);
% r_vmec = rho_vmec * a0;

s_local = ((1:radial_surfs)-1)/(radial_surfs-1);
rho_local = sqrt(s_local);
r_local = rho_local * a0;

% set up a 'half mesh', between the nodes of the rho_vmec mesh, along with
% the endpoints 0 and 1.
rho_halfmesh_vmec = [0 sqrt(.5*(s_vmec(2:end)+s_vmec(1:(end-1)))) 1];
 
% iota_vmec(1) = iota_vmec(2);
iota_spline_d = csaps(rho_vmec, iota_vmec);
iota_spline_d = fnxtr(iota_spline_d);
% [a,b,c] = v_spline(ns_vmec, rho_vmec, iota_vmec);
% iota_spline_a = a;
% iota_spline_b = b;
% iota_spline_c = c;

for ii = 1:radial_surfs
    iota_local(ii) = ppval(iota_spline_d, rho_local(ii));
    %     iota_local(ii) = w_spline(ns_vmec, rho_local(ii), ...
    %         rho_vmec, iota_vmec, iota_spline_a, ...
    %         iota_spline_b, iota_spline_c, 0);
end

disp('Radial grids calculated and iota spline coeffs calculated');
disp('Normalizing lambda, rmn, and zmn terms and generating spline cooefficients');

% normalize lambda to rho^m for better axial resolution
for ii = 1:mnmax_vmec
    %     disp(['Surface # ' num2str(ii)]);
    rho_halfmesh_m_norm = rho_halfmesh_vmec.^abs(xm_vmec(ii));
    rho_m_norm = rho_vmec.^abs(xm_vmec(ii));
    
    % normalize the bulk of the lambda terms
    lmns_normalized(ii,2:ns_vmec) = lmns_vmec(2:ns_vmec,ii) ./ ...
        rho_halfmesh_m_norm(2:(end-1))';
    %     lmns_normalized(2:ns_vmec) = lmns_vmec(2:ns_vmec,ii);
    
    % extrapolate to the edges of the halfmesh radial grid (the 1.5 & .5
    % factors are because of the grid space differences
    lmns_normalized(ii,1) = 1.5 * lmns_normalized(ii,2) - ...
        .5 * lmns_normalized(ii,3);
    lmns_normalized(ii,ns_vmec+1) = 1.5 * lmns_normalized(ii,ns_vmec) - ...
        .5 * lmns_normalized(ii,ns_vmec-1);
    
    % interpolate the lambda terms back to the full grid
    lmns_local(ii,:) = interp1(rho_halfmesh_vmec, lmns_normalized(ii,:), rho_local);
    
    % normalize the bulk of the rmnc and zmns terms
    rmnc_normalized(ii, 2:ns_vmec) = rmnc_vmec(2:ns_vmec, ii)' ./ ...
        rho_m_norm(2:end);
    zmns_normalized(ii, 2:ns_vmec) = zmns_vmec(2:ns_vmec, ii)' ./ ...
        rho_m_norm(2:end);
    
    % extrapolate to the axis for the normalized rmnc and zmns terms
    % Note by JCS: Here, the extrapolation is 'normal'.
    rmnc_normalized(ii, 1) = 2*rmnc_normalized(ii, 2) - rmnc_normalized(ii, 3);
    zmns_normalized(ii, 1) = 2*zmns_normalized(ii, 2) - zmns_normalized(ii, 3);
    % and get the spline coeffiecients for r,z, and radial deriviatives
    %     rmnc_spline_d(ii) = csaps(rho_vmec, rmnc_normalized(ii, :), csap_val);
    rmnc_spline_d(ii) = csaps(rho_vmec, rmnc_normalized(ii, :));
    %     zmns_spline_d(ii) = csaps(rho_vmec, zmns_normalized(ii, :), csap_val);
    zmns_spline_d(ii) = csaps(rho_vmec, zmns_normalized(ii, :));
    
    for jj = 1:radial_surfs
        rmnc_local(ii,jj) = ppval(rmnc_spline_d(ii), rho_local(jj));
        zmns_local(ii,jj) = ppval(zmns_spline_d(ii), rho_local(jj));
                
        if jj == 1
            drmnc_local(ii,jj) = 0;
            dzmns_local(ii,jj) = 0;
        else
            drmnc_local(ii,jj) = (ppval(rmnc_spline_d(ii), rho_local(jj)+1e-3) - ...
                ppval(rmnc_spline_d(ii), rho_local(jj)-1e-3)) / 2e-3;
            dzmns_local(ii,jj) = (ppval(zmns_spline_d(ii), rho_local(jj)+1e-3) - ...
                ppval(zmns_spline_d(ii), rho_local(jj)-1e-3)) / 2e-3;
        end
        
    end
end
disp('Done normalizing and finding spline coefficients')

if PLOTMORESTUFF
    % find the index of the m=1, n=0 component
    ind_m1n0_term = intersect(find(xm_vmec == -1), find(xn_vmec==0));
    ind_m0n1_term = intersect(find(xm_vmec == 0), find(xn_vmec==1));
    ind_m1n1_term = intersect(find(xm_vmec == -1), find(xn_vmec==1));
    ind_m0n0_term = intersect(find(xm_vmec == 0), find(xn_vmec==0));
    
    x_dense_1 = linspace(0, .25, 1000);
    
    
    for ind2plot = [ind_m1n0_term ind_m0n1_term ind_m1n1_term ind_m0n0_term]
        
        for kk = 1:length(x_dense_1)
            %             rmnc_dense_1(1,kk) = w_spline(ns_vmec, x_dense_1(kk), ...
            %                 rho_vmec, rmnc_normalized(ind2plot,:), rmnc_spline_a(ind2plot,:), ...
            %                 rmnc_spline_b(ind2plot,:), rmnc_spline_c(ind2plot,:), 0);
            %             zmns_dense_1(1,kk) = w_spline(ns_vmec, x_dense_1(kk), ...
            %                 rho_vmec, zmns_normalized(ind2plot,:), zmns_spline_a(ind2plot,:), ...
            %                 zmns_spline_b(ind2plot,:), zmns_spline_c(ind2plot,:), 0);
            %             drmnc_dense_1(1,kk) = w_spline(ns_vmec, x_dense_1(kk), ...
            %                 rho_vmec, rmnc_normalized(ind2plot, :), rmnc_spline_a(ind2plot,:), ...
            %                 rmnc_spline_b(ind2plot,:), rmnc_spline_c(ind2plot,:), 1);
            %             dzmns_dense_1(1,kk) = w_spline(ns_vmec, x_dense_1(kk), ...
            %                 rho_vmec, zmns_normalized(ind2plot, :), zmns_spline_a(ind2plot,:), ...
            %                 zmns_spline_b(ind2plot,:), zmns_spline_c(ind2plot,:), 1);
            rmnc_dense_1(1,kk) = ppval(rmnc_spline_d(ind2plot), x_dense_1(kk));
            zmns_dense_1(1,kk) = ppval(zmns_spline_d(ind2plot), x_dense_1(kk));
            drmnc_dense_1(1,kk) =(ppval(rmnc_spline_d(ind2plot), x_dense_1(kk)+1e-3) - ...
                ppval(rmnc_spline_d(ind2plot), x_dense_1(kk)-1e-3)) / 2e-3;
            dzmns_dense_1(1,kk) = (ppval(zmns_spline_d(ind2plot), x_dense_1(kk)+1e-3) - ...
                ppval(zmns_spline_d(ind2plot), x_dense_1(kk)-1e-3)) / 2e-3;
        end
        
        norm_power = abs(xm_vmec(ind2plot));
        figure
        subplot(3,2,1)
        plot(rho_vmec, rmnc_vmec(:,ind2plot)'./(rho_vmec.^norm_power) ,'o', rho_vmec, rmnc_normalized(ind2plot,:)', 'rx');
        hold on
        plot(rho_local, rmnc_local(ind2plot,:), 'ks')
        
        legend('rmnc', 'normalized', 'local');
        title(['mode# ' num2str(ind2plot)]);
        subplot(3,2,3)
        plot(rho_vmec, zmns_vmec(:,ind2plot)'./(rho_vmec.^norm_power) ,'o', rho_vmec, zmns_normalized(ind2plot,:)', 'rx');
        hold on
        plot(rho_local, zmns_local(ind2plot,:), 'ks')
        legend('zmns', 'normalized', 'local');
        subplot(3,2,5)
        plot(rho_halfmesh_vmec(2:end-1), lmns_vmec(2:ns_vmec,ind2plot)'./(rho_halfmesh_vmec(2:end-1).^norm_power) ,'o');
        hold on
        plot(rho_halfmesh_vmec, lmns_normalized(ind2plot,:)', 'rx');
        hold on
        plot(rho_local, lmns_local(ind2plot,:), 'ks')
        legend('lmns', 'normalized', 'local');
        
        subplot(3,2,2)
        plot(rho_vmec, rmnc_vmec(:,ind2plot)' ,'o', rho_vmec, rmnc_normalized(ind2plot,:)', 'rx');
        hold on
        plot(rho_local, rmnc_local(ind2plot,:), 'ks')
        legend('rmnc', 'normalized', 'local');
        subplot(3,2,4)
        plot(rho_vmec, zmns_vmec(:,ind2plot)' ,'o', rho_vmec, zmns_normalized(ind2plot,:)', 'rx');
        hold on
        plot(rho_local, zmns_local(ind2plot,:), 'ks')
        legend('zmns', 'normalized', 'local');
        subplot(3,2,6)
        plot(rho_halfmesh_vmec(2:end-1), lmns_vmec(2:ns_vmec,ind2plot)' ,'o');
        hold on
        plot(rho_halfmesh_vmec, lmns_normalized(ind2plot,:)', 'rx');
        hold on
        plot(rho_local, lmns_local(ind2plot,:), 'ks')
        legend('lmns', 'normalized', 'local');
        
        figure
        subplot(2,2,1)
        plot(rho_vmec, rmnc_vmec(:,ind2plot)'./(rho_vmec.^norm_power) ,'o', rho_vmec, rmnc_normalized(ind2plot,:)', 'rx');
        hold on
        plot(rho_local, rmnc_local(ind2plot,:), 'ks')
        plot(x_dense_1, rmnc_dense_1, 'y.');
        legend('rmnc', 'normalized', 'local');
        title(['mode# ' num2str(ind2plot)]);
        
        subplot(2,2,3)
        plot(rho_local, drmnc_local(ind2plot,:), 'ks')
        hold on
        plot(x_dense_1, drmnc_dense_1, 'y.');
        legend('drmnc');
        
        subplot(2,2,2)
        plot(rho_vmec, zmns_vmec(:,ind2plot)'./(rho_vmec.^norm_power) ,'o', rho_vmec, zmns_normalized(ind2plot,:)', 'rx');
        hold on
        plot(rho_local, zmns_local(ind2plot,:), 'ks')
        plot(x_dense_1, zmns_dense_1, 'y.');
        legend('zmns', 'normalized', 'local');
        
        subplot(2,2,4)
        plot(rho_local, dzmns_local(ind2plot,:), 'ks')
        hold on
        plot(x_dense_1, dzmns_dense_1, 'y.');
        
        legend('dzmns');
        
    end
end

disp('Beginning the metric coefficient calculation on the MHD grid.')

% set up angle deltas
d_phi = 2 * pi / (nfp_vmec * (tor_nodes_per_period-1));
d_theta = 2 * pi / (pol_nodes-1);

% integration (weighting) factors for theta and some averaging
% theta_weights = [1 4 2 4 2 ... 4 2 4 (2) 4 1];
% integration (weighting) factors for phi
% phi_weights = [1 4 2 4 2 ... 4 2 4 (2) 4 1];
theta_weights = ones(1,pol_nodes);
phi_weights = ones(1,tor_nodes_per_period);
if 0
    theta_weights(3:2:(end-2)) = 2;
    theta_weights(2:2:(end-1)) = 4;
    theta_weights(end-1) = 4;  % just in case there is an even # of nodes/period
    theta_norm = 3 * (pol_nodes-1);
    phi_weights(3:2:(end-2)) = 2;
    phi_weights(2:2:(end-1)) = 4;
    phi_weights(end-1) = 4;  % just in case there is an even # of nodes/period
    phi_norm = 3 * (tor_nodes_per_period-1);
else
    theta_weights(end) = 0;
    theta_norm = pol_nodes-1;
    phi_weights(end) = 0;
    phi_norm = tor_nodes_per_period-1;
end

for ii = [2:radial_surfs] % loop over each surface, calculate stuff
    disp(['Surface # ' num2str(ii) '/' num2str(radial_surfs)]);
    % loop over toroidal angle
    for jj = 1:tor_nodes_per_period
        %         disp(['Tor node # ' num2str(jj) '/' num2str(tor_nodes_per_period) ]);
        %         pause(.01);
        phi = (jj-1) * d_phi;
        % loop over poloidal angle
        for kk = 1:pol_nodes
            %             disp(['Pol node # ' num2str(kk) '/' num2str(pol_nodes) ]);
            theta = (kk-1) * d_theta;
            % find r, z, and deriviatisve w.r.t. rho, theta and phi
            % find cylindrical components of magnetic field, and jacobian
            [r(ii,jj,kk), z(ii,jj,kk), r_rho(ii,jj,kk), ...
                r_theta(ii,jj,kk), r_phi(ii,jj,kk), z_rho(ii,jj,kk), ...
                z_theta(ii,jj,kk), z_phi(ii,jj,kk), ...
                lambda_phi(ii,jj,kk), lambda_theta(ii,jj,kk), ...
                b_r(ii,jj,kk), b_z(ii,jj,kk), b_tor(ii,jj,kk), ...
                sqrt_g(ii,jj,kk)] = ...
                tragrz(rho_local(ii), theta, phi, mnmax_vmec, ...
                xm_vmec, xn_vmec, rmnc_local(:,ii), zmns_local(:,ii), ...
                drmnc_local(:,ii), dzmns_local(:,ii), lmns_local(:,ii), ...
                iota_local(ii), phitot_vmec(end));

            % some geometric quantities that I need
            % this is V' = dV/drho ??
            %             geom_1(ii,jj,kk) = r * (r_theta * z_rho - r_rho * z_theta) / a0;
            geom_1_a(ii,jj,kk) = r(ii,jj,kk) * ...
                (r_theta(ii,jj,kk) * z_rho(ii,jj,kk) - ...
                r_rho(ii,jj,kk) * z_theta(ii,jj,kk));
            
            % this is what THRIFT used
            g1_thrift(ii,jj,kk) = geom_1_a(ii,jj,kk) / a0;
            g_inv_thrift(ii,jj,kk) = 1/g1_thrift(ii,jj,kk);
            %             g2_thrift(ii,jj,kk) = r(ii,jj,kk).^2 * ...
            %                 (r_theta(ii,jj,kk).^2 + z_theta.^2) / g1_thrift(ii,jj,kk);
            
            % These are my Jacobians in rho- and r- space
            g1_rho_tlite(ii,jj,kk) = geom_1_a(ii,jj,kk) *4*pi*pi;
            g1_r_tlite(ii,jj,kk) = geom_1_a(ii,jj,kk) *4*pi*pi / a0;
            
            g1_method_3(ii,jj,kk) = geom_1_a(ii,jj,kk);

            % quantities related to the susceptance matrix
            s11_a_thrift(ii,jj,kk) = (r_theta(ii,jj,kk).^2 + z_theta(ii,jj,kk).^2) / ...
                g1_thrift(ii,jj,kk);
            s12_a_thrift(ii,jj,kk) = (r_theta(ii,jj,kk) * r_phi(ii,jj,kk) +...
                z_theta(ii,jj,kk) * z_phi(ii,jj,kk)) * ...
                (1 + lambda_theta(ii,jj,kk)) / g1_thrift(ii,jj,kk) - ...
                s11_a_thrift(ii,jj,kk) * lambda_phi(ii,jj,kk);
            s21_a_thrift(ii,jj,kk) = (r_theta(ii,jj,kk) * r_phi(ii,jj,kk) + ...
                z_theta(ii,jj,kk) * z_phi(ii,jj,kk)) / g1_thrift(ii,jj,kk);
            s22_a_thrift(ii,jj,kk) = (r(ii,jj,kk)^2 + r_phi(ii,jj,kk)^2 + ...
                z_phi(ii,jj,kk)^2) * (1 + lambda_theta(ii,jj,kk)) / ...
                g1_thrift(ii,jj,kk) - s21_a_thrift(ii,jj,kk) * lambda_phi(ii,jj,kk);
            
            s11_a_rho_tlite(ii,jj,kk) = (r_theta(ii,jj,kk).^2 + z_theta(ii,jj,kk).^2) / ...
                g1_rho_tlite(ii,jj,kk);
            s12_a_rho_tlite(ii,jj,kk) = (r_theta(ii,jj,kk) * r_phi(ii,jj,kk) +...
                z_theta(ii,jj,kk) * z_phi(ii,jj,kk)) * ...
                (1 + lambda_theta(ii,jj,kk)) / g1_rho_tlite(ii,jj,kk) - ...
                s11_a_rho_tlite(ii,jj,kk) * lambda_phi(ii,jj,kk);
            s21_a_rho_tlite(ii,jj,kk) = (r_theta(ii,jj,kk) * r_phi(ii,jj,kk) + ...
                z_theta(ii,jj,kk) * z_phi(ii,jj,kk)) / g1_rho_tlite(ii,jj,kk);
            s22_a_rho_tlite(ii,jj,kk) = (r(ii,jj,kk)^2 + r_phi(ii,jj,kk)^2 + ...
                z_phi(ii,jj,kk)^2) * (1 + lambda_theta(ii,jj,kk)) / ...
                g1_rho_tlite(ii,jj,kk) - s21_a_rho_tlite(ii,jj,kk) * lambda_phi(ii,jj,kk);
            
            s11_a_r_tlite(ii,jj,kk) = (r_theta(ii,jj,kk).^2 + z_theta(ii,jj,kk).^2) / ...
                g1_r_tlite(ii,jj,kk);
            s12_a_r_tlite(ii,jj,kk) = (r_theta(ii,jj,kk) * r_phi(ii,jj,kk) +...
                z_theta(ii,jj,kk) * z_phi(ii,jj,kk)) * ...
                (1 + lambda_theta(ii,jj,kk)) / g1_r_tlite(ii,jj,kk) - ...
                s11_a_r_tlite(ii,jj,kk) * lambda_phi(ii,jj,kk);
            s21_a_r_tlite(ii,jj,kk) = (r_theta(ii,jj,kk) * r_phi(ii,jj,kk) + ...
                z_theta(ii,jj,kk) * z_phi(ii,jj,kk)) / g1_r_tlite(ii,jj,kk);
            s22_a_r_tlite(ii,jj,kk) = (r(ii,jj,kk)^2 + r_phi(ii,jj,kk)^2 + ...
                z_phi(ii,jj,kk)^2) * (1 + lambda_theta(ii,jj,kk)) / ...
                g1_r_tlite(ii,jj,kk) - s21_a_r_tlite(ii,jj,kk) * lambda_phi(ii,jj,kk);
            

            s11_a_method_3(ii,jj,kk) = (r_theta(ii,jj,kk).^2 + z_theta(ii,jj,kk).^2) / ...
                g1_method_3(ii,jj,kk).^2;
            s12_a_method_3(ii,jj,kk) = (r_theta(ii,jj,kk) * r_phi(ii,jj,kk) +...
                z_theta(ii,jj,kk) * z_phi(ii,jj,kk)) * ...
                (1 + lambda_theta(ii,jj,kk)) / g1_method_3(ii,jj,kk).^2 - ...
                s11_a_thrift(ii,jj,kk) * lambda_phi(ii,jj,kk);
            s21_a_method_3(ii,jj,kk) = (r_theta(ii,jj,kk) * r_phi(ii,jj,kk) + ...
                z_theta(ii,jj,kk) * z_phi(ii,jj,kk)) / g1_method_3(ii,jj,kk).^2;
            s22_a_method_3(ii,jj,kk) = (r(ii,jj,kk)^2 + r_phi(ii,jj,kk)^2 + ...
                z_phi(ii,jj,kk)^2) * (1 + lambda_theta(ii,jj,kk)) / ...
                g1_method_3(ii,jj,kk).^2 - s21_a_thrift(ii,jj,kk) * lambda_phi(ii,jj,kk);

            % for <B^2> and <B>
            bsq(ii,jj,kk) = b_r(ii,jj,kk)^2 + b_z(ii,jj,kk)^2 + ...
                b_tor(ii,jj,kk)^2;
            b1(ii,jj,kk) = sqrt(bsq(ii,jj,kk));
            g6_thrift(ii,jj,kk) = g1_thrift(ii,jj,kk) * ...
                bsq(ii,jj,kk);
            g5_thrift(ii,jj,kk) = g1_thrift(ii,jj,kk) * ...
                b1(ii,jj,kk);
        end % loop over poloidal nodes
        
        
        % average over the poloidal angle.
        for kk = 1:pol_nodes
            g1_pol_ave_thrift(kk) = g1_thrift(ii,jj,kk);
            g_inv_pol_ave_inv(kk) = g_inv_thrift(ii,jj,kk);
            g1_pol_ave_rho_tlite(kk) = g1_rho_tlite(ii,jj,kk);
            g1_pol_ave_r_tlite(kk) = g1_r_tlite(ii,jj,kk);
            g5_pol_ave_thrift(kk) = g5_thrift(ii,jj,kk);
            g6_pol_ave_thrift(kk) = g6_thrift(ii,jj,kk);
            s11_pol_ave_thrift(kk) = s11_a_thrift(ii,jj,kk);
            s12_pol_ave_thrift(kk) = s12_a_thrift(ii,jj,kk);
            s21_pol_ave_thrift(kk) = s21_a_thrift(ii,jj,kk);
            s22_pol_ave_thrift(kk) = s22_a_thrift(ii,jj,kk);
            s11_pol_ave_rho_tlite(kk) = s11_a_rho_tlite(ii,jj,kk);
            s12_pol_ave_rho_tlite(kk) = s12_a_rho_tlite(ii,jj,kk);
            s21_pol_ave_rho_tlite(kk) = s21_a_rho_tlite(ii,jj,kk);
            s22_pol_ave_rho_tlite(kk) = s22_a_rho_tlite(ii,jj,kk);
            s11_pol_ave_r_tlite(kk) = s11_a_r_tlite(ii,jj,kk);
            s12_pol_ave_r_tlite(kk) = s12_a_r_tlite(ii,jj,kk);
            s21_pol_ave_r_tlite(kk) = s21_a_r_tlite(ii,jj,kk);
            s22_pol_ave_r_tlite(kk) = s22_a_r_tlite(ii,jj,kk);
            s11_pol_ave_method_3(kk) = s11_a_method_3(ii,jj,kk);
            s12_pol_ave_method_3(kk) = s12_a_method_3(ii,jj,kk);
            s21_pol_ave_method_3(kk) = s21_a_method_3(ii,jj,kk);
            s22_pol_ave_method_3(kk) = s22_a_method_3(ii,jj,kk);
        end
        
        g1_theta_avg(ii,jj) = sum(g1_pol_ave_thrift .* theta_weights) / theta_norm;
        g_inv_theta_avg(ii,jj) = sum(g_inv_pol_ave_inv .* theta_weights) / theta_norm;
        
        g5_theta_avg(ii,jj) = sum(g5_pol_ave_thrift .* theta_weights) / theta_norm;
        g6_theta_avg(ii,jj) = sum(g6_pol_ave_thrift .* theta_weights) / theta_norm;
        Jac_rho_theta_avg(ii,jj) = sum(g1_pol_ave_rho_tlite .* theta_weights) / theta_norm;
        Jac_r_theta_avg(ii,jj) = sum(g1_pol_ave_r_tlite .* theta_weights) / theta_norm;
        
        s11_theta_avg_thrift(ii,jj) = sum(s11_pol_ave_thrift .* theta_weights) / theta_norm;
        s12_theta_avg_thrift(ii,jj) = sum(s12_pol_ave_thrift .* theta_weights) / theta_norm;
        s21_theta_avg_thrift(ii,jj) = sum(s21_pol_ave_thrift .* theta_weights) / theta_norm;
        s22_theta_avg_thrift(ii,jj) = sum(s22_pol_ave_thrift .* theta_weights) / theta_norm;
        s11_theta_avg_rho_tlite(ii,jj) = sum(s11_pol_ave_rho_tlite .* theta_weights) / theta_norm;
        s12_theta_avg_rho_tlite(ii,jj) = sum(s12_pol_ave_rho_tlite .* theta_weights) / theta_norm;
        s21_theta_avg_rho_tlite(ii,jj) = sum(s21_pol_ave_rho_tlite .* theta_weights) / theta_norm;
        s22_theta_avg_rho_tlite(ii,jj) = sum(s22_pol_ave_rho_tlite .* theta_weights) / theta_norm;
        s11_theta_avg_r_tlite(ii,jj) = sum(s11_pol_ave_r_tlite .* theta_weights) / theta_norm;
        s12_theta_avg_r_tlite(ii,jj) = sum(s12_pol_ave_r_tlite .* theta_weights) / theta_norm;
        s21_theta_avg_r_tlite(ii,jj) = sum(s21_pol_ave_r_tlite .* theta_weights) / theta_norm;
        s22_theta_avg_r_tlite(ii,jj) = sum(s22_pol_ave_r_tlite .* theta_weights) / theta_norm;
        s11_theta_avg_method_3(ii,jj) = sum(s11_pol_ave_method_3 .* theta_weights) / theta_norm;
        s12_theta_avg_method_3(ii,jj) = sum(s12_pol_ave_method_3 .* theta_weights) / theta_norm;
        s21_theta_avg_method_3(ii,jj) = sum(s21_pol_ave_method_3 .* theta_weights) / theta_norm;
        s22_theta_avg_method_3(ii,jj) = sum(s22_pol_ave_method_3 .* theta_weights) / theta_norm;
    end % looping of tor nodes
    
    % average over the toroidal angle.
    g1_phitheta_avg(ii) = sum(g1_theta_avg(ii,:) .* phi_weights) / phi_norm;
    g_inv_phitheta_avg(ii) = sum(g_inv_theta_avg(ii,:) .* phi_weights) / phi_norm;
    g5_phitheta_avg(ii) = sum(g5_theta_avg(ii,:) .* phi_weights) / phi_norm;
    g6_phitheta_avg(ii) = sum(g6_theta_avg(ii,:) .* phi_weights) / phi_norm;
    Jac_rho_phitheta_avg(ii) = sum(Jac_rho_theta_avg(ii,:) .* phi_weights) / phi_norm;
    Jac_r_phitheta_avg(ii) = sum(Jac_r_theta_avg(ii,:) .* phi_weights) / phi_norm;
    s11_phitheta_avg_thrift(ii) = sum(s11_theta_avg_thrift(ii,:) .* phi_weights) / phi_norm;
    s12_phitheta_avg_thrift(ii) = sum(s12_theta_avg_thrift(ii,:) .* phi_weights) / phi_norm;
    s21_phitheta_avg_thrift(ii) = sum(s21_theta_avg_thrift(ii,:) .* phi_weights) / phi_norm;
    s22_phitheta_avg_thrift(ii) = sum(s22_theta_avg_thrift(ii,:) .* phi_weights) / phi_norm;
    s11_phitheta_avg_rho_tlite(ii) = sum(s11_theta_avg_rho_tlite(ii,:) .* phi_weights) / phi_norm;
    s12_phitheta_avg_rho_tlite(ii) = sum(s12_theta_avg_rho_tlite(ii,:) .* phi_weights) / phi_norm;
    s21_phitheta_avg_rho_tlite(ii) = sum(s21_theta_avg_rho_tlite(ii,:) .* phi_weights) / phi_norm;
    s22_phitheta_avg_rho_tlite(ii) = sum(s22_theta_avg_rho_tlite(ii,:) .* phi_weights) / phi_norm;
    s11_phitheta_avg_r_tlite(ii) = sum(s11_theta_avg_r_tlite(ii,:) .* phi_weights) / phi_norm;
    s12_phitheta_avg_r_tlite(ii) = sum(s12_theta_avg_r_tlite(ii,:) .* phi_weights) / phi_norm;
    s21_phitheta_avg_r_tlite(ii) = sum(s21_theta_avg_r_tlite(ii,:) .* phi_weights) / phi_norm;
    s22_phitheta_avg_r_tlite(ii) = sum(s22_theta_avg_r_tlite(ii,:) .* phi_weights) / phi_norm;
    s11_phitheta_avg_method_3(ii) = sum(s11_theta_avg_method_3(ii,:) .* phi_weights) / phi_norm;
    s12_phitheta_avg_method_3(ii) = sum(s12_theta_avg_method_3(ii,:) .* phi_weights) / phi_norm;
    s21_phitheta_avg_method_3(ii) = sum(s21_theta_avg_method_3(ii,:) .* phi_weights) / phi_norm;
    s22_phitheta_avg_method_3(ii) = sum(s22_theta_avg_method_3(ii,:) .* phi_weights) / phi_norm;
    
    g1_b_thrift(ii) = 1 * g1_phitheta_avg(ii);
    g_inv_b_thrift(ii) = 1 * g_inv_phitheta_avg(ii);
    g5_b_thrift(ii) = 1 * g5_phitheta_avg(ii);
    g6_b_thrift(ii) = 1 * g6_phitheta_avg(ii);
    Vp_r_thrift(ii) = g1_b_thrift(ii) * (4*pi*pi);
    Jac_rho_b_tlite(ii) = 1 * Jac_rho_phitheta_avg(ii);
    Jac_r_b_tlite(ii) = 1 * Jac_r_phitheta_avg(ii);
    
    % normalized the values according to their radial behavior
    s11_b_thrift(ii) = 1 * s11_phitheta_avg_thrift(ii) / (r_local(ii));
    s12_b_thrift(ii) = 1 * s12_phitheta_avg_thrift(ii) / (r_local(ii));
    s21_b_thrift(ii) = 1 * s21_phitheta_avg_thrift(ii) / (r_local(ii));
    s22_b_thrift(ii) = 1 * s22_phitheta_avg_thrift(ii) * (r_local(ii));
    s11_b_rho_tlite(ii) = 1 * s11_phitheta_avg_rho_tlite(ii) / (rho_local(ii));
    s12_b_rho_tlite(ii) = 1 * s12_phitheta_avg_rho_tlite(ii) / (rho_local(ii));
    s21_b_rho_tlite(ii) = 1 * s21_phitheta_avg_rho_tlite(ii) / (rho_local(ii));
    s22_b_rho_tlite(ii) = 1 * s22_phitheta_avg_rho_tlite(ii) * (rho_local(ii));
    s11_b_r_tlite(ii) = 1 * s11_phitheta_avg_r_tlite(ii) / (r_local(ii));
    s12_b_r_tlite(ii) = 1 * s12_phitheta_avg_r_tlite(ii) / (r_local(ii));
    s21_b_r_tlite(ii) = 1 * s21_phitheta_avg_r_tlite(ii) / (r_local(ii));
    s22_b_r_tlite(ii) = 1 * s22_phitheta_avg_r_tlite(ii) * (r_local(ii));
    s11_b_method_3(ii) = 1 * s11_phitheta_avg_method_3(ii);
    s12_b_method_3(ii) = 1 * s12_phitheta_avg_method_3(ii);
    s21_b_method_3(ii) = 1 * s21_phitheta_avg_method_3(ii);
    s22_b_method_3(ii) = 1 * s22_phitheta_avg_method_3(ii);
end

disp('phi =')
phi

% extrapoloate to the axis
% geom_1_b(1) = geom_1_b(2) - rho_local(2) * (geom_1_b(3) - geom_1_b(2)) / ...
%     (rho_local(3) - rho_local(2));

g1_b_thrift(1) = 0;
g5_b_thrift(1) = g5_b_thrift(2);
g6_b_thrift(1) = g6_b_thrift(2);
Jac_rho_b_tlite(1) = 0;
Jac_r_b_tlite(1) = 0;
s11_b_thrift(1) = s11_b_thrift(2) - rho_local(2) * (s11_b_thrift(3) - s11_b_thrift(2)) / ...
    (rho_local(3) - rho_local(2));
s12_b_thrift(1) = s12_b_thrift(2) - rho_local(2) * (s12_b_thrift(3) - s12_b_thrift(2)) / ...
    (rho_local(3) - rho_local(2));
s21_b_thrift(1) = s21_b_thrift(2) - rho_local(2) * (s21_b_thrift(3) - s21_b_thrift(2)) / ...
    (rho_local(3) - rho_local(2));
s22_b_thrift(1) = s22_b_thrift(2) - rho_vmec(2) * (s22_b_thrift(3) - s22_b_thrift(2)) / ...
    (rho_local(3) - rho_local(2));
s11_b_rho_tlite(1) = s11_b_rho_tlite(2) - rho_local(2) * (s11_b_rho_tlite(3) - s11_b_rho_tlite(2)) / ...
    (rho_local(3) - rho_local(2));
s12_b_rho_tlite(1) = s12_b_rho_tlite(2) - rho_local(2) * (s12_b_rho_tlite(3) - s12_b_rho_tlite(2)) / ...
    (rho_local(3) - rho_local(2));
s21_b_rho_tlite(1) = s21_b_rho_tlite(2) - rho_local(2) * (s21_b_rho_tlite(3) - s21_b_rho_tlite(2)) / ...
    (rho_local(3) - rho_local(2));
s22_b_rho_tlite(1) = s22_b_rho_tlite(2) - rho_vmec(2) * (s22_b_rho_tlite(3) - s22_b_rho_tlite(2)) / ...
    (rho_local(3) - rho_local(2));
s11_b_r_tlite(1) = s11_b_r_tlite(2) - rho_local(2) * (s11_b_r_tlite(3) - s11_b_r_tlite(2)) / ...
    (rho_local(3) - rho_local(2));
s12_b_r_tlite(1) = s12_b_r_tlite(2) - rho_local(2) * (s12_b_r_tlite(3) - s12_b_r_tlite(2)) / ...
    (rho_local(3) - rho_local(2));
s21_b_r_tlite(1) = s21_b_r_tlite(2) - rho_local(2) * (s21_b_r_tlite(3) - s21_b_r_tlite(2)) / ...
    (rho_local(3) - rho_local(2));
s22_b_r_tlite(1) = s22_b_r_tlite(2) - rho_vmec(2) * (s22_b_r_tlite(3) - s22_b_r_tlite(2)) / ...
    (rho_local(3) - rho_local(2));
s11_b_method_3(1) = NaN;
s12_b_method_3(1) = NaN;
s21_b_method_3(1) = NaN;
s22_b_method_3(1) = NaN;

% dvol_r_thrift = volume of cell 'i' in m^3
dvol_r_thrift(1) = Vp_r_thrift(2) * r_local(2) * 0.5;
for ii = 2:(radial_surfs - 1)
    dvol_r_thrift(ii) = 0.5 * (Vp_r_thrift(ii) / r_local(ii) + ...
        Vp_r_thrift(ii+1) / r_local(ii+1) ) * ...
        0.5 * (r_local(ii+1).^2 - r_local(ii).^2);
end

for ii = 2:radial_surfs
    Bsquared(ii) = g6_b_thrift(ii) / g1_b_thrift(ii);
    B(ii) = g5_b_thrift(ii) / g1_b_thrift(ii);
end
Bsquared(1) = Bsquared(2);
B(1) = B(2);
% Bsquared(1) = NaN;
% B(1) = NaN;

% the volume of the plasma
sum(dvol_r_thrift)
trapz(rho_local, Jac_rho_b_tlite)
trapz(r_local, Jac_r_b_tlite)
trapz(r_local, Vp_r_thrift)

% these should all be '1' if the volume, R0, and a0 are the 'effective'
% volume, major and minor radius, respectively

mean(Vp_r_thrift(2:end) ./r_local(2:end)) / (R0 * 4*pi*pi)
mean(Jac_r_b_tlite(2:end) ./r_local(2:end)) / (R0 * 4*pi*pi)
mean(Jac_rho_b_tlite(2:end) ./rho_local(2:end)) / (R0 * 4*pi*pi*a0*a0)

g_inv_b_thrift;

S11_thrift = s11_b_thrift .* r_local;
S12_thrift = s12_b_thrift .* r_local;
S21_thrift = s21_b_thrift .* r_local;
S22_thrift = (s22_b_thrift ./ r_local);

S11_rho_tlite = s11_b_rho_tlite .* rho_local / (4*pi*pi);
S12_rho_tlite = s12_b_rho_tlite .* rho_local / (4*pi*pi);
S21_rho_tlite = s21_b_rho_tlite .* rho_local / (4*pi*pi);
S22_rho_tlite = (s22_b_rho_tlite ./ rho_local) / (4*pi*pi);

S11_r_tlite = s11_b_r_tlite .* r_local / (4*pi*pi);
S12_r_tlite = s12_b_r_tlite .* r_local / (4*pi*pi);
S21_r_tlite = s21_b_r_tlite .* r_local / (4*pi*pi);
S22_r_tlite = (s22_b_r_tlite ./ r_local) / (4*pi*pi);

S11_method_3 = Vp_r_thrift .* s11_b_method_3 / (4*pi^2);
S12_method_3 = Vp_r_thrift .* s12_b_method_3 / (4*pi^2);
S21_method_3 = Vp_r_thrift .* s21_b_method_3 / (4*pi^2);
S22_method_3 = Vp_r_thrift .* s22_b_method_3 / (4*pi^2);


if PLOTSTUFF
    disp('Making plots')
    title_w2 = texit(vmec_filename);
    figure
    surf(squeeze(g_inv_thrift(15,:,:))*Vp_r_thrift(15)/(R0*4*pi*pi))
    
    figure
    subplot(3,4,1);box on; hold on;%axis([0 1 0 .3]);
    plot(rho_local, S11_thrift, '.');
    subplot(3,4,2);box on; hold on;%axis([0 1 0 .3]);
    title(title_w2);
    plot(rho_local, S12_thrift, '.');
    subplot(3,4,3);box on; hold on;%axis([0 1 0 .3]);
    plot(rho_local, S21_thrift, '.');
    subplot(3,4,4);box on; hold on;%axis([0 1 0 40]);
    plot(rho_local, S22_thrift, '.');
    subplot(3,4,5);box on; hold on
    plot(rho_local, S11_rho_tlite, '.');
    plot(rho_local, S11_thrift/(16*pi^4*a0), 'ro');
    legend('S11_{\rho tlite}', 'S11_{thrift} / (16*pi^4a_0)');
    subplot(3,4,6);box on; hold on
    plot(rho_local, S12_rho_tlite, '.');
    subplot(3,4,7);box on; hold on
    plot(rho_local, S21_rho_tlite, '.');
    subplot(3,4,8);box on; hold on
    plot(rho_local, S22_rho_tlite, '.');
    subplot(3,4,9);box on; hold on
    plot(rho_local, S11_r_tlite, '.');
    plot(rho_local, S11_thrift/(16*pi^4), 'ro');
    legend('S11_{r tlite}', 'S11_{thrift} / (16*pi^4)');
    subplot(3,4,10);box on; hold on
    plot(rho_local, S12_r_tlite, '.');
    subplot(3,4,11);box on; hold on
    plot(rho_local, S21_r_tlite, '.');
    subplot(3,4,12);box on; hold on
    plot(rho_local, S22_r_tlite, '.');
    
    
    figure
    subplot(2,2,1);box on;hold on
    plot(rho_local,S11_thrift, 'b.', rho_local,S11_thrift./rho_local , 'b.--');
    legend('S11')
    subplot(2,2,2);box on;hold on
    plot(rho_local,S12_thrift, 'b.', rho_local,S12_thrift./rho_local , 'b.--');
    legend('S12')
    subplot(2,2,3);box on;hold on
    plot(rho_local,S21_thrift, 'b.', rho_local,S21_thrift./rho_local , 'b.--');
    legend('S21')
    xlabel('\rho')
    subplot(2,2,4);box on;hold on
    plot(rho_local,S22_thrift, 'b.', rho_local,S22_thrift.*rho_local , 'b.--');
    legend('S22')
    xlabel('\rho')

    figure
    subplot(2,2,1);box on;hold on
    plot(rho_local,S11_method_3, 'b.', rho_local,S11_method_3./rho_local , 'b.--');
    legend('S11 method 3')
    subplot(2,2,2);box on;hold on
    plot(rho_local,S12_method_3, 'b.', rho_local,S12_method_3./rho_local , 'b.--');
    legend('S12 method 3')
    subplot(2,2,3);box on;hold on
    plot(rho_local,S21_method_3, 'b.', rho_local,S21_method_3./rho_local , 'b.--');
    legend('S21 method 3')
    xlabel('\rho')
    subplot(2,2,4);box on;hold on
    plot(rho_local,S22_method_3, 'b.', rho_local,S22_method_3.*rho_local , 'b.--');
    legend('S22 method 3')
    xlabel('\rho')
    
    figure
    subplot(2,2,1);box on;hold on
    plot(rho_local,S11_method_3, 'b.', rho_local,S11_method_3./rho_local , 'b.--');
    legend('S11 method 3')
    subplot(2,2,2);box on;hold on
    plot(rho_local,S12_method_3, 'b.', rho_local,S12_method_3./rho_local , 'b.--');
    legend('S12 method 3')
    subplot(2,2,3);box on;hold on
    plot(rho_local,S21_method_3, 'b.', rho_local,S21_method_3./rho_local , 'b.--');
    legend('S21 method 3')
    xlabel('\rho')
    subplot(2,2,4);box on;hold on
    plot(rho_local,S22_method_3, 'b.', rho_local,S22_method_3.*rho_local , 'b.--');
    legend('S22 method 3')
    xlabel('\rho')

    figure;
    subplot(2,2,1);box on;hold on;
    title(title_w2);
    plot(rho_local, -S12_thrift./S11_thrift,'b.');
    plot(rho_local, -S12_rho_tlite./S11_rho_tlite,'ko');
    plot(rho_local, -S12_r_tlite./S11_r_tlite,'rx');
    plot(rho_local, iota_local, 'g+');
    legend('iota_{thrift}', 'iota_{tlite, rho}', 'iota_{tlite, r}', 'iota_{vmec}')
    subplot(2,2,2);box on;hold on
    plot(rho_local, iota_local - -S12_thrift./S11_thrift, '.');
    legend('\Delta iota (VMEC-THRIFT)')
    subplot(2,1,2);box on; hold on;
    plot(rho_local, Bsquared, '.');
    legend('\langleB^2\rangle');
    
    
    
    
    for ii = [radial_surfs]
        figure
        box on;hold on;title(['Surf #' num2str(ii)]);
        for jj = 1:tor_nodes_per_period
            
            for kk = 1:pol_nodes
                phi(jj,kk) = (jj-1) * d_phi;
                theta(jj,kk) = (kk-1) * d_theta;
                my_x(jj,kk) = r(ii,jj,kk) * cos(phi(jj,kk));
                my_y(jj,kk) = r(ii,jj,kk) * sin(phi(jj,kk));
                my_z(jj,kk) = z(ii,jj,kk);
                my_r_rho(jj,kk) = r_rho(ii,jj,kk);
                my_r_theta(jj,kk) = r_theta(ii,jj,kk);
                my_r_phi(jj,kk) = r_phi(ii,jj,kk);
                my_z_rho(jj,kk) = z_rho(ii,jj,kk);
                my_z_theta(jj,kk) = z_theta(ii,jj,kk);
                my_z_phi(jj,kk) = z_phi(ii,jj,kk);
                my_lambda_phi(jj,kk) = lambda_phi(ii,jj,kk);
                my_lambda_theta(jj,kk) = lambda_theta(ii,jj,kk);
                %                 my_geom_1(jj,kk) = geom_1(ii,jj,kk);
                %                 my_S11(jj,kk) = s11_a(ii,jj,kk);
                %                 my_S12(jj,kk) = s12_a(ii,jj,kk);
                %                 my_S21(jj,kk) = s21_a(ii,jj,kk);
                %                 my_S22(jj,kk) = s22_a(ii,jj,kk);
                
            end
            plot3(my_x(jj,:),my_y(jj,:),my_z(jj,:),'.--');
            
        end
        view([115 25]);
        axis equal;
        
    end
    
    if DEBUG_MODE
        disp('Returning to original directory')
        cd(my_pwd)
    end
end

disp('Done generating susceptance matrix');

if nargout == 0
    r_eff = 0;
else
    r_eff = r_local;
    Vprime = Vp_r_thrift;
    S11 = S11_thrift;
    S12 = S12_thrift;
    S21 = S21_thrift;
    S22 = S22_thrift;
end



function [r, z, r_rho, r_theta, r_phi, z_rho, z_theta, z_phi, ...
    lambda_phi, lambda_theta, b_r, b_z, b_tor, sqrt_g] = ...
    tragrz(rho, theta, phi, mnmax, xm, xn, ...
    rmnc_full, zmns_full, drmnc_full, dzmns_full, lmns_full, ...
    iota, torflux_LCFS)
%***********************************************************************
%tragrz finds r and z and their derivatives with respect to rho, theta,
%   and phi at the given rho, theta, and phi.  It also returns the
%   cylindrical components of the magnetic field and the jacobian.
%References:
%  Attenberger, Houlberg, Hirshman J Comp Phys 72 (1987) 435
%  W.A.Houlberg 3/98
%Input:
%  rho-radial flux coordinate (-)
%  theta-poloidal angle coordinate (radians)
%  phi-toroidal angle coordinate (radians)
%  mnmax-total number of modes
%  xn-toroidal mode #s
%  xm-poloidal mode #s
%  rmnc_spline_pp - rmnc spline coeffiecinets for each mode
%  zmns_spline_pp - zmns spline coeffiecinets for each mode
%  drmnc_full, dzmns_full, lmns_full, iota, torflux_LCFS
%Output:
%  r-radius of point from major axis (m)
%  z-distance of point from midplane (m)
%  r_rho-dr/dx (m)
%  r_theta-dr/dtheta (m/radians)
%  r_phi-dr/dphi (m/radians)
%  z_rho-dz/dx (m)
%  z_theta-dz/dtheta (m/radians)
%  z_phi-dz/dphi (m/radians)
%  lambda_phi  d(lambda)/d(phi)
%  lambda_theta d(lambda)/d(theta)
%  br-b dot grad(r) (T)
%  bz-b dot grad(z) (T)
%  btor-toroidal field (T)
%  sqrt_g
%Comments:
%  A factor of rho**m is factored out of the Fourier coefficients prior
%    to spline fitting for increased accuracy. For rho less than
%    rho_thrift(1), all terms are assumed to vary as rho**m
%
% Adapted from THRIFT fortran source code and converted to MATLAB
% by JC Schmitt, 2009
%***********************************************************************

z_precision = 2.0e-9;

r=0.0;
z=0.0;
r_rho=0.0;
r_theta=0.0;
r_phi=0.0;
z_rho=0.0;
z_theta=0.0;
z_phi=0.0;
lambda_phi = 0.0;
lambda_theta = 0.0;

% loop over each mode
for ii = 1:mnmax
    angle_ = xm(ii) * theta - xn(ii) * phi;
    cosangle_ = cos(angle_);
    sinangle_ = sin(angle_);
    
    if rho == 0
        rho_m_norm = 1;
    else
        rho_m_norm = rho^abs(xm(ii));
    end
    
    r = r + rmnc_full(ii) * cosangle_ * rho_m_norm;
    z = z + zmns_full(ii) * sinangle_ * rho_m_norm;
    
    r_theta = r_theta - xm(ii) * rmnc_full(ii) * sinangle_ * rho_m_norm;
    r_phi = r_phi + xn(ii) * rmnc_full(ii) * sinangle_ * rho_m_norm;
    z_theta = z_theta + xm(ii) * zmns_full(ii) * cosangle_ * rho_m_norm;
    z_phi = z_phi - xn(ii) * zmns_full(ii) * cosangle_ * rho_m_norm;
    
    drmn_rho = drmnc_full(ii);
    dzmn_rho = dzmns_full(ii);
    
    %     r_rho = r_rho + drmnc_full(ii) * cosangle_;
    %     z_rho = z_rho + dzmns_full(ii) * sinangle_;
    rho_dumm = max([rho 10*z_precision]);
    %         rho_dumm = rho;
    drho_m = abs(xm(ii)) * rho_dumm^(abs(xm(ii))-1);
    r_rho = r_rho + drmn_rho * cosangle_ * rho_m_norm + ...
        rmnc_full(ii) * cosangle_* drho_m;
    z_rho = z_rho + dzmn_rho * sinangle_ * rho_m_norm + ...
        zmns_full(ii) * sinangle_ * drho_m;
    
    lambda_0 = lmns_full(ii);
    % unnormalize lambda_0
    lambda_0 = lambda_0 * rho_m_norm;
    lambda_theta = lambda_theta + lambda_0 * xm(ii) * cosangle_;
    lambda_phi = lambda_phi - lambda_0 * xn(ii) * cosangle_;
    
end

% Find jacobian and magnetic field components
tau = z_rho * r_theta - r_rho * z_theta;
sqrt_g = r * tau;
phiprm = 2 * rho * torflux_LCFS;
b_r = phiprm * ((iota - lambda_phi) * r_theta + ...
    (1 + lambda_theta) * r_phi) / (2*pi * sqrt_g);
b_z = phiprm * ((iota - lambda_phi) * z_theta + ...
    (1 + lambda_theta) * z_phi) / (2*pi * sqrt_g);
b_tor = r * phiprm * (1 + lambda_theta) / (2*pi * sqrt_g);

