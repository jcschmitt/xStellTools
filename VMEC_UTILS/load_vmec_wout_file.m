function [mnmax_vmec, xm_vmec, xn_vmec, spline_vmec, spline_d_vmec, ...
    phitot_vmec] = load_vmec_wout_file(file_vmec_output)

% get matlab version; version >= 7.6 has native netcdf support
matlabVersion = version;

% check for matlab version supporting mexnc/getnc
if (7.2 <= str2num(matlabVersion(1:3))) & (str2num(matlabVersion(1:3)) < 7.6)
    % load the variables
    if strfind(file_vmec_output, '.nc') % is it an .nc file?
        disp('Checking file');
        input_extension_vmec = getnc(file_vmec_output,'input_extension');
        disp('Opening .nc file');
        
        %USED
        disp('Loading data from VMEC wout file');
        ns_vmec = getnc(file_vmec_output,'ns')
        nfp_vmec = getnc(file_vmec_output, 'nfp')  % # of field periods
        mpol_vmec = getnc(file_vmec_output,'mpol')
        ntor_vmec = getnc(file_vmec_output,'ntor')
        mnmax_vmec = getnc(file_vmec_output,'mnmax')
        xn_vmec = getnc(file_vmec_output,'xn');      %(mn_mode)  toridal mode numbers
        xm_vmec = getnc(file_vmec_output,'xm');      %(mn_mode)  polidal mode numbers
        %         xm_vmec = -xm_vmec;  % this flips the VMEC Left handed coordinate system to a right handed system
        xm_vmec = xm_vmec;  %
        rmnc_vmec = getnc(file_vmec_output,'rmnc');        %(radius mn_mode ) cosmn component of cylindrical R, full mesh in [m]
        zmns_vmec = getnc(file_vmec_output,'zmns');        %(radius mn_mode ) sinmn component of cylindrical Z, full mesh in [m]
        lmns_vmec = getnc(file_vmec_output,'lmns');        %(radius mn_mode ) sinmn component of lambda, half mesh
        %         lmns_vmec = -lmns_vmec; % flipping the coordinate system from left- to right-handed
        %         iota_vmec = -getnc(file_vmec_output, 'iotaf'); % Flipped for RHS system
        iota_vmec = getnc(file_vmec_output, 'iotaf'); %
        phitot_vmec = getnc(file_vmec_output, 'phi');
    end; % gnetnc load variables
elseif (str2num(matlabVersion(1:3)) >= 7.6 || ...
        str2num(matlabVersion(1:4)) >= 7.10 )
    if strfind(file_vmec_output, '.nc') % is it an .nc file?
        disp('Checking file');
        ncID = netcdf.open(file_vmec_output,'nc_nowrite');
        input_extension_vmecID = netcdf.inqVarID(ncID,'input_extension');
        input_extension_vmec = netcdf.getVar(ncID,input_extension_vmecID);
        clear input_extension_vmecID;
        disp('Opening .nc file');
        
        %USED
        disp('Loading data from VMEC wout file');
        %        ns_vmec = getnc(file_vmec_output,'ns')
        ns_vmecID = netcdf.inqVarID(ncID,'ns');
        ns_vmec = double(netcdf.getVar(ncID,ns_vmecID));
        clear ns_vmecID;
        
        % # of field periods
        nfp_vmecID = netcdf.inqVarID(ncID,'nfp');
        nfp_vmec = double(netcdf.getVar(ncID,nfp_vmecID));
        clear nfp_vmecID;
        
        mpol_vmecID = netcdf.inqVarID(ncID,'mpol');
        mpol_vmec = double(netcdf.getVar(ncID,mpol_vmecID));
        clear mpol_vmecID;
        
        ntor_vmecID = netcdf.inqVarID(ncID,'ntor');
        ntor_vmec = double(netcdf.getVar(ncID, ntor_vmecID));
        clear ntor_vmec;
        
        mnmax_vmecID = netcdf.inqVarID(ncID, 'mnmax');
        mnmax_vmec = double(netcdf.getVar(ncID, mnmax_vmecID));
        clear mnmax_vmecID;
        
        %(mn_mode)  toridal mode numbers
        xn_vmecID = netcdf.inqVarID(ncID, 'xn');
        xn_vmec = double(netcdf.getVar(ncID, xn_vmecID));
        clear xn_vmecID;
        
        %(mn_mode)  polidal mode numbers
        xm_vmecID = netcdf.inqVarID(ncID, 'xm');
        xm_vmec = double(netcdf.getVar(ncID, xm_vmecID));
        clear xm_vmecID;
        
        %         xm_vmec = -xm_vmec;  % this flips the VMEC Left handed coordinate system to a right handed system
        
        %(radius mn_mode ) cosmn component of cylindrical R, full mesh in [m]
        rmnc_vmecID = netcdf.inqVarID(ncID, 'rmnc');
        rmnc_vmec = transpose(netcdf.getVar(ncID, rmnc_vmecID));
        clear rmnc_vmecID;
        
        %(radius mn_mode ) sinmn component of cylindrical Z, full mesh in [m]
        zmns_vmecID = netcdf.inqVarID(ncID, 'zmns');
        zmns_vmec = transpose(netcdf.getVar(ncID, zmns_vmecID));
        clear zmns_vmecID;
        
        %(radius mn_mode ) sinmn component of lambda, half mesh
        lmns_vmecID = netcdf.inqVarID(ncID, 'lmns');
        lmns_vmec = transpose(netcdf.getVar(ncID, lmns_vmecID));
        clear lmns_vmecID;
        %         lmns_vmec = -lmns_vmec; % flip the cood sys from left- to right-handed
        size(lmns_vmec)
        ns_vmec
        
        % NOT
        % Flipped for RHS system
        iota_vmecID = netcdf.inqVarID(ncID, 'iotaf');
        %         iota_vmec = -1.0*netcdf.getVar(ncID, iota_vmecID);
        iota_vmec = 1*netcdf.getVar(ncID, iota_vmecID);
        clear iota_vmecID;
        
        phitot_vmecID = netcdf.inqVarID(ncID, 'phi');
        phitot_vmec_radial = netcdf.getVar(ncID, phitot_vmecID);
        phitot_vmec = phitot_vmec_radial(end);
        clear phitot_vmecID;
        netcdf.close(ncID);
    end;
end;  % check matlab version
disp('Done loading data from VMEC wout file');



% set up radial array
s_vmec = ((1:ns_vmec)-1)/(ns_vmec-1);
rho_vmec = sqrt(s_vmec);

% set up a 'half mesh', between the nodes of the rho_vmec mesh, along with
% the endpoints 0 and 1.
rho_halfmesh_vmec = [0 sqrt(.5*(s_vmec(2:end)+s_vmec(1:(end-1)))) 1];

disp('Radial grids calculated');

disp('Normalizing lambda, rmn, and zmn terms');
% normalize lambda, r, and z coefficients to rho^m for better axial resolution
for ii = 1:mnmax_vmec
    %     disp(['Surface # ' num2str(ii)]);
    
    % calcuate normalizing factors
    rho_halfmesh_m_norm = rho_halfmesh_vmec.^abs(xm_vmec(ii));
    rho_m_norm = rho_vmec.^abs(xm_vmec(ii));
    
    % normalize the bulk of the lambda terms
    lmns_normalized(ii,2:ns_vmec) = lmns_vmec(2:ns_vmec,ii) ./ ...
        rho_halfmesh_m_norm(2:(end-1))';
    
    % extrapolate to the edges of the halfmesh radial grid
    % I use the 1.5/.5 weighting here  -JC Schmitt
    lmns_normalized(ii,1) = 1.5 * lmns_normalized(ii,2) - ...
        .5 * lmns_normalized(ii,3);
    lmns_normalized(ii,ns_vmec+1) = 1.5 * lmns_normalized(ii,ns_vmec) - ...
        .5 * lmns_normalized(ii,ns_vmec-1);
    
    % interpolate the lambda terms back to the full grid
    % lmns_on_rho_normalized(ii,:) = linterp(rho_halfmesh_vmec, lmns_normalized(ii,:), rho_vmec);
    lmns_on_rho_normalized(ii,:) = interp1(rho_halfmesh_vmec, lmns_normalized(ii,:), rho_vmec,'linear');
    
    % normalize the bulk of the rmnc and zmns terms
    rmnc_normalized(ii, 2:ns_vmec) = rmnc_vmec(2:ns_vmec, ii)' ./ ...
        rho_m_norm(2:end);
    zmns_normalized(ii, 2:ns_vmec) = zmns_vmec(2:ns_vmec, ii)' ./ ...
        rho_m_norm(2:end);
    % extrapolate to the axis for the normalized rmnc and zmns terms
    % Note by JCS: Here, the extrapolation is 'normal'. 2/1 weighting
    rmnc_normalized(ii, 1) = 2*rmnc_normalized(ii, 2) - rmnc_normalized(ii, 3);
    zmns_normalized(ii, 1) = 2*zmns_normalized(ii, 2) - zmns_normalized(ii, 3);
end
disp('Done normalizing.')

disp('Generating spline cooefficients');
% Going to make spline fits, but with 'double-sided' profiles
% get the spline coefficients for each of the lambda terms
% spline coefficients describe the radial variation of each harmonic.
% and get the spline coeffiecients for r,z, and radial deriviatives
% build the x-inputs
spline_vmec_x_val = [fliplr(-rho_vmec) rho_vmec];
% build the y-inputs
splines_vmec_y_vals(1,:) = [fliplr(iota_vmec') iota_vmec'];
ind_spline_start = 2;
ind_spline_end = ind_spline_start + mnmax_vmec - 1;
splines_vmec_y_vals((ind_spline_start:ind_spline_end),:) = ...
    [fliplr(lmns_on_rho_normalized) lmns_on_rho_normalized];
ind_spline_start = ind_spline_end + 1;
ind_spline_end = ind_spline_start + mnmax_vmec - 1;
splines_vmec_y_vals((ind_spline_start:ind_spline_end),:) = ...
    [fliplr(rmnc_normalized) rmnc_normalized];
ind_spline_start = ind_spline_end + 1;
ind_spline_end = ind_spline_start + mnmax_vmec - 1;
splines_vmec_y_vals((ind_spline_start:ind_spline_end),:) = ...
    [fliplr(zmns_normalized) zmns_normalized];

% get the spline coefficients
spline_vmec = csaps(spline_vmec_x_val, splines_vmec_y_vals);

% need the expressions for the derivates or r and z, too...
% build the y-inputs (note, I could use Matlab's built-in routines to build
% the spline for the derivates, but this is more transparent.
ind_spline_start = 1;
ind_spline_end = ind_spline_start + mnmax_vmec - 1;
splines_vmec_dy_vals((ind_spline_start:ind_spline_end),:) = ...
    [fliplr(rmnc_normalized) rmnc_normalized];
ind_spline_start = ind_spline_end + 1;
ind_spline_end = ind_spline_start + mnmax_vmec - 1;
splines_vmec_dy_vals((ind_spline_start:ind_spline_end),:) = ...
    [fliplr(zmns_normalized) zmns_normalized];
%
spline_d_vmec_1 = csaps(spline_vmec_x_val, splines_vmec_dy_vals);
spline_d_vmec = fnder(spline_d_vmec_1);
VMEC_DATA_LOADED = 1;
disp('Splines are generated. Data is persistent between function calls.')