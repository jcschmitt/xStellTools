function boozer_data = read_Boozer_output(fname_ext)


% load basic information.  see booz_xform (part of stellopt package) for details
try
    % matlab 2009a and newer
    disp('<---Using Matlab''s built-in netcdf support');
    ncid = netcdf.open(fname_ext, 'NC_NOWRITE');
    
    ns_bID = netcdf.inqVarID(ncid, 'ns_b');
    boozer_data.ns_b = netcdf.getVar(ncid, ns_bID); % number of surfs
    
    mboz_bID = netcdf.inqVarID(ncid, 'mboz_b');
    boozer_data.mboz_b = netcdf.getVar(ncid, 6); % something related to number of poloidal modes
    
    nboz_bID = netcdf.inqVarID(ncid, 'nboz_b');
    boozer_data.nboz_b = netcdf.getVar(ncid, nboz_bID); % " " of toroidal modes
    
    mnboz_bID = netcdf.inqVarID(ncid, 'mnboz_b');
    boozer_data.mnboz_b = netcdf.getVar(ncid, mnboz_bID); % " " of total modes = (nboz_b * 2 - 1) * (mboz_b) - nboz_b
    
    ixm_bID = netcdf.inqVarID(ncid, 'ixm_b');
    boozer_data.xm_b = double(netcdf.getVar(ncid, ixm_bID)); % poloidal mode numbers
    
    ixn_bID = netcdf.inqVarID(ncid, 'ixn_b');
    boozer_data.xn_b = double(netcdf.getVar(ncid, ixn_bID)); % toroidal mode numbers
        
    try
        bmnc_bID = netcdf.inqVarID(ncid, 'bmnc_b');
        boozer_data.bmnc_b = netcdf.getVar(ncid, bmnc_bID)'; % bmn mode magnitudes (signed)
    catch
        disp('<---Did not find bmnc_b.  Looking for bmn_b')
        boozer_data.bmnc_bID = netcdf.inqVarID(ncid, 'bmn_b');
        boozer_data.bmnc_b = netcdf.getVar(ncid, bmnc_bID)'; % bmn mode magnitudes (signed)
    end
    phi_bID = netcdf.inqVarID(ncid, 'phi_b');
    boozer_data.phi_b = netcdf.getVar(ncid, phi_bID)'; % enclosed toroidal flux
    
    gmn_bID = netcdf.inqVarID(ncid, 'gmn_b');
    boozer_data.gmn_b = netcdf.getVar(ncid, gmn_bID)'; % boozer g factor
    
    iota_bID = netcdf.inqVarID(ncid, 'iota_b');
    boozer_data.iota_b = netcdf.getVar(ncid, iota_bID)'; % boozer iota
    
    rmnc_bID = netcdf.inqVarID(ncid, 'rmnc_b');
    boozer_data.rmnc_b = netcdf.getVar(ncid, rmnc_bID)'; % cos terms of r-expansion
    
    zmns_bID = netcdf.inqVarID(ncid, 'zmns_b');
    boozer_data.zmns_b = netcdf.getVar(ncid, zmns_bID)'; % sin terms of z-expansion

    bvco_bID = netcdf.inqVarID(ncid, 'bvco_b');
    boozer_data.bvco_b = netcdf.getVar(ncid, bvco_bID)'; % The Boozer g factor
    
    buco_bID = netcdf.inqVarID(ncid, 'buco_b');
    boozer_data.buco_b = netcdf.getVar(ncid, buco_bID)'; % The Boozer I factor

    netcdf.close(ncid)
catch
    % matlab 2008b and older
    disp('<----Warning: Using 3rd party netcdf support.');
    disp('<----Warning: Some variables may need to have a transpose applied.');
    boozer_data.ns_b = getnc(fname_ext, 'ns_b'); % number of surfs
    boozer_data.mboz_b = getnc(fname_ext, 'mboz_b'); % something related to number of poloidal modes
    boozer_data.nboz_b = getnc(fname_ext, 'nboz_b'); % " " of toroidal modes
    boozer_data.mnboz_b = getnc(fname_ext, 'mnboz_b'); % " " of total modes = (nboz_b * 2 - 1) * (mboz_b) - nboz_b
    boozer_data.ixm_b = getnc(fname_ext, 'ixm_b'); % poloidal mode numbers
    boozer_data.ixn_b = getnc(fname_ext, 'ixn_b'); % toroidal mode numbers
    boozer_data.bmnc_b = getnc(fname_ext, 'bmnc_b'); % mode magnitudes (signed)
    boozer_data.phi_b = getnc(fname_ext, 'phi_b'); % enclosed toroidal flux
    boozer_data.gmn_b = getnc(fname_ext, 'gmn_b'); % boozer g factor
    boozer_data.iota_b = getnc(fname_ext, 'iota_b'); % boozer iota
    boozer_data.rmnc_b = getnc(fname_ext, 'rmnc_b'); % cos terms of r-expansion
    boozer_data.zmns_b = getnc(fname_ext, 'zmns_b'); % sin terms of z-expansion
    boozer_data.bvco_b = getnc(fname_ext, 'bvco_b'); % The Boozer g factor
    boozer_data.buco_b = getnc(fname_ext, 'buco_b'); % The Boozer I factor
end
