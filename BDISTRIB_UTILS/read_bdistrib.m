function boozer_data = read_Boozer_output(fname_ext)


try
    % matlab 2009a and newer
    disp('<---Using Matlab''s built-in netcdf support');
    ncid = netcdf.open(fname_ext, 'NC_NOWRITE');
    
    ns_bID = netcdf.inqVarID(ncid, 'ns_b');
    boozer_data.ns_b = netcdf.getVar(ncid, ns_bID); % number of surfs

    nfp_bID = netcdf.inqVarID(ncid, 'nfp_b');
    boozer_data.nfp_b = double(netcdf.getVar(ncid, nfp_bID)); % number of surfs

    lasym_ID = netcdf.inqVarID(ncid, 'lasym__logical__');
    boozer_data.lasym = netcdf.getVar(ncid, lasym_ID); % number of surfs

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
        error('<---where did this output file come from?')
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

    if (boozer_data.lasym)
        rmns_bID = netcdf.inqVarID(ncid, 'rmns_b');
        boozer_data.rmns_b = netcdf.getVar(ncid, rmns_bID)'; % cos terms of r-expansion
        
        bmns_bID = netcdf.inqVarID(ncid, 'bmns_b');
        boozer_data.bmns_b = netcdf.getVar(ncid, bmns_bID)'; % cos terms of r-expansion
        
        zmnc_bID = netcdf.inqVarID(ncid, 'zmnc_b');
        boozer_data.zmnc_b = netcdf.getVar(ncid, zmnc_bID)'; % sin terms of z-expansion
        
    end
    
    netcdf.close(ncid)
catch