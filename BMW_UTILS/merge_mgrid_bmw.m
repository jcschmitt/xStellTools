function merge_mgrid_bmw(filename_MGRID_Original, filename_MGRID_to_add, ...
    filename_MGRID_out)
% function merge_mgrid_bmw(filename_MGRID_Original, filename_MGRID_to_add, ...
%     filename_MGRID_out)
% This function makes a copy of the Original MGRID file
% 'filename_MGRID_Original' using the name given in 'filename_MGRID_out'.
% Next, it merges the mgrid data from 'filename_MGRID_to_add' into
% 'filename_MGRID_out', replacing the contents of the mgrid due to the
% last field coil in the original file.
%
% Generally, the original file should be generated with an extra 'dummy'
% coil that isn't actually in the experiment. The sole purpose of this
% dummy coiil is so that the mgrid output .nc file contains the proper
% number of dimensions, ordering, etc., to maintain compatibility with vmec
% and other related codes. It is easeier to 'swap' data in this file than
% it is to 'add' data to this .nc file.
%
% To Do: Add the option to create new files, rather than merging files.
%
% MGRID_Original = 'mgrid_w7x_with_BMW_space.nc';
% MGRID_to_add = 'bmw_CJV8I13_flux.nc';
% MGRID_out = 'merged_mgrid_out.nc';
%
% This should work with internal matlab netcdf support.

% make copy of original .nc file
disp(['<----Making a copy of ' filename_MGRID_Original ' with name ' ...
    filename_MGRID_out]);
copyfile(filename_MGRID_Original, filename_MGRID_out);

% open copy for writing
disp('<----Opening output file for writing');
ncid_out = netcdf.open(filename_MGRID_out, 'NC_WRITE');

disp('<----Opening bmw output to merge into output file')
% open mgrid to add for reading
ncid_merge = netcdf.open(filename_MGRID_to_add, 'NC_NOWRITE');

% Get grid sizes and confirm they are the same
disp('<----Checking grid sizes')
radInd_out = netcdf.inqDimID(ncid_out, 'rad');
[~, radDim_out] = netcdf.inqDim(ncid_out, radInd_out);
zeeInd_out = netcdf.inqDimID(ncid_out, 'zee');
[~, zeeDim_out] = netcdf.inqDim(ncid_out, zeeInd_out);
phiInd_out = netcdf.inqDimID(ncid_out, 'phi');
[~, phiDim_out] = netcdf.inqDim(ncid_out, phiInd_out);

radInd_merge = netcdf.inqDimID(ncid_merge, 'r');
[~, radDim_merge] = netcdf.inqDim(ncid_merge, radInd_merge);
zeeInd_merge = netcdf.inqDimID(ncid_merge, 'z');
[~, zeeDim_merge] = netcdf.inqDim(ncid_merge, zeeInd_merge);
phiInd_merge = netcdf.inqDimID(ncid_merge, 'phi');
[~, phiDim_merge] = netcdf.inqDim(ncid_merge, phiInd_merge);


if ~(radDim_out == radDim_merge)
    error('! Radial grid dimensions differ!')
end
if ~(zeeDim_out == zeeDim_merge)
    error('! Z grid dimensions differ!')
end
if ~(phiDim_out == phiDim_merge)
    error('! Phi grid dimensions differ!')
end

% Get the info for the original mgrid

ext_coils_Dim = netcdf.inqDimID(ncid_out, 'external_coils');
[~, ext_coils] = netcdf.inqDim(ncid_out, ext_coils_Dim);

disp(['<---- Number of external coils = ' num2str(ext_coils)]);
disp('<---- Will be replacing the LAST coil with the incoming Merge grid');

disp('<----Reading in original data Indexes (not the data itself)');
% Get the correct 'extension'
if ext_coils < 10
    field_extension = ['00' num2str(ext_coils)];
elseif ext_coils < 100
    field_extension = ['0' num2str(ext_coils)];
else
    error('<---- code not ready to handle more than 100 coils');
end

orig_br = ['br_' field_extension];
orig_bphi = ['bp_' field_extension];
orig_bz = ['bz_' field_extension];

brOrigIndex = netcdf.inqVarID(ncid_out, orig_br);
bphiOrigIndex = netcdf.inqVarID(ncid_out, orig_bphi);
bzOrigIndex = netcdf.inqVarID(ncid_out, orig_bz);

% br_orig_val = netcdf.getVar(ncid_out, brOrigIndex);
% bphi_orig_val = netcdf.getVar(ncid_out, bphiOrigIndex);
% bz_orig_val = netcdf.getVar(ncid_out, bzOrigIndex);

disp('<----Reading in new data to be merged');

brNewIndex = netcdf.inqVarID(ncid_merge, 'br_grid');
bphiNewIndex = netcdf.inqVarID(ncid_merge, 'bp_grid');
bzNewIndex = netcdf.inqVarID(ncid_merge, 'bz_grid');

br_merge_val = netcdf.getVar(ncid_merge, brNewIndex);
bphi_merge_val = netcdf.getVar(ncid_merge, bphiNewIndex);
bz_merge_val = netcdf.getVar(ncid_merge, bzNewIndex);

% write the data back to the file
disp('<----Writing merged data into output file');
% br
netcdf.putVar(ncid_out, brOrigIndex, br_merge_val);
% bp
netcdf.putVar(ncid_out, bphiOrigIndex, bphi_merge_val);
% bz
netcdf.putVar(ncid_out, bzOrigIndex, bz_merge_val);

disp('<----Closing files');

netcdf.close(ncid_out);
netcdf.close(ncid_merge);

disp('<----Finished. Check your work.');
