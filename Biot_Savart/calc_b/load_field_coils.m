function [coils, coilnames] = load_field_coils(field_coilset_id)
% function QHS46_coils = load_QHS46_coils()
%
% Loads the information for the details of the field coils into the
% coils structure.
%
% JC Schmitt, 2018

if nargin < 1
    error('k=====NO Field Coilset ID specified!. Stopping');
end

% this can be modified based on the local directory structre
in_foldername = '/p/stellopt/ANALYSIS/jschmitt/src/xStellTools/Biot_Savart/coilsets/';

in_filename = [field_coilset_id '.mat'];   

in_filename_complete = [in_foldername in_filename];

disp(['Loading field coils from: ' in_filename_complete]);

coils = load(in_filename_complete);

all_fields = fieldnames(coils);
coilnames = {};
num_coils = 0;
for ii = 1:length(all_fields)
    fieldname = all_fields{ii};
    if ~(strcmpi(fieldname, 'coil_order') || strcmpi(fieldname, 'winding_factors') )        
        disp(['<----Found coil: ' fieldname]);
        num_coils = num_coils + 1;
        coilnames{num_coils} = fieldname;
    end
end
