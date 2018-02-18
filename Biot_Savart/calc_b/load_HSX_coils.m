function HSX_coils = load_HSX_coils()
% function LTX_coils = load_HSX_coils()
%
% Loads the information for the details of the HSX field coils into the
% HSX_coils structure.
%

% The variable 'in_foldername' is a string that points to the location
% where this Biot-Savart code is located.  For instance, the following line
% works on a Macbook with the repository extracted in a directory in the
% MATLAB dicrectory.
% in_foldername = '/Users/user123/Documents/MATLAB/Biot_Savart/calc_b/';

% Initially, in_foldername will be empty, forcing you to fix it.
%in_foldername = '';
in_foldername = '/Users/schmittj/Documents/MATLAB/Biot_Savart/calc_b/';

in_filename_1 = 'HSX_main_coils.mat';
in_filename_2 = 'HSX_aux_coils.mat';

in_filename_complete_1 = [in_foldername in_filename_1];
in_filename_complete_2 = [in_foldername in_filename_2];

disp(['Loading HSX coils from: ' in_filename_complete_1 ' and ' in_filename_complete_2]);

% See comments below.
HSX_main_coils = load(in_filename_complete_1);
HSX_aux_coils = load(in_filename_complete_2);


% HSX_main_coils has the coordinates for the main field
% coils in coils_x, coils_y, and coils_z
% number of main coils: 48
% number of points in each main coil: 842
%
% HSX_aux_coils has a single 3x2096 structure, which needs to be split up
% into coil groups 1 - 6

% Now the hard part. Convert these files to the format that calc_b_HSX and
% calc_B_BiotSavart are excpecting

% This section builds the HSX coil set

HSX_coils.main.num_turns = 48;
for ii = 1:48
    HSX_coils.main.turn_number(ii).num_vertices = 842;
    ind_start = 842 * (ii - 1) + 1;
    ind_end = ind_start + 841;
    HSX_coils.main.turn_number(ii).x = HSX_main_coils.coils_x(ind_start:ind_end);
    HSX_coils.main.turn_number(ii).y = HSX_main_coils.coils_y(ind_start:ind_end);
    HSX_coils.main.turn_number(ii).z = HSX_main_coils.coils_z(ind_start:ind_end);
end

HSX_coils.aux1.num_turns = 8;
HSX_coils.aux2.num_turns = 8;
HSX_coils.aux3.num_turns = 8;
HSX_coils.aux4.num_turns = 8;
HSX_coils.aux5.num_turns = 8;
HSX_coils.aux6.num_turns = 8;

for jj = 1:8
    HSX_coils.aux1.turn_number(jj).num_vertices = 45;
    HSX_coils.aux2.turn_number(jj).num_vertices = 43;
    HSX_coils.aux3.turn_number(jj).num_vertices = 45;
    HSX_coils.aux4.turn_number(jj).num_vertices = 41;
    HSX_coils.aux5.turn_number(jj).num_vertices = 45;
    HSX_coils.aux6.turn_number(jj).num_vertices = 43;
    
    ind1_start = (jj - 1) * 45 + 1;
    ind1_end = ind1_start + 44;
    ind2_start = (jj - 1) * 43 + 361;
    ind2_end = ind2_start + 42;
    ind3_start = (jj - 1) * 45 + 705;
    ind3_end = ind3_start + 44;
    ind4_start = (jj - 1) * 41 + 1065;
    ind4_end = ind4_start + 40;
    ind5_start = (jj - 1) * 45 + 1393;
    ind5_end = ind5_start + 44;
    ind6_start = (jj - 1) * 43 + 1753;
    ind6_end = ind6_start + 42;
    
    HSX_coils.aux1.turn_number(jj).x = HSX_aux_coils.aux(1, ind1_start:ind1_end);
    HSX_coils.aux1.turn_number(jj).y = HSX_aux_coils.aux(2, ind1_start:ind1_end);
    HSX_coils.aux1.turn_number(jj).z = HSX_aux_coils.aux(3, ind1_start:ind1_end);
    
    HSX_coils.aux2.turn_number(jj).x = HSX_aux_coils.aux(1, ind2_start:ind2_end);
    HSX_coils.aux2.turn_number(jj).y = HSX_aux_coils.aux(2, ind2_start:ind2_end);
    HSX_coils.aux2.turn_number(jj).z = HSX_aux_coils.aux(3, ind2_start:ind2_end);
    
    HSX_coils.aux3.turn_number(jj).x = HSX_aux_coils.aux(1, ind3_start:ind3_end);
    HSX_coils.aux3.turn_number(jj).y = HSX_aux_coils.aux(2, ind3_start:ind3_end);
    HSX_coils.aux3.turn_number(jj).z = HSX_aux_coils.aux(3, ind3_start:ind3_end);
    
    HSX_coils.aux4.turn_number(jj).x = HSX_aux_coils.aux(1, ind4_start:ind4_end);
    HSX_coils.aux4.turn_number(jj).y = HSX_aux_coils.aux(2, ind4_start:ind4_end);
    HSX_coils.aux4.turn_number(jj).z = HSX_aux_coils.aux(3, ind4_start:ind4_end);
    
    HSX_coils.aux5.turn_number(jj).x = HSX_aux_coils.aux(1, ind5_start:ind5_end);
    HSX_coils.aux5.turn_number(jj).y = HSX_aux_coils.aux(2, ind5_start:ind5_end);
    HSX_coils.aux5.turn_number(jj).z = HSX_aux_coils.aux(3, ind5_start:ind5_end);
    
    HSX_coils.aux6.turn_number(jj).x = HSX_aux_coils.aux(1, ind6_start:ind6_end);
    HSX_coils.aux6.turn_number(jj).y = HSX_aux_coils.aux(2, ind6_start:ind6_end);
    HSX_coils.aux6.turn_number(jj).z = HSX_aux_coils.aux(3, ind6_start:ind6_end);
    
end



