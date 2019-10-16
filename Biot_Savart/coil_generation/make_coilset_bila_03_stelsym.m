function make_coilset_bila_03_stelsym 


import_coils_from_mgrid_input_format('make_plots', 1, 'coils_file_in', ...
     'coils.bila_03', 'save_mat_data', 1, 'output_mat_filename', ...
     'coilset_bila_03_stelsym', 'winding_factors', [1], ...
     'coil_order', {'Mod'},'make_stellarator_symmetric', 1, ...
     'coil_format', 'FOCUS', 'number_field_periods', 5,...
     'coils_per_period', 12)
 
 
 coil_data = load('coilset_bila_03_stelsym.mat');
 
 figure; box on; hold on; axis equal
 plot_fieldcoils(coil_data.Mod, 'r');
 
 figure; box on; hold on; axis equal
 plot_coils(read_coils('coils.bila_03'))
 
 