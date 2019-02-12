function import_coils_from_mgrid_input_format(varargin)
% Import and (optionally) save and plot field coils.
% Import format is that which is compatible with mgrid/xgrid
% Saved format is compatible with the Biot-Savart and Matlab-based
% field-line following tools included in xStellTools.
%
% version_number: 1 = coils.extension (original)
% coil: structure defining the coil
% coil.num_turns : number of turns, or loops in the coil
% coil.turn_number(:).num_vertices : number of vertices for a turn. Does
%           not need to be equal for all turn_numbers
% coil.turn_number(:).x = x coordinates of the vertics, in meters
% coil.turn_number(:).y = y coordinates of the vertics, in meters
% coil.turn_number(:).z = z coordinates of the vertics, in meters
%
% Add example here:
% import_coils_from_mgrid_input_format('coils_file_in',
%   'coils.wistell_a_004', 'make_plots', 1, ...
%   'save_mat_data', 1, 'output_mat_filename', 'coilset_wistell_a_004.mat') 

% Parsing the input. Solution from
% https://stackoverflow.com/questions/2775263/how-to-deal-with-name-value-pairs-of-function-arguments-in-matlab
%# define defaults at the beginning of the code so that you do not need to
%# scroll way down in case you want to change something or if the help is
%# incomplete
options = struct('output_mat_filename', '','flf_file_prefix', '', ...
    'save_mat_data', 0, 'save_flf_data', 0, 'make_plots', 0, ...
    'coils_file_in', '', 'debug_level', 0, 'coil_order', 0, ...
    'winding_factors', 1);

%# read the acceptable names
option_names = fieldnames(options);

%# count arguments
num_args = length(varargin);
if round(num_args/2)~=num_args/2
    error('import_coils_from_mgrid_input_format needs propertyName/propertyValue pairs')
end

for pair = reshape(varargin,2,[]) %# pair is {propName;propValue}
    input_name = lower(pair{1}); %# make case insensitive
    
    if any(strcmp(input_name,option_names))
        %# overwrite options. If you want you can test for the right class here
        %# Also, if you find out that there is an option you keep getting wrong,
        %# you can use "if strcmp(inpName,'problemOption'),testMore,end"-statements
        options.(input_name) = pair{2};
    else
        error('%s is not a recognized parameter name',input_name)
    end
end

options

if isempty(options.coils_file_in)
    error('<----No coil file found. Check input struct');
end
if (options.save_mat_data && isempty(options.output_mat_filename))
    error('<---If saving the Matlab coil data, specify ''output_mat_filename''.');
end

fid = fopen(options.coils_file_in);
% Start reading file
% Skip first 4 lines
for ii = 1:3
    skipped = fgetl(fid);
    if options.debug_level >= 1000
        disp(skipped)
    end
end

coil_struct_temp_out = [];

% Start reading coords (and currents and labels)
num_coilturns = 0;
coilturn_started = 0;

% Read in line
A = fgetl(fid);
while (ischar(A) && ~strcmpi(A, 'end'))
    A_split = strsplit(A);
    len_A = length(A_split);
    if isempty(A_split{1})  % A leading space get through the strsplit
        ind_start = 2;
    else
        ind_start = 1;
    end
    
    
    % Determine if it is the last line of a coil, then append to correct coil
    % structure
    try
        if str2num(A_split{ind_start + 3}) == 0 % end of coil
            num_vertices = num_vertices + 1;
            x(num_vertices) = str2num(A_split{ind_start});
            y(num_vertices) = str2num(A_split{ind_start+1});
            z(num_vertices) = str2num(A_split{ind_start+2});
            coil_number = str2num(A_split{ind_start+4});
            coil_name = A_split{ind_start+5};
            
            coil_struct_temp_out{num_coilturns}.num_vertices = num_vertices;
            coil_struct_temp_out{num_coilturns}.x = x;
            coil_struct_temp_out{num_coilturns}.y = y;
            coil_struct_temp_out{num_coilturns}.z = z;
            coil_struct_temp_out{num_coilturns}.coil_number = coil_number;
            coil_struct_temp_out{num_coilturns}.coil_name = coil_name;
            
            coilturn_started = 0;
            % Else, if first, create new coil turn structure, add the coords from the line
            % to that structure
        elseif coilturn_started == 0
            coilturn_started = 1;
            num_coilturns = num_coilturns + 1;
            num_vertices = 1;
            x = str2num(A_split{ind_start});
            y = str2num(A_split{ind_start+1});
            z = str2num(A_split{ind_start+2});
            % Else, add the coords from the lineto that structure
        else
            num_vertices = num_vertices + 1;
            x(num_vertices) = str2num(A_split{ind_start});
            y(num_vertices) = str2num(A_split{ind_start+1});
            z(num_vertices) = str2num(A_split{ind_start+2});
            
        end
    catch
        keyboard
    end
    % read in next line
    A = fgetl(fid);
end


% Process coil_struct_temp_out to create coils structure

% coil: structure defining the coil
% coil.num_turns : number of turns, or loops in the coil
% coil.turn_number(:).num_vertices : number of vertices for a turn. Does
%           not need to be equal for all turn_numbers
% coil.turn_number(:).x = x coordinates of the vertics, in meters
% coil.turn_number(:).y = y coordinates of the vertics, in meters
% coil.turn_number(:).z = z coordinates of the vertics, in meters
% coil_name_list: a cell-list of the (unique) names of all of the coils
coil_name_list = {};
num_coils = 0;

for ii = 1:num_coilturns
    coil_name = coil_struct_temp_out{ii}.coil_name;
    % Check to see if coil with name 'coil_name' exists
    coil_exists = exist(coil_name);
    % If not, create it, and put the coilturn into the coil
    if ~coil_exists
        eval([coil_name, '.num_turns = 1;']);
        eval([coil_name, '.turn_number(1).x = coil_struct_temp_out{ii}.x;'])
        eval([coil_name, '.turn_number(1).y = coil_struct_temp_out{ii}.y;'])
        eval([coil_name, '.turn_number(1).z = coil_struct_temp_out{ii}.z;'])
        eval([coil_name, '.turn_number(1).num_vertices = coil_struct_temp_out{ii}.num_vertices;'])
        num_coils = num_coils + 1;
        coil_name_list{num_coils} = coil_name;
    else % Otherwise, add the coilturn to the existing coil structure
        eval(['next_index = ', coil_name, '.num_turns + 1'])
        eval([coil_name, '.num_turns = next_index']);
        eval([coil_name, '.turn_number(next_index).x = coil_struct_temp_out{ii}.x;'])
        eval([coil_name, '.turn_number(next_index).y = coil_struct_temp_out{ii}.y;'])
        eval([coil_name, '.turn_number(next_index).z = coil_struct_temp_out{ii}.z;'])
        eval([coil_name, '.turn_number(next_index).num_vertices = coil_struct_temp_out{ii}.num_vertices;'])
    end
end


if options.save_mat_data
    disp('<----Saving matlab data.');
    
    for ii = 1:num_coils
        % loop over each unique coil name and save it to the file
        this_coil_name = coil_name_list{ii};
        if ii == 1
            coil_order = options.coil_order;
            winding_factors = options.winding_factors;
            eval(['save(options.output_mat_filename, ''', this_coil_name, ''')']);
            eval(['save(options.output_mat_filename, ''coil_order'', ''-append'')']);
            eval(['save(options.output_mat_filename, ''winding_factors'', ''-append'')']);
        else
            eval(['save(options.output_mat_filename, ''', this_coil_name, ''', ''-append'')']);
        end
    end
end

if options.save_flf_data
    disp('<-----ENtering an untested section of the code!!!!')
    disp('<----Writing coils to individual files for FLF');
    for ii = 1:num_coilturns
        coil_name = coil_struct_temp_out{ii}.coil_name;
        filename = [options.flf_file_prefix 'index_' num2str(ii) '_' ...
            coil_name];
        filehandle = fopen(filename, 'w');
        num_vertices = coil_struct_temp_out{ii}.num_vertices;
        fprintf(filehandle, [num2str(num_vertices), '\n']);
        x = coil_struct_temp_out{ii}.x;
        y = coil_struct_temp_out{ii}.y;
        z = coil_struct_temp_out{ii}.z;
        for jj = 1:num_vertices
            fprintf(filehandle, '%.12e %.12e %.12e\n', x(jj), y(jj), z(jj));
        end
        fclose(filehandle);
    end
end

if options.make_plots
    
    % Step 6:  Pretty pictures
    figure;
    box on;hold on; axis equal; view(3)
    xlabel('x');ylabel('y');zlabel('z')
    
    %plot_fieldcoils(Mod, 'm');
    for ii = 1:num_coils
        % loop over each unique coil name and save it to the file
        this_coil_name = coil_name_list{ii};
        eval(['plot_fieldcoils(', this_coil_name, ', ''k'' )']);
    end
    
    
    
    % plot_coils(Coil_2, 'k');
    % plot_coils(Coil_3, 'r');
    % plot_coils(Coil_4, 'g');
    % plot_coils(Coil_5, [.2 .5 0]);
    % plot_coils(Coil_A, [.8 .5 0]);
    % plot_coils(Coil_B, 'y');
    % plot_coils(Sweep_Coil_1, 'r');
    % plot_coils(Sweep_Coil_2, 'g');
    % plot_coils(TRIM_A1, 'k');
    % plot_coils(TRIM_A2, 'b');
    % plot_coils(TRIM_A3, 'r');
    % plot_coils(TRIM_A4, 'c');
    % plot_coils(TRIM_B1, 'g');
    %
end
