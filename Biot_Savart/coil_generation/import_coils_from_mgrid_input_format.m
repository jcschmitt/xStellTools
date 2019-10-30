function import_coils_from_mgrid_input_format(varargin)
% Import and (optionally) save and plot field coils.
% Import format is that which is compatible with mgrid/xgrid
% Saved format is compatible with the Matlab-based Biot-Savart and
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
% import_coils_from_mgrid_input_format('make_plots', 1, 'coils_file_in',  ...
%   'coils.hsx_complete', 'save_mat_data', 1, 'output_mat_filename', ...
%   'coilset_hsx_complete', 'winding_factors', [1 1 1 1 1 1 1], ...
%   'coil_order', {'MainCoil', 'AuxCoil1', 'AuxCoil2', 'AuxCoil3',  ...
%   'AuxCoil4', 'AuxCoil5', 'AuxCoil6'})


% import_coils_from_mgrid_input_format('coils_file_in', ...
%   'coils.wistell_a_004', 'make_plots', 1, ...
%   'save_mat_data', 1, 'output_mat_filename', 'coilset_wistell_a_004.mat')
% import_coils_from_mgrid_input_format('make_plots', 1, 'coils_file_in',
%   'coils.w7x_asbuilt_16_v3fit_nfp_5', 'save_mat_data', 1,
%   'output_mat_filename', 'coilset_w7x_asbuilt16', 'winding_factors', [108
%   108 108 108 108 36 36 8 8 8 8 8 8 8 8 8 8 1 1 1 1 1 108 108 108 108 36
%   36], 'coil_order', {'AAE10_SC', 'AAE29_SC', 'AAE38_SC', 'AAE47_SC',
%   'AAE56_SC', 'AAE14_SC', 'AAE23_SC', 'CC1L', 'CC1U', 'CC5L', 'CC5U',
%   'CC4L', 'CC4U', 'CC3L', 'CC3U', 'CC2L', 'CC2U', 'AAQ11', 'AAQ22',
%   'AAQ31', 'AAQ41', 'AAQ51', 'AAE10', 'AAE29', 'AAE38', 'AAE47', 'AAE56',
%   'AAE14', 'AAE23'})
%
% import_coils_from_mgrid_input_format('make_plots', 1, 'coils_file_in',
%   'coils.bila_03', 'save_mat_data', 1, 'output_mat_filename',
%   'coilset_bila_03', 'winding_factors', [1], 'coil_order', {'MainCoil'})
%
% Parsing the input. Solution from
% https://stackoverflow.com/questions/2775263/how-to-deal-with-name-value-pairs-of-function-arguments-in-matlab
%# define defaults at the beginning of the code so that you do not need to
%# scroll way down in case you want to change something or if the help is
%# incomplete
options = struct('output_mat_filename', '', ...
    'flf_file_prefix', '', ...
    'save_mat_data', 0, ...
    'save_flf_data', 0,...
    'make_plots', 0, ...
    'coils_file_in', '', ...
    'debug_level', 0, ...
    'coil_order', 0, ...
    'winding_factors', 1, ...
    'make_stellarator_symmetric', 0, ...
    'coil_format', '', ...
    'number_field_periods', 0, ...
    'coils_per_period', 0);

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
end_found = 0;
while ~feof(fid);
    skipped = fgetl(fid);
    if options.debug_level >= 1000
        disp(skipped)
    end
    if (strcmpi(strtrim(skipped), '** coils_dot_starts_below **'))
        end_found = 1;
        break;
    end
end

if end_found
    % Skip first 4 lines
else
    % close it, repoen it, read three lines
    fclose(fid);
    fid = fopen(options.coils_file_in);
end
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

% Read in line:1
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
            % Else, add the coords from the line to that structure
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

% Check to see if we need to 'symmeterize' the coil (Coils aren't
% necessarily stellarator symmetric)
if (options.make_stellarator_symmetric)
    if (strcmpi(options.coil_format, 'focus'))
        disp('<----Making stellarator symmetric coils from FOCUS output');
        
        NFP = options.number_field_periods;
        CPP = options.coils_per_period;
        %keyboard
        % According to discusions with Tom K., field periods are already
        % identical. We just need to fix a single FP and rotate/flip
        % Construct an 'average' coil by flip/rotate. Structure starts out as:
        % coil_struct_temp_out{:}
        %     num_vertices: XYZ
        %                x: [1xXYZ double]
        %                y: [1xXYZ double]
        %                z: [1xXYZ double]
        %      coil_number: 1 ... NFP*CPP
        %        coil_name: 'Mod'
        if (mod(CPP,2) == 0) % even # of coils / period
            DEBUG_STELSYM = 1;
            CPHP = CPP / 2; % coils per half period
            if DEBUG_STELSYM
                figure(10001); hold on; box on; grid on;
                axis equal;
            end
            for ii = 1:CPHP
                index_coil_A = ii;
                index_coil_B = CPP + 1 - ii;
                % check that # of vertices match
                if (coil_struct_temp_out{index_coil_A}.num_vertices ~= ...
                        coil_struct_temp_out{index_coil_B}.num_vertices)
                    error('<----Vertices count mismatch in import-coils');
                end
                % assign number of vertices to stel sym structure
                coil_struct_stelsym_out{index_coil_A}.num_vertices = ...
                    coil_struct_temp_out{index_coil_A}.num_vertices;
                coil_struct_stelsym_out{index_coil_B}.num_vertices = ...
                    coil_struct_temp_out{index_coil_A}.num_vertices;
                
                % check that names match
                if ~(strcmpi(coil_struct_temp_out{index_coil_A}.coil_name, ...
                        coil_struct_temp_out{index_coil_B}.coil_name))
                    error('<----Coil name mismatch in import-coils');
                end
                % assign name to stel sym structure
                coil_struct_stelsym_out{index_coil_A}.coil_name = ...
                    coil_struct_temp_out{index_coil_A}.coil_name;
                coil_struct_stelsym_out{index_coil_B}.coil_name = ...
                    coil_struct_temp_out{index_coil_A}.coil_name;
                
                % check that numbers match
                if (coil_struct_temp_out{index_coil_A}.coil_number ~= ...
                        coil_struct_temp_out{index_coil_B}.coil_number)
                    error('<----Coil name mismatch in import-coils');
                end
                % assign numbers to stel sym structure
                coil_struct_stelsym_out{index_coil_A}.coil_number = ...
                    coil_struct_temp_out{index_coil_A}.coil_number;
                coil_struct_stelsym_out{index_coil_B}.coil_number = ...
                    coil_struct_temp_out{index_coil_A}.coil_number;
                
                % now, get an 'average x y z' by flipping and reversing
                % keyboard;
                x_p1 = coil_struct_temp_out{index_coil_A}.x;
                y_p1 = coil_struct_temp_out{index_coil_A}.y;
                z_p1 = coil_struct_temp_out{index_coil_A}.z;
                x_p2 = coil_struct_temp_out{index_coil_B}.x;
                y_p2 = coil_struct_temp_out{index_coil_B}.y;
                z_p2 = coil_struct_temp_out{index_coil_B}.z;
                % the _p2 portions need to be 'mirrored' (and their indices
                % are opposite, too)
                fpangle = (2*pi) / (NFP * 2);
                z_p2_fix = -z_p2(end:-1:1);
                phi_p2_fix = atan2(y_p2(end:-1:1), x_p2(end:-1:1));
                rho_p2_fix = sqrt(x_p2(end:-1:1).^2 + y_p2(end:-1:1).^2);
                %phi_p2_fix = fpangle - phi_p2;
                phi_p2_flip = fpangle - (phi_p2_fix - fpangle);
                x_p2_fix = rho_p2_fix .* cos(phi_p2_flip);
                y_p2_fix = rho_p2_fix .* sin(phi_p2_flip);
                x_new = 0 * x_p1;
                y_new = 0 * y_p1;
                z_new = 0 * z_p1;
                
                for jj = 1:coil_struct_stelsym_out{ii}.num_vertices
                    x_new(jj) = 0.5 * (x_p1(jj) + x_p2_fix(jj));
                    y_new(jj) = 0.5 * (y_p1(jj) + y_p2_fix(jj));
                    z_new(jj) = 0.5 * (z_p1(jj) + z_p2_fix(jj));
                    
                end
                
                coil_struct_stelsym_out{index_coil_A}.x = x_new;
                coil_struct_stelsym_out{index_coil_A}.y = y_new;
                coil_struct_stelsym_out{index_coil_A}.z = z_new;
                
                phi_new_fix = atan2(y_new(end:-1:1), x_new(end:-1:1));
                rho_new_fix = sqrt(x_new(end:-1:1).^2 + y_new(end:-1:1).^2);
                phi_new_flip = fpangle - (phi_new_fix - fpangle);
                x_new_fix = rho_new_fix .* cos(phi_new_flip);
                y_new_fix = rho_new_fix .* sin(phi_new_flip);
                
                coil_struct_stelsym_out{index_coil_B}.x = x_new_fix;
                coil_struct_stelsym_out{index_coil_B}.y = y_new_fix;
                coil_struct_stelsym_out{index_coil_B}.z = -z_new(end:-1:1);
                
                if DEBUG_STELSYM
                    figure(10001);
                    plot3(x_p1, y_p1, z_p1, 'ro-');
                    plot3(x_p2, y_p2, z_p2, 'g+-');
                    plot3(x_p2_fix, y_p2_fix, z_p2_fix, 'b+-');
                end
                
                if DEBUG_STELSYM
                    figure(10002);
                    subplot(3,2,1);box on; hold on;
                    plot(x_p1, '.-');
                    plot(x_p2_fix, 'o:');
                    plot(x_new, 's--');
                    legend('x p1', 'x p2 fix', 'x new');
                    
                    subplot(3,2,2);box on; hold on;
                    plot(x_p1 - x_p2_fix, '.-');
                    legend('x p1 - x p2 fix');
                    
                    subplot(3,2,3);box on; hold on;
                    plot(y_p1, '.-');
                    plot(y_p2_fix, 'o:');
                    plot(y_new, 's--');
                    legend('y p1', 'y p2 fix', 'y new');
                    
                    subplot(3,2,4);box on; hold on;
                    plot(y_p1 - y_p2_fix, '.-');
                    legend('y p1 - y p2 fix');
                    
                    subplot(3,2,5);box on; hold on;
                    plot(z_p1, '.-');
                    plot(z_p2_fix, 'o:');
                    plot(z_new, 's--');
                    legend('z p1', 'z p2 fix', 'z new');
                    
                    subplot(3,2,6);box on; hold on;
                    plot(z_p1 - z_p2_fix, '.-');
                    legend('z p1 - z p2 fix');
                    
                end
            end % for ii = 1:CPHP
            
            
            if DEBUG_STELSYM
                figure(10001);
                plot3(x_p1, y_p1, z_p1, 'ro-');
                plot3(x_p2, y_p2, z_p2, 'g+-');
                plot3(x_p2_fix, y_p2_fix, z_p2_fix, 'b+-');
            end
            
            if DEBUG_STELSYM
                figure(10002);
                subplot(3,2,1);box on; hold on;
                plot(x_p1, '.-');
                plot(x_p2_fix, 'o:');
                plot(x_new, 's--');
                legend('x p1', 'x p2 fix', 'x new');
                
                subplot(3,2,2);box on; hold on;
                plot(x_p1 - x_p2_fix, '.-');
                legend('x p1 - x p2 fix');
                
                subplot(3,2,3);box on; hold on;
                plot(y_p1, '.-');
                plot(y_p2_fix, 'o:');
                plot(y_new, 's--');
                legend('y p1', 'y p2 fix', 'y new');
                
                subplot(3,2,4);box on; hold on;
                plot(y_p1 - y_p2_fix, '.-');
                legend('y p1 - y p2 fix');
                
                subplot(3,2,5);box on; hold on;
                plot(z_p1, '.-');
                plot(z_p2_fix, 'o:');
                plot(z_new, 's--');
                legend('z p1', 'z p2 fix', 'z new');
                
                subplot(3,2,6);box on; hold on;
                plot(z_p1 - z_p2_fix, '.-');
                legend('z p1 - z p2 fix');
                
            end
            
            % finished with one field period now dow the rest via rotation
            for ii = 2:NFP
                % repeat the above steps, but skip all of the safety
                % checks.
                for jj = 1:CPP
                    index_coil_Afp = jj; % the index of the coil in the first period
                    index_coil_A = (ii-1)*CPP + jj; % the index of the coil to generate
                    
                    coil_struct_stelsym_out{index_coil_A}.num_vertices = ...
                        coil_struct_temp_out{index_coil_Afp}.num_vertices;
                    
                    coil_struct_stelsym_out{index_coil_A}.coil_name = ...
                        coil_struct_temp_out{index_coil_Afp}.coil_name;
                    
                    coil_struct_stelsym_out{index_coil_A}.coil_number = ...
                        coil_struct_temp_out{index_coil_Afp}.coil_number;
                    
                    offset_angle = (2*pi*(ii-1)) / (NFP);
                    
                    phi_new_rotate = atan2(...
                        coil_struct_temp_out{index_coil_Afp}.y, ...
                        coil_struct_temp_out{index_coil_Afp}.x);
                    rho_new_rotate = sqrt( ...
                        coil_struct_temp_out{index_coil_Afp}.x.^2 + ...
                        coil_struct_temp_out{index_coil_Afp}.y.^2);
                    phi_new_rotate = offset_angle + phi_new_rotate;
                    x_new_rotate = rho_new_rotate .* cos(phi_new_rotate);
                    y_new_rotate = rho_new_rotate .* sin(phi_new_rotate);
                    z_new_rotate = coil_struct_temp_out{index_coil_Afp}.z;
                    
                    coil_struct_stelsym_out{index_coil_A}.x = x_new_rotate;
                    coil_struct_stelsym_out{index_coil_A}.y = y_new_rotate;
                    coil_struct_stelsym_out{index_coil_A}.z = z_new_rotate;
                    
                   
                    figure(10003);hold on; box on; grid on;axis equal;
                    plot3(coil_struct_stelsym_out{index_coil_A}.x, ...
                        coil_struct_stelsym_out{index_coil_A}.y, ...
                        coil_struct_stelsym_out{index_coil_A}.z, 'o');
                    
                end
                
                
            end
        else
            error('<---An odd number of coils / period.  This section of code is not yet developed.');
        end
        
        
        
        
    end
    
end
% on exit:
%   coil: structure defining the coil
%   coil.num_turns : number of turns, or loops in the coil
%   coil.turn_number(:).num_vertices : number of vertices for a turn. Does
%             not need to be equal for all turn_numbers
%   coil.turn_number(:).x = x coordinates of the vertics, in meters
%   coil.turn_number(:).y = y coordinates of the vertics, in meters
%   coil.turn_number(:).z = z coordinates of the vertics, in meters
%   coil_name_list: a cell-list of the (unique) names of all of the coils
coil_name_list = {};
num_coils = 0;

if (options.make_stellarator_symmetric)
    disp('<----Using stellarator symmetric coils');
    coil_struct_temp_out = coil_struct_stelsym_out;
end

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
        eval(['next_index = ', coil_name, '.num_turns + 1;'])
        eval([coil_name, '.num_turns = next_index;']);
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
    colorlist = {'b', 'r', 'g', 'c', 'm', 'k', 'b', 'r', 'g', 'c', 'm', 'k', ...
        'b', 'r', 'g', 'c', 'm', 'k', 'b', 'r', 'g', 'c', 'm', 'k', ...
        'b', 'r', 'g', 'c', 'm'}
    % Step 6:  Pretty pictures
    figure;
    box on;hold on; axis equal; view(3)
    xlabel('x');ylabel('y');zlabel('z')
    
    %plot_fieldcoils(Mod, 'm');
    for ii = 1:num_coils
        % loop over each unique coil name and save it to the file
        this_coil_name = coil_name_list{ii};
        this_color = colorlist{ii};
        %eval(['plot_fieldcoils(', this_coil_name, ', ''k'' )']);
        eval(['plot_fieldcoils(', this_coil_name, ', ''', this_color, ''' )']);
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
