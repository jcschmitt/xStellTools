function [data_out, err_msg] = load_v3fit_result(path_to_file)


try
    data_out = [];
    
    %disp('<----Using Matlab''s built-in netcdf support');
    ncid = netcdf.open(path_to_file, 'NC_NOWRITE');
    
    
    nsignal_ID = netcdf.inqDimID(ncid, 'nsignal');
    [~, data_out.nsignal] = netcdf.inqDim(ncid, nsignal_ID); % number of signals
    nparam_ID = netcdf.inqDimID(ncid, 'nparam');
    [~, data_out.nparam] = netcdf.inqDim(ncid, nparam_ID); % number of signals
    
    var_list = {'derived_param_name', 'derived_param_index', 'derived_param_value', ...
        'derived_param_sigma', 'derived_param_corr', 'param_name', 'param_index', ...
        'param_value', 'param_sigma', 'param_corr', 'signal_eff_matrix', 'signal_name', ...
        'signal_type', 'signal_weight', 'signal_observed_value', 'signal_model_value', ...
        'signal_sigma', 'signal_model_value', 'signal_type', 'g2', 'nsteps'};
    
    for this_var = var_list
        the_expression = ['get_ncdata_helper(ncid, ''', this_var{1}, ''');'];
        this_data = eval(the_expression);
        next_expression = ['data_out.', this_var{1}, ' = this_data;'];
        try
            eval(next_expression)
        catch
            disp(['<----load_v3fit_result failed to load ' this_var{1} ' or evaluate: ' next_expression]);
        end
        
    end
    
    try
       data_out.signal_name_clean = data_out.signal_name';
       data_out.signal_name_clean2 = data_out.signal_name_clean(:,1:10);
       data_out.param_name_clean = data_out.param_name';
       
    catch
    end
    % if we got here, yea!
    %     C_S_1d = data_out.signal_sigma(:,end).^2;
    %     data_out.C_S_2d = zeros(data_out.nsignal, data_out.nsignal);
    %     for ii = 1:data_out.nsignal
    %         data_out.CS_2d(ii,ii) = C_S_1d(ii);
    %     end
    %     data_out.J_2d = zeros(data_out.nsignal, data_out.nparam);
    %     for ii = 1:data_out.nsignal
    %         delta_e_i = data_out.signal_model_value(1,ii, end) - ...
    %             data_out.signal_model_value(1,ii, end-1);
    %         for jj = 1:data_out.nparam
    %             delta_p_j = data_out.param_value(end) - data_out.param_value(end-1);
    %             data_out.J_2d(ii,jj) = delta_e_i / delta_p_j;
    %         end
    %     end
    %     data_out.J_2d.' * inv(data_out.C_S_2d) * data_out.J_2d
    %     data_out.J_2d.' * data_out.C_S_2d \ data_out.J_2d
    %
    %     data_out.K_2d = zeros(data_out.nsignal, data_out.nparam);
    %     for ii = 1:data_out.nsignal
    %         delta_e_i = data_out.signal_model_value(1,ii, end) - ...
    %             data_out.signal_model_value(1,ii, end-1);
    %         for jj = 1:data_out.nparam
    %             delta_p_j = data_out.param_value(end) - data_out.param_value(end-1);
    %             data_out.J_2d(ii,jj) = delta_e_i / delta_p_j;
    %         end
    %     end
    %
    %     C_P = ;
    %     J = ;
    %
    
catch
    disp(['<----Error reading ' this_var ' in load_v3fit_result with file: ' path_to_file]);
    err_msg = -1;
    keyboard
end

err_msg = 0;
try
    
    netcdf.close(ncid);
catch
    disp(['<----Error closing load_v3fit_result with file: ' path_to_file]);
    %data_out = [];
    err_msg = -1;
end


function data_out = get_ncdata_helper(ncid, id)

this_ID = netcdf.inqVarID(ncid, id);
data_out = netcdf.getVar(ncid, this_ID); % number of signals


