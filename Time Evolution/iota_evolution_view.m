function [] = ...
    iota_evolution_view(evolution_parameter_file_in)
% function [] = ...
%     iota_evolution(evolution_parameter_file_in, SHOW_PLOTS, RUN_EVOLUTION)
%
% evolution_parameter_file_in = filename of specification

% FIX THIS PATH FOR YOUR FILE STRUCTURE
%output_path = 'D:\JohnS\Plots and Notes\Time Evolution\simulation_data\output'
%output_file = [output_path '\' evolution_parameter_file_in];

output_path = '/Users/schmittj/Documents/MATLAB/TIme Evolution/simulation_data/output';
%output_path = '/Users/jschmitt/Documents/HSX Material/Work PC/D drive/JohnS/Plots and Notes/Time Evolution/simulation_data/output'
output_file = [output_path '/' evolution_parameter_file_in];

% options
MAKE_3D_PLOTS = 1;
MAKE_STANDARD_PLOTS = 1;
PLOT_ALL_IOTA_TIMES = 0;
SURF_PLOT = 0;
PLOT_R = 1; % if not r, then plot in s!
% FIT_TYPE = '';  % default, no fit
%FIT_TYPE = 'sum_atan';  %
 FIT_TYPE = '';  %
% FIT_TYPE = '';  %

t_end = 1.001; % [sec]  % end of simulation time
delta_t = .001; % msec timestep in simulation time

mu_0 = pi * 4e-7; % permeability of free space

%TIMES_TO_PLOT_IOTA = [0.01 0.05 .1 .2 .4 1 2 4 8];%.01 .05 .10 .2 .4 .6 .8 1];
TIMES_TO_PLOT_IOTA = [1 2.5 4 7 9 ];%.01 .05 .10 .2 .4 .6 .8 1];
%TIMES_TO_PLOT_IOTA = [1 2.5 5 10 25 50 90 ];%.01 .05 .10 .2 .4 .6 .8 1];
%TIMES_TO_PLOT_IOTA = [1:4:29 ];%.01 .05 .10 .2 .4 .6 .8 1];
%TIMES_TO_PLOT_IOTA = [.01:.01:.07];%.01 .05 .10 .2 .4 .6 .8 1];
% TIMES_TO_PLOT_IOTA = [.002:.002:.05];
% TIMES_TO_PLOT_IOTA = [.2 2 25  ];
% TIMES_TO_PLOT_IOTA = [.01 .015 0.02  ];
% TIMES_TO_PLOT_IOTA = [.01 .025 .05  ];
% TIMES_TO_PLOT_IOTA = [.025 .047  ];
%    TIMES_TO_PLOT_IOTA = [ 0.001 0.002 0.003 0.004 .005 0.006 0.007  ];
% TIMES_TO_PLOT_IOTA = [.0 .005  .025 .125 1  ];
%TIMES_TO_PLOT_IOTA = [0.015  0.030 0.05 0.075 ...
%     0.1 0.2  ];
% TIMES_TO_PLOT_IOTA = [.001 .005 0.01 .025 0.033 0.05 0.066 0.075 ...
%     .0875 0.095 0.1:0.01:0.250  ];

COLORS_OF_TIMES = 'brmcygkbmrcygkbmrcygkbmrcygkbmrcygkbmrcygkbmrcyg';

PROFILE_LW = 2;
PROFILE_FS = 14;
MarkerSize1 = 6;
Linetype1 = '-';
Symboltype1 = '';


disp(['Loading data from ' output_file]);
load(output_file)

if MAKE_STANDARD_PLOTS
    disp('Making plots')
    
    for ii = 1:length(TIMES_TO_PLOT_IOTA)
        % INDICES_OF_TIMES_TO_PLOT_IOTA(ii) = find( abs(tspan - TIMES_TO_PLOT_IOTA(ii)) < 1e-5);
        INDICES_OF_TIMES_TO_PLOT_IOTA(ii) = ...
            round( (TIMES_TO_PLOT_IOTA(ii) - tspan(1)) / ...
            (tspan(2) - tspan(1)));
    end
    
    if SURF_PLOT
        figure
        subplot(1,2,1); box on;
        surf(s_local,tspan,iota, 'LineStyle', 'None')
        subplot(1,2,2); box on;
        surf(s_local,tspan,delta_iota_thing, 'LineStyle', 'None')
        xlabel('s');
        ylabel('time');
        %         alpha(.6)
    end
    
    if PLOT_R
        xlabels = '\rho';
    else
        xlabels = 's';
    end
    
    if PLOT_ALL_IOTA_TIMES
        
        figure
        %     subplot(2,1,1)
        hold on;box on
        for ii = 1:ceil(t_end/delta_t)
            if PLOT_R
                plot(rho_local(1:end), iota(ii, (1:end)), 'b:', 'LineWidth', 2);
            else
                plot(s_local(1:end), iota(ii, (1:end)), 'b:', 'LineWidth', 2);
            end
        end
        xlabel(xlabels);
        ylabel('\iota');
        %     axis([0 1 .8 1.2]);
    end
    
    figure;hold on;box on
    plot(tspan, I_tor_net, 'LineWidth', PROFILE_LW);
    plot(tspan, I_tor(:, end), 'LineWidth', PROFILE_LW);
    xlabel('sec', 'FontSize', PROFILE_FS);
    title('I_{tor, net}', 'FontSize', PROFILE_FS);
    title('I_{tor, simul}', 'FontSize', PROFILE_FS);
    ylabel('A', 'FontSize', PROFILE_FS);
    %     cur_axis = axis;
    %     axis([0 TIMES_TO_PLOT_IOTA(end) cur_axis(3) cur_axis(4)]);
    %xlim([-0.005 0.075])
    
    figure(7003);hold on;box on
    plot(tspan, V_loop, 'LineWidth', PROFILE_LW);
    xlabel('sec', 'FontSize', PROFILE_FS);
    title('V_{loop}', 'FontSize', PROFILE_FS);
    ylabel('A', 'FontSize', PROFILE_FS);
    %xlim([-0.005 0.075])
    
    
    figure;hold on;box on
    for ii = 1:length(INDICES_OF_TIMES_TO_PLOT_IOTA)
        index_to_plot = INDICES_OF_TIMES_TO_PLOT_IOTA(ii);
        if PLOT_R
            plot(rho_local(2:end), iota(index_to_plot, ([ 2:end])), [Linetype1 Symboltype1 COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW, 'MarkerSize', MarkerSize1);
            %             plot(mid_rho_local, iota_halfmesh(index_to_plot, :), [Linetype1 'o' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW, 'MarkerSize', MarkerSize1);
        else
            plot(s_local(2:end), iota(index_to_plot, (2:end)), ['.' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW, 'MarkerSize', MarkerSize1); %#ok<*COLND>
            %             plot(mid_s_local, iota_halfmesh(index_to_plot, :), [Linetype1 'o' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW, 'MarkerSize', MarkerSize1);
        end
    end
    %         xlabel('\rho', 'FontSize', PROFILE_FS);
    xlabel(xlabels, 'FontSize', PROFILE_FS);
    ylabel('\iota', 'FontSize', PROFILE_FS);
    cur_axis = axis;
    xlim([0 1])
    %     axis([0 1 cur_axis(3) cur_axis(4)]);
    legend(num2str(TIMES_TO_PLOT_IOTA'),'Location','best');

    figure;hold on;box on
    for ii = 1:length(INDICES_OF_TIMES_TO_PLOT_IOTA)
        index_to_plot = INDICES_OF_TIMES_TO_PLOT_IOTA(ii);
        if PLOT_R
            plot(rho_local(2:end), iota(index_to_plot, 2:end)-...
                iota(1, 2:end), [Linetype1 Symboltype1 COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW, 'MarkerSize', MarkerSize1);
            %             plot(mid_rho_local, iota_halfmesh(index_to_plot, :), [Linetype1 'o' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW, 'MarkerSize', MarkerSize1);
        else
            plot(s_local(2:end), iota(index_to_plot, 2:end)-...
                iota(index_to_plot, 2:end), ['.' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW, 'MarkerSize', MarkerSize1);
            %             plot(mid_s_local, iota_halfmesh(index_to_plot, :), [Linetype1 'o' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW, 'MarkerSize', MarkerSize1);
        end
    end
    %         xlabel('\rho', 'FontSize', PROFILE_FS);
    xlabel(xlabels, 'FontSize', PROFILE_FS);
    ylabel('\delta\iota_{t=0}', 'FontSize', PROFILE_FS);
    cur_axis = axis;
    xlim([0 1])
    %     axis([0 1 cur_axis(3) cur_axis(4)]);
    legend(num2str(TIMES_TO_PLOT_IOTA'),'Location','best');

    
    %     figure;hold on;box on
    %     for ii = 1:length(INDICES_OF_TIMES_TO_PLOT_IOTA)
    %         index_to_plot = INDICES_OF_TIMES_TO_PLOT_IOTA(ii);
    %         if PLOT_R
    %             plot(rho_local(2:end), -1./(4+1*iota(index_to_plot, (2:end))), [Linetype1 Symboltype1 COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW, 'MarkerSize', MarkerSize1);
    %         else
    %             plot(s_local(2:end), -1./(4+1*iota(index_to_plot, (2:end))), [Linetype1 Symboltype1 COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW, 'MarkerSize', MarkerSize1);
    %         end
    %     end
    %     %     xlabel('\rho', 'FontSize', PROFILE_FS);
    %     xlabel(xlabels, 'FontSize', PROFILE_FS);
    %     ylabel('-q', 'FontSize', PROFILE_FS);
    %     legend(num2str(TIMES_TO_PLOT_IOTA'),'Location','best');
    
    figure;hold on; box on;
    %     subplot(3,1,2);hold on;box on
    for ii = 1:length(INDICES_OF_TIMES_TO_PLOT_IOTA)
        index_to_plot = INDICES_OF_TIMES_TO_PLOT_IOTA(ii);
        
        %         make fit
        fit_options = optimset('MaxFunEvals', 1e4, 'TolFun', 1e-8, 'TolX', 1e-6, ...
            'MaxIter', 1e3);
        
        if 1
            pcurr_type = 'sum_atan, 3 terms, ac(0)=0, ac(1) = 1, ';
            pcurr_guess = [1.  1.5 1. ];
            pcurr_UB = [ 10 1.51 1.0002 ];
            %                 pcurr_UB = [ 1 1.0001   1.0001  ];
            %         pcurr_LB = [ 14.999  eps 0.9999 ];
            %                 pcurr_LB = [ 1  eps 0.9999 ];
            pcurr_LB = [ .25  1.49 .999 ];
            %                 pcurr_type = 'sum_atan, 7 terms, ac(0)=0, ac(1) = 1';
            %                 pcurr_guess = [ 8  1.0535  1 -.5  1 .1 5];
            %                 pcurr_guess = [ 10  0.91387  3.4016 -0.066685  9.9992 0.0019433  6.6819];
            %                 pcurr_UB = [ 10 10   10  10 10 4   10];
            %                 pcurr_LB = [ 0  eps eps -10 0  eps eps];
            try
                pcurr_fit_params{ii} = lsqcurvefit(@pcurr_sa01_fit, pcurr_guess, ...
                    s_local(1:end), ...
                    I_tor(index_to_plot, (1:end))/I_tor(index_to_plot, end), ...
                    pcurr_LB, pcurr_UB, fit_options);
            catch
                pcurr_fit_params{ii} = [0 1 1];
            end
            pcurr_fitline = I_tor(index_to_plot, end) * ...
                pcurr_sa01_fit(pcurr_fit_params{ii}, s_local(1:end));
            dpcurr_fitline{ii} = diff(pcurr_fitline) ./ diff(s_local(1:end));
        elseif 0
        pcurr_type = 'sum_atan, 4 terms';
        pcurr_guess = [-0 1  5 1 1 ];
        pcurr_UB = [   1   4.0001   10.0001 5.51   1.0001  ];
        pcurr_LB = [   -1  0       .25       .5  0.9999 ];
        pcurr_fit_params{ii} = lsqcurvefit(@pcurr_sa_fit, pcurr_guess, ...
            s_local(1:end), ...
            I_tor(index_to_plot, (1:end))/I_tor(index_to_plot, end), ...
            pcurr_LB, pcurr_UB, fit_options);
        pcurr_fitline = I_tor(index_to_plot, end) * ...
            pcurr_sa_fit(pcurr_fit_params{ii}, s_local(1:end));
        dpcurr_fitline{ii} = diff(pcurr_fitline) ./ diff(s_local(1:end));
        end
        
        %         pcurr_type = 'power series, ac(0..4)';
        %         pcurr_guess = [1 1  1 1 ];
        %         pcurr_UB = [ 100  100  100  100];
        %         pcurr_LB = [-100 -100 -100 -100];
        %         pcurr_fit_params = lsqcurvefit(@pcurr_ps_fit, pcurr_guess, ...
        %             s_local(1:end), ...
        %             I_tor(index_to_plot, (1:end))/I_tor(index_to_plot, end), ...
        %             pcurr_LB, pcurr_UB, fit_options);
        %         pcurr_fitline = I_tor(index_to_plot, end) * ...
        %             pcurr_ps_fit(pcurr_fit_params, s_local(1:end));
        
        
        %         disp(['V3FIT Initial guess @ t= ' num2str(TIMES_TO_PLOT_IOTA(ii))]);
        %         disp(pcurr_type);
        %         disp(['cutor:  based on I_tor(time_index, end): ' ...
        %             num2str(I_tor(index_to_plot, end))]);
        % %         disp(['ac:  ' num2str([0 pcurr_fit_params])]);
        %         disp(['ac:  ' num2str([0 1 pcurr_fit_params])]);
        
        if PLOT_R
            plot(rho_local(1:end), I_tor(index_to_plot, (1:end)), [Linetype1 Symboltype1 COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
            plot(rho_local(1:end), pcurr_fitline, ['' '.' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
%             plot(mid_rho_local, I_tor_halfmesh(index_to_plot, :), [Linetype1 'o' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
        else
            plot(s_local(2:end), I_tor(index_to_plot, (2:end)), [Linetype1 Symboltype1 COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
            plot(s_local(2:end), pcurr_fitline(2:end), ['' '.' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
            %             plot(mid_s_local, I_tor_halfmesh(index_to_plot, :), [Linetype1 'o' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
        end
    end
    %     xlabel('\rho', 'FontSize', PROFILE_FS);
    xlabel(xlabels, 'FontSize', PROFILE_FS);
    title('I_{tor}', 'FontSize', PROFILE_FS);
    ylabel('A', 'FontSize', PROFILE_FS);
    legend(num2str(TIMES_TO_PLOT_IOTA'),'Location','best');
    %     axis([0 1 .8 1.2]);

    for ii = 1:length(INDICES_OF_TIMES_TO_PLOT_IOTA)
        index_to_plot = INDICES_OF_TIMES_TO_PLOT_IOTA(ii);
        disp(['V3FIT Initial guess @ t= ' num2str(TIMES_TO_PLOT_IOTA(ii))]);
        disp(pcurr_type);
        disp(['cutor:  based on I_tor(time_index, end): ' ...
            num2str(I_tor(index_to_plot, end))]);
%         disp(['ac = ' num2str([0 1 pcurr_fit_params{ii}])]);
        disp(['ac = ' num2str([pcurr_fit_params{ii}])]);
    end

    figure;hold on; box on;
    %     subplot(3,1,2);hold on;box on
    for ii = 1:length(INDICES_OF_TIMES_TO_PLOT_IOTA)
        index_to_plot = INDICES_OF_TIMES_TO_PLOT_IOTA(ii);
        if PLOT_R
            plot(rho_local(1:end), F_pol(index_to_plot, (1:end)), [Linetype1 Symboltype1 COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
            plot(mid_rho_local, F_pol_halfmesh(index_to_plot, :), [Linetype1 'o' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
        else
            plot(s_local(1:end), F_pol(index_to_plot, (1:end)), [Linetype1 Symboltype1 COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
            plot(mid_s_local, F_pol_halfmesh(index_to_plot, :), [Linetype1 'o' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
        end
    end
    %     xlabel('\rho', 'FontSize', PROFILE_FS);
    xlabel(xlabels, 'FontSize', PROFILE_FS);
    title('F_{pol}', 'FontSize', PROFILE_FS);
    ylabel('A', 'FontSize', PROFILE_FS);
    legend(num2str(TIMES_TO_PLOT_IOTA'),'Location','best');
    %     axis([0 1 .8 1.2]);
    xlim([0 1])
    
    figure;hold on; box on;
    %     subplot(3,1,2);hold on;box on
    for ii = 1:length(INDICES_OF_TIMES_TO_PLOT_IOTA)
        index_to_plot = INDICES_OF_TIMES_TO_PLOT_IOTA(ii);
        if PLOT_R
            plot(rho_local(1:end), F_pol(index_to_plot, (1:end))-F_pol(1, (1:end)), [Linetype1 Symboltype1 COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
            plot(mid_rho_local(1:end), F_pol_halfmesh(index_to_plot, (1:end))-F_pol_halfmesh(1, (1:end)), [Linetype1 'o' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
        else
            plot(s_local(1:end), F_pol(index_to_plot, (1:end))-F_pol(1, (1:end)), [Linetype1 Symboltype1 COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
            plot(mid_s_local(1:end), F_pol_halfmesh(index_to_plot, (1:end))-F_pol_halfmesh(1, (1:end)), [Linetype1 'o' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
        end
    end
    %     xlabel('\rho', 'FontSize', PROFILE_FS);
    xlabel(xlabels, 'FontSize', PROFILE_FS);
    title('F_{pol} - F_{pol vacuum}', 'FontSize', PROFILE_FS);
    ylabel('A', 'FontSize', PROFILE_FS);
    legend(num2str(TIMES_TO_PLOT_IOTA'),'Location','best');
    %     axis([0 1 .8 1.2]);
    xlim([0 1])
    
    figure;hold on; box on;
    %     subplot(3,1,3);hold on;box on
    for ii = 1:length(INDICES_OF_TIMES_TO_PLOT_IOTA)
        index_to_plot = INDICES_OF_TIMES_TO_PLOT_IOTA(ii);
        if PLOT_R
            %             plot(rho_local(1:end), dI_tor(index_to_plot, (1:end))/Phi_LCFS_in, [Linetype1 Symboltype1 COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
            plot(mid_rho_local(1:end), sqrt(B2_halfmesh) .* dI_tor_halfmesh(index_to_plot, (1:end))/Phi_LCFS_in/1e3, [ '-' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
%            plot(mid_rho_local(2:end), sqrt(B2_halfmesh(2:end)) .* dpcurr_fitline{ii}((2:end))/Phi_LCFS_in/1e3, [ '--' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
        else
            %             plot(s_local(1:end), dI_tor(index_to_plot, (1:end))/Phi_LCFS_in, [Linetype1 Symboltype1 COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
%             disp('applying quick fix of -1  (JC SCHMITT)') % see the next line
%             plot(mid_s_local(2:end), (1) * sqrt(B2_halfmesh(2:end)) .*  dI_tor_halfmesh(index_to_plot, (2:end))/Phi_LCFS_in/1e3, [Linetype1 'o' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
            plot(mid_s_local(1:end), (1) * sqrt(B2_halfmesh(1:end)) .*  dI_tor_halfmesh(index_to_plot, (1:end))/Phi_LCFS_in/1e3, [Linetype1 'o' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
        end
    end
    %     xlabel('\rho', 'FontSize', PROFILE_FS);
    xlabel(xlabels, 'FontSize', PROFILE_FS);
    title('J_{tor}', 'FontSize', PROFILE_FS);
    ylabel('kA / m^2', 'FontSize', PROFILE_FS);
    legend(num2str(TIMES_TO_PLOT_IOTA'),'Location','best');
    xlim([0 1])

        figure;hold on; box on;
    %     subplot(3,1,3);hold on;box on
    my_r_mat = zeros(length(tspan), length(B2_halfmesh));
    my_t_mat = zeros(length(tspan), length(B2_halfmesh));

    my_jdotb_mat = zeros(length(tspan), length(B2_halfmesh));
        
    for ii = 1:length(tspan)
        my_t_mat(ii,:) = tspan(ii);
        my_jdotb_mat(ii,:) = (-sqrt(B2_halfmesh) .* dI_tor_halfmesh(ii, :) / (1e3*Phi_LCFS_in));            
        if PLOT_R
            my_r_mat(ii,:) = (mid_rho_local);
        else
            my_r_mat(ii,:) = (mid_s_local);
        end
    end
    surf(my_t_mat(100:end,:), -1*my_r_mat(100:end,:), my_jdotb_mat(100:end,:), 'EdgeColor', 'Interp');
    
    %     xlabel('\rho', 'FontSize', PROFILE_FS);
    ylabel(xlabels, 'FontSize', PROFILE_FS);
    xlabel('time', 'FontSize', PROFILE_FS);
    title('J_{tor}', 'FontSize', PROFILE_FS);
    zlabel('kA / m^2', 'FontSize', PROFILE_FS);
    legend(num2str(TIMES_TO_PLOT_IOTA'),'Location','best');
    %xlim([0 1])

    
    figure;hold on; box on;
    %     subplot(3,1,3);hold on;box on
    for ii = 1:length(INDICES_OF_TIMES_TO_PLOT_IOTA)
        index_to_plot = INDICES_OF_TIMES_TO_PLOT_IOTA(ii);
        if PLOT_R
            %             plot(rho_local(1:end), J_inductive(index_to_plot, (1:end)), [Linetype1 Symboltype1 COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
            plot(mid_rho_local(1:end), 1* sqrt(B2_halfmesh(1:end)) .* J_inductive_halfmesh(index_to_plot, (1:end))/1e3, [Linetype1 '' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
        else
            %             plot(s_local(1:end), J_inductive(index_to_plot, (1:end)), [Linetype1 Symboltype1 COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
%                 disp('applying quick fix of -1  (JC SCHMITT)') % see the next line
            plot(mid_s_local(1:end), 1* sqrt(B2_halfmesh(1:end)) .* J_inductive_halfmesh(index_to_plot, (1:end))/1e3, [Linetype1 'o' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
        end
    end
    %     xlabel('\rho', 'FontSize', PROFILE_FS);
    xlabel(xlabels, 'FontSize', PROFILE_FS);
    title('J_{ind}', 'FontSize', PROFILE_FS);
    ylabel('kA / m^2', 'FontSize', PROFILE_FS);
    legend(num2str(TIMES_TO_PLOT_IOTA'),'Location','best');
    xlim([0 1])

    figure;hold on; box on;
    %     subplot(3,1,3);hold on;box on
    for ii = 1:length(INDICES_OF_TIMES_TO_PLOT_IOTA)
        index_to_plot = INDICES_OF_TIMES_TO_PLOT_IOTA(ii);
        if PLOT_R
            plot(mid_rho_local(1:end), ...
                J_inductive_halfmesh(index_to_plot, (1:end)) / ...
                J_inductive_halfmesh(end, (1)), ...
                [Linetype1 '' COLORS_OF_TIMES(ii)], ...
                'LineWidth', PROFILE_LW);
        else
            plot(mid_s_local(1:end), ...
                J_inductive_halfmesh(index_to_plot, (1:end)) / ...
                J_inductive_halfmesh(end, (1)), ...
                [Linetype1 'o' COLORS_OF_TIMES(ii)], ...
                'LineWidth', PROFILE_LW);
        end
    end
        if PLOT_R
            plot(mid_rho_local(1:end), ...
                J_inductive_halfmesh(index_to_plot, (1:end)) / ...
                J_inductive_halfmesh(end, (1)), ...
                'k--', ...
                'LineWidth', PROFILE_LW);
        else
            plot(mid_s_local(1:end), ...
                J_inductive_halfmesh(index_to_plot, (1:end)) / ...
                J_inductive_halfmesh(end, (1)), ...
                'k--', ...
                'LineWidth', PROFILE_LW);
        end
    
    %     xlabel('\rho', 'FontSize', PROFILE_FS);
    xlabel(xlabels, 'FontSize', PROFILE_FS);
    title('J_{ind} / J_{ind, s=0, t=tend}', 'FontSize', PROFILE_FS);
    ylabel('Normalized Current Density', 'FontSize', PROFILE_FS);
    legend([num2str(TIMES_TO_PLOT_IOTA')],'Location','best');
    xlim([0 1])

    figure;hold on;box on;
    plot(tspan, J_inductive_halfmesh(:, 1) / ...
        J_inductive_halfmesh(end, 1));
    xlabel('sec', 'FontSize', PROFILE_FS);
    title('J_{ind, s=0, t} / J_{ind, s=0, t=tend}', 'FontSize', PROFILE_FS);
    ylabel('Normalized Near-Axis Current Density', 'FontSize', PROFILE_FS);
    
    
    figure;hold on; box on;
    %     subplot(3,1,3);hold on;box on
    for ii = 1:length(INDICES_OF_TIMES_TO_PLOT_IOTA)
        index_to_plot = INDICES_OF_TIMES_TO_PLOT_IOTA(ii);
        if PLOT_R
            %             plot(rho_local(1:end), dF_pol(index_to_plot, (1:end))/Phi_LCFS_in, [Linetype1 Symboltype1 COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
            %             plot(mid_rho_local(1:end), dF_pol_halfmesh(index_to_plot, (1:end))/Phi_LCFS_in, [Linetype1 'o' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
            plot(mid_rho_local(1:end), J_pol_halfmesh(index_to_plot, (1:end)), [Linetype1 'x' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
        else
            %             plot(s_local(1:end), dF_pol(index_to_plot, (1:end))/Phi_LCFS_in, [Linetype1 Symboltype1 COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
            %             plot(mid_s_local(1:end), dF_pol_halfmesh(index_to_plot, (1:end))/Phi_LCFS_in, [Linetype1 'o' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
            plot(mid_s_local(1:end), J_pol_halfmesh(index_to_plot, (1:end)), [Linetype1 'x' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
        end
    end
    %     xlabel('\rho', 'FontSize', PROFILE_FS);
    xlabel(xlabels, 'FontSize', PROFILE_FS);
    title('J_{pol}', 'FontSize', PROFILE_FS);
    ylabel('A / Wb', 'FontSize', PROFILE_FS);
    legend(num2str(TIMES_TO_PLOT_IOTA'),'Location','best');
    xlim([0 1])
    
    figure;hold on; box on;
    %     subplot(3,1,2);hold on;box on
    for ii = 1:length(INDICES_OF_TIMES_TO_PLOT_IOTA)
        index_to_plot = INDICES_OF_TIMES_TO_PLOT_IOTA(ii);
        if PLOT_R
            plot(rho_local(1:end), PsiPrime_local(index_to_plot, (1:end)), [Linetype1 Symboltype1 COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
            %             plot(mid_rho_local(1:end), PsiPrime_halfmesh(index_to_plot, (1:end))/Phi_LCFS_in, [Linetype1 'o' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
        else
            plot(s_local(1:end), PsiPrime_local(index_to_plot, (1:end)), [Linetype1 Symboltype1 COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
            %             plot(mid_s_local(1:end), PsiPrime_halfmesh(index_to_plot, (1:end))/Phi_LCFS_in, [Linetype1 'o' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
        end
    end
    %     xlabel('\rho', 'FontSize', PROFILE_FS);
    xlabel(xlabels, 'FontSize', PROFILE_FS);
    title('Radial Derivatives of +(Pol Flux) and -(Tor Flux)', 'FontSize', PROFILE_FS);
    ylabel('Wb / s    s = \Psi / \Psi_{LCFS}', 'FontSize', PROFILE_FS);
    legend(num2str(TIMES_TO_PLOT_IOTA'),'Location','best');
    %     axis([0 1 .8 1.2]);
    ii = 1;
    if PLOT_R
        plot(rho_local(1:end), -PhiPrime_local((1:end)), [Linetype1 Symboltype1 COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
        %         plot(mid_rho_local(1:end), -PhiPrime_halfmesh((1:end)), [Linetype1 'o' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
    else
        plot(s_local(1:end), -PhiPrime_local((1:end)), [Linetype1 Symboltype1 COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
        %         plot(mid_s_local(1:end), -PhiPrime_halfmesh((1:end)), [Linetype1 'o' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
    end
    xlim([0 1])
    
    
    figure;box on;hold on;
    for ii = 1:length(INDICES_OF_TIMES_TO_PLOT_IOTA)
        index_to_plot = INDICES_OF_TIMES_TO_PLOT_IOTA(ii);
        if PLOT_R
            %             plot(rho_local(1:end), EdotB_local(index_to_plot, (1:end)), ...
            %                 [Linetype1 Symboltype1 COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
            plot(mid_rho_local(1:end), EdotB_halfmesh(index_to_plot,(1:end)), ...
                [Linetype1 'o' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
        else
            %             plot(s_local(1:end), EdotB_local(index_to_plot, (1:end)), ...
            %                 [Linetype1 Symboltype1 COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
            plot(mid_s_local(1:end), EdotB_halfmesh(index_to_plot,(1:end)), ...
                [Linetype1 'o' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
        end
    end
    %     xlabel('\rho', 'FontSize', PROFILE_FS);
    xlabel(xlabels, 'FontSize', PROFILE_FS);
    title('\langleE\cdotB\rangle', 'FontSize', PROFILE_FS)
    ylabel('T-V / m', 'FontSize', PROFILE_FS)
    %     axis([0 1 -1e-2 1e-2])
    xlim([0 1])
    legend(num2str(TIMES_TO_PLOT_IOTA'),'Location','best');
    
    figure;% hold on; box on;
    
    subplot(3,1,1);hold on;box on
    for ii = 1:length(INDICES_OF_TIMES_TO_PLOT_IOTA)
        index_to_plot = INDICES_OF_TIMES_TO_PLOT_IOTA(ii);
        if PLOT_R
            plot(mid_rho_local(1:end), QL1_halfmesh(index_to_plot, (1:end)), ['.' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
            plot(mid_rho_local(1:end), QR1_halfmesh(index_to_plot, (1:end)), ['o' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
        else
            plot(mid_s_local(1:end), QL1_halfmesh(index_to_plot, (1:end)), ['.' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
            plot(mid_s_local(1:end), QR1_halfmesh(index_to_plot, (1:end)), ['o' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
        end
    end
    %     xlabel('\rho', 'FontSize', PROFILE_FS);
    xlabel(xlabels, 'FontSize', PROFILE_FS);
    legend('\langleB^2\rangle V'' / \mu_0', '(I\Psi'' + F\Phi'')');
    
    subplot(3,2,3);hold on;box on
    for ii = 1:length(INDICES_OF_TIMES_TO_PLOT_IOTA)
        index_to_plot = INDICES_OF_TIMES_TO_PLOT_IOTA(ii);
        if PLOT_R
            plot(mid_rho_local(1:end), QL2_halfmesh(index_to_plot, (1:end)), ['.' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
            plot(mid_rho_local(1:end), QR2_halfmesh(index_to_plot, (1:end)), ['o' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
        else
            plot(mid_s_local(1:end), QL2_halfmesh(index_to_plot, (1:end)), ['.' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
            plot(mid_s_local(1:end), QR2_halfmesh(index_to_plot, (1:end)), ['o' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
        end
    end
    %     xlabel('\rho', 'FontSize', PROFILE_FS);
    xlabel(xlabels, 'FontSize', PROFILE_FS);
    legend('\langleJ_{BS}\cdotB\rangle V''', '\mu_0 (F I'' - I F'') == \langleJ\cdotB\rangle V''');
    
    subplot(3,2,4);hold on;box on
    for ii = 1:length(INDICES_OF_TIMES_TO_PLOT_IOTA)
        index_to_plot = INDICES_OF_TIMES_TO_PLOT_IOTA(ii);
        if PLOT_R
            plot(mid_rho_local(1:end), EdotB_halfmesh(index_to_plot, (1:end)), ...
                ['o' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
        else
            plot(mid_s_local(1:end), EdotB_halfmesh(index_to_plot, (1:end)), ...
                ['o' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
        end
    end
    %     xlabel('\rho', 'FontSize', PROFILE_FS);
    xlabel(xlabels, 'FontSize', PROFILE_FS);
    legend('\langleE\cdotB\rangle');
    
    subplot(3,1,3);hold on;box on
    for ii = 1:length(INDICES_OF_TIMES_TO_PLOT_IOTA)
        index_to_plot = INDICES_OF_TIMES_TO_PLOT_IOTA(ii);
        if PLOT_R
            plot(mid_rho_local(1:end), QL3_halfmesh(index_to_plot, (1:end)), ['.' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
            plot(mid_rho_local(1:end), QR3_halfmesh(index_to_plot, (1:end)), ['o' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
        else
            plot(mid_s_local(1:end), QL3_halfmesh(index_to_plot, (1:end)), ['.' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
            plot(mid_s_local(1:end), QR3_halfmesh(index_to_plot, (1:end)), ['o' COLORS_OF_TIMES(ii)], 'LineWidth', PROFILE_LW);
        end
    end
    %     xlabel('\rho', 'FontSize', PROFILE_FS);
    xlabel(xlabels, 'FontSize', PROFILE_FS);
    legend('-p'' V''', 'I''\Psi'' + F''\Phi''');
    
end
disp('Done with iota evolution view');

end



