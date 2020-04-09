function vmec_data=make_well_plots(VMEC_FILENAME)

% ==========================================================
% ==========================================================
% ==========================================================

num_configs = length(VMEC_FILENAME);
plot_data_beta = [];
plot_data_max_gammatau = [];
plot_data_well_depth_mid_s = [];

plot_index = [];
for ii = 1:num_configs
    ii
    try
        vmec_data{ii} = load_vmec(VMEC_FILENAME{ii});
        
        
        s{ii} = vmec_data{ii}.phi ./ vmec_data{ii}.phi(end);
        rho{ii} = sqrt(s{ii});
        a0{ii} = vmec_data{ii}.Aminor;
        ns = length(s{ii});
        % [~,nsm1]=size(cobradata{ii}.grate);
        sall{ii} = linspace(0,1,ns);
        splot{ii} = sall{ii}(2:end);
        rhoplot{ii} = sqrt(splot{ii});
        plot_index = [plot_index ii];
        
        COMMENTS{ii} = ['\beta = ' num2str(vmec_data{ii}.betatot * 100) ' %'];
        jdotbCOMMENTS{ii} = ['\beta = ' num2str(vmec_data{ii}.betatot * 100) ' %, ' ...
            'Itor = ' num2str(vmec_data{ii}.ctor/1e3) ' kA'];
        
        plot_data_beta = [plot_data_beta vmec_data{ii}.betatot];
        
        ind_mid_s = round(ns * 0.5);
        
        plot_data_well_depth_mid_s = [plot_data_well_depth_mid_s vmec_data{ii}.welldepth(ind_mid_s)];
%         
%         try
%             cobradata{ii} = read_cobra(COBRA_FILENAME{ii});
%             plot_data_max_gammatau = [plot_data_max_gammatau max(max(cobradata{ii}.grate))];
%         catch
%             disp(['<----Did not find some cobra data for index #' num2str(ii)]);
%         end
        
        
    catch
        disp(['<----Did not find some data for index #' num2str(ii)]);
    end
    
end
scale_yaxis = 10;


linewidth = 2;
colors{1} = 'k';
colors{2} = [20 20 235 ]/255;
colors{3} = [235 20 20]/255;
colors{4} = [20 220 220]/255;
colors{5} = [150 150 215]/255;


%plot_index = 1:num_configs;
%plot_index = 1:16;





% ==========================================================
% ==========================================================
% ==========================================================

disp('<----Here come the figures!')

% fh1 = figure;
%
% for ii = plot_index
%     ii
%     %keyboard
%     plot_indices = [2:(length(rho{ii})-1)];
%
%     subplot(2,3,1);box on;hold on; grid on;
%     plot(rho{ii}(plot_indices), vmec_data{ii}.DCurr(plot_indices),  'Linewidth', linewidth);
%     ylim(scale_yaxis*[-1.5 1]);
%     title('DCurr');
%     xlabel('$\rho = \sqrt{\psi_{tor} / \psi_{tor,LCFS}}$', 'Interpreter', 'latex')
%
%     subplot(2,3,2);box on;hold on;
%     plot(rho{ii}(plot_indices), vmec_data{ii}.DGeod(plot_indices),  'Linewidth', linewidth);
%     ylim(scale_yaxis*[-1.5 1]);
%     title('DGeod');
%     xlabel('$\rho = \sqrt{\psi_{tor} / \psi_{tor,LCFS}}$', 'Interpreter', 'latex')
%
%     subplot(2,3,4);box on;hold on;
%     plot(rho{ii}(plot_indices), vmec_data{ii}.DShear(plot_indices), 'Linewidth', linewidth);
%     ylim(scale_yaxis*[-1.5 1]);
%     title('DShear');
%     xlabel('$\rho = \sqrt{\psi_{tor} / \psi_{tor,LCFS}}$', 'Interpreter', 'latex')
%
%     subplot(2,3,5);box on;hold on;
%     plot(rho{ii}(plot_indices), vmec_data{ii}.DWell(plot_indices),  'Linewidth', linewidth);
%     ylim(scale_yaxis*[-1.5 1]);
%     title('DWell');
%     xlabel('$\rho = \sqrt{\psi_{tor} / \psi_{tor,LCFS}}$', 'Interpreter', 'latex')
%
%     subplot(1,3,3);box on;hold on;
%     plot(rho{ii}(plot_indices), vmec_data{ii}.DMerc(plot_indices), 'Linewidth', linewidth);
%     ylim(scale_yaxis*[-1.5 1]);
%     title('DMerc');
%     xlabel('$\rho = \sqrt{\psi_{tor} / \psi_{tor,LCFS}}$', 'Interpreter', 'latex')
%
%     %     subplot(2,3,6);box on;hold on;
%     %     plot(rhoplot{ii}, cobradata{ii}.grate',  'Marker', 'o', 'Linewidth', linewidth);
%     %     %ylim(scale_yaxis*[-1.5 1]);
%     %     title('\gamma \tau_A');
%     %     xlabel('$\rho = \sqrt{\psi_{tor} / \psi_{tor,LCFS}}$', 'Interpreter', 'latex')
%
% end
% %
% % subplot(2,3,1)
% % legend(LEGEND_LIST{plot_index})
% % grid on
% % subplot(2,3,2)
% % grid on
% % legend(LEGEND_LIST{plot_index})
% % grid on
% % subplot(2,3,4)
% % grid on
% % legend(LEGEND_LIST{plot_index})
% % grid on
% % subplot(2,3,5)
% % grid on
% % legend(LEGEND_LIST{plot_index})
% % subplot(1,3,3)
% % grid on
% % legend(COMMENTS{plot_index})
% %
%
% subplot(2,3,1)
% legend(COMMENTS{plot_index})
% grid on
% subplot(2,3,2)
% grid on
% legend(COMMENTS{plot_index})
% grid on
% subplot(2,3,4)
% grid on
% legend(COMMENTS{plot_index})
% grid on
% subplot(2,3,5)
% grid on
% legend(COMMENTS{plot_index})
% subplot(1,3,3)
% grid on
% legend(COMMENTS{plot_index})
% % subplot(2,3,6)
% % grid on
% % legend(COMMENTS{plot_index})
%
%
% fh2 = figure;
%
% for ii = plot_index
%     %keyboard
%     plot_indices = [2:(length(s{ii})-1)];
%
%     subplot(1,3,3);box on;hold on;
%     plot(s{ii}(plot_indices), vmec_data{ii}.DMerc(plot_indices),  'Linewidth', linewidth);
%     ylim(scale_yaxis*[-1.5 1]);
%     title('DMerc');
%
%     subplot(2,3,1);box on;hold on;
%     plot(s{ii}(plot_indices), vmec_data{ii}.DCurr(plot_indices),  'Linewidth', linewidth);
%     ylim(scale_yaxis*[-1.5 1]);
%     title('DCurr');
%
%     subplot(2,3,2);box on;hold on;
%     plot(s{ii}(plot_indices), vmec_data{ii}.DGeod(plot_indices),   'Linewidth', linewidth);
%     ylim(scale_yaxis*[-1.5 1]);
%     title('DGeod');
%
%     subplot(2,3,4);box on;hold on;
%     plot(s{ii}(plot_indices), vmec_data{ii}.DShear(plot_indices),   'Linewidth', linewidth);
%     ylim(scale_yaxis*[-1.5 1]);
%     title('DShear');
%
%     subplot(2,3,5);box on;hold on;
%     plot(s{ii}(plot_indices), vmec_data{ii}.DWell(plot_indices),   'Linewidth', linewidth);
%     ylim(scale_yaxis*[-1.5 1]);
%     title('DWell');
%
%     %     subplot(2,3,6);box on;hold on;
%     %     plot(splot{ii}, cobradata{ii}.grate', 'Marker', 'o', 'Linewidth', linewidth);
%     %     %ylim(scale_yaxis*[-1.5 1]);
%     %     title('\gamma \tau_A');
%
%
%     xlabel('$s = \psi_{tor} / \psi_{tor,LCFS}$', 'Interpreter', 'latex')
%
% end
%
% % subplot(2,3,1)
% % legend(LEGEND_LIST{plot_index})
% % grid on
% % subplot(2,3,2)
% % grid on
% % legend(LEGEND_LIST{plot_index})
% % grid on
% % subplot(2,3,4)
% % grid on
% % legend(LEGEND_LIST{plot_index})
% % grid on
% % subplot(2,3,5)
% % grid on
% % legend(LEGEND_LIST{plot_index})
% % subplot(1,3,3)
% % grid on
% % legend(COMMENTS{plot_index})
%
% subplot(2,3,1)
% legend(COMMENTS{plot_index})
% grid on
% subplot(2,3,2)
% grid on
% legend(COMMENTS{plot_index})
% grid on
% subplot(2,3,4)
% grid on
% legend(COMMENTS{plot_index})
% grid on
% subplot(2,3,5)
% grid on
% legend(COMMENTS{plot_index})
% subplot(1,3,3)
% grid on
% legend(COMMENTS{plot_index})
% % subplot(2,3,6)
% % grid on
% % legend(COMMENTS{plot_index})
%
%
% fh3 = figure;
% sqr_grid_dim = ceil(sqrt(length(plot_index)));
% cols = sqr_grid_dim;
% rows = ceil(length(plot_index) / cols);
% %keyboard; % modify cols, rows, as neccessary.
%
% counter = 0;
% legend_str = {};
% for ii = plot_index
%     counter = counter +1;
%     plot_indices = [1:(length(rhoplot{ii})-1)];
%
%     subplot(rows, cols, counter)
%
%     try
%         plot(rhoplot{ii}(plot_indices), cobradata{ii}.grate(:,plot_indices)', '--', 'Marker', 'None', 'Linewidth', linewidth);
%         ylabel('\gamma \tau_A');
%         title(COMMENTS{ii})
%         xlabel('$\rho = \sqrt{\psi_{tor} / \psi_{tor,LCFS}}$', 'Interpreter', 'latex')
%         grid on;
%         for jj = 1:length(cobradata{ii}.grate(:,1))
%             if jj == 1
%                 legend_str{jj} = ['(\zeta = ' num2str(cobradata{ii}.zeta(jj)) ...
%                     ', \theta = ' num2str(cobradata{ii}.theta(jj)) ')'];
%             else
%                 legend_str{jj} = ['(' num2str(cobradata{ii}.zeta(jj)) ...
%                     ', ' num2str(cobradata{ii}.theta(jj)) ')'];
%             end
%         end
%     catch
%     end
%
%     end
%
%     legend(legend_str);
%
% fh4 = figure;
%
% counter = 0;
% for ii = plot_index
%     counter = counter +1;
%     plot_indices = [1:(length(splot{ii})-1)];
%
%     subplot(rows, cols, counter)
%     try
%
%     plot(splot{ii}(plot_indices), cobradata{ii}.grate(:,plot_indices)', '--', 'Marker', 'None', 'Linewidth', linewidth);
%     catch
%     end
%     ylabel('\gamma \tau_A');
%     title(COMMENTS{ii})
%     grid on;
%     xlabel('$s = \psi_{tor} / \psi_{tor,LCFS}$', 'Interpreter', 'latex')
% end
% legend(legend_str);
%
%
% figure;
%
%
% for ii = plot_index
%     hold on; box on; grid on;
%     plot(rho{ii}, vmec_data{ii}.iotaf  , 'Marker', 'None', 'Linewidth', linewidth);
%     ylabel('\iota_f');
%     xlabel('$\rho = \sqrt{\psi_{tor} / \psi_{tor,LCFS}}$', 'Interpreter', 'latex')
%
% end
% legend(COMMENTS{plot_index})
%
% figure;
%
%
% for ii = plot_index
%     hold on; box on; grid on;
%     plot(rho{ii}, vmec_data{ii}.presf/1e3  , 'Marker', 'None', 'Linewidth', linewidth);
%     ylabel('Pressure [kPa]');
%     xlabel('$\rho = \sqrt{\psi_{tor} / \psi_{tor,LCFS}}$', 'Interpreter', 'latex')
%
% end
%
% legend(COMMENTS{plot_index})
% figure;
%
%
% for ii = plot_index
%     hold on; box on; grid on;
%     plot(rho{ii}, vmec_data{ii}.jdotb/1e3  , 'Marker', 'None', 'Linewidth', linewidth);
%     ylabel('<j . b>, kA');
%     xlabel('$\rho = \sqrt{\psi_{tor} / \psi_{tor,LCFS}}$', 'Interpreter', 'latex')
%
% end
% legend(jdotbCOMMENTS{plot_index})
%
% figure
%
% for ii = plot_index
%     hold on; box on; grid on;
%     plot(s{ii}, vmec_data{ii}.iotaf  , 'Marker', 'None', 'Linewidth', linewidth);
%     ylabel('\iota_f');
%     xlabel('$s = \psi_{tor} / \psi_{tor,LCFS}$', 'Interpreter', 'latex')
%
% end
% legend(COMMENTS{plot_index})
%
% figure;
%
%
% for ii = plot_index
%     hold on; box on; grid on;
%     plot(s{ii}, vmec_data{ii}.presf/1e3  , 'Marker', 'None', 'Linewidth', linewidth);
%     ylabel('Pressure [kPa]');
%     xlabel('$s = \psi_{tor} / \psi_{tor,LCFS}$', 'Interpreter', 'latex')
%
% end
% legend(COMMENTS{plot_index})
%
% figure;
%
%
% for ii = plot_index
%     hold on; box on; grid on;
%     plot(s{ii}, vmec_data{ii}.jdotb/1e3  , 'Marker', 'None', 'Linewidth', linewidth);
%     ylabel('<j . b>, kA');
%     xlabel('$s = \psi_{tor} / \psi_{tor,LCFS}$', 'Interpreter', 'latex')
%
% end
% legend(jdotbCOMMENTS{plot_index})
%
%
%
% figure
% box on;hold on
% plot(100*plot_data_beta,plot_data_max_gammatau, 'o');
% xlabel('Beta %')
% ylabel('max(\gamma \tau)');
%


figure

well_depth_mid_s = [];

for ii = plot_index
    hold on; box on; grid on;
    plot(s{ii}, vmec_data{ii}.welldepth  , 'Marker', 'None', 'Linewidth', linewidth);
    ylabel('Well depth');
    xlabel('$s = \psi_{tor} / \psi_{tor,LCFS}$', 'Interpreter', 'latex')
    
    
end
legend(COMMENTS{plot_index})

figure

for ii = plot_index
    hold on; box on; grid on;
    plot(a0{ii} * sqrt(s{ii}), vmec_data{ii}.welldepth  , 'Marker', 'None', 'Linewidth', linewidth);
    ylabel('Well depth');
    xlabel('$r_{eff}  m$', 'Interpreter', 'latex')
    
    
end
legend(COMMENTS{plot_index})


figure
box on;hold on
plot(100*plot_data_beta,plot_data_well_depth_mid_s, 'o');
xlabel('Beta %')
ylabel('Well Depth @ s=0.5');


