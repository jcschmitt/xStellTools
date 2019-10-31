function make_merc_balloon_plots(CONFIG_LIST, VMEC_FILENAME, ...
    COBRA_FILENAME)

% ==========================================================
% ==========================================================
% ==========================================================

num_configs = length(CONFIG_LIST);
plot_index = [];
for ii = 1:num_configs
    ii
    try
        vmec_data{ii} = load_mercier(VMEC_FILENAME{ii});
        COMMENTS{ii} = ['\beta = ' num2str(vmec_data{ii}.betatot * 100) ' %'];
        jdotbCOMMENTS{ii} = ['\beta = ' num2str(vmec_data{ii}.betatot * 100) ' %, ' ...
                             'Itor = ' num2str(vmec_data{ii}.ctor/1e3) ' kA'];
        
        cobradata{ii} = read_cobra(COBRA_FILENAME{ii});
        
        s{ii} = vmec_data{ii}.phi ./ vmec_data{ii}.phi(end);
        rho{ii} = sqrt(s{ii});
        
        [~,nsm1]=size(cobradata{ii}.grate);
        ns = nsm1+1;
        sall{ii} = linspace(0,1,ns);
        splot{ii} = sall{ii}(2:end);
        rhoplot{ii} = sqrt(splot{ii});
        plot_index = [plot_index ii];
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

fh1 = figure;

for ii = plot_index
    ii
    %keyboard
    plot_indices = [2:(length(rho{ii})-1)];
    
    subplot(2,3,1);box on;hold on;
    plot(rho{ii}(plot_indices), vmec_data{ii}.DCurr(plot_indices),  'Linewidth', linewidth);
    ylim(scale_yaxis*[-1.5 1]);
    title('DCurr');
    xlabel('$\rho = \sqrt{\psi_{tor} / \psi_{tor,LCFS}}$', 'Interpreter', 'latex')
    
    subplot(2,3,2);box on;hold on;
    plot(rho{ii}(plot_indices), vmec_data{ii}.DGeod(plot_indices),  'Linewidth', linewidth);
    ylim(scale_yaxis*[-1.5 1]);
    title('DGeod');
    xlabel('$\rho = \sqrt{\psi_{tor} / \psi_{tor,LCFS}}$', 'Interpreter', 'latex')
    
    subplot(2,3,4);box on;hold on;
    plot(rho{ii}(plot_indices), vmec_data{ii}.DShear(plot_indices), 'Linewidth', linewidth);
    ylim(scale_yaxis*[-1.5 1]);
    title('DShear');
    xlabel('$\rho = \sqrt{\psi_{tor} / \psi_{tor,LCFS}}$', 'Interpreter', 'latex')
    
    subplot(2,3,5);box on;hold on;
    plot(rho{ii}(plot_indices), vmec_data{ii}.DWell(plot_indices),  'Linewidth', linewidth);
    ylim(scale_yaxis*[-1.5 1]);
    title('DWell');
    xlabel('$\rho = \sqrt{\psi_{tor} / \psi_{tor,LCFS}}$', 'Interpreter', 'latex')
    
    subplot(1,3,3);box on;hold on;
    plot(rho{ii}(plot_indices), vmec_data{ii}.DMerc(plot_indices), 'Linewidth', linewidth);
    ylim(scale_yaxis*[-1.5 1]);
    title('DMerc');
    xlabel('$\rho = \sqrt{\psi_{tor} / \psi_{tor,LCFS}}$', 'Interpreter', 'latex')
    
    %     subplot(2,3,6);box on;hold on;
    %     plot(rhoplot{ii}, cobradata{ii}.grate',  'Marker', 'o', 'Linewidth', linewidth);
    %     %ylim(scale_yaxis*[-1.5 1]);
    %     title('\gamma \tau_A');
    %     xlabel('$\rho = \sqrt{\psi_{tor} / \psi_{tor,LCFS}}$', 'Interpreter', 'latex')
    
end
%
% subplot(2,3,1)
% legend(LEGEND_LIST{plot_index})
% grid on
% subplot(2,3,2)
% grid on
% legend(LEGEND_LIST{plot_index})
% grid on
% subplot(2,3,4)
% grid on
% legend(LEGEND_LIST{plot_index})
% grid on
% subplot(2,3,5)
% grid on
% legend(LEGEND_LIST{plot_index})
% subplot(1,3,3)
% grid on
% legend(COMMENTS{plot_index})
%

subplot(2,3,1)
legend(COMMENTS{plot_index})
grid on
subplot(2,3,2)
grid on
legend(COMMENTS{plot_index})
grid on
subplot(2,3,4)
grid on
legend(COMMENTS{plot_index})
grid on
subplot(2,3,5)
grid on
legend(COMMENTS{plot_index})
subplot(1,3,3)
grid on
legend(COMMENTS{plot_index})
% subplot(2,3,6)
% grid on
% legend(COMMENTS{plot_index})


fh2 = figure;

for ii = plot_index
    %keyboard
    plot_indices = [2:(length(s{ii})-1)];
    
    subplot(1,3,3);box on;hold on;
    plot(s{ii}(plot_indices), vmec_data{ii}.DMerc(plot_indices),  'Linewidth', linewidth);
    ylim(scale_yaxis*[-1.5 1]);
    title('DMerc');
    
    subplot(2,3,1);box on;hold on;
    plot(s{ii}(plot_indices), vmec_data{ii}.DCurr(plot_indices),  'Linewidth', linewidth);
    ylim(scale_yaxis*[-1.5 1]);
    title('DCurr');
    
    subplot(2,3,2);box on;hold on;
    plot(s{ii}(plot_indices), vmec_data{ii}.DGeod(plot_indices),   'Linewidth', linewidth);
    ylim(scale_yaxis*[-1.5 1]);
    title('DGeod');
    
    subplot(2,3,4);box on;hold on;
    plot(s{ii}(plot_indices), vmec_data{ii}.DShear(plot_indices),   'Linewidth', linewidth);
    ylim(scale_yaxis*[-1.5 1]);
    title('DShear');
    
    subplot(2,3,5);box on;hold on;
    plot(s{ii}(plot_indices), vmec_data{ii}.DWell(plot_indices),   'Linewidth', linewidth);
    ylim(scale_yaxis*[-1.5 1]);
    title('DWell');
    
    %     subplot(2,3,6);box on;hold on;
    %     plot(splot{ii}, cobradata{ii}.grate', 'Marker', 'o', 'Linewidth', linewidth);
    %     %ylim(scale_yaxis*[-1.5 1]);
    %     title('\gamma \tau_A');
    
    
    xlabel('$s = \psi_{tor} / \psi_{tor,LCFS}$', 'Interpreter', 'latex')
    
end

% subplot(2,3,1)
% legend(LEGEND_LIST{plot_index})
% grid on
% subplot(2,3,2)
% grid on
% legend(LEGEND_LIST{plot_index})
% grid on
% subplot(2,3,4)
% grid on
% legend(LEGEND_LIST{plot_index})
% grid on
% subplot(2,3,5)
% grid on
% legend(LEGEND_LIST{plot_index})
% subplot(1,3,3)
% grid on
% legend(COMMENTS{plot_index})

subplot(2,3,1)
legend(COMMENTS{plot_index})
grid on
subplot(2,3,2)
grid on
legend(COMMENTS{plot_index})
grid on
subplot(2,3,4)
grid on
legend(COMMENTS{plot_index})
grid on
subplot(2,3,5)
grid on
legend(COMMENTS{plot_index})
subplot(1,3,3)
grid on
legend(COMMENTS{plot_index})
% subplot(2,3,6)
% grid on
% legend(COMMENTS{plot_index})


fh3 = figure;
sqr_grid_dim = ceil(sqrt(length(plot_index)));
cols = sqr_grid_dim;
rows = ceil(length(plot_index) / cols);

counter = 0;
for ii = plot_index
    counter = counter +1;
    plot_indices = [1:(length(rhoplot{ii})-1)];
    
    subplot(rows, cols, counter)
    
    
    plot(rhoplot{ii}(plot_indices), cobradata{ii}.grate(:,plot_indices)', 'Marker', 'o', 'Linewidth', linewidth);
    ylabel('\gamma \tau_A');
    title(COMMENTS{ii})
    xlabel('$\rho = \sqrt{\psi_{tor} / \psi_{tor,LCFS}}$', 'Interpreter', 'latex')
    
end


fh4 = figure;

counter = 0;
for ii = plot_index
    counter = counter +1;
    plot_indices = [1:(length(splot{ii})-1)];
    
    subplot(rows, cols, counter)
    
    
    plot(splot{ii}(plot_indices), cobradata{ii}.grate(:,plot_indices)', 'Marker', 'o', 'Linewidth', linewidth);
    ylabel('\gamma \tau_A');
    title(COMMENTS{ii})
    
    xlabel('$s = \psi_{tor} / \psi_{tor,LCFS}$', 'Interpreter', 'latex')
end


figure;


for ii = plot_index
    hold on; box on; grid on;
    plot(rho{ii}, vmec_data{ii}.iotaf  , 'Marker', 'None', 'Linewidth', linewidth);
    ylabel('\iota_f');
    xlabel('$\rho = \sqrt{\psi_{tor} / \psi_{tor,LCFS}}$', 'Interpreter', 'latex')
    
end
    legend(COMMENTS{plot_index})

figure;


for ii = plot_index
    hold on; box on; grid on;
    plot(rho{ii}, vmec_data{ii}.presf/1e3  , 'Marker', 'None', 'Linewidth', linewidth);
    ylabel('Pressure [kPa]');
    xlabel('$\rho = \sqrt{\psi_{tor} / \psi_{tor,LCFS}}$', 'Interpreter', 'latex')
    
end

    legend(COMMENTS{plot_index})
figure;


for ii = plot_index
    hold on; box on; grid on;
    plot(rho{ii}, vmec_data{ii}.jdotb/1e3  , 'Marker', 'None', 'Linewidth', linewidth);
    ylabel('<j . b>, kA');
    xlabel('$\rho = \sqrt{\psi_{tor} / \psi_{tor,LCFS}}$', 'Interpreter', 'latex')
    
end
    legend(jdotbCOMMENTS{plot_index})

figure

for ii = plot_index
    hold on; box on; grid on;
    plot(s{ii}, vmec_data{ii}.iotaf  , 'Marker', 'None', 'Linewidth', linewidth);
    ylabel('\iota_f');
    xlabel('$s = \psi_{tor} / \psi_{tor,LCFS}$', 'Interpreter', 'latex')
    
end
    legend(COMMENTS{plot_index})

figure;


for ii = plot_index
    hold on; box on; grid on;
    plot(s{ii}, vmec_data{ii}.presf/1e3  , 'Marker', 'None', 'Linewidth', linewidth);
    ylabel('Pressure [kPa]');
    xlabel('$s = \psi_{tor} / \psi_{tor,LCFS}$', 'Interpreter', 'latex')
    
end
    legend(COMMENTS{plot_index})

figure;


for ii = plot_index
    hold on; box on; grid on;
    plot(s{ii}, vmec_data{ii}.jdotb/1e3  , 'Marker', 'None', 'Linewidth', linewidth);
    ylabel('<j . b>, kA');
    xlabel('$s = \psi_{tor} / \psi_{tor,LCFS}$', 'Interpreter', 'latex')
    
end
    legend(jdotbCOMMENTS{plot_index})




