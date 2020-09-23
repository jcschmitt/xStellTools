function plot_cobra(cobra_extensions, colorcode, figure_number, title_prefix)
%plot_cobra(cobra_extensions, title_prefix)
% Plot COBRA output
if nargin < 2
    colorcode = 'k';
end
if nargin < 3
    figure_number = figure;
end
if nargin < 4
    title_prefix = '';
end

figure(figure_number);box on; hold on;
for ii = 1:length(cobra_extensions)
    cobradata{ii} = read_cobra(['cobra_grate.' cobra_extensions{ii}]);
    
    [~,nsm1]=size(cobradata{ii}.grate);
    ns = nsm1+1;
    sall = linspace(0,1,ns);
    splot = sall(2:end);
    delta_splot = sall(2) - sall(1);
    %rhoplot=sqrt(splot);
    [grate_rows, grate_cols] = size(cobradata{ii}.grate)
    num_inner_pts = grate_rows;
    %splot_new = zeros(1, num_inner_pts);
    inner_span = linspace(0,1, num_inner_pts+1);
    inner_span = delta_splot * inner_span(1:(end-1));
    for jj_use = 1:nsm1
        %splot_new(((jj-1)*num_inner_pts + 1):(jj*num_inner_pts)) = splot(jj) + inner_span;
        splot_new = splot(jj_use) + inner_span;
        cobra_new = cobradata{ii}.grate(:,jj_use);
        
        %plot(rhoplot, cobradata{ii}.grate', 'o');
        plot(splot_new, cobra_new, '.', 'MarkerSize', 14, 'Color', colorcode);
        %plot(rhoplot, cobradata{ii}.grate', '--', 'Linewidth', 3);
    end
    title(strrep( ([title_prefix, ' ', cobra_extensions{ii}]), '_', '\_'));
    
end
xlabel('s ', 'Interpreter', 'Latex')
ylabel('\gamma\tau');
grid on



%figure;box on; hold on;
for ii = 1:length(cobra_extensions)
    cobradata{ii} = read_cobra(['cobra_grate.' cobra_extensions{ii}]);
    
    [~,nsm1]=size(cobradata{ii}.grate);
    ns = nsm1+1;
    sall = linspace(0,1,ns);
    splot = sall(2:end);
    %rhoplot=sqrt(splot);
    %     contourf(cobradata{ii}.theta, cobradata{ii}.zeta, cobradata{ii}.grate', '--', 'Linewidth', 3);
    %[   [tgrid, zgrid] = meshgrid(cobradata{ii}.theta, cobradata{ii}.zeta)
    num_zeta = length(unique(cobradata{ii}.zeta));
    num_theta = length(unique(cobradata{ii}.theta));
    tgrid = reshape(cobradata{ii}.theta, num_zeta,num_theta);
    zgrid = reshape(cobradata{ii}.zeta, num_zeta,num_theta);
    for jj = ns*[0.3, 0.5, 0.7]
        jj_use = round(jj)
        figure
        grate_plot= reshape(cobradata{ii}.grate(:,jj_use), num_zeta,num_theta);
        %plot(rhoplot, cobradata{ii}.grate', 'o');
        contourf(zgrid, tgrid, grate_plot, '--', 'Linewidth', 3);
        title(strrep( ([title_prefix, ' ', cobra_extensions{ii} '  ns=', num2str(jj_use)]), '_', '\_'));
        colorbar
        xlabel('$\zeta$', 'Interpreter', 'Latex')
        ylabel('\theta');
        grid on
        make_my_plot_pretty3
        caxis([-.2, .2])
    end
    
end


