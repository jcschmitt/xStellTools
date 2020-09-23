function plot_cobra(cobra_extensions, title_prefix)
%plot_cobra(cobra_extensions, title_prefix)
% Plot COBRA output
if nargin < 2
    title_prefix = '';
end

figure;box on; hold on; 
for ii = 1:length(cobra_extensions)
    cobradata{ii} = read_cobra(['cobra_grate.' cobra_extensions{ii}]);

    [~,nsm1]=size(cobradata{ii}.grate);
    ns = nsm1+1;
    sall = linspace(0,1,ns);
    splot = sall(2:end);
    delta_splot = sall(2) - sall(1);
    rhoplot=sqrt(splot);
    [grate_rows, grate_cols] = size(cobradata{ii}.grate)
    num_inner_pts = grate_rows;
    %splot_new = zeros(1, num_inner_pts);
    inner_span = linspace(0,1, num_inner_pts+1);
    inner_span = delta_splot * inner_span(1:(end-1));
    for jj = 1:nsm1
        %splot_new(((jj-1)*num_inner_pts + 1):(jj*num_inner_pts)) = splot(jj) + inner_span;
        splot_new = splot(jj) + inner_span;
        cobra_new = cobradata{ii}.grate(:,jj);
        
        %plot(rhoplot, cobradata{ii}.grate', 'o');
        plot(splot_new, cobra_new, 'b.', 'MarkerSize', 10);
        %plot(rhoplot, cobradata{ii}.grate', '--', 'Linewidth', 3);
    end
    title(strrep( ([title_prefix, ' ', cobra_extensions{ii}]), '_', '\_'));
    
end
xlabel('$\rho$', 'Interpreter', 'Latex')
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
    tgrid = reshape(cobradata{ii}.theta, 5,5)
    zgrid = reshape(cobradata{ii}.zeta, 5,5)
    for jj = [1:5:nsm1]
        figure
        grate_plot= reshape(cobradata{ii}.grate(:,jj), 5,5);
        %plot(rhoplot, cobradata{ii}.grate', 'o');
        contourf(tgrid, zgrid, grate_plot, '--', 'Linewidth', 3);
        title(strrep( ([title_prefix, ' ', cobra_extensions{ii} '  ns=', num2str(jj)]), '_', '\_'));
        colorbar
    end
    xlabel('$\theta$', 'Interpreter', 'Latex')
    ylabel('\zeta');
    grid on

end


