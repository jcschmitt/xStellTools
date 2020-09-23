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
    rhoplot=sqrt(splot);

    %plot(rhoplot, cobradata{ii}.grate', 'o');
    plot(rhoplot, cobradata{ii}.grate', '--', 'Linewidth', 3);
    title(strrep( ([title_prefix, ' ', cobra_extensions{ii}]), '_', '\_'));
    
end
xlabel('$\rho$', 'Interpreter', 'Latex')
ylabel('\gamma\tau');
grid on


