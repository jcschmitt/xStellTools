function make_my_plot_pretty

mycurrent_fig = gcf;
set(mycurrent_fig,'Units', 'points')
set(mycurrent_fig, 'Position', [200 200 300 187.5])
% set(mycurrent_fig, 'Position', [200 200 400 250])

textobj = findobj(mycurrent_fig, 'type', 'text');
% set(textobj, 'fontunits', 'points');
set(textobj, 'FontSize', 14);
set(textobj, 'FontWeight', 'bold');
lineobj = findobj(mycurrent_fig, 'type', 'line');
set(lineobj, 'linewidth', 2);

% figure children
child_list = get(mycurrent_fig, 'Children');

for ii = child_list'
    try
        mycurrent_axes = ii;
        mycurrent_yaxis = get(mycurrent_axes, 'YLabel');
        mycurrent_xaxis = get(mycurrent_axes, 'XLabel');
        
        % set(mycurrent_axes, 'FontSize', 12)
        % set(mycurrent_axes, 'FontWeight', 'bold')
        
        
        
        set(mycurrent_axes, 'FontSize', 12);
        set(mycurrent_axes, 'FontWeight', 'bold');
        set(mycurrent_xaxis, 'FontSize', 12);
        set(mycurrent_xaxis, 'FontWeight', 'bold');
        set(mycurrent_yaxis, 'FontSize', 12);
        set(mycurrent_yaxis, 'FontWeight', 'bold');
        % set(mycurrent_fig, 'Units', 'points')
    catch
    end
    
end
