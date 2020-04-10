function make_my_plot_pretty3

mycurrent_fig = gcf;
set(mycurrent_fig,'Units', 'points')
set(mycurrent_fig, 'Position', [200 200 400 350])
set(legend, 'location', 'best');


textobj = findobj(mycurrent_fig, 'type', 'text');
set(textobj, 'FontSize', 14);
set(textobj, 'FontWeight', 'normal');
lineobj = findobj(mycurrent_fig, 'type', 'line');
set(lineobj, 'linewidth', 2);

% figure children
child_list = get(mycurrent_fig, 'Children');

for ii = child_list'
    try    mycurrent_axes = ii;
        mycurrent_yaxis = get(mycurrent_axes, 'YLabel');
        mycurrent_xaxis = get(mycurrent_axes, 'XLabel');
        set(mycurrent_axes, 'FontSize', 14);
        set(mycurrent_axes, 'FontWeight', 'normal');
        set(mycurrent_xaxis, 'FontSize', 14);
        set(mycurrent_xaxis, 'FontWeight', 'normal');
        set(mycurrent_yaxis, 'FontSize', 14);
        set(mycurrent_yaxis, 'FontWeight', 'normal');
    catch
    end
    
end
