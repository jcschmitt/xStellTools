function compare_Poincare(filename_list)

if nargin < 2
    color_list = 'ykmgbr';
end

fh = figure;
hold on;
box on;

for ii = 1:length(filename_list)
    filename = filename_list{ii}
    [R, Z, T, PolFlux] = read_Poincare(filename);
    the_color = color_list(mod(ii,6)+1);
    figure(fh);
    plot(R, Z, 'Marker', '.', 'color', the_color, 'LineStyle', 'None');
end

axis equal
legend(filename_list)



