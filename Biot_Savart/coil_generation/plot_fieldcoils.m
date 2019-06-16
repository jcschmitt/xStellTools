function plot_coils(Coil_in, colorcode)
% plot a coil on a figure that is assumed to have the focus
for ii = 1:Coil_in.num_turns
    if length(colorcode) > 1
        plot3(Coil_in.turn_number(ii).x, Coil_in.turn_number(ii).y, ...
            Coil_in.turn_number(ii).z, 'Color', colorcode{ii});
    else
        plot3(Coil_in.turn_number(ii).x, Coil_in.turn_number(ii).y, ...
            Coil_in.turn_number(ii).z, 'Color', colorcode);
    end
end