function blah

close all

% Scan of am(17)
am{1} = [0  2  1    1];
am{2} = [0  1  1    1];
am{3} = [0  .5  1   1];
am{4} = [0  .1  1   1];
make_plots(am, 'am(17)')

clear am
% Scan of am(18)
am{1} = [0  1  -1    .5];
am{2} = [0  1  -.5   .5];
am{3} = [0  1  -.1  .5];
am{4} = [0  1  0   .5];
am{5} = [0  1  .1   .5];
am{6} = [0  1  .5   .5];
am{7} = [0  1  1   .5];
make_plots(am, 'am(18)')

clear am
% Scan of am(19)
am{1} = [0  1  1    2];
am{2} = [0  1  1    1];
am{3} = [0  1  1    .5];
am{4} = [0  1  1    .1];
make_plots(am, 'am(19)')

clear am
% Scan of am(18) and am(19)
am{1} = [0  1  .9    .1];
am{2} = [0  1  .95   .1];
am{3} = [0  1  .99  .1];
am{4} = [0  1  .999   .1];
am{5} = [0  1  .9999   .1];
am{6} = [0  1  1   .1];
am{7} = [0  1  .9    .2];
am{8} = [0  1  .95   .2];
am{9} = [0  1  .99  .2];
am{10} = [0  1  .999   .2];
am{11} = [0  1  .9999   .2];
am{12} = [0  1  1   .2];
am{13} = [0  1  .9    .5];
am{14} = [0  1  .95   .5];
am{15} = [0  1  .99  .5];
am{16} = [0  1  .999   .5];
am{17} = [0  1  .9999   .5];
am{18} = [0  1  1   .5];
make_plots(am, 'am(18) and am(19)')

function make_plots(am, title_in)
ss = linspace(0,1,51);

for ii = 1:length(am)
    pres_out{ii} = pedestal_pressure(am{ii}, ss);
end

figure;box on;hold on;title(title_in)
for ii = 1:length(am)
    plot(ss, pres_out{ii})
end
axis tight
legend(num2str([1:length(am)]'))
make_my_plot_pretty3

function pressure = pedestal_pressure(am, ss_in)
if am(4) <= 0
    am(1:5) = 0;
    am(4) = 1e30;
else
    am(5) = 1.0 / (tanh(2*am(3)/am(4)) - tanh(2*(am(3)-1) / am(4)));
end

pressure = 0 * ss_in;
pressure = pressure + am(5) * am(2) * (tanh(2*(am(3) - sqrt(ss_in)) / ...
    am(4)) - tanh(2*(am(3) - 1.0) / am(4) ) );

