s = linspace(0,1,51);

P1{1} = pp_two_power_vmec(s, [1 .2 1]);
P1{2} = pp_two_power_vmec(s, [1 .5 1]);
P1{3} = pp_two_power_vmec(s, [1 1 1]);
P1{4} = pp_two_power_vmec(s, [1 2 1]);
P1{5} = pp_two_power_vmec(s, [1 5 1]);

P2{1} = pp_two_power_vmec(s, [1 1 .2]);
P2{2} = pp_two_power_vmec(s, [1 1 .5]);
P2{3} = pp_two_power_vmec(s, [1 1 1]);
P2{4} = pp_two_power_vmec(s, [1 1 2]);
P2{5} = pp_two_power_vmec(s, [1 1 5]);

for ii = 1:length(P1)
    dP1{ii} = [0 diff(P1{ii})];
end
for ii = 1:length(P2)
    dP2{ii} = [0 diff(P2{ii})];
end

figure;
subplot(2,1,1);box on;hold on;
plot(s, P1{1}, s, P1{2}, s, P1{3}, s, P1{4}, s, P1{5});
title('P1')
legend('1', '2', '3', '4', '5')
subplot(2,1,2);box on;hold on;
plot(s, dP1{1}, s, dP1{2}, s, dP1{3}, s, dP1{4}, s, dP1{5});
legend('1', '2', '3', '4', '5')
title('dP1')

figure;
subplot(2,1,1);box on;hold on;
plot(s, P2{1}, s, P2{2}, s, P2{3}, s, P2{4}, s, P2{5});
title('P2')
legend('1', '2', '3', '4', '5')
subplot(2,1,2);box on;hold on;
plot(s, dP2{1}, s, dP2{2}, s, dP2{3}, s, dP2{4}, s, dP2{5});
legend('1', '2', '3', '4', '5')
title('dP2')



