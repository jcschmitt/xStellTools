s = linspace(0,1,51);


theI3 = 1;
I2{1} = sum_atan([0 1 .05 theI3 1], s);
I2{2} = sum_atan([0 1 1 theI3 1], s);
I2{3} = sum_atan([0 1 1.5 theI3 1], s);
I2{4} = sum_atan([0 1 2 theI3 1], s);
I2{5} = sum_atan([0 1 50 theI3 1], s);

I3{1} = sum_atan([0 1 1 .5 1], s);
I3{2} = sum_atan([0 1 1 1  1], s);
I3{3} = sum_atan([0 1 1 1.5 1], s);
I3{4} = sum_atan([0 1 1 2 1], s);
I3{5} = sum_atan([0 1 1 2.5 1], s);

IG1{1} = gauss_trunc([1 .01], s);
IG1{2} = gauss_trunc([1 .1], s);
IG1{3} = gauss_trunc([1 .2], s);
IG1{4} = gauss_trunc([1 .3], s);
IG1{5} = gauss_trunc([1 .4], s);

IG2{1} = gauss_trunc([1 .2], s);
IG2{2} = gauss_trunc([1 .5], s);
IG2{3} = gauss_trunc([1 1], s);
IG2{4} = gauss_trunc([1 2], s);
IG2{5} = gauss_trunc([1 5], s);

IDAT{1} = sum_atan([0 1 5 2 1 -.5 1 3 1 ], s) ;
IDAT{2} = sum_atan([0 1 5 2 1 -.5 2 3 1 ], s) ;
IDAT{3} = sum_atan([0 1 5 2 1 -.5 3 3 1 ], s) ;
IDAT{4} = sum_atan([0 1 5 2 1 -.5 4 3 1 ], s) ;
IDAT{5} = sum_atan([0 1 5 2 1 -.5 5 3 1 ], s) ;  % 5
IDAT{6} = sum_atan([0 1 5 2 1 -.5 1 1 1 ], s) ;
IDAT{7} = sum_atan([0 1 5 2 1 -.5 2 1 1 ], s) ;  % 7
IDAT{8} = sum_atan([0 1 5 2 1 -.5 3 1 1 ], s) ; %  8
IDAT{9} = sum_atan([0 1 5 2 1 -.5 4 1 1 ], s) ;
IDAT{10} = sum_atan([0 1 5 2 1 -.5 5 1 1 ], s) ; % 10
IDAT{11} = sum_atan([0 1 5 2 1 0 0 0 0 ], s) ;  % 11
IDAT{12} = sum_atan([0 1 5 2 1 -.5 2 .5 1 ], s) ;  % 12
IDAT{13} = sum_atan([0 1 5 2 1 -.5 2 2 1 ], s) ;  % 13
IDAT{14} = sum_atan([0 1 5 2 1 -.5 2 14 1 ], s) ;  % 14


for ii = 1:length(I2)
    J2{ii} = [nan diff(I2{ii})];
    J2{ii} = J2{ii} / max(J2{ii});
end
for ii = 1:length(I3)
    J3{ii} = [nan diff(I3{ii})];
    J3{ii} = J3{ii} / max(J3{ii});
end
for ii = 1:length(IG1)
    JG1{ii} = [nan diff(IG1{ii})];
    JG1{ii} = JG1{ii} / max(JG1{ii});
end
for ii = 1:length(IG2)
    JG2{ii} = [nan diff(IG2{ii})];
    JG2{ii} = JG2{ii} / max(JG2{ii});
end
for ii = 1:length(IDAT)
    JDAT{ii} = [nan diff(IDAT{ii})];
    JDAT{ii} = JDAT{ii} / max(abs(JDAT{ii}));
    % I want to normalize the current at the edge, not the density.
    IDAT2{ii} = IDAT{ii} / IDAT{ii}(end);
    JDAT2{ii} = [nan diff(IDAT2{ii})];
    
end



%=======================================
%=======================================
%=======================================

s=sqrt(s);

figure;
subplot(2,1,1);box on;hold on;
plot(s, I2{1}, s, I2{2}, s, I2{3}, s, I2{4}, s, I2{5});
title('I2')
legend('1', '2', '3', '4', '5')
subplot(2,1,2);box on;hold on;
plot(s, J2{1}, s, J2{2}, s, J2{3}, s, J2{4}, s, J2{5});
legend('1', '2', '3', '4', '5')
title('J2')

figure;
subplot(2,1,1);box on;hold on;
plot(s, I3{1}, s, I3{2}, s, I3{3}, s, I3{4}, s, I3{5});
title('I3')
legend('1', '2', '3', '4', '5')
subplot(2,1,2);box on;hold on;
plot(s, J3{1}, s, J3{2}, s, J3{3}, s, J3{4}, s, J3{5});
legend('1', '2', '3', '4', '5')
title('J3')

figure;
subplot(2,1,1);box on;hold on;
plot(s, IG1{1}, s, IG1{2}, s, IG1{3}, s, IG1{4}, s, IG1{5});
title('IG1')
legend('1', '2', '3', '4', '5')
subplot(2,1,2);box on;hold on;
plot(s, JG1{1}, s, JG1{2}, s, JG1{3}, s, JG1{4}, s, JG1{5});
legend('1', '2', '3', '4', '5')
title('JG1')

figure;
subplot(2,1,1);box on;hold on;
plot(s, IG2{1}, s, IG2{2}, s, IG2{3}, s, IG2{4}, s, IG2{5});
title('IG2')
legend('1', '2', '3', '4', '5')
subplot(2,1,2);box on;hold on;
plot(s, JG2{1}, s, JG2{2}, s, JG2{3}, s, JG2{4}, s, JG2{5});
legend('1', '2', '3', '4', '5')
title('JG2')

figure;
subplot(2,1,1);box on;hold on;
plot(s, IDAT{1}, s, IDAT{2}, s, IDAT{3}, s, IDAT{4}, s, IDAT{5}, ...
    s, IDAT{6}, s, IDAT{7}, s, IDAT{8}, s, IDAT{9}, s, IDAT{10}, ...
    s, IDAT{11}, s, IDAT{12}, s, IDAT{13}, s, IDAT{14});
title('IDAT')
legend('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14')
subplot(2,1,2);box on;hold on;
plot(s, JDAT{1}, s, JDAT{2}, s, JDAT{3}, s, JDAT{4}, s, JDAT{5}, ...
    s, JDAT{6}, s, JDAT{7}, s, JDAT{8}, s, JDAT{9}, s, JDAT{10}, s, JDAT{11}, ...
    s, JDAT{12}, s, JDAT{13}, s, JDAT{14});
legend('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14')
title('JDAT')

figure;
subplot(2,1,1);box on;hold on;
plot(s, IDAT2{1}, s, IDAT2{2}, s, IDAT2{3}, s, IDAT2{4}, s, IDAT2{5}, ...
    s, IDAT2{6}, s, IDAT2{7}, s, IDAT2{8}, s, IDAT2{9}, s, IDAT2{10}, ...
    s, IDAT2{11}, s, IDAT2{2}, s, IDAT2{13}, s, IDAT2{14});
title('IDAT2')

legend('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14')
subplot(2,1,2);box on;hold on;
plot(s, JDAT2{1}, s, JDAT2{2}, s, JDAT2{3}, s, JDAT2{4}, s, JDAT2{5}, ...
    s, JDAT2{6}, s, JDAT2{7}, s, JDAT2{8}, s, JDAT2{9}, s, JDAT2{10}, ...
    s, JDAT2{11}, s, JDAT2{12}, s, JDAT2{13}, s, JDAT2{14});
legend('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14')
title('JDAT2')


figure;
subplot(2,1,1);box on;hold on;
plot(... %s, IDAT2{1}, s, IDAT2{2}, s, IDAT2{3}, s, IDAT2{4}, s, IDAT2{5}, ...
    ... %s, IDAT2{6}, s, IDAT2{7}, s, IDAT2{8}, s, IDAT2{9}, s, IDAT2{10}, ...
    s, IDAT2{7}, s, IDAT2{2}, s, IDAT2{13}, s, IDAT2{14});
title('IDAT2')
legend('7', '12', '13', '14')

subplot(2,1,2);box on;hold on;
plot(... %s, JDAT{1}, s, JDAT{2}, s, JDAT{3}, s, JDAT{4}, s, JDAT{5}, ...
    ... %s, JDAT{6}, s, JDAT{7}, s, JDAT{8}, s, JDAT{9}, s, JDAT{10}, s, JDAT{11}, ...
    s, JDAT2{7}, s, JDAT{12}, s, JDAT{13}, s, JDAT{14});
legend('7', '12', '13', '14')
title('JDAT')

