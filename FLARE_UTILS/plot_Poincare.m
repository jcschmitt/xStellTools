function plot_Poincare(filename)

[R, Z, T, PolFlux] = read_Poincare(filename);
figure;
plot(R, Z, 'r.');
axis equal
title(filename);



