function [] = test_calc_b_LTX_yellow()
% function [] = test_calc_b_LTX_yellow()
%


% Make a grid
x_min = 0;
x_max = .7;
y_min = 0;
y_max = 0;
z_min = 0;
z_max = 0.4;

num_x_pts = 100;
num_y_pts = 1;
num_z_pts = 5;

color_list = 'brgkc';

xs = linspace(x_min, x_max, num_x_pts);
ys = linspace(y_min, y_max, num_y_pts);
zs = linspace(z_min, z_max, num_z_pts);

% Yellow
coil_currents = [0 0 0 0 0 0 1 1 0 0 0 0 0 0];
for ii = 1:num_x_pts
    for jj = 1:num_y_pts
        for kk = 1:num_z_pts
            [Bx(ii,jj,kk), By(ii,jj,kk), Bz(ii,jj,kk)] = ...
            calc_b_LTX(xs(ii), ys(jj), zs(kk), coil_currents);            
        end
    end
end

figure;box on;hold on

for ii = 1:num_z_pts
    plot(xs, 1e4*Bz(:,1,ii), color_list(ii));
end
legend('Z = 0', 'Z = 0.1', 'Z = 0.2', 'Z = 0.3', 'Z = 0.4');
xlabel('x (or R) in meters')
ylabel('Gauss')
title('Vertical field (Gauss) due to 1 Amp in Yellow coils');

