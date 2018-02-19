function [] = test_calc_b_LTX_maxB()

separate_plots = 1;
% Make a grid
x_min = 0.14;
x_max = 0.68;
y_min = 0;
y_max = 0;
z_min = -0.4;
z_max = 0.4;

num_x_pts = 16;
num_y_pts = 1;
num_z_pts = 16;


[X, Y, Z] = ndgrid(linspace(x_min, x_max, num_x_pts), ...
    linspace(y_min, y_max, num_y_pts), ...
    linspace(z_min, z_max, num_z_pts));
Bx = zeros(size(X));
By = zeros(size(Y));
Bz = zeros(size(Z));

tic
% TF
coil_currents = [2765 -40000 5000 5000 5000 5000 70000 70000 4000 4000 10000 10000 5000 5000];
for ii = 1:num_x_pts
    for jj = 1:num_y_pts
        for kk = 1:num_z_pts
            [Bx(ii,jj,kk), By(ii,jj,kk), Bz(ii,jj,kk)] = ...
            calc_b_LTX(X(ii,jj,kk), Y(ii,jj,kk), Z(ii,jj,kk), coil_currents); 
        end
    end
end
% B_R_mag = sqrt(Bx.^2 + Bz.^2);

figure;subplot(2,3,1);box on;hold on
quiver3(X, Y, Z, Bx, By, Bz, 'k');
title('B vectors');
axis tight;view(0,0);axis equal; 

subplot(2,3,2);box on;hold on
contourf(squeeze(X), squeeze(Z), squeeze(sqrt(Bz.^2 + By.^2 + Bx.^2)))
axis equal; axis tight
colorbar;
title('|B|');

% subplot(2,3,3);box on;hold on
% contourf(squeeze(X), squeeze(Z), squeeze(B_R_mag))
% axis equal; axis tight
% colorbar;
% title('|BR|');

subplot(2,3,4);box on;hold on
contourf(squeeze(X), squeeze(Z), squeeze(Bx))
axis equal; axis tight
colorbar;
title('Bx');

subplot(2,3,5);box on;hold on
contourf(squeeze(X), squeeze(Z), squeeze(By))
axis equal; axis tight
colorbar;
title('By = Btoroidal');

subplot(2,3,6);box on;hold on
contourf(squeeze(X), squeeze(Z), squeeze(Bz))
axis equal; axis tight
colorbar;
title('Bz')

