function [] = test_calc_b_LTX_OH()
% function [] = test_calc_b_LTX_OH()
%

% Make a grid
x_min = 0.1;
x_max = 0.7;
y_min = 0;
y_max = 0;
z_min = -0.45;
z_max = 0.45;

% z_min = -0.2;
% z_max = 0.2;
quiver_scale = 5;

% % Make a grid
% x_min = 0.1;
% x_max = 0.7;
% y_min = 0;
% y_max = 0;
% z_min = -0.1;
% z_max = 0.1;

num_x_pts = 21;
num_y_pts = 1;
num_z_pts = 21;

[X, Y, Z] = ndgrid(linspace(x_min, x_max, num_x_pts), ...
    linspace(y_min, y_max, num_y_pts), ...
    linspace(z_min, z_max, num_z_pts));
Bx = zeros(size(X));
By = zeros(size(Y));
Bz = zeros(size(Z));
Bxm = zeros(size(X));
Bym = zeros(size(Y));
Bzm = zeros(size(Z));

tic
% OH
coil_currents = 1000*[0 1 0 0 0 0 0 0 0 0 0 0 0 0];
for ii = 1:num_x_pts
    for jj = 1:num_y_pts
        for kk = 1:num_z_pts
            [Bx(ii,jj,kk), By(ii,jj,kk), Bz(ii,jj,kk)] = ...
            calc_b_LTX(X(ii,jj,kk), Y(ii,jj,kk), Z(ii,jj,kk), coil_currents);            
            [Bxm(ii,jj,kk), Bym(ii,jj,kk), Bzm(ii,jj,kk)] = ...
            calc_b_LTX(X(ii,jj,kk), Y(ii,jj,kk), Z(ii,jj,kk), -coil_currents);            
        end
    end
end


figure;box on;hold on
quiver3(X, Y, Z, Bx, By, Bz, quiver_scale, 'm');
title('OH');
axis tight;view(0,0);axis equal; 

num_contours = 45;

figure;
subplot(2,2,1);box on;hold on
[C,h]=contourf(squeeze(X), squeeze(Z), squeeze(sqrt(Bz.^2 + By.^2 + Bx.^2)))
clabel(C,h);
axis tight;axis equal; 
title('OH - |B|');
colorbar

subplot(2,2,2);box on;hold on
[C,h]=contourf(squeeze(X), squeeze(Z), squeeze(By),num_contours);
clabel(C,h);
axis equal; axis tight
title('B_y');
colorbar

subplot(2,2,3);box on;hold on
[C,h]=contourf(squeeze(X), squeeze(Z), squeeze(Bz),num_contours);
clabel(C,h);
axis equal; axis tight
title('B_z')
colorbar

subplot(2,2,4);box on;hold on
[C,h]=contourf(squeeze(X), squeeze(Z), squeeze(Bx),num_contours);
clabel(C,h);
axis equal; axis tight
title('B_x');
colorbar

figure;
subplot(2,2,1);box on;hold on
[C,h]=contourf(squeeze(X), squeeze(Z), squeeze(sqrt(Bzm.^2 + Bym.^2 + Bxm.^2)))
clabel(C,h);
axis tight;axis equal; 
title('OH - |B|');
colorbar

subplot(2,2,2);box on;hold on
[C,h]=contourf(squeeze(X), squeeze(Z), squeeze(Bym),num_contours);
clabel(C,h);
axis equal; axis tight
title('B_y');
colorbar

subplot(2,2,3);box on;hold on
[C,h]=contourf(squeeze(X), squeeze(Z), squeeze(Bzm),num_contours);
clabel(C,h);
axis equal; axis tight
title('B_z')
colorbar

subplot(2,2,4);box on;hold on
[C,h]=contourf(squeeze(X), squeeze(Z), squeeze(Bxm),num_contours);
clabel(C,h);
axis equal; axis tight
title('B_x');
colorbar
