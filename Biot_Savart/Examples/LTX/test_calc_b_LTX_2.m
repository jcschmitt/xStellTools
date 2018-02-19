function [] = test_calc_b_LTX_2()
% function [] = test_calc_b_LTX_2()
%

% Make a grid
x_min = 0;
x_max = 0.8;
y_min = 0;
y_max = 0;
z_min = -0.5;
z_max = 0.5;

num_x_pts = 20;
num_y_pts = 1;
num_z_pts = 20;

num_contours = 20;

% % Make a grid
% x_min = -1;
% x_max = 1;
% y_min = -1;
% y_max = 1;
% z_min = -1.2;
% z_max = 1.2;
% 
% num_x_pts = 10;
% num_y_pts = 10;
% num_z_pts = 10;

[X, Y, Z] = ndgrid(linspace(x_min, x_max, num_x_pts), ...
    linspace(y_min, y_max, num_y_pts), ...
    linspace(z_min, z_max, num_z_pts));
Bx = zeros(size(X));
By = zeros(size(Y));
Bz = zeros(size(Z));

tic
% TF
coil_currents = [1 0 0 0 0 0 0 0 0 0 0 0 0 0];
for ii = 1:num_x_pts
    for jj = 1:num_y_pts
        for kk = 1:num_z_pts
            [Bx(ii,jj,kk), By(ii,jj,kk), Bz(ii,jj,kk)] = ...
            calc_b_LTX(X(ii,jj,kk), Y(ii,jj,kk), Z(ii,jj,kk), coil_currents);            
        end
    end
end
figure;subplot(2,2,1);box on;hold on
quiver3(X, Y, Z, Bx, By, Bz, 'k');
title('TF');
axis tight;view(0,0);axis equal; 

subplot(2,2,2);box on;hold on
contourf(squeeze(X), squeeze(Z), squeeze(sqrt(Bz.^2 + By.^2 + Bx.^2)))
colorbar;
axis equal; axis tight
title('|B|');

subplot(2,2,3);box on;hold on
contourf(squeeze(X), squeeze(Z), squeeze(Bz), num_contours); 
colorbar;
axis equal; axis tight
title('B_z')

subplot(2,2,4);box on;hold on
contourf(squeeze(X), squeeze(Z), squeeze(Bx))
colorbar;
axis equal; axis tight
title('B_x');

toc
% OH
coil_currents = [0 1 0 0 0 0 0 0 0 0 0 0 0 0];
for ii = 1:num_x_pts
    for jj = 1:num_y_pts
        for kk = 1:num_z_pts
            [Bx(ii,jj,kk), By(ii,jj,kk), Bz(ii,jj,kk)] = ...
            calc_b_LTX(X(ii,jj,kk), Y(ii,jj,kk), Z(ii,jj,kk), coil_currents);            
        end
    end
end

figure;subplot(2,2,1);box on;hold on
quiver3(X, Y, Z, Bx, By, Bz, 'm');
title('OH');
axis tight;view(0,0);axis equal; 

subplot(2,2,2);box on;hold on
contourf(squeeze(X), squeeze(Z), squeeze(sqrt(Bz.^2 + By.^2 + Bx.^2)))
colorbar;
axis equal; axis tight
title('|B|');

subplot(2,2,3);box on;hold on
contourf(squeeze(X), squeeze(Z), squeeze(Bz), num_contours);
colorbar;
axis equal; axis tight
title('B_z')

subplot(2,2,4);box on;hold on
contourf(squeeze(X), squeeze(Z), squeeze(Bx))
colorbar;
axis equal; axis tight
title('B_x');


toc
% Red, both
coil_currents = [0 0 1 1 0 0 0 0 0 0 0 0 0 0];
for ii = 1:num_x_pts
    for jj = 1:num_y_pts
        for kk = 1:num_z_pts
            [Bx(ii,jj,kk), By(ii,jj,kk), Bz(ii,jj,kk)] = ...
            calc_b_LTX(X(ii,jj,kk), Y(ii,jj,kk), Z(ii,jj,kk), coil_currents);            
        end
    end
end

figure;subplot(2,2,1);box on;hold on
quiver3(X, Y, Z, Bx, By, Bz, 'r');
title('Red, Both');
axis tight;view(0,0);axis equal; 

subplot(2,2,2);box on;hold on
contourf(squeeze(X), squeeze(Z), squeeze(sqrt(Bz.^2 + By.^2 + Bx.^2)))
axis equal; axis tight
title('|B|');

subplot(2,2,3);box on;hold on
contourf(squeeze(X), squeeze(Z), squeeze(Bz), num_contours);  colorbar; 
axis equal; axis tight
title('B_z')

subplot(2,2,4);box on;hold on
contourf(squeeze(X), squeeze(Z), squeeze(Bx))
axis equal; axis tight
title('B_x');

toc
% Orange
coil_currents = [0 0 0 0 1 1 0 0 0 0 0 0 0 0];
for ii = 1:num_x_pts
    for jj = 1:num_y_pts
        for kk = 1:num_z_pts
            [Bx(ii,jj,kk), By(ii,jj,kk), Bz(ii,jj,kk)] = ...
            calc_b_LTX(X(ii,jj,kk), Y(ii,jj,kk), Z(ii,jj,kk), coil_currents);            
        end
    end
end

figure;subplot(2,2,1);box on;hold on
quiver3(X, Y, Z, Bx, By, Bz, 'Color', [1 .5 0]);
title('Orange');
axis tight;view(0,0);axis equal; 

subplot(2,2,2);box on;hold on
contourf(squeeze(X), squeeze(Z), squeeze(sqrt(Bz.^2 + By.^2 + Bx.^2)))
colorbar;
axis equal; axis tight
title('|B|');

subplot(2,2,3);box on;hold on
contourf(squeeze(X), squeeze(Z), squeeze(Bz), num_contours); 
colorbar; 
axis equal; axis tight
title('B_z')

subplot(2,2,4);box on;hold on
contourf(squeeze(X), squeeze(Z), squeeze(Bx))
colorbar;
axis equal; axis tight
title('B_x');

toc
% Yellow
coil_currents = [0 0 0 0 0 0 1 1 0 0 0 0 0 0];
for ii = 1:num_x_pts
    for jj = 1:num_y_pts
        for kk = 1:num_z_pts
            [Bx(ii,jj,kk), By(ii,jj,kk), Bz(ii,jj,kk)] = ...
            calc_b_LTX(X(ii,jj,kk), Y(ii,jj,kk), Z(ii,jj,kk), coil_currents);            
        end
    end
end

figure;subplot(2,2,1);box on;hold on
quiver3(X, Y, Z, Bx, By, Bz, 'y');
title('Yellow');
axis tight;view(0,0);axis equal; 

subplot(2,2,2);box on;hold on
contourf(squeeze(X), squeeze(Z), squeeze(sqrt(Bz.^2 + By.^2 + Bx.^2)))
colorbar;
axis equal; axis tight
title('|B|');

subplot(2,2,3);box on;hold on
contourf(squeeze(X), squeeze(Z), squeeze(Bz), num_contours);  
colorbar; 
axis equal; axis tight
title('B_z')

subplot(2,2,4);box on;hold on
contourf(squeeze(X), squeeze(Z), squeeze(Bx))
colorbar;
axis equal; axis tight
title('B_x');

toc
% Green
coil_currents = [0 0 0 0 0 0 0 0 1 1 0 0 0 0];
for ii = 1:num_x_pts
    for jj = 1:num_y_pts
        for kk = 1:num_z_pts
            [Bx(ii,jj,kk), By(ii,jj,kk), Bz(ii,jj,kk)] = ...
            calc_b_LTX(X(ii,jj,kk), Y(ii,jj,kk), Z(ii,jj,kk), coil_currents);            
        end
    end
end

figure;subplot(2,2,1);box on;hold on
quiver3(X, Y, Z, Bx, By, Bz, 'g');
title('Green');
axis tight;view(0,0);axis equal; 
title('|B|');

subplot(2,2,2);box on;hold on
contourf(squeeze(X), squeeze(Z), squeeze(sqrt(Bz.^2 + By.^2 + Bx.^2)))
colorbar;
axis equal; axis tight
title('B_z')

subplot(2,2,3);box on;hold on
contourf(squeeze(X), squeeze(Z), squeeze(Bz), num_contours); 
colorbar; 
axis equal; axis tight

subplot(2,2,4);box on;hold on
contourf(squeeze(X), squeeze(Z), squeeze(Bx))
colorbar;
axis equal; axis tight
title('B_x');

toc
% Blue
coil_currents = [0 0 0 0 0 0 0 0 0 0 1 1 0 0];
for ii = 1:num_x_pts
    for jj = 1:num_y_pts
        for kk = 1:num_z_pts
            [Bx(ii,jj,kk), By(ii,jj,kk), Bz(ii,jj,kk)] = ...
            calc_b_LTX(X(ii,jj,kk), Y(ii,jj,kk), Z(ii,jj,kk), coil_currents);            
        end
    end
end

figure;subplot(2,2,1);box on;hold on
quiver3(X, Y, Z, Bx, By, Bz, 'b');
title('Blue, B vector');
axis tight;view(0,0);axis equal; 

subplot(2,2,2);box on;hold on
contourf(squeeze(X), squeeze(Z), squeeze(sqrt(Bz.^2 + By.^2 + Bx.^2)))
colorbar;
axis equal; axis tight
title('|B|');

subplot(2,2,3);box on;hold on
contourf(squeeze(X), squeeze(Z), squeeze(Bz), num_contours); 
colorbar;
axis equal; axis tight
title('B_z');

subplot(2,2,4);box on;hold on
contourf(squeeze(X), squeeze(Z), squeeze(Bx))
colorbar;
axis equal; axis tight
title('B_x');

toc
% Internal
coil_currents = [0 0 0 0 0 0 0 0 0 0 0 0 1 1];
for ii = 1:num_x_pts
    for jj = 1:num_y_pts
        for kk = 1:num_z_pts
            [Bx(ii,jj,kk), By(ii,jj,kk), Bz(ii,jj,kk)] = ...
            calc_b_LTX(X(ii,jj,kk), Y(ii,jj,kk), Z(ii,jj,kk), coil_currents);            
        end
    end
end

figure;subplot(2,2,1);box on;hold on
quiver3(X, Y, Z, Bx, By, Bz, 'Color', [.5 .5 .5]);
title('Internal, B vector');
axis tight;view(0,0);axis equal; 

subplot(2,2,2);box on;hold on
contourf(squeeze(X), squeeze(Z), squeeze(sqrt(Bz.^2 + By.^2 + Bx.^2)))
colorbar;
axis equal; axis tight
title('|B|');

subplot(2,2,3);box on;hold on
contourf(squeeze(X), squeeze(Z), squeeze(Bz), num_contours);  
colorbar; 
axis equal; axis tight
title('B_z')

subplot(2,2,4);box on;hold on
contourf(squeeze(X), squeeze(Z), squeeze(Bx))
colorbar;
axis equal; axis tight
title('B_x');

toc
