function [] = test_calc_b_LTX_3()
% function [] = test_calc_b_LTX_3()
%

separate_plots = 1;
% Make a grid
x_min = -1;
x_max = 1;
y_min = -1;
y_max = 1;
z_min = -1.2;
z_max = 1.2;

num_x_pts = 10;
num_y_pts = 10;
num_z_pts = 10;

[X, Y, Z] = ndgrid(linspace(x_min, x_max, num_x_pts), ...
    linspace(y_min, y_max, num_y_pts), ...
    linspace(z_min, z_max, num_z_pts));
Bx = zeros(size(X));
By = zeros(size(Y));
Bz = zeros(size(Z));

tic
% TF
coil_currents = [-1 0 0 0 0 0 0 0 0 0 0 0 0 0];
for ii = 1:num_x_pts
    for jj = 1:num_y_pts
        for kk = 1:num_z_pts
            [Bx(ii,jj,kk), By(ii,jj,kk), Bz(ii,jj,kk)] = ...
            calc_b_LTX(X(ii,jj,kk), Y(ii,jj,kk), Z(ii,jj,kk), coil_currents);            
        end
    end
end
figure;box on;hold on
quiver3(X, Y, Z, Bx, By, Bz, 'k');
title('TF');
legend('B');xlabel('m');ylabel('m');zlabel('m');
view(3);axis equal; axis equal;

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
if separate_plots 
    figure;box on;hold on
end
quiver3(X, Y, Z, Bx, By, Bz, 'm');
title('OH');
legend('B');xlabel('m');ylabel('m');zlabel('m');
view(3);axis equal; 

toc
% Red, upper
coil_currents = [0 0 1 0 0 0 0 0 0 0 0 0 0 0];
for ii = 1:num_x_pts
    for jj = 1:num_y_pts
        for kk = 1:num_z_pts
            [Bx(ii,jj,kk), By(ii,jj,kk), Bz(ii,jj,kk)] = ...
            calc_b_LTX(X(ii,jj,kk), Y(ii,jj,kk), Z(ii,jj,kk), coil_currents);            
        end
    end
end
if separate_plots 
    figure;box on;hold on
end
quiver3(X, Y, Z, Bx, By, Bz, 'r');
title('Red, Upper');
legend('B');xlabel('m');ylabel('m');zlabel('m');
view(3);axis equal; 

toc
% Red, lower
coil_currents = [0 0 0 1 0 0 0 0 0 0 0 0 0 0];
for ii = 1:num_x_pts
    for jj = 1:num_y_pts
        for kk = 1:num_z_pts
            [Bx(ii,jj,kk), By(ii,jj,kk), Bz(ii,jj,kk)] = ...
            calc_b_LTX(X(ii,jj,kk), Y(ii,jj,kk), Z(ii,jj,kk), coil_currents);            
        end
    end
end
if separate_plots 
    figure;box on;hold on
end
quiver3(X, Y, Z, Bx, By, Bz, 'r');
title('Red, Lower');
legend('B');xlabel('m');ylabel('m');zlabel('m');
view(3);axis equal; 

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
if separate_plots 
    figure;box on;hold on
end
quiver3(X, Y, Z, Bx, By, Bz, 'Color', [1 .5 0]);
title('Orange');
legend('B');xlabel('m');ylabel('m');zlabel('m');
view(3);axis equal; 

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
if separate_plots 
    figure;box on;hold on
end
quiver3(X, Y, Z, Bx, By, Bz, 'y');
title('Yellow');
legend('B');xlabel('m');ylabel('m');zlabel('m');
view(3);axis equal; 

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
if separate_plots 
    figure;box on;hold on
end
quiver3(X, Y, Z, Bx, By, Bz, 'g');
title('Green');
legend('B');xlabel('m');ylabel('m');zlabel('m');
view(3);axis equal; 

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
if separate_plots 
    figure;box on;hold on
end
quiver3(X, Y, Z, Bx, By, Bz, 'b');
title('Blue');
legend('B');xlabel('m');ylabel('m');zlabel('m');
view(3);axis equal; 

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
if separate_plots 
    figure;box on;hold on
end
quiver3(X, Y, Z, Bx, By, Bz, 'Color', [.5 .5 .5]);
legend('B');xlabel('m');ylabel('m');zlabel('m');
title('Internal');
view(3);axis equal; 

toc

if ~separate_plots 
    title('All fields')
end

