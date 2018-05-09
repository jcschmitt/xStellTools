function [fig_handle, R, Z, theta_wrt_axis] = ...
    tangent_map(IOTA_DIRECTORY_PATH, surfaceIndices, indexInc)

if (nargin < 3)
    indexInc = 1;
end

if surfaceIndices(1) ~= 1
    disp('<---- Prepending 1 to the surface index list');
    surfaceIndices = [1 surfaceIndices]    
end

[fig_handle, R, Z] = plotPuncturePlot(IOTA_DIRECTORY_PATH, surfaceIndices, indexInc);

% Build tangent map information

% Will need to calculate the poloidal angle.  
% Estimate magnetic axis
R_axis = mean(R{1});
Z_axis = mean(Z{1});


for ii = 1:length(surfaceIndices)
    % do the actual theta calculation (arctan)
    theta_wrt_axis{ii} = atan2(Z{ii} - Z_axis, R{ii} - R_axis);    
end

fig_handle = figure;
for ii = 1:length(surfaceIndices)
    plot(theta_wrt_axis{ii}, 'o:');
    hold on;
end
legend(num2str(surfaceIndices'))

close(fig_handle);

fig_handle2 = figure;
sample_freq = 4;
for ii = 1:length(surfaceIndices)
    mydtheta{ii} = diff(theta_wrt_axis{ii});
    L = length(mydtheta{ii});
    L_lf(ii) = L;
    n = 2^nextpow2(L);
    my_fft = fft(mydtheta{ii}, n);
    f{ii} = sample_freq * (0:(n/2))/n;
    P{ii} = abs(my_fft / n)';
    P_plot{ii} = P{ii}(1:n/2+1);
    
    plot3(surfaceIndices(ii)*ones(size(f{ii})), f{ii}, P_plot{ii}, 'bo-');
    hold on;
end
xlabel('Index');
ylabel('Eigenvalue / Full Transit');
zlabel('Power');

[f_grid, index_grid] = meshgrid(f{1}, surfaceIndices(2:end));
power_surf_data = zeros(size(f_grid));
L = length(mydtheta{1});
n = 2^nextpow2(L);
for ii = 2:length(surfaceIndices)
    try
    power_surf_data(ii-1, :) = interp1(f{ii-1}, P_plot{ii-1}, f_grid(ii-1,:));
    catch
        %keyboard
        
    end
end

figure
plot(surfaceIndices, L_lf, 'o')
xlabel('Index');
ylabel('Map Length');

figure_handle3 = figure;
surf(index_grid, f_grid, power_surf_data, 'EdgeColor', 'None', 'FaceColor', ...
    'interp');

xlabel('Index');
ylabel('Poloidal Transform');
zlabel('Power');
view(2);