function iota = calculate_iota_nonsymplanes(coordsAxis, coords, ...
        phi(end), phiIncInDegrees);
%[iota] = calculate_iota_nonstelsym(rAxis, coords, phiIncInDegrees);
        
% construct (R,Z) maps at each plane
surf_R{num_surfs} = 0;
map_R{num_surfs, num_toroidal_planes} = 0;
map_Z{num_surfs, num_toroidal_planes} = 0;

num_pts_per_surface = length(surf_x{1});

disp('<----Generating R, Z maps');
for ii = 1:num_surfs
    % construct R data from x,y data
    surf_R{ii} = sqrt(surf_x{ii}.^2 + surf_y{ii}.^2);
    for jj = 1:num_toroidal_planes
        map_R{ii, jj} = surf_R{ii}(jj:num_toroidal_planes:end);
        map_Z{ii, jj} = surf_z{ii}(jj:num_toroidal_planes:end);
    end
end


% find 'center' of innermost surface at each plane
disp('<----Finding axis of innermost surface map');
axis_R(num_toroidal_planes) = 0;
axis_Z(num_toroidal_planes) = 0;
for jj = 1:num_toroidal_planes
    axis_R(jj) = mean(map_R{1, jj});
    axis_Z(jj) = mean(map_Z{1, jj});
    
end

% use center of inner most surface as raxis and zaxis at each plane -
% calculate effective 'poloidal' angle for each point at each plane
% map_theta = theta 'map' at local toroidal angle
% map_theta_2 = alternative thata 'map'
% the two theta maps are used in conjunction to determine the 'real' delta
% theta between two toroidal slics, allowing for the possiblity that the
% theta angles may cross the -pi/pi boundary (which is accounted for by
% using the 2nd map to find the 'true' delta theta.

disp('<----Finding effecting poloidal angle for points at each plane');
for jj = 1:num_toroidal_planes
    axis_R(jj) = mean(map_R{1, jj});
    axis_Z(jj) = mean(map_Z{1, jj});
    for ii = 1:num_surfs
        map_theta_1{ii,jj} = atan2(map_Z{ii, jj} - axis_Z(jj), map_R{ii,jj} - axis_R(jj));
        map_theta_2{ii,jj} = map_theta_1{ii,jj};
        ind_fix = find(map_theta_2{ii,jj} < 0);
        map_theta_2{ii,jj}(ind_fix) = map_theta_2{ii,jj}(ind_fix) + 2*pi;        
    end
end

% loop over all surfaces. find the differential poloidal angle between each
% toroidal plane 
disp('<----Finding differential poloidal angle for points at each plane');
delta_theta{num_surfs} = 0;
delta_phi{num_surfs} = 0;

for ii = 1:num_surfs
    delta_theta{ii} = zeros(1, num_pts_per_surface-1);
    delta_phi{ii} = (2*pi/num_toroidal_planes) * ones(1, num_pts_per_surface-1);
    for jj = 1:(num_pts_per_surface-1)
        ind_plane_1 = mod(jj, num_toroidal_planes);
        ind_point_1 = floor(jj/num_toroidal_planes)+1;
        ind_plane_2 = mod(jj+1, num_toroidal_planes);
        ind_point_2 = floor((jj+1)/num_toroidal_planes)+1;
        if (ind_plane_1 == 0)
            ind_plane_1 = num_toroidal_planes;
            ind_point_1 = ind_point_1 - 1;
        elseif (ind_plane_2 == 0)
            ind_plane_2 = num_toroidal_planes;
            ind_point_2 = ind_point_2 - 1;
        end
        delta_theta_1 = map_theta_1{ii, ind_plane_2}(ind_point_2) - ...
            map_theta_1{ii, ind_plane_1}(ind_point_1);
        delta_theta_2 = map_theta_2{ii, ind_plane_2}(ind_point_2) - ...
            map_theta_2{ii, ind_plane_1}(ind_point_1);
        if (abs(delta_theta_1) <= abs(delta_theta_2))
            delta_theta{ii}(jj) = delta_theta_1;
        else
            delta_theta{ii}(jj) = delta_theta_2;
        end
    end    
end

% are there any corrections or fixes necessary? (differences less than 0?)

% calculate total poloida angle
disp('<----Finding total differential poloidal and toroidal angles');
total_theta(num_surfs) = 0;
total_phi(num_surfs) = 0;
iota(num_surfs) = 0;

for ii = 1:num_surfs
    total_theta(ii) = sum(delta_theta{ii});
    total_phi(ii) = sum(delta_phi{ii});
    % calculate iota
    %iota(ii) = total_theta(ii) / total_phi(ii);
    iota(ii) = total_theta(ii) / surf_phi(end);
end


%keyboard

