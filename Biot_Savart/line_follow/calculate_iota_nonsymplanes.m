function iota = calculate_iota_nonsymplanes(coordsAxis, coordsSurf, ...
        phiEnd, phiIncInDegrees)

num_toroidal_planes = 360 / phiIncInDegrees;
    
% construct (R,Z) maps at each plane
%surf_R{num_surfs} = 0;
map_R{num_toroidal_planes} = 0;
map_Z{num_toroidal_planes} = 0;

num_pts_per_surface = length(coordsSurf(:,1));

disp('<----Generating R, Z maps');
% construct R data from x,y data
for jj = 1:num_toroidal_planes
    map_R{jj} = coordsSurf(jj:num_toroidal_planes:end, 1);
    map_Z{jj} = coordsSurf(jj:num_toroidal_planes:end, 2);
end

% find 'center' of innermost surface at each plane
axis_R(num_toroidal_planes) = 0;
axis_Z(num_toroidal_planes) = 0;
for jj = 1:num_toroidal_planes
    axis_R(jj) = coordsAxis(jj,1);
    axis_Z(jj) = coordsAxis(jj,2);    
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
        map_theta_1{jj} = atan2(map_Z{jj} - axis_Z(jj), map_R{jj} - axis_R(jj));
        map_theta_2{jj} = map_theta_1{jj};
        ind_fix = find(map_theta_2{jj} < 0);
        map_theta_2{jj}(ind_fix) = map_theta_2{jj}(ind_fix) + 2*pi;        
end

% loop over all surfaces. find the differential poloidal angle between each
% toroidal plane 
disp('<----Finding differential poloidal angle for points at each plane');

delta_theta = zeros(1, num_pts_per_surface-1);
delta_phi = (2*pi/num_toroidal_planes) * ones(1, num_pts_per_surface-1);
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
    delta_theta_1 = map_theta_1{ind_plane_2}(ind_point_2) - ...
        map_theta_1{ind_plane_1}(ind_point_1);
    delta_theta_2 = map_theta_2{ind_plane_2}(ind_point_2) - ...
        map_theta_2{ind_plane_1}(ind_point_1);
    if (abs(delta_theta_1) <= abs(delta_theta_2))
        delta_theta(jj) = delta_theta_1;
    else
        delta_theta(jj) = delta_theta_2;
    end
end

% are there any corrections or fixes necessary? (differences less than 0?)

% calculate total poloida angle
disp('<----Finding total differential poloidal and toroidal angles');

total_theta = sum(delta_theta);
total_phi = sum(delta_phi);
% calculate iota
%iota(ii) = total_theta(ii) / total_phi(ii);
iota = total_theta / total_phi;
%total_phi
%phiEnd

    

%keyboard

