function [coils, MVALS, NVALS, RMNC, ZMNS, RMNS, ZMNC] = ...
    convert_coils_to_winding_surface(coil_filename, M_max, N_max)

% DEBUG = 0
DEBUG = 1;
%coil_filename = 'coils.wistell_a_0007';
%M_max = 16; N_max = 24;

coils = read_coils(coil_filename);
MVALS = [];
NVALS = [];
RMNC = [];
RMNS = [];
ZMNC = [];
ZMNS = [];

vert_x = coils.vert(1,:);
vert_y = coils.vert(2,:);
vert_z = coils.vert(3,:);
num_vertices = length(vert_x);
nfp = coils.periods;

if DEBUG
    fh1 = figure;
    plot3(vert_x, vert_y, vert_z, '.');
    hold on;
end

% Got a bunch of points in space.  Assume that they are supposed to
% represent some sort of simple toroidal winding surface for a stellarator
% with stellarator-symmetry.  Each point is given by the vertices in
% coils.vert(1:3,:)

% Going to use lsqnonlin fto find the solution.
% The function will be an anonymous function based on 
% The inital guess
for mm = 0:M_max
    for nn = (0:-1:-N_max)*nfp
        MVALS = [MVALS mm];
        NVALS = [NVALS nn];
    end
end
num_modes = length(MVALS);

% inital guess and bounds
% Input x0 contains all of the RMNC modes followed by all of the ZMNS
% modes.  They will be extracted in the function below.

x0 = 0 * [MVALS MVALS];
x0(1) = 2; % The R00 cos component is the 1st element
%x0((num_modes+1) = 0.4;  % The ???
lb = -3*ones(size(x0));
ub =  3*ones(size(x0));
%options = optimoptions('lsqnonlin', 'MaxIterations', 4, 'UseParallel', true); 
options = optimoptions('lsqnonlin', 'UseParallel', true); 

fun = @(x)tor_surf_diff_fun(x, MVALS, NVALS, vert_x, vert_y, vert_z);
x = lsqnonlin(fun,x0,lb,ub,options);

[x_surf, y_surf, z_surf] = make_tor_surf(x, MVALS, NVALS);

if DEBUG
    figure(fh1);
    surf(x_surf, y_surf, z_surf);
    hold on;
end
RMNC = x(1:num_modes);
ZMNS = x((num_modes+1):end);
% if DEBUG
%     figure(fh1);
%     plot3(x_reb, y_reb, z_reb, 'ro')
%     plot3(x_surf, y_surf, z_surf, 'bx')
% end
% 
% [MVALS' NVALS' RMNC' ZMNS']
% 
keyboard



function deltas = tor_surf_diff_fun(x_in, MVALS, NVALS, vert_x, ...
    vert_y, vert_z)

[x_surf, y_surf, z_surf] = make_tor_surf(x_in, MVALS, NVALS);
% Find detlas
deltas = 0 * length(vert_x);
%disp('<----Deltas');
for ii = 1:length(vert_x)
    delta_mesh = sqrt((x_surf - vert_x(ii)).^2 + ...
                      (y_surf - vert_y(ii)).^2 + ...
                      (z_surf - vert_z(ii)).^2);
    deltas(ii) = min(min(delta_mesh));
end       
%disp(['<----sqrt(sum(delta^2))) = ' num2str(sqrt(sum(deltas.^2)))]);


function [x_surf, y_surf, z_surf] = make_tor_surf(x_in, MVALS, NVALS)
tor_grid_pts = 360;
pol_grid_ts = 64;
%tor_grid_pts = 10*2*4;
%pol_grid_ts = 10*2;

num_modes = length(MVALS);

% Assume RMNC modes are in the first half of x_in, ZMNS modes are in 2nd
% half of x_in  (this is the stellarator symmetric fit)
RMNC = x_in(1:num_modes);
ZMNS = x_in((num_modes+1):(2*num_modes));

% Build the winding surface
%disp('<----Building');

[zeta_surf, theta_surf] = meshgrid(linspace(0,2*pi,tor_grid_pts), ...
    linspace(0,2*pi,pol_grid_ts));
% initialize to 0
r_surf = 0 * zeta_surf;
z_surf = 0 * zeta_surf;

for jj = 1:num_modes
    mm = MVALS(jj);
    nn = NVALS(jj);
    mtheta_plus_nzeta = mm * theta_surf + nn * zeta_surf;
    r_surf  = r_surf + RMNC(jj) * cos(mtheta_plus_nzeta);% + ...
        %RMNS(jj) * sin(mtheta_plus_nzeta);
    z_surf = z_surf  + ZMNS(jj) * sin(mtheta_plus_nzeta);% + ...
        %ZMNC(jj) * cos(mtheta_plus_nzeta);
end

x_surf = r_surf .* cos(zeta_surf);
y_surf = r_surf .* sin(zeta_surf);

% Surface is generated


