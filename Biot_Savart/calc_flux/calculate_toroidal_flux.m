function [flux] = calc_toroidal_flux(coilsetID, current, R, Z, MAKE_PLOT)
% Returns the flux through a surface defined by [R, Z].  Assumes the
% surface is at the boxport of HSX (at the 'D' shaped toroidal plane).  

% If MAKE_PLOT was not specified, set to 0.
if nargin < 5
    MAKE_PLOT = 0; % Default is 'no'
end

% 1e-9 seems to be a good value for convergence. Should be checked for each
% configuration.
% 1e-9 works in most cases
%TOL_FLUX = 1e-4;
%TOL_FLUX = 1e-9;
%TOL_FLUX = 1e-14;

% trying to improve the EIM edge island situation
%TOL_FLUX = 1e-4;
% A 40th order polynomial usually works, unless the surface is very
% 'corrugated'.
POLYFIT_ORDER = 40;


% Step 1: Separate the points into two regions of integration.
%       1a. Identify the apex of the surface, in terms of z.  
%       1b. Identify the points to the each side of the apex (inboard and
%           outboard w.r.t. the apex). 
%       1c. Each set of points is fit with a high-order polynomial, and
%       this polynomial is used to define the flux boundary during surface
%       integration.
% Old way: Rely on the assumption that the surfaces are NOT concave/convex
% New way: Permit concavity/conexity.  Sort point according to some sort of
% poloidal angle and separate the 'inboard' and 'outboard' regions with
% this additional knowledge.
% New(er) way: Sort points manuallly.  Start with 1st index.  Search
% counter-clockwise for the next nearest point and append it to the
% 'sorted' list. Remove that point from future searches

% Define R_center as the mean of the radius of the 4 points nearest

% Old way
R_center = mean(R);
Z_center = mean(Z);

poloidal_angle = atan2( (Z - Z_center), (R - R_center) );
poloidal_angle_2 = atan2( (Z - Z_center), (R) );

% Old way ends here.

% New way starts here
% if the min and max of poloidal angle aren't close to -/+ pi, then it is
% likely that the 'center' found above is not enclosed by the locus of
% points.
max_pa = max(poloidal_angle);
min_pa = min(poloidal_angle);
% tolerance of +/1 degree
angle_tolerance = pi/180;
if ( ( (pi - max_pa) > angle_tolerance) && ...
     ( (pi + min_pa) > angle_tolerance) )
    % 'center' not enclosed'
    % try to find points near the poloidal_angle = 0 that have different radii
    % and use them to determine a better 'center'.
    angle_tolerance_2 = 1*pi/180; % use degree 
    ind_near_0 = find(abs(poloidal_angle_2) < angle_tolerance_2);
    min_R = min(R(ind_near_0));
    max_R = max(R(ind_near_0));
    R_center = 0.5*(max_R + min_R);
    % recalculate the poloidal angle based on this new center
    poloidal_angle = atan2( (Z - Z_center), (R - R_center) );
end

[~, ind_sorted] = sort(poloidal_angle);
Z_sorted = Z(ind_sorted);
R_sorted = R(ind_sorted);
% New way ends here

% New(er) way starts here
% Assume 1st point has z=0
if (Z(1) ~= 0)
    disp('<----Problem with new(er) method.')
else
    % Create two arrays.  One will be appended wth the 'nearest neighbor',
    % and the other one will have the nearest neighbor removed from it
    % Assumes the first point is at the 'largest' R value - Surfaces that
    % are constructed differently will need to  be 'fixed' somehow (not
    % handeled here)
    ind_ordered = 1; 
    ind_available = 2:length(Z);
    ind_original = 1:length(Z);
    for ii = 2:length(Z)
        % find the index of the nearest neighborin the ORIGINAL list
        ind_nearest = find_nearest(R(ind_ordered(ii-1)), ...
            Z(ind_ordered(ii-1)), R(ind_available), Z(ind_available), ...
            R, Z, ii-1);
        ind_ordered = [ind_ordered ind_nearest];
        ind_available = setdiff(ind_available, ind_nearest);
    end
    R_sorted = R(ind_ordered);
    Z_sorted = Z(ind_ordered);
end


% New(er) way ends here

%[Z_apex, index_apex] = max(Z);
%[Z_bottom, ~] = min(Z);
[Z_apex, index_apex] = max(Z_sorted);
[Z_bottom, index_bottom] = min(Z_sorted);

R_apex = R(index_apex);

% for the 'Old way'
%index_inboard = find(R <= R_apex);
%index_outboard = find(R >= R_apex);

%index_inboard = index_apex:index_bottom;
%index_outboard = [1:index_apex index_bottom:length(Z_sorted)];
% for the 'New way'
%index_inboard = [1:index_bottom index_apex:length(Z_sorted)];
%index_outboard = index_bottom:index_apex;

% for the 'New(er) way'
% Check to make sure sets are ordered the correct way
if (index_apex < index_bottom) 
    index_inboard = [1:index_apex index_bottom:length(Z_sorted)];
    index_outboard = index_apex:index_bottom;
else 
    index_inboard = index_bottom:index_apex;
    index_outboard = [index_apex:length(Z_sorted) 1:index_bottom];
end

% use the mean value of the two sets to determine inboard/outboard
if (mean(R_sorted(index_inboard)) < mean(R_sorted(index_outboard)))
    R_inboard = R_sorted(index_inboard);
    R_outboard = R_sorted(index_outboard);
    Z_inboard = Z_sorted(index_inboard);
    Z_outboard = Z_sorted(index_outboard);
else
    R_inboard = R_sorted(index_outboard);
    R_outboard = R_sorted(index_inboard);
    Z_inboard = Z_sorted(index_outboard);
    Z_outboard = Z_sorted(index_inboard);    
end



% Perform the polyfit
poly_inboard = polyfit(Z_inboard, R_inboard, POLYFIT_ORDER);
poly_outboard = polyfit(Z_outboard, R_outboard, POLYFIT_ORDER);

% check the fit properties
if ( min([length(R_inboard) length(R_outboard)]) < POLYFIT_ORDER)
    warning('Polynomial order is too high in calculateFlux');
    flux = -1;
    return;
end

% Make a diagnostic plot, if desired.
if MAKE_PLOT
    z_dense = linspace(Z_bottom - 0.01, Z_apex + 0.01, 1e4);
    figure(1001); hold on;
    plot(R_inboard, Z_inboard, 'r.');
    plot(R_outboard, Z_outboard, 'b.');
    plot(R_center, Z_center, 'k+');
    plot(polyval(poly_inboard, z_dense), z_dense,'r');
    plot(polyval(poly_outboard, z_dense), z_dense,'b');
    axis([min(R_inboard)*0.9 max(R_outboard)*1.1 ...
        min([Z_inboard' Z_outboard'])*1.1 max([Z_inboard' Z_outboard'])*1.1]);
    %axis([1.3 1.55 -.3 .3]);
    pause(1);  % Allows time to update the figure on the screen.
end

% Step 2: Calculate the surface integral to determine the flux through the
% surface. A buffer of 5e-4 meters seems to work well.
%flux = dblquad(@flux_integrand, Z_bottom - 5e-4, Z_apex + 5e-4, ...
%    min(R) - 5e-4, max(R) + 5e-4, TOL_FLUX, [], ...
%    poly_inboard, poly_outboard, R_apex, coilsetID, current);
flux_helper = @(Z, R) flux_integrand2(Z, R, poly_inboard, poly_outboard, R_apex, coilsetID, current);
flux = integral2(flux_helper, Z_bottom - 5e-4, Z_apex + 5e-4, ...
    min(R) - 5e-4, max(R) + 5e-4, 'AbsTol', 1e-8, 'RelTol', 1e-4);


% ++++++++++++++++
% Helper functions
% ++++++++++++++++

function dflux = flux_integrand(Z, R, poly_inboard, poly_outboard, ...
    R_apex, coilsetID, current)
% This function calculates the field at points given by the array 'Z' and
% the scalar 'R'. The two polynomial coefficent arrays, poly_inboard and
% poly_outboard, define the 'inboard' and 'outboard' boundaries of the
% surface integral. Points outside of these boundaries do not contribute to
% the integral.
% R_apex is the radius of the 'peak' or 'apex' of the flux surface.

% Initialize values
dflux = zeros(size(Z)); 

% Assume the points are at the box-port 'D' shaped region.
Phi = 0;

if (R < R_apex)  % Consider the inboard points
    R_fitline = polyval(poly_inboard, Z);    
    index_in = find(R >= R_fitline);   % points inside the flux surface
else  % Consider the outboard points
    R_fitline = polyval(poly_outboard, Z);    
    index_in = find(R <= R_fitline); % points inside the flux surface
end

% Points (R, Z(index_in)) are 'inside' the flux surface.
% Loop over these points and calculate the field at these points.
if ( ~isempty(index_in) )
    for ii = 1:length(index_in)
        % [~, B_y, ~] = calc_b_HSX_RPhiZ(R, ...
        %     Phi, Z(index_in(ii)), current, taper);
        % dflux(index_in(ii)) = B_y;
        %[~, dflux(index_in(ii)), ~] = calc_b_HSX_RPhiZ(R, ...
        %    Phi, Z(index_in(ii)), current, taper);
        [~, dflux(index_in(ii)), ~, ~, ~] = calc_b_RPhiZ(coilsetID, R, Phi, Z(index_in(ii)), current);
    end
end

function dflux = flux_integrand2(Z, R, poly_inboard, poly_outboard, ...
    R_apex, coilsetID, current)
% This function calculates the field at points given by the matrices 'Z'
% and 'R'. The two polynomial coefficent arrays, poly_inboard and
% poly_outboard, define the 'inboard' and 'outboard' boundaries of the
% surface integral. Points outside of these boundaries do not contribute to
% the integral.
% R_apex is the radius of the 'peak' or 'apex' of the flux surface.

% Initialize values
dflux = zeros(size(Z)); 

% Assume the points are at the box-port 'D' shaped region.
Phi = 0;

% if (R < R_apex)  % Consider the inboard points
R_fitline_ib = polyval(poly_inboard, Z);    
index_in1 = find(R >= R_fitline_ib);   % points inside the flux surface
% else  % Consider the outboard points
R_fitline_ob = polyval(poly_outboard, Z);    
index_in2 = find(R <= R_fitline_ob); % points inside the flux surface
index_in = intersect(index_in1, index_in2);
%end

% Points (R, Z(index_in)) are 'inside' the flux surface.
% Loop over these points and calculate the field at these points.
if ( ~isempty(index_in) )
    for ii = 1:length(index_in)
        % [~, B_y, ~] = calc_b_HSX_RPhiZ(R, ...
        %     Phi, Z(index_in(ii)), current, taper);
        % dflux(index_in(ii)) = B_y;
        %[~, dflux(index_in(ii)), ~] = calc_b_HSX_RPhiZ(R, ...
        %    Phi, Z(index_in(ii)), current, taper);
        [~, dflux(index_in(ii)), ~, ~, ~] = calc_b_RPhiZ(coilsetID, ...
            R(index_in(ii)), Phi, Z(index_in(ii)), current);
    end
end


function ind_nearest_orig = find_nearest(R, Z, R_avail, Z_avail, ...
    R_orig, Z_orig, ii_last)

if (ii_last == 1) % handle the first point carefully- force the search to be in the +Z direction
    distance = sqrt( (R - R_avail).^2 + (Z - Z_avail).^2);
    ind_Z_negative = find(Z_avail < 0);
    distance(ind_Z_negative) = inf;    
else
    distance = sqrt( (R - R_avail).^2 + (Z - Z_avail).^2);
end

[~, ind_nearest] = min(distance);
% now, find the index of that point in the ORIGINAL list
ind_nearest_orig = 0;
for ii = 1:length(R_orig)
    if ( (R_avail(ind_nearest) == R_orig(ii)) & ...
            (Z_avail(ind_nearest) == Z_orig(ii)) )
        ind_nearest_orig = ii;
    end
end




