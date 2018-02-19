function [flux] = calc_flux_HSX(R, Z, current, taper, MAKE_PLOT)
% Returns the flux through a surface defined by [R, Z].  Assumes the
% surface is at the boxport of HSX (at the 'D' shaped toroidal plane).  

% If MAKE_PLOT was not specified, set to 0.
if nargin < 5
    MAKE_PLOT = 0; % Default is 'no'
end

% 1e-9 seems to be a good value for convergence.
TOL_FLUX = 1e-9;
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

[Z_apex, index_apex] = max(Z);
[Z_bottom, ~] = min(Z);

R_apex = R(index_apex);
index_inboard = find(R <= R_apex);
index_outboard = find(R >= R_apex);
R_inboard = R(index_inboard);
R_outboard = R(index_outboard);
Z_inboard = Z(index_inboard);
Z_outboard = Z(index_outboard);

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
    z_dense = linspace(-.25, .25, 1e4);
    figure; hold on;
    plot(R_inboard, Z_inboard, 'r.');
    plot(R_outboard, Z_outboard, 'b.');
    plot(polyval(poly_inboard, z_dense), z_dense,'r');
    plot(polyval(poly_outboard, z_dense), z_dense,'b');
    axis([1.3 1.55 -.3 .3]);
    pause(1);  % Allows time to update the figure on the screen.
end

% Step 2: Calculate the surface integral to determine the flux through the
% surface. A buffer of 5e-4 meters seems to work well.
flux = dblquad(@flux_integrand, Z_bottom - 5e-4, Z_apex + 5e-4, ...
    min(R) - 5e-4, max(R) + 5e-4, TOL_FLUX, [], ...
    poly_inboard, poly_outboard, R_apex, current, taper);



% +++++++++++++++
% Helper function
% +++++++++++++++

function dflux = flux_integrand(Z, R, poly_inboard, poly_outboard, ...
    R_apex, current, taper)
% This function calculates the field at points given by the array 'Z' and
% the scalar 'R'. The two polynomial coefficent arrays, poly_inboard and
% poly_outboard, define the 'inboard' and 'outboard' boundaries of the
% surface integral. Points outside of these boundaries do not contribute to
% the integral.
% R_apex is the radius of the 'peak' or 'apex' of the flux surface.
% This function is intended or HSX.

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
        [~, dflux(index_in(ii)), ~] = calc_b_HSX_RPhiZ(R, ...
            Phi, Z(index_in(ii)), current, taper);
    end
end

