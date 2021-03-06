function [Bx, By, Bz] = calc_B_BiotSavart(P_x, P_y, P_z, ...
    coil, coil_current)
% function [Bx, By, Bz] = calc_B_BiotSavart(P_x, P_y, P_z, ...
%     coil, coil_current)
% Input: Px, Py, Pz:  Cartesian coordinates of the observation point
%        coil_current: the coil current
%
% Output:  Bx, By, Bz:  Magnetic flux density (T/m) at observation point
%
% Calculates cartesian B field from given cartesian lab frame
% coordinates and coil current.  If the observation point lies
% on a current filament (or extremely close) the contribution from that
% filament is ignored.
%
% Note that you can comment out the on filament check for speed, make sure
% this has not been done if you are looking near the coils!
%
% Method from [J. Hanson, S. Hirshman PoP 9, 4410 (2002)]
% Original Matlab version by J.Lore, HSX, UW-Madison (2006)
% Modified by JC Schmitt, LTX, W7X, HSX Auburn (2011 - 2019)

% coil: structure defining the coil
% coil.num_turns : number of turns, or loops in the coil
% coil.turn_number(:).num_vertices : number of vertices for a turn. Does
%           not need to be equal for all turn_numbers
% coil.turn_number(:).x = x coordinates of the vertics, in meters
% coil.turn_number(:).y = y coordinates of the vertics, in meters
% coil.turn_number(:).z = z coordinates of the vertics, in meters
%

% tolerance for checking if the test point is on a current filament
onfil_tol = 1e-10;
I = coil_current * 1e-7; % this is mu0 * I / 4*pi

% Initialize the values
if 0
    % this section is under development for handling of non-standard input
    % sizes
    [ndims_x, dims_x, total_x] = find_matrix_parameters(P_x);
    [ndims_y, dims_y, total_y] = find_matrix_parameters(P_y);
    [ndims_z, dims_z, total_z] = find_matrix_parameters(P_z);
    
    
    % not validated for all types of size combinations, but it does handle
    % points, lines, planes and volumes
    
    % Determine the case, initalize the B components and determine the final
    % output structure dimension and sizes
    % Modify any P_x, P_y, P_z matrices as neccessary.
    total_point_count = max([total_x, total_y, total_z]);
    if ((total_x == 1) && (total_y == 1) && (total_z == 1))
        % Individual point
        dimensions_out = [1 1];
    elseif ((ndims_x == 1) && (ndims_y == 1) && (ndims_z == 1))
        % Handle lines
        if ( (total_x == total_y) && (total_x == total_z) )
            dimensions_out = dims_x;
        elseif ( (total_x == total_y) && (total_z == 1) )
            dimensions_out = dims_x;
            P_z = P_z * ones(dims_x);
        elseif ( (total_x == total_z) && (total_y == 1) )
            dimensions_out = dims_x;
            P_y = P_y * ones(dims_x);
        elseif ( (total_y == total_z) && (total_x == 1) )
            dimensions_out = dims_y;
            P_x = P_x * ones(dims_x);
        else
            disp('K=====Error in calc_B_BiotSavart');
        end
        
    elseif ((ndims_x == ndims_y) && (ndims_z == 1))
        % Handle X-Y planes
        dimensions_out = dims_x;
        P_z = P_z * ones(dims_x);
    elseif ((ndims_x == ndims_z) && (ndims_y == 1))
        % Handle X-Z planes
        dimensions_out = dims_x;
        P_y = P_y * ones(dims_x);
    elseif ((ndims_y == ndims_z) && (ndims_x == 1))
        % Handle Y-Z planes
        dimensions_out = dims_y;
        P_x = P_x * ones(dims_y);
    elseif ((ndims_x == ndims_y) && (ndims_x == ndims_z))
        % Handle X-Y-Z volumes
        if ( (total_x == total_y) && (total_x == total_z) )
            dimensions_out = dims_x;
        elseif ( (total_x == total_y) && (total_z == 1) )
            dimensions_out = dims_x;
            P_z = P_z * ones(dims_x);
        elseif ( (total_x == total_z) && (total_y == 1) )
            dimensions_out = dims_x;
            P_y = P_y * ones(dims_x);
        elseif ( (total_y == total_z) && (total_x == 1) )
            dimensions_out = dims_y;
            P_x = P_x * ones(dims_x);
        else
            disp('K=====Error in calc_B_BiotSavart');
        end
    end
    
    % Reshape the (modified) matrices for the loop
    newP_x = reshape(P_x, total_point_count, 1);
    newP_y = reshape(P_y, total_point_count, 1);
    newP_z = reshape(P_z, total_point_count, 1);
    
    % These have the same dimension as the reshaped matrices.  Will need to
    % undo the reshaping at the end.
    
    Bx = zeros(total_point_count, 1);
    By = zeros(total_point_count, 1);
    Bz = zeros(total_point_count, 1);
else
    Bx = 0;
    By = 0;
    Bz = 0;
end

% save some time by checking for zero current in the coil
if coil_current == 0
    return;
end

% Loop over each point in space and find the Bx, By, Bz for that point
% for jj = 1:total_point_count
%     thisP_x = newP_x(jj);
%     thisP_y = newP_y(jj);
%     thisP_z = newP_z(jj);
thisP_x = P_x;
thisP_y = P_y;
thisP_z = P_z;

% Loop over each turn of the coil - assumes coils are 'closed', i.e. that
% the first and last points are the same physical point
for ii = 1:coil.num_turns
    % npts = coil.turn_number(ii).num_vertices;
    
    R_x = thisP_x - coil.turn_number(ii).x;
    R_y = thisP_y - coil.turn_number(ii).y;
    R_z = thisP_z - coil.turn_number(ii).z;
    R_mag = sqrt(R_x.^2 + R_y.^2 + R_z.^2);
    
    % R_i is from beginning of filament to point x
    Ri_x = R_x(1:end-1);
    Ri_y = R_y(1:end-1);
    Ri_z = R_z(1:end-1);
    Ri_mag = R_mag(1:end-1);
    
    % R_f is from end of filament to x
    Rf_x = R_x(2:end);
    Rf_y = R_y(2:end);
    Rf_z = R_z(2:end);
    Rf_mag = R_mag(2:end);
    
    % this is the magnitude of the vector along the filament squared (for speed)
    e_mag_sq = (Ri_x - Rf_x).^2 + (Ri_y - Rf_y).^2 + (Ri_z - Rf_z).^2;
    
    % cross product: RfxRi
    RfcrRi_x = Rf_y .* Ri_z - Rf_z .* Ri_y;
    RfcrRi_y = Rf_z .* Ri_x - Rf_x .* Ri_z;
    RfcrRi_z = Rf_x .* Ri_y - Rf_y .* Ri_x;
    
    eps_sq = e_mag_sq ./ (Ri_mag + Rf_mag).^2;
    
    % The common part of B equation
    back = -2 * I ./ (Rf_mag .* Ri_mag .* (Ri_mag + Rf_mag));
    
    % Uncomment this if you are going to calculate things near the coils
    % check for test point on a filament, within tolerance 'onfil_tol'
    while find(abs(eps_sq-1) < onfil_tol);
        A = find(abs(eps_sq-1) < onfil_tol);
        % Uncomment this for debugging.
        % if ~isempty(A)
        %     disp(['<----Close to a coil: (x,y,z) = (' ...
        %         num2str(P_x) ', ' num2str(P_y) ', ' ...
        %         num2str(P_z) ')']);
        % end
        % this is a more stringent test, you can reduce the
        % tolerance if this is used, but it adds two extra sqrt
        % calculations
        %   while find(abs(sqrt(eps_sq)-1)<onfil_tol);
        %       A=find(abs(sqrt(eps_sq)-1)<onfil_tol);
        % if such a filament found, set contribution (via back)
        % to 0 and change eps_sq to a 'large number' so search can continue
        back(A) = 0;
        eps_sq(A) = 2;
    end
    
    back = back ./ (1 - eps_sq);
    
    %         Bx(jj) = Bx(jj) + sum(RfcrRi_x .* back);
    %         By(jj) = By(jj) + sum(RfcrRi_y .* back);
    %         Bz(jj) = Bz(jj) + sum(RfcrRi_z .* back);
    Bx = Bx + sum(RfcrRi_x .* back);
    By = By + sum(RfcrRi_y .* back);
    Bz = Bz + sum(RfcrRi_z .* back);
end
% end

% use the 'correct' shape
% Bx = reshape(Bx, dimensions_out);
% By = reshape(By, dimensions_out);
% Bz = reshape(Bz, dimensions_out);
