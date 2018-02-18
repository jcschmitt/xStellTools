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
% Modified by JC Schmitt, LTX, W7X, HSX Auburn (2011 - 2017)

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
Bx = 0; By = 0; Bz = 0;

% save some time by checking for zero current in the coil
if coil_current == 0
    return;
end

% Loop over each turn of the coil - assumes coils are 'closed', i.e. that
% the first and last points are the same physical point
for ii = 1:coil.num_turns
    % npts = coil.turn_number(ii).num_vertices;
    
    R_x = P_x - coil.turn_number(ii).x;
    R_y = P_y - coil.turn_number(ii).y;
    R_z = P_z - coil.turn_number(ii).z;
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
        if ~isempty(A)
            disp(['<----Close to a coil: (x,y,z) = (' ...
                num2str(P_x) ', ' num2str(P_y) ', ' ...
                num2str(P_z) ')']);
        end
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
    
    Bx = Bx + sum(RfcrRi_x .* back);
    By = By + sum(RfcrRi_y .* back);
    Bz = Bz + sum(RfcrRi_z .* back);
end

