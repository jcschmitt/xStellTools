function [Ax, Ay, Az, Bx, By, Bz] = ...
    calc_AB_BiotSavart(P_x, P_y,P_z, coil, coil_current, CALC_B)
% function [Ax, Ay, Az, Bx, By, Bz] = ...
%     calc_AB_BiotSavart(P_x, P_y,P_z, coil, coil_current, CALC_B)
% This function calculates the magnetic vector potential, 'A', and, if
% requested, it will also calculate the magnetic field strength, 'B'

muI_div_4pi = coil_current * 1e-7; % this is mu0 * I / 4*pi

Ax = 0; Ay = 0; Az = 0; Bx = 0; By = 0; Bz = 0;  % Initialize the values

% save some time by checking for zero current in the coil
if coil_current == 0
    return;
end

% Loop over each turn of the coil - assumes coils are 'closed', i.e. that
% the first and last points are the same physical point
for ii = 1:coil.num_turns
    % Points per turn
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
    
    Lx = Ri_x - Rf_x;
    Ly = Ri_y - Rf_y;
    Lz = Ri_z - Rf_z;
    L = sqrt(Lx.^2 + Ly.^2 + Lz.^2);
    
    ex = Lx ./ L;
    ey = Ly ./ L;
    ez = Lz ./ L;
    
    Ri_p_Rf = (Ri_mag + Rf_mag);
    epsilon = L ./ Ri_p_Rf;
    ln_arg = log(( 1 + epsilon ) ./ ( 1 - epsilon));
    
    Ax = Ax + sum(ex .* ln_arg);
    Ay = Ay + sum(ey .* ln_arg);
    Az = Az + sum(ez .* ln_arg);
    
    if CALC_B
        % cp is 'cross product'
        % L * (e cross Ri) is same as
        %    (Lx, Ly, Lz) cross (Ri_x, Ri_y, Ri_z)
        cp_x = Ly .* Ri_z - Lz .* Ri_y;
        cp_y = Lz .* Ri_x - Lx .* Ri_z;
        cp_z = Lx .* Ri_y - Ly .* Ri_x;
        
        Bx = Bx + sum( cp_x .* Ri_p_Rf ./ ...
            ( (Ri_mag .* Rf_mag) .* (Ri_p_Rf .^2 - L .^2) ) );
        By = By + sum( cp_y .* Ri_p_Rf ./ ...
            ( (Ri_mag .* Rf_mag) .* (Ri_p_Rf .^2 - L .^2) ) );
        Bz = Bz + sum( cp_z .* Ri_p_Rf ./ ...
            ( (Ri_mag .* Rf_mag) .* (Ri_p_Rf .^2 - L .^2) ) );
    end
    
end

Ax = muI_div_4pi * Ax;
Ay = muI_div_4pi * Ay;
Az = muI_div_4pi * Az;

if CALC_B
    Bx = 2 * muI_div_4pi * Bx;
    By = 2 * muI_div_4pi * By;
    Bz = 2 * muI_div_4pi * Bz;
end
