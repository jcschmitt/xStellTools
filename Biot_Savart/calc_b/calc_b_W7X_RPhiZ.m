function [Bx, By, Bz, BR, BPhi] = ...
    calc_b_HSX_RPhiZ(P_R, P_Phi, P_Z, coil_current_array)
% function [Bx,By,Bz,Br,Bphi] = calc_b_W7X(P_R, P_Phi, P_Z, main_current,
% taper)
%
% Input: P_R, P_Phi, P_Z:  Cylindrical coordinates of the observation point
%        main_current:  The current in the main field coils (per turn)
%        taper: 'Taper Array' to define the ratio of the coil courrent in
%        the auxilliarry coilset: 
%               [Coil_1 Coil_2 Coil_3 Coil_4 Coil_5 Coil_6]
%  
%
% Output:  Bx, By, Bz, BR, BPhi:  Magnetic flux density (T/m) at
%          observation point (aka Magnetic field strength)
%
% Calculates cartesian B field from given cylindrical lab frame
% coordinates and coil current array.  If the observation point lies
% on a current filament (or extremely close to it) then the contribution
% from that filament is ignored.
%
% The HSX coil fileset (coils.mat aux_coils.mat) must be available to this
% function
%

% Convert from cylindrical to cartesian coordinates, and call the
% 'cartesian' version for HSX
P_x = P_R .* cos(P_Phi);
P_y = P_R .* sin(P_Phi);
P_z = P_Z;

% Could have used Matlab internal:
%[P_x,P_y,P_z]=pol2cart(P_Phi,P_R,P_Z);    %convert to cartesian lab space

[Bx, By, Bz, BR, BPhi] = calc_b_W7X(P_x, P_y, P_z, coil_current_array);


