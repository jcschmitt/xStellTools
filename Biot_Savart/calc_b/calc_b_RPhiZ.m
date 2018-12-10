function [Bx, By, Bz, BR, BPhi] = ...
    calc_b_RPhiZ(coilset, P_R, P_Phi, P_Z, coil_current_array)
% function [Bx,By,Bz,Br,Bphi] = calc_b_QHS46_RPhiZ(P_R, P_Phi, P_Z,
% coil_current_array)
%
% Input: P_R, P_Phi, P_Z:  Cylindrical coordinates of the observation point
%        coil_current_array:  The current in the field coils (per turn)
%
% Output:  Bx, By, Bz, BR, BPhi:  Magnetic flux density (T/m) at
%          observation point (aka Magnetic field strength)
%
% Calculates cartesian B field from given cylindrical lab frame
% coordinates and coil current array.  If the observation point lies
% on a current filament (or extremely close to it) then the contribution
% from that filament is ignored.
%
% The QHS46 coil fileset (coils.QHS46_xxx) must be available to this
% function
%

% Convert from cylindrical to cartesian coordinates, and call the
% 'cartesian' version for QSH46
P_x = P_R .* cos(P_Phi);
P_y = P_R .* sin(P_Phi);
P_z = P_Z;

% Could have used Matlab internal:
%[P_x,P_y,P_z]=pol2cart(P_Phi,P_R,P_Z);    %convert to cartesian lab space

[Bx, By, Bz, BR, BPhi] = calc_b(coilset, P_x, P_y, P_z, coil_current_array);


