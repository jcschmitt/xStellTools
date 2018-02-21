function [dX_dchi] = Boozer_chi_derivatives_HSX(chi, Position, current, taper)
% function [dX_dchi] = Boozer_chi_derivatives(chi, Position, current, taper)
%
% This returns the ode for a field line as a function of chi
% (chi = the integral of |B|*dl along the line)r
% 
% This is inteded to be solved with one of matlab's ode solvers.
%
% Input:
%   Chi: the independent variable to be integrated over (not used)
%   Position: column vector containing r, phi, z (Position = [r phi z])
%
% Output:
%   dX_dchi: row vector containg the derivatives of r and z w/ respect to
%   chi : (dX_dchi = [drdchi dphidchi dzdchi]')
%
% Requirements: 
%   The calculation of the magnetic field is required. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = Position(1);
phi = Position(2);
Z = Position(3);

[B_X, B_Y, B_Z, B_R, B_phi] = calc_b_HSX_RPhiZ(R, phi, Z, current, taper);

% Calculate |B|^2
%B=sqrt(bx.*bx+by.*by+bz.*bz);
Bsq = (B_X.*B_X + B_Y.*B_Y + B_Z.*B_Z);

% These are the 'dX / dchi' components
dR_dchi = B_R ./ Bsq;
dPhi_dchi = B_phi ./ (R .* Bsq);
dZ_dchi = B_Z ./ Bsq;

% Arrange them as needed for the ode solver.
dX_dchi = [dR_dchi dPhi_dchi dZ_dchi]';    
