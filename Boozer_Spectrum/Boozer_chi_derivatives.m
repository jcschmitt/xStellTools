function [dxdz]=chi_Boozer_derivs(Chi, Position)
% This returns the ode for a field line as a function of chi
% the integral of B*dl
% which can be solved with one of matlab's ode solvers.
%
% input:
%   Chi: the independent variable to be integrate over
%   Position: column vector containing r, phi, z (Position = [r phi z])
%
% output:
%   dxdz: row vector containg the derivatives of r and z w/ respect to phi
%       (dxdz = [drdchi dphidchi dzdphi]')
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global current taper
R = Position(1);
phi = Position(2);
Z = Position(3);

[bx, by, bz] = calc_b_bs(R, phi, Z, current, taper);

sinphi=sin(phi);
cosphi=cos(phi);

% Convert to cyl coords
br = bx*cosphi + by*sinphi;
bphi = -bx*sinphi + by*cosphi;

B=sqrt(bx.*bx+by.*by+bz.*bz);

dRdChi=br./(B.*B);
dPhidChi=bphi./(R.*B.*B);
dZdChi=bz./(B.*B);

dxdz=[dRdChi dPhidChi dZdChi]';    
