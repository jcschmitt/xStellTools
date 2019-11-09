% See calc_descur_RZ_points.m

% function [surf_x,surf_y,surf_z,phitemp] = make_lcfs(rboundary,nphi,ntheta, coilsetID, coilCurrents)
% % inputs:
% %   rboundary- R coord of boundary at the boxport
% %   nphi- number of toridal cuts
% %   ntheta- number of points on each cut
% %   current- coil current in amps (- for counter-clockwise)
% %   taper- taper array for relative magnitudes of the aux coil ampturns
% %
% % outputs:
% %   surf_x,surf_y,surf_z- x,y,z coords of the points on the surface
% %   phitemp-the phi coordinate of each point
% %
% %                                             1/07 JL mod of JC file
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Keep track of how long this takes
% tic;
% 
% % Starting position of field line to follow [r0 z0]
% x0=[rboundary 0];
% 
% % Values of phi at which B is to be calculated. 
% phi_points=0 : 2*pi/nphi : ntheta*2*pi;
% 
% 
% % Integrate the ode from field_line_derivs.
% % NOTE: the tolerances here make a big difference, so be careful!
% options = odeset('RelTol',1e-10,'AbsTol',1e-10);
% 
% %x is [r z]
% [phitemp,xtemp] = ode45(@field_line_derivs,phi_points,x0,options, coilsetID, coilCurrents); 
% rsurf=xtemp(:,1);
% zsurf=xtemp(:,2);
% surf_x = [rsurf.*cos(phitemp)]';
% surf_y = [rsurf.*sin(phitemp)]';
% surf_z = [zsurf]';
% 
% 
% if 0
% % plot the surface
% 
% figure
% hold on
% %ph = plot3(rsurf.*cos(phi),rsurf.*sin(phi),zsurf,'.');
% plot(rsurf,zsurf,'.')
% axis equal
% hold off
% end
% 
% toc
% 
% 
% 
% function [dxdz] = field_line_derivs(phi, x, coilsetID, coilCurrents)
% % function [dxdz] = field_line_derivs(phi,x)
% %
% % This returns the ode for a field line, which can be solved with one 
% % of matlab's ode solvers.
% %
% % input:
% %   phi: the independent variable to be integrate over
% %   x: column vector containing r and z (x = [r z])
% %
% % output:
% %   dxdz: row vector containg the derivatives of r and z w/ respect to 
% %         phi. (dxdz = [drdphi dzdphi]')
% %
% %
% %                                               5/3/02 John Canik
% %                                   Updated 2015-2019, John Schmitt
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% r = x(1);
% z = x(2);
% 
% [bx, by, bz] = calc_b_RPhiZ(coilsetID, r, phi, z, coilCurrents);
% 
% 
% br = bx*cos(phi) + by*sin(phi);
% bphi = -bx*sin(phi) + by*cos(phi);
% 
% % Return derivs
% drdphi = (r.*br)./bphi;
% dzdphi = (r.*bz)./bphi;
% dxdz(1) = drdphi;
% dxdz(2) = dzdphi;
% dxdz=dxdz';