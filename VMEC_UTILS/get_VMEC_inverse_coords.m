function [r, z, r_rho, r_theta, r_phi, z_rho, z_theta, z_phi, ...
    lambda_phi, lambda_theta, b_r, b_z, b_tor, sqrt_g] = ...
    get_VMEC_inverse_coords(rho, theta, phi, mnmax, xm, xn, ...
    spline_vmec, spline_d_vmec, torflux_LCFS)
%***********************************************************************
%get_VMEC_inverse_coords finds r and z and their derivatives with respect to rho, theta,
%   and phi at the given rho, theta, and phi.  It also returns the
%   cylindrical components of the magnetic field and the jacobian.
%References:
%  Attenberger, Houlberg, Hirshman J Comp Phys 72 (1987) 435
%  W.A.Houlberg 3/98
%  Original Fortran (tragrz) by W.A. Houlberg 12/97
%  Adapted by JC Schmitt for Matlab  11/2009
%Input:
%  rho-radial flux coordinate (-)
%  theta-poloidal angle coordinate (radians)
%  phi-toroidal angle coordinate (radians)
%  mnmax-total number of modes
%  xm-poloidal mode #s
%  xn-toroidal mode #s
%  spline_vmec - spline coeffiecinets structure for each mode iota, lmns,
%       rmnc, zmns
%  spline_d_vmec - spline coeffiecinets structure for radial derivatives
%       of each mode of rmnc, zmns
%  torflux_LCFS-The enclosed toroidal flux at the last closed flux surface.
%Output:
%  r-radius of point from major axis (m)
%  z-distance of point from midplane (m)
%  r_rho-dr/dx (m)
%  r_theta-dr/dtheta (m/radians)
%  r_phi-dr/dphi (m/radians)
%  z_rho-dz/dx (m)
%  z_theta-dz/dtheta (m/radians)
%  z_phi-dz/dphi (m/radians)
%  lambda_phi  d(lambda)/d(phi)
%  lambda_theta d(lambda)/d(theta)
%  br-b dot grad(r) (T)
%  bz-b dot grad(z) (T)
%  btor-toroidal field (T)
%  sqrt_g
%Comments:
%  A factor of rho**m is factored out of the Fourier coefficients prior
%    to spline fitting for increased accuracy. For rho less than
%    rho_thrift(1), all terms are assumed to vary as rho**m
% For points eyond the LCFS, the m=1,n=0  (m=-1,n=0) term is used to
%  extend the values of r and z, while the rest of the term are held to
%  their values on the LCFS.  The allows a unique flux representation of
%  all space for convex, not bean-shaped plasmas.  This does hold exactly
%  for HSX.
%***********************************************************************

z_precision = 2.0e-9;
DELTA_RHO = 1e-7;

r=0.0;
z=0.0;
r_rho=0.0;
r_theta=0.0;
r_phi=0.0;
z_rho=0.0;
z_theta=0.0;
z_phi=0.0;
lambda_phi = 0.0;
lambda_theta = 0.0;

% get the spline-fitted values
if rho > 1
    rho_to_use = 1;
else
    rho_to_use = rho;
end

spline_vals_vmec = ppval(spline_vmec, rho_to_use);
spline_dvals_vmec = ppval(spline_d_vmec, rho_to_use);

% extract the values from the output matrix
iota = spline_vals_vmec(1);
ind_spline_start = 2;
ind_spline_end = ind_spline_start + mnmax - 1;
lmns = spline_vals_vmec(ind_spline_start:ind_spline_end);
ind_spline_start = ind_spline_end + 1;
ind_spline_end = ind_spline_start + mnmax - 1;
rmnc = spline_vals_vmec(ind_spline_start:ind_spline_end);
ind_spline_start = ind_spline_end + 1;
ind_spline_end = ind_spline_start + mnmax - 1;
zmns = spline_vals_vmec(ind_spline_start:ind_spline_end);

% extract the values from the output matrix
ind_spline_start = 1;
ind_spline_end = ind_spline_start + mnmax - 1;
drmnc = spline_dvals_vmec(ind_spline_start:ind_spline_end);
ind_spline_start = ind_spline_end + 1;
ind_spline_end = ind_spline_start + mnmax - 1;
dzmns = spline_dvals_vmec(ind_spline_start:ind_spline_end);

ind_m_eq_1_n_eq_0 = NaN;
% loop over each mode
for ii = 1:mnmax
    % find the m=1 or -1, n=0 term for later
    if (xm(ii) == 1 && xn(ii) == 0)
        ind_m_eq_1_n_eq_0 = ii;
    end
    
    % find angle matrix, get cosin and sin
    angle_ = xm(ii) * theta - xn(ii) * phi;
    cosangle_ = cos(angle_);
    sinangle_ = sin(angle_);
    
    % get normalization factor
    if rho == 0 || rho > 1;
        rho_m_norm = 1;
    else
        rho_m_norm = rho^abs(xm(ii));
    end
    
    % add the mode's contribution to r and z
    r = r + rmnc(ii) * cosangle_ * rho_m_norm;
    z = z + zmns(ii) * sinangle_ * rho_m_norm;
    
    % add terms for derivaties
    r_theta = r_theta - xm(ii) * rmnc(ii) * sinangle_ * rho_m_norm;
    r_phi = r_phi + xn(ii) * rmnc(ii) * sinangle_ * rho_m_norm;
    z_theta = z_theta + xm(ii) * zmns(ii) * cosangle_ * rho_m_norm;
    z_phi = z_phi - xn(ii) * zmns(ii) * cosangle_ * rho_m_norm;
    
    % the radial terms require more precision.  use the spliens
    drmn_rho = drmnc(ii);
    dzmn_rho = dzmns(ii);
    
    if rho > (10*z_precision)
        rho_dumm = rho;
    else
        rho_dumm = 10*z_precision;
    end
    drho_m = abs(xm(ii)) * rho_dumm^(abs(xm(ii))-1);
    r_rho = r_rho + drmn_rho * cosangle_ * rho_m_norm + ...
        rmnc(ii) * cosangle_* drho_m;
    z_rho = z_rho + dzmn_rho * sinangle_ * rho_m_norm + ...
        zmns(ii) * sinangle_ * drho_m;
    
    lambda_0 = lmns(ii);
    % unnormalize lambda_0 and add the mode's contribution
    lambda_0 = lambda_0 * rho_m_norm;
    lambda_theta = lambda_theta + lambda_0 * xm(ii) * cosangle_;
    lambda_phi = lambda_phi - lambda_0 * xn(ii) * cosangle_;
end

if rho > 1
    % add on a contribution from the extension of ONLY the m=1, n=0 term
    angle_ =  theta;
    cosangle_ = cos(angle_);
    sinangle_ = sin(angle_);
    r_rho = rmnc(ind_m_eq_1_n_eq_0) * cosangle_;
    z_rho = zmns(ind_m_eq_1_n_eq_0) * sinangle_;
    r = r + (rho - 1) * rmnc(ind_m_eq_1_n_eq_0) * cosangle_;
    z = z + (rho - 1) * zmns(ind_m_eq_1_n_eq_0) * sinangle_;
    r_theta = r_theta - (rho - 1) * rmnc(ind_m_eq_1_n_eq_0) * ...
        sinangle_;
    z_theta = z_theta + (rho - 1) * zmns(ind_m_eq_1_n_eq_0) * ...
        cosangle_;
end

% Find jacobian and magnetic field components
tau = z_rho * r_theta - r_rho * z_theta;
sqrt_g = r * tau;
phiprm = 2 * rho * torflux_LCFS;
b_r = phiprm * ((iota - lambda_phi) * r_theta + ...
    (1 + lambda_theta) * r_phi) / (2*pi * sqrt_g);
b_z = phiprm * ((iota - lambda_phi) * z_theta + ...
    (1 + lambda_theta) * z_phi) / (2*pi * sqrt_g);
b_tor = r * phiprm * (1 + lambda_theta) / (2*pi * sqrt_g);

