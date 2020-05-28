function [r_out, z_out, rho_out, theta_out, r_rho_out, ...
    r_theta_out, r_phi_out, z_rho_out, ...
    z_theta_out, z_phi_out, ...
    lambda_phi_out, lambda_theta_out, ...
    b_r_out, b_z_out, b_tor_out, ...
    sqrt_g_out, bsq_out, b1_out, b_x_out, b_y_out, ...
    Psi_Phi_out, Psi_R_out, Psi_Z_out, Psi_x_out, Psi_y_out] = ...
    get_VMEC_derived_values(file_vmec_output, lab_x, lab_y, lab_z, TOL_INV, ...
    rho_guess, theta_guess)
% function [r_out, z_out, rho_out, theta_out, r_rho_out, ...
%     r_theta_out, r_phi_out, z_rho_out, ...
%     z_theta_out, z_phi_out, ...
%     lambda_phi_out, lambda_theta_out, ...
%     b_r_out, b_z_out, b_tor_out, ...
%     sqrt_g_out, bsq_out, b1_out, b_x_out, b_y_out, ...
%     Psi_Phi_out, Psi_R_out, Psi_Z_out, Psi_x_out, Psi_y_out] = ...
%     get_VMEC_derived_values(file_vmec_output, lab_x, lab_y, lab_z, TOL_INV, ...
%     rho_guess, theta_guess)
% Input
%
% The input list is:
%
%     * file_vmec_output: The complete name of the VMEC output file
%     * lab_x, lab_y, lab_z: Each is an array of points in lab coordinate space.
%     * TOL_INV: (optional) A tolerance on the inversion. Default = 1e-6;
%     * rho_guess, theta_guess: (optional) Initial guesses on the inversion. default = .5 and pi/2
%
% Output
%
% Each output is an array the same size as lab_x, lab_y, etc, and correspond to the following:
%
%     * r_out: The real space R coordinate, after inversion
%     * z_out: The real space Z coordinate, after inversion
%     * rho_out, theta_out: The vmec rho coordinate (radial) and theta (poloidal) coordinate
%     * r_rho_out, r_theta_out, r_phi_out: The derivatives of 'R' w.r.t. VMEC rho, theta, and phi.
%     * z_rho_out, z_theta_out, z_phi_out: The derivatives of 'Z' w.r.t. VMEC rho, theta, and phi
%     * lambda_phi_out, lambda_theta_out: The derivatives of the VMEC poloidal angle stream function, w.r.t. phi and theta
%     * b_r_out, b_z_out, b_tor_out: Magnetic field components from VMEC inversion point.
%     * sqrt_g_out, bsq_out, b1_out: Jacobian, B^2, and |B| at the VMEC inversion point
%     * b_x_out, b_y_out: The lab x and y projections of the VMEC output
%     * Psi_Phi_out, Psi_R_out, Psi_Z_out, Psi_x_out, Psi_y_out: The components of grad(Psi) w.r.t. Phi (lab toroidal angle), R (major radius), Z (vertical coordinate), x (laboratory x coordinate), and y (lab y coordinate).
%
% The code
%
% Requirements:
%
%    1. Although this VMEC code does require it, the DEMOs do require you to have hsx_grid_interp in your matlab path.
%    2. You need netcdf support enabled in matlab. If you can not get this to work, please see me and I'll help you get it working on your machine.
%    3. You need a vmec output file in your path. There are a few in Y:\CODES\Auxiliary VMEC and Boozer codes
%
%=============================================
% Updates
%=============================================
%
% 2009/12/02 - Added native netcdf support for Matlab 7.6+.  Replaced depracated linterp funtion with interp1. - JW Radder
% 2009/11/18 - version 1.0a - JC Schmitt
%=============================================


DO_TIME_TEST = 0;

if nargin < 7
    rho_guess = .5; theta_guess = pi/2;
end
if nargin < 5
    TOL_INV = 1e-6;
end

if nargin < 4
    disp('Not enuf inputs')
end

if isvector(lab_x) ~= 1
    disp('Sorry.  This only works on 1-d arrays');
end

persistent VMEC_DATA_LOADED
persistent VMEC_OUTPUT_FILENAME
persistent mnmax_vmec
persistent xm_vmec
persistent xn_vmec
persistent spline_vmec
persistent spline_d_vmec
persistent phitot_vmec

if ( (isempty(VMEC_DATA_LOADED)) || (VMEC_DATA_LOADED == 0) || ...
        ~strcmp(VMEC_OUTPUT_FILENAME, file_vmec_output) )
    LOAD_VMEC_DATA = 1;
    VMEC_OUTPUT_FILENAME = file_vmec_output;
    VMEC_DATA_LOADED = 0;
else
    LOAD_VMEC_DATA = 0;
end

% load VMEC_LOADED
if LOAD_VMEC_DATA
    [mnmax_vmec, xm_vmec, xn_vmec, spline_vmec, spline_d_vmec, ...
        phitot_vmec] = load_vmec_wout_file(file_vmec_output);
else
    % debug space
end


num_pts = length(lab_x);

if DO_TIME_TEST
    tic
end

% First convert lab (x,y,z) to cylindrical
lab_R = sqrt(lab_x.^2 + lab_y.^2);
lab_Phi = atan2(lab_y, lab_x);

for ii = 1:num_pts
    % find the the vmec coordinates of the point
    [rho_out(ii), theta_out(ii), r_rho(ii), ...
        r_theta(ii), r_phi(ii), z_rho(ii), ...
        z_theta(ii), z_phi(ii), output_flag(ii)] = ...
        find_VMEC_coords(rho_guess, theta_guess, ...
        lab_R(ii), lab_z(ii), lab_Phi(ii), TOL_INV, ...
        mnmax_vmec, xm_vmec, xn_vmec, spline_vmec, spline_d_vmec, ...
        phitot_vmec(end));
    if (output_flag(ii) < 100)  % normal solution found
    else
        % try a whole bunch of inital guess
        disp('<---Did not invert correctly. Looping over search grid, keeping 1st ''good'' inverted solution')
        num_attempts = 0;
        for iii = 0.9:-0.2:0.1
            for jjj = (2*pi) * linspace(0,1,5);
                num_attempts  = num_attempts +1 ;
                rho_guess = iii;
                theta_guess = jjj;
                
                [rho_out(ii), theta_out(ii), r_rho(ii), ...
                    r_theta(ii), r_phi(ii), z_rho(ii), ...
                    z_theta(ii), z_phi(ii), output_flag(ii)] = ...
                    find_VMEC_coords(rho_guess, theta_guess, ...
                    lab_R(ii), lab_z(ii), lab_Phi(ii), TOL_INV, ...
                    mnmax_vmec, xm_vmec, xn_vmec, spline_vmec, spline_d_vmec, ...
                    phitot_vmec(end));
                
                if output_flag(ii) < 100 % normal solution found
                    disp(['<-recovered #1 after ' num2str(num_attempts) ...
                        ' attempts'])
                    break
                end
                
            end
            if output_flag(ii) < 100 % normal solution found
                %disp('<-recovered #2')
                break
            end
        end
        if output_flag(ii) >= 100 % no solution found
            disp('<---Did not invert correctly. Setting rho==Nan, theta=Nan')
            rho_out(ii) = NaN;
            theta_out(ii) = NaN;
        end
        
    end
    
    
    if output_flag(ii) < 100 % normal solution found
        % find cylindrical components of magnetic field, and jacobian
        [r_out(ii), z_out(ii), r_rho_out(ii), ...
            r_theta_out(ii), r_phi_out(ii), z_rho_out(ii), ...
            z_theta_out(ii), z_phi_out(ii), ...
            lambda_phi_out(ii), lambda_theta_out(ii), ...
            b_r_out(ii), b_z_out(ii), b_tor_out(ii), ...
            sqrt_g_out(ii)] = ...
            get_VMEC_inverse_coords(rho_out(ii), theta_out(ii), lab_Phi(ii), ...
            mnmax_vmec, ...
            xm_vmec, xn_vmec, spline_vmec, spline_d_vmec, phitot_vmec(end));
    else
        % either way, don't bother w/ the derivatives or derived quantites
        % out there.
        r_out(ii) = NaN;
        z_out(ii) = NaN;
        rho_out(ii) = NaN;
        theta_out(ii) = NaN;
    end
    if output_flag(ii) ~= 0
        r_rho_out(ii) = NaN;
        r_theta_out(ii) = NaN;
        r_phi_out(ii) = NaN;
        z_rho_out(ii) = NaN;
        z_theta_out(ii) = NaN;
        z_phi_out(ii) = NaN;
        lambda_phi_out(ii) = NaN;
        lambda_theta_out(ii) = NaN;
        b_r_out(ii) = NaN;
        b_z_out(ii) = NaN;
        b_tor_out(ii) = NaN;
        sqrt_g_out(ii) = NaN;
    end
end

% for <B^2> and <B>
bsq_out = b_r_out.^2 + b_z_out.^2 + b_tor_out.^2;
b1_out = sqrt(bsq_out);

% construct grad Psi in cylindrical coords.
% gradPsi = cross product of co-variant basis vector components,
% e_Phi and e_Theta, divded by the Jacobian
Psi_Phi_out = (r_phi_out .* z_theta_out - r_theta_out .* z_phi_out) ./ ...
    sqrt_g_out;
Psi_R_out = -z_theta_out ./ sqrt_g_out;
Psi_Z_out = r_theta_out ./ sqrt_g_out;

% calculate gradPsi and B in cartesian coords.
b_x_out = -b_tor_out .* sin(lab_Phi) + b_r_out .* cos(lab_Phi);
b_y_out = b_tor_out .* cos(lab_Phi) + b_r_out .* sin(lab_Phi);
% use chain rule for gradPsi
Psi_x_out = (lab_x .* Psi_R_out ./ lab_R) - ...
    (lab_y .* Psi_Phi_out ./ lab_R.^2);
Psi_y_out = (lab_y .* Psi_R_out ./ lab_R) + ...
    (lab_x .* Psi_Phi_out ./ lab_R.^2);


if DO_TIME_TEST
    wall_time = toc;
    disp(['Total time in main loop of get_VMEC_derived_values: ' num2str(wall_time)]);
    disp(['Avg time/point in main loop of get_VMEC_derived_values: ' num2str(wall_time/num_pts)]);
end





