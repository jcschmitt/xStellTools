function [] = make_descur_file_W7X(LCFS_filename, LCFS_pathname, N_phi, N_theta)
% function [] = make_descur_file_W7X(LCFS_filename, LCFS_pathname)
% This file takes converts the data points from a 'calc_descur_RZ_points'
% run into a format that is expected by DESCUR.
%
% If no input arguments are included, the program will query the user to
% select a file.
%

try
    if nargin < 1
        %get the cart. coords. of the lcfs (from make_lcfs_jl.m)
        [LCFS_filename, LCFS_pathname] = uigetfile('*.mat', ...
            'Indicate the file containing the x,y,z coordinates of the LCFS');
    elseif nargin < 2
        if isunix
            LCFS_pathname = './';
        elseif iswindowx
            LCFS_pathname = '.\';
        else
            LCFS_pathname = './';  % Taking a guess here
        end
    elseif isempty(LCFS_pathname);
        if isunix
            LCFS_pathname = './';
        elseif iswindowx
            LCFS_pathname = '.\';
        else
            LCFS_pathname = './';  % Taking a guess here
        end
    end
    
    cur_pwd = pwd;
    cd(LCFS_pathname);
    load(LCFS_filename, ...
        'surf_x', ...
        'surf_y', ...
        'surf_z');
    outfilename = ['rzdata_' LCFS_filename(6:(end-4)) '.dat'];
    
    
    if nargin < 4
        % Try to determine ntheta and nphi from the file
        [theta, ~, ~] = cart2pol(surf_x, surf_y, surf_z);
        dtheta = theta(2) - theta(1);
        
        % W7X has four field periods
        NFP = 5;
        nphi = round((2 * pi / NFP) / dtheta);
        ntheta = round((length(theta) - 1) / nphi);
        
        disp(['<----Found file. Looks like nphi=' num2str(nphi) ...
            ' and ntheta=' num2str(ntheta)]);
    else
        nphi = N_phi;
        ntheta = N_theta;
        disp(['<----Using inputs. nphi=' num2str(nphi) ...
            ' and ntheta=' num2str(ntheta)]);
    end
    
    
    %open output file for writing
    fid = fopen(outfilename,'w');
    fprintf(fid,'%s %s %s \n',num2str(ntheta),num2str(nphi),'4');
    % R = sqrt(surf_x.^2 + surf_y.^2);
    % PHI = atan2(surf_y,surf_x);
    % Z = surf_z;
    
    for phi_ind = 1:nphi %loop over cuts and calc R,phi,Z
        Rsurf = sqrt(surf_x(phi_ind:nphi:end).^2 + surf_y(phi_ind:nphi:end).^2);
        PHI = atan2(surf_y(phi_ind:nphi:end),surf_x(phi_ind:nphi:end));
        Zsurf = surf_z(phi_ind:nphi:end);
        figure(400+phi_ind);hold on;box on;hold on;
        plot(Rsurf, Zsurf, 'k.');
    end
    
    for phi_ind = 1:nphi %loop over cuts and calc R,phi,Z
        Rsurf = sqrt(surf_x(phi_ind:nphi:end).^2 + surf_y(phi_ind:nphi:end).^2);
        PHI = atan2(surf_y(phi_ind:nphi:end),surf_x(phi_ind:nphi:end));
        Zsurf = surf_z(phi_ind:nphi:end);
        for theta_ind = 1:ntheta %print each point on cut
            fprintf(fid,'%f %f %f \n',Rsurf(theta_ind),PHI(theta_ind),Zsurf(theta_ind));
        end
    end
        
    fclose(fid);
    cd(cur_pwd);
catch
    warning('error in make_descur_file_W7X.m');
    LCFS_filename
    LCFS_pathname
    cd(cur_pwd);
    fclose(fid);
end