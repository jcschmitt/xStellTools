function [] = check_descur_output(lcfs_file, outcurve_file, num_field_periods)
% A slight modificaiton of a file written and modified by J. Lore. (HSX)
% If the input files are not entered, the user will be queried.

if nargin < 1
    %get the R,Z coords from field line following
    [filename,pathname] = uigetfile('*.mat', ...
        'Indicate the file containing the x,y,z coordinates of the LCFS');
    load([pathname '/' filename], 'surf_x', 'surf_y', 'surf_z');
else
    load(lcfs_file, 'surf_x', 'surf_y', 'surf_z');
end
if nargin < 3
    num_field_periods = input('----> How many field periods? ');
end

rsurf = sqrt(surf_x.^2 + surf_y.^2);
zsurf = surf_z;

% Try to determine ntheta and nphi from the file
[theta, ~, ~] = cart2pol(surf_x, surf_y, surf_z);
dtheta = theta(2) - theta(1);

% 
numphi = round((2 * pi / num_field_periods) / dtheta);
ntheta = round((length(theta) - 1) / numphi);

disp(['<----Found file. Looks like nphi=' num2str(numphi) ...
    ' and ntheta=' num2str(ntheta)]);

% This is probalby related to numphi, but I haven't checked. I know this
% works for the current version of DESCUR and 50 slices/toroidal field
% period. In the end, I don't think it matters because all modes are
% eventually read in.
%numtormodes_first_pass = 25;
numtormodes_first_pass = 30;

%phi angle of each cut
%phi = 0:2*pi/(4*numphi):pi/2;
phi = [0:numphi] *  ( 2*pi / ( num_field_periods ) ) / numphi;
theta = linspace(0,2*pi,1000);

if nargin < 2
    [fid_outcurve]=fopen(uigetfile('*.*', ...
        'Indicate the descur outcurve file'),'r');
else
    fid_outcurve = fopen(outcurve_file);
end

%parse the file to find the RBC lines
try
    tline = fgetl(fid_outcurve)
    if length(tline) < 5
        cmpSTR = ' ';
    else
        cmpSTR = tline(4:5);
    end
    while((~strcmp(cmpSTR,'MB')))
        tline = fgetl(fid_outcurve);
        if length(tline) < 5
            cmpSTR = ' ';
        else
            cmpSTR = tline(4:5);
        end
        if (tline == -1)
            break;
            disp('<----');
            
        end
    end
    
    for ii = 1:numtormodes_first_pass
        temp_data = fgetl(fid_outcurve);
        data(ii,:) = temp_data(1:57);
    end
    
    keep_going = 1;
    counter = ii;
    while keep_going % check the outcurve file, this ??? changes.
        temp_data = fgetl(fid_outcurve)
        if max(size(temp_data)) > 1
            counter = counter + 1;
            data(counter,:) = temp_data;
        else
            keep_going = 0;
        end
    end
    fclose(fid_outcurve);

    % Convert all those strings into numbers at one time
    data = str2num(data);
    numtotalmodes = counter;
    
    disp(['<----Total number of modes found: ' num2str(numtotalmodes)]);
    
    % Pull out values
    m = data(:,1);
    n = data(:,2);
    rbc = data(:,3); % cos component of R
    rbs = data(:,4); % sin component of R
    zbc = data(:,5); % cos component of Z
    zbs = data(:,6); % sin component of Z
    
    ERR_OUTCURVE = false;
catch
    ERR_OUTCURVE = true;
    warning('Something messed up when trying to read outcurve file')
end

for kk = 1:length(phi) % for each cut
% for kk = 1 % just the 1st toroidal cut
    r_rec = zeros(size(theta));
    z_rec = zeros(size(theta));
    
    if ~ERR_OUTCURVE
        for ii = 1:length(theta) %sum up the fourier components of R,Z for each m,n for all angles theta
            r_rec(ii) = sum(rbc.*cos(m*theta(ii) - num_field_periods*n*phi(kk)))+sum(rbs.*sin(m*theta(ii) - num_field_periods*n*phi(kk)));
            z_rec(ii) = sum(zbc.*cos(m*theta(ii) - num_field_periods*n*phi(kk)))+sum(zbs.*sin(m*theta(ii) - num_field_periods*n*phi(kk)));
        end
    end
    
    figure(kk);box on; hold on;
    % plot field line follow points
    plot(rsurf(kk:numphi:end), zsurf(kk:numphi:end), 'k.')
    hold on
    % plot fourier points
    if ~ERR_OUTCURVE
        plot(r_rec,z_rec,'r-','LineWidth',1)
    end
    axis equal
end

