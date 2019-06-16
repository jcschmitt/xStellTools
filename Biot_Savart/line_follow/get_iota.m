function [rStart, iota, Psi_toroidal] = getIotaPsi_tor(mode, surfaceIndices, VIEW_PLOTS)
% function [rStart, iotaProfile, Psi_toroidal] = getIotaPsi_tor(mode, surfaceIndices, VIEW_PLOTS);
% view the results of the calculations of HSC Iota Profile Toroidal FLux

% determine OS.  If *nix, set the path to include Biot-Savart code and grid
% interpolation paths.  Useful for Condor jobs 

%IOTA_DIRECTORY_PATH = 'Y:\CODES\HSX Iota Profile Toroidal Flux'
%IOTA_DIRECTORY_PATH = '/Users/schmittj/Documents/MATLAB/HSX Iota Profile Toroidal Flux'


if (nargin < 3)
    VIEW_PLOTS = 0;
end

%foldername = [IOTA_DIRECTORY_PATH '/' mode];
foldername = pwd; 

surfaceIndices = sort(surfaceIndices);
numFiles = length(surfaceIndices);
for ii = 1:numFiles
    filenames{ii} = ['LineFollow_Coord_Data_surface_' num2str(surfaceIndices(ii))];
end

% loop through each surface (one surface per file) and load the radial start location, iota, and toroidal flux
% from each file
for ii = 1:numFiles
    try
        filedata = load([foldername '/' filenames{ii}]);
        disp(['ii:' num2str(ii)])
        rStart(ii) = filedata.rStart(surfaceIndices(ii));
        iota(ii) = filedata.iota;
        Psi_toroidal(ii) = filedata.flux_startLoc;
    catch
        warning(['Unable to load data for ' filenames{ii}]);
        %try
            %disp('Gonna try archive directory')
            %IOTA_DIRECTORY_PATH2 = ['\\Hsx4f_backup\Backup_1\JohnS\HSX Iota Profile Toroidal Flux archive\' mode]
            %filedata = load([IOTA_DIRECTORY_PATH2 '/' filenames{ii}]);
            %disp(['ii:' num2str(ii)])
            %rStart(ii) = filedata.rStart(surfaceIndices(ii));
            %iota(ii) = filedata.iota;
            %Psi_toroidal(ii) = filedata.flux_startLoc;
        %catch
            %disp('It bombed');
            
            rStart(ii) = 0;
            iota(ii) = 0.99;
            Psi_toroidal = .01;
        %end
    end
end

%r_eff = sqrt( abs(Psi_toroidal) / (0.5 * pi) );  % Assuming 0.5 Tesla on axis
% r_eff = sqrt( abs(Psi_toroidal) / (1.0 * pi) );  % Assuming 1 Tesla on axis
r_eff = sqrt( abs(Psi_toroidal) / (2.5 * pi) );  % Assuming 2.5 Tesla on axis

if VIEW_PLOTS
    h = figure;
    subplot(3,1,1);
    plot(surfaceIndices, iota, '+');
    hold on
    plot(surfaceIndices, -iota, 'r+');
    ylabel('\iota/2\pi');
    title(mode);
    subplot(3,1,2);
    plot(surfaceIndices, Psi_toroidal, '+');
    ylabel('\psi_{tor}');
    subplot(3,1,3);
    plot(surfaceIndices, r_eff, '+');
    ylabel('r_{eff}');
    h = figure;
    subplot(2,1,1);
    plot(r_eff, iota, '+');
    ylabel('\iota/2\pi');
    title(mode);
    subplot(2,1,2);
    plot(r_eff, Psi_toroidal, '+');
    ylabel('\psi_{tor}');
    xlabel('r_{eff}');
    h = figure;
    subplot(2,1,1);
    plot(surfaceIndices, r_eff/r_eff(end), '+');
    ylabel('\rho');
    xlabel('Index #');
    title(mode);
    subplot(2,1,2);
    plot(r_eff/r_eff(end), Psi_toroidal, '+');
    ylabel('\psi_{tor}');
    xlabel('\rho');

    figure;
    plot(r_eff/r_eff(end), iota, '+');
    ylabel('Rotational Transform');
    xlabel('\rho');
    title(mode);

end
    
    
