function [rStart, iota, Psi_toroidal] = get_iota(mode, surfaceIndices, VIEW_PLOTS, LCFS_index)
% function [rStart, iotaProfile, Psi_toroidal] = getIotaPsi_tor(mode, surfaceIndices, VIEW_PLOTS);
% view the results of the calculations of HSC Iota Profile Toroidal FLux

% determine OS.  If *nix, set the path to include Biot-Savart code and grid
% interpolation paths.  Useful for Condor jobs 

%IOTA_DIRECTORY_PATH = 'Y:\CODES\HSX Iota Profile Toroidal Flux'
%IOTA_DIRECTORY_PATH = '/Users/schmittj/Documents/MATLAB/HSX Iota Profile Toroidal Flux'

if (nargin < 4)
    LCFS_index = length(surfaceIndices);
end


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
        disp(['ii:' num2str(ii)]);
        rStart(ii) = filedata.rStart(surfaceIndices(ii));
        iota(ii) = filedata.iota;
        Psi_toroidal(ii) = filedata.flux_startLoc;
        dV_dPsi_toroidal(ii) = (filedata.IntdldB(end)) / filedata.transits;
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
            
            rStart(ii) = NaN;
            iota(ii) = NaN;
            Psi_toroidal(ii) = NaN;
            dV_dPsi_toroidal(ii) = NaN;
        %end
    end
end

%r_eff = sqrt( abs(Psi_toroidal) / (0.5 * pi) );  % Assuming 0.5 Tesla on axis
% r_eff = sqrt( abs(Psi_toroidal) / (1.0 * pi) );  % Assuming 1 Tesla on axis
r_eff = sqrt( abs(Psi_toroidal) / (2.5 * pi) );  % Assuming 2.5 Tesla on axis

try
Volume = cumtrapz(Psi_toroidal, dV_dPsi_toroidal);
catch
Volume = 0;
end
%Volume = cumtrapz(Psi_toroidal, dV_dPsi_toroidal);

if VIEW_PLOTS
    h = figure;
    subplot(3,1,1);
    plot(surfaceIndices, iota, '+');
    hold on
    %plot(surfaceIndices, -iota, 'r+');
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
    %plot(surfaceIndices, r_eff/r_eff(end), '+');
    plot(surfaceIndices, r_eff/r_eff(LCFS_index), '+');
    ylabel('\rho');
    xlabel('Index #');
    title(mode);
    subplot(2,1,2);
    plot(r_eff/r_eff(end), Psi_toroidal, '+');
    ylabel('\psi_{tor}');
    xlabel('\rho');

    figure;
    plot(r_eff/r_eff(end), iota, '+');
    %plot(r_eff/r_eff(end), iota, '--');
    %ylabel('Rotational Transform');
    ylabel('\iota / 2\pi');
    xlabel('\rho');
    title(mode);

    figure;
    plot(r_eff/r_eff(end), dV_dPsi_toroidal, '+');
    %plot(r_eff/r_eff(end), iota, '--');
    %ylabel('Rotational Transform');
    ylabel('dV / d\psi');
    xlabel('\rho');
    title(mode);

    figure;
    plot(r_eff/r_eff(end), 100* (dV_dPsi_toroidal(1) - dV_dPsi_toroidal)/dV_dPsi_toroidal(1), '+');
    %plot(r_eff/r_eff(end), iota, '--');
    %ylabel('Rotational Transform');
    ylabel('Well Depth, %');
    xlabel('\rho');
    title(mode);


    figure;
    plot(r_eff/r_eff(end), Volume, '+');
    %plot(r_eff/r_eff(end), iota, '--');
    %ylabel('Rotational Transform');
    ylabel('Volume');
    xlabel('\rho');
    title(mode);


end
    
    
