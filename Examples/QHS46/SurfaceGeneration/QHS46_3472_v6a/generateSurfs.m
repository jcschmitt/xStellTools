function [] = generateSurfs(surfaceToGenerate)
% sample code that uses line following and Biot-Savart to perform line
% following and generate puncture plot data.  


% This is incremented because the Condor system starts counting from 0, not
% from 1.
surfaceToGenerate = surfaceToGenerate + 1;

% determine OS.  If *nix, set the path to include Biot-Savart code and grid
% interpolation paths--especially useful for Condor jobs
% if (isunix)
%      path(path,'~/HSX/bs_mex');
%      path(path,'~/HSX/grid_interp');
%      disp('Paths updated to include bs_mex and grid_interp');
% end


%----------------------------------------
% Option settings
%
% All options are read in from an ASCII file stored within the starting
% directory
%----------------------------------------

optionsFile = ['Line_Follow_Options.txt'];

try
    fid_input=fopen(optionsFile,'r');
catch
    disp(['Error opening options file: ', optionsFile]);
    return;
end

if (fid_input == -1)
    disp(['Error opening options file: ', optionsFile]);
    return
end
    
disp(['Loading options from ' optionsFile]);

fillerLine = fgetl(fid_input);
magneticConfiguration = fgetl(fid_input);

fillerLine = fgetl(fid_input);
fillerLine = fgetl(fid_input);
current = fscanf(fid_input, '%f');

% fillerLine = fgetl(fid_input)
% It seems that a fscacf reads in the next white space, removing the need
% for two of these reads...

%fillerLine = fgetl(fid_input);
%taper = fscanf(fid_input, '%f');
%taper = taper';

%fillerLine = fgetl(fid_input);
%earthField = fscanf(fid_input, '%i');

fillerLine = fgetl(fid_input);
fillerLine = fgetl(fid_input);
relTol = fscanf(fid_input, '%f');

fillerLine = fgetl(fid_input);
absTol = fscanf(fid_input, '%f');

fillerLine = fgetl(fid_input);
numSurfaces = fscanf(fid_input,'%i');

fillerLine = fgetl(fid_input);
transits = fscanf(fid_input,'%f');

fillerLine = fgetl(fid_input);
phiIncInDegrees = fscanf(fid_input,'%f');
phiInc = phiIncInDegrees * pi / 180;

fillerLine = fgetl(fid_input);
fillerLine = fgetl(fid_input);
rStart = fscanf(fid_input,'%f');

fillerLine = fgetl(fid_input);
zStart = fscanf(fid_input,'%f');

fillerLine = fgetl(fid_input);
phiStart = fscanf(fid_input,'%f');

fclose(fid_input);

%---------------------
% Options done
%---------------------


for ii = surfaceToGenerate
    
    transitSkip = inf;  % This sets how many transits the line following code will follow before saving the data to a file. 
                       % Useful if the process may be interrupted (such as with Condor, or a reboot)
    if (transits > transitSkip)
        phiEnd = phiStart(ii) + (transitSkip:transitSkip:transits) * 2 * pi;
    else
        phiEnd = phiStart(ii) + transits * 2 * pi;
    end
    transitBlocks = length(phiEnd);
    
    disp(['Generating surface #', num2str(ii)]);
    filename = ['LineFollow_Coord_Data_surface_' num2str(ii)];
    try
        load(filename);
        disp('Found previous data.')
        startingBlock = currentBlock + 1;
    catch
        disp('Found no previous data.  Initialized surface data');
        startingBlock = 1;
    end
    
    for jj = startingBlock:transitBlocks
        currentBlock = jj;
        disp(['Starting transit block ', num2str(currentBlock), ' of ', num2str(transitBlocks)]);
        
        % Set up the line following parameters (starting, ending location) for
        % current transit block
        if (jj == 1)
            rStartBlock = rStart(ii);
            zStartBlock = zStart(ii);
            phiStartBlock = phiStart(ii);
        else
            rStartBlock = coords(end,1);
            zStartBlock = coords(end,2);
            phiStartBlock = phi(end);
        end
        phiEndBlock = phiEnd(jj);
        
        startTime(jj,:) = clock;
        
        % Events are enabled and recorded
        [phiBlock, coordsBlock] = ... 
            LineFollow_QHS46_34742_v6a(current, rStartBlock, zStartBlock, phiStartBlock, ...
            phiInc, phiEndBlock, relTol, absTol);
         
        stopTime(jj,:) = clock;
        
        if ~(exist('coords'))
            phi = phiBlock;
            coords = coordsBlock;
        else
            phi = [phi ; phiBlock(2:end) ];
            % fix the dl/dB integral
            coordsBlock(:,3) = coordsBlock(:,3) + coords(end,3);
            coords = [coords ; coordsBlock(2:end,:) ];
        end
        
        % Save the data to file
        save(filename, ...
            'magneticConfiguration', ...
            'current', ...
            'relTol', ...
            'absTol', ...
            'numSurfaces', ...
            'transits', ...
            'phiIncInDegrees', ...
            'startTime', ...
            'stopTime', ...
            'rStart', ...
            'zStart', ...
            'phiStart', ...
            'surfaceToGenerate', ...
            'currentBlock', ...
            'phi', ...
            'coords');
%             '-v6');
        
    end
    if (jj > 1)
        totalTime = sum(stopTime - startTime);    
    else
        totalTime = (stopTime - startTime);    
    end
        
    
    disp(['Total time for line following is: ',  num2str(totalTime(3)), ':', num2str(totalTime(4)), ':', ...
            num2str(totalTime(5)), ':', num2str(totalTime(6))]);
    
    disp(['Now calculating iota based on puncture plot']);
    % Use the puncture plot of the surface at the starting location (usually a boxport)
    phiIndexInc = round(360 / phiIncInDegrees);
    end
        
    % if this is the first surface, assume it it close to the axis and determine the magnetic axis from
    % the surface.
    if (ii == 1)
        rAxis = mean( coords(1:phiIndexInc:end, 1) );
    else  % otherwise, just use the first value of the r-coords array and assume it was the axis
        rAxis = rStart(1);
    end
    % Calculate the rotational transform.
    iota = calculateIota(rAxis, coords(1:phiIndexInc:end, 1:2), phiIncInDegrees);
    
    % 	break; % break out of this for testing purposes
    if ~exist('flux_startLoc')
        disp(['Now calculating the enclosed flux based on puncture plot ']);
        %flux_startLoc = calculateFlux(current, phi(1:phiIndexInc:end), coords(1:phiIndexInc:end,1:2));
        flux_startLoc = calculateFlux(current, phi(1:phiIndexInc:end), coords(1:phiIndexInc:end,1:2), 40, 1);
        IntdldB = coords(:, 3);
        coords = coords(:, 1:2);
    else
	    disp('<----Found previous enclosed flux.');
    end
    
    save(filename, ...
        'magneticConfiguration', ...
        'current', ...
        'relTol', ...
        'absTol', ...
        'numSurfaces', ...
        'transits', ...
        'iota', ...
        'flux_startLoc', ...
        'phiIncInDegrees', ...
        'startTime', ...
        'stopTime', ...
        'totalTime', ...
        'rStart', ...
        'zStart', ...
        'phiStart', ...
        'surfaceToGenerate', ...
        'currentBlock', ...
        'phi', ...
        'coords', ...
        'IntdldB');
%             '-v6');
end
