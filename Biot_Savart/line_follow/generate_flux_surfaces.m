function [] = generate_flux_surfaces(surfaceToGenerate)
% sample code that uses line following and Biot-Savart to perform line
% following and generate puncture plot data.


% To do: Add a check to see if it is the Condor system.  If so, 
% then this is incremented. The Condor system starts counting from 0, not
% from 1.
%surfaceToGenerate = surfaceToGenerate + 1;

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
coilsetID = fgetl(fid_input);

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
rAxis = fscanf(fid_input, '%f');

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
    
    transitSkip = 10;  % This sets how many transits the line following code will follow before saving the data to a file.
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
            field_line_follow(coilsetID, current, rStartBlock, zStartBlock, phiStartBlock, ...
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
            'coilsetID', ...
            'current', ...
            'relTol', ...
            'absTol', ...
            'numSurfaces', ...
            'transits', ...
            'phiIncInDegrees', ...
            'startTime', ...
            'stopTime', ...
            'rAxis', ...
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
    
end

% Calculate the rotational transform.
disp(['Now calculating iota based on puncture plot']);
% Use the puncture plot of the surface at the starting location (usually the phi=0 plane)
phiIndexInc = round(360 / phiIncInDegrees);

iota = calculate_iota(rAxis, coords(:, 1:2), phiIncInDegrees);

% 	break; % break out of this for testing purposes
disp(['Now calculating the enclosed flux based on puncture plot ']);
%if ~(exist('flux_startLoc'))
    
    %flux_startLoc = calculateFlux(current, phi(1:phiIndexInc:end), coords(1:phiIndexInc:end,1:2));
    %flux_startLoc = calculateFlux_W7X(current, phi(1:phiIndexInc:end), coords(1:phiIndexInc:end,1:2), 40, 1);
    %flux_startLoc = -1;
    flux_startLoc = calculate_toroidal_flux(coilsetID, current,  coords(1:phiIndexInc:end,1),  coords(1:phiIndexInc:end,2), 1);
%else
%    disp('<----Flux found. Skipping the flux calculation.')
%end
IntdldB = coords(:, 3); % The line-integral of (dPhi r/Bphi)
coordsRZ = coords(:, 1:2);
save(filename, ...
    'magneticConfiguration', ...
    'coilsetID', ...
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
    'rAxis', ...
    'rStart', ...
    'zStart', ...
    'phiStart', ...
    'surfaceToGenerate', ...
    'currentBlock', ...
    'phi', ...
    'coords', ...
    'coordsRZ', ...
    'IntdldB');
%             '-v6');
end