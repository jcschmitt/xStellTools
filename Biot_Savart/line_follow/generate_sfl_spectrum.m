function [] = generate_sfl_spectrum(surfaceToGenerate)
% sample code that uses HSXLineFollow and bs_derics_aux_me to perform line
% following and generate puncture plot data.  Events are enabled in this
% demo.
DEBUG=0
% This is incremented because the Condor system starts counting from 0, not
% from 1.
% surfaceToGenerate = surfaceToGenerate + 1;

% determine OS.  If *nix, set the path to include Biot-Savart code and grid
% interpolation paths.  Useful for Condor jobs
% if (isunix)
%     path(path,'~/HSX/');
%end

%----------------------------------------
% Option settings
%

% All options are read in from an ASCII file stored within the starting
% directory
%----------------------------------------

optionsFile = ['SFL_Spectrum_Options.txt'];

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
coilCurrents = fscanf(fid_input, '%f');

fillerLine = fgetl(fid_input);
fillerLine = fgetl(fid_input);
stellaratorSymmetry = fscanf(fid_input, '%f');

fillerLine = fgetl(fid_input);
spectrumType = lower(fgetl(fid_input));

% fillerLine = fgetl(fid_input)
% It seems that a fscacf reads in the next white space, removing the need
% for two of these reads...

fillerLine = fgetl(fid_input);
fillerLine = fgetl(fid_input);
fillerLine = fgetl(fid_input);
relTol = fscanf(fid_input, '%f');

fillerLine = fgetl(fid_input);
absTol = fscanf(fid_input, '%f');

fillerLine = fgetl(fid_input);
numSurfaces = fscanf(fid_input,'%i');

fillerLine = fgetl(fid_input);
chiEnd = fscanf(fid_input,'%f');

fillerLine = fgetl(fid_input);
chiDensity = fscanf(fid_input,'%f');
chiInc = 1/2^chiDensity;

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
    
    chiSkip = 10;  % This sets a length (in chi) the line following code will follow before saving the data to a file.
    % Useful if the process may be interrupted (such as with Condor, or a reboot)
    if (chiEnd > chiSkip)
        chiEnd = (chiSkip:chiSkip:chiEnd);
    end
    transitBlocks = length(chiEnd);
    
    disp(['Generating surface #', num2str(ii)]);
    filename = ['SFL_SpectrumFollow_' spectrumType '_Data_surface_' num2str(ii)];
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
            chiStartBlock = 0;
            rStartBlock2 = rStart(ii);
            zStartBlock2 = zStart(ii);
            phiStartBlock2 = phiStart(ii);
            chiStartBlock2 = 0;
        else
            rStartBlock = coords(end,1);
            zStartBlock = coords(end,3);
            phiStartBlock = coords(end,2);
            chiStartBlock = chi(end);
            rStartBlock2 = coords2(end,1);
            zStartBlock2 = coords2(end,3);
            phiStartBlock2 = coords2(end,2);
            chiStartBlock2 = chi2(end);
        end
        chiEndBlock = chiEnd(jj);
        chiEndBlock2 = chiEnd(jj);
        
        startTime(jj,:) = clock;
        
        % Aiight.  Let's do the calculation!
        [chiBlock, coordsBlock] = ...
            field_line_follow_sfl_spectrum(coilsetID, coilCurrents, ...
            spectrumType, ...
            rStartBlock, phiStartBlock, ...
            zStartBlock, chiStartBlock, chiInc, chiEndBlock, relTol, absTol);
        
        if (stellaratorSymmetry)
            chiBlock2 = chiBlock;
            coordsBlock2 = coordsBlock;
        else
            [chiBlock2, coordsBlock2] = ...
                field_line_follow_sfl_spectrum(coilsetID, -1*coilCurrents, ...
                spectrumType, ...
                rStartBlock2, phiStartBlock2, ...
                zStartBlock2, chiStartBlock2, chiInc, chiEndBlock2, relTol, absTol);
        end
        
        stopTime(jj,:) = clock;
        
        if ~(exist('coords'))
            chi = chiBlock;
            coords = coordsBlock;
            chi2 = chiBlock2;
            coords2 = coordsBlock2;
        else
            chi = [chi ; chiBlock(2:end) ];
            coords = [coords ; coordsBlock(2:end,:) ];
            chi2 = [chi2 ; chiBlock2(2:end) ];
            coords2 = [coords2 ; coordsBlock2(2:end,:) ];
        end
        
        % Save the data to file
        save(filename, ...
            'magneticConfiguration', ...
            'coilsetID', ...
            'coilCurrents', ...
            'stellaratorSymmetry', ...
            'spectrumType', ...
            'relTol', ...
            'absTol', ...
            'numSurfaces', ...
            'startTime', ...
            'stopTime', ...
            'rStart', ...
            'zStart', ...
            'phiStart', ...
            'chiEnd', ...
            'chiDensity', ...
            'surfaceToGenerate', ...
            'currentBlock', ...
            'chi', ...
            'coords', ...
            'chi2', ...
            'coords2');
        %             '-v6');
    end  % for jj = startingBlock:transitBlocks
    
    totalTime = sum(stopTime - startTime);
    if (length(totalTime) == 1)
        totalTime = stopTime - startTime;
    end
    
    if strcmpi(spectrumType, 'hamada')
        disp('<---entering an un-verirfied part of the code!!!!');
        disp('<---entering an un-verirfied part of the code!!!!');
        disp('<---entering an un-verirfied part of the code!!!!');
        len_chi = length(chi);
        len_chi_90 = round( 0.9 * len_chi );
        dV_dPsi_all = (chi - chi(1)) ./ ( (coords(:, 2) - coords(1, 2)) / (2*pi) );
        dV_dPsi = mean( dV_dPsi_all(len_chi_90:end) );
        g_Boozer = -1; % The Boozer g factor is not needed for Hamada spectrum calculations
    elseif strcmpi(spectrumType, 'boozer')
        dV_dPsi = -1;  % Does it matter, if you're not doing the Hamada spectrum calculation?
        if ~(exist('g_Boozer'))
            g_Boozer = calculate_sfl_gFactor(coilsetID, coilCurrents, rStart(1), relTol, absTol);  % rStart(1) should be the magnetic axis
        else
            disp('value of g_Boozer already found--going with it.')
        end
    else error('Unknown spectrum request');
    end
    
    disp(['Times around the machine: ' num2str(( (coords(end, 2) - coords(1, 2)) / (2*pi) )) ]);
    
    disp(['Total time for line following is: ',  num2str(totalTime(3)), ':', num2str(totalTime(4)), ':', ...
        num2str(totalTime(5)), ':', num2str(totalTime(6))]);
    
    
    if ~(exist('modB'))
        disp(['Doubling up the chi and coordinate data.']);
        chi_L = -chi2(1:end)';  % the first half of the doubled data array
        chi_L = fliplr(chi_L);
        r_L = coords2(1:end, 1)';      r_L = fliplr(r_L);
        phi_L = coords2(1:end, 2)';    phi_L = fliplr(phi_L);
        z_L = coords2(1:end, 3)';      z_L = fliplr(z_L);
        
        % skipping the 1st point since it is shared with chi_L
        chi_R = chi(2:(end-1))'; % the second half of the doubled data array
        r_R = coords((2:(end-1)), 1)';
        phi_2 = coords((2:(end-1)), 2)';
        z_R = coords((2:(end-1)), 3)';
        
        chi_FFT = [chi_L chi_R];    % Build the double length arrays
        r_FFT = [r_L r_R];
        phi_FFT = [phi_L phi_2];
        z_FFT = [z_L z_R];
        
        disp(['Now calculating |B| along the entire chi-curve']);
        % Calculate the values of the |B| for each of the points specified by
        % by *_FFT
        tic
        numPoints_alongChi = length(r_FFT);
        % preallocate memory for B components
        bx = zeros(1,numPoints_alongChi);
        by = bx; bz = bx;
        for kk = 1:numPoints_alongChi
            [bx(kk), by(kk), bz(kk)] = calc_b_RPhiZ(coilsetID, r_FFT(kk), phi_FFT(kk), z_FFT(kk), coilCurrents);
        end
        
        modB = sqrt(bx.^2 + by.^2 + bz.^2);
        toc
        
        % Save the data to file
        save(filename, ...
            'magneticConfiguration', ...
            'coilsetID', ...
            'coilCurrents', ...
            'stellaratorSymmetry', ...
            'spectrumType', ...
            'relTol', ...
            'absTol', ...
            'numSurfaces', ...
            'startTime', ...
            'stopTime', ...
            'rStart', ...
            'zStart', ...
            'phiStart', ...
            'chiEnd', ...
            'chiDensity', ...
            'surfaceToGenerate', ...
            'currentBlock', ...
            'chi', ...
            'coords', ...
            'chi2', ...
            'coords2', ...
            'dV_dPsi', ...
            'g_Boozer', ...
            'chi_FFT', ...
            'modB');
        %             '-v6');
        
        disp('Done calculating |B|');
    else
        disp('|B| already calculated.  Moving on to the spectrum calculation.');
    end  % ~(exist('modB'))
    
    % load the iota from the iota/toroidal flux calculation
    [rStart, iota_pp, psi_tor] = get_iota(magneticConfiguration, ii);
    
    % determine the fourier spectrum from modB and chi_FFT
    % %     [nm_amp, pk_n_sorted, pk_m_sorted, nm_avail, nm_error, nm_next_best_error] = calculateSpectrum(chi_FFT, modB, spectrumType, dV_dPsi, g_Boozer, iota_pp);
    if 0
        [nm_amp, pk_n_sorted, pk_m_sorted, pk_pos_sorted, nm_avail, nm_error, nm_next_best_error, n_values, m_values, iota_best] = ...
            calculate_sfl_spectrum2(chi_FFT, modB, spectrumType, dV_dPsi, g_Boozer, iota_pp, DEBUG);
        
        % Save the data to file
        save(filename, ...
            'magneticConfiguration', ...
            'coilsetID', ...
            'coilCurrents', ...
            'stellaratorSymmetry', ...
            'spectrumType', ...
            'relTol', ...
            'absTol', ...
            'numSurfaces', ...
            'startTime', ...
            'stopTime', ...
            'rStart', ...
            'zStart', ...
            'phiStart', ...
            'chiEnd', ...
            'chiDensity', ...
            'surfaceToGenerate', ...
            'currentBlock', ...
            'chi', ...
            'coords', ...
            'dV_dPsi', ...
            'g_Boozer', ...
            'chi_FFT', ...
            'modB', ...
            'nm_amp', ...
            'pk_n_sorted', ...
            'pk_m_sorted', ...
            'pk_pos_sorted', ...
            'nm_avail', ...
            'nm_error', ...
            'nm_next_best_error', ...
            'n_values', ...
            'm_values', ...
            'iota_best');
        %         '-v6');
    end
end  % for ii = surfaceToGenerate



