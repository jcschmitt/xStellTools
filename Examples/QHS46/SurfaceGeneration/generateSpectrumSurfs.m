function [] = generateSpectrumSurfs(surfaceToGenerate, DEBUG)
% sample code that uses HSXLineFollow and bs_derics_aux_mex to perform line
% following and generate puncture plot data.  Events are enabled in this
% demo.

if (nargin < 2)
    DEBUG = 0;
end

% This is incremented because the Condor system starts counting from 0, not
% from 1.
surfaceToGenerate = surfaceToGenerate + 1;

%----------------------------------------
% Option settings
%
% All options are read in from an ASCII file stored within the starting
% directory
%----------------------------------------

optionsFile = ['Spectrum_Options.txt'];

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
spectrumType = fgetl(fid_input);

fillerLine = fgetl(fid_input);
fillerLine = fgetl(fid_input);
current = fscanf(fid_input, '%f');

% fillerLine = fgetl(fid_input)
% It seems that a fscacf reads in the next white space, removing the need
% for two of these reads...

fillerLine = fgetl(fid_input);
taper = fscanf(fid_input, '%f');
taper = taper';

fillerLine = fgetl(fid_input);
earthField = fscanf(fid_input, '%i');

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
    filename = ['HSXLineFollow_' spectrumType '_Spectrum_Data_surface_' num2str(ii)];
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
        else
            rStartBlock = coords(end,1);
            zStartBlock = coords(end,3);
            phiStartBlock = coords(end,2);
            chiStartBlock = chi(end);
        end
        chiEndBlock = chiEnd(jj);

        startTime(jj,:) = clock;

        % Aiight.  Let's do the calculation!
        [chiBlock, coordsBlock] = ...
            HSXLineFollow_Spectrum(spectrumType, earthField, current, taper, rStartBlock, phiStartBlock, ...
            zStartBlock, chiStartBlock, chiInc, chiEndBlock, relTol, absTol);

        stopTime(jj,:) = clock;

        if ~(exist('coords'))
            chi = chiBlock;
            coords = coordsBlock;
        else
            chi = [chi ; chiBlock(2:end) ];
            coords = [coords ; coordsBlock(2:end,:) ];
        end

        % Save the data to file
        save(filename, ...
            'magneticConfiguration', ...
            'spectrumType', ...
            'current', ...
            'taper', ...
            'earthField', ...
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
            'coords');
        %             '-v6');
    end  % for jj = startingBlock:transitBlocks

    totalTime = sum(stopTime - startTime);

    if strcmp(lower(spectrumType), 'hamada')
        len_chi = length(chi);
        len_chi_90 = round( 0.9 * len_chi );
        dV_dPsi_all = (chi - chi(1)) ./ ( (coords(:, 2) - coords(1, 2)) / (2*pi) );
        dV_dPsi = mean( dV_dPsi_all(len_chi_90:end) );
        g_Boozer = -1; % The Boozer g factor is not needed for Hamada spectrum calculations
    elseif strcmp(lower(spectrumType), 'boozer')
        dV_dPsi = -1;  % Does it matter, if you're not doing the Hamada spectrum calculation?
        if ~(exist('g_Boozer'))
            g_Boozer = calculateBoozer_gFactor(earthField, current, taper, rStart(1), relTol, absTol);  % rStart(1) should be the magnetic axis
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

        chi_1 = -chi(1:end)';  % the first half of the doubled data array
        chi_1 = fliplr(chi_1);
        r_1 = coords(1:end, 1)';      r_1 = fliplr(r_1);
        phi_1 = coords(1:end, 2)';    phi_1 = fliplr(phi_1);
        z_1 = coords(1:end, 3)';      z_1 = fliplr(z_1);

        chi_2 = chi(2:(end-1))'; % the second half of the doubled data array
        r_2 = coords((2:(end-1)), 1)';
        phi_2 = coords((2:(end-1)), 2)';
        z_2 = coords((2:(end-1)), 3)';

        chi_FFT = [chi_1 chi_2];    % Build the double length arrays
        r_FFT = [r_1 r_2];
        phi_FFT = [phi_1 phi_2];
        z_FFT = [z_1 z_2];

        disp(['Now calculating |B| along the entire chi-curve']);
        % Calculate the values of the |B| for each of the points specified by
        % by *_FFT
        tic
        numPoints_alongChi = length(r_FFT);
        % preallocate memory for B components
        bx = zeros(1,numPoints_alongChi);
        by = bx; bz = bx;
        if 0
            FIELD_SOURCE = 1;
            for kk = 1:numPoints_alongChi
                [bx(kk), by(kk), bz(kk)] = calc_b_bs_JL(r_FFT(kk), phi_FFT(kk), z_FFT(kk), FIELD_SOURCE);
            end
        else
            disp(['total # of points: ' num2str(numPoints_alongChi)])
            %factor = 0.001;
            %for kk = 1:numPoints_alongChi
            %    if (mod(kk,round(factor*numPoints_alongChi)) == 0)
            %        disp(['Finished with ' num2str(100*factor*round((1/factor)*kk/numPoints_alongChi)) '% of the points']);
            %    end
            parfor kk = 1:numPoints_alongChi
                [bx(kk), by(kk), bz(kk)] = calc_b_bs_JL(r_FFT(kk), phi_FFT(kk), z_FFT(kk), current, taper);
            end
        end

        modB = sqrt(bx.^2 + by.^2 + bz.^2);
        toc

        % Save the data to file
        save(filename, ...
            'magneticConfiguration', ...
            'spectrumType', ...
            'current', ...
            'taper', ...
            'earthField', ...
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
            'modB');
        %             '-v6');

        disp('Done calculating |B|');
    else
        disp('|B| already calculated.  Moving on to the spectrum calculation.');
    end  % ~(exist('modB'))
    
    if 1
        % load the iota from the iota/toroidal flux calculation
        [rStart, iota_pp, psi_tor] = getIotaPsi_tor(magneticConfiguration, ii);
        
        % determine the fourier spectrum from modB and chi_FFT
        % %     [nm_amp, pk_n_sorted, pk_m_sorted, nm_avail, nm_error, nm_next_best_error] = calculateSpectrum(chi_FFT, modB, spectrumType, dV_dPsi, g_Boozer, iota_pp);
        [nm_amp, pk_n_sorted, pk_m_sorted, pk_pos_sorted, nm_avail, nm_error, nm_next_best_error, n_values, m_values, iota_best] = ...
            calculateSpectrum(chi_FFT, modB, spectrumType, dV_dPsi, g_Boozer, iota_pp, DEBUG);
        
        % Save the data to file
        save(filename, ...
            'magneticConfiguration', ...
            'spectrumType', ...
            'current', ...
            'taper', ...
            'earthField', ...
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
        
        length_expected_chi = (chiEnd*2^chiDensity) + 1;
        length_chi = length(chi);
        if (length_expected_chi ~= length_chi)
            disp('length of chi does not seem correct.');
        else
            disp('length of chi seems correct.');
        end
    end
end  % for ii = surfaceToGenerate



