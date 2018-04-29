function [] = displaySpectrum(folder, spectrumType, surfaceNumber, maxMode, minMode, recon_B, figNum);
%
if (nargin < 6)
    recon_B = 0
end
if (nargin < 7)
    figNum = 0;
end
if (nargin < 5)
    minMode = Inf;
end
if (nargin < 4)
    maxMode = 1;
end


disp(['Displaying ' spectrumType ' spectrum of surface #' num2str(surfaceNumber) ' of ' folder]);
filename = ['LineFollow_' spectrumType '_Spectrum_Data_surface_' num2str(surfaceNumber)];
filenameComplete = [folder '/' filename];

load(filenameComplete, ...
    'chi', ...
    'coords', ...
    'modB', ...
    'spectrumType', ...
    'dV_dPsi', ...
    'g_Boozer', ...
    'n_values', ...
    'm_values', ...
    'pk_pos_sorted', ...
    'iota_best', ...
    'nm_amp', ...
    'pk_n_sorted', ...
    'pk_m_sorted', ...
    'nm_avail', ...
    'nm_error', ...
    'nm_next_best_error', ...
    'n_values', ...
    'm_values');

if recon_B
    reconstructModB(chi, modB, spectrumType, dV_dPsi, g_Boozer, n_values, m_values, ...
        pk_n_sorted, pk_m_sorted, pk_pos_sorted, nm_amp, iota_best, 1);
    % -----------------------------------
    % the following section is broken!!!
    % -----------------------------------
    % % % %     % Reorder the r, phi, z arrays into matrices
    % % % %     % Reorder the modB array into a matrix
        r = coords(:,1); phi = coords(:,2); z = coords(:,3);
        BS_B_x = r.*cos(phi);
        BS_B_y = r.*sin(phi);
        BS_B_z = z;
        modB_plot = modB(1:(length(modB)/2 + 1));  % grab first part of the modB contour-the part that lines up w/ these coords
    % % % %     figure;
    % % % % %     surface(BS_B_x, BS_B_y, BS_B_z,  'linestyle', 'none');
end

% % zero out the (0,0)) and (4,1) modes
% m_zero = find(m_values == 0);
% n_zero = find(n_values == 0);
% nm_amp(n_zero, m_zero) = 0;
% m_one = find(m_values == 1);
% n_four = find(n_values == 4);
% nm_amp(n_four, m_one) = 0;

% zero out the modes not between maxMode and minMode
minMode = min(minMode, length(pk_pos_sorted));
disp(['Displaying largest modes: ' num2str(maxMode) ' through ' num2str(minMode)]);
allModes = 1:length(pk_pos_sorted);
modesToZero = setdiff(allModes,maxMode:minMode);
for ii = modesToZero
    %     nm_amp(pk_n_sorted(ii), pk_m_sorted(ii));
    nm_amp(pk_n_sorted(ii), pk_m_sorted(ii)) = 0;
end    

% if recon_B
%     reconstructModB(chi, spectrumType, dV_dPsi, g_Boozer, n_values, m_values, ...
%         pk_n_sorted, pk_m_sorted, pk_pos_sorted, nm_amp, iota_best, 1);
% end

% re-do the arrays to help make the plots look "neat-o!"
expandFactor = 10;
delta_n = (n_values(2) - n_values(1));
delta_m = (m_values(2) - m_values(1));
last_n = n_values(end) + delta_n * (expandFactor - 1) / expandFactor;
last_m = m_values(end) + delta_m * (expandFactor - 1) / expandFactor;
exp_n_values = n_values(1):(delta_n/expandFactor):(last_n);
exp_m_values = m_values(1):(delta_m/expandFactor):(last_m);
for ii = ( [1:length(n_values)] - 1)
    for jj = ( [1:length(m_values)] - 1)
        exp_nm_amp( expandFactor*ii + [1:expandFactor], ...
            expandFactor*jj + [1:expandFactor] ) = nm_amp(ii+1, jj+1);
    end
end

if (figNum < 1)
    figure
else
    figure(figNum);
end

surf(exp_m_values, exp_n_values, abs(exp_nm_amp) )
view(2)
colormap('autumn')
colorbar
shading flat
xlabel('m numbers');
ylabel('n numbers');
title(['Surface #' num2str(surfaceNumber)]);
%  axis([-10 10 0 80])