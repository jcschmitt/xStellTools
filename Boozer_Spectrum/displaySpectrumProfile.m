function [] = displaySpectrumProfile(folder, text_title, spectrumType, ...
    surfs, ModesToDisplay, COLOR_SPEC, STYLE_SPEC, XAXIS, fig_handle_in)
% function [] = displaySpectrumProfile(folder, text_title, spectrumType, ...
%     surfs, ModesToDisplay, COLOR_SPEC, STYLE_SPEC, XAXIS, fig_handle_in)
% folder: absolute or relative folderpath for the data files
% text_title: The title string that will be added to the figure.
% spectrumType: 'boozer' or 'hamada'  (only boozer has been tested so far)
% surfs: the indices of the surfaces that are to be ploted.
% the following are optional:
%   ModesToDisplay: the indices of the modes to display (sorted by magnitude)
%   COLOR_SPEC: the color of the lines/symbols
%   STYLE_SPEC: the style of the plotted data
%   XAXIS: The x-axis to be used, 
%   fig_handle_in: For overplotting on the figure #(fig_handle_in);
<<<<<<< HEAD
%  displaySpectrumProfile('QHS', 'QHS', 'boozer', 1:26, [9 11 15 18 26 27
%           29 40 62 74 80 91  119], 'k', '+', [1:26]/29) 
=======

>>>>>>> de3a1a81b199526aa12d96643625b3a77d223cc6

if nargin < 9
    fig_handle_in = figure;
end

if nargin < 8
    XAXIS = surfs;
elseif length(XAXIS) ~= length(surfs)
    error('<----Number of surfaces does not match specs of XAXIS');
end
if nargin < 7
    STYLE_SPEC = '-';
end

if nargin < 6
    COLOR_SPEC = 'b';
end
if (nargin < 5)
    ModesToDisplay = 1:10;
end
if (nargin < 4)
    error('<----Not enough inputs');
end

SHOW_LABELS = 1;

disp(['<----Displaying radial profile of ' spectrumType ' spectrum in folder: ' folder]);

numSurfs = length(surfs);

for ii = 1:numSurfs
    filename = ['SFL_SpectrumFollow_' spectrumType '_Data_surface_' num2str(surfs(ii)) '.mat'];
    filenameComplete = [folder '/' filename];

    load(filenameComplete, ...
        'spectrumType', ...
        'dV_dPsi', ...
        'g_Boozer', ...
        'pk_pos_sorted', ...
        'iota_best', ...
        'nm_amp', ...
        'nm_avail', ...
        'n_values', ...
        'm_values');
     nm_amp_all{ii}.nm_amp = nm_amp;
    %nm_amp_all{ii}.nm_amp = nm_amp*2;  % for comparing 1/2 T to 1 T calculations
    %disp('Mysterious factor of 2..Are you sure???')
    iota_best_all(ii) = iota_best;
    disp(['<----File # ' num2str(surfs(ii)) ' loaded']);
end

% Re-sort the amplitude data into nice arrays
num_n = length(n_values);
num_m = length(m_values);
for ii = 1:num_n
    for jj = 1:num_m
        for kk = 1:numSurfs
            spectrumAmp_versus_surfNum( (ii - 1) * num_m + jj, kk ) = nm_amp_all{kk}.nm_amp(ii, jj);

        end
        mode_label{ (ii - 1) * num_m + jj} = [ '(' num2str(n_values(ii)) ',' num2str(m_values(jj)) ')' ];
    end
end

% all mode except (0,0) should be x2
% spectrumAmp_versus_surfNum = 4* spectrumAmp_versus_surfNum;  % for when you want 1T spectrum from 0.5 T line follow data
%spectrumAmp_versus_surfNum = 2* spectrumAmp_versus_surfNum; % for normal
% display

%m_zero = find(m_values == 0);
%n_zero = find(n_values == 0);
%spectrumAmp_versus_surfNum( (n_zero-1) * num_m + m_zero, :) = ...
%    0.5 * spectrumAmp_versus_surfNum( (n_zero-1) * num_m + m_zero, :);


for ii = 1:num_n*num_m
    max_of_mode(ii) = max(abs(spectrumAmp_versus_surfNum(ii,:)));
end

[~, sorted_indices] = sort(max_of_mode, 'descend');

figure(fig_handle_in)
legend_text = [];
for ii = ModesToDisplay
    jj = mod(ii-1, length(COLOR_SPEC)) + 1;
    plot(XAXIS, spectrumAmp_versus_surfNum(sorted_indices(ii), :), [COLOR_SPEC(jj) STYLE_SPEC(jj) '']);
    hold on
    if SHOW_LABELS
        %text( XAXIS(end)+.01*(XAXIS(end) - XAXIS(1)), spectrumAmp_versus_surfNum(sorted_indices(ii), end), ...
        %    mode_label(sorted_indices(ii)), 'Color', COLOR_SPEC(jj));
        [XAXIS_max, ind_max] = max(XAXIS);
        [XAXIS_min, ind_min] = min(XAXIS);
        text( XAXIS_max+.01*(XAXIS_max - XAXIS_min), spectrumAmp_versus_surfNum(sorted_indices(ii), ind_max), ...
            mode_label(sorted_indices(ii)), 'Color', COLOR_SPEC(jj));
    end
    legend_text = [legend_text ; mode_label(sorted_indices(ii))];
end
xlabel('Surf #');
ylabel('Spectrum amplitude');
cur_axis = axis;
axis([cur_axis(1) cur_axis(2)*1.015 cur_axis(3) cur_axis(4)]);
title(text_title);
%title([folder ' ' spectrumType]);
legend(legend_text);
