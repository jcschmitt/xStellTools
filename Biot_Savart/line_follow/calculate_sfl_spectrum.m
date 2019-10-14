function [nm_amp, pk_n_sorted, pk_m_sorted, pk_pos_sorted, nm_avail, nm_error, nm_next_best_error, n_values, m_values, iota_best] = ...
    calculateSpectrum(chi, modB, spectrumType, dV_dPsi, g_Boozer, iota_pp, DEBUG);
%               Apply a Gaussian window function to the data
%               Rearrange the data to guarantee that the modB curve has even symmetry
%               (so that the FFT has only real components)
%               FFT the curve
%               find the peaks of the FFT
%               determine best value for iota
%               go find n, m
%               done!

if (nargin < 7)
    DEBUG = 0;
end
PAUSE_TIME_FACTOR = 0.00;  % 0 is fastest, 1 is fast, 100 is sloow.
% make the gaussian window to smooth the ends of the B-curve.
eta = 4;  % Sets the width of the applied window
chi_length = max(chi);  % the length of the chi (in each direction, almost)
% gauss_window = (  eta / ( chi_length * sqrt(2*pi) ) ) * exp( - eta^2 * chi.^2 / (2 * chi_length^2) );
gauss_window = (  eta / ( 1 * sqrt(2*pi) ) ) * exp( - eta^2 * chi.^2 / (2 * chi_length^2) );
%               Apply a Gaussian window function to the data
modB_window = modB .* gauss_window;

%               Rearrange the data to guarantee that the modB curve has even symmetry
%               (so that the FFT has only real components)
midpoint = length(modB_window) / 2;
modB_even = [modB_window((midpoint+1):end) modB_window(1:midpoint)];
modB_orig_even = [modB((midpoint+1):end) modB(1:midpoint)];
gw_even = [gauss_window((midpoint+1):end) gauss_window(1:midpoint)];
%               FFT the curve
B_FFT_complex = eta*fft(modB_even) / length(modB_even); % fft w normalization
B_FFT_orig_complex = fft(modB_orig_even) / length(modB_orig_even); % fft w/o normalization
B_FFT_gw_complex = fft(gw_even) / length(gw_even); % fft w/o normalization
B_FFT_real = real(B_FFT_complex);  % the imaginary portion is pretty small-plot it and check to make sure, tho.
B_FFT_orig_real = abs(B_FFT_orig_complex);  % the imaginary portion is pretty small-plot it and check to make sure, tho.
B_FFT_gw_real = real(B_FFT_gw_complex);  % the imaginary portion is pretty small-plot it and check to make sure, tho.
B_FFT = B_FFT_real(1:midpoint);
B_FFT_orig = B_FFT_orig_real(1:midpoint) * B_FFT(1)/B_FFT_orig_real(1);
B_FFT_gw = B_FFT_gw_real(1:midpoint);

% discrepency between SPGs code and this code in the following section---please fix!!!
% create frequency range of fft
dchi = chi(2) - chi(1);
omega = 2 * pi * ([1:midpoint] - 1) / (2*midpoint * dchi);

if strcmp(lower(spectrumType), 'hamada')
    n_minus_miota = omega * dV_dPsi / (2 * pi);
elseif strcmp(lower(spectrumType), 'boozer')
    n_minus_miota = omega * g_Boozer;
else
    error('Unknown spectrum request');
end
n_minus_miota = abs(n_minus_miota);  % this needs to be fixed.
disp(['Delta in n_minus_miota: ' num2str(n_minus_miota(2)-n_minus_miota(1))]);
% find the peaks of the FFT
% To avoid the (0,0) peak, start the search a few indices into the
% arrays
% =====loop====
% make a seek window, and slide it until a peak above the minimum threshold is reached
% Once found, create a fit window around the peak, and try to slide the window to
%           center the peak (slide up to 1/2 the fit window width).
%           After sliding try to fit a gaussian curve to the window.
%           Check the validity of the fit.  Store if it looks good
%           If peak is repeated, jump ahead by (1/2 fit window width + a
%           skip factor)
%           If it doesn't look good, slide seek window ahead by 1/2 fit
%           window width
% ==end loop==
hyper = .5e-4; %  Threshold on minimum peak amplitude to look for.  
seek_window_width = 4;  % Size of the window (in array slots) used during the peek search
seek_window_start = 30; % Where to start, to avoid the (0, 0)-peak
max_centroid_shifts = 10; 
fit_window_width = 4;  % size of window for Gaussian fits attempts, in each direction. total size = 2*fit_window_width
fit_window_skip_size = 4; % SPG used 8.  Must be even #
num_peaks_found = 1;  % The (0,0) peak is the first position of the array
num_failed_fits = 0;  % Keep track of the number of 'failed fits'
gauss_fit_options = optimset('TolX', 1e-15, 'TolFun', 1e-15, 'MaxFunEvals', 50000, 'Display', 'Off');

tic
disp('Seeking peaks in FFT');
toc

peak_position(1) = n_minus_miota(1);
peak_amplitude(1) = B_FFT(1)/2; % The (0,0)-peak 
ii = seek_window_start;

while ii <= midpoint - seek_window_width * 2;  % Set the last start point of the seek window
    seek_window = ii:(ii + seek_window_width - 1);  % Set the seek window
    x_seek = n_minus_miota(seek_window);
    y_seek = B_FFT(seek_window);
    [max_val, index] = max(abs(y_seek));
    if max_val < hyper * peak_amplitude(1)  % Slightly different than SPG's.  My threshold is lower-SPG's is a factor of two higher
        ii = ii + 2;  % slide forward if threshold isn't met
    else
        % Make a fit window and and try to slide it so that that peak is centered
        fit_window = [(ii + index - fit_window_width):(ii + index + fit_window_width - 1 )];
        x_fit = n_minus_miota(fit_window);
        y_fit = B_FFT(fit_window);
        if DEBUG*PAUSE_TIME_FACTOR
            % Do something debuggy
            if (~(exist('plotWindowOpen')))
                h = figure;
                plot(n_minus_miota, B_FFT_orig, 'c+');
                hold on
                plot(n_minus_miota, B_FFT, 'r+');
                plot(n_minus_miota, B_FFT_gw, 'g:');
                plotWindowOpen = 1;
            else
                figure(h);
            end
            widnow_h = plot(x_fit, zeros(size(x_fit)), 'ko');
            x_width_seek = max(x_fit)-min(x_fit);
            cur_axis = axis;
            axis([min(x_fit)-5*x_width_seek max(x_fit)+5*x_width_seek cur_axis(3) cur_axis(4) ]);
        end
        % find 'indices' of the the centroid and mean of the window--these
        % 'indices' may not be integers
        mean_ind = mean(x_fit);
        centroid_ind = sum(x_fit .* abs(y_fit)) / sum(abs(y_fit));

        % slide window to put the centroid close to the mean
        centroid_shifts = 0;
        while ( (centroid_shifts <= max_centroid_shifts) & (abs(mean_ind - centroid_ind) > (x_fit(2)-x_fit(1)) ) )
            if (centroid_ind > mean_ind)  % if centroid is to the right of center, slide window right
                fit_window = fit_window + 1;
            else
                fit_window = fit_window - 1; % otherwise, slide window to the left
            end
            % update the values and repeat the loop
            x_fit = n_minus_miota(fit_window);
            y_fit = B_FFT(fit_window);
            if DEBUG*PAUSE_TIME_FACTOR
                delete(widnow_h);
                widnow_h = plot(x_fit, zeros(size(x_fit)), 'ko');
                pause(.1*PAUSE_TIME_FACTOR);
            end
            mean_ind = mean(x_fit);
            centroid_ind = sum(x_fit .* abs(y_fit)) / sum(abs(y_fit));
            centroid_shifts = centroid_shifts + 1;
        end
        % Try to fit a Gaussian curve to the peak
        if strcmp(lower(spectrumType), 'hamada')
            gauss_fit_width = 0.004;
        elseif strcmp(lower(spectrumType), 'boozer')
            gauss_fit_width = 0.001;
        end
        gauss_init_guess = [sign(mean(y_fit)) * max(abs(y_fit)), mean_ind, gauss_fit_width]; % the initial guess
        [fit_results, resnorm, residual, exitflag] = ...
            lsqcurvefit(@gaussCurveFit, gauss_init_guess, x_fit, y_fit, [-Inf 0 0], [+Inf 600 1], gauss_fit_options);
        if DEBUG*PAUSE_TIME_FACTOR
            x_dense = linspace(x_fit(1), x_fit(end), 100);
            y_dense = gaussCurveFit(fit_results, x_dense);
            y_guess_dense = gaussCurveFit(gauss_init_guess, x_dense);
            plot(x_dense, y_guess_dense, 'b:');
            ph = plot(x_dense, y_dense, 'b');
            x_width = max(x_dense)-min(x_dense); y_height = max(y_dense)-min(y_dense);
            axis([min(x_dense)-5*x_width max(x_dense)+5*x_width min(y_dense)-1.5*y_height max(y_dense)+1.5*y_height ]);
            pause(.2*PAUSE_TIME_FACTOR)
        end

        %         if (exitflag ~= 1)
        %             disp(['Fit not quite perfect near n-m*iota of ' num2str(mean(x_fit))]);
        %         end
        if ( (exitflag > 0) & (fit_results(1) < peak_amplitude(1)) ) % check to see if peak makes sense
            
            % check to see if it is a repeated peak
            if abs( (fit_results(2) - peak_position(num_peaks_found)) ) < 0.01 % repeated peak
                ii =  ii + fit_window_skip_size;
                disp('Repeated fit found');
                if DEBUG*PAUSE_TIME_FACTOR
                    set(ph, 'Color', 'Black');
                    pause(5*PAUSE_TIME_FACTOR);
                end
            else % new peak-store it
                num_peaks_found = num_peaks_found + 1;
                peak_position(num_peaks_found) = fit_results(2);
                peak_amplitude(num_peaks_found) = fit_results(1);
                ii = floor(mean(fit_window)) + floor(fit_window_width/2) + fit_window_skip_size/2;
                disp(['Fit #' num2str(num_peaks_found) ' at n-m(iota) = ' num2str(peak_position(num_peaks_found))]); 
            end
        else
            num_failed_fits = num_failed_fits + 1;
            disp('Another fit bit the dust');
            ii = ii + floor(fit_window_width/2);
            if DEBUG*PAUSE_TIME_FACTOR
                set(ph, 'Color','g');
                pause(5*PAUSE_TIME_FACTOR)
%                 set(ph, 'Color', 'k');
            end
        end % validity check
        if DEBUG*PAUSE_TIME_FACTOR
%             pause(1)
            delete(widnow_h);
        end
    end  % threshold check
end % seek window slide loop
toc

% Try to adjust the peak positions so that the (48,0)-mode is exactly on
% 48.
[peak48, peak48_index] = min( abs( abs(peak_position) - 48) ) ;
renorm = 48 / peak_position(peak48_index);
if ( abs(renorm - 1) <= 0.01)
    peak_position = peak_position * renorm;
    disp(['Correcting spectrum by factor of: ' num2str(renorm)]);
end

disp('Peaks found.  Assigning n, m numbers');

% In HSX, we can cheat by looking for the(n,m) = (4,1) mode to determine
% iota
[values, nm41index] = find( (peak_position > 2.5) & (peak_position < 3.1) );
if (length(nm41index) == 0)
    iota_best = iota_pp;  % no candidates found, use puncture plot value
    disp('No (4,1) candidates found.  Iota from puncture plot is used.');
else
    peaks_41_possible = peak_position(nm41index);
    values_41_possible = peak_amplitude(nm41index);
    [peak_41, ind_41] = max( abs(values_41_possible));
    iota_41 = 4 - peaks_41_possible(ind_41); % got the best guess-compare it to puncture plot
    if abs(iota_41 - iota_pp) / iota_pp > 0.3
        iota_best = iota_pp;    % baaad iota from iota_41: use iota_pp
        disp('Iota from (4,1) list is not close to puncture plot.  Iota from puncture plot is used.');
    else
        iota_best = iota_41;  % iota_41 looks ok.
        disp('Iota from spectrum is used.');
    end
end
toc

% Find n, m.
% sort the peak values and positions from largest to smallest
disp('Finding (n, m) values and magnitudes.');
toc
[bad_pk_amp_sorted, pk_sort_order] = sort(-1 * abs(peak_amplitude));
pk_pos_sorted = peak_position(pk_sort_order);
pk_amp_sorted = peak_amplitude(pk_sort_order);
% For each peak, find the best match of n, m numbers
% Create 5 n by m arrays-one indicates that the (n,m) combination is
% still available, while another keeps track of the peak amplitude that
% goes with it.  A third m by n array keeps track of the mismatch of the
% (n,m)-choice in percentage.  The fourth is the value of (n-m*iota) which
% is compared to each value of pk_pos_sorted.  The fifth is the error of
% the 'next best' (n, m) fit.
% A one-dim array keeps track of (n, m) indices from largest to smallest.
% n_fit_max = 23*4; m_fit_max = 15;  % Maximum on n & m numbers-are these long enough?
n_fit_max = 33*4; m_fit_max = 25;  % Maximum on n & m numbers-are these long enough?
num_fl_per = 4;  % Number of field periods in HSX
% this section for n numbers that are scrictly non-negative
nm_avail = ones(n_fit_max/num_fl_per+1,2*m_fit_max+1); % an 'Inf' value makes the (n, m) combo invalid--see below.
nm_amp = zeros(n_fit_max/num_fl_per+1,2*m_fit_max+1);
nm_error = zeros(n_fit_max/num_fl_per+1,2*m_fit_max+1);
nm_next_best_error = zeros(n_fit_max/num_fl_per+1,2*m_fit_max+1);
nm_pos_value = [0:num_fl_per:n_fit_max]' * ones(1, 2*m_fit_max+1) - iota_best * (ones(n_fit_max/num_fl_per+1, 1)) * [m_fit_max:-1:-m_fit_max];
disp(['Maximum (n - m*iota) =' num2str(max(max(nm_pos_value))) ]);
n_values = 0:num_fl_per:n_fit_max; m_values = m_fit_max:-1:-m_fit_max;
pk_n_sorted = zeros(size(pk_pos_sorted));
pk_m_sorted = zeros(size(pk_pos_sorted));
% this section for n numbers to go from -n_fit_max:n_fit_max
% nm_avail = ones(2*n_fit_max/num_fl_per+1,2*m_fit_max+1); % an 'Inf' value makes the (n, m) combo invalid--see below.
% nm_amp = zeros(2*n_fit_max/num_fl_per+1,2*m_fit_max+1);
% nm_error = zeros(2*n_fit_max/num_fl_per+1,2*m_fit_max+1);
% nm_next_best_error = zeros(2*n_fit_max/num_fl_per+1,2*m_fit_max+1);
% nm_pos_value = [n_fit_max:-num_fl_per:-n_fit_max]' * ones(1, 2*m_fit_max+1) - iota_best * (ones(2*n_fit_max/num_fl_per+1, 1)) * [m_fit_max:-1:-m_fit_max];
% disp(['Maximum (n - m*iota) =' num2str(max(max(nm_pos_value))) ]);
% n_values = n_fit_max:-num_fl_per:-n_fit_max; m_values = m_fit_max:-1:-m_fit_max;
% pk_n_sorted = zeros(size(pk_pos_sorted));
% pk_m_sorted = zeros(size(pk_pos_sorted));

disp('Beginning (n, m) search');
toc
% to match up n, m
% (0, 0) component is automatic
% this section for n numbers that are scrictly non-negative
pk_n_sorted(1) = 1;
pk_m_sorted(1) = m_fit_max+1;
nm_error(1,m_fit_max+1) = 0;
nm_next_best_error(1,m_fit_max+1) = -1;
nm_amp(1,m_fit_max+1) = pk_amp_sorted(1);
nm_avail(n_fit_max/num_fl_per+1,m_fit_max+1) = Inf;
% this section for n numbers to go from -n_fit_max:n_fit_max
% pk_n_sorted(1) = n_fit_max/num_fl_per+1;
% pk_m_sorted(1) = m_fit_max+1;
% nm_error(n_fit_max/num_fl_per+1,m_fit_max+1) = 0;
% nm_next_best_error(n_fit_max/num_fl_per+1,m_fit_max+1) = -1;
% nm_amp(n_fit_max/num_fl_per+1,m_fit_max+1) = pk_amp_sorted(1);
% nm_avail(n_fit_max/num_fl_per+1,m_fit_max+1) = Inf;

for jj = 2:length(pk_pos_sorted)  % skip first one
    % determine mismatch value for all combinations
    nm_mismatch = abs(abs(nm_pos_value) - abs(pk_pos_sorted(jj))) .* nm_avail;  % If nm_avail(n,m) == Inf, then, mismatch = Inf
    % find the (n, m)-indices of the minimum mismatch value
    [min_per_each_m, ind_min_per_n] = min(nm_mismatch, [], 1);
    [min_of_all, ind_min_m] = min(min_per_each_m, [], 2);
    ind_min_n = ind_min_per_n(ind_min_m);
    % record value.  update both the (n by m) arrays and the sorted array
    nm_amp(ind_min_n, ind_min_m) = pk_amp_sorted(jj);
    nm_error(ind_min_n, ind_min_m) = min_of_all / pk_pos_sorted(jj);
    nm_avail(ind_min_n, ind_min_m) = Inf;
    pk_n_sorted(jj) = ind_min_n;
    pk_m_sorted(jj) = ind_min_m;
    % Find 'next best' match
    nm_mismatch_next_best = abs(nm_pos_value - pk_amp_sorted(jj)) .* nm_avail;  % If nm_avail(n,m) == Inf, then, mismatch = Inf
    % find the (n, m)-indices of the minimum mismatch value
    [min_per_each_m_2, ind_min_per_n_2] = min(nm_mismatch_next_best, [], 1);
    [min_of_all_2, ind_min_m_2] = min(min_per_each_m_2, [], 2);
    ind_min_n_2 = ind_min_per_n_2(ind_min_m_2);
    nm_next_best_error(ind_min_n, ind_min_m) = min_of_all_2 /pk_pos_sorted(jj);
end  % and repeat

toc    
disp('All peaks are matched to (n, m) values.  Results are: ');
disp('Peak Amplitude  n   m   error   next_best_error');
disp('====|=========|===|====|=======|================');
for kk = 1:length(pk_pos_sorted)
    disp([num2str( nm_amp(pk_n_sorted(kk), pk_m_sorted(kk)) ) ' === ' ...
        num2str( n_values(pk_n_sorted(kk)) ) ' === ' ...
        num2str( m_values(pk_m_sorted(kk)) ) ' === ' ...
        num2str( nm_error(pk_n_sorted(kk), pk_m_sorted(kk)) ) ' === ' ...
        num2str( nm_next_best_error(pk_n_sorted(kk), pk_m_sorted(kk)) )]);
end
num_failed_fits

if DEBUG
    % reconstitute the B-curve from its spactral components
    if strcmp(lower(spectrumType), 'hamada')
        %         n_minus_miota = omega * dV_dPsi / (2 * pi);
        omega_factor = (2*pi)/dV_dPsi;
    elseif strcmp(lower(spectrumType), 'boozer')
        %         n_minus_miota = omega * g_Boozer;
        omega_factor = 1/g_Boozer;
    end
    
    modB_recon = zeros(size(chi));
    for mm = 1:length(pk_pos_sorted)
        omega_value = (n_values(pk_n_sorted(mm)) - iota_best * m_values(pk_m_sorted(mm))) * omega_factor;
        modB_recon = modB_recon + nm_amp( pk_n_sorted(mm), pk_m_sorted(mm) ) * cos(chi*omega_value);
    end
    %     off_factor = modB ./ modB_recon;
    %     figure
    %     plot(chi, off_factor);

    figure
    subplot(2,1,1);
    plot(chi, modB, 'k');
    hold on;
    plot(chi, modB_recon, 'b');
    subplot(2,1,2);
    plot(chi, (modB-modB_recon)./modB, 'r');

end
% returning the following variables:
% nm_avail, nm_amp, nm_error, nm_next_best_error, pk_n_sorted, pk_m_sorted

%========================================================
%========================================================
%               auxilliary functions
%========================================================
%========================================================

function y = gaussCurveFit(params, x)
% function y = gaussCurveFit(params, x)
%
% y = params(1) * exp( - ( (x - params(2)).^2) ./ (2*(params(3))^2);
%

y = params(1) * exp( -( (x - params(2)).^2) ./ (2*(params(3))^2) );


