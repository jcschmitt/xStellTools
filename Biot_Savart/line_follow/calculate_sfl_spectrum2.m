function [nm_amp, pk_n_sorted, pk_m_sorted, pk_pos_sorted, nm_avail, nm_error, nm_next_best_error, n_values, m_values, iota_best] = ...
    calculate_sfl_spectrum(chi, modB, spectrumType, dV_dPsi, g_Boozer, iota_pp, DEBUG)
%               Apply a Gaussian window function to the data
%               Rearrange the data to guarantee that the modB curve has
%               even symmetry. This usually means that the FFT has only
%               real components, but if there are non-stellarator symmetric
%               terms, these will also be captured.
%               FFT the curve
%               find the peaks of the FFT
%               determine best value for iota
%               go find n, m
%               done!

if (nargin < 7)
    DEBUG = 0;
end
DEBUG = 1

min_sym_mode = 1e-12;
min_asym_mode = 1e-12;

min_sym_mode = 4e-3;
min_asym_mode = 1e-11;

PAUSE_TIME_FACTOR = 1.00;  % 0 is fastest, 1 is fast, 100 is sloow.
% make the gaussian window to smooth the ends of the B-curve.
%eta = 4;  % Sets the width of the applied window
eta = 6;  % Sets the width of the applied window
mask_mult = 3;
chi_length = max(chi);  % the length of the chi (in each direction, almost)
% gauss_window = (  eta / ( chi_length * sqrt(2*pi) ) ) * exp( - eta^2 * chi.^2 / (2 * chi_length^2) );
gauss_window = (  eta / ( 1 * sqrt(2*pi) ) ) * exp( - eta^2 * chi.^2 / (2 * chi_length^2) );
%               Apply a Gaussian window function to the data
modB_gw_window_filt = modB .* gauss_window;

%               Rearrange the data to guarantee that the modB curve has even symmetry
%               (so that the FFT has primarily real components)
%    Question: Doesn't this simply undo the way the data was constructed
%    prior to this function?  Implementation (re)consideration required here.

%midpoint = length(modB_gw_window_filt) / 2;
midpoint = find(chi == 0);

if 1
    disp('<----Applying the even symmetry stuff')
    modB_gw_filt_even = [modB_gw_window_filt((midpoint+1):end) modB_gw_window_filt(1:midpoint)];
    modB_orig_even = [modB((midpoint+1):end) modB(1:midpoint)];
    gw_even = [gauss_window((midpoint+1):end) gauss_window(1:midpoint)];
    
else
    disp('<----Undoing the even symmetry stuff')
    modB_gw_filt_even = modB_gw_window_filt;
    modB_orig_even = modB;
    gw_even = gauss_window;
end

if DEBUG
    fh0 = figure;
    figure
    subplot(2,1,1);
    plot(chi, modB, chi, modB_gw_window_filt)
    xlabel('chi'); ylabel('|B|');
    subplot(2,1,2);
    plot(chi, modB_orig_even, chi, modB_gw_filt_even)
    xlabel('chi'); ylabel('|B|');
    legend('orig, even?', 'gw filt, even?');
end

%               FFT the curve

B_FFT_gw_filt_complex = eta*fft(modB_gw_filt_even) / length(modB_gw_filt_even); % fft w normalization
B_FFT_orig_complex = fft(modB_orig_even) / length(modB_orig_even); % fft w/o normalization

B_FFT_gw_even_complex = fft(gw_even) / length(gw_even); % fft w/o normalization
B_FFT_gw_filt_real = real(B_FFT_gw_filt_complex);  % the imaginary portion is pretty small-plot it and check to make sure, tho.
B_FFT_gw_filt_imag = imag(B_FFT_gw_filt_complex);  % the imaginary portion is pretty small-plot it and check to make sure, tho.

B_FFT_orig_real = real(B_FFT_orig_complex);  % the imaginary portion is pretty small-plot it and check to make sure, tho.
B_FFT_orig_imag = imag(B_FFT_orig_complex);  % the imaginary portion is pretty small-plot it and check to make sure, tho.

B_FFT_gw_even_real = real(B_FFT_gw_even_complex);  % the imaginary portion is pretty small-plot it and check to make sure, tho.
B_FFT_gw_filt_real_half = B_FFT_gw_filt_real(1:midpoint);
B_FFT_gw_filt_imag_half = B_FFT_gw_filt_imag(1:midpoint);

B_FFT_orig_half = B_FFT_orig_real(1:midpoint) * B_FFT_gw_filt_real_half(1)/B_FFT_orig_real(1);
B_FFT_real_half = B_FFT_orig_real(1:midpoint);

B_FFT_imag_half = B_FFT_orig_imag(1:midpoint);

B_FFT_gw_even_real_half = B_FFT_gw_even_real(1:midpoint);


B_max_gw_real = max(B_FFT_gw_filt_real);
B_max_real = max(B_FFT_orig_real);

scale_gw2orig = B_max_real / B_max_gw_real;


BFFT_gw_real_scaled = B_FFT_gw_filt_real_half * scale_gw2orig;
BFFT_gw_imag_scaled =  B_FFT_gw_filt_imag_half * scale_gw2orig;

% discrepency between SPGs code and this code in the following section---please fix!!!
% create frequency range of fft
dchi = chi(2) - chi(1);
omega_s = 1/dchi;
numchi = length(chi)
omega_full = omega_s * (0:(numchi/2) / numchi);
omega = 2 * pi * ([1:midpoint] - 1) / (2*midpoint * dchi);
%omega = pi * ([1:midpoint] - 1) / (midpoint * dchi);

if strcmpi(spectrumType, 'hamada')
    n_minus_miota = omega * dV_dPsi / (2 * pi);
elseif strcmpi(spectrumType, 'boozer')
    n_minus_miota = omega * g_Boozer;
else
    error('Unknown spectrum request');
end

% n_minus_miota = abs(n_minus_miota);  % this needs to be fixed.  maybe not?
delta_nmmi = abs(n_minus_miota(2)-n_minus_miota(1));
disp(['Delta in n_minus_miota: ' num2str(delta_nmmi)]);

% find the peaks of the FFT
% fit the largest modes first, remove from spectrum, and research.
% If a repeat peak is found, alert, add new contribution to old, and
% continue
if DEBUG
    axis_current = 0; axis_max = 4; rows_max = 2; cols_max = 2;
    fh_current = figure;
    
    figure
    subplot(2,2,1)
    plot(fftshift(B_FFT_gw_filt_real) * scale_gw2orig);
    hold on;
    plot(fftshift(B_FFT_orig_real));
    legend('gw filt', 'orig')
    
    %ylim([-1e-7, 1e-7])
    %xlim(4200+[-250 250])
    title('shifted FFT of spectrum, real, full');
    grid on
    hold on
    subplot(2,2,3)
    plot(fftshift(B_FFT_gw_filt_imag) * scale_gw2orig);
    hold on;
    plot(fftshift(B_FFT_orig_imag));
    %ylim([-1e-9, 1e-9])
    %xlim(4200+[-250 250])
    title('shifted FFT of spectrum, imag, full');
    grid on
    
    subplot(2,2,2)
    plot(omega, B_FFT_real_half, omega, BFFT_gw_real_scaled);
    %ylim([-1e-7, 1e-7])
    %xlim(4200+[-250 250])
    legend('real', 'gw')
    title('FFT of spectrum, half');
    grid on
    hold on
    subplot(2,2,4)
    plot(omega, B_FFT_imag_half, omega, BFFT_gw_imag_scaled);
    %ylim([-1e-9, 1e-9])
    %xlim(4200+[-250 250])
    title('imag half')
    legend('real', 'gw')
    grid on
    
    %keyboard
end

% find the 0,0 peak
next_task = 'find00';
keep_going = 1;
gauss_fit_options = optimset('TolX', 1e-18, 'TolFun', 1e-18, 'MaxFunEvals', 50000, 'Display', 'Off');

data_new = BFFT_gw_real_scaled;
search_mask = ones(size(data_new));
loop_index = 0;

while keep_going
    task = next_task;
    switch task
        case 'find00'
            data_current = data_new;
            fit_window = 1:(20*eta);
            x_fit = n_minus_miota(fit_window);
            y_fit = data_current(fit_window);
            
            gauss_init_guess = [y_fit(1), 0, .1]; % the initial guess
            loop_index = loop_index + 1;
            [fit_results_sym{loop_index}, resnorm, residual, exitflag] = ...
                lsqcurvefit(@gaussCurveFit, gauss_init_guess, x_fit, y_fit, [-Inf -1 0], [+Inf 1 1], gauss_fit_options);
            data_fit = gaussCurveFit(fit_results_sym{loop_index}, n_minus_miota);
            print_fit(loop_index, fit_results_sym{loop_index});
            
            data_new = data_current - data_fit;
            search_mask_min = 1;
            search_mask_max = ceil(mask_mult*(fit_results_sym{loop_index}(3) / delta_nmmi));
            search_mask(search_mask_min:search_mask_max) = 0;
            
            next_task = 'findsym';
            
            
        case 'findsym'
            data_current = data_new;
            [~, index_max] = max(abs(data_current .* search_mask));
            this_max = data_current(index_max);
            if ( (index_max-10*eta) < 1)
                fit_window = 1:(20*eta);
            elseif ( (index_max+10*eta) > length(data_current) )
                fit_window = ((-20*eta+1):0) + length(data_current);
            else
                fit_window = index_max + [(-10*eta):(10*eta)];
            end
            
            try
                x_fit = n_minus_miota(fit_window);
                y_fit = data_current(fit_window);
            catch
                keyboard
            end
            gauss_init_guess = [this_max, n_minus_miota(index_max), fit_results_sym{1}(3)]; % the initial guess
            min_nmmi = min(n_minus_miota(fit_window));
            max_nmmi = max(n_minus_miota(fit_window));
            
            loop_index = loop_index + 1;
            [fit_results_sym{loop_index}, resnorm, residual, exitflag] = ...
                lsqcurvefit(@gaussCurveFit, gauss_init_guess, x_fit, y_fit, [-Inf min_nmmi 0], [+Inf max_nmmi 1], gauss_fit_options);
            data_fit = gaussCurveFit(fit_results_sym{loop_index}, n_minus_miota);
            print_fit(loop_index, fit_results_sym{loop_index});
            
            data_new = data_current - data_fit;
            
            search_mask_min = index_max - ceil(mask_mult*(fit_results_sym{loop_index}(3) / delta_nmmi));
            search_mask_max = index_max + ceil(mask_mult*(fit_results_sym{loop_index}(3) / delta_nmmi));
            search_mask_min = max(search_mask_min, 1);  % to make sure the lower bound is not 0 or negative
            search_mask_max = min(search_mask_max, length(search_mask));  % to make sure the upper bound is not larger than the total length
            
            search_mask(search_mask_min:search_mask_max) = 0;
            
            if ( abs(fit_results_sym{loop_index}(1) / fit_results_sym{1}(1)) >= min_sym_mode)            
                next_task = 'findsym';
            else
                % prepare to search the stellarator asymmetric terms
                next_task = 'findasym';
                asym_first_loop = 1;
            end
            
        case 'findasym'
            if asym_first_loop
                % initialize the data for the loop.
                data_new = BFFT_gw_imag_scaled;
                search_mask = ones(size(data_new));
                loop_index = 0;
                asym_first_loop = 0;  % otherwise, this section will be run again the next time through
            end
            
            data_current = data_new;
            [~, index_max] = max(abs(data_current .* search_mask));
            this_max = data_current(index_max);
            if ( (index_max-10*eta) < 1)
                fit_window = 1:(20*eta);
            elseif ( (index_max+10*eta) > length(data_current) )
                fit_window = ((-20*eta+1):0) + length(data_current);
            else
                fit_window = index_max + [(-10*eta):(10*eta)];
            end
            
            try
                x_fit = n_minus_miota(fit_window);
                y_fit = data_current(fit_window);
            catch
                keyboard
            end
            gauss_init_guess = [this_max, n_minus_miota(index_max), fit_results_sym{1}(3)]; % the initial guess
            min_nmmi = min(n_minus_miota(fit_window));
            max_nmmi = max(n_minus_miota(fit_window));
            
            loop_index = loop_index + 1;
            [fit_results_asym{loop_index}, resnorm, residual, exitflag] = ...
                lsqcurvefit(@gaussCurveFit, gauss_init_guess, x_fit, y_fit, [-Inf min_nmmi 0], [+Inf max_nmmi 1], gauss_fit_options);
            data_fit = gaussCurveFit(fit_results_asym{loop_index}, n_minus_miota);
            print_fit(loop_index, fit_results_asym{loop_index});
            
            data_new = data_current - data_fit;
            
            search_mask_min = index_max - ceil(mask_mult*(fit_results_asym{loop_index}(3) / delta_nmmi));
            search_mask_max = index_max + ceil(mask_mult*(fit_results_asym{loop_index}(3) / delta_nmmi));
            search_mask_min = max(search_mask_min, 1);  % to make sure the lower bound is not 0 or negative
            search_mask_max = min(search_mask_max, length(search_mask));  % to make sure the upper bound is not larger than the total length
            
            search_mask(search_mask_min:search_mask_max) = 0;
            
            if ( abs(fit_results_asym{loop_index}(1) / fit_results_asym{1}(1)) >= min_asym_mode)            
                next_task = 'findasym';
            else
                next_task = 'finished';
            end
            
        case 'finished'
            break;
            
            
    end
    
    if DEBUG
        x_dense = linspace(x_fit(1), x_fit(end), 400);
        switch task
            case {'find00', 'findsym'}
                this_result = fit_results_sym{loop_index};
            case 'findasym'
                this_result = fit_results_asym{loop_index};
        end
        
        y_dense = gaussCurveFit(this_result, x_dense);
        y_guess_dense = gaussCurveFit(gauss_init_guess, x_dense);
        x_width = max(x_dense)-min(x_dense); y_height = max(y_dense)-min(y_dense);
        
        if (strcmp(task, 'find00'))
            figure(fh0);
            box on; hold on;
            plot(n_minus_miota, B_FFT_real_half, n_minus_miota, BFFT_gw_real_scaled);
            
            plot(x_dense, y_guess_dense, 'g:');
            ph = plot(x_dense, y_dense, 'g');
            plot(n_minus_miota, data_new, 'm:');
            legend('B_{FFT, real, half}', 'B_{FFT, gw, real, scaled}', 'y_{guess,dense}', 'y_{dense}', 'new data')
            axis([min(x_dense)-5*x_width max(x_dense)+5*x_width min((y_dense))-1.5*y_height max(((y_dense)))+1.5*y_height ]);
            pause(.2*PAUSE_TIME_FACTOR)
        end
        
        if 1
            figure(fh_current);
            axis_current = axis_current + 1;
            if (axis_current > axis_max)
                axis_current = 1;
                %keyboard
            end
            subplot(rows_max, cols_max, axis_current);
            hold on;
            hold off;
            plot(n_minus_miota, data_current, n_minus_miota, data_new);
            hold on;
            grid on;
            plot(x_dense, y_guess_dense, 'g:');
            ph = plot(x_dense, y_dense, 'g');
            plot(n_minus_miota(search_mask_min:search_mask_max), 0*n_minus_miota(search_mask_min:search_mask_max), 'ko')
            plot(n_minus_miota, data_current.*search_mask, 'mo')
            legend('data current', ' data new', 'y_{guess,dense}', 'y_{dense}', 'mask off', 'mask_on')
            try
                axis([min(x_dense)-5*x_width max(x_dense)+5*x_width min((y_dense))-1.5*y_height max(((y_dense)))+1.5*y_height ]);
                ylim([min(data_current.*search_mask), max(data_current.*search_mask)])
            catch
                keyboard
            end
            title(['Loop index # ' num2str(loop_index)]);
            pause(.2*PAUSE_TIME_FACTOR)
        end
    end
    
    
end

% 
% if 0
%     
%     % To avoid the (0,0) peak, start the search a few indices into the
%     % arrays
%     % =====loop====
%     % make a seek window, and slide it until a peak above the minimum threshold is reached
%     % Once found, create a fit window around the peak, and try to slide the window to
%     %           center the peak (slide up to 1/2 the fit window width).
%     %           After sliding try to fit a gaussian curve to the window.
%     %           Check the validity of the fit.  Store if it looks good
%     %           If peak is repeated, jump ahead by (1/2 fit window width + a
%     %           skip factor)
%     %           If it doesn't look good, slide seek window ahead by 1/2 fit
%     %           window width
%     % ==end loop==
%     hyper = .5e-4; %  Threshold on minimum peak amplitude to look for.
%     seek_window_width = 4;  % Size of the window (in array slots) used during the peek search
%     seek_window_start = 30; % Where to start, to avoid the (0, 0)-peak
%     max_centroid_shifts = 10;
%     fit_window_width = 4;  % size of window for Gaussian fits attempts, in each direction. total size = 2*fit_window_width
%     fit_window_skip_size = 4; % SPG used 8.  Must be even #
%     num_peaks_found = 1;  % The (0,0) peak is the first position of the array
%     num_failed_fits = 0;  % Keep track of the number of 'failed fits'
%     gauss_fit_options = optimset('TolX', 1e-15, 'TolFun', 1e-15, 'MaxFunEvals', 50000, 'Display', 'Off');
%     
%     tic
%     disp('Seeking peaks in FFT');
%     toc
%     
%     peak_position(1) = n_minus_miota(1);
%     peak_amplitude(1) = B_FFT_gw_filt_real_half(1)/2; % The (0,0)-peak
%     ii = seek_window_start;
%     
%     while ii <= midpoint - seek_window_width * 2;  % Set the last start point of the seek window
%         seek_window = ii:(ii + seek_window_width - 1);  % Set the seek window
%         x_seek = n_minus_miota(seek_window);
%         y_seek = B_FFT_gw_filt_real_half(seek_window);
%         %figure
%         %subplot(2,1,1)
%         %plot(fftshift(B_FFT_gw_filt_real));
%         %ylim([-1e-7, 1e-7])
%         %xlim(4200+[-250 250])
%         %legend('real')
%         %grid on
%         %hold on
%         %subplot(2,1,2)
%         %plot(fftshift(B_FFT_gw_filt_imag));
%         %ylim([-1e-9, 1e-9])
%         %xlim(4200+[-250 250])
%         %legend('imag')
%         %grid on
%         
%         
%         [max_val, index] = max(abs(y_seek));
%         if max_val < hyper * peak_amplitude(1)  % Slightly different than SPG's.  My threshold is lower-SPG's is a factor of two higher
%             ii = ii + 2;  % slide forward if threshold isn't met
%         else
%             % Make a fit window and and try to slide it so that that peak is centered
%             fit_window = [(ii + index - fit_window_width):(ii + index + fit_window_width - 1 )];
%             x_fit = n_minus_miota(fit_window);
%             y_fit = B_FFT_gw_filt_real_half(fit_window);
%             if DEBUG*PAUSE_TIME_FACTOR
%                 % Do something debuggy
%                 if (~(exist('plotWindowOpen')))
%                     h = figure;
%                     plot(n_minus_miota, B_FFT_orig_half, 'c+');
%                     hold on
%                     plot(n_minus_miota, B_FFT_gw_filt_real_half, 'r+');
%                     plot(n_minus_miota, B_FFT_gw_even_real_half, 'g:');
%                     plotWindowOpen = 1;
%                     legend('B_{FFT, orig}', 'B_{FFT}', 'B_{FFT,gw}')
%                 else
%                     figure(h);
%                 end
%                 widnow_h = plot(x_fit, zeros(size(x_fit)), 'ko');
%                 legend('B_{FFT, orig}', 'B_{FFT}', 'B_{FFT,gw}', 'zeros')
%                 x_width_seek = max(x_fit)-min(x_fit);
%                 cur_axis = axis;
%                 axis([min(x_fit)-5*x_width_seek max(x_fit)+5*x_width_seek cur_axis(3) cur_axis(4) ]);
%             end
%             % find 'indices' of the the centroid and mean of the window--these
%             % 'indices' may not be integers
%             mean_ind = mean(x_fit);
%             centroid_ind = sum(x_fit .* abs(y_fit)) / sum(abs(y_fit));
%             
%             % slide window to put the centroid close to the mean
%             centroid_shifts = 0;
%             while ( (centroid_shifts <= max_centroid_shifts) & (abs(mean_ind - centroid_ind) > (x_fit(2)-x_fit(1)) ) )
%                 if (centroid_ind > mean_ind)  % if centroid is to the right of center, slide window right
%                     fit_window = fit_window + 1;
%                 else
%                     fit_window = fit_window - 1; % otherwise, slide window to the left
%                 end
%                 % update the values and repeat the loop
%                 x_fit = n_minus_miota(fit_window);
%                 y_fit = B_FFT_gw_filt_real_half(fit_window);
%                 if DEBUG*PAUSE_TIME_FACTOR
%                     delete(widnow_h);
%                     widnow_h = plot(x_fit, zeros(size(x_fit)), 'ko');
%                     legend('zeros vs x\_fit')
%                     pause(.1*PAUSE_TIME_FACTOR);
%                 end
%                 mean_ind = mean(x_fit);
%                 centroid_ind = sum(x_fit .* abs(y_fit)) / sum(abs(y_fit));
%                 centroid_shifts = centroid_shifts + 1;
%             end
%             % Try to fit a Gaussian curve to the peak
%             if strcmp(lower(spectrumType), 'hamada')
%                 gauss_fit_width = 0.004;
%             elseif strcmp(lower(spectrumType), 'boozer')
%                 gauss_fit_width = 0.001;
%             end
%             gauss_init_guess = [sign(mean(y_fit)) * max(abs(y_fit)), mean_ind, gauss_fit_width]; % the initial guess
%             [fit_results, resnorm, residual, exitflag] = ...
%                 lsqcurvefit(@gaussCurveFit, gauss_init_guess, x_fit, y_fit, [-Inf 0 0], [+Inf 600 1], gauss_fit_options);
%             if DEBUG*PAUSE_TIME_FACTOR
%                 x_dense = linspace(x_fit(1), x_fit(end), 100);
%                 y_dense = gaussCurveFit(fit_results, x_dense);
%                 y_guess_dense = gaussCurveFit(gauss_init_guess, x_dense);
%                 plot(x_dense, y_guess_dense, 'b:');
%                 ph = plot(x_dense, y_dense, 'b');
%                 legend('B_{FFT, orig}', 'B_{FFT}', 'B_{FFT,gw}', 'y_{guess,dense}', 'y_{dense}')
%                 x_width = max(x_dense)-min(x_dense); y_height = max(y_dense)-min(y_dense);
%                 axis([min(x_dense)-5*x_width max(x_dense)+5*x_width min(y_dense)-1.5*y_height max(y_dense)+1.5*y_height ]);
%                 pause(.2*PAUSE_TIME_FACTOR)
%             end
%             
%             %         if (exitflag ~= 1)
%             %             disp(['Fit not quite perfect near n-m*iota of ' num2str(mean(x_fit))]);
%             %         end
%             if ( (exitflag > 0) & (fit_results(1) < peak_amplitude(1)) ) % check to see if peak makes sense
%                 
%                 % check to see if it is a repeated peak
%                 if abs( (fit_results(2) - peak_position(num_peaks_found)) ) < 0.01 % repeated peak
%                     ii =  ii + fit_window_skip_size;
%                     disp('Repeated fit found');
%                     if DEBUG*PAUSE_TIME_FACTOR
%                         set(ph, 'Color', 'Black');
%                         pause(5*PAUSE_TIME_FACTOR);
%                     end
%                 else % new peak-store it
%                     num_peaks_found = num_peaks_found + 1;
%                     peak_position(num_peaks_found) = fit_results(2);
%                     peak_amplitude(num_peaks_found) = fit_results(1);
%                     ii = floor(mean(fit_window)) + floor(fit_window_width/2) + fit_window_skip_size/2;
%                     disp(['Fit #' num2str(num_peaks_found) ' at n-m(iota) = ' num2str(peak_position(num_peaks_found))]);
%                 end
%             else
%                 num_failed_fits = num_failed_fits + 1;
%                 disp('Another fit bit the dust');
%                 ii = ii + floor(fit_window_width/2);
%                 if DEBUG*PAUSE_TIME_FACTOR
%                     set(ph, 'Color','g');
%                     pause(5*PAUSE_TIME_FACTOR)
%                     %                 set(ph, 'Color', 'k');
%                 end
%             end % validity check
%             if DEBUG*PAUSE_TIME_FACTOR
%                 %             pause(1)
%                 delete(widnow_h);
%             end
%         end  % threshold check
%     end % seek window slide loop
%     toc
%     
% end

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
    legend('modB vs chi', 'modB(recon) vs chi', 'err')
    
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


function print_fit(loop_index, data_fit)
disp(['<----Fit #' num2str(loop_index) ' Amplitude: ' num2str(data_fit(1)) ...
    '  @ m-n*iota=' num2str(data_fit(2))]);





