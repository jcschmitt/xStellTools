function [nm_cos_amp, pk_cos_n_sorted, pk_cos_m_sorted, ...
    pk_pos_sym_sorted, nm_cos_avail, nm_cos_error, ...
    nm_cos_next_best_error, n_values, m_values, iota_best] = ...
    calculate_sfl_spectrum2(chi, modB, spectrumType, dV_dPsi, ...
    g_Boozer, iota_pp, numFieldPeriods, DEBUG_PLOTS)
%               Apply a Gaussian window function to the data
%               Rearrange the data to guarantee that the modB curve has
%               even symmetry. This usually means that the FFT has only
%               real components, but if there are non-stellarator symmetric
%               terms, these will also be captured.
%               FFT the curve
%               find the peaks of the FFT
%               determine best value for iota
%               Find (n, m) values that best match each peak
%               done!

if (nargin < 8)
    DEBUG_PLOTS = 0;
end

DEBUG_PLOTS = 0;
WHEN2STARTPLOTTING = 1;
MAKE_RECON_PLOTS = 1;

%min_sym_mode = 1e-12;
%min_asym_mode = 1e-12;

min_sym_mode = 1e-8;
min_asym_mode = 1e-8 ;

PAUSE_TIME_FACTOR = 1.00;  % 0 is fastest, 1 is fast, 100 is sloow.
% make the gaussian window to smooth the ends of the B-curve.
%eta = 4;  % Sets the width of the applied window
eta = 6;  % Sets the width of the applied window
mask_mult = 3;
chi_length = (max(chi) - min(chi)) / 2;  % the length of the chi (in each direction)
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

if DEBUG_PLOTS
    % open figures for later
    fh0 = figure;
    fh1 = figure;
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
numchi = length(chi);
omega_full = omega_s * (0:(numchi/2) / numchi);
omega = 2 * pi * ([1:midpoint] - 1) / (2*midpoint * dchi);
%omega = pi * ([1:midpoint] - 1) / (midpoint * dchi);

if strcmpi(spectrumType, 'hamada')
    n_minus_miota = omega * dV_dPsi / (2 * pi);
elseif strcmpi(spectrumType, 'boozer')
    n_minus_miota = omega * g_Boozer;
else
    error('<----Unknown spectrum request');
end

% n_minus_miota = abs(n_minus_miota);  % this needs to be fixed.  maybe not?
delta_nmmi = abs(n_minus_miota(2)-n_minus_miota(1));
disp(['<----Delta in n_minus_miota: ' num2str(delta_nmmi)]);

% find the peaks of the FFT
% fit the largest modes first, remove from spectrum, and research.
% If a repeat peak is found, alert, add new contribution to old, and
% continue
if DEBUG_PLOTS
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
                lsqcurvefit(@gaussCurveFit, gauss_init_guess, x_fit, y_fit, [-Inf -1 0], [+Inf 1 1], gauss_fit_options); %#ok<*ASGLU>
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
                disp('<----Issues with fit_window')
                %keyboard
            end
            gauss_init_guess = [this_max, n_minus_miota(index_max), fit_results_sym{1}(3)]; % the initial guess
            min_nmmi = min(n_minus_miota(fit_window));
            max_nmmi = max(n_minus_miota(fit_window));
            
            loop_index = loop_index + 1;
            [fit_results_sym{loop_index}, resnorm, residual, exitflag] = ...
                lsqcurvefit(@gaussCosFit, gauss_init_guess, x_fit, y_fit, [-Inf min_nmmi 0], [+Inf max_nmmi 1], gauss_fit_options);
            data_fit = gaussCosFit(fit_results_sym{loop_index}, n_minus_miota);
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
                lsqcurvefit(@gaussSinFit, gauss_init_guess, x_fit, y_fit, [-Inf min_nmmi 0], [+Inf max_nmmi 1], gauss_fit_options);
            data_fit = gaussSinFit(fit_results_asym{loop_index}, n_minus_miota);
            print_fit(loop_index, fit_results_asym{loop_index});
            
            data_new = data_current - data_fit;
            
            search_mask_min = index_max - ceil(mask_mult*(fit_results_asym{loop_index}(3) / delta_nmmi));
            search_mask_max = index_max + ceil(mask_mult*(fit_results_asym{loop_index}(3) / delta_nmmi));
            search_mask_min = max(search_mask_min, 1);  % to make sure the lower bound is not 0 or negative
            search_mask_max = min(search_mask_max, length(search_mask));  % to make sure the upper bound is not larger than the total length
            
            search_mask(search_mask_min:search_mask_max) = 0;
            
            if ( abs(fit_results_asym{loop_index}(1) / fit_results_sym{1}(1)) >= min_asym_mode)            
                next_task = 'findasym';
            else
                next_task = 'finished';
            end
            
        case 'finished'
            break;
            
            
    end
    
    if DEBUG_PLOTS
        x_dense = linspace(x_fit(1), x_fit(end), 400);
        switch task 
            case {'find00'}
                this_result = fit_results_sym{loop_index};
                y_dense = gaussCurveFit(this_result, x_dense);
                y_guess_dense = gaussCurveFit(gauss_init_guess, x_dense);
            case {'findsym'}
                this_result = fit_results_sym{loop_index};
                y_dense = gaussCosFit(this_result, x_dense);
                y_guess_dense = gaussCosFit(gauss_init_guess, x_dense);
            case 'findasym'
                this_result = fit_results_asym{loop_index};
                y_dense = gaussSinFit(this_result, x_dense);
                y_guess_dense = gaussSinFit(gauss_init_guess, x_dense);
        end
        
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
        
        if (loop_index >= WHEN2STARTPLOTTING)

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
            title(['<----Loop index # ' num2str(loop_index)]);
            pause(.2*PAUSE_TIME_FACTOR)
        end
    end
    
    
end

% put the fits data into individual arrays for easier handling
num_sym_fits = length(fit_results_sym);
num_asym_fits = length(fit_results_asym);
peak_position_sym = zeros(1, num_sym_fits) ;
peak_amplitude_sym = peak_position_sym;
peak_position_asym = zeros(1, num_asym_fits);
peak_amplitude_asym = peak_position_asym;

peak_position_sym(1) = fit_results_sym{1}(2);
peak_amplitude_sym(1) = fit_results_sym{1}(1) / 1.0; % the (0,0) component

for ii = 2:num_sym_fits
    peak_position_sym(ii) = fit_results_sym{ii}(2);
    peak_amplitude_sym(ii) = 2*fit_results_sym{ii}(1);
end
for ii = (1:num_asym_fits)
    peak_position_asym(ii) = fit_results_asym{ii}(2);
    peak_amplitude_asym(ii) = 2*fit_results_asym{ii}(1);
end

% Try to adjust the peak positions so that the {(N, 0), (2*N, 0),
% (3*N, 0), ...} modes are exactly on integers
% 'max_mode_multiple' is motivated by the fact that the max # of
% coils/field period will probably be 12 (or less)  
max_mode_multiple = 12;

peak_fit = [];
peak_mode = [];
for ii = 1:max_mode_multiple
    this_mode = numFieldPeriods * ii;
    [~, this_peak_index] = min( abs( abs(peak_position_sym) - this_mode) );
    this_renorm = this_mode / peak_position_sym(this_peak_index);
    if ( abs( abs(this_renorm) - 1) <= 0.01)
        peak_fit = [peak_fit peak_position_sym(this_peak_index)];
        peak_mode = [peak_mode this_mode];
    end;
end

if length(peak_mode) < 1
    disp('<----Found no peaks that were multiples of the number of field periods.')
    best_renorm = 1;
else
    disp('<----Found peaks that were multiples of the number of field perios');
    peak_mode
    peak_fit
    best_renorm = peak_fit' \ peak_mode';
end

disp(['<----Renormalizing by : ' num2str(best_renorm)]);
peak_position_sym_renormed = peak_position_sym * best_renorm;
peak_position_asym_renormed = peak_position_asym * best_renorm;

iota_best = iota_pp; % why? (answer: in case there is better information from *somewhere* else)

disp('<----Peaks found.  Assigning n, m numbers');

% Find n, m.
% sort the peak values and positions from largest to smallest
disp('<----Finding (n, m) values and magnitudes.');

[~, pk_sym_sort_order] = sort(-1 * abs(peak_amplitude_sym));
pk_pos_sym_sorted = peak_position_sym_renormed(pk_sym_sort_order);
pk_amp_sym_sorted = peak_amplitude_sym(pk_sym_sort_order);

[~, pk_asym_sort_order] = sort(-1 * abs(peak_amplitude_asym));
pk_pos_asym_sorted = peak_position_asym_renormed(pk_asym_sort_order);
pk_amp_asym_sorted = peak_amplitude_asym(pk_asym_sort_order);

% Note: The following range of acceptable n,m combinations (same as
% BOOZ_XFORM) is:
% 0 <= m <= m_max, increments of 1
% m = (0, 1, 2, ..., m_max)
% -n_max <= n <= n_max, inncrements of numFieldPeriods
% n = (-n_max, -n_max+numFieldPeriods, -n_max+numFieldsPeriods*2, ..., n_max)
% If m = 0, n must be non-negative

% For each peak, find the best match of n, m numbers
% Create 5 n by m matrices for each of the sin and cos components.
% For each:
% _avail: A mask that indicates that the (n,m) combination is still
%         available. '1' means the (n,m) component is available, 'inf'
%         means it is already used (or unavailable as a legitimate (n,m)
%         value)  Note: if m = 0, n must be non-negative. So, if m=0, n<0,
%         then the corresponding _avail component should be marked as 'inf'
% _amp: keeps track of the peak amplitudes that goes with it.
% _error: keeps track of the mismatch of the (n,m)-choice in percentage.
% _pos_value: the value of (n-m*iota) which is compared to each value of
% pk_pos_sorted.
% _next_best_error: The fifth is the error of the 'next best' (n, m) fit.
% Also, the following are stored
% pk_*_sorted: A one-dim array to track the (n, m) indices from largest to
% smallest. This points to the INDEX in pk_pos_sym_sorted and pk_s


% n_fit_max = 23*4; m_fit_max = 15;  % Maximum on n & m numbers-are these long enough?
% Maximum on n & m numbers-are these long enough?
num_n_harmonics = 24;
num_m_harmonics = 24;

n_fit_max = num_n_harmonics * numFieldPeriods;
m_fit_max = num_m_harmonics; 

% Iniitialize the arrays.
n_values = -n_fit_max:numFieldPeriods:n_fit_max;
m_values = 0:1:m_fit_max;
n_matrix = n_values' * ones(1, m_fit_max+1);
m_matrix = (ones(2*n_fit_max/numFieldPeriods+1, 1)) * m_values;
nm_pos_value = n_matrix - iota_pp * m_matrix;

nm_cos_avail = ones(2*n_fit_max/numFieldPeriods+1, m_fit_max+1);
% mark the (n<0, m=0) components as unavailable (set them to inf)
nm_cos_avail(1:num_n_harmonics, 1) = Inf;
nm_cos_amp = zeros(2 * num_n_harmonics + 1, m_fit_max + 1);
nm_cos_error = zeros(2 * num_n_harmonics + 1, m_fit_max + 1);
nm_cos_next_best_error = zeros(2 * num_n_harmonics + 1, m_fit_max + 1);
pk_cos_n_sorted = zeros(size(pk_pos_sym_sorted));
pk_cos_m_sorted = zeros(size(pk_pos_sym_sorted));

nm_sin_avail = ones(2*n_fit_max/numFieldPeriods+1, m_fit_max+1);
% mark the (n<=0, m=0) components as unavailable (set them to inf)
nm_sin_avail(1:(num_n_harmonics+1), 1) = Inf;
nm_sin_amp = zeros(2 * num_n_harmonics + 1, m_fit_max + 1);
nm_sin_error = zeros(2 * num_n_harmonics + 1, m_fit_max + 1);
nm_sin_next_best_error = zeros(2 * num_n_harmonics + 1, m_fit_max + 1);
pk_sin_n_sorted = zeros(size(pk_pos_asym_sorted));
pk_sin_m_sorted = zeros(size(pk_pos_asym_sorted));

disp(['<----Minimum (n - m*iota) = ' num2str(min(min(nm_pos_value))) ]);
disp(['<----Maximum (n - m*iota) = ' num2str(max(max(nm_pos_value))) ]);

disp('<----Beginning (n, m) search');

% handle the (0, 0) component
% the index for the (0,0) mode is (num_n_harmonics + 1, 1)
pk_cos_n_sorted(1) = num_n_harmonics + 1;
pk_cos_m_sorted(1) = 1;
nm_cos_error(num_n_harmonics + 1, 1) = 0;
nm_cos_next_best_error(num_n_harmonics + 1, 1) = -1;
nm_cos_amp(num_n_harmonics + 1, 1) = pk_amp_sym_sorted(1);
% and, handle the mask
nm_cos_avail(num_n_harmonics + 1, 1) = Inf; 

% handle the rest of the symmetric components
for jj = 2:length(pk_pos_sym_sorted)  % skip first one
    % determine mismatch value for all combinations
%     nm_mismatch = abs(abs(nm_pos_value) - ...
%                       abs(pk_pos_sym_sorted(jj))) .* nm_cos_avail; 
    %nm_mismatch = abs(nm_pos_value - abs(pk_pos_sym_sorted(jj))) .* nm_cos_avail; 
    nm_mismatch = abs(nm_pos_value - pk_pos_sym_sorted(jj)) .* nm_cos_avail; 
    % find the (n, m)-indices of the minimum mismatch value
    [min_per_each_m, ind_min_per_n] = min(nm_mismatch, [], 1);
    [min_of_all, ind_min_m] = min(min_per_each_m, [], 2);
    ind_min_n = ind_min_per_n(ind_min_m);
    % record value.  update both the (n by m) arrays and the sorted array
    nm_cos_amp(ind_min_n, ind_min_m) = pk_amp_sym_sorted(jj);
    nm_cos_error(ind_min_n, ind_min_m) = min_of_all / pk_pos_sym_sorted(jj);
    nm_cos_avail(ind_min_n, ind_min_m) = Inf;
    pk_cos_n_sorted(jj) = ind_min_n;
    pk_cos_m_sorted(jj) = ind_min_m;
    % Find 'next best' match
    nm_mismatch_next_best = abs(nm_pos_value - pk_amp_sym_sorted(jj)) .* nm_cos_avail;  % If nm_cos_avail(n,m) == Inf, then, mismatch = Inf
    % find the (n, m)-indices of the minimum mismatch value
    [min_per_each_m_2, ind_min_per_n_2] = min(nm_mismatch_next_best, [], 1);
    [min_of_all_2, ind_min_m_2] = min(min_per_each_m_2, [], 2);
    ind_min_n_2 = ind_min_per_n_2(ind_min_m_2);
    nm_cos_next_best_error(ind_min_n, ind_min_m) = min_of_all_2 / pk_pos_sym_sorted(jj);
end  % and repeat

% handle the asymmetric components
for jj = 1:length(pk_pos_asym_sorted)
    % determine mismatch value for all combinations
    nm_mismatch = abs(nm_pos_value - pk_pos_asym_sorted(jj)) .* nm_sin_avail; 
    % find the (n, m)-indices of the minimum mismatch value
    [min_per_each_m, ind_min_per_n] = min(nm_mismatch, [], 1);
    [min_of_all, ind_min_m] = min(min_per_each_m, [], 2);
    ind_min_n = ind_min_per_n(ind_min_m);
    % record value.  update both the (n by m) arrays and the sorted array
    nm_sin_amp(ind_min_n, ind_min_m) = pk_amp_asym_sorted(jj);
    nm_sin_error(ind_min_n, ind_min_m) = min_of_all / pk_pos_asym_sorted(jj);
    nm_sin_avail(ind_min_n, ind_min_m) = Inf;
    pk_sin_n_sorted(jj) = ind_min_n;
    pk_sin_m_sorted(jj) = ind_min_m;
    % Find 'next best' match
    nm_mismatch_next_best = abs(nm_pos_value - pk_amp_asym_sorted(jj)) .* nm_sin_avail;  % If nm_sin_avail(n,m) == Inf, then, mismatch = Inf
    % find the (n, m)-indices of the minimum mismatch value
    [min_per_each_m_2, ind_min_per_n_2] = min(nm_mismatch_next_best, [], 1);
    [min_of_all_2, ind_min_m_2] = min(min_per_each_m_2, [], 2);
    ind_min_n_2 = ind_min_per_n_2(ind_min_m_2);
    nm_sin_next_best_error(ind_min_n, ind_min_m) = min_of_all_2 / pk_pos_asym_sorted(jj);
end  % and repeat


disp('<----All peaks are matched to (n, m) values.  Results are: ');
disp('Index Cos Peak Amplitude   n   m   error   next_best_error');
disp('=====|========|=========|====|====|=======|================');
for kk = 1:length(pk_pos_sym_sorted)
    disp([num2str(kk) ' === ' ...
        num2str( nm_cos_amp(pk_cos_n_sorted(kk), pk_cos_m_sorted(kk)) ) ' === ' ...
        num2str( n_values(pk_cos_n_sorted(kk)) ) ' === ' ...
        num2str( m_values(pk_cos_m_sorted(kk)) ) ' === ' ...
        num2str( nm_cos_error(pk_cos_n_sorted(kk), pk_cos_m_sorted(kk)) ) ' === ' ...
        num2str( nm_cos_next_best_error(pk_cos_n_sorted(kk), pk_cos_m_sorted(kk)) )]);
end

disp('Index  Sin Peak Amplitude   n   m   error   next_best_error');
disp('======|========|=========|====|====|=======|================');
for kk = 1:length(pk_pos_asym_sorted)
    disp([num2str(kk) ' === ' ...
        num2str( nm_sin_amp(pk_sin_n_sorted(kk), pk_sin_m_sorted(kk)) ) ' === ' ...
        num2str( n_values(pk_sin_n_sorted(kk)) ) ' === ' ...
        num2str( m_values(pk_sin_m_sorted(kk)) ) ' === ' ...
        num2str( nm_sin_error(pk_sin_n_sorted(kk), pk_sin_m_sorted(kk)) ) ' === ' ...
        num2str( nm_sin_next_best_error(pk_sin_n_sorted(kk), pk_sin_m_sorted(kk)) )]);
end


if MAKE_RECON_PLOTS
    % reconstitute the B-curve from its spactral components
    if strcmp(lower(spectrumType), 'hamada')
        %         n_minus_miota = omega * dV_dPsi / (2 * pi);
        omega_factor = (2*pi)/dV_dPsi;
    elseif strcmp(lower(spectrumType), 'boozer')
        %         n_minus_miota = omega * g_Boozer;
        omega_factor = 1/g_Boozer;
    end
    
    modB_sym_recon = zeros(size(chi));
    modB_asym_recon = zeros(size(chi));
    %  Handle cos-terms
    for mm = 1:length(pk_pos_sym_sorted)
        omega_value = (n_values(pk_cos_n_sorted(mm)) - iota_pp * m_values(pk_cos_m_sorted(mm))) * omega_factor;
        if mm == 1
            modB_sym_recon = modB_sym_recon + nm_cos_amp( pk_cos_n_sorted(mm), pk_cos_m_sorted(mm) ) * cos(chi*omega_value);
        else
            modB_sym_recon = modB_sym_recon + nm_cos_amp( pk_cos_n_sorted(mm), pk_cos_m_sorted(mm) ) * cos(chi*omega_value);
        end
    end
    
    %  Handle sin-terms
    for mm = 1:length(pk_pos_asym_sorted)
        omega_value = (n_values(pk_sin_n_sorted(mm)) - iota_pp * m_values(pk_sin_m_sorted(mm))) * omega_factor;
        modB_asym_recon = modB_asym_recon + nm_sin_amp( pk_sin_n_sorted(mm), pk_sin_m_sorted(mm) ) * cos(chi*omega_value);
    end

    modB_recon = modB_sym_recon + modB_asym_recon;
    
    figure
    subplot(2,1,1);
    plot(chi, modB, 'k');
    hold on;
    plot(chi, modB_sym_recon, 'r');
    plot(chi, modB_recon, 'b');
    legend('modB', 'modB(recon) sym terms only', 'modB(recon) total')
    subplot(2,1,2);
    plot(chi, 100*(modB-modB_sym_recon)./modB, 'r');
    hold on
    plot(chi, 100*(modB-modB_recon)./modB, 'b');
    legend('Error (%) - Sym Terms', 'Error (%) - Total')
    
end
% returning the following variables:
% nm_cos_avail, nm_cos_amp, nm_cos_error, nm_cos_next_best_error, pk_cos_n_sorted, pk_cos_m_sorted

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


function y = gaussCosFit(params, x)
% function y = gaussCurveFit(params, x)
%
% y = params(1) * exp( - ( (x - params(2)).^2) ./ (2*(params(3))^2);
%

y = params(1) * (exp( -( (x - params(2)).^2) ./ (2*(params(3))^2) ) + ...
    exp( -( (x + params(2)).^2) ./ (2*(params(3))^2) ));


function y = gaussSinFit(params, x)
% function y = gaussCurveFit(params, x)
%
% y = params(1) * exp( - ( (x - params(2)).^2) ./ (2*(params(3))^2);
%

y = params(1) * (exp( -( (x - params(2)).^2) ./ (2*(params(3))^2) ) - ...
    exp( -( (x + params(2)).^2) ./ (2*(params(3))^2) ));


function print_fit(loop_index, data_fit)
disp(['<----Fit #' num2str(loop_index) ' Amplitude: ' num2str(data_fit(1)) ...
    '  @ m-n*iota=' num2str(data_fit(2))]);





