function [] = ...
    iota_evolution_generate(evolution_parameter_file_in)
% function [] = ...
%     iota_evolution(evolution_parameter_file_in, SHOW_PLOTS, RUN_EVOLUTION)
%
% evolution_parameter_file_in = filename of specification

% FIX THIS PATH FOR YOUR FILE STRUCTURE
% data_path='D:\JohnS\Plots and Notes\Time Evolution\simulation_data\input';
% output_path = 'D:\JohnS\Plots and Notes\Time Evolution\simulation_data\output'
% evolution_parameter_file = [data_path '\' evolution_parameter_file_in];
% % attempts to whack the '.txt'off of the input filename
% output_file = [output_path '\' evolution_parameter_file_in(1:end-4)];
% options

% data_path= ['/Users/jschmitt/Documents/HSX Material/Work PC/D drive/' ...
%     'JohnS/Plots and Notes/Time Evolution/simulation_data/input'];
% output_path = ['/Users/jschmitt/Documents/HSX Material/Work PC/D drive/' ...
%     'JohnS/Plots and Notes/Time Evolution/simulation_data/output'];
data_path= '/Users/schmittj/Documents/MATLAB/Time Evolution/simulation_data/input';
output_path = '/Users/schmittj/Documents/MATLAB/Time Evolution/simulation_data/output';
evolution_parameter_file = [data_path '/' evolution_parameter_file_in];
output_file = [output_path '/' evolution_parameter_file_in(1:end-4)];

t_end = 9.5;
% t_end = 0.051; % [sec]  % end of simulation time
% t_end = 0.069; % [sec]  % end of simulation time
% t_end = 1.001; % [sec]  % end of simulation time
delta_t = .01; % msec timestep in simulation time

mu_0 = pi * 4e-7; % permeability of free space

% load the variables
disp(['Opening ' evolution_parameter_file]);
[fid_input, message] = fopen(evolution_parameter_file,'r');
if isempty(message)
    disp(['Loading data from ' evolution_parameter_file]);
    
    fillerLine = fgetl(fid_input);
    fillerLine = fgetl(fid_input);
    fillerLine = fgetl(fid_input);
    susceptance_matrix_datafile = fgetl(fid_input);
    
    fillerLine = fgetl(fid_input);
    decay_type_in = (fgetl(fid_input));
    
    fillerLine = fgetl(fid_input);
    pressure_poly_form = str2num(fgetl(fid_input));
    
    fillerLine = fgetl(fid_input);
    pressure_poly_terms = str2num(fgetl(fid_input));
    
    fillerLine = fgetl(fid_input);
    pressure_time_evolve_form = str2num(fgetl(fid_input));
    
    fillerLine = fgetl(fid_input);
    pressure_tau_on = str2num(fgetl(fid_input));
    
    fillerLine = fgetl(fid_input);
    pressure_breakdown_delay = str2num(fgetl(fid_input));

    fillerLine = fgetl(fid_input);
    JdotB_poly_form = str2num(fgetl(fid_input));
    
    fillerLine = fgetl(fid_input);
    JdotB_poly_terms_x = str2num(fgetl(fid_input));

    fillerLine = fgetl(fid_input);
    JdotB_poly_terms_y = str2num(fgetl(fid_input));
    
    fillerLine = fgetl(fid_input);
    JdotB_tau_on = str2num(fgetl(fid_input));
    
    fillerLine = fgetl(fid_input);
    eta_parallel_poly_form = str2num(fgetl(fid_input));
    
    fillerLine = fgetl(fid_input);
    eta_parallel_poly_terms_x = str2num(fgetl(fid_input));

    fillerLine = fgetl(fid_input);
    eta_parallel_poly_terms_y = str2num(fgetl(fid_input));
    
    fillerLine = fgetl(fid_input);
    Itor_form = str2num(fgetl(fid_input));
    
    fillerLine = fgetl(fid_input);
    Itor_ss_value = str2num(fgetl(fid_input));
    
    fillerLine = fgetl(fid_input);
    Itor_decay_rate = str2num(fgetl(fid_input));
    
    fillerLine = fgetl(fid_input);
    Phi_LCFS_in = str2num(fgetl(fid_input));
end
fclose(fid_input);
disp('Done reading evolution parameter file');

disp('Loading susceptance matrix data file');
sus_mat_data = load(susceptance_matrix_datafile);

disp('Susceptance matrix data loaded');

disp('Generating static variables and derivatives in s-space')
disp('s-space means "normalized toroidal flux"');
% For the matlab pdepe formulation incorporated here, the variables are:
%  x = spatial = s;  t = time; u = S11*iota+S12; dudx
%  m = 0:  implies a 'slab-like' geometry
% Need the following terms for the evolution solution, as a function of s
% Phi_a : (constant) flux enclosed by LCFS (last closed flux surface)
% eta_parallel: Parallel conductivity
% dV_ds: d(Volume) / d(s)
% dP_ds: d(plasma pressure) / d(s)
% fa_B2: flux surface average of B^2  <B^2>
% S11, S12, :  Susceptance matrix components and their
%              radial derivitives
% Need the following terms as a function of time, t
% I_tor_total:  Total toroidal current enclosed by the LCFS
% Need the following terms as a function of time, t, ans space, s
% fa_JdotB:  Flux surface average of J_noninductive dot B  <J_NI . B>

num_surfs_vmec = length(sus_mat_data.r_eff);
s_local_in = ((1:num_surfs_vmec)-1)/(num_surfs_vmec-1);

% core_rad_grid_add_pts = 10; % for fine resolution on axis
core_rad_grid_add_pts = 0; % coarse resoluiont.  maybe closer to experiment? 
disp(['Adding ' num2str(core_rad_grid_add_pts) ...
    ' pts to radial (s) grid for simulation']);
% the points are spaced evenly in rho-space
s_dense_core_pts = linspace(0, sqrt(s_local_in(2)), 2+core_rad_grid_add_pts).^2;
s_local = [s_dense_core_pts s_local_in(3:end)];
disp(['<----Length of s_local: ' num2str(length(s_local))])
disp('overriding vmec s density.');
s_local = linspace(0,1,31);
%s_local = s_local_in;
disp(['<----' num2str(length(s_local)) ' s points']);
%
mid_s_local = 0.5*(s_local(2:end)+s_local(1:end-1));
rho_local = sqrt(s_local);
rho_local_in = sqrt(s_local_in);
mid_rho_local = sqrt(mid_s_local);
a0 = sus_mat_data.r_eff(end);
r_local = rho_local * a0;
r_local_in = rho_local_in * a0;
mid_r_local = mid_rho_local * a0;

% chain rule stuff:  from (r_eff) r- to s- space
dr_ds = a0^2 ./ (2*r_local);
ds_dr = (2*r_local) ./ (a0^2);
dr_ds_in = a0^2 ./ (2*r_local_in);
ds_dr_in = (2*r_local_in) ./ (a0^2);
ds_drho = 2*rho_local;
drho_ds = 1./(2*rho_local);

disp('Applying fix to <B^2> at edge')
% sus_mat_data.Bsquared(end) = csaps(s_local(1:end-1), ...
%     sus_mat_data.Bsquared(1:end-1), [], 1);
sus_mat_data.Bsquared(end) = interp1(s_local_in(2:end-1), ...
    sus_mat_data.Bsquared(2:end-1), 1, 'spline', 'extrap');
disp('Applying fix to S22 at edge')
sus_mat_data.S22(end) = interp1(s_local_in(2:end-1), ...
    sus_mat_data.S22(2:end-1), 1, 'spline', 'extrap');
disp('Applying fix to iota at edge')
sus_mat_data.S22(end) = interp1(s_local_in(2:end-1), ...
    sus_mat_data.S22(2:end-1), 1, 'spline', 'extrap');

% set each of the quantities for s-space
get_set_Phi_LCFS(Phi_LCFS_in);
get_set_DECAY_TYPE(decay_type_in);

figure; subplot(4,5,1); box on;
get_set_Phiprime(s_local_in, Phi_LCFS_in*ones(size(s_local_in)));
subplot(4,5,2); box on;
get_set_conductivity_profile(eta_parallel_poly_form, ...
    [eta_parallel_poly_terms_x eta_parallel_poly_terms_y]);
subplot(4,5,3); box on;
get_set_Vprime(s_local_in, sus_mat_data.Vprime .* dr_ds_in);
subplot(4,5,4); box on;
get_set_Pprime(pressure_tau_on, polyder(pressure_poly_terms), ...
    pressure_time_evolve_form);
subplot(4,5,5); box on;
get_set_faB(s_local_in, sus_mat_data.B);
subplot(3,4,5); box on;hold on;
get_set_faB2(s_local_in, sus_mat_data.Bsquared);
subplot(3,4,6); box on;hold on;
% the ds_dr term is an absorbed term from the change in derivative/chain rule stuff
get_set_S11(s_local_in, sus_mat_data.S11 .* ds_dr_in);
subplot(3,4,7); box on;
% the ds_dr term is an absorbed term from the change in derivative/chain rule stuff
get_set_S12(s_local_in, sus_mat_data.S12 .* ds_dr_in);
subplot(3,4,8); box on;
plot(s_local_in, -sus_mat_data.S12 ./ sus_mat_data.S11, '.');
title('\iota_{vac}');
subplot(3,4,9); box on;
% the ds_dr term is an absorbed term from the change in derivative/chain rule stuff
get_set_S21(s_local_in, sus_mat_data.S21 .* ds_dr_in);
subplot(3,4,10); box on;
% the ds_dr term is an absorbed term from the change in derivative/chain rule stuff
get_set_S22(s_local_in, sus_mat_data.S22 .* ds_dr_in);

subplot(3,4,11); box on;
get_set_JdotB(JdotB_tau_on, [JdotB_poly_terms_x JdotB_poly_terms_y], JdotB_poly_form);
subplot(3,4,12); box on;
get_set_I_tor(Itor_decay_rate, Itor_ss_value, Itor_form);

get_set_pressure_delay_time(pressure_breakdown_delay);

pause(1)

tic
disp('<----Beginning evolution')
m = 0;
pdefun = @iotaevol_pdefun1;
icfun = @iotaevol_icfun1;
bcfun = @iotaevol_bcfun1;
tspan = 0:delta_t:t_end;
sol = pdepe(m, pdefun, icfun, bcfun, s_local, tspan);
disp('<----Evolution done');
toc
disp('<----Re-evaluation');
for ii = 1:length(sol(:,1,1))
    [uout(ii,:), duoutdx(ii,:)] = pdeval(m, s_local, sol(ii, :), s_local);
end
disp('<----Re-evaluation done');

delta_iota_thing = sol(:,:,1);

num_surfs_simulation = length(s_local);
disp(['<----number of surfs in simulation: ' ...
    num2str(num_surfs_simulation)]);

for ii = 1:num_surfs_simulation
    S11_local(ii) = get_set_S11(s_local(ii)); %#ok<AGROW>
    S12_local(ii) = get_set_S12(s_local(ii)); %#ok<AGROW>
    S21_local(ii) = get_set_S21(s_local(ii)); %#ok<AGROW>
    S22_local(ii) = get_set_S22(s_local(ii)); %#ok<AGROW>
    sigma_local(ii) = get_set_conductivity_profile(s_local(ii)); %#ok<AGROW>
    eta_local(ii) = 1./sigma_local(ii); %#ok<AGROW>
    PhiPrime_local(ii) = get_set_Phiprime(s_local(ii)); %#ok<AGROW>
    B2_local(ii) = get_set_faB2(s_local(ii)); %#ok<AGROW>
    VPrime_local(ii) = get_set_Vprime(s_local(ii)); %#ok<AGROW>
end

% QL1 = <B2>V'/mu0  QR1 = (I Psi' + F Phi')
% QL2 = <J.B>V' QR2 = mu0(F I' - I F')
% QL3 = -p'V'   QR3 = I' Psi' + F' Phi'
iota = zeros(length(tspan), num_surfs_simulation);
PPrime_local = zeros(length(tspan), num_surfs_simulation);
PsiPrime_local = zeros(length(tspan), num_surfs_simulation);
I_tor = zeros(length(tspan), num_surfs_simulation);
F_pol = zeros(length(tspan), num_surfs_simulation);
faJdotB_local = zeros(length(tspan), num_surfs_simulation);
dI_tor = zeros(length(tspan), num_surfs_simulation);
dF_pol = zeros(length(tspan), num_surfs_simulation);
QL1 = zeros(length(tspan), num_surfs_simulation);
QR1 = zeros(length(tspan), num_surfs_simulation);
QL2 = zeros(length(tspan), num_surfs_simulation);
QR2 = zeros(length(tspan), num_surfs_simulation);
QL3 = zeros(length(tspan), num_surfs_simulation);
QR3 = zeros(length(tspan), num_surfs_simulation);
EdotB_local = zeros(length(tspan), num_surfs_simulation);
V_loop = zeros(size(tspan));
% J_inductive = EdotB_local;
I_tor_net = 0*tspan;
dI_tor_halfmesh = zeros(length(tspan), num_surfs_simulation - 1);
dF_pol_halfmesh = zeros(length(tspan), num_surfs_simulation - 1);
I_tor_halfmesh = zeros(length(tspan), num_surfs_simulation - 1);
F_pol_halfmesh = zeros(length(tspan), num_surfs_simulation - 1);
dF_pol_dt = zeros(length(tspan), num_surfs_simulation);
J_pol_halfmesh = zeros(length(tspan), num_surfs_simulation - 1);
iota_halfmesh = zeros(length(tspan), num_surfs_simulation - 1);
PsiPrime_halfmesh = zeros(length(tspan), num_surfs_simulation - 1);

QL1_halfmesh = zeros(length(tspan), num_surfs_simulation - 1);
QR1_halfmesh = zeros(length(tspan), num_surfs_simulation - 1);
QL2_halfmesh = zeros(length(tspan), num_surfs_simulation - 1);
QR2_halfmesh = zeros(length(tspan), num_surfs_simulation - 1);
QL3_halfmesh = zeros(length(tspan), num_surfs_simulation - 1);
QR3_halfmesh = zeros(length(tspan), num_surfs_simulation - 1);
J_inductive_halfmesh = zeros(length(tspan), num_surfs_simulation - 1);
EdotB_halfmesh = zeros(length(tspan), num_surfs_simulation - 1);
    
% Variable to help monitor how close we are to done
status_perc = 0;
status_inc = 1;
for ii = 1:length(tspan)
    if ( (100*ii/(length(tspan))) > status_perc)
        disp(['<---- Percent complete: ' num2str(status_perc)])
        status_perc = status_perc + status_inc;
    end
    I_tor_net(ii) = get_set_I_tor(tspan(ii));
    for jj = 1:num_surfs_simulation
        if jj == 1
            iota(ii,jj) = (delta_iota_thing(ii,1) - S12_local(2)) / S11_local(2);
        else
            iota(ii,jj) = (delta_iota_thing(ii,jj) - S12_local(jj)) / S11_local(jj);
        end
        PsiPrime_local(ii,jj) = iota(ii,jj) * PhiPrime_local(jj);
        I_tor(ii,jj) =  PhiPrime_local(jj) * delta_iota_thing(ii,jj) / mu_0;
        F_pol(ii,jj) = (S21_local(jj) * PsiPrime_local(ii,jj) + S22_local(jj) * PhiPrime_local(jj)) / mu_0;
        faJdotB_local(ii,jj) = get_set_JdotB(s_local(jj), tspan(ii));
        PPrime_local(ii,jj) = get_set_Pprime(s_local(jj), tspan(ii));
    end
    %     for jj = 1:num_surfs_simulation
    %         if jj == 1
    %             dI_tor(ii,jj) = diff(I_tor(ii,1:2)) ./ diff(s_local(1:2));
    %             dF_pol(ii,jj) = diff(F_pol(ii,1:2)-F_pol(1,1:2)) ./ diff(s_local(1:2));
    %             %                         dF_pol(ii,jj) = diff(F_pol(ii,1:2)) ./ diff(s_local(1:2));
    %         elseif jj == num_surfs_simulation
    %             dI_tor(ii,jj) = diff(I_tor(ii,(jj-1):jj )) ./ diff(s_local((jj-1):jj ));
    %             dF_pol(ii,jj) = diff(F_pol(ii,(jj-1):jj )-F_pol(1,(jj-1):jj )) ./ diff(s_local((jj-1):jj ));
    %             %                         dF_pol(ii,jj) = diff(F_pol(ii,(jj-1):jj )) ./ diff(s_local((jj-1):jj ));
    %         else
    %             dI_tor(ii,jj) = diff(I_tor(ii,(jj-1):2:(jj+1) )) ./ diff(s_local((jj-1):2:(jj+1) ));
    %             dF_pol(ii,jj) = diff(F_pol(ii,(jj-1):2:(jj+1) )-F_pol(1,(jj-1):2:(jj+1) )) ./ diff(s_local((jj-1):2:(jj+1) ));
    %             %                         dF_pol(ii,jj) = diff(F_pol(ii,(jj-1):2:(jj+1) )) ./ diff(s_local((jj-1):2:(jj+1) ));
    %         end
    %     end
    dI_tor_halfmesh(ii,:) = diff(I_tor(ii,:)) ./ diff(s_local);
    dF_pol_halfmesh(ii,:) = diff(F_pol(ii,:)) ./ diff(s_local);
    I_tor_halfmesh(ii,:) = interp1(s_local, I_tor(ii,:), mid_s_local);
    F_pol_halfmesh(ii,:) = interp1(s_local, F_pol(ii,:), mid_s_local);
    dF_pol_dt(ii,:) = F_pol(ii, :) - F_pol(1, :);
    J_pol_halfmesh(ii,:) = ( diff(dF_pol_dt(ii,:)) ./ diff(s_local) ) / ...
        Phi_LCFS_in;
    
    for jj = 1:num_surfs_simulation-1
        sigma_halfmesh(jj) = get_set_conductivity_profile(mid_s_local(jj)); %#ok<AGROW>
        eta_halfmesh(jj) = 1./sigma_halfmesh(jj); %#ok<AGROW>
        PhiPrime_halfmesh(jj) = get_set_Phiprime(mid_s_local(jj)); %#ok<AGROW>
        B2_halfmesh(jj) = get_set_faB2(mid_s_local(jj)); %#ok<AGROW>
        VPrime_halfmesh(jj) = get_set_Vprime(mid_s_local(jj)); %#ok<AGROW>
        PPrime_halfmesh(ii,jj) = get_set_Pprime(mid_s_local(jj), tspan(ii)); %#ok<AGROW>
        faJdotB_halfmesh(ii,jj) = get_set_JdotB(mid_s_local(jj), tspan(ii)); %#ok<AGROW>
    end
    iota_halfmesh(ii,:) = interp1(s_local, iota(ii,:), mid_s_local);
    PsiPrime_halfmesh(ii,:) = iota_halfmesh(ii,:) .* PhiPrime_halfmesh;
    
    QL1_halfmesh(ii,:) = B2_halfmesh .* VPrime_halfmesh / mu_0;
    QR1_halfmesh(ii,:) = (I_tor_halfmesh(ii,:) .* PsiPrime_halfmesh(ii,:) + ...
        F_pol_halfmesh(ii,:) .* PhiPrime_halfmesh(:)');
    QL2_halfmesh(ii,:) = faJdotB_halfmesh(ii,:) .* VPrime_halfmesh;
    
    QR2_halfmesh(ii,:) = mu_0 * ...
        (F_pol_halfmesh(ii,:) .* dI_tor_halfmesh(ii,:) - ...
        I_tor_halfmesh(ii,:) .* dF_pol_halfmesh(ii,:));
    J_inductive_halfmesh(ii,:) = ( QR2_halfmesh(ii,:) - ...
        QL2_halfmesh(ii,:) ) ./ VPrime_halfmesh;
    EdotB_halfmesh(ii,:) = eta_halfmesh .* ( QR2_halfmesh(ii,:) - ...
        QL2_halfmesh(ii,:) ) ./ VPrime_halfmesh;
    
    dir_B = sign(Phi_LCFS_in); % based on convention of VMEC and jcschmitt
    % cylindrical approximation using geometric values from vmec, vacuum
    % run.  The dir_B term takes care of the 'dotB' part.
    %V_loop(ii) = dir_B*1.2*2*pi*EdotB_halfmesh(ii,end);
    V_loop(ii) = dir_B*2*pi*EdotB_halfmesh(ii,end);
    
    QL3_halfmesh(ii,:) = -PPrime_halfmesh(ii,:) .* VPrime_halfmesh;
    QR3_halfmesh(ii,:) = dI_tor_halfmesh(ii,:) .* PsiPrime_halfmesh(ii,:) + ...
        dF_pol_halfmesh(ii,:) .* PhiPrime_halfmesh;
    
    %     QR2(ii,:) = mu_0 * (F_pol(ii,:) .* dI_tor(ii,:) - I_tor(ii,:) .* dF_pol(ii,:));
    %     J_inductive(ii,:) = ( QR2(ii,:) - QL2(ii,:) ) ./ VPrime_local;
    %     EdotB_local(ii,:) = eta_local .* ( QR2(ii,:) - QL2(ii,:) ) ./ VPrime_local;
    %     QL3(ii,:) = -PPrime_local(ii,:) .* VPrime_local;
    %     QR3(ii,:) = dI_tor(ii,:) .* PsiPrime_local(ii,:) + dF_pol(ii,:) .* PhiPrime_local;
end

disp(['Saving data to ' output_file]);
save(output_file)


disp('Done with iota evolution');
end


function [c, f, s] = iotaevol_pdefun1(x, t, u, dudx)
% global a_MinorRadius
% sigma_0 = 1e7;
% x is the normalized flux. Normally this is 's', but that variable
% is taken
mu_0 = pi * 4e-7;
Phi_a = get_set_Phi_LCFS;
S11 = get_set_S11(x);
%S12 = get_set_S12(x);
%S21 = get_set_S21(x);
%S22 = get_set_S22(x);
sigma_p = get_set_conductivity_profile(x);
eta_p = 1./sigma_p;
Vprime = get_set_Vprime(x);
Pprime = get_set_Pprime(x, t);
faB2 = get_set_faB2(x);
% tau_ee = 0.001;
JdotB = get_set_JdotB(x, t);

c = (Phi_a^2 / S11);
s = zeros(size(u));

f = eta_p * Vprime * (Pprime * u + ...
    faB2 * dudx / mu_0 - JdotB);
% x
% if x == 1
%     x
%     c = 19e-6;
%     f = eta_p * Vprime * faB2 * dudx / mu_0;
% end

% chunk_1 = (Phi_a^2 * dudx / mu_0) * ( S22 + (S21 + u) * (u - S12) / S11 );
% f = eta_p * (Pprime * u * Vprime + ...
%     chunk_1 - JdotB * Vprime);

% disp(['pdefun1   x:' num2str(x)  '  t:' num2str(t) ...
%     '  u:' num2str(u) '  dudx:' num2str(dudx) '  c:' num2str(c) ...
%     '  f:' num2str(f) '  s:' num2str(s)]);
% if x == 0.99
%     pause
% end

end

function u_init = iotaevol_icfun1(s)
% S11 = get_S11(s);
% S12 = get_S12(s);
% iota_init = -S12/S11;
% if s == 0
%     iota_init = iotaevol_icfun1(.02);
% end
u_init = 0;
% s
% iota_init

end

function [pl, ql, pr, qr] = iotaevol_bcfun1(xl, ul, xr, ur, t)

DECAY_TYPE = get_set_DECAY_TYPE;
% DECAY_TYPE = 1;  % Set by measured net toroidal current
% DECAY_TYPE = 2;  % L/R inductive-resistive relaxation

% W7X Numbers, Configuration Ref 1
mu_0 = 4 * pi * 1e-7;
R_major = 5.507988; %1.2;
R_minor = 0.538789; %0.112;
F_shape = .25;
L_ext = mu_0 * R_major * (log ( 8*R_major/R_minor ) - 2 + F_shape);
%L_ext = mu_0 * R_major * 2.32; % Experimental 9/24/07 27 (cal 5)

S11 = get_set_S11(xr);
S12 = get_set_S12(xr);
%S21 = get_set_S21(x);
%S22 = get_set_S22(x);
Vprime = get_set_Vprime(xr);
Pprime = get_set_Pprime(xr, t);
JdotB = get_set_JdotB(xr, t);
% Phi_a = get_set_Phi_LCFS;
sigma_p = get_set_conductivity_profile(xr);
eta_p = 1./sigma_p;
faB2 = get_set_faB2(xr);
Phi_prime = get_set_Phiprime(xr);

pl = ul;
ql = 0;

switch lower(DECAY_TYPE)
    case 'netcurrent'
        % DECAY_TYPE == 1  % set by observed (measured) net current
        I_tor = get_set_I_tor(t);
        pr = ur - mu_0 * I_tor / Phi_prime;
        qr = 0;
    case 'inductive'
        % DECAY_TYPE == 2 % Inductive decay
        t_delay_override = 2;
        I_Inf = get_set_I_tor(Inf);
        tau_LR = L_ext * R_minor^2 / (2 * eta_p * R_major);
        pr = 2 * Phi_prime * eta_p * R_major * ...
            I_Inf * exp(-(t-t_delay_override) / tau_LR) / R_minor^2;
        qr = 1;
    case 'decaytest'
        % f = eta_p * Vprime * (Pprime * u + ...
        %     faB2 * dudx / mu_0 - JdotB)
        % want dudx = 0; have pr + qr * f = 0
        ur * Phi_prime / mu_0
        %Pprime
        %JdotB
        %pr = -eta_p * Vprime * (Pprime * ur - JdotB);
        %qr = 1;
        %pr = -eta_p * Vprime * (Pprime * ur - JdotB) + ...
        %    eta_p * Vprime * faB2 * ur * R_major / (R_minor * L_ext);
        %qr = 1;
        pr = -( + ...
            ur * (1 + R_major * S11 / R_minor) - S12) * ...
            (faB2 * eta_p * Vprime) / (S11 * L_ext);
        qr =  1;
    otherwise
        error('what?');
end



% disp(['bcfun1   xl:' num2str(xl)  '  ul:' num2str(ul) ...
%     '  xr:' num2str(xr) '  ur:' num2str(ur) '  t:' num2str(t) ...
%      '  pl:' num2str(pl) '  ql:' num2str(ql) '  pr:' num2str(pr) ...
%       '  qr:' num2str(qr)]);

end


function sigma_plasma = get_set_conductivity_profile(s, Y_IN)
persistent sigma_sub sigma_x sigma_y sigma_type
if nargin == 2
    sigma_type = s;
    if sigma_type == 1 % poly form
        sigma_sub = Y_IN;
        ss = linspace(0,1,100);
        semilogy(ss, polyval(sigma_sub,ss),'.');
    elseif sigma_type == 2 % ([x] [y])
        length_input = length(Y_IN);
        sigma_x = Y_IN(1:length_input/2);
        sigma_y = Y_IN(length_input/2+1:end);
        % disp('Reducing sigma by 1/3')
        % sigma_y = Y_IN(length_input/2+1:end)/3;
        plot(sigma_x, sigma_y, '.');
    end
    title('\sigma');
    %disp('Skin time based on central conductivity')
    %a1=trapz(sigma_x, sigma_y)
    %a2=trapz(sigma_x, 1./sigma_y)
    %4*pi*1e-7*.12^2 * a1
    %4*pi*1e-7*.12^2 /a2
elseif nargin == 1
    if sigma_type == 1 % poly form
        sigma_plasma = polyval(sigma_sub, s);
        % what is the factor of .5 for? reduced conductivity
    elseif sigma_type == 2 % ([x] [y])
        sigma_plasma = interp1(sigma_x, sigma_y, s);
    end
    
    ind_toosmall = find(sigma_plasma <= 0);
    if ~isempty(ind_toosmall)
        sigma_plasma(ind_toosmall) = 1e2;
    end
else
    error('blah');
end

end

function JdotB = get_set_JdotB(s, JdotB_poly_IN, JdotB_form_IN)
persistent JdotB_sub JdotB_form tau_on
persistent JdotB_x JdotB_y

if nargin == 3
    tau_on = s;
    JdotB_form = JdotB_form_IN;
    if JdotB_form == 1 || JdotB_form == 2 || JdotB_form == 5
        JdotB_sub = JdotB_poly_IN;
    elseif JdotB_form == 4
        length_input = length(JdotB_poly_IN);
        JdotB_x = JdotB_poly_IN(1:length_input/2);
        JdotB_y = JdotB_poly_IN(length_input/2+1:end);
    else
        error('what?');
    end
    ss = linspace(0,1,100);
    if JdotB_form == 1
        plot(ss, polyval(JdotB_sub,ss),'.');
    elseif JdotB_form == 2
        plot(ss, double_gauss_curve([JdotB_sub JdotB_form], ss),'.');
    elseif JdotB_form == 5
        plot(ss, W7X_ECCD_gauss_curve([JdotB_sub], ss),'.');
    elseif JdotB_form == 4
        plot(JdotB_x, JdotB_y, '.');
    end
    title('JdotB');
elseif nargin == 2
    t = JdotB_poly_IN;
    if JdotB_form == 1
        JdotB = polyval(JdotB_sub, s);
    elseif JdotB_form == 2
        JdotB = double_gauss_curve([JdotB_sub JdotB_form], s);
    elseif JdotB_form == 5
        JdotB = W7X_ECCD_gauss_curve([JdotB_sub], s);
    elseif JdotB_form == 4
        % use interp1 for now.  You should specify s=0 and s=1
        JdotB = interp1(JdotB_x, JdotB_y, s);
    end
    T_ECHOFF = 60.0;
    %     T_ECHON = 0.002;  % about a 2.ms delay b4 anything happens.
    % T_ZERO = 0.06; 
    T_ZERO = 0.0; 
    %     tau_off = tau_on;
    tau_off = 0.001;
    
    pressurce_delay_time = ...
        get_set_pressure_delay_time();

    
    if t <= pressurce_delay_time
        JdotB = 0 * JdotB;
    elseif t <= T_ECHOFF
        exp_decay = 1 - exp(-(t-pressurce_delay_time) / tau_on);
        JdotB = JdotB * exp_decay;
    elseif t <= T_ZERO
        exp_decay = 1 - exp(-(T_ECHOFF-pressurce_delay_time) / tau_on);
        JdotB = JdotB * exp_decay;
        exp_decay_2 = exp( - (t-T_ECHOFF) / tau_off);
        JdotB = JdotB * exp_decay_2;
    else
        JdotB = 0 * JdotB;
    end
    
    
    %     Phi_LCFS = get_set_Phi_LCFS;
    %     JdotB = JdotB / Phi_LCFS;
else
    error('blah');
end

end

function S11 = get_set_S11(s, Y_IN)
persistent s_sub y_sub
if nargin == 2
    s_sub = s;
    y_sub = Y_IN;
    %     figure;box on;hold on;
    plot(s_sub,y_sub,'.');
    title('S11');
elseif nargin == 1
    S11 = interp1(s_sub, y_sub, s);
else
    error('blah');
end

end

function S12 = get_set_S12(s, Y_IN)
persistent s_sub y_sub
if nargin == 2
    s_sub = s;
    y_sub = Y_IN;
    %     figure;box on;hold on;
    plot(s_sub,y_sub,'.');
    title('S12');
elseif nargin == 1
    S12 = interp1(s_sub, y_sub, s);
else
    error('blah');
end
end

function S21 = get_set_S21(s, Y_IN)
persistent s_sub y_sub
if nargin == 2
    s_sub = s;
    y_sub = Y_IN;
    %     figure;box on;hold on;
    plot(s_sub,y_sub,'.');
    title('S21');
elseif nargin == 1
    S21 = interp1(s_sub, y_sub, s);
else
    error('blah');
end
end

function S22 = get_set_S22(s, Y_IN)
persistent s_sub y_sub
if nargin == 2
    s_sub = s;
    y_sub = Y_IN;
    if isnan(y_sub(1))
        y_sub(1) = y_sub(3) - (s_sub(3)- s_sub(1)) * ...
            ( ( y_sub(3) - y_sub(2) ) / (s_sub(3) - s_sub(2) ) );
    end
    
    %     figure;box on;hold on;
    plot(s_sub,y_sub,'.');
    title('S22');
elseif nargin == 1
    S22 = interp1(s_sub, y_sub, s);
else
    error('blah');
end
end

function Vprime = get_set_Vprime(s, Y_IN)
persistent s_sub y_sub
if nargin == 2
    s_sub = s;
    y_sub = Y_IN;
    if isnan(y_sub(1))
        y_sub(1) = y_sub(3) - (s_sub(3)- s_sub(1)) * ...
            ( ( y_sub(3) - y_sub(2) ) / (s_sub(3) - s_sub(2) ) );
    end
    %     figure;box on;hold on;
    plot(s_sub,y_sub,'.');
    title('V''');
elseif nargin == 1
    Vprime = interp1(s_sub, y_sub, s);
else
    error('blah');
end

end

function faB2 = get_set_faB2(s, Y_IN)
persistent s_sub y_sub
if nargin == 2
    s_sub = s;
    y_sub = Y_IN;
    %     figure;box on;hold on;
    plot(s_sub,y_sub,'.');
    title('\langleB^2\rangle');
elseif nargin == 1
    faB2 = interp1(s_sub, y_sub, s);
else
    error('blah');
end
end

function faB = get_set_faB(s, Y_IN)
persistent s_sub y_sub
if nargin == 2
    s_sub = s;
    y_sub = Y_IN;
    %     figure;box on;hold on;
    plot(s_sub,y_sub,'.');
    title('\langleB\rangle');
elseif nargin == 1
    faB2 = interp1(s_sub, y_sub, s);
else
    error('blah');
end

end

function Pprime = get_set_Pprime(s, poly_sub_in, poly_form_in)
persistent tau_pres poly_sub_inf poly_form
if nargin == 3
    tau_pres = s;
    poly_sub_inf = poly_sub_in;
    poly_form = poly_form_in;
    ss = linspace(0,1,100);
    plot(ss, polyval(poly_sub_inf,ss),'.');
    title('P''');
elseif nargin == 2
    t = poly_sub_in;
    Pprime = polyval(poly_sub_inf, s);
    pressurce_delay_time = get_set_pressure_delay_time;
    if poly_form == 1
        Pprime = polyval(poly_sub_inf, s) * (1 - ...
            exp(- (t-pressurce_delay_time) / tau_pres));
    elseif poly_form == 2
        if t <= tau_pres
            Pprime = polyval(poly_sub_inf, s) * (1 + (t - tau_pres - pressurce_delay_time)/tau_pres);
        else
            Pprime = polyval(poly_sub_inf, s);
        end
    end
    if t <= pressurce_delay_time
        Pprime = 0 * Pprime;
    end
else
    error('blah');
end

end

function pressurce_delay_time = get_set_pressure_delay_time(VAL_in)
persistent delay_time
if nargin == 1
    delay_time = VAL_in;
elseif nargin == 0
    pressurce_delay_time = delay_time;
else
    error('blah');
end
    
end

function Phiprime = get_set_Phiprime(s, Y_IN)
persistent s_sub y_sub
if nargin == 2
    s_sub = s;
    y_sub = Y_IN;
    %     figure;box on;hold on;
    plot(s_sub,y_sub,'.');
    title('\Phi''');
elseif nargin == 1
    Phiprime = interp1(s_sub, y_sub, s);
else
    error('blah');
end

end

function I_tor = get_set_I_tor(t, I_tor_inf_IN, I_tor_form_IN)
persistent tau_ee I_tor_inf I_tval I_Ival I_tor_form
if nargin == 3
    tau_ee = t;
    I_tor_inf = I_tor_inf_IN;
    I_tor_form = I_tor_form_IN;
    t_plot = linspace(0,.1,100);
    %     figure;box on;hold on;
    if I_tor_form == 1
        ptemp = I_tor_inf * (1 - exp(- t_plot / tau_ee));
    elseif I_tor_form == 2
        for ii = 1:length(t_plot)
            if t_plot(ii) < tau_ee
                ptemp(ii) = I_tor_inf * (1 + (t_plot(ii) - tau_ee)/tau_ee);
            else
                ptemp(ii) = I_tor_inf;
            end
        end
    elseif I_tor_form == 3
        I_tval = t;
        I_Ival = I_tor_inf_IN;
        t_plot = I_tval;
        ptemp = I_Ival;
    end
    plot(t_plot, ptemp,'.');
    title('Itor(t)');
elseif nargin == 1
    if I_tor_form == 1
        t_delay = 2.0;
        if t <= t_delay
            I_tor = 0;
        else
            I_tor = I_tor_inf * (1 - exp(- (t - t_delay) / tau_ee));
        end
    elseif I_tor_form == 2
        if t <= tau_ee
            I_tor = I_tor_inf * (1 + (t - tau_ee)/tau_ee);
        else
            I_tor = I_tor_inf;
        end
    elseif I_tor_form == 3
        if t < 0
            I_tor = 0;
        elseif t < max(I_tval)
            I_tor = interp1(I_tval, I_Ival, t, 'linear');
        else
            I_tor = 0;
        end
    else
        error('blah');
    end
    % disp(['get_I_tor   t:' num2str(t)  '  I_tor:' num2str(I_tor)]);
    
end
end

function Phi_a = get_set_Phi_LCFS(Val_IN)
persistent Phi_LCFS
if nargin == 1
    Phi_LCFS = Val_IN;
elseif nargin == 0
    Phi_a = Phi_LCFS;
else
    error('blah');
end
end

function DECAY_TYPE_out = get_set_DECAY_TYPE(Val_IN)
persistent DECAY_TYPE
if nargin == 1
    DECAY_TYPE = Val_IN;
elseif nargin == 0
    DECAY_TYPE_out = DECAY_TYPE;
else
    error('blah');
end
end

