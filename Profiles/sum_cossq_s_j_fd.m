function [ dJ_dparam ] = sum_cossq_s_j_fd(s, scsq_params, ind_active )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% ttp_params : pres_scale, A0, A1, A2, A3, A4, A5
curtor = scsq_params(1);
AC_0 = scsq_params(2);
N_cssq = AC_0;
AC_1 = scsq_params(3);
AC_2 = scsq_params(4);
AC_3 = scsq_params(5);
AC_4 = scsq_params(6);
AC_5 = scsq_params(7);

% turn 's' into a single colomn
[rows, cols] = size(s);
if (cols > 1)
    s = s';
end

sum_cossq_profile = [AC_0 AC_1 AC_2 AC_3 AC_4 AC_5];
eqn_1 = sum_cossq_s_j(s, sum_cossq_profile);
eqn_2 = sum(sum_cossq_profile(2:(N_cssq-1)));
J_Norm = curtor * (N_cssq - 1) / eqn_2;
current_solution = J_Norm * eqn_1;
delta_x = 1.0 / (N_cssq - 1);
xi = ([1:N_cssq]-1) * delta_x;


dJ_dparam = zeros(length(s), length(ind_active));
ind_current = 0;
%dJ / d CURTOR
if ismember(1, ind_active)
    ind_current = ind_current + 1;
    dJ_dparam(:, ind_current) = current_solution / curtor;
end

%dJ / d A_0
if ismember(2, ind_active)
    error('<- AC0 should not vary in the sum_cossq model')
end

%dJ / d A_1
if ismember(3, ind_active)
    ind_current = ind_current + 1;
    dJ_dparam(:, ind_current) = -curtor * (N_cssq - 1) * eqn_1 / eqn_2^2 + ...
        J_Norm * 2 * (heaviside(s) - heaviside(s - delta_x)) .* ...
        (cos(pi * (s - xi(1)) / (2 * delta_x))).^2;
end

%dJ / d A_2
if ismember(4, ind_active)
    ind_current = ind_current + 1;
    if (N_cssq == 2)
        dJ_dparam(:, ind_current) = -curtor * (N_cssq - 1) * eqn_1 / eqn_2^2 + ...
            J_Norm * 2 * (heaviside(s - (xi(N_cssq) - delta_x)) - ...
            heaviside(s - xi(N_cssq))) .* ...
            (cos(pi * (s - xi(N_cssq)) / (2 * delta_x))).^2;
    else
        dJ_dparam(:, ind_current) = -curtor * (N_cssq - 1) * eqn_1 / eqn_2^2 + ...
            J_Norm * (heaviside(s - (xi(2) - delta_x)) - ...
            heaviside(s - (xi(2) + delta_x))) .* ...
            (cos(pi * (s - xi(2)) / (2 * delta_x))).^2;
    end
end

%dJ / d A_3
if ismember(5, ind_active)
    ind_current = ind_current + 1;
    if (N_cssq == 3)
        dJ_dparam(:, ind_current) = -curtor * (N_cssq - 1) * eqn_1 / eqn_2^2 + ...
            J_Norm * 2 * (heaviside(s - (xi(N_cssq) - delta_x)) - ...
            heaviside(s - xi(N_cssq))) .* ...
            (cos(pi * (s - xi(N_cssq)) / (2 * delta_x))).^2;
    else
        dJ_dparam(:, ind_current) = -curtor * (N_cssq - 1) * eqn_1 / eqn_2^2 + ...
            J_Norm * (heaviside(s - (xi(3) - delta_x)) - ...
            heaviside(s - (xi(3) + delta_x))) .* ...
            (cos(pi * (s - xi(3)) / (2 * delta_x))).^2;
    end
end

%dJ / d A_4
if ismember(6, ind_active)
    ind_current = ind_current + 1;
    if (N_cssq == 4)
        dJ_dparam(:, ind_current) = -curtor * (N_cssq - 1) * eqn_1 / eqn_2^2 + ...
            J_Norm * 2 * (heaviside(s - (xi(N_cssq) - delta_x)) - ...
            heaviside(s - xi(N_cssq))) .* ...
            (cos(pi * (s - xi(N_cssq)) / (2 * delta_x))).^2;
    else
        dJ_dparam(:, ind_current) = -curtor * (N_cssq - 1) * eqn_1 / eqn_2^2 + ...
            J_Norm * (heaviside(s - (xi(4) - delta_x)) - ...
            heaviside(s - (xi(4) + delta_x))) .* ...
            (cos(pi * (s - xi(4)) / (2 * delta_x))).^2;
    end
end

%dJ / d A_5
if ismember(7, ind_active)
    ind_current = ind_current + 1;
    if (N_cssq == 5)
        dJ_dparam(:, ind_current) = -curtor * (N_cssq - 1) * eqn_1 / eqn_2^2 + ...
            J_Norm * 2 * (heaviside(s - (xi(N_cssq) - delta_x)) - ...
            heaviside(s - xi(N_cssq))) .* ...
            (cos(pi * (s - xi(N_cssq)) / (2 * delta_x))).^2;
    else
        dJ_dparam(:, ind_current) = -curtor * (N_cssq - 1) * eqn_1 / eqn_2^2 + ...
            J_Norm * (heaviside(s - (xi(5) - delta_x)) - ...
            heaviside(s - (xi(5) + delta_x))) .* ...
            (cos(pi * (s - xi(5)) / (2 * delta_x))).^2;
    end
end

%keyboard
%ind_0 = find(s == 0);
% dP_dparam(ind_0,3) = 0;
% dP_dparam(ind_0,6) = 0;
%
% ind_1 = find(s == 1.0);
% dP_dparam(ind_1,4) = 0;
% dP_dparam(ind_1,7) = 0;



