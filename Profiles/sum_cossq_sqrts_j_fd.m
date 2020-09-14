function [ dJ_dparam ] = sum_cossq_sqrts_j_fd(s, scsq_params, ind_active )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% scsq_params : curtor, A0, A1, A2, A3, A4, A5
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
eqn_1 = sum_cossq_sqrts_j(s, sum_cossq_profile);
switch N_cssq
    case 2
        error('not developed yet')
    case 3
        warning('not developed yet')
        F1 = (pi^2-4)/(8*pi^2);
        F2 = 1/2;
        F3 = 3/8 + 1/(2 *pi^2);
        F4 = 0;
        F5 = 0;
        eqn_2 = AC_1*F1 + AC_2*F2 + AC_3*F3+ AC_4*F4;
        
    case 4
        F1 = (pi^2-4)/(18*pi^2);
        F2 = 2/9;
        F3 = 4/9;
        F4 = 5/18 + 2/(9 * pi^2);
        F5 = 0;
        eqn_2 = AC_1*F1 + AC_2*F2 + AC_3*F3+ AC_4*F4;
    case 5
        F1 = (pi^2-4)/(32*pi^2);
        F2 = 1/8;
        F3 = 1/4;
        F4 = 3/8;
        F5 = 7/32 + 1/(8*pi^2);
        eqn_2 = AC_1*F1 + AC_2*F2 + AC_3*F3+ AC_4*F4 + ...
            AC_5*F5;
end
J_Norm = curtor / eqn_2;
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
    dJ_dparam(:, ind_current) = -curtor * eqn_1 * F1 / eqn_2^2 + ...
        J_Norm * 2 * (heaviside(sqrt(s)) - heaviside(sqrt(s) - delta_x)) .* ...
        (cos(pi * (sqrt(s) - xi(1)) / (2 * delta_x))).^2;
end

%dJ / d A_2
if ismember(4, ind_active)
    ind_current = ind_current + 1;
    if (N_cssq == 2)
        error('Not devel')
    else
        dJ_dparam(:, ind_current) = -curtor * eqn_1 * F2 / eqn_2^2 + ...
            J_Norm * (heaviside(sqrt(s) - (xi(2) - delta_x)) - ...
            heaviside(sqrt(s) - (xi(2) + delta_x))) .* ...
            (cos(pi * (sqrt(s) - xi(2)) / (2 * delta_x))).^2;
    end
end

%dJ / d A_3
if ismember(5, ind_active)
    ind_current = ind_current + 1;
    if (N_cssq == 3)
        error('Not devel')
    else
        dJ_dparam(:, ind_current) = -curtor * eqn_1 * F3 / eqn_2^2 + ...
            J_Norm * (heaviside(sqrt(s) - (xi(3) - delta_x)) - ...
            heaviside(sqrt(s) - (xi(3) + delta_x))) .* ...
            (cos(pi * (sqrt(s) - xi(3)) / (2 * delta_x))).^2;
    end
end

%dJ / d A_4
if ismember(6, ind_active)
    ind_current = ind_current + 1;
    if (N_cssq == 4)
        dJ_dparam(:, ind_current) = -curtor * eqn_1 * F4 / eqn_2^2 + ...
            J_Norm * (heaviside(sqrt(s) - (xi(4) - delta_x)) - ...
            heaviside(sqrt(s) - (xi(4) + delta_x))) .* ...
            (cos(pi * (sqrt(s) - xi(4)) / (2 * delta_x))).^2;
    else
        dJ_dparam(:, ind_current) = -curtor * eqn_1 * F4 / eqn_2^2 + ...
            J_Norm * (heaviside(sqrt(s) - (xi(4) - delta_x)) - ...
            heaviside(sqrt(s) - (xi(4) + delta_x))) .* ...
            (cos(pi * (sqrt(s) - xi(4)) / (2 * delta_x))).^2;
    end
end

%dJ / d A_5
if ismember(7, ind_active)
    ind_current = ind_current + 1;
    if (N_cssq == 5)
        dJ_dparam(:, ind_current) = -curtor * eqn_1 * F5 / eqn_2^2 + ...
            J_Norm * 2 * (heaviside(sqrt(s) - (xi(N_cssq) - delta_x)) - ...
            heaviside(sqrt(s) - xi(N_cssq))) .* ...
            (cos(pi * (sqrt(s) - xi(N_cssq)) / (2 * delta_x))).^2;
    else
        dJ_dparam(:, ind_current) = -curtor * eqn_1 * F5 / eqn_2^2 + ...
            J_Norm * (heaviside(sqrt(s) - (xi(5) - delta_x)) - ...
            heaviside(sqrt(s) - (xi(5) + delta_x))) .* ...
            (cos(pi * (sqrt(s) - xi(5)) / (2 * delta_x))).^2;
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



