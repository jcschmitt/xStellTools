function [ dP_dparam ] = two_two_power_fd(s, ttp_params, ind_active )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% ttp_params : pres_scale, A0, A1, A2, A3, A4, A5
pres_scale = ttp_params(1);
A_0 = ttp_params(2);
A_1 = ttp_params(3);
A_2 = ttp_params(4);
A_3 = ttp_params(5);
A_4 = ttp_params(6);
A_5 = ttp_params(7);

% turn 's' into a single colomn
[rows, cols] = size(s);
if (cols > 1)
    s = s';
end


dP_dparam = zeros(length(s), length(ind_active));
ind_current = 0;
%dP / d pres_scale
if ismember(1, ind_active)
    ind_current = ind_current + 1;
    dP_dparam(:, ind_current) = A_0 * (  A_3      * (1 - s.^A_1).^A_2 + ...
        (1 - A_3) * (1 - s.^A_4).^A_5  );
end

%dP / d A_0
if ismember(2, ind_active)
    ind_current = ind_current + 1;
    dP_dparam(:, ind_current) = pres_scale * ( A_3       * (1 - s.^A_1).^A_2 + ...
        ( 1 - A_3) * (1 - s.^A_4).^A_5  );
end

%dP / d A_1
if ismember(3, ind_active)
    ind_current = ind_current + 1;
    dP_dparam(:, ind_current) = -pres_scale * A_0 * A_2 * A_3 * s.^A_1 .* log(s) .* ...
        (1 - s.^A_1).^(A_2 - 1);
end

%dP / d A_2
if ismember(4, ind_active)
    ind_current = ind_current + 1;
    dP_dparam(:, ind_current) = pres_scale * A_0 * A_3 * ...
        log(1 - s.^A_1) .* (1 - s.^A_1)  .^A_2;
end

%dP / d A_3
if ismember(5, ind_active)
    ind_current = ind_current + 1;
    dP_dparam(:, ind_current) = pres_scale * A_0 * ( ( 1 - s.^A_1).^A_2 - ...
        ( 1 - s.^A_4).^A_5 );
end

%dP / d A_4
if ismember(6, ind_active)
    ind_current = ind_current + 1;
    dP_dparam(:, ind_current) = -pres_scale * A_0 * (1 - A_3) * A_5 * s.^A_4 .* log(s) .* ...
        (1 - s.^A_4).^(A_5 - 1);
end

%dP / d A_5
if ismember(7, ind_active)
    ind_current = ind_current + 1;
    dP_dparam(:, ind_current) = pres_scale * A_0 * (1 - A_3) * ...
        log(1 - s.^A_4) .* (1 - s.^A_4)  .^A_5;
end

%keyboard
%ind_0 = find(s == 0);
% dP_dparam(ind_0,3) = 0;
% dP_dparam(ind_0,6) = 0;
%
% ind_1 = find(s == 1.0);
% dP_dparam(ind_1,4) = 0;
% dP_dparam(ind_1,7) = 0;



