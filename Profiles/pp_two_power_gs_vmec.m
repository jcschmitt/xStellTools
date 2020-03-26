function data_out = pp_two_power_gs_vmec(s, pp_in)
% ac(0) * ( (1 - s ** ac(1)) ** ac(2) ) * 
%    ( 1 + Sum[ ac(i) * Exp( -( (s - ac(i+1)) / ac(i+2))** 2 )])
a0 = pp_in(1);
a1 = pp_in(2);
a2 = pp_in(3);
gs_params = pp_in(4:end);
num_gs_params = length(gs_params);
mod_num_gs_params = mod(num_gs_params, 3);
if (mod_num_gs_params ~= 0)
    disp('<----Wrong number of paramters sent into pp_two_power_gs_vmec');
    disp('<----Padding with zeros');
    gs_params = [gs_params zeros(1, 3 - mod_num_gs_params)];
end

num_gs_peaks = length(gs_params) / 3;
gs_mult_factor = ones(size(s));
for ii = 1:num_gs_peaks
    a3 = gs_params(1 + (ii-1)*3);
    a4 = gs_params(2 + (ii-1)*3);
    a5 = gs_params(3 + (ii-1)*3);
    gs_mult_factor = gs_mult_factor + a3 * exp(-( (s - a4) / a5).^2);
end
data_out = a0 * (1 - s.^a1).^a2 .* gs_mult_factor;
%keyboard