function data_out = pp_two_power_v3fit(s, pp_in)
% 
a0 = pp_in(1);
a1 = pp_in(2);
a2 = pp_in(3);
a3 = pp_in(4);

data_out = a0 + a1 * (1 - s.^a2).^a3;
