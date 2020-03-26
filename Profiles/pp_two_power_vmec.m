function data_out = pp_two_power_vmec(s, pp_in)
% 
a0 = pp_in(1);
a1 = pp_in(2);
a2 = pp_in(3);

data_out = a0 * (1 - s.^a1).^a2;
