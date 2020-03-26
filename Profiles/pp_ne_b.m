function data_out = pp_ne_b(s, pp_in, ne_min)
% 
a0 = pp_in(1);
a1 = pp_in(2);
a2 = pp_in(3);
a3 = pp_in(4);

if nargin < 3
    ne_min = 0;
end

data_out = a0 + a1 * (1 - s.^a2).^a3;

ind = find(data_out<ne_min);
data_out(ind) = ne_min;

