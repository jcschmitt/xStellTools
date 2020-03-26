function data_out = pp_two_lorentz(s, AM)
% 
if length(AM) == 5
    AM(6:8) = [1 1 1];
end

C0 = (1 + (1 / (AM(3)^2)) ^AM(4)) ^ -AM(5);
C1 = (1 + (1 / (AM(6)^2)) ^AM(7)) ^ -AM(8);
N0 = 1 - C0;
N1 = 1 - C1;

data_out = AM(1) * (AM(2) / N0) * (1 ./ ...
    ((1 + (s / (AM(3)^2)).^AM(4)).^AM(5)) - C0) + ...
    AM(1) * ((1 - AM(2)) / N1) * (1 ./ ...
    ((1 + (s / (AM(6)^2)).^AM(7)).^AM(8)) - C1);

