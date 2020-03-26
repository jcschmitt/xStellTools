function pcurr = pp_sum_atan(s_in, ac_in)
% returns I, not dI/ds

% buffer the input array with zeros
len_ac_in = length(ac_in);
ac = [ac_in zeros(1,(21-len_ac_in))];
%       IF (x .ge. one) THEN
%          pcurr = ac(0) + ac(1) + ac(5) + ac(9) + ac(13) + ac(17)
%       ELSE

% Scale factor to ensure integral is 1
scale_factor = 1 / (ac(1) + ac(2) + ac(6) + ac(14) + ac(18));

for ii = 1:length(s_in)
    s = s_in(ii);
    if s >= 1
        pcurr(ii) = ac(1) + ac(2) + ac(6) + ac(14) + ac(18);
    else
        pcurr(ii) = ac(1) +   ...
            ac(2) * (2/pi) * atan(ac(3)*s.^ac(4)./(1-s).^ac(5)) +  ...
            ac(6) * (2/pi) * atan(ac(7)*s.^ac(8)./(1-s).^ac(9)) +  ...
            ac(10) * (2/pi) * atan(ac(11)*s.^ac(12)./(1-s).^ac(13)) +  ...
            ac(14) * (2/pi) * atan(ac(15)*s.^ac(16)./(1-s).^ac(17)) + ...
            ac(18) * (2/pi) * atan(ac(19)*s.^ac(20)./(1-s).^ac(21));
        %       pcurr = ac(0) +                                                         &
        %      &         ac(1) * (2/pi) * atan(ac(2)*x**ac(3)/(1-x)**ac(4)) +            &
        %      &         ac(5) * (2/pi) * atan(ac(6)*x**ac(7)/(1-x)**ac(8)) +            &
        %      &         ac(9) * (2/pi) * atan(ac(10)*x**ac(11)/(1-x)**ac(12)) +         &
        %      &         ac(13) * (2/pi) * atan(ac(14)*x**ac(15)/(1-x)**ac(16)) +        &
        %      &         ac(17) * (2/pi) * atan(ac(18)*x**ac(19)/(1-x)**ac(20))
    end
end
pcurr = scale_factor * pcurr;
