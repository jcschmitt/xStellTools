function pcurr = gauss_trunc(ac_in, s_in)
% returns I, not dI/ds

% copied, almost verbatim from
% VMEC/Sources/Initialiazation_cleanup/profile_functions.f
%!  Truncated Gaussian
%!  I-prime(s) = ac(0) * (exp(-(s/ac(1))**2) - exp(-(1/ac(1))**2)
%!  Gauss-Legendre Quadrature to get I(s)

% convert from matlab indexing to fortran code indexing
ac0 = ac_in(1);
ac1 = ac_in(2);

gln = 10;

glx = [0.01304673574141414, 0.06746831665550774, 0.1602952158504878, ...
       0.2833023029353764, 0.4255628305091844, 0.5744371694908156, ...
       0.7166976970646236, 0.8397047841495122, 0.9325316833444923, ...
       0.9869532642585859];
glw = [0.03333567215434407, 0.0747256745752903, 0.1095431812579910, ...
       0.1346333596549982, 0.1477621123573764, 0.1477621123573764, ...
       0.1346333596549982, 0.1095431812579910, 0.0747256745752903, ...
       0.03333567215434407];
   
pcurr = 0 * s_in;

for ii = 1:length(s_in)
    s = s_in(ii);
    if s >= 1
        s = 1;
    end
    
    for gli = 1:gln
        sp = s * glx(gli);
        pcurr(ii) = pcurr(ii) + glw(gli) * ac0 * (exp(-(sp / ac1)^2) - ...
            exp(-(1 / ac1)^2));
    end
    pcurr(ii) = pcurr(ii) * s;     %! correct for x interval
end

