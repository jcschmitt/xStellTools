function pcurr = pp_pedestal(x, ac)
x=linspace(0, 1, 51);
ac = [0 1 -1 0 0 0 0 0 .5 .75 .8 .1 0 0 0 0 0 0 0 0 0];
pcurr = 0 * x;
% CASE ('pedestal')
% !  Parameterization for I(s)
% !  Added by SPH 2010-05-26
% ioff is 'index offset', and should be 0 (see profile_functions.f) for
% details
    for ii = 8:-1:1
        pcurr = x .* pcurr + ac(ii) / (ii - 0 + 1);
    end
    pcurr = x .* pcurr;

    ii = 9;
    if ac(ii+3) <= 0.0 
        ac(ii:ii+4) = 0;
        ac(ii+3) = 1.0e30;
    else
         ac(ii+4)= 1.0 / ...
             (tanh(2*ac(ii+2)/ac(ii+3)) - tanh(2*(ac(ii+2)-1) / ac(ii+3)) );
    end
    
    a8 = max(ac(ii+8), 0.01);
    a12 = max(ac(ii+12), 0.01);
    
    g1 = (x - ac(ii+7)) / a8;
    g3 = (-ac(ii+7)) / a8;
    g2 = (x - ac(ii+11)) / a12;
    g4 = (-ac(ii+11)) / a12;
    
    pcurr = pcurr + ac(ii+4) * ac(ii+0) * ...
        ( tanh( 2 * ac(ii+2) / ac(ii+3))  - ...
        tanh( 2 * (ac(ii+2) - sqrt(x)) / ac(ii+3) ) ) + ...
        ac(ii+5) * ( tanh(g1) - tanh(g3) ) + ...
        ac(ii+9)*( tanh(g2) - tanh(g4) );

figure;plot(x, pcurr, 'o'); axis tight
