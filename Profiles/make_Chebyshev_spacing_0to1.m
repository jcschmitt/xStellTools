function pts = make_Chebyshev_spacing_0to1(num_chebypts)
chebypts = num_chebypts;
% chebypts  = 	51 ;
x_k = (fliplr(cos( (2*(1:(chebypts+2))-1) / (2*(chebypts+2)) * pi)) + 1 )' /2;
%s_ntss = x_k;
pts = x_k;