function pcurr_out = pcurr(xx, ac_vmec, ac_aux_s, ac_aux_f, pcurr_type)

%  if = 0, then plots are made, otherwise this function runs silent
TEST_THIS_FUNCTION = 0;
% TEST_THIS_FUNCTION = 1;

% ac_matlab = fliplr(ac_vmec);

%  xp               variable for Gauss_legendre quadrature
%  gli              index for Gauss-Legendre quadrature loop
%  gln              order of Gauss-Legendre quadrature
%  glx              array of abscissa values for Gauss-Legendre quadrature
%  glw              array of wieghts for Gauss-Legendre quadrature

%!  Note that the profile that is parameterized is often I-prime, whereas
%!  I(x) (= Integral_from_0_to_x I-prime(s) ds) is the function that pcurr
%!  returns. For the default case of a power series, the integral can be
%!  computed analytically. For other cases, a numerical quadrature is done,
%!  using a 10-point Gauss-Legendre quadrature.


gln = 10;
glx = [0.01304673574141414, 0.06746831665550774, 0.1602952158504878, ...
    0.2833023029353764, 0.4255628305091844, 0.5744371694908156, ...
    0.7166976970646236, 0.8397047841495122, 0.9325316833444923, ...
    0.9869532642585859];
glw = [0.03333567215434407, 0.0747256745752903, 0.1095431812579910, ...
    0.1346333596549982, 0.1477621123573764, 0.1477621123573764, ...
    0.1346333596549982, 0.1095431812579910, 0.0747256745752903, ...
    0.03333567215434407];

pcurr_out = [];
% check_out = [];

for jj = xx
    pcurr_ = 0;
    
    switch lower(regexprep(pcurr_type, ' ', ''))
        %     switch pcurr_type
        %         case (2)
        case 'gauss_trunc'
            % Truncated Gaussian
            % I-prime(s) = ac(0) * (exp(-(s/ac(1))**2) - exp(-(1/ac(1))**2)
            for gli = 1:gln
                xp = jj * glx(gli);
                % JC Schmitt, 2010/05/03
                % I add this check, becuase the routine behaves poorly for
                % ac_vmec(2) < 0.1;
                %                 if abs(ac_vmec(2)) < 0.1
                %                     ac_vmec(2) = 0.1;%  the sign check doesn't matter
                %                     % since the function is even, anyway...
                %                 end
                pcurr_ = pcurr_ + glw(gli) * ac_vmec(1) * ( ...
                    exp(-(xp / ac_vmec(2)).^2) - exp(-(1 / ac_vmec(2)).^2)  );
            end
            pcurr_ = pcurr_ * jj;     % correct for x interval
            %             if isempty(check_out)
            %                 check = NaN;
            %             else
            %                 check = pcurr - pcurr_out(end);
            %             end
            %         case (1)
        case 'power_series'
            % polynomial
            pcurr_ = 0;            
            for ii = length(ac_vmec):-1:1
                pcurr_ = jj*pcurr_ + ac_vmec(ii)/(ii);
            end
            pcurr_ = jj*pcurr_;
            %             check = polyval(polyint(ac_matlab), jj);
        case 'power_series_i'
            %!  I(s)  - not I-prime(s)
            %!  I(s) = Sum(i,0,-)[ac(i) * s ** i]
            %             DO i = UBOUND(ac,1), ioff, -1
            %             pcurr = (pcurr + ac(i))*x
            %             END DO
            
            pcurr_ = 0;
            for ii = length(ac_vmec):-1:1
                pcurr_ = (pcurr_ + ac_vmec(ii)) * jj;
            end
            %             check = polyval(polyint(ac_matlab), jj);
        case 'akima_spline_i'
            % the akima spline routine doesn't like the trailing zeros  in
            % the ac_aux_* arrays, so trim them off
            len_ac_aux_arrays = 1;
            for ii = 2:length(ac_aux_s)
                if ac_aux_s(ii) > ac_aux_s(ii-1)
                    len_ac_aux_arrays = ii;
                end
            end
            pcurr_ = akima(ac_aux_s(1:len_ac_aux_arrays), ...
                ac_aux_f(1:len_ac_aux_arrays), jj);
        case 'sum_atan'
            pcurr_ = sum_atan(ac_vmec, jj);
        otherwise
            error('Unknown current profile specification');
            
    end
    
    pcurr_out = [pcurr_out pcurr_];
    %     check_out = [check_out check];
end

if TEST_THIS_FUNCTION
    figure;box on;hold on;
    plot(xx, pcurr_out, 'b:', 'Linewidth', 2);
    %     plot(xx, check_out, 'r');
end



