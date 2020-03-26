function y = double_gauss_curve(fit_params, x)
A1 = fit_params(1);
sigma1 = fit_params(2);
mu1 = fit_params(3);
A2 = fit_params(4);
sigma2 = fit_params(5);
mu2 = 0;

if length(fit_params) == 6
    if fit_params(6) == 3
        % take derivative
%         y = A1 .* -(x - mu1) ./ (sigma1^2) .* ( 1 ./ (sigma1 * sqrt(2*pi)) ) .* exp( -(x - mu1).^2 ./ (2 * sigma1^2) ) + ...
%             A2 .* -(x - mu2) ./ (sigma2^2) .* ( 1 ./ (sigma2 * sqrt(2*pi)) ) .* exp( -(x - mu2).^2 ./ (2 * sigma2^2) );
%         disp('blah')
    elseif fit_params(6) == 2
        % return normal gaussian        
        y = A1 .* ( 1 ./ (sigma1 .* sqrt(2*pi)) ) .* exp( -(x - mu1).^2 ./ (2 * sigma1^2) ) + ...
            A2 .* ( 1 ./ (sigma2 .* sqrt(2*pi)) ) .* exp( -(x - mu2).^2 ./ (2 * sigma2^2) );
        %         disp('normal')
    end
else  % return normal gaussian
    y = A1 .* ( 1 ./ (sigma1 .* sqrt(2*pi)) ) .* exp( -(x - mu1).^2 ./ (2 * sigma1^2) ) + ...
        A2 .* ( 1 ./ (sigma2 .* sqrt(2*pi)) ) .* exp( -(x - mu2).^2 ./ (2 * sigma2^2) );
end