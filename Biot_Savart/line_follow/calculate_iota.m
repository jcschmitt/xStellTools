function [iota] = calculate_iota(rAxis, coords, phiIncInDegrees);

% Using all boxport coords.  Don't do this for if the configuration is not
% stellarator symmetric.
R = coords(:, 1); z = coords(:, 2);  
polAngle = atan2(z, R - rAxis);
polAngle_fix = polAngle;
for ii = 1:length(polAngle_fix)
    if ( polAngle_fix(ii) < 0 )
        polAngle_fix(ii) = polAngle_fix(ii) + 2*pi;
    end
end
diffPolAngle = diff(polAngle_fix);
ind = find(diffPolAngle < 0);
diffPolAngle(ind) = diffPolAngle(ind) + 2*pi;
diffTorAngle = phiIncInDegrees * (pi / 180) * ones(size(diffPolAngle));
%diffTorAngle = 2 * pi * ones(size(diffPolAngle));
totPolAngle = sum(diffPolAngle);
totTorAngle = sum(diffTorAngle);
iota = totPolAngle / totTorAngle;
