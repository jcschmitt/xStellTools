function [iota] = calculateIota(rAxis, coords, phiIncInDegrees)

R = coords(:, 1); z = coords(:, 2);  % Using all boxport coords.  Don't do this for if earthField is on!!!
polAngle = atan2(z, R - rAxis);
for ii = 1:length(polAngle)
    if ( polAngle(ii) < 0 )
        polAngle(ii) = polAngle(ii) + 2*pi;
    end
end
diffPolAngle = diff(polAngle);
ind = find(diffPolAngle < 0);
diffPolAngle(ind) = diffPolAngle(ind) + 2*pi;
diffTorAngle = phiIncInDegrees * (pi / 180) * ones(size(diffPolAngle));
totPolAngle = sum(diffPolAngle);
totTorAngle = sum(diffTorAngle);
iota = totPolAngle / totTorAngle;

