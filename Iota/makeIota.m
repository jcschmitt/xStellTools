% function makeIota(fieldlines_data)
%fieldlines_data = read_fieldlines('fieldlines_vmec.h5');

phiIncInDegrees = 72; % W7-X
rAxis = 5.948;
for ii = 1:fieldlines_data.ns
    coords(:,1) = fieldlines_data.R_lines(ii, 1:fieldlines_data.npoinc:end);
    coords(:,2) = fieldlines_data.Z_lines(ii, 1:fieldlines_data.npoinc:end);
    iota(ii) = calculateIota(rAxis, coords, phiIncInDegrees);
end
    
figure;box on;hold on;
plot(iota, 'o');