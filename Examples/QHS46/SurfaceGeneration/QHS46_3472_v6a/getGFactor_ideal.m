function fa_G_BS = getGFactor_ideal(mode, surfaceIndexes, VIEW_PLOTS);

[myiota, myJacobian, mypsi, diota_dpsi_new, dJacobian_dpsi_new] = ...
    getSurfaceQuantities(mode, surfaceIndexes, VIEW_PLOTS);

L = 1
M = 4
I = 0
J = 5361 * 48 * 14 / (2*pi)

dV_dPsi = myJacobian * (2*pi)^2;

fa_G_BS = (L * J + M * I) * dV_dPsi ./ (L * myiota - M);

if VIEW_PLOTS
    figure    
    plot(mypsi, fa_G_BS, 'o');
    xlabel('\Psi');
    ylabel('<G_{BS}>_{1/\nu}');
end