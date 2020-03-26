function [] = dfdf()


s = linspace(0, 1, 51);
ds = (s(2) - s(1));
s_mid = s + ds / 2;
s_mid = s_mid([1:(end-1)]);

rho = sqrt(s);
rho_mid = sqrt(s_mid);

Te_0 = 4e3;

% Te ~ (1 - s)
Te_poly = [-1 1];

Te = Te_0 * polyval(Te_poly, s);
Ti = Te;
dTe = diff(Te);

Ne_0 = 2e20;
% Ne ~ (1 - s^5)
Ne_poly = [-1 0 0 0 0 1];
Ne = Ne_0 * polyval(Ne_poly, s);
Ni = Ne;
dNe = diff(Ne);

% Pe ~ (1-s^5) * (1-s) = 1 -s -s^5 - s^6
P_poly = [1 -1 0 0 0 -1 1];
P_calc = 1.6e-19 * 2 * Te_0 * Ne_0 * polyval(P_poly, s);

Pe = Te .* Ne * 1.6e-19;
Pi = Ti .* Ni * 1.6e-19;
P = Pe + Pi;
P_mean = mean(P)
dPe = diff(Pe);





figure
subplot(2,3,1);box on; hold on
plot(rho, Te/1e3, 'o');
title('T_{e} = T_i');
xlabel('\rho');
ylabel('keV');
grid on;

subplot(2,3,2);box on; hold on
plot(rho, Ne/1e20, 'o');
title('N_{e} = N_{i}');
xlabel('\rho');
ylabel('1e20 / m^3');
grid on;

subplot(2,3,3);box on; hold on
plot(rho, P, 'o', rho, P_calc, 'x');
title('P_{total}');
xlabel('\rho');
ylabel('Pa');
grid on;


subplot(2,3,4);box on; hold on
plot(rho_mid, dTe/1e3 / ds, 'o');
title('d T_{e} / ds');
xlabel('\rho');
ylabel('keV / s');
grid on;

subplot(2,3,5);box on; hold on
plot(rho_mid, dNe/1e20 / ds, 'o');
title('dN_{e} / ds');
xlabel('\rho');
ylabel('1e20 / m^3  / s');
grid on;

subplot(2,3,6);box on; hold on
plot(rho_mid, dPe / ds, 'o');
title('dP_{e} / ds');
xlabel('\rho');
ylabel('Pa / s');
grid on;
