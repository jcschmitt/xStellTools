function vmec_data = load_vmec(vmec_woutfile_ext)

len_wout = length(vmec_woutfile_ext);
last3 = [(len_wout-2):len_wout];

if (strcmpi(vmec_woutfile_ext(last3), '.nc') )   
    vmec_woutfile = vmec_woutfile_ext;
else
    vmec_woutfile = ['wout_' vmec_woutfile_ext '.nc'];
end

vmec_data = read_vmec(vmec_woutfile);


% calculate the 'real' Vp, Volume and Well depth (see mercier.f)
vmec_data.phip_real = 2*pi*2*pi * vmec_data.phips;
vmec_data.vp_real = NaN * vmec_data.phip_real;

vmec_data.vp_real(2:end) = double(vmec_data.signgs) * ...
    (4*pi*pi)*vmec_data.vp(2:end) ./ vmec_data.phip_real(2:end);

vmec_data.nsm1 = vmec_data.ns - 1;
vmec_data.hs = 1.0 / vmec_data.nsm1;

% interpoloate/compute on full radial mesh
vmec_data.phip_real2 = NaN;
for jj = 2:(vmec_data.ns-1)
    vmec_data.phip_real2(jj) = 0.5 * ...
        (vmec_data.phip_real(jj+1) + vmec_data.phip_real(jj));
    denom = 1.0 / (vmec_data.hs * vmec_data.phip_real2(jj));
    vmec_data.vpp(jj) = (vmec_data.vp_real(jj+1) - vmec_data.vp_real(jj)) ...
        * denom;
end
vmec_data.phip_real2(vmec_data.ns ) = NaN;
vmec_data.vpp(vmec_data.ns ) = NaN;

        

s = vmec_data.phi ./ vmec_data.phi(end);
rho = sqrt(s);

        
figure
ii = 0;
maxii = 7;

ii = ii + 1;
subplot(maxii,1,ii);
plot(rho,  vmec_data.phip_real, '+--');
xlabel('rho'); ylabel('phip_real');
axis tight;

ii = ii + 1;
subplot(maxii,1,ii);
plot(rho,  vmec_data.vp_real, '+--');
xlabel('rho'); ylabel('vp_real');
axis tight;

ii = ii + 1;
subplot(maxii,1,ii);
plot(rho,  vmec_data.phip_real2, '+--');
xlabel('rho'); ylabel('phip_real2');
axis tight;

ii = ii + 1;
subplot(maxii,1,ii);
plot(rho,  vmec_data.vpp, '+--');
xlabel('rho'); ylabel('vpp');
axis tight;

ii = ii + 1;
subplot(maxii,1,ii);
plot(rho,  vmec_data.vp, '+--');
xlabel('rho'); ylabel('vp');
axis tight;

ii = ii + 1;
subplot(maxii,1,ii);
plot(rho,  (vmec_data.vp(2) - vmec_data.vp) / vmec_data.vp(2), '+--');
xlabel('rho'); ylabel('well from vp');
axis tight;

ii = ii + 1;
subplot(maxii,1,ii);
plot(rho,  (vmec_data.vp_real(2) - vmec_data.vp_real) / vmec_data.vp_real(2), '+--');
xlabel('rho'); ylabel('well from vp_real');
axis tight;






