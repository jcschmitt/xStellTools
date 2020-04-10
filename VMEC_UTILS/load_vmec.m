function vmec_data = load_vmec(vmec_woutfile_ext)

make_well_plot = false;
make_pb_plot = false;

len_wout = length(vmec_woutfile_ext);
last3 = [(len_wout-2):len_wout];

if (strcmpi(vmec_woutfile_ext(last3), '.nc') )
    vmec_woutfile = vmec_woutfile_ext;
else
    vmec_woutfile = ['wout_' vmec_woutfile_ext '.nc'];
end

vmec_data = read_vmec(vmec_woutfile);



% Use vmec info to extrapolate the magnetic well depth

% calculate the 'real' Vp, Volume and Well depth (see mercier.f)
% for 'typical' vmec run, the phips array is on the 'half mesh', and its
% first value is '0' because the value on-axis has not been computed by
% VMEC.
vmec_data.phip_real = 2*pi*2*pi * vmec_data.phips;
vmec_data.vp_real = NaN * vmec_data.phip_real;

vmec_data.vp_real(2:end) = double(vmec_data.signgs) * ...
    (4*pi*pi)*vmec_data.vp(2:end) ./ vmec_data.phip_real(2:end);


% interpoloate/compute on full radial mesh
% vmec_data.phip_real_fullmesh(vmec_data.ns) = NaN;
% vmec_data.vpp(vmec_data.ns) = NaN;
vmec_data.nsm1 = vmec_data.ns - 1;
vmec_data.hs = 1.0 / vmec_data.nsm1;

% interpoloate/compute on full radial mesh
vmec_data.vpp = NaN;
vmec_data.phip_real_fullmesh = NaN;
for jj = 2:(vmec_data.ns-1)
    vmec_data.phip_real_fullmesh(jj) = 0.5 * ...
        (vmec_data.phip_real(jj+1) + vmec_data.phip_real(jj));
    denom = 1.0 / (vmec_data.hs * vmec_data.phip_real_fullmesh(jj));
    vmec_data.vpp(jj) = (vmec_data.vp_real(jj+1) - vmec_data.vp_real(jj)) ...
        * denom;
end
%vmec_data.vpp(1) = 2*vmec_data.vpp(2) - vmec_data.vpp(3);

% Now we assing the last
vmec_data.phip_real_fullmesh(vmec_data.ns ) = NaN;
vmec_data.vpp(vmec_data.ns ) = NaN;
%vmec_data.phip_real_fullmesh(vmec_data.ns ) = NaN;
%vmec_data.vpp(vmec_data.ns ) = NaN;





s = vmec_data.phi ./ vmec_data.phi(end);
rho = sqrt(s);


% make a parabolic fit to the first few indices to extrapoloate to s=0
% for some variables
temp_indices = [2:10];
disp(['<----Extrapolatiing VP_REAL to axis based on surfaces ' num2str(temp_indices(1)) ' to ' num2str(temp_indices(end))]);
x_temp = s(temp_indices);
y_temp = vmec_data.vp_real(temp_indices);
fit_temp = polyfit(x_temp, y_temp, 2);
x0_temp = polyval(fit_temp, 0);
vmec_data.vp_real_fit = vmec_data.vp_real;
vmec_data.vp_real_fit(1) = x0_temp;

%vmec_data.welldepth =  (vmec_data.vp_real(2) - vmec_data.vp_real) / vmec_data.vp_real(2);
vmec_data.welldepth =  (vmec_data.vp_real_fit(1) - vmec_data.vp_real_fit) / vmec_data.vp_real_fit(1);


% Peeling-Ballooning
% (-p')(-V'') + (iota') <J_{||} B> / <B^2> > 0
% ' == derivative w.r.t. toroidal flux (psi)
% V'' = 1/(4*pi^2)  d^2 (V) / d (psi)^2
% V'' and p' are not available at s=0, s=1, because phip_real_fullmesh is
% not on the full grid.
mVpp_pb = -vmec_data.vpp;

% interpoloate/compute on full radial mesh
vmec_data.pprime = NaN;
vmec_data.iotaprime = NaN; % (Going to use the halfmesh iota)

%vmec_data.phip_real_fullmesh = NaN;
for jj = 2:(vmec_data.ns-1)
    %vmec_data.phip_real_fullmesh(jj) = 0.5 * ...
    %    (vmec_data.phip_real(jj+1) + vmec_data.phip_real(jj));
    denom = 1.0 / (vmec_data.hs * vmec_data.phip_real_fullmesh(jj));
    vmec_data.pprime(jj) = (vmec_data.pres(jj+1) - vmec_data.pres(jj)) ...
        * denom;
    vmec_data.iotaprime(jj) = (vmec_data.iotas(jj+1) - vmec_data.iotas(jj)) ...
        * denom;
end
%vmec_data.vpp(1) = 2*vmec_data.vpp(2) - vmec_data.vpp(3);

vmec_data.pprime(vmec_data.ns ) = NaN;
vmec_data.iotaprime(vmec_data.ns ) = NaN;

mpp_pb = -vmec_data.pprime;
iotap_pb = vmec_data.iotaprime;

fa_jdotb_pb = vmec_data.jdotb;
fa_bdotb_pb = vmec_data.bdotb;

metric1a_pb = mpp_pb .* mVpp_pb;
metric1b_pb = iotap_pb .* fa_jdotb_pb ./ fa_bdotb_pb;
metric1_pb = metric1a_pb + metric1b_pb;



if make_well_plot
    figure
    ii = 0;
    maxii = 8;
    sqr_layout_size = ceil(sqrt(maxii));
    
    
    ii = ii + 1;
    subplot(sqr_layout_size, sqr_layout_size,ii);
    plot(s,  vmec_data.phip_real, '+--');
    xlabel('s'); ylabel('phip\_real');
    axis tight;
    
    ii = ii + 1;
    subplot(sqr_layout_size, sqr_layout_size,ii);
    plot(s,  vmec_data.vp_real, '+');
    hold on;
    plot(s, vmec_data.vp_real_fit, '-');
    xlabel('s'); ylabel('vp\_real');
    axis tight;
    
    ii = ii + 1;
    subplot(sqr_layout_size, sqr_layout_size,ii);
    plot(s,  vmec_data.phip_real_fullmesh, '+--');
    xlabel('s'); ylabel('phip\_real\_fullmesh');
    axis tight;
    
    ii = ii + 1;
    subplot(sqr_layout_size, sqr_layout_size,ii);
    plot(s,  vmec_data.vpp, '+--');
    xlabel('s'); ylabel('vpp');
    axis tight;
    
    ii = ii + 1;
    subplot(sqr_layout_size, sqr_layout_size,ii);
    plot(s,  vmec_data.vp, '+--');
    xlabel('s'); ylabel('vp');
    axis tight;
    
    ii = ii + 1;
    subplot(sqr_layout_size, sqr_layout_size,ii);
    plot(s,  (vmec_data.vp(2) - vmec_data.vp) / vmec_data.vp(2), '+--');
    xlabel('s'); ylabel('well from vp');
    axis tight;
    
    ii = ii + 1;
    subplot(sqr_layout_size, sqr_layout_size,ii);
    plot(s, vmec_data.welldepth, '+');
    xlabel('s'); ylabel('well from vp\_real (extrap to 0)');
    axis tight;
    
    ii = ii + 1;
    subplot(sqr_layout_size, sqr_layout_size,ii);
    plot(s, vmec_data.welldepth, '+');
    xlabel('s'); ylabel('Well');
    axis tight;
    this_axis = axis;
    hold on
    
    plot(s,  (vmec_data.vp(2) - vmec_data.vp) / vmec_data.vp(2), '+--');
    axis(this_axis);
    legend('from vp\_real (extrap to 0)', 'from vp');
    
    figure
    ii = 0;
    maxii = 8;
    
    ii = ii + 1;
    subplot(sqr_layout_size, sqr_layout_size,ii);
    plot(rho,  vmec_data.phip_real, '+--');
    xlabel('rho'); ylabel('phip\_real');
    axis tight;
    
    ii = ii + 1;
    subplot(sqr_layout_size, sqr_layout_size,ii);
    plot(rho,  vmec_data.vp_real, '+');
    hold on;
    plot(rho, vmec_data.vp_real_fit, '-');
    xlabel('rho'); ylabel('vp\_real');
    axis tight;
    
    ii = ii + 1;
    subplot(sqr_layout_size, sqr_layout_size,ii);
    plot(rho,  vmec_data.phip_real_fullmesh, '+--');
    xlabel('rho'); ylabel('phip\_real\_fullmesh');
    axis tight;
    
    ii = ii + 1;
    subplot(sqr_layout_size, sqr_layout_size,ii);
    plot(rho,  vmec_data.vpp, '+--');
    xlabel('rho'); ylabel('vpp');
    axis tight;
    
    ii = ii + 1;
    subplot(sqr_layout_size, sqr_layout_size,ii);
    plot(rho,  vmec_data.vp, '+--');
    xlabel('rho'); ylabel('vp');
    axis tight;
    
    ii = ii + 1;
    subplot(sqr_layout_size, sqr_layout_size,ii);
    plot(rho,  (vmec_data.vp(2) - vmec_data.vp) / vmec_data.vp(2), '+--');
    xlabel('rho'); ylabel('well from vp');
    axis tight;
    
    ii = ii + 1;
    subplot(sqr_layout_size, sqr_layout_size,ii);
    plot(rho, vmec_data.welldepth, '+');
    xlabel('rho'); ylabel('well from vp\_real (extrap to 0)');
    axis tight;
    
    ii = ii + 1;
    subplot(sqr_layout_size, sqr_layout_size,ii);
    plot(rho, vmec_data.welldepth, '+');
    xlabel('rho'); ylabel('Well');
    axis tight;
    this_axis = axis;
    hold on
    plot(rho,  (vmec_data.vp(2) - vmec_data.vp) / vmec_data.vp(2), '+--');
    axis(this_axis);
    legend('from vp\_real (extrap to 0)', 'from vp');
end

if make_pb_plot
    
    figure
    ii = 0;
    maxii = 9;
    
    ii = ii + 1;
    subplot(sqr_layout_size, sqr_layout_size,ii);
    
    plot(s,  vmec_data.phi, '+--');
    xlabel('s'); title('tor flux');
    axis tight;
    grid on
    
    
    ii = ii + 1;
    subplot(sqr_layout_size, sqr_layout_size,ii);
    
    plot(s,  vmec_data.presf, '+--');
    xlabel('s'); title('presf');
    axis tight;
    grid on
    
    ii = ii + 1;
    subplot(sqr_layout_size, sqr_layout_size,ii);
    
    plot(s,  vmec_data.iotaf, '+--');
    xlabel('s'); title('iotaf');
    axis tight;
    grid on
    
    
    ii = ii + 1;
    subplot(sqr_layout_size, sqr_layout_size,ii);
    
    plot(s,  vmec_data.pprime, '+--');
    xlabel('s'); title('pprime');
    axis tight;
    grid on
    
    ii = ii + 1;
    subplot(sqr_layout_size, sqr_layout_size,ii);
    plot(s,  vmec_data.vpp, '+--');
    xlabel('s'); title('vpp');
    axis tight;
    grid on
    
    ii = ii + 1;
    subplot(sqr_layout_size, sqr_layout_size,ii);
    plot(s,  vmec_data.iotaprime, '+--');
    xlabel('s'); title('iotaprime');
    axis tight;
    grid on
    
    ii = ii + 1;
    subplot(sqr_layout_size, sqr_layout_size,ii);
    plot(s,  vmec_data.jdotb, '+--');
    xlabel('s'); title('jdotb');
    axis tight;
    grid on
    
    ii = ii + 1;
    subplot(sqr_layout_size, sqr_layout_size,ii);
    plot(s,  vmec_data.bdotb, '+--');
    xlabel('s'); title('bdotb');
    axis tight;
    grid on
    
    ii = ii + 1;
    subplot(sqr_layout_size, sqr_layout_size,ii);
    plot(s,  metric1_pb, 'o-');
    hold on;
    plot(s,  5*metric1a_pb, '+:');
    plot(s,  metric1b_pb, 'x--');
    legend('Left Side', '5x Part 1', 'Part 2');
    xlabel('s');
    title('Equation 2');
    axis tight;
    grid on
    
    
    figure
    ii = 0;
    maxii = 9;
    
    ii = ii + 1;
    subplot(sqr_layout_size, sqr_layout_size,ii);
    
    plot(rho,  vmec_data.phi, '+--');
    xlabel('rho'); title('tor flux');
    axis tight;
    grid on
    
    ii = ii + 1;
    subplot(sqr_layout_size, sqr_layout_size,ii);
    
    plot(rho,  vmec_data.presf, '+--');
    xlabel('rho'); title('presf');
    axis tight;
    grid on
    
    ii = ii + 1;
    subplot(sqr_layout_size, sqr_layout_size,ii);
    
    plot(rho,  vmec_data.iotaf, '+--');
    xlabel('rho'); title('iotaf');
    axis tight;
    grid on
    
    ii = ii + 1;
    subplot(sqr_layout_size, sqr_layout_size,ii);
    
    plot(rho,  vmec_data.pprime, '+--');
    xlabel('rho'); title('pprime');
    axis tight;
    grid on
    
    ii = ii + 1;
    subplot(sqr_layout_size, sqr_layout_size,ii);
    plot(rho,  vmec_data.vpp, '+--');
    xlabel('rho'); title('vpp');
    axis tight;
    grid on
    
    ii = ii + 1;
    subplot(sqr_layout_size, sqr_layout_size,ii);
    plot(rho,  vmec_data.iotaprime, '+--');
    xlabel('rho'); title('iotaprime');
    axis tight;
    grid on
    
    ii = ii + 1;
    subplot(sqr_layout_size, sqr_layout_size,ii);
    plot(rho,  vmec_data.jdotb, '+--');
    xlabel('rho'); title('jdotb');
    axis tight;
    grid on
    
    ii = ii + 1;
    subplot(sqr_layout_size, sqr_layout_size,ii);
    plot(rho,  vmec_data.bdotb, '+--');
    xlabel('rho'); title('bdotb');
    axis tight;
    grid on
    
    ii = ii + 1;
    subplot(sqr_layout_size, sqr_layout_size,ii);
    plot(rho,  metric1_pb, 'o-');
    hold on;
    plot(rho,  5*metric1a_pb, '+:');
    plot(rho,  metric1b_pb, 'x--');
    legend('Left Side', '5x Part 1', 'Part 2');
    xlabel('rho'); 
    title('Equation 2');
    axis tight;
    grid on
    
end

%keyboard







%%%   
%%%   
%%%   function vmec_data = load_vmec(varargin)
%%%   
%%%   <<<<<<< HEAD
%%%   
%%%   options = struct('vmec_woutfile_ext', '', ...
%%%       'make_plots', 0, ...
%%%       'debug_level', 0);
%%%   
%%%   if (~isempty(varargin))
%%%       options = parse_input_options(options, varargin{:});
%%%   end
%%%   
%%%   
%%%   
%%%   len_wout = length(options.vmec_woutfile_ext);
%%%   if (len_wout >= 3)
%%%       last3 = [(len_wout-2):len_wout];
%%%       
%%%       if (strcmpi(options.vmec_woutfile_ext(last3), '.nc') )
%%%           vmec_woutfile = options.vmec_woutfile_ext;
%%%       else
%%%           vmec_woutfile = ['wout_' options.vmec_woutfile_ext '.nc'];
%%%       end
%%%   =======
%%%   make_well_plot = false;
%%%   make_pb_plot = false;
%%%   
%%%   len_wout = length(vmec_woutfile_ext);
%%%   last3 = [(len_wout-2):len_wout];
%%%   
%%%   if (strcmpi(vmec_woutfile_ext(last3), '.nc') )
%%%       vmec_woutfile = vmec_woutfile_ext;
%%%   >>>>>>> 90fc100ac0685d04da2c0e3f33904ba2f34925be
%%%   else
%%%       vmec_woutfile = ['wout_' options.vmec_woutfile_ext '.nc'];
%%%   end
%%%   
%%%   vmec_data = read_vmec(vmec_woutfile);
%%%   
%%%   
%%%   
%%%   % Use vmec info to extrapolate the magnetic well depth
%%%   
%%%   % calculate the 'real' Vp, Volume and Well depth (see mercier.f)
%%%   % for 'typical' vmec run, the phips array is on the 'half mesh', and its
%%%   % first value is '0' because the value on-axis has not been computed by
%%%   % VMEC.
%%%   vmec_data.phip_real = 2*pi*2*pi * vmec_data.phips;
%%%   vmec_data.vp_real = NaN * vmec_data.phip_real;
%%%   
%%%   vmec_data.vp_real(2:end) = double(vmec_data.signgs) * ...
%%%       (4*pi*pi)*vmec_data.vp(2:end) ./ vmec_data.phip_real(2:end);
%%%   
%%%   
%%%   % interpoloate/compute on full radial mesh
%%%   % vmec_data.phip_real_fullmesh(vmec_data.ns) = NaN;
%%%   % vmec_data.vpp(vmec_data.ns) = NaN;
%%%   vmec_data.nsm1 = vmec_data.ns - 1;
%%%   vmec_data.hs = 1.0 / vmec_data.nsm1;
%%%   
%%%   % interpoloate/compute on full radial mesh
%%%   vmec_data.vpp = NaN;
%%%   vmec_data.phip_real_fullmesh = NaN;
%%%   for jj = 2:(vmec_data.ns-1)
%%%       vmec_data.phip_real_fullmesh(jj) = 0.5 * ...
%%%           (vmec_data.phip_real(jj+1) + vmec_data.phip_real(jj));
%%%       denom = 1.0 / (vmec_data.hs * vmec_data.phip_real_fullmesh(jj));
%%%       vmec_data.vpp(jj) = (vmec_data.vp_real(jj+1) - vmec_data.vp_real(jj)) ...
%%%           * denom;
%%%   end
%%%   %vmec_data.vpp(1) = 2*vmec_data.vpp(2) - vmec_data.vpp(3);
%%%   
%%%   % Now we assing the last
%%%   vmec_data.phip_real_fullmesh(vmec_data.ns ) = NaN;
%%%   vmec_data.vpp(vmec_data.ns ) = NaN;
%%%   %vmec_data.phip_real_fullmesh(vmec_data.ns ) = NaN;
%%%   %vmec_data.vpp(vmec_data.ns ) = NaN;
%%%   
%%%   
%%%   
%%%   
%%%   <<<<<<< HEAD
%%%   
%%%   =======
%%%   >>>>>>> 90fc100ac0685d04da2c0e3f33904ba2f34925be
%%%   
%%%   s = vmec_data.phi ./ vmec_data.phi(end);
%%%   rho = sqrt(s);
%%%   
%%%   <<<<<<< HEAD
%%%   vmec_data.welldepth =  (vmec_data.vp_real(2) - vmec_data.vp_real) / vmec_data.vp_real(2);
%%%   
%%%   
%%%   if options.make_plots
%%%       figure
%%%       ii = 0;
%%%       maxii = 7;
%%%       
%%%       ii = ii + 1;
%%%       subplot(maxii,1,ii);
%%%       plot(rho,  vmec_data.phip_real, '+--');
%%%       xlabel('rho'); ylabel('phip_real');
%%%       axis tight;
%%%       
%%%       ii = ii + 1;
%%%       subplot(maxii,1,ii);
%%%       plot(rho,  vmec_data.vp_real, '+--');
%%%       xlabel('rho'); ylabel('vp_real');
%%%       axis tight;
%%%       
%%%       ii = ii + 1;
%%%       subplot(maxii,1,ii);
%%%       plot(rho,  vmec_data.phip_real2, '+--');
%%%       xlabel('rho'); ylabel('phip_real2');
%%%       axis tight;
%%%       
%%%       ii = ii + 1;
%%%       subplot(maxii,1,ii);
%%%   =======
%%%   
%%%   % make a parabolic fit to the first few indices to extrapoloate to s=0
%%%   % for some variables
%%%   temp_indices = [2:10];
%%%   disp(['<----Extrapolatiing VP_REAL to axis based on surfaces ' num2str(temp_indices(1)) ' to ' num2str(temp_indices(end))]);
%%%   x_temp = s(temp_indices);
%%%   y_temp = vmec_data.vp_real(temp_indices);
%%%   fit_temp = polyfit(x_temp, y_temp, 2);
%%%   x0_temp = polyval(fit_temp, 0);
%%%   vmec_data.vp_real_fit = vmec_data.vp_real;
%%%   vmec_data.vp_real_fit(1) = x0_temp;
%%%   
%%%   %vmec_data.welldepth =  (vmec_data.vp_real(2) - vmec_data.vp_real) / vmec_data.vp_real(2);
%%%   vmec_data.welldepth =  (vmec_data.vp_real_fit(1) - vmec_data.vp_real_fit) / vmec_data.vp_real_fit(1);
%%%   
%%%   
%%%   % Peeling-Ballooning
%%%   % (-p')(-V'') + (iota') <J_{||} B> / <B^2> > 0
%%%   % ' == derivative w.r.t. toroidal flux (psi)
%%%   % V'' = 1/(4*pi^2)  d^2 (V) / d (psi)^2
%%%   % V'' and p' are not available at s=0, s=1, because phip_real_fullmesh is
%%%   % not on the full grid.
%%%   mVpp_pb = -vmec_data.vpp;
%%%   
%%%   % interpoloate/compute on full radial mesh
%%%   vmec_data.pprime = NaN;
%%%   vmec_data.iotaprime = NaN; % (Going to use the halfmesh iota)
%%%   
%%%   %vmec_data.phip_real_fullmesh = NaN;
%%%   for jj = 2:(vmec_data.ns-1)
%%%       %vmec_data.phip_real_fullmesh(jj) = 0.5 * ...
%%%       %    (vmec_data.phip_real(jj+1) + vmec_data.phip_real(jj));
%%%       denom = 1.0 / (vmec_data.hs * vmec_data.phip_real_fullmesh(jj));
%%%       vmec_data.pprime(jj) = (vmec_data.pres(jj+1) - vmec_data.pres(jj)) ...
%%%           * denom;
%%%       vmec_data.iotaprime(jj) = (vmec_data.iotas(jj+1) - vmec_data.iotas(jj)) ...
%%%           * denom;
%%%   end
%%%   %vmec_data.vpp(1) = 2*vmec_data.vpp(2) - vmec_data.vpp(3);
%%%   
%%%   vmec_data.pprime(vmec_data.ns ) = NaN;
%%%   vmec_data.iotaprime(vmec_data.ns ) = NaN;
%%%   
%%%   mpp_pb = -vmec_data.pprime;
%%%   iotap_pb = vmec_data.iotaprime;
%%%   
%%%   fa_jdotb_pb = vmec_data.jdotb;
%%%   fa_bdotb_pb = vmec_data.bdotb;
%%%   
%%%   metric1a_pb = mpp_pb .* mVpp_pb;
%%%   metric1b_pb = iotap_pb .* fa_jdotb_pb ./ fa_bdotb_pb;
%%%   metric1_pb = metric1a_pb + metric1b_pb;
%%%   
%%%   
%%%   
%%%   if make_well_plot
%%%       figure
%%%       ii = 0;
%%%       maxii = 8;
%%%       sqr_layout_size = ceil(sqrt(maxii));
%%%       
%%%       
%%%       ii = ii + 1;
%%%       subplot(sqr_layout_size, sqr_layout_size,ii);
%%%       plot(s,  vmec_data.phip_real, '+--');
%%%       xlabel('s'); ylabel('phip\_real');
%%%       axis tight;
%%%       
%%%       ii = ii + 1;
%%%       subplot(sqr_layout_size, sqr_layout_size,ii);
%%%       plot(s,  vmec_data.vp_real, '+');
%%%       hold on;
%%%       plot(s, vmec_data.vp_real_fit, '-');
%%%       xlabel('s'); ylabel('vp\_real');
%%%       axis tight;
%%%       
%%%       ii = ii + 1;
%%%       subplot(sqr_layout_size, sqr_layout_size,ii);
%%%       plot(s,  vmec_data.phip_real_fullmesh, '+--');
%%%       xlabel('s'); ylabel('phip\_real\_fullmesh');
%%%       axis tight;
%%%       
%%%       ii = ii + 1;
%%%       subplot(sqr_layout_size, sqr_layout_size,ii);
%%%       plot(s,  vmec_data.vpp, '+--');
%%%       xlabel('s'); ylabel('vpp');
%%%       axis tight;
%%%       
%%%       ii = ii + 1;
%%%       subplot(sqr_layout_size, sqr_layout_size,ii);
%%%       plot(s,  vmec_data.vp, '+--');
%%%       xlabel('s'); ylabel('vp');
%%%       axis tight;
%%%       
%%%       ii = ii + 1;
%%%       subplot(sqr_layout_size, sqr_layout_size,ii);
%%%       plot(s,  (vmec_data.vp(2) - vmec_data.vp) / vmec_data.vp(2), '+--');
%%%       xlabel('s'); ylabel('well from vp');
%%%       axis tight;
%%%       
%%%       ii = ii + 1;
%%%       subplot(sqr_layout_size, sqr_layout_size,ii);
%%%       plot(s, vmec_data.welldepth, '+');
%%%       xlabel('s'); ylabel('well from vp\_real (extrap to 0)');
%%%       axis tight;
%%%       
%%%       ii = ii + 1;
%%%       subplot(sqr_layout_size, sqr_layout_size,ii);
%%%       plot(s, vmec_data.welldepth, '+');
%%%       xlabel('s'); ylabel('Well');
%%%       axis tight;
%%%       this_axis = axis;
%%%       hold on
%%%       
%%%       plot(s,  (vmec_data.vp(2) - vmec_data.vp) / vmec_data.vp(2), '+--');
%%%       axis(this_axis);
%%%       legend('from vp\_real (extrap to 0)', 'from vp');
%%%       
%%%       figure
%%%       ii = 0;
%%%       maxii = 8;
%%%       
%%%       ii = ii + 1;
%%%       subplot(sqr_layout_size, sqr_layout_size,ii);
%%%       plot(rho,  vmec_data.phip_real, '+--');
%%%       xlabel('rho'); ylabel('phip\_real');
%%%       axis tight;
%%%       
%%%       ii = ii + 1;
%%%       subplot(sqr_layout_size, sqr_layout_size,ii);
%%%       plot(rho,  vmec_data.vp_real, '+');
%%%       hold on;
%%%       plot(rho, vmec_data.vp_real_fit, '-');
%%%       xlabel('rho'); ylabel('vp\_real');
%%%       axis tight;
%%%       
%%%       ii = ii + 1;
%%%       subplot(sqr_layout_size, sqr_layout_size,ii);
%%%       plot(rho,  vmec_data.phip_real_fullmesh, '+--');
%%%       xlabel('rho'); ylabel('phip\_real\_fullmesh');
%%%       axis tight;
%%%       
%%%       ii = ii + 1;
%%%       subplot(sqr_layout_size, sqr_layout_size,ii);
%%%   >>>>>>> 90fc100ac0685d04da2c0e3f33904ba2f34925be
%%%       plot(rho,  vmec_data.vpp, '+--');
%%%       xlabel('rho'); ylabel('vpp');
%%%       axis tight;
%%%       
%%%       ii = ii + 1;
%%%   <<<<<<< HEAD
%%%       subplot(maxii,1,ii);
%%%   =======
%%%       subplot(sqr_layout_size, sqr_layout_size,ii);
%%%   >>>>>>> 90fc100ac0685d04da2c0e3f33904ba2f34925be
%%%       plot(rho,  vmec_data.vp, '+--');
%%%       xlabel('rho'); ylabel('vp');
%%%       axis tight;
%%%       
%%%       ii = ii + 1;
%%%   <<<<<<< HEAD
%%%       subplot(maxii,1,ii);
%%%   =======
%%%       subplot(sqr_layout_size, sqr_layout_size,ii);
%%%   >>>>>>> 90fc100ac0685d04da2c0e3f33904ba2f34925be
%%%       plot(rho,  (vmec_data.vp(2) - vmec_data.vp) / vmec_data.vp(2), '+--');
%%%       xlabel('rho'); ylabel('well from vp');
%%%       axis tight;
%%%       
%%%       ii = ii + 1;
%%%   <<<<<<< HEAD
%%%       subplot(maxii,1,ii);
%%%       plot(rho,  vmec_data.welldepth, '+--');
%%%       xlabel('rho'); ylabel('well from vp_real');
%%%       axis tight;
%%%       
%%%   end
%%%   =======
%%%       subplot(sqr_layout_size, sqr_layout_size,ii);
%%%       plot(rho, vmec_data.welldepth, '+');
%%%       xlabel('rho'); ylabel('well from vp\_real (extrap to 0)');
%%%       axis tight;
%%%       
%%%       ii = ii + 1;
%%%       subplot(sqr_layout_size, sqr_layout_size,ii);
%%%       plot(rho, vmec_data.welldepth, '+');
%%%       xlabel('rho'); ylabel('Well');
%%%       axis tight;
%%%       this_axis = axis;
%%%       hold on
%%%       plot(rho,  (vmec_data.vp(2) - vmec_data.vp) / vmec_data.vp(2), '+--');
%%%       axis(this_axis);
%%%       legend('from vp\_real (extrap to 0)', 'from vp');
%%%   end
%%%   
%%%   if make_pb_plot
%%%       
%%%       figure
%%%       ii = 0;
%%%       maxii = 9;
%%%       
%%%       ii = ii + 1;
%%%       subplot(sqr_layout_size, sqr_layout_size,ii);
%%%       
%%%       plot(s,  vmec_data.phi, '+--');
%%%       xlabel('s'); title('tor flux');
%%%       axis tight;
%%%       grid on
%%%       
%%%       
%%%       ii = ii + 1;
%%%       subplot(sqr_layout_size, sqr_layout_size,ii);
%%%       
%%%       plot(s,  vmec_data.presf, '+--');
%%%       xlabel('s'); title('presf');
%%%       axis tight;
%%%       grid on
%%%       
%%%       ii = ii + 1;
%%%       subplot(sqr_layout_size, sqr_layout_size,ii);
%%%       
%%%       plot(s,  vmec_data.iotaf, '+--');
%%%       xlabel('s'); title('iotaf');
%%%       axis tight;
%%%       grid on
%%%       
%%%       
%%%       ii = ii + 1;
%%%       subplot(sqr_layout_size, sqr_layout_size,ii);
%%%       
%%%       plot(s,  vmec_data.pprime, '+--');
%%%       xlabel('s'); title('pprime');
%%%       axis tight;
%%%       grid on
%%%       
%%%       ii = ii + 1;
%%%       subplot(sqr_layout_size, sqr_layout_size,ii);
%%%       plot(s,  vmec_data.vpp, '+--');
%%%       xlabel('s'); title('vpp');
%%%       axis tight;
%%%       grid on
%%%       
%%%       ii = ii + 1;
%%%       subplot(sqr_layout_size, sqr_layout_size,ii);
%%%       plot(s,  vmec_data.iotaprime, '+--');
%%%       xlabel('s'); title('iotaprime');
%%%       axis tight;
%%%       grid on
%%%       
%%%       ii = ii + 1;
%%%       subplot(sqr_layout_size, sqr_layout_size,ii);
%%%       plot(s,  vmec_data.jdotb, '+--');
%%%       xlabel('s'); title('jdotb');
%%%       axis tight;
%%%       grid on
%%%       
%%%       ii = ii + 1;
%%%       subplot(sqr_layout_size, sqr_layout_size,ii);
%%%       plot(s,  vmec_data.bdotb, '+--');
%%%       xlabel('s'); title('bdotb');
%%%       axis tight;
%%%       grid on
%%%       
%%%       ii = ii + 1;
%%%       subplot(sqr_layout_size, sqr_layout_size,ii);
%%%       plot(s,  metric1_pb, 'o-');
%%%       hold on;
%%%       plot(s,  5*metric1a_pb, '+:');
%%%       plot(s,  metric1b_pb, 'x--');
%%%       legend('Left Side', '5x Part 1', 'Part 2');
%%%       xlabel('s');
%%%       title('Equation 2');
%%%       axis tight;
%%%       grid on
%%%       
%%%       
%%%       figure
%%%       ii = 0;
%%%       maxii = 9;
%%%       
%%%       ii = ii + 1;
%%%       subplot(sqr_layout_size, sqr_layout_size,ii);
%%%       
%%%       plot(rho,  vmec_data.phi, '+--');
%%%       xlabel('rho'); title('tor flux');
%%%       axis tight;
%%%       grid on
%%%       
%%%       ii = ii + 1;
%%%       subplot(sqr_layout_size, sqr_layout_size,ii);
%%%       
%%%       plot(rho,  vmec_data.presf, '+--');
%%%       xlabel('rho'); title('presf');
%%%       axis tight;
%%%       grid on
%%%       
%%%       ii = ii + 1;
%%%       subplot(sqr_layout_size, sqr_layout_size,ii);
%%%       
%%%       plot(rho,  vmec_data.iotaf, '+--');
%%%       xlabel('rho'); title('iotaf');
%%%       axis tight;
%%%       grid on
%%%       
%%%       ii = ii + 1;
%%%       subplot(sqr_layout_size, sqr_layout_size,ii);
%%%       
%%%       plot(rho,  vmec_data.pprime, '+--');
%%%       xlabel('rho'); title('pprime');
%%%       axis tight;
%%%       grid on
%%%       
%%%       ii = ii + 1;
%%%       subplot(sqr_layout_size, sqr_layout_size,ii);
%%%       plot(rho,  vmec_data.vpp, '+--');
%%%       xlabel('rho'); title('vpp');
%%%       axis tight;
%%%       grid on
%%%       
%%%       ii = ii + 1;
%%%       subplot(sqr_layout_size, sqr_layout_size,ii);
%%%       plot(rho,  vmec_data.iotaprime, '+--');
%%%       xlabel('rho'); title('iotaprime');
%%%       axis tight;
%%%       grid on
%%%       
%%%       ii = ii + 1;
%%%       subplot(sqr_layout_size, sqr_layout_size,ii);
%%%       plot(rho,  vmec_data.jdotb, '+--');
%%%       xlabel('rho'); title('jdotb');
%%%       axis tight;
%%%       grid on
%%%       
%%%       ii = ii + 1;
%%%       subplot(sqr_layout_size, sqr_layout_size,ii);
%%%       plot(rho,  vmec_data.bdotb, '+--');
%%%       xlabel('rho'); title('bdotb');
%%%       axis tight;
%%%       grid on
%%%       
%%%       ii = ii + 1;
%%%       subplot(sqr_layout_size, sqr_layout_size,ii);
%%%       plot(rho,  metric1_pb, 'o-');
%%%       hold on;
%%%       plot(rho,  5*metric1a_pb, '+:');
%%%       plot(rho,  metric1b_pb, 'x--');
%%%       legend('Left Side', '5x Part 1', 'Part 2');
%%%       xlabel('rho'); 
%%%       title('Equation 2');
%%%       axis tight;
%%%       grid on
%%%       
%%%   end
%%%   
%%%   %keyboard
%%%   >>>>>>> 90fc100ac0685d04da2c0e3f33904ba2f34925be
%%%   
%%%   
%%%   


