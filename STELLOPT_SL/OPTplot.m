function varargout = OPTplot(varargin)
% OPTPLOT GUI for interactive STELLOPT plots.
%   OPTPLOT is a MATLAB Graphical User Interface (GUI) for
%   visualizing the data stored in the STELLOPT output file (stellopt.*).
%   It allows the user to plot various quantities in the output file.
%
%   See also read_vmec, read_stellopt.
%
%   Example:
%      OPTplot
%
%   Written by: Samuel Lazerson (lazerson@pppl.gov)
%   Version:    2.5
%   Date:       03/16/2016

% Edit the above text to modify the response to help OPTplot

% Last Modified by GUIDE v2.5 08-Mar-2018 16:50:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @OPTplot_OpeningFcn, ...
                   'gui_OutputFcn',  @OPTplot_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before OPTplot is made visible.
function OPTplot_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to OPTplot (see VARARGIN)

% Choose default command line output for OPTplot
handles.phi0 = 0.0;
handles.output = hObject;
filename = dir('stellopt.*');
if (isempty(filename))
    return;
end
handles.chi_types={'ASPECT','BETA','CURTOR','EXTCUR','SEPARATRIX',...
    'PHIEDGE','RBTOR','R0','Z0','VOLUME','WP','KAPPA',...
    'B_PROBES','FARADAY','FLUXLOOPS','SEGROG','MSE',...
    'NE','NELINE','TE','TELINE','TI','TILINE',...
    'XICS','XICS_BRIGHT','SXR','VPHI',...
    'IOTA','BALLOON','BOOTSTRAP','DKES','HELICITY','HELICITY_FULL',...
    'KINK','ORBIT','JDOTB','J_STAR','NEO','TXPORT','ECEREFLECT',...
    'REGCOIL_CHI2_B', 'PTSM3D'};
handles.data = read_stellopt(filename(1).name);  % Read the file
handles.plot_types={'Chi-Squared';'Chi-Stacked'};
for i=1:length(handles.chi_types)
    if isfield(handles.data,[handles.chi_types{i} '_chisq'])
        handles.plot_types = [handles.plot_types; [handles.chi_types{i} '_chisq']];
    end
end
handles.plot_types = [handles.plot_types; '-----SPECIAL-----'];
if isfield(handles.data,'BALLOON_chisq')
    handles.plot_types = [handles.plot_types; 'BALLOON_evolution'];
end
if isfield(handles.data,'KINK_chisq')
    handles.plot_types = [handles.plot_types; 'KINK_evolution'];
end
if isfield(handles.data,'ORBIT_chisq')
    handles.plot_types = [handles.plot_types; 'ORBIT_evolution'];
end
if isfield(handles.data,'NEO_chisq')
    handles.plot_types = [handles.plot_types; 'NEO_evolution'];
end
if isfield(handles.data,'HELICITY_chisq')
    handles.plot_types = [handles.plot_types; 'HELICITY_evolution'];
end
if isfield(handles.data,'HELICITY_FULL_chisq')
    handles.plot_types = [handles.plot_types; 'HELICITY_evolution'];
    handles.plot_types = [handles.plot_types; 'HELICITY_B'];
end
if isfield(handles.data,'TXPORT_chisq')
    handles.plot_types = [handles.plot_types; 'TXPORT_evolution'];
end
if isfield(handles.data,'B_PROBES_chisq')
    handles.plot_types = [handles.plot_types; 'B_PROBE_evolution'];
end
if isfield(handles.data,'FLUXLOOPS_chisq')
    handles.plot_types = [handles.plot_types; 'FLUXLOOP_evolution'];
end
if isfield(handles.data,'SEGROG_chisq')
    handles.plot_types = [handles.plot_types; 'SEGROG_evolution'];
end
if isfield(handles.data,'NE_chisq')
    handles.plot_types = [handles.plot_types; 'NE_evolution'];
end
if isfield(handles.data,'PRESS_chisq')
    handles.plot_types = [handles.plot_types; 'PRESS_evolution'];
end
if isfield(handles.data,'TE_chisq')
    handles.plot_types = [handles.plot_types; 'TE_evolution'];
end
if isfield(handles.data,'TI_chisq')
    handles.plot_types = [handles.plot_types; 'TI_evolution'];
end
if isfield(handles.data,'NELINE_chisq')
    handles.plot_types = [handles.plot_types; 'NELINE_evolution'];
end
if isfield(handles.data,'TELINE_chisq')
    handles.plot_types = [handles.plot_types; 'TELINE_evolution'];
end
if isfield(handles.data,'TILINE_chisq')
    handles.plot_types = [handles.plot_types; 'TILINE_evolution'];
end
if isfield(handles.data,'XICS_chisq')
    handles.plot_types = [handles.plot_types; 'XICS_evolution'];
end
if isfield(handles.data,'XICS_BRIGHT_chisq')
    handles.plot_types = [handles.plot_types; 'XICS_BRIGHT_evolution'];
end
if isfield(handles.data,'IOTA_chisq')
    handles.plot_types = [handles.plot_types; 'IOTA_evolution'];
end
if isfield(handles.data,'ECEREFLECT_chisq')
    handles.plot_types = [handles.plot_types; 'ECE_evolution'];
    handles.plot_types = [handles.plot_types; 'ECE_xmode'];
    handles.plot_types = [handles.plot_types; 'ECE_omode'];
end
if isfield(handles.data,'MSE_chisq')
    handles.plot_types = [handles.plot_types; 'MSE_evolution'];
end
if isfield(handles.data,'SXR_chisq')
    handles.plot_types = [handles.plot_types; 'SXR_evolution'];
end
% Here we create the ability to plot equilibria
filename = dir('wout*.*.nc');
if ~isempty(filename)
    handles.plot_types = [handles.plot_types; '-----WOUT-----'];
    handles.plot_types = [handles.plot_types; 'FLUX0'];
    handles.plot_types = [handles.plot_types; 'FLUXPI'];
    handles.plot_types = [handles.plot_types; 'PRESSURE'];
    handles.plot_types = [handles.plot_types; 'CURRENT'];
    handles.plot_types = [handles.plot_types; 'IOTA'];
    handles.plot_types = [handles.plot_types; '<J*B>'];
    handles.plot_types = [handles.plot_types; 'CURTOR'];
    handles.plot_types = [handles.plot_types; 'PHIEDGE'];
    handles.plot_types = [handles.plot_types; 'EXTCUR'];
end
filename_prof = dir('tprof.*.*');
if ~isempty(filename_prof)
    handles.plot_types = [handles.plot_types; '-----PROF-----'];
    handles.plot_types = [handles.plot_types; 'TE'];
    handles.plot_types = [handles.plot_types; 'TI'];
    handles.plot_types = [handles.plot_types; 'NE'];
    handles.plot_types = [handles.plot_types; 'ZEFF'];
end
filename_prof = dir('dprof.*.*');
if ~isempty(filename_prof)
    handles.plot_types = [handles.plot_types; 'EMIS_XICS'];
end
filename_prof = dir('gist_genet_*.*');
if ~isempty(filename_prof)
    handles.plot_types = [handles.plot_types; '-----GIST-----'];
    handles.plot_types = [handles.plot_types; 'G11'];
    handles.plot_types = [handles.plot_types; 'G12'];
    handles.plot_types = [handles.plot_types; 'G22'];
    handles.plot_types = [handles.plot_types; 'Bhat'];
    handles.plot_types = [handles.plot_types; '|JAC|'];
    handles.plot_types = [handles.plot_types; 'L1'];
    handles.plot_types = [handles.plot_types; 'L2'];
    handles.plot_types = [handles.plot_types; 'DBDT'];
    handles.plot_types = [handles.plot_types; 'KP1'];
end
filename_prof = dir('txport_out.*.*');
if ~isempty(filename_prof)
    handles.plot_types = [handles.plot_types; '-----TXPORT-----'];
    handles.plot_types = [handles.plot_types; 'PROXY'];
end
filename_prof = dir('terpsichore_16.*');
if ~isempty(filename_prof)
    handles.plot_types = [handles.plot_types; '-----TERPSICHORE-----'];
    handles.plot_types = [handles.plot_types; 'XI_00'];
    handles.plot_types = [handles.plot_types; 'XI_01'];
    handles.plot_types = [handles.plot_types; 'ETA_00'];
    handles.plot_types = [handles.plot_types; 'ETA_01'];
end
        
% Update the UI
set(handles.pulldownmenu,'String',handles.plot_types)
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes OPTplot wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = OPTplot_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in pulldownmenu.
function pulldownmenu_Callback(hObject, eventdata, handles)
% hObject    handle to pulldownmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pulldownmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pulldownmenu
update_plots(handles)

function update_plots(handles)
contents = cellstr(get(handles.pulldownmenu,'String'));
stemp=contents{get(handles.pulldownmenu,'Value')};
        set(gca,'YScale','linear');
cinitial='b';
cfinal='g';
marker='+';
legend off;
switch stemp
    case{'Chi-Stacked'}
        data_array=[];
        lab_array={};
        for i=1:length(handles.chi_types)
            if isfield(handles.data,[handles.chi_types{i} '_chisq'])
               chi_temp = handles.data.([handles.chi_types{i} '_chisq']);
               data_array = [data_array sum(chi_temp,2)];
               lab_array = [lab_array; strrep(handles.chi_types{i},'_','\_')];
            end
        end
        area(handles.data.iter,data_array);
        legend(lab_array);
        xlabel('ITERATION');
        ylabel('Total \chi^2');
        title('TOTAL CHISQ');
    case{'Chi-Squared'}
        data_array=[];
        lab_array={};
        for i=1:length(handles.chi_types)
            if isfield(handles.data,[handles.chi_types{i} '_chisq'])
               chi_temp = handles.data.([handles.chi_types{i} '_chisq']);
               data_array = [data_array sum(chi_temp,2)];
               lab_array = [lab_array; strrep(handles.chi_types{i},'_','\_')];
            end
        end
        f = (sum( ( abs(handles.data.TARGETS'-handles.data.VALS') ./ ...
                   abs(handles.data.SIGMAS')).^2 ));
        %data_array = [f' data_array];
        lab_array = ['TOTAL'; lab_array];
        plot(handles.data.iter,f','or','MarkerSize',18,'LineWidth',4.0);
        hold on; 
        plot(handles.data.iter,data_array,marker,'MarkerSize',18,'LineWidth',4.0); hold off;
        legend(lab_array);
        set(gca,'YScale','log');
        xlabel('ITERATION');
        ylabel('Total \chi^2');
        title('TOTAL CHISQ');
    case{'CURTOR_chisq','VOLUME_chisq','BETA_chisq','ASPECT_chisq',...
            'KAPPA_chisq','REGCOIL_CHI2_B_chisq', 'PTSM3D_chisq'}
        f = handles.data.(stemp);
        plot(handles.data.iter,f,marker,'MarkerSize',18,'LineWidth',4.0);
        set(gca,'YScale','log');
        xlabel('ITERATION');
        ylabel('\chi^2');
        tstr = stemp; tstr(strfind(tstr,'_'))=' '; tstr=upper(tstr);
        title(tstr);
    case{'TXPORT_chisq','B_PROBES_chisq','FLUXLOOPS_chisq',...
            'BALLOON_chisq','NEO_chisq','HELICITY_chisq',...
            'HELICITY_FULL_chisq','NELINE_chisq',...
            'TELINE_chisq','TILINE_chisq','SXR_chisq','ECEREFLECT_chisq',...
            'KINK_chisq','XICS_chisq','XICS_BRIGHT_chisq'}
        f = sum(handles.data.(stemp),2);
        plot(handles.data.iter,f,marker,'MarkerSize',18,'LineWidth',4.0);
        set(gca,'YScale','log');
        xlabel('ITERATION');
        ylabel('\chi^2');
        tstr = stemp; tstr(strfind(tstr,'_'))=' '; tstr=upper(tstr);
        title(tstr);
    case{'BALLOON_evolution'}
        plot(handles.data.BALLOON_k(:,:)',handles.data.BALLOON_grate(:,:)','k');
        hold on;
        plot(handles.data.BALLOON_k(1,:),handles.data.BALLOON_grate(1,:),cinitial,'LineWidth',2.0);
        plot(handles.data.BALLOON_k(end,:),handles.data.BALLOON_grate(end,:),cfinal,'LineWidth',2.0);
        plot(handles.data.BALLOON_k(1,:),handles.data.BALLOON_target(end,:),'r','LineWidth',2.0);
        hold off;
        xlabel('Surface Number');
        ylabel('Ballooning Stability');
        title('Ballooning Stability Evolution');
    case{'TXPORT_evolution'}
        plot(handles.data.TXPORT_s(:,:)',handles.data.TXPORT_equil(:,:)','k');
        hold on;
        plot(handles.data.TXPORT_s(1,:),handles.data.TXPORT_equil(1,:),cinitial,'LineWidth',2.0);
        plot(handles.data.TXPORT_s(end,:),handles.data.TXPORT_equil(end,:),cfinal,'LineWidth',2.0);
        hold off;
        xlim([0 1]);
        xlabel('Normalized Flux');
        ylabel('Proxy Function');
        title('Turbulent Transport Evolution');
    case{'ORBIT_evolution'}
        if (size(handles.data.ORBIT_s,2)==1)
            plot(handles.data.ORBIT_s(:,:),handles.data.ORBIT_equil(:,:),'k.');
            cinitial = [cinitial 'o'];
            cfinal = [cfinal 'o'];
        else
            plot(handles.data.ORBIT_s(:,:)',handles.data.ORBIT_equil(:,:)','k');
        end
        hold on;
        plot(handles.data.ORBIT_s(1,:),handles.data.ORBIT_equil(1,:),cinitial,'LineWidth',2.0);
        plot(handles.data.ORBIT_s(end,:),handles.data.ORBIT_equil(end,:),cfinal,'LineWidth',2.0);
        hold off;
        xlim([0 1]);
        xlabel('Normalized Flux');
        ylabel('Orbit Losses');
        title('Energetic Particle Losses');
    case{'NEO_evolution'}
        plot(handles.data.NEO_k(:,:)',handles.data.NEO_equil(:,:)','k');
        hold on;
        plot(handles.data.NEO_k(1,:),handles.data.NEO_equil(1,:),cinitial,'LineWidth',2.0);
        plot(handles.data.NEO_k(end,:),handles.data.NEO_equil(end,:),cfinal,'LineWidth',2.0);
        hold off;
        %xlim([0 1]);
        xlabel('Flux Surface');
        ylabel('Epsilon Effective');
        title('Neoclassical (NEO) Transport Evolution');
    case{'HELICITY_evolution'}
        if isfield(handles.data,'HELICITY2_target')
            plot(handles.data.HELICITY2_equil(:,:)','k');
            hold on;
            plot(handles.data.HELICITY2_equil(1,:),cinitial,'LineWidth',2.0);
            plot(handles.data.HELICITY2_equil(end,:),cfinal,'LineWidth',2.0);
        else
            plot(handles.data.HELICITY_equil(:,:)','k');
            hold on;
            plot(handles.data.HELICITY_equil(1,:),cinitial,'LineWidth',2.0);
            plot(handles.data.HELICITY_equil(end,:),cfinal,'LineWidth',2.0);
        end
        hold off;
        %xlim([0 1]);
        xlabel('');
        ylabel('Helicity');
        title('Boozer Spectrum Helicity Evolution');
    case{'HELICITY_B'}
        mmax = max(handles.data.HELICITY_FULL_m(1,:));
        nmax = max(handles.data.HELICITY_FULL_n(1,:));
        kmax = max(handles.data.HELICITY_FULL_k(1,:));
        nstep = size(handles.data.HELICITY_FULL_equil,1);
        for i=1:size(handles.data.HELICITY_FULL_equil,2)
            m=handles.data.HELICITY_FULL_m(1,i);
            n=handles.data.HELICITY_FULL_n(1,i);
            k=handles.data.HELICITY_FULL_k(1,i);
            b(:,k,m+1,n+nmax+1)=log10(abs(handles.data.HELICITY_FULL_equil(:,i)));
        end
        pixplot(0:mmax,-nmax:nmax,squeeze(b(nstep,kmax,:,:)))
        %plot(handles.data.HELICITY_equil(:,:)','k');
        %hold on;
        %plot(handles.data.HELICITY_equil(1,:),cinitial,'LineWidth',2.0);
        %plot(handles.data.HELICITY_equil(end,:),cfinal,'LineWidth',2.0);
        %hold off;
        %xlim([0 1]);
        xlabel('m');
        ylabel('n');
        title('Boozer Spectrum Helicity Evolution');
        colormap hot;
    case{'B_PROBE_evolution'}
        errorbar(handles.data.B_PROBES_target(1,:),handles.data.B_PROBES_sigma(1,:),'ok');
        hold on;
        plot(handles.data.B_PROBES_equil(:,:)','.k');
        plot(handles.data.B_PROBES_equil(1,:),'+b');
        plot(handles.data.B_PROBES_equil(end,:),'+r');
        hold off;
        xlabel('B-Probe Index');
        ylabel('Signal');
        title('B-Probe Reconstruction');
    case{'KINK_evolution'}
        errorbar(handles.data.KINK_target(1,:),handles.data.KINK_sigma(1,:),'ok');
        hold on;
        plot(handles.data.KINK_equil(:,:)','.k');
        plot(handles.data.KINK_equil(1,:),'+b');
        plot(handles.data.KINK_equil(end,:),'+r');
        hold off;
        xlabel('KINK Family');
        ylabel('Stability');
        title('KINK Evolution');
    case{'FLUXLOOP_evolution'}
        x_loop=1:size(handles.data.FLUXLOOPS_target,2);
        dex = handles.data.FLUXLOOPS_sigma(1,:) < 1.0E10;
        errorbar(x_loop(dex),handles.data.FLUXLOOPS_target(1,dex),handles.data.FLUXLOOPS_sigma(1,dex),'ok');
        hold on;
        plot(x_loop(dex),handles.data.FLUXLOOPS_equil(:,dex)','.k');
        plot(x_loop(dex),handles.data.FLUXLOOPS_equil(1,dex),'+b');
        plot(x_loop(dex),handles.data.FLUXLOOPS_equil(end,dex),'+g');
        hold off;
        xlabel('Fluxloop Index');
        ylabel('Signal');
        title('Fluxloop Reconstruction');
        ylim([min(ylim)-0.2*diff(ylim) max(ylim)+0.2*diff(ylim)]);
        xlim([min(xlim)-1 max(xlim)+1]);
    case{'SEGROG_evolution'}
        x_loop=1:size(handles.data.SEGROG_target,2);
        dex = handles.data.SEGROG_sigma(1,:) < 1.0E10;
        errorbar(x_loop(dex),handles.data.SEGROG_target(1,dex),handles.data.SEGROG_sigma(1,dex),'ok');
        hold on;
        plot(x_loop(dex),handles.data.SEGROG_equil(:,dex)','.k');
        plot(x_loop(dex),handles.data.SEGROG_equil(1,dex),'+b');
        plot(x_loop(dex),handles.data.SEGROG_equil(end,dex),'+g');
        hold off;
        xlabel('Rogowski Index');
        ylabel('Signal');
        title('Rogowski Reconstruction');
        ylim([min(ylim)-0.2*diff(ylim) max(ylim)+0.2*diff(ylim)]);
        xlim([min(xlim)-1 max(xlim)+1]);
    case{'PRESS_evolution'}
        %errorbar(handles.data.PRESS_S(1,:),handles.data.PRESS_target(1,:),handles.data.PRESS_sigma(1,:),'ok');
        plot(handles.data.PRESS_S(:,:)',handles.data.PRESS_equil(:,:)','k');
        hold on;
        plot(handles.data.PRESS_S(1,:),handles.data.PRESS_equil(1,:),cinitial,'LineWidth',2.0);
        plot(handles.data.PRESS_S(end,:),handles.data.PRESS_equil(end,:),cfinal,'LineWidth',2.0);
        hold off;
        xlim([0 1]);
        xlabel('Normalized Flux');
        ylabel('Pressure');
        title('Pressure Reconstruction');
    case{'NE_evolution'}
        errorbar(handles.data.NE_S(1,:),handles.data.NE_target(1,:),handles.data.NE_sigma(1,:),'ok');
        hold on;
        plot(handles.data.NE_S(:,:)',handles.data.NE_equil(:,:)','.k');
        plot(handles.data.NE_S(1,:),handles.data.NE_equil(1,:),'+b','LineWidth',2.0);
        plot(handles.data.NE_S(end,:),handles.data.NE_equil(end,:),'+g','LineWidth',2.0);
        hold off;
        xlim([0 1]);
        xlabel('Normalized Flux');
        ylabel('Electron Density');
        title('n_e Reconstruction');
    case{'TE_evolution'}
        errorbar(handles.data.TE_S(1,:),handles.data.TE_target(1,:),handles.data.TE_sigma(1,:),'ok');
        hold on;
        plot(handles.data.TE_S(:,:)',handles.data.TE_equil(:,:)','.k');
        plot(handles.data.TE_S(1,:),handles.data.TE_equil(1,:),'+b','LineWidth',2.0);
        plot(handles.data.TE_S(end,:),handles.data.TE_equil(end,:),'+g','LineWidth',2.0);
        hold off;
        xlim([0 1]);
        xlabel('Normalized Flux');
        ylabel('Electron Temperature');
        title('T_e Reconstruction');
    case{'ECE_evolution'}
        nchan = size(handles.data.ECEREFLECT_target,2);
        errorbar(1:nchan,handles.data.ECEREFLECT_target(1,:),handles.data.ECEREFLECT_sigma(1,:),'ok');
        hold on;
        plot(1:nchan,handles.data.ECEREFLECT_equil(:,:)','xk');
        plot(1:nchan,handles.data.ECEREFLECT_equil(1,:),['x' cinitial],'LineWidth',2.0);
        plot(1:nchan,handles.data.ECEREFLECT_equil(end,:),['x' cfinal],'LineWidth',2.0);
        hold off;
        xlabel('Channel');
        ylabel('Radiative Temperature');
        title('ECE Reflectrometry');
        xlim([1 nchan+1]);
    case{'ECE_xmode'}
        nchan = size(handles.data.ECEREFLECT_target,2);
        plot(1:nchan,handles.data.ECEREFLECT_tradx(:,:)','xk');
        hold on;
        plot(1:nchan,handles.data.ECEREFLECT_tradx(1,:),['x' cinitial],'LineWidth',2.0);
        plot(1:nchan,handles.data.ECEREFLECT_tradx(end,:),['x' cfinal],'LineWidth',2.0);
        hold off;
        xlabel('Channel');
        ylabel('Radiative Temperature');
        title('ECE (X-mode) Reflectrometry');
        xlim([1 nchan+1]);
    case{'ECE_omode'}
        nchan = size(handles.data.ECEREFLECT_target,2);
        plot(1:nchan,handles.data.ECEREFLECT_trado(:,:)','xk');
        hold on;
        plot(1:nchan,handles.data.ECEREFLECT_trado(1,:),['x' cinitial],'LineWidth',2.0);
        plot(1:nchan,handles.data.ECEREFLECT_trado(end,:),['x' cfinal],'LineWidth',2.0);
        hold off;
        xlabel('Channel');
        ylabel('Radiative Temperature');
        title('ECE (O-mode) Reflectrometry');
        xlim([1 nchan+1]);
    case{'TI_evolution'}
        errorbar(handles.data.TI_S(1,:),handles.data.TI_target(1,:),handles.data.TI_sigma(1,:),'ok');
        hold on;
        plot(handles.data.TI_S(:,:)',handles.data.TI_equil(:,:)','.k');
        plot(handles.data.TI_S(1,:),handles.data.TI_equil(1,:),'+b','LineWidth',2.0);
        plot(handles.data.TI_S(end,:),handles.data.TI_equil(end,:),'+g','LineWidth',2.0);
        hold off;
        xlim([0 1]);
        xlabel('Normalized Flux');
        ylabel('Ion Temperature');
        title('T_i Reconstruction');
    case{'IOTA_evolution'}
        errorbar(handles.data.IOTA_S(1,:),handles.data.IOTA_target(1,:),handles.data.IOTA_sigma(1,:),'ok');
        hold on;
        plot(handles.data.IOTA_S(:,:)',handles.data.IOTA_equil(:,:)','.k');
        plot(handles.data.IOTA_S(1,:),handles.data.IOTA_equil(1,:),'+b','LineWidth',2.0);
        plot(handles.data.IOTA_S(end,:),handles.data.IOTA_equil(end,:),'+g','LineWidth',2.0);
        hold off;
        xlim([0 1]);
        xlabel('Normalized Flux');
        ylabel('Iota');
        title('Rotational Transform');
    case{'NELINE_evolution'}
        nchan = size(handles.data.NELINE_target,2);
        errorbar(1:nchan,handles.data.NELINE_target(1,:),handles.data.NELINE_sigma(1,:),'ok');
        hold on;
        plot(1:nchan,handles.data.NELINE_equil(:,:)','xk');
        plot(1:nchan,handles.data.NELINE_equil(1,:),['x' cinitial],'LineWidth',2.0);
        plot(1:nchan,handles.data.NELINE_equil(end,:),['+' cfinal],'LineWidth',2.0);
        xlim([0 nchan+1]);
        hold off;
        xlabel('Channel');
        ylabel('Electron Density');
        title('Line Int. n_e Reconstruction');
    case{'TELINE_evolution'}
        nchan = size(handles.data.TELINE_target,2);
        errorbar(1:nchan,handles.data.TELINE_target(1,:),handles.data.TELINE_sigma(1,:),'ok');
        hold on;
        plot(1:nchan,handles.data.TELINE_equil(:,:)','xk');
        plot(1:nchan,handles.data.TELINE_equil(1,:),['x' cinitial],'LineWidth',2.0);
        plot(1:nchan,handles.data.TELINE_equil(end,:),['+' cfinal],'LineWidth',2.0);
        xlim([0 nchan+1]);
        hold off;
        xlabel('Channel');
        ylabel('Electron Temperature');
        title('Line Int. T_e Reconstruction');
    case{'TILINE_evolution'}
        nchan = size(handles.data.TILINE_target,2);
        errorbar(1:nchan,handles.data.TILINE_target(1,:),handles.data.TILINE_sigma(1,:),'ok');
        hold on;
        plot(1:nchan,handles.data.TILINE_equil(:,:)','xk');
        plot(1:nchan,handles.data.TILINE_equil(1,:),['x' cinitial],'LineWidth',2.0);
        plot(1:nchan,handles.data.TILINE_equil(end,:),['+' cfinal],'LineWidth',2.0);
        xlim([0 nchan+1]);
        hold off;
        xlabel('Channel');
        ylabel('Ion Temperature');
        title('Line Int. T_i Reconstruction');
    case{'XICS_evolution'}
        nchan = size(handles.data.XICS_target,2);
        errorbar(1:nchan,handles.data.XICS_target(1,:),handles.data.XICS_sigma(1,:),'ok');
        hold on;
        plot(1:nchan,handles.data.XICS_equil(:,:)','xk');
        plot(1:nchan,handles.data.XICS_equil(1,:),['x' cinitial],'LineWidth',2.0);
        plot(1:nchan,handles.data.XICS_equil(end,:),['+' cfinal],'LineWidth',2.0);
        xlim([0 nchan+1]);
        hold off;
        xlabel('Channel');
        ylabel('XICS Signal');
        title('Line Int. XICS T_i*Emis Reconstruction');
    case{'XICS_BRIGHT_evolution'}
        nchan = size(handles.data.XICS_BRIGHT_target,2);
        errorbar(1:nchan,handles.data.XICS_BRIGHT_target(1,:),handles.data.XICS_BRIGHT_sigma(1,:),'ok');
        hold on;
        plot(1:nchan,handles.data.XICS_BRIGHT_equil(:,:)','xk');
        plot(1:nchan,handles.data.XICS_BRIGHT_equil(1,:),['x' cinitial],'LineWidth',2.0);
        plot(1:nchan,handles.data.XICS_BRIGHT_equil(end,:),['+' cfinal],'LineWidth',2.0);
        xlim([0 nchan+1]);
        hold off;
        xlabel('Channel');
        ylabel('XICS Brightness');
        title('Line Int. XICS Emis Reconstruction');
    case{'SXR_evolution'}
        nchan = size(handles.data.SXR_target,2);
        errorbar(1:nchan,handles.data.SXR_target(1,:),handles.data.SXR_sigma(1,:),'ok');
        hold on;
        plot(1:nchan,handles.data.SXR_equil(:,:)','xk');
        plot(1:nchan,handles.data.SXR_equil(1,:),['x' cinitial],'LineWidth',2.0);
        plot(1:nchan,handles.data.SXR_equil(end,:),['+' cfinal],'LineWidth',2.0);
        xlim([0 nchan+1]);
        hold off;
        xlabel('Channel');
        ylabel('Soft X-Ray');
        title('SXR Reconstruction');
    case{'MSE_evolution'}
        errorbar(handles.data.MSE_R(1,:),handles.data.MSE_target(1,:),handles.data.MSE_sigma(1,:),'ok');
        hold on;
        plot(handles.data.MSE_R(:,:)',handles.data.MSE_equil(:,:)','.k');
        plot(handles.data.MSE_R(1,:),handles.data.MSE_equil(1,:),'+b','LineWidth',2.0);
        plot(handles.data.MSE_R(end,:),handles.data.MSE_equil(end,:),'+g','LineWidth',2.0);
        hold off;
        xlabel('R [m]');
        ylabel('Pitch Angle');
        title('MSE Reconstruction');
    case{'FLUX0'}
        files = dir('wout*.*.nc');
        theta=0:2*pi/180:2*pi;
        zeta0 =handles.phi0;
        cla;
        xlim([0 1]);
        ylim([0 1]);
        hold on;
        for i=2:length(files)-1
            vmec_data=read_vmec(files(i).name);
            r = cfunct(theta,zeta0,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
            z = sfunct(theta,zeta0,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
            plot(r(end,:,1),z(end,:,1),'k');
            plot(r(2,:,1),z(2,:,1),'k');
            plot(r(1,1,1),z(1,1,1),'k+');
            title(['Loading ' num2str(i)]);
            axis equal
            pause(0.1);
        end
        title('Loading Last two');
        vmec_data=read_vmec(files(1).name);
        r = cfunct(theta,zeta0,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
        z = sfunct(theta,zeta0,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
        plot(r(end,:,1),z(end,:,1),'b','LineWidth',2.0);
        plot(r(2,:,1),z(2,:,1),'b','LineWidth',2.0);
        plot(r(1,1,1),z(1,1,1),'b+','LineWidth',2.0);
        vmec_data=read_vmec(files(end).name);
        r = cfunct(theta,zeta0,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
        z = sfunct(theta,zeta0,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
        plot(r(end,:,1),z(end,:,1),'g','LineWidth',2.0);
        plot(r(2,:,1),z(2,:,1),'g','LineWidth',2.0);
        plot(r(1,1,1),z(1,1,1),'g+','LineWidth',2.0);
        hold off;
        axis tight;
        title('Flux surfaces');
        ylabel('Z [m]');
        xlabel('R [m]');
        axis equal;
    case{'FLUXPI'}
        files = dir('wout*.*.nc');
        theta=0:2*pi/180:2*pi;
        %zeta =handles.phi0;
        cla;
        xlim([0 1]);
        ylim([0 1]);
        hold on;
        for i=2:length(files)-1
            vmec_data=read_vmec(files(i).name);
            zeta0 = pi/vmec_data.nfp;
            r = cfunct(theta,zeta0,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
            z = sfunct(theta,zeta0,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
            plot(r(end,:,1),z(end,:,1),'k');
            plot(r(2,:,1),z(2,:,1),'k');
            plot(r(1,1,1),z(1,1,1),'k+');
            title(['Loading ' num2str(i)]);
            axis equal
            pause(0.1);
        end
        title('Loading Last two');
        vmec_data=read_vmec(files(1).name);
        zeta0 = pi/vmec_data.nfp;
        r = cfunct(theta,zeta0,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
        z = sfunct(theta,zeta0,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
        plot(r(end,:,1),z(end,:,1),'b','LineWidth',2.0);
        plot(r(2,:,1),z(2,:,1),'b','LineWidth',2.0);
        plot(r(1,1,1),z(1,1,1),'b+','LineWidth',2.0);
        vmec_data=read_vmec(files(end).name);
        zeta0 = pi/vmec_data.nfp;
        r = cfunct(theta,zeta0,vmec_data.rmnc,vmec_data.xm,vmec_data.xn);
        z = sfunct(theta,zeta0,vmec_data.zmns,vmec_data.xm,vmec_data.xn);
        plot(r(end,:,1),z(end,:,1),'g','LineWidth',2.0);
        plot(r(2,:,1),z(2,:,1),'g','LineWidth',2.0);
        plot(r(1,1,1),z(1,1,1),'g+','LineWidth',2.0);
        hold off;
        axis tight;
        title('Flux surfaces');
        ylabel('Z [m]');
        xlabel('R [m]');
        axis equal;
    case{'CURRENT'}
        files = dir('wout*.*.nc');
        cla;
        xlim([0 1]);
        ylim([0 1]);
        vmec_data=read_vmec(files(1).name);
        plot(vmec_data.phi./vmec_data.phi(end),vmec_data.jcurv,cinitial,'LineWidth',2.0);
        hold on;
        for i=2:length(files)-1
            vmec_data=read_vmec(files(i).name);
            plot(vmec_data.phi./vmec_data.phi(end),vmec_data.jcurv,'k');
        end
        vmec_data=read_vmec(files(end).name);
        plot(vmec_data.phi./vmec_data.phi(end),vmec_data.jcurv,cfinal,'LineWidth',2.0);
        %xlim([0 1]);
        axis tight;
        hold off;
        title('Current Profile');
        xlabel('Normalized Flux');
        ylabel('Current Density A/m^2');
    case{'<J*B>'}
        files = dir('wout*.*.nc');
        cla;
        xlim([0 1]);
        ylim([0 1]);
        hold on;
        for i=2:length(files)-1
            vmec_data=read_vmec(files(i).name);
            plot(vmec_data.phi./vmec_data.phi(end),vmec_data.jdotb,'k');
        end
        vmec_data=read_vmec(files(1).name);
        plot(vmec_data.phi./vmec_data.phi(end),vmec_data.jdotb,cinitial,'LineWidth',2.0);
        vmec_data=read_vmec(files(end).name);
        plot(vmec_data.phi./vmec_data.phi(end),vmec_data.jdotb,cfinal,'LineWidth',2.0);
        %xlim([0 1]);
        axis tight;
        hold off;
        title('Flux surface Averaged J\cdot B');
        xlabel('Normalized Flux');
        ylabel('<J*B>');
    case{'IOTA'}
        files = dir('wout*.*.nc');
        cla;
        xlim([0 1]);
        ylim([0 1]);
        vmec_data=read_vmec(files(1).name);
        plot(vmec_data.phi./vmec_data.phi(end),vmec_data.iotaf,cinitial,'LineWidth',2.0);
        hold on;
        for i=2:length(files)-1
            vmec_data=read_vmec(files(i).name);
            plot(vmec_data.phi./vmec_data.phi(end),vmec_data.iotaf,'k');
        end
        vmec_data=read_vmec(files(end).name);
        plot(vmec_data.phi./vmec_data.phi(end),vmec_data.iotaf,cfinal,'LineWidth',2.0);
        %xlim([0 1]);
        axis tight;
        hold off;
        title('Rotational Transform Profile');
        xlabel('Normalized Flux');
        ylabel('\iota');
    case{'PRESSURE'}
        files = dir('wout*.*.nc');
        cla;
        xlim([0 1]);
        ylim([0 1]);
        vmec_data=read_vmec(files(1).name);
        plot(vmec_data.phi./vmec_data.phi(end),vmec_data.presf,cinitial,'LineWidth',2.0);
        hold on;
        for i=2:length(files)-1
            vmec_data=read_vmec(files(i).name);
            plot(vmec_data.phi./vmec_data.phi(end),vmec_data.presf,'k');
        end
        vmec_data=read_vmec(files(end).name);
        plot(vmec_data.phi./vmec_data.phi(end),vmec_data.presf,cfinal,'LineWidth',2.0);
        %xlim([0 1]);
        axis tight;
        hold off;
        title('Pressure Profile');
        xlabel('Normalized Flux');
        ylabel('Pressure [Pa]');
    case{'EXTCUR'}
        files = dir('wout*.*.nc');
        cla;
        extcur=[];
        filename={};
        for i=1:length(files)
            vmec_data=read_vmec(files(i).name);
            extcur=[extcur; vmec_data.extcur];
            filename=[filename; files(i).name];
        end
        bar3(extcur');
        ylabel('Current Group');
        set(gca,'XTick',1:length(files));
        set(gca,'XTickLabel',filename);
        zlabel('Vaccum Field Currents');
        view(3);
    case{'CURTOR'}
        files = dir('wout*.*.nc');
        cla;
        xlim([0 1]);
        ylim([0 1]);
        for i=1:length(files)
            vmec_data=read_vmec(files(i).name);
            f(i) = vmec_data.Itor;
        end
        plot(handles.data.iter,f,'+k','MarkerSize',18,'LineWidth',4.0);
        if isfield(handles.data,'CURTOR_target')
            hold on;
            plot(xlim,handles.data.CURTOR_target(1).*[1 1],'r');
        end
        hold off;
        axis tight;
        title('Total Toroidal Current');
        xlabel('ITERATION');
        ylabel('I [kA]');
    case{'PHIEDGE'}
        files = dir('wout*.*.nc');
        cla;
        xlim([0 1]);
        ylim([0 1]);
        for i=1:length(files)
            vmec_data=read_vmec(files(i).name);
            f(i) = vmec_data.phi(end);
        end
        plot(handles.data.iter,f,'+k','MarkerSize',18,'LineWidth',4.0);
        if isfield(handles.data,'PHIEDGE_target')
            hold on;
            plot(xlim,handles.data.PHIEDGE_target(1).*[1 1],'r');
        end
        hold off;
        axis tight;
        title('Enclosed Toroidal Flux');
        xlabel('ITERATION');
        ylabel('Flux [Wb]');
    case{'NE'}
        files = dir('tprof.*.*');
        cla;
        xlim([0 1]);
        ylim([0 1]);
        for i=1:length(files)
            prof_data=importdata(files(i).name);
            f(:,i) = prof_data.data(:,2)./1E19;
        end
        s = prof_data.data(:,1);
        plot(s,f,'k');
        hold on;
        plot(s,f(:,1),cinitial,'LineWidth',2.0);
        plot(s,f(:,end),cfinal,'LineWidth',2.0);
        hold off;
        axis tight;
        title('Electron Density Profile');
        xlabel('Normalized Toroidal Flux');
        ylabel('N_e x10^{19} [m^3]');
    case{'TE'}
        files = dir('tprof.*.*');
        cla;
        xlim([0 1]);
        ylim([0 1]);
        for i=1:length(files)
            prof_data=importdata(files(i).name);
            f(:,i) = prof_data.data(:,3)./1E3;
        end
        s = prof_data.data(:,1);
        plot(s,f,'k');
        hold on;
        plot(s,f(:,1),cinitial,'LineWidth',2.0);
        plot(s,f(:,end),cfinal,'LineWidth',2.0);
        hold off;
        axis tight;
        title('Electron Temperature Profile');
        xlabel('Normalized Toroidal Flux');
        ylabel('T_e [keV]');
    case{'TI'}
        files = dir('tprof.*.*');
        cla;
        xlim([0 1]);
        ylim([0 1]);
        for i=1:length(files)
            prof_data=importdata(files(i).name);
            f(:,i) = prof_data.data(:,4)./1E3;
        end
        s = prof_data.data(:,1);
        plot(s,f,'k');
        hold on;
        plot(s,f(:,1),cinitial,'LineWidth',2.0);
        plot(s,f(:,end),cfinal,'LineWidth',2.0);
        hold off;
        axis tight;
        title('Ion Temperature Profile');
        xlabel('Normalized Toroidal Flux');
        ylabel('T_i [keV]');
    case{'ZEFF'}
        files = dir('tprof.*.*');
        cla;
        xlim([0 1]);
        ylim([0 1]);
        for i=1:length(files)
            prof_data=importdata(files(i).name);
            f(:,i) = prof_data.data(:,5);
        end
        s = prof_data.data(:,1);
        plot(s,f,'k');
        hold on;
        plot(s,f(:,1),cinitial,'LineWidth',2.0);
        plot(s,f(:,end),cfinal,'LineWidth',2.0);
        hold off;
        axis tight;
        title('Zeff Profile');
        xlabel('Normalized Toroidal Flux');
        ylabel('Zeff');
    case{'EMIS_XICS'}
        files = dir('dprof.*.*');
        cla;
        xlim([0 1]);
        ylim([0 1]);
        for i=1:length(files)
            prof_data=importdata(files(i).name);
            f(:,i) = prof_data.data(:,2);
        end
        s = prof_data.data(:,1);
        plot(s,f,'k');
        hold on;
        plot(s,f(:,1),cinitial,'LineWidth',2.0);
        plot(s,f(:,end),cfinal,'LineWidth',2.0);
        hold off;
        axis tight;
        title('Emissivity (XICS) Profile');
        xlabel('Normalized Toroidal Flux');
        ylabel('Emissivity');
    case{'G11'}
        files = dir('gist_genet_*.*');
        cla;
        xlim([0 1]);
        ylim([0 1]);
        for i=1:length(files)
            prof_data=importdata(files(i).name,' ',9);
            f(:,i) = prof_data.data(:,1);
        end
        ns = size(f,1);
        s  = -pi:2*pi/(ns-1):pi;
        plot(s,f,'k');
        hold on;
        plot(s,f(:,1),cinitial,'LineWidth',2.0);
        plot(s,f(:,end),cfinal,'LineWidth',2.0);
        hold off;
        axis tight;
        title('GIST G11');
        xlabel('\alpha');
        ylabel('G11');
    case{'G12'}
        files = dir('gist_genet_*.*');
        cla;
        xlim([0 1]);
        ylim([0 1]);
        for i=1:length(files)
            prof_data=importdata(files(i).name,' ',9);
            f(:,i) = prof_data.data(:,2);
        end
        ns = size(f,1);
        s  = -pi:2*pi/(ns-1):pi;
        plot(s,f,'k');
        hold on;
        plot(s,f(:,1),cinitial,'LineWidth',2.0);
        plot(s,f(:,end),cfinal,'LineWidth',2.0);
        hold off;
        axis tight;
        title('GIST G12');
        xlabel('\alpha');
        ylabel('G12');
    case{'G22'}
        files = dir('gist_genet_*.*');
        cla;
        xlim([0 1]);
        ylim([0 1]);
        for i=1:length(files)
            prof_data=importdata(files(i).name,' ',9);
            f(:,i) = prof_data.data(:,3);
        end
        ns = size(f,1);
        s  = -pi:2*pi/(ns-1):pi;
        plot(s,f,'k');
        hold on;
        plot(s,f(:,1),cinitial,'LineWidth',2.0);
        plot(s,f(:,end),cfinal,'LineWidth',2.0);
        hold off;
        axis tight;
        title('GIST G22');
        xlabel('\alpha');
        ylabel('G22');
    case{'Bhat'}
        files = dir('gist_genet_*.*');
        cla;
        xlim([0 1]);
        ylim([0 1]);
        for i=1:length(files)
            prof_data=importdata(files(i).name,' ',9);
            f(:,i) = prof_data.data(:,4);
        end
        ns = size(f,1);
        s  = -pi:2*pi/(ns-1):pi;
        plot(s,f,'k');
        hold on;
        plot(s,f(:,1),cinitial,'LineWidth',2.0);
        plot(s,f(:,end),cfinal,'LineWidth',2.0);
        hold off;
        axis tight;
        title('GIST Bhat');
        xlabel('\alpha');
        ylabel('Bhat');
    case{'|JAC|'}
        files = dir('gist_genet_*.*');
        cla;
        xlim([0 1]);
        ylim([0 1]);
        for i=1:length(files)
            prof_data=importdata(files(i).name,' ',9);
            f(:,i) = prof_data.data(:,5);
        end
        ns = size(f,1);
        s  = -pi:2*pi/(ns-1):pi;
        plot(s,f,'k');
        hold on;
        plot(s,f(:,1),cinitial,'LineWidth',2.0);
        plot(s,f(:,end),cfinal,'LineWidth',2.0);
        hold off;
        axis tight;
        title('GIST |JAC|');
        xlabel('\alpha');
        ylabel('JAC');
    case{'L1'}
        files = dir('gist_genet_*.*');
        cla;
        xlim([0 1]);
        ylim([0 1]);
        for i=1:length(files)
            prof_data=importdata(files(i).name,' ',9);
            f(:,i) = prof_data.data(:,7);
        end
        ns = size(f,1);
        s  = -pi:2*pi/(ns-1):pi;
        plot(s,f,'k');
        hold on;
        plot(s,f(:,1),cinitial,'LineWidth',2.0);
        plot(s,f(:,end),cfinal,'LineWidth',2.0);
        hold off;
        axis tight;
        title('GIST L1');
        xlabel('\alpha');
        ylabel('L1');
    case{'L2'}
        files = dir('gist_genet_*.*');
        cla;
        xlim([0 1]);
        ylim([0 1]);
        for i=1:length(files)
            prof_data=importdata(files(i).name,' ',9);
            f(:,i) = prof_data.data(:,6);
        end
        ns = size(f,1);
        s  = -pi:2*pi/(ns-1):pi;
        plot(s,f,'k');
        hold on;
        plot(s,f(:,1),cinitial,'LineWidth',2.0);
        plot(s,f(:,end),cfinal,'LineWidth',2.0);
        hold off;
        axis tight;
        title('GIST L2');
        xlabel('\alpha');
        ylabel('L2');
    case{'DBDT'}
        files = dir('gist_genet_*.*');
        cla;
        xlim([0 1]);
        ylim([0 1]);
        for i=1:length(files)
            prof_data=importdata(files(i).name,' ',9);
            f(:,i) = prof_data.data(:,8);
        end
        ns = size(f,1);
        s  = -pi:2*pi/(ns-1):pi;
        plot(s,f,'k');
        hold on;
        plot(s,f(:,1),cinitial,'LineWidth',2.0);
        plot(s,f(:,end),cfinal,'LineWidth',2.0);
        hold off;
        axis tight;
        title('GIST DBDT');
        xlabel('\alpha');
        ylabel('DBDT');
    case{'KP1'}
        files = dir('curv_*.*');
        cla;
        xlim([0 1]);
        ylim([0 1]);
        for i=1:length(files)
            prof_data=importdata(files(i).name);
            f(:,i) = prof_data(:,1);
        end
        ns = size(f,1);
        s  = -pi:2*pi/(ns-1):pi;
        plot(s,f,'k');
        hold on;
        plot(s,f(:,1),cinitial,'LineWidth',2.0);
        plot(s,f(:,end),cfinal,'LineWidth',2.0);
        hold off;
        axis tight;
        title('GIST KP1');
        xlabel('\alpha');
        ylabel('KP1');
    case{'PROXY'}
        files = dir('txport_out.*.*');
        cla;
        xlim([0 1]);
        ylim([0 1]);
        for i=1:length(files)
            prof_data=importdata(files(i).name);
            f(:,i) = prof_data(2:end,1);
        end
        ns = size(f,1);
        s  = -pi:2*pi/(ns-1):pi;
        plot(s,f,'k');
        hold on;
        plot(s,f(:,1),cinitial,'LineWidth',2.0);
        plot(s,f(:,end),cfinal,'LineWidth',2.0);
        hold off;
        axis tight;
        title([ 'TXPORT PROXY ' handles.data.TXPORT_prox]);
        xlabel('\alpha');
        ylabel('PROXY');
    case{'XI_00','XI_01'}
        files = dir(['terpsichore_16.*' stemp(3:end)]);
        cla;
        xlim([0 1]);
        ylim([0 1]);
        data=read_terpsichore(files(1).name);
        ximax = data.ximax;
        for j= 1:5
            [val idex] = max(ximax);
            hold on; plot(data.xi(idex,:),'Color',cinitial);
            ximax(idex)=-1E30;
        end
        data=[]; ximax=[];
        data=read_terpsichore(files(end).name);
        ximax = data.ximax;
        for j= 1:5
            [val idex] = max(ximax);
            hold on; plot(data.xi(idex,:),'Color',cfinal);
            ximax(idex)=-1E30;
        end
        axis tight;
        title(['\xi (' stemp(4:end) ')']);
        xlabel('radial index');
        ylabel('\xi');
    case{'ETA_00','ETA_01'}
        files = dir(['terpsichore_16.*' stemp(4:end)]);
        cla;
        xlim([0 1]);
        ylim([0 1]);
        data=read_terpsichore(files(1).name);
        etamax = data.etamax;
        for j= 1:5
            [val idex] = max(etamax);
            hold on; plot(data.eta(idex,:),'Color',cinitial);
            etamax(idex)=-1E30;
        end
        data=[]; etamax=[];
        data=read_terpsichore(files(end).name);
        etamax = data.etamax;
        for j= 1:5
            [val idex] = max(etamax);
            hold on; plot(data.eta(idex,:),'Color',cfinal);
            etamax(idex)=-1E30;
        end
        axis tight;
        title(['\eta (' stemp(5:end) ')']);
        xlabel('radial index');
        ylabel('\eta');
end
hold off;
% Add labels
if isempty(strfind(stemp,'Chi-')) && isempty(strfind(stemp,'_chisq'))
    x1=xlim;
    y1=ylim;
    text(min(x1)+0.05*diff(x1),max(y1)-0.05*diff(y1),'Initial','Color',cinitial,'FontSize',18)
    text(min(x1)+0.05*diff(x1),max(y1)-0.10*diff(y1),'Intermediate','Color','k','FontSize',18)
    text(min(x1)+0.05*diff(x1),max(y1)-0.15*diff(y1),'Final','Color',cfinal,'FontSize',18)
end


% --- Executes during object creation, after setting all properties.
function pulldownmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pulldownmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function save_name_Callback(hObject, eventdata, handles)
% hObject    handle to save_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of save_name as text
%        str2double(get(hObject,'String')) returns contents of save_name as a double


% --- Executes during object creation, after setting all properties.
function save_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to save_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in save_button.
function save_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filename=get(handles.save_name,'String');
c = strfind(filename,'.');
c = max(c);
ddex=c+1;
%[S c]=textscan(filename,'%[^.]');
%ddex=c+2;
nlen=size(filename,2);
lvrml=0;
if ddex < nlen
    filename(ddex:nlen)=lower(filename(ddex:nlen));
    switch filename(ddex:nlen)
        case {'ai','bmp','emf','eps','fig','jpg','m','pbm','pcx','pdf',...
                'pgm','png','ppm','tif'}
        case {'wrl'}
            lvrml=1;
        otherwise
            filename=strcat(filename(1:c),'.fig');
            %set(handles.statustext,'String','Saving as matlab fig!','ForegroundColor','red');
            %pause(1.0);
    end
else
    filename=strcat(filename,'.fig');
end
h=figure('Visible','on');
update_plots(handles);
pause(.1);
if lvrml
    vrml(h,filename);
else
    saveas(h,filename);
end
%set(handles.statustext,'String','Ready','ForegroundColor','black');
%pause(.1);
close(h);
