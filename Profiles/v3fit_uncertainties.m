function varargout = v3fit_uncertainties(varargin)
% V3FIT_UNCERTAINTIES MATLAB code for v3fit_uncertainties.fig
%      V3FIT_UNCERTAINTIES, by itself, creates a new V3FIT_UNCERTAINTIES or raises the existing
%      singleton*.
%
%      H = V3FIT_UNCERTAINTIES returns the handle to a new V3FIT_UNCERTAINTIES or the handle to
%      the existing singleton*.
%
%      V3FIT_UNCERTAINTIES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in V3FIT_UNCERTAINTIES.M with the given input arguments.
%
%      V3FIT_UNCERTAINTIES('Property','Value',...) creates a new V3FIT_UNCERTAINTIES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before v3fit_uncertainties_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to v3fit_uncertainties_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help v3fit_uncertainties

% Last Modified by GUIDE v2.5 09-Sep-2020 19:46:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @v3fit_uncertainties_OpeningFcn, ...
    'gui_OutputFcn',  @v3fit_uncertainties_OutputFcn, ...
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


% --- Executes just before v3fit_uncertainties is made visible.
function v3fit_uncertainties_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to v3fit_uncertainties (see VARARGIN)

% Choose default command line output for v3fit_uncertainties
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes v3fit_uncertainties wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = v3fit_uncertainties_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3


% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton4


% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton5


% --- Executes on button press in radiobutton6.
function radiobutton6_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton6


% --- Executes on button press in radiobutton7.
function radiobutton7_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton7

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the parameter values for the pressure and current profiles
pres_scale = str2num(get(handles.edit1, 'String'));
AM0 = str2num(get(handles.edit2, 'String'));
AM1 = str2num(get(handles.edit3, 'String'));
AM2 = str2num(get(handles.edit4, 'String'));
AM3 = str2num(get(handles.edit5, 'String'));
AM4 = str2num(get(handles.edit6, 'String'));
AM5 = str2num(get(handles.edit7, 'String'));

curtor = str2num(get(handles.edit10, 'String'));
AC0 = str2num(get(handles.edit11, 'String'));
AC1 = str2num(get(handles.edit12, 'String'));
AC2 = str2num(get(handles.edit13, 'String'));
AC3 = str2num(get(handles.edit14, 'String'));
AC4 = str2num(get(handles.edit15, 'String'));
AC5 = str2num(get(handles.edit16, 'String'));

% Check the radio buttons and find the 'active' variables
num_active_P = 0;
ind_active_P = [];
num_active_I = 0;
ind_active_I = [];

rb_ps = get(handles.radiobutton1, 'Value');
if (rb_ps) num_active_P = num_active_P + 1; ind_active_P = [1]; end
rb_AM0 = get(handles.radiobutton2, 'Value');
if (rb_AM0) num_active_P = num_active_P + 1; ind_active_P = [ind_active_P 2];end
rb_AM1 = get(handles.radiobutton3, 'Value');
if (rb_AM1) num_active_P = num_active_P + 1; ind_active_P = [ind_active_P 3]; end
rb_AM2 = get(handles.radiobutton4, 'Value');
if (rb_AM2) num_active_P = num_active_P + 1; ind_active_P = [ind_active_P 4]; end
rb_AM3 = get(handles.radiobutton5, 'Value');
if (rb_AM3) num_active_P = num_active_P + 1; ind_active_P = [ind_active_P 5]; end
rb_AM4 = get(handles.radiobutton6, 'Value');
if (rb_AM4) num_active_P = num_active_P + 1; ind_active_P = [ind_active_P 6]; end
rb_AM5 = get(handles.radiobutton7, 'Value');
if (rb_AM5) num_active_P = num_active_P + 1; ind_active_P = [ind_active_P 7]; end

rb_ct = get(handles.radiobutton8, 'Value');
if (rb_ct) num_active_I = num_active_I + 1; ind_active_I = [1]; end
rb_AC0 = get(handles.radiobutton9, 'Value');
if (rb_AC0) num_active_I = num_active_I + 1; ind_active_I = [ind_active_I 2];end
rb_AC1 = get(handles.radiobutton10, 'Value');
if (rb_AC1) num_active_I = num_active_I + 1; ind_active_I = [ind_active_I 3]; end
rb_AC2 = get(handles.radiobutton11, 'Value');
if (rb_AC2) num_active_I = num_active_I + 1; ind_active_I = [ind_active_I 4]; end
rb_AC3 = get(handles.radiobutton12, 'Value');
if (rb_AC3) num_active_I = num_active_I + 1; ind_active_I = [ind_active_I 5]; end
rb_AC4 = get(handles.radiobutton13, 'Value');
if (rb_AC4) num_active_I = num_active_I + 1; ind_active_I = [ind_active_I 6]; end
rb_AC5 = get(handles.radiobutton14, 'Value');
if (rb_AC5) num_active_I = num_active_I + 1; ind_active_I = [ind_active_I 7]; end

C_P_all = get(handles.uitable1, 'Data')
C_I_all = get(handles.uitable2, 'Data')

C_P2 = zeros(num_active_P, num_active_P);
C_P2 = C_P_all([ind_active_P], [ind_active_P])

C_I2 = zeros(num_active_I, num_active_I);
C_I2 = C_I_all([ind_active_I], [ind_active_I])


sdata = linspace(0,1,51);
rhodata = sqrt(sdata);

pressure_profile_params = [AM0 AM1 AM2 AM3 AM4 AM5];
pressure_solution = pres_scale * two_two_power(sdata, pressure_profile_params);


ttp_params = [pres_scale AM0 AM1 AM2 AM3 AM4 AM5];
[ dP_dparam ] = two_two_power_fd(sdata, ttp_params, ind_active_P  );

current_profile_params = [AC0 AC1 AC2 AC3 AC4 AC5];
if get(handles.radiobutton15, 'Value')
    J_Norm = curtor * (AC0 - 1) / sum(current_profile_params(2:end));
    current_solution = J_Norm * sum_cossq_s_j(sdata, current_profile_params);
    scsq_params = [curtor AC0 AC1 AC2 AC3 AC4 AC5];
    [ dI_dparam ] = sum_cossq_s_j_fd(sdata, scsq_params, ind_active_I );
elseif get(handles.radiobutton16, 'Value')
    J_Norm = curtor / ( AC1*(pi^2-4)/(32*pi^2) + AC2/8 + AC3/4 + AC4*3/8 + ...
        AC5*(7/32 + 1/(8*pi^2)));
    current_solution = J_Norm * sum_cossq_sqrts_j(sdata, current_profile_params);
    scsq_params = [curtor AC0 AC1 AC2 AC3 AC4 AC5];
    [ dI_dparam ] = sum_cossq_sqrts_j_fd(sdata, scsq_params, ind_active_I );    
end

%
% if 0
%     C_M = zeros(length(sdata), 1);
%     % Now project
%     for ii = 1:length(sdata)
%         this_K = dP_dparam(ii,:);
%         C_M(ii) = this_K * C_p * this_K';
%     end
% else
C_M_P = dP_dparam * C_P2 * dP_dparam';
C_M_I = dI_dparam * C_I2 * dI_dparam';
% end

sigma_MP = zeros(51, 1);
for ii = 1:length(sdata)
    sigma_MP(ii) = sqrt(C_M_P(ii,ii));
end
sigma_MP

sigma_MI = zeros(51, 1);
for ii = 1:length(sdata)
    sigma_MI(ii) = sqrt(C_M_I(ii,ii));
end
sigma_MI



% Figures

% Pressure

axes(handles.axes2a);
box on; hold off
for ii = 1:num_active_P
    plot(sdata, dP_dparam(:,ii), '-');
    hold on
end
title('d Pressure / d param ')
%legend('P_{scale}', 'A_0', 'A_1', 'A_2', 'A_3', 'A_4', 'A_5' );
legend([strtrim(handles.v3data.param_name_clean),num2str( handles.v3data.param_index(1,:)')])

handles.v3data.param_sigma
%da_sum = 0*dP_dparam(:,1);


axes(handles.axes2b);
box on; hold off
for ii = 1:num_active_P
    plot(sdata, handles.v3data.param_sigma(handles.map_P(ii),end) * dP_dparam(:,ii), '-');
    %da_sum = sqrt(da_sum.^2 + (1.0/num_active_P) * (handles.v3data.param_sigma(ii,end) * dP_dparam(:,ii)).^2);
    hold on
end
title('d Pressure / d param ')
%legend('P_{scale}', 'A_0', 'A_1', 'A_2', 'A_3', 'A_4', 'A_5' );

legend([strtrim(handles.v3data.param_name_clean(handles.map_P,:)),num2str( handles.v3data.param_index(1,(handles.map_P))')])


axes(handles.axes3);
box on; hold off
%plot(sdata, sqrt(C_M), '-');
plot(sdata, sigma_MP, '-');
hold on;
%plot(sdata, da_sum, '-');
hold on
legend('\sigma_M', 'da sum' );


axes(handles.axes1);
box on; hold off
plot(sdata, pressure_solution, '-');
hold on
legend('Total Pressure');

%errorbar(sdata, pressure_solution, sqrt(C_M), 'r:');
%errorbar(sdata, pressure_solution, da_sum, 'g--', 'Linewidth', 2);
errorbar(sdata, pressure_solution, sigma_MP, 'r:', 'Linewidth', 1);
hold on
% legend('Total Pressure');

%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Current

axes(handles.axes5a);
box on; hold off
for ii = 1:num_active_I
    plot(sdata, dI_dparam(:,ii), '-');
    hold on
end
title('d J / d param ')
%legend('P_{scale}', 'A_0', 'A_1', 'A_2', 'A_3', 'A_4', 'A_5' );
legend([strtrim(handles.v3data.param_name_clean(handles.map_I,:)),num2str( handles.v3data.param_index(1,(handles.map_I))')])

handles.v3data.param_sigma
%da_sum = 0*dI_dparam(:,1);


axes(handles.axes5b);
box on; hold off
for ii = 1:num_active_I
    plot(sdata, handles.v3data.param_sigma(handles.map_I(ii),end) * dI_dparam(:,ii), '-');
    %da_sum = sqrt(da_sum.^2 + (1.0/num_active_I) * (handles.v3data.param_sigma(ii,end) * dI_dparam(:,ii)).^2);
    hold on
end
title('d J / d param ')
%legend('P_{scale}', 'A_0', 'A_1', 'A_2', 'A_3', 'A_4', 'A_5' );



axes(handles.axes6);
box on; hold off
%plot(sdata, sqrt(C_M), '-');
plot(sdata, sigma_MI, '-');
hold on;
%plot(sdata, da_sum, '-');
hold on
legend('\sigma_M', 'da sum' );


axes(handles.axes4);
box on; hold off
plot(sdata, current_solution, '-');
hold on
legend('Total J');

%errorbar(sdata, pressure_solution, sqrt(C_M), 'r:');
%errorbar(sdata, current_solution, da_sum, 'g--', 'Linewidth', 2);
errorbar(sdata, current_solution, sigma_MI, 'r:', 'Linewidth', 1);
hold on
% legend('Total Pressure');







% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

vmec_result_ext = get(handles.edit9, 'String');
handles.vmdata = load_vmec(vmec_result_ext);

v3fit_result_filename = get(handles.edit8, 'String');
handles.v3data = load_v3fit_result(v3fit_result_filename);

guidata(hObject, handles);



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
v3data = handles.v3data;

set(handles.radiobutton1, 'Value', 0.0);
set(handles.radiobutton2, 'Value', 0.0);
set(handles.radiobutton3, 'Value', 0.0);
set(handles.radiobutton4, 'Value', 0.0);
set(handles.radiobutton5, 'Value', 0.0);
set(handles.radiobutton6, 'Value', 0.0);
set(handles.radiobutton7, 'Value', 0.0);

set(handles.radiobutton8, 'Value', 0.0);
set(handles.radiobutton9, 'Value', 0.0);
set(handles.radiobutton10, 'Value', 0.0);
set(handles.radiobutton11, 'Value', 0.0);
set(handles.radiobutton12, 'Value', 0.0);
set(handles.radiobutton13, 'Value', 0.0);
set(handles.radiobutton14, 'Value', 0.0);

pind_active = [];
cind_active = [];
pmap = [];
cmap = [];

% import vmec data.  Can't handle pres_scale, so will set to '1' and hope
% it gets set later
set(handles.edit1, 'String', num2str(1));
set(handles.edit2, 'String', num2str(handles.vmdata.am(1)));
set(handles.edit3, 'String', num2str(handles.vmdata.am(2)));
set(handles.edit4, 'String', num2str(handles.vmdata.am(3)));
set(handles.edit5, 'String', num2str(handles.vmdata.am(4)));
set(handles.edit6, 'String', num2str(handles.vmdata.am(5)));
set(handles.edit7, 'String', num2str(handles.vmdata.am(6)));

set(handles.edit10, 'String', num2str(handles.vmdata.ctor));
set(handles.edit11, 'String', num2str(handles.vmdata.ac(1)));
set(handles.edit12, 'String', num2str(handles.vmdata.ac(2)));
set(handles.edit13, 'String', num2str(handles.vmdata.ac(3)));
set(handles.edit14, 'String', num2str(handles.vmdata.ac(4)));
set(handles.edit15, 'String', num2str(handles.vmdata.ac(5)));
set(handles.edit16, 'String', num2str(handles.vmdata.ac(6)));


% apply v3fit reconstruction data
% import parameter covariance matrix - update table on screen
handles.C_p = v3data.param_corr(:,:,end);

new_data_P = zeros(7,7);
new_data_I = zeros(7,7);

the_nparam = v3data.nparam;
P_map = zeros(the_nparam, the_nparam);
I_map = zeros(the_nparam, the_nparam);

for ii = 1:the_nparam
    this_param = strtrim(lower(v3data.param_name_clean(ii,:)))
    switch this_param
        case 'pres_scale'
            set(handles.edit1, 'String', num2str(v3data.param_value(ii,end)));
            set(handles.radiobutton1, 'Value', 1.0);
            pind_active = [pind_active 1];
            pmap = [pmap ii];
        case 'am'
            this_ind = v3data.param_index(1,ii)
            switch this_ind
                case 0
                    set(handles.edit2, 'String', num2str(v3data.param_value(ii,end)));
                    set(handles.radiobutton2, 'Value', 1.0);
                    pind_active = [pind_active 2];
                    pmap = [pmap ii];
                case 1
                    set(handles.edit3, 'String', num2str(v3data.param_value(ii,end)));
                    set(handles.radiobutton3, 'Value', 1.0);
                    pind_active = [pind_active 3];
                    pmap = [pmap ii];
                case 2
                    set(handles.edit4, 'String', num2str(v3data.param_value(ii,end)));
                    set(handles.radiobutton4, 'Value', 1.0);
                    pind_active = [pind_active 4];
                    pmap = [pmap ii];
                case 3
                    set(handles.edit5, 'String', num2str(v3data.param_value(ii,end)));
                    set(handles.radiobutton5, 'Value', 1.0);
                    pind_active = [pind_active 5];
                    pmap = [pmap ii];
                case 4
                    set(handles.edit6, 'String', num2str(v3data.param_value(ii,end)));
                    set(handles.radiobutton6, 'Value', 1.0);
                    pind_active = [pind_active 6];
                    pmap = [pmap ii];
                case 5
                    set(handles.edit7, 'String', num2str(v3data.param_value(ii,end)));
                    set(handles.radiobutton7, 'Value', 1.0);
                    pind_active = [pind_active 7];
                    pmap = [pmap ii];
            end
        case 'curtor'
            set(handles.edit10, 'String', num2str(v3data.param_value(ii,end)));
            set(handles.radiobutton8, 'Value', 1.0);
            cind_active = [cind_active 1];
            cmap = [cmap ii];
        case 'ac'
            this_ind = v3data.param_index(1,ii)
            switch this_ind
                case 0
                    set(handles.edit11, 'String', num2str(v3data.param_value(ii,end)));
                    set(handles.radiobutton9, 'Value', 1.0);
                    cind_active = [cind_active 2];
                    cmap = [cmap ii];
                case 1
                    set(handles.edit12, 'String', num2str(v3data.param_value(ii,end)));
                    set(handles.radiobutton10, 'Value', 1.0);
                    cind_active = [cind_active 3];
                    cmap = [cmap ii];
                case 2
                    set(handles.edit13, 'String', num2str(v3data.param_value(ii,end)));
                    set(handles.radiobutton11, 'Value', 1.0);
                    cind_active = [cind_active 4];
                    cmap = [cmap ii];
                case 3
                    set(handles.edit14, 'String', num2str(v3data.param_value(ii,end)));
                    set(handles.radiobutton12, 'Value', 1.0);
                    cind_active = [cind_active 5];
                    cmap = [cmap ii];
                case 4
                    set(handles.edit15, 'String', num2str(v3data.param_value(ii,end)));
                    set(handles.radiobutton13, 'Value', 1.0);
                    cind_active = [cind_active 6];
                    cmap = [cmap ii];
                case 5
                    set(handles.edit16, 'String', num2str(v3data.param_value(ii,end)));
                    set(handles.radiobutton14, 'Value', 1.0);
                    cind_active = [cind_active 7];
                    cmap = [cmap ii];
            end
        otherwise
            %keyboard
    end
end


num_active_P = length(pind_active);
num_active_I = length(cind_active);


for ii = 1:num_active_P
    for jj = 1:num_active_P
        new_data_P(pind_active(ii), pind_active(jj)) = handles.C_p(pmap(ii), pmap(jj));
    end
end

for ii = 1:num_active_I
    for jj = 1:num_active_I
        new_data_I(cind_active(ii), cind_active(jj)) = handles.C_p(cmap(ii), cmap(jj));
    end
end

set(handles.uitable1, 'Data', new_data_P);
set(handles.uitable2, 'Data', new_data_I);

handles.ind_active_P = pind_active;
handles.ind_active_I = cind_active;
handles.map_P = pmap
handles.map_I = cmap;

handles.v3data.param_sigma

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function pushbutton3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton8.
function radiobutton8_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton8


% --- Executes on button press in radiobutton9.
function radiobutton9_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton9


% --- Executes on button press in radiobutton10.
function radiobutton10_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton10


% --- Executes on button press in radiobutton11.
function radiobutton11_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton11


% --- Executes on button press in radiobutton12.
function radiobutton12_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton12


% --- Executes on button press in radiobutton13.
function radiobutton13_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton13


% --- Executes on button press in radiobutton14.
function radiobutton14_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton14
