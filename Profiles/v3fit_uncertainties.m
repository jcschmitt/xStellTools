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

% Last Modified by GUIDE v2.5 09-Sep-2020 17:29:46

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

pres_scale = str2num(get(handles.edit1, 'String'));
A0 = str2num(get(handles.edit2, 'String'));
A1 = str2num(get(handles.edit3, 'String'));
A2 = str2num(get(handles.edit4, 'String'));
A3 = str2num(get(handles.edit5, 'String'));
A4 = str2num(get(handles.edit6, 'String'));
A5 = str2num(get(handles.edit7, 'String'));

num_active = 0;
ind_active = [];
rb_ps = get(handles.radiobutton1, 'Value');
if (rb_ps) num_active = num_active + 1; ind_active = [1]; end
rb_A0 = get(handles.radiobutton2, 'Value');
if (rb_A0) num_active = num_active + 1; ind_active = [ind_active 2];end
rb_A1 = get(handles.radiobutton3, 'Value');
if (rb_A1) num_active = num_active + 1; ind_active = [ind_active 3]; end
rb_A2 = get(handles.radiobutton4, 'Value');
if (rb_A2) num_active = num_active + 1; ind_active = [ind_active 4]; end
rb_A3 = get(handles.radiobutton5, 'Value');
if (rb_A3) num_active = num_active + 1; ind_active = [ind_active 5]; end
rb_A4 = get(handles.radiobutton6, 'Value');
if (rb_A4) num_active = num_active + 1; ind_active = [ind_active 6]; end
rb_A5 = get(handles.radiobutton7, 'Value');
if (rb_A5) num_active = num_active + 1; ind_active = [ind_active 7]; end

C_p_all = get(handles.uitable1, 'Data')

C_p = zeros(num_active, num_active);
C_p = C_p_all([ind_active], [ind_active])



sdata = linspace(0,1,51);
rhodata = sqrt(sdata);
two_two_power_profile = [A0 A1 A2 A3 A4 A5];
pressure_solution = pres_scale * two_two_power(sdata, two_two_power_profile);

axes(handles.axes1);
box on; hold off
plot(sdata, pressure_solution, '-');
hold on
legend('Total Pressure');

ttp_params = [pres_scale A0 A1 A2 A3 A4 A5];
[ dP_dparam ] = two_two_power_fd(sdata, ttp_params, ind_active  );

axes(handles.axes2a);
box on; hold off
for ii = 1:num_active
    plot(sdata, dP_dparam(:,ii), '-');
    hold on
end
title('d Pressure / d param ')
%legend('P_{scale}', 'A_0', 'A_1', 'A_2', 'A_3', 'A_4', 'A_5' );

legend([strtrim(handles.v3data.param_name_clean),num2str( handles.v3data.param_index(1,:)')])

handles.v3data.param_sigma

da_sum = 0*dP_dparam(:,1);
axes(handles.axes2b);
box on; hold off
for ii = 1:num_active
    plot(sdata, handles.v3data.param_sigma(ii,end) * dP_dparam(:,ii), '-');
    da_sum = sqrt(da_sum.^2 + (1.0/num_active) * (handles.v3data.param_sigma(ii,end) * dP_dparam(:,ii)).^2);
    hold on
end
title('d Pressure / d param ')
%legend('P_{scale}', 'A_0', 'A_1', 'A_2', 'A_3', 'A_4', 'A_5' );

legend([strtrim(handles.v3data.param_name_clean),num2str( handles.v3data.param_index(1,:)')])


%
% if 0
%     C_M = zeros(length(sdata), 1);
%     % Now project
%     for ii = 1:length(sdata)
%         this_K = dP_dparam(ii,:);
%         C_M(ii) = this_K * C_p * this_K';
%     end
% else
C_M = dP_dparam * C_p * dP_dparam';
% end

sigma_M = zeros(51, 1);
for ii = 1:length(sdata)
    sigma_M(ii) = sqrt(C_M(ii,ii));
end
sigma_M
axes(handles.axes3);
box on; hold off
%plot(sdata, sqrt(C_M), '-');
plot(sdata, sigma_M, '-');
hold on;
plot(sdata, da_sum, '-');

hold on
legend('\sigma_M', 'da sum' );


axes(handles.axes1);
%errorbar(sdata, pressure_solution, sqrt(C_M), 'r:');
errorbar(sdata, pressure_solution, da_sum, 'g--', 'Linewidth', 2);
errorbar(sdata, pressure_solution, sigma_M, 'r:', 'Linewidth', 2);
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
ind_active = [];

% import vmec data.  Can't handle pres_scale, so will set to '1' and hope
% it gets set later
set(handles.edit1, 'String', num2str(1));
set(handles.edit2, 'String', num2str(handles.vmdata.am(1)));
set(handles.edit3, 'String', num2str(handles.vmdata.am(2)));
set(handles.edit4, 'String', num2str(handles.vmdata.am(3)));
set(handles.edit5, 'String', num2str(handles.vmdata.am(4)));
set(handles.edit6, 'String', num2str(handles.vmdata.am(5)));
set(handles.edit7, 'String', num2str(handles.vmdata.am(6)));


% apply v3fit reconstruction data
the_nparam = v3data.nparam
for ii = 1:the_nparam
    this_param = strtrim(lower(v3data.param_name_clean(ii,:)))
    switch this_param
        case 'blah'
            %keyboard
        case 'pres_scale'
            set(handles.edit1, 'String', num2str(v3data.param_value(ii,end)));
            set(handles.radiobutton1, 'Value', 1.0);
            ind_active = [ind_active 1];
            
        case 'am'
            this_ind = v3data.param_index(1,ii)
            switch this_ind
                case 0
                    set(handles.edit2, 'String', num2str(v3data.param_value(ii,end)));
                    set(handles.radiobutton2, 'Value', 1.0);
                    ind_active = [ind_active 2];
                case 1
                    set(handles.edit3, 'String', num2str(v3data.param_value(ii,end)));
                    set(handles.radiobutton3, 'Value', 1.0);
                    ind_active = [ind_active 3];
                case 2
                    set(handles.edit4, 'String', num2str(v3data.param_value(ii,end)));
                    set(handles.radiobutton4, 'Value', 1.0);
                    ind_active = [ind_active 4];
                case 3
                    set(handles.edit5, 'String', num2str(v3data.param_value(ii,end)));
                    set(handles.radiobutton5, 'Value', 1.0);
                    ind_active = [ind_active 5];
                case 4
                    set(handles.edit6, 'String', num2str(v3data.param_value(ii,end)));
                    set(handles.radiobutton6, 'Value', 1.0);
                    ind_active = [ind_active 6];
                case 5
                    set(handles.edit7, 'String', num2str(v3data.param_value(ii,end)));
                    set(handles.radiobutton7, 'Value', 1.0);
                    ind_active = [ind_active 7];
            end
        otherwise
            %keyboard
    end
end

% import parameter covariance matrix - update table on screen
handles.C_p = v3data.param_corr(:,:,end);
num_active = length(ind_active);
new_data = zeros(7,7);
for ii = 1:num_active
    for jj = 1:num_active
        new_data(ind_active(ii), ind_active(jj)) = handles.C_p(ii, jj);
    end
end

set(handles.uitable1, 'Data', new_data);

handles.ind_active = ind_active;
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
