function varargout = pressure_fit(varargin)
% PRESSURE_FIT MATLAB code for pressure_fit.fig
%      PRESSURE_FIT, by itself, creates a new PRESSURE_FIT or raises the existing
%      singleton*.
%
%      H = PRESSURE_FIT returns the handle to a new PRESSURE_FIT or the handle to
%      the existing singleton*.
%
%      PRESSURE_FIT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PRESSURE_FIT.M with the given input arguments.
%
%      PRESSURE_FIT('Property','Value',...) creates a new PRESSURE_FIT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pressure_fit_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pressure_fit_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pressure_fit

% Last Modified by GUIDE v2.5 17-Jun-2020 09:59:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pressure_fit_OpeningFcn, ...
                   'gui_OutputFcn',  @pressure_fit_OutputFcn, ...
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


% --- Executes just before pressure_fit is made visible.
function pressure_fit_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pressure_fit (see VARARGIN)

% Choose default command line output for pressure_fit
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pressure_fit wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = pressure_fit_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in loadAnalysisData.
function loadAnalysisData_Callback(hObject, eventdata, handles)
% hObject    handle to loadAnalysisData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


handles.analysisData = load_analysis_data;

% Update handles structure
guidata(hObject, handles);



function ProfileIndex_Callback(hObject, eventdata, handles)
% hObject    handle to ProfileIndex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ProfileIndex as text
%        str2double(get(hObject,'String')) returns contents of ProfileIndex as a double


% --- Executes during object creation, after setting all properties.
function ProfileIndex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ProfileIndex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listFitTypes.
function listFitTypes_Callback(hObject, eventdata, handles)
% hObject    handle to listFitTypes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listFitTypes contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listFitTypes


% --- Executes during object creation, after setting all properties.
function listFitTypes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listFitTypes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in FindFit.
function FindFit_Callback(hObject, eventdata, handles)
% hObject    handle to FindFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


index_req = str2num(get(handles.ProfileIndex, 'String'));
MinorRadius = str2num(get(handles.MinorRadius, 'String')) * 0.01;
xdata = handles.analysisData{index_req}(:,1);
ydata = handles.analysisData{index_req}(:,4);
sigdata = handles.analysisData{index_req}(:,5);


ind_in = find(xdata<=MinorRadius);

sdata_fit = (xdata(ind_in) / MinorRadius).^2;
ydata_fit = ydata(ind_in);
two_power_guess = [sdata_fit(1) 2  2  .5 10 10];
two_power_min =   [0            1  1  0  1  1 ];
two_power_max =   [1e3          30 30 1  30 30]

pressure_fit = @(am_array)two_two_power(sdata_fit, am_array) - ydata_fit;

two_power_solution = lsqnonlin(pressure_fit, two_power_guess, two_power_min, two_power_max)

pressure_solution = two_two_power(sdata_fit, two_power_solution);





axis(handles.axes1);


box on; hold off
errorbar(xdata, ydata, sigdata);
hold on
plot(xdata(ind_in), ydata(ind_in), 'x');
plot(xdata(ind_in), pressure_solution, '-');
legend('A.R.', 'In Plasma', 'Fit');

text(0., two_power_solution(1) /2, num2str(two_power_solution))


function MinorRadius_Callback(hObject, eventdata, handles)
% hObject    handle to MinorRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MinorRadius as text
%        str2double(get(hObject,'String')) returns contents of MinorRadius as a double


% --- Executes during object creation, after setting all properties.
function MinorRadius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MinorRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
