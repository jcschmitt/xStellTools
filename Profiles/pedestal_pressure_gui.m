function varargout = pedestal_pressure_gui(varargin)
% PEDESTAL_PRESSURE_GUI MATLAB code for pedestal_pressure_gui.fig
%      PEDESTAL_PRESSURE_GUI, by itself, creates a new PEDESTAL_PRESSURE_GUI or raises the existing
%      singleton*.
%
%      H = PEDESTAL_PRESSURE_GUI returns the handle to a new PEDESTAL_PRESSURE_GUI or the handle to
%      the existing singleton*.
%
%      PEDESTAL_PRESSURE_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PEDESTAL_PRESSURE_GUI.M with the given input arguments.
%
%      PEDESTAL_PRESSURE_GUI('Property','Value',...) creates a new PEDESTAL_PRESSURE_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pedestal_pressure_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pedestal_pressure_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pedestal_pressure_gui

% Last Modified by GUIDE v2.5 22-Aug-2018 20:50:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pedestal_pressure_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @pedestal_pressure_gui_OutputFcn, ...
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


% --- Executes just before pedestal_pressure_gui is made visible.
function pedestal_pressure_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pedestal_pressure_gui (see VARARGIN)

% Choose default command line output for pedestal_pressure_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pedestal_pressure_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = pedestal_pressure_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function pressure = pedestal_pressure(am, ss_in)
if am(4) <= 0
    am(1:5) = 0;
    am(4) = 1e30;
else
    am(5) = 1.0 / (tanh(2*am(3)/am(4)) - tanh(2*(am(3)-1) / am(4)));
end

pressure = 0 * ss_in;
pressure = pressure + am(5) * am(2) * (tanh(2*(am(3) - sqrt(ss_in)) / ...
    am(4)) - tanh(2*(am(3) - 1.0) / am(4) ) );


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
doit(handles);


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
doit(handles);


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
doit(handles);


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

doit(handles);

function doit(handles)

ss = linspace(0,1,51);

am16 = 0;
am17 = get(handles.slider1,'Value');
am18 = get(handles.slider2,'Value');
am19 = get(handles.slider3,'Value');
am = [am16 am17 am18 am19];
pressure = pedestal_pressure(am, ss);
axis(handles.axes1);
plot(ss, pressure);
plot(sqrt(ss), pressure);
legend(num2str([am17 am18 am19]))
