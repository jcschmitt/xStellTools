function varargout = current_gui(varargin)
% CURRENT_GUI MATLAB code for current_gui.fig
%      CURRENT_GUI, by itself, creates a new CURRENT_GUI or raises the existing
%      singleton*.
%
%      H = CURRENT_GUI returns the handle to a new CURRENT_GUI or the handle to
%      the existing singleton*.
%
%      CURRENT_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CURRENT_GUI.M with the given input arguments.
%
%      CURRENT_GUI('Property','Value',...) creates a new CURRENT_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before current_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to current_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help current_gui

% Last Modified by GUIDE v2.5 13-Feb-2018 19:42:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @current_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @current_gui_OutputFcn, ...
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


% --- Executes just before current_gui is made visible.
function current_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to current_gui (see VARARGIN)

% Choose default command line output for current_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes current_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = current_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% After an update, Call the 'Update Plots' button.
set(handles.edit1, 'String', num2str(get(hObject, 'Value')));
pushbutton1_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% After an update, Call the 'Update Plots' button.
set(handles.edit2, 'String', num2str(get(hObject, 'Value')));
pushbutton1_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% After an update, Call the 'Update Plots' button.
set(handles.edit3, 'String', num2str(get(hObject, 'Value')));
pushbutton1_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% After an update, Call the 'Update Plots' button.
set(handles.edit4, 'String', num2str(get(hObject, 'Value')));
pushbutton1_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider5_Callback(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% After an update, Call the 'Update Plots' button.
set(handles.edit5, 'String', num2str(get(hObject, 'Value')));
pushbutton1_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function slider5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider6_Callback(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% After an update, Call the 'Update Plots' button.
set(handles.edit6, 'String', num2str(get(hObject, 'Value')));
pushbutton1_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function slider6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider8_Callback(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% After an update, Call the 'Update Plots' button.
set(handles.edit8, 'String', num2str(get(hObject, 'Value')));
pushbutton1_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function slider8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
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

% Get the coeeffiencient of the current profile from the AC(2) ... AC(8)
% buttons
% AC0 and AC1 are always 0 and 1 (although, maybe 1 could be modified...)
AC0 = get(handles.slider0, 'Value');
AC1 = get(handles.slider1, 'Value');
AC2 = get(handles.slider2, 'Value');
AC3 = get(handles.slider3, 'Value');
AC4 = get(handles.slider4, 'Value');
AC5 = get(handles.slider5, 'Value');
AC6 = get(handles.slider6, 'Value');
AC7 = get(handles.slider7, 'Value');
AC8 = get(handles.slider8, 'Value');

ns = get(handles.slider_ns, 'Value');
s = linspace(1/ns, 1, ns);
s_ext = 1+s(1:(floor(ns/5)));

profile_type = handles.listbox1.String{handles.listbox1.Value};

switch profile_type
    % Pressure profiles
    %
    
    case 'pp_ne_b'
        PP_NE_B_array = [AC0 AC1 AC2 AC3];
        PP_NE_min = AC4;
        legend_str = num2str([AC0 AC1 AC2 AC3 AC4]);
        pp_ne_b_profile = pp_ne_b(s, PP_NE_B_array, PP_NE_min);
        pp_ne_b_profile_ext = pp_ne_b(s_ext, PP_NE_B_array, PP_NE_min);
        ne_profile_plot = pp_ne_b_profile;
        ne_profile_ext_plot = pp_ne_b_profile_ext;
        grad_ne_profile = [nan diff(ne_profile_plot)];
        grad_ne_profile_ext = [nan diff(ne_profile_ext_plot)];
        grad_ne_profile_plot = grad_ne_profile / max(abs(grad_ne_profile));
        grad_ne_profile_ext_plot = grad_ne_profile_ext / max(abs(grad_ne_profile_ext));
        
        % set axes1 as active axes
        axes(handles.axes1);
        cla;
        box on;hold on;
        plot(s, ne_profile_plot, 'k-');
        plot(s_ext, ne_profile_ext_plot, 'k:');
        title('N_e(s)')
        legend(legend_str)
        
        % set axes2 as active axes
        axes(handles.axes2);
        cla;
        box on;hold on;
        plot(s, grad_ne_profile_plot, 'k-');
        plot(s_ext, grad_ne_profile_ext_plot, 'k:');
        title('grad N_e(s)')
        xlabel('s');
        
        % set axes3 as active axes
        axes(handles.axes3);
        cla;
        box on;hold on;
        plot(sqrt(s), ne_profile_plot, 'k-');
        plot(sqrt(s_ext), ne_profile_ext_plot, 'k:');
        title('N_e(s)')
        legend(legend_str)
        
        % set axes4 as active axes
        axes(handles.axes4);
        cla;
        box on;hold on;
        plot(sqrt(s), grad_ne_profile_plot, 'k-');
        plot(sqrt(s_ext), grad_ne_profile_ext_plot, 'k:');
        title('grad N_e(s)')
        %legend_str = num2str(AC_array);
        %legend(legend_str)
        xlabel('\rho');
    case 'pedestal'
        AM = zeros(1, 20);
        AM(15) = AC0;
        AM(14) = AC1;
        AM(13) = AC2;
        AM(12) = AC3;
        AM(11) = AC4;
        AM(16) = AC5;
        AM(17) = AC6;
        AM(18) = AC7;
        AM(19) = AC8;
        
            %DO i = 15, LBOUND(am,1), -1
            %  pmass = x*pmass + am(i)
            %END DO
        pmass = 0;
        for ii = 15:-1:11
            pmass = s .* pmass + AM(ii);
        end
        
            %i = 16
            %IF (am(i+3) .le. 0._dp ) THEN
            %   am(i:i+4) = 0
            %   am(i+3) = 1.0e30_dp
            % ELSE
            %   am(i+4) = 1.0_dp/(TANH(2*am(i+2)/am(i+3))-TANH(2*(am(i+2)-1)/am(i+3)))
            %END IF

        if AM(16)<= 0
            AM(20) = 0;
            AM(19) = 1e30;
        else
            AM(20) = 1 ./ (tanh(2*AM(18) / AM(19) ) - tanh(2*(AM(18)-1)/AM(19)));
        end
        
            %pmass = pmass + am(i+4) * am(i+1) * ( TANH( 2*(am(i+2)-SQRT(x))/am(i+3) ) - TANH( 2*(am(i+2)-1._dp)  /am(i+3) ) )
        edge_part = AM(20) * AM(17) * ( tanh( 2*(AM(18) - sqrt(s))/AM(19) ) - tanh( 2*(AM(18) -1) / AM(19) ) );
        pmass_profile = pmass + edge_part;
        
        legend_str = num2str([AC0 AC1 AC2 AC3 AC4 AC5 AC6 AC7 AC8]);
        ptotal_profile_plot = pmass_profile / pmass_profile(1);
        ppoly_profile_plot = pmass / pmass_profile(1);
        pedge_profile_plot = edge_part / pmass_profile(1);
        grad_ptotal_profile = [nan diff(pmass_profile)];
        grad_ppoly_profile = [nan diff(pmass)];
        grad_pedge_profile = [nan diff(edge_part)];
        grad_ptotal_profile_plot = grad_ptotal_profile / max(abs(grad_ptotal_profile));
        grad_ppoly_profile_plot = grad_ppoly_profile / max(abs(grad_ptotal_profile));
        grad_pedge_profile_plot = grad_pedge_profile / max(abs(grad_ptotal_profile));
        
        % set axes1 as active axes
        axes(handles.axes1);
        cla;
        box on;hold on;
        plot(s, ptotal_profile_plot, 'k-');
        plot(s, ppoly_profile_plot, 'k--');
        plot(s, pedge_profile_plot, 'k:');
        title('P_{total}(s)')
        ylim([-.1 1.1])
        legend(legend_str)
        grid;
        
        % set axes2 as active axes
        axes(handles.axes2);
        cla;
        box on;hold on;
        plot(s, grad_ptotal_profile_plot, 'k-');
        plot(s, grad_ppoly_profile_plot, 'k--');
        plot(s, grad_pedge_profile_plot, 'k:');
        title('grad P_{total}(s)')
        xlabel('s');
        grid;
        
        % set axes3 as active axes
        axes(handles.axes3);
        cla;
        box on;hold on;
        plot(sqrt(s), ptotal_profile_plot, 'k-');
        plot(sqrt(s), ppoly_profile_plot, 'k--');
        plot(sqrt(s), pedge_profile_plot, 'k:');
        title('P_{total}(\rho)')
        ylim([-.1 1.1])
        legend(legend_str)
        grid;
        
        % set axes4 as active axes
        axes(handles.axes4);
        cla;
        box on;hold on;
        plot(sqrt(s), grad_ptotal_profile_plot, 'k-');
        plot(sqrt(s), grad_ppoly_profile_plot, 'k--');
        plot(sqrt(s), grad_pedge_profile_plot, 'k:');
        title('grad P_{total}(\rho)')
        %legend_str = num2str(AC_array);
        %legend(legend_str)
        xlabel('\rho');
        grid;

        
    % Current profiles
    %
    
    case 'sum_atan'
        AC_array = [AC0 AC1 AC2 AC3 AC4 AC5 AC6 AC7 AC8];
        AC_array_p1 = [0 AC1 AC2 AC3 AC4 ];
        AC_array_p2 = [0 AC5 AC6 AC7 AC8];
        I_profile = sum_atan(AC_array, s) ;  % 14
        I_profile_p1 = sum_atan(AC_array_p1, s) ;  % 14
        I_profile_p2 = sum_atan(AC_array_p2, s) ;  % 14
        
        % I want to normalize the current at the edge
        I_profile_plot = I_profile / I_profile(end);
        I_profile_plot_p1 = I_profile_p1 / I_profile(end);
        I_profile_plot_p2 = I_profile_p2 / I_profile(end);
        
        J_profile = [nan diff(I_profile)];
        J_profile_p1 = [nan diff(I_profile_p1)];
        J_profile_p2 = [nan diff(I_profile_p2)];
        J_profile_plot = J_profile / max(abs(J_profile));
        J_profile_plot_p1 = J_profile_p1 / max(abs(J_profile));
        J_profile_plot_p2 = J_profile_p2 / max(abs(J_profile));
        
        
        % set axes1 as active axes
        axes(handles.axes1);
        cla;
        box on;hold on;
        plot(s, I_profile_plot, 'k-');
        plot(s, I_profile_plot_p1, 'r^--');
        plot(s, I_profile_plot_p2, 'bv--');
        title('I(s)')
        legend_str = num2str(AC_array);
        legend(legend_str)
        
        % set axes2 as active axes
        axes(handles.axes2);
        cla;
        box on;hold on;
        plot(s, J_profile_plot, 'k-');
        plot(s, J_profile_plot_p1, 'r^--');
        plot(s, J_profile_plot_p2, 'bv--');
        title('J(s)')
        %legend_str = num2str(AC_array);
        %legend(legend_str)
        
        xlabel('s');
        
        % set axes3 as active axes
        axes(handles.axes3);
        cla;
        box on;hold on;
        plot(sqrt(s), I_profile_plot, 'o');
        plot(sqrt(s), I_profile_plot_p1, 'r^--');
        plot(sqrt(s), I_profile_plot_p2, 'bv--');
        title('I(\rho)')
        legend_str = num2str(AC_array);
        legend(legend_str)
        
        % set axes4 as active axes
        axes(handles.axes4);
        cla;
        box on;hold on;
        plot(sqrt(s), J_profile_plot, 'o');
        plot(sqrt(s), J_profile_plot_p1, 'r^--');
        plot(sqrt(s), J_profile_plot_p2, 'bv--');
        title('J(\rho)')
        %legend_str = num2str(AC_array);
        %legend(legend_str)
        xlabel('\rho');
end

% --- Executes on slider movement.
function slider_NS_Callback(hObject, eventdata, handles)
% hObject    handle to slider_NS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_NS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_NS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_ns_Callback(hObject, eventdata, handles)
% hObject    handle to slider_ns (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_ns_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_ns (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit0 as text
%        str2double(get(hObject,'String')) returns contents of edit0 as a double

%If edited, set the corresponding slider bar value, and call the callback
%for the update plots push button
set(handles.slider1, 'Value', str2num(get(hObject, 'String')));
pushbutton1_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit0 as text
%        str2double(get(hObject,'String')) returns contents of edit0 as a double
set(handles.slider2, 'Value', str2num(get(hObject, 'String')));
pushbutton1_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit0 as text
%        str2double(get(hObject,'String')) returns contents of edit0 as a double

set(handles.slider3, 'Value', str2num(get(hObject, 'String')));
pushbutton1_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit0 as text
%        str2double(get(hObject,'String')) returns contents of edit0 as a double
set(handles.slider4, 'Value', str2num(get(hObject, 'String')));
pushbutton1_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit0 as text
%        str2double(get(hObject,'String')) returns contents of edit0 as a double

set(handles.slider5, 'Value', str2num(get(hObject, 'String')));
pushbutton1_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit0 as text
%        str2double(get(hObject,'String')) returns contents of edit0 as a double

set(handles.slider6, 'Value', str2num(get(hObject, 'String')));
pushbutton1_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit0 as text
%        str2double(get(hObject,'String')) returns contents of edit0 as a double

set(handles.slider7, 'Value', str2num(get(hObject, 'String')));
pushbutton1_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on slider movement.
function slider0_Callback(hObject, eventdata, handles)
% hObject    handle to slider0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% After an update, Call the 'Update Plots' button.
set(handles.edit0, 'String', num2str(get(hObject, 'Value')));
pushbutton1_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function slider0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end






function edit0_Callback(hObject, eventdata, handles)
% hObject    handle to edit0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit0 as text
%        str2double(get(hObject,'String')) returns contents of edit0 as a double

set(handles.slider0, 'Value', str2num(get(hObject, 'String')));
pushbutton1_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider7_Callback(hObject, eventdata, handles)
% hObject    handle to slider7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
% After an update, Call the 'Update Plots' button.
set(handles.edit7, 'String', num2str(get(hObject, 'Value')));
pushbutton1_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function slider7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes during object deletion, before destroying properties.
function edit7_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
