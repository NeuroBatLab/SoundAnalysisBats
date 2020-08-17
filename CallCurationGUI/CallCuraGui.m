function varargout = CallCuraGui(varargin)
% CALLCURAGUI MATLAB code for CallCuraGui.fig
%      CALLCURAGUI, by itself, creates a new CALLCURAGUI or raises the existing
%      singleton*.
%
%      H = CALLCURAGUI returns the handle to a new CALLCURAGUI or the handle to
%      the existing singleton*.
%
%      CALLCURAGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CALLCURAGUI.M with the given input arguments.
%
%      CALLCURAGUI('Property','Value',...) creates a new CALLCURAGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CallCuraGui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CallCuraGui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CallCuraGui

% Last Modified by GUIDE v2.5 17-Aug-2020 10:21:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CallCuraGui_OpeningFcn, ...
                   'gui_OutputFcn',  @CallCuraGui_OutputFcn, ...
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


% --- Executes just before CallCuraGui is made visible.
function CallCuraGui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CallCuraGui (see VARARGIN)

% Choose default command line output for CallCuraGui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CallCuraGui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CallCuraGui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function redoEditVoc_Callback(hObject, eventdata, handles)
% hObject    handle to redoEditVoc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of redoEditVoc as text
%        str2double(get(hObject,'String')) returns contents of redoEditVoc as a double


% --- Executes during object creation, after setting all properties.
function redoEditVoc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to redoEditVoc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in micRadioh.
function micRadioh_Callback(hObject, eventdata, handles)
% hObject    handle to micRadioh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of micRadioh


% --- Executes on button press in logRadioh.
function logRadioh_Callback(hObject, eventdata, handles)
% hObject    handle to logRadioh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of logRadioh



function redoEditSet_Callback(hObject, eventdata, handles)
% hObject    handle to redoEditSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of redoEditSet as text
%        str2double(get(hObject,'String')) returns contents of redoEditSet as a double


% --- Executes during object creation, after setting all properties.
function redoEditSet_CreateFcn(hObject, eventdata, handles)
% hObject    handle to redoEditSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function sliderLeft_Callback(hObject, eventdata, handles)
% hObject    handle to sliderLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function sliderLeft_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliderRight_Callback(hObject, eventdata, handles)
% hObject    handle to sliderRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function sliderRight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in start.
function start_Callback(hObject, eventdata, handles)
% hObject    handle to start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in playLogEval.
function playLogEval_Callback(hObject, eventdata, handles)
% hObject    handle to playLogEval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in playMicEval.
function playMicEval_Callback(hObject, eventdata, handles)
% hObject    handle to playMicEval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
