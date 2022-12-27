function varargout = GUI_Example(varargin)
%GUI_EXAMPLE M-file for GUI_Example.fig
%      GUI_EXAMPLE, by itself, creates a new GUI_EXAMPLE or raises the existing
%      singleton*.
%
%      H = GUI_EXAMPLE returns the handle to a new GUI_EXAMPLE or the handle to
%      the existing singleton*.
%
%      GUI_EXAMPLE('Property','Value',...) creates a new GUI_EXAMPLE using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to GUI_Example_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      GUI_EXAMPLE('CALLBACK') and GUI_EXAMPLE('CALLBACK',hObject,...) call the
%      local function named CALLBACK in GUI_EXAMPLE.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_Example

% Last Modified by GUIDE v2.5 26-Dec-2022 17:04:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_Example_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_Example_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before GUI_Example is made visible.
function GUI_Example_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for GUI_Example
handles.output = hObject;

set(hObject,'toolbar','figure');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_Example wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_Example_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

varargout{1} = handles.output;


% --- Executes on button press in open_file.
function open_file_Callback(hObject, eventdata, handles)
% hObject    handle to open_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname]= uigetfile({'*.jpg;*.tif;*.png;*.gif','All Image Files'; '*.*','All Files' },'Choose an image file ...');

if isequal(filename,0)
   warndlg('File not selected');
   return
end

image=imread(strcat(pathname,filename));

if not(size(image,3)==3)
   errordlg('The image must be colored');
   return
end

axes(handles.imag_orig);
imshow(image);

set(handles.imag_orig,'UserData',image);

% --- Executes on button press in saveas.
function saveas_Callback(hObject, eventdata, handles)
% hObject    handle to saveas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname]= uiputfile({'*.jpg;*.tif;*.png;*.gif','All Image Files'; '*.*','All Files' },'Save image file as ...');

if isequal(filename,0)
   warndlg('File not selected');
   return
end

image=get(handles.imag_processed,'UserData');

imwrite(image,strcat(pathname,filename));

% --- Executes on button press in new_window.
function new_window_Callback(hObject, eventdata, handles)
% hObject    handle to new_window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

image=get(handles.imag_processed,'UserData');

figure
imshow(image);

% --- Executes on button press in copy_orig.
function copy_orig_Callback(hObject, eventdata, handles)
% hObject    handle to copy_orig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

image=get(handles.imag_processed,'UserData');

axes(handles.imag_orig);
imshow(image);

set(handles.imag_orig,'UserData',image);

% --- Executes on button press in apply.
function apply_Callback(hObject, eventdata, handles)
% hObject    handle to apply (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

image=get(handles.imag_orig,'UserData');
image_d=double(image)/255;

image_proc=(image_d(:,:,1)+image_d(:,:,2)+image_d(:,:,3))/3; % intensity component

axes(handles.imag_processed);
imshow(image_proc);

set(handles.imag_processed,'UserData',image_proc);


% --- Executes on button press in capture.
function capture_Callback(hObject, eventdata, handles)
% hObject    handle to capture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

res_values=cellstr(get(handles.capture_resolution,'String'));
res=res_values{get(handles.capture_resolution,'Value')};

image=imread (['http://visionartificia.gnd.upv.es/axis-cgi/jpg/image.cgi?resolution=' res]);

axes(handles.imag_orig);
imshow(image);

set(handles.imag_orig,'UserData',image);


% --- Executes on selection change in capture_resolution.
function capture_resolution_Callback(hObject, eventdata, handles)
% hObject    handle to capture_resolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns capture_resolution contents as cell array
%        contents{get(hObject,'Value')} returns selected item from capture_resolution

% --- Executes during object creation, after setting all properties.
function capture_resolution_CreateFcn(hObject, eventdata, handles)
% hObject    handle to capture_resolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
