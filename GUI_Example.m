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
%imagePrimeraVista=imread(strcat(pathname,filename));

if not(size(image,3)==3)
   errordlg('The image must be colored');
   return
end

%
axes(handles.imag_orig);
imshow(image);
set(handles.imag_orig,'UserData',image);


axes(handles.imag_processed);
imshow(image);
set(handles.imag_processed,'UserData',image);


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



% --- Executes on button press in apply.
function apply_Callback(hObject, eventdata, handles)

%% Important parameters
% Cargar la imagen del hueso
bone_image = get(handles.imag_processed,'UserData');
% Convertir a escala de grises
%bone_image = imread('30.jpg');

% Convertir a escala de grises
im_gray = rgb2gray(bone_image);

% Umbralizacion con umbral de 0.9 (va con la mayoria de los imagenes)
threshold = graythresh(im_gray);
im_bw = im2bw(im_gray, 0.9);


% Usar la función bwlabel para identificar la conectividad
[labeled, num_objects] = bwlabel(im_bw, 8);
% Determinar si hay al menos una fractura
fracture_detected = num_objects > 0;

% Mostrar la imagen original con un círculo rojo o azul sobre la fractura

imshow(bone_image); hold on;
if fracture_detected
    % Encontrar el centro de la primera área conectada
    stats = regionprops(labeled, 'Centroid');
    center = stats(1).Centroid;
    
    % Dibujar un círculo rojo en el centro de la fractura
    radius = 50;
    viscircles(center, radius, 'Color', 'r');
else
   msgbox('No hay fracturas', 'Pop-Up Message');
  
end

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
