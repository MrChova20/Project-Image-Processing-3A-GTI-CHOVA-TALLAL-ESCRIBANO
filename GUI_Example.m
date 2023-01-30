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
% hObject    handle to apply (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

image1=get(handles.imag_processed,'UserData');

%% Important parameters

imageBlurSigma = 2; %cantidad de suavizado que se aplicará a la imagen de entrada.
MinHoughPeakDistance = 5; % es la distancia mínima permitida entre los picos en la detección de ángulos en la transformada de Hough.
HoughConvolutionLength = 40; %  tamaño de la línea que se utilizará para detectar regiones óseas en la imagen.
HoughConvolutionDilate = 2; % cantidad de dilatación que se aplicará al kernel para la detección de huesos.
BreakLineTolerance = 0.25; % la tolerancia para la detección del final de los huesos.
breakPointDilate = 6; %es la cantidad de dilatación que se aplicará a los puntos detectados en el final de los huesos.

%%%

image_gray = (rgb2gray(image1)); % 1º se convierte a escala de grises con la función "rgb2gray"



image_filtered = imfilter(image_gray, fspecial('gaussian', 10, imageBlurSigma), 'symmetric'); %aplica un filtro gaussiano a la imagen en escala de grises "image_gray" y guarda el resultado en la variable "image_filtered". El filtro gaussiano se crea con la función "fspecial" y tiene un tamaño de 10x10 y un sigma (desviación estándar) especificado por la variable "imageBlurSigma". El parámetro "symmetric" indica que se debe
                                                                                                %utilizar una forma de simetría en el filtrado para mejorar la conservación de los bordes.


boneEdges = edge(image_filtered, 'canny');%  detecta los bordes en la imagen "image_filtered" utilizando el algoritmo de Canny. El resultado es almacenado en la variable "boneEdges".



boneEdges1 = bwmorph(boneEdges, 'close');%Esta linea de codigo utiliza la funcion "bwmorph" de Matlab, la cual permite realizar operaciones morfologicas en imágenes binarias. En este caso, se utiliza la operación "close" para cerrar pequeñas brechas o agujeros en los bordes detectados previamente en la imagen y unificar los contornos conectados. La salida se almacena en la variable "boneEdges1".



edgeRegs = regionprops(boneEdges1, 'Area', 'PixelIdxList');%  Esta línea de código utiliza la función "regionprops" para identificar las regiones conectadas en la imagen "boneEdges1" y calcular propiedades de estas regiones. La propiedad "Area" se refiere al número de píxeles en cada región, mientras que la propiedad "PixelIdxList" es una lista de índices de los píxeles en cada región. Es decir, "edgeRegs" es una estructura que contiene información sobre las regiones conectadas en la imagen "boneEdges1".
AreaList = sort(vertcat(edgeRegs.Area), 'descend'); %oncatena y ordena en orden descendente los valores de la propiedad "Area" de la estructura "edgeRegs" y los almacena en la matriz "AreaList".
edgeRegs(~ismember(vertcat(edgeRegs.Area), AreaList(1:2))) = [];
edgeimage = zeros(size(image_filtered, 1), size(image_filtered,2));%    GRAFICA 6 EN PROCESOOOOOOOOOOOOOO
edgeimage(vertcat(edgeRegs.PixelIdxList)) = 1;

% Do hough transform on edge image to find angles at which bone pieces are
% found
% Use max value of Hough transform vs angle to find angles at which lines
% are oriented.  If there is more than one major angle contribution there
% will be two peaks detected but only one peak if there is only one major
% angle contribution (ie peaks here = number of located bones = Number of
% breaks + 1)
[H,T,R] = hough(edgeimage,'RhoResolution',1,'Theta',-90:2:89.5);
maxHough = max(H, [], 1);
HoughThresh = (max(maxHough) - min(maxHough))/2 + min(maxHough);
[~, HoughPeaks] = findpeaks(maxHough,'MINPEAKHEIGHT',HoughThresh, 'MinPeakDistance', MinHoughPeakDistance);

% Plot Hough detection results
%figure(6)% MOSTRAR GRAFICA
%plot(T, maxHough);
hold on
%plot([min(T) max(T)], [HoughThresh, HoughThresh], 'g');
%plot(T(HoughPeaks), maxHough(HoughPeaks), 'rx', 'MarkerSize', 12, 'LineWidth', 2);
hold off
%xlabel('Theta Value'); ylabel('Max Hough Transform');
%legend({'Max Hough Transform', 'Hough Peak Threshold', 'Detected Peak'});
%title('Hough Detection Plot : Max Hough transform vs Theta');% MOSTRAR GRAFICA



%%2nd CT 8th sem

% Locate site of break
if numel(HoughPeaks) > 1;
    BreakStack = zeros(size(image_filtered, 1), size(image_filtered, 2), numel(HoughPeaks));
    % Convolute edge image with line of detected angle from hough transform
    for m = 1:numel(HoughPeaks);

        boneKernel = strel('line', HoughConvolutionLength, T(HoughPeaks(m)));
        kern = double(bwmorph(boneKernel.getnhood(), 'dilate', HoughConvolutionDilate));
        BreakStack(:,:,m) = imfilter(edgeimage, kern).*edgeimage;
        %figure(7) %MOSTRAR  IMAGEN
        %imshow(BreakStack(:,:,m));
        
    end

    % Take difference between convolution images.  Where this crosses zero
    % (within tolerance) should be where the break is.  Have to filter out
    % regions elsewhere where the bone simply ends.
    
    
    brimage = abs(diff(BreakStack, 1, 3)) < BreakLineTolerance*max(BreakStack(:)) & edgeimage > 0;
    [BpY, BpX] = find(abs(diff(BreakStack, 1, 3)) < BreakLineTolerance*max(BreakStack(:)) & edgeimage > 0);
    brimage = bwmorph(brimage, 'dilate', breakPointDilate);
    %figure(8);% MOSTRAR IMAGE
    %imshow(brimage);% MOSTRAR IMAGE
    brReg = regionprops(brimage, 'Area', 'MajorAxisLength', 'MinorAxisLength', ...
        'Orientation', 'Centroid');
    brReg(vertcat(brReg.Area) ~= max(vertcat(brReg.Area))) = [];

    % Calculate bounding ellipse
    % CIRCULOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    brReg.EllipseCoords = zeros(100, 2);
    t = linspace(0, 2*pi, 100);
    brReg.EllipseCoords(:,1) = brReg.Centroid(1) + brReg.MajorAxisLength/2*cos(t - brReg.Orientation);
    brReg.EllipseCoords(:,2) = brReg.Centroid(2) + brReg.MinorAxisLength/2*sin(t - brReg.Orientation);

else
    brReg = [];      %% No Fracture points are there

end

% Draw ellipse around break location
%figure(9)% MOSTRAR IMAGE
hold on
colormap('gray')
if ~isempty(brReg)
    plot(brReg.EllipseCoords(:,1), brReg.EllipseCoords(:,2), 'r');
end
hold off
    
axes(handles.imag_orig);
imshow(image1);
set(handles.imag_orig,'UserData',image1);


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
