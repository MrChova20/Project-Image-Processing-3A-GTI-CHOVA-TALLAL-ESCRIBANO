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

%--------------------------------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------------------------

% --- Outputs from this function are returned to the command line.
function varargout = GUI_Example_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

varargout{1} = handles.output;
%--------------------------------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------------------------



%--------------------------------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------------------------

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

%
axes(handles.imag_orig);
imshow(image);
set(handles.imag_orig,'UserData',image);


axes(handles.imag_processed);
imshow(image);
set(handles.imag_processed,'UserData',image);

%--------------------------------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------------------------

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

%--------------------------------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------------------------



%--------------------------------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------------------------
% --- Executes on button press in apply.
function apply_Callback(hObject, eventdata, handles)
% hObject    handle to apply (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

image1=get(handles.imag_processed,'UserData');

%% Important parameters

imageBlurSigma = 2; % Amount to denoise input image
MinHoughPeakDistance = 5; % Distance between peaks in Hough transform angle detection
HoughConvolutionLength = 40; % Length of line to use to detect bone regions
HoughConvolutionDilate = 2; % Amount to dilate kernel for bone detection
BreakLineTolerance = 0.25; % Tolerance for bone end detection
breakPointDilate = 6; % Amount to dilate detected bone end points

%%%

image_gray = (rgb2gray(image1)); % 1º FILTRADOOOOOOOOOOOOOOOOOOOOOOOOO

%figure(2)% MOSTRAR IMAGE
%imshow(image_gray);% MOSTRAR IMAGE
%title('Gray Scale X Ray Image');% MOSTRAR IMAGE

image_filtered = imfilter(image_gray, fspecial('gaussian', 10, imageBlurSigma), 'symmetric'); % Denoise  2º FILTRADOOOOOOOOOOOOOOOOOOOOOOOOO

%figure(3)% MOSTRAR IMAGE
%imshow(image_filtered);% MOSTRAR IMAGE
%title('denoised Gray Scale X Ray image');% MOSTRAR IMAGE

% Do edge detection to find bone edges in image
% Filter out all but the two longest lines
% This feature may need to be changed if break is not in middle of bone 
% ESO SIGNIFICA Q SI ES UNA FISURA HAY Q CAMBAIR ESTO
boneEdges = edge(image_filtered, 'canny');% 3º FILTRADOOOOOOOOOOOOOOOOOOOOOOOOO

%figure(4)% MOSTRAR IMAGE
%imshow(boneEdges);% MOSTRAR IMAGE
%title('Edges of the bones');% MOSTRAR IMAGE

boneEdges1 = bwmorph(boneEdges, 'close');% 4º FILTRADOOOOOOOOOOOOOOOOOOOOOOOOO

%figure(5)% MOSTRAR IMAGE
%imshow(boneEdges1);% MOSTRAR IMAGE
%title('Morphological operation on Edges of the bones ');% MOSTRAR IMAGE

edgeRegs = regionprops(boneEdges1, 'Area', 'PixelIdxList');%  CREAMOSSSSSSSSSSSSSSSS GRAFICA 6 EN PROCESOOOOOOOOOOOOOO
AreaList = sort(vertcat(edgeRegs.Area), 'descend');
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

%--------------------------------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------------------------
%--------------------------------------------------------------------------------------------------------------------------------

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
