
clear all;
close all;
%1st CT 8th sem
%image = imread('http://i.stack.imageur.com/mHo7s.jpg');
%image = imread('frac3.jpg');
image = imread('frac3.jpg');%LLEGA IMAGEN

figure(1)% MOSTRAR IMAGE
imshow(image);% MOSTRAR IMAGE
title('Input X Ray Image');% MOSTRAR IMAGE

%% Important parameters

imageBlurSigma = 2; % Amount to denoise input image
MinHoughPeakDistance = 5; % Distance between peaks in Hough transform angle detection
HoughConvolutionLength = 40; % Length of line to use to detect bone regions
HoughConvolutionDilate = 2; % Amount to dilate kernel for bone detection
BreakLineTolerance = 0.25; % Tolerance for bone end detection
breakPointDilate = 6; % Amount to dilate detected bone end points

%%%

image_gray = (rgb2gray(image)); % 1º FILTRADOOOOOOOOOOOOOOOOOOOOOOOOO

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
    brReg.EllipseCoords = zeros(100, 2);
    t = linspace(0, 2*pi, 100);
    brReg.EllipseCoords(:,1) = brReg.Centroid(1) + brReg.MajorAxisLength/2*cos(t - brReg.Orientation);
    brReg.EllipseCoords(:,2) = brReg.Centroid(2) + brReg.MinorAxisLength/2*sin(t - brReg.Orientation);

else
    brReg = [];      %% No Fracture points are there

end

% Draw ellipse around break location
%figure(9)% MOSTRAR IMAGE
imshow(image)% MOSTRAR IMAGE
hold on
colormap('gray')
if ~isempty(brReg)
    plot(brReg.EllipseCoords(:,1), brReg.EllipseCoords(:,2), 'r');
end
hold off
    