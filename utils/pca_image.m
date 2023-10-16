clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
clear;  % Erase all existing variables. Or clearvars if you want.
workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;
fontSize = 20;

% Check that user has the Image Processing Toolbox installed.
% This is only needed for supplying demo images.
hasIPT = license('test', 'image_toolbox');
if ~hasIPT
	% User does not have the toolbox installed.
	message = sprintf('Sorry, but you do not seem to have the Image Processing Toolbox.\nDo you want to try to continue anyway?');
	reply = questdlg(message, 'Toolbox missing', 'Yes', 'No', 'Yes');
	if strcmpi(reply, 'No')
		% User said No, so exit.
		return;
	end
end

% Check that user has the Statistics Toolbox installed.
% This is needed for the PCA function.
hasStatsToolbox = license('test', 'statistics_toolbox');
if ~hasStatsToolbox
	% User does not have the toolbox installed.
	message = sprintf('Sorry, but you do not seem to have the Statistics and Machine Learning Toolbox.\nDo you want to try to continue anyway?');
	reply = questdlg(message, 'Toolbox missing', 'Yes', 'No', 'Yes');
	if strcmpi(reply, 'No')
		% User said No, so exit.
		return;
	end
end

% Read in a standard MATLAB color demo image.
folder = fileparts(which('peppers.png')); % Determine where demo folder is (works with all versions).
button = menu('Use which demo image?', 'Peppers', 'Onion', 'Hands1', 'Colored Chips', 'He stain', 'Football', 'Fabric', 'Pears', 'Yellow Lily');
if button == 1
	baseFileName = 'peppers.png';
elseif button == 2
	baseFileName = 'onion.png';
elseif button == 3
	baseFileName = 'hands1.jpg';
elseif button == 4
	baseFileName = 'coloredChips.png';
elseif button == 5
	baseFileName = 'hestain.png';
elseif button == 6
	baseFileName = 'Football.jpg';
elseif button == 7
	baseFileName = 'fabric.png';
elseif button == 8
	baseFileName = 'pears.png';
elseif button == 9
	baseFileName = 'yellowlily.jpg';
end

%===============================================================================
% Read in a standard MATLAB color demo image.
folder = fileparts(which('peppers.png')); % Determine where demo folder is (works with all versions).
% Get the full filename, with path prepended.
fullFileName = fullfile(folder, baseFileName);
if ~exist(fullFileName, 'file')
	% Didn't find it there.  Check the search path for it.
	fullFileName = baseFileName; % No path this time.
	if ~exist(fullFileName, 'file')
		% Still didn't find it.  Alert user.
		errorMessage = sprintf('Error: %s does not exist.', fullFileName);
		uiwait(warndlg(errorMessage));
		return;
	end
end
rgbImage = imread(fullFileName);
% Get the dimensions of the image.  numberOfColorBands should be = 3.
[rows, columns, numberOfColorBands] = size(rgbImage);
% Display the original color image.
subplot(1, 2, 1);
imshow(rgbImage);
caption = sprintf('Original Color Image : %s', baseFileName);
title(caption, 'FontSize', fontSize, 'Interpreter', 'None');
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'Outerposition', [0, 0, 1, 1]);
drawnow;

% Extract the individual red, green, and blue color channels.
redChannel = rgbImage(:, :, 1);
greenChannel = rgbImage(:, :, 2);
blueChannel = rgbImage(:, :, 3);

% Display the red channel image.
fontSize = 14;
subplot(2, 6, 4);
imshow(redChannel);
title('Red Channel Image', 'FontSize', fontSize, 'Interpreter', 'None');
% Display the green channel image.
subplot(2, 6, 5);
imshow(greenChannel);
title('Green Channel Image', 'FontSize', fontSize, 'Interpreter', 'None');
% Display the blue channel image.
subplot(2, 6, 6);
imshow(blueChannel);
title('Blue Channel Image', 'FontSize', fontSize, 'Interpreter', 'None');

X = double(reshape(rgbImage, rows * columns, 3));

coeff = pca(X);

Itransformed = X * coeff;

pca1Image = reshape(Itransformed(:,1), rows, columns);
pca2Image = reshape(Itransformed(:,2), rows, columns);
pca3Image = reshape(Itransformed(:,3), rows, columns);

% Display the PCA 1 image.
subplot(2, 6, 10);
imshow(pca1Image, []);
title('PCA 1 Image', 'FontSize', fontSize, 'Interpreter', 'None');

% Display the PCA 2 image.
subplot(2, 6, 11);
imshow(pca2Image, []);
title('PCA 2 Image', 'FontSize', fontSize, 'Interpreter', 'None');

% Display the PCA 3 image.
subplot(2, 6, 12);
imshow(pca3Image, []);
title('PCA 3 Image', 'FontSize', fontSize, 'Interpreter', 'None');

