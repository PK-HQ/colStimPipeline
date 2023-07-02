function [labeledImage,areasThresholded]=getBlobAreas(varargin)
binaryImage=varargin{1};
if ~isempty(varargin{2})
    minAreaThreshold=varargin{2};
else
    minAreaThreshold=0;
end
% label img blobs
labeledImage = bwlabel(binaryImage);
% Measure the area
areas = struct2mat(regionprops(labeledImage, 'Area'));
areasThresholded=areas(areas>minAreaThreshold); %remove pesky small dots
end