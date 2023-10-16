function transformedImage=transformImage(image,filenameStructCurrent)
%% This applies the transformation that morphs the reference orientation map to the current imaging window view of the cortex
% load the transformation, estimated between the vasculature image of reference to current imaging window view
load(filenameStructCurrent.transformParams) 
% Apply transformation to all images in image
for nImage=1:size(image,3)
    transformedImage(:,:,nImage)=double(imwarp(image(:,:,nImage),transformParams,'OutputView',imref2d(size(image,[1 2]))));
end
end