function [angle, amplitude, angleAverage, amplitudeAverage]=calculateProjectionVector(image, referenceImage, orientation)
%% Load ort map, but separated into each color
% init
nContrasts=size(image,3);
orientations=repmat(orientation,1,nContrasts);
angle=nan(1, nContrasts);
amplitude=nan(1, nContrasts);
x=nan(1, nContrasts);
y=nan(1, nContrasts);

%% Calculate angle and amplitude of image projected into reference image space

% Rescale
normRange=[-1 1];
referenceImageRescaled=rescale(referenceImage,normRange(1),normRange(2));

% Mask reference map to match region used in image
mask=double(~isnan(image(:,:,1)));
mask(mask==0)=NaN;
% Resize if necessary
if size(referenceImage,1) > size(mask,1) % if binned image
    mask=imresize(mask,size(referenceImage,[1 2]) ,'bilinear');
end
referenceImageRescaled=referenceImageRescaled.*mask;

for ort=1:nContrasts
    % Dot product of image with reference image
    projectedMap=image(:,:,ort).*referenceImageRescaled;
    % Sum all pixels, and normalize by converting the summed value into a percentage of ort map
    imageSum=sum(projectedMap,'all','omitnan'); %exclude nans
    referenceSum=abs(sum(referenceImageRescaled,'all','omitnan')); %exclude nans
    normSum=imageSum*100./referenceSum;
    % Get angle and amplitude of vector
    angle(ort)=orientations(ort); amplitude(ort)=normSum;
    % Convert to 2D vector coordinates
    x(ort)=cos(angle(ort)) * amplitude(ort); %x
    y(ort)=sin(angle(ort)) * amplitude(ort); %y
end
% Get vector average (angle and amplitude)
angleAverage=atan(mean(y)/mean(x));
amplitudeAverage=mean(y)/sin(angleAverage);
end