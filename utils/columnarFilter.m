function images=columnarFilter(TS,images,trialOutcomeType)
% Bandpass filter for columns
bandpassSF=[.8 3];
switch trialOutcomeType
    case 'average'
        %images=images.average;
    case 'averageColumn'
        images.(trialOutcomeType)=FilterFermi3D(images.average, bandpassSF(1), bandpassSF(2), TS.Header.Imaging.SizePxl);
end
%{
tic
for cond=1:size(images,3)
    imagesColumnar(:,:,cond)=FilterFermi2D(images(:,:,cond), bandpassSF(1), bandpassSF(2), TS.Header.Imaging.SizePxl);
end
toc
%}
end