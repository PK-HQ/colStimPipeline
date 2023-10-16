function images=columnarfilter(TS,images)
% Bandpass filter for columns
bandpassSF=[.8 3];
for cond=1:size(images,3)
    images(:,:,cond)=FilterFermi2D(images(:,:,cond), bandpassSF(1), bandpassSF(2), TS.Header.Imaging.SizePxl);
end
end