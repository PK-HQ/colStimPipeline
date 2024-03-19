function imagesColumnar=columnarFilter(TS,images)
% Bandpass filter for columns
bandpassSF=[.8 3];
imagesColumnar=FilterFermi2D(images(:,:,:), bandpassSF(1), bandpassSF(2), TS.Header.Imaging.SizePxl);
end