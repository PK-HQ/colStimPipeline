function [similarityScore]=showGreenImage(dsMovingImg,dsFixedImg, plotOrNot)

method='hessian';
similarityScore=NaN;
%% Plotting
if plotOrNot==0
    set(groot,'defaultFigureVisible','off')
else
    set(groot,'defaultFigureVisible','on')
end

%% Get filenames
structMovingImg=generateFilenames(dsMovingImg);
structFixedImg=generateFilenames(dsFixedImg);

pdfFilename=[structFixedImg.plotPath 'M' dsFixedImg.monkeyNo 'R' dsMovingImg.date 'S' dsFixedImg.date '.pdf'];

%% Get coregistration params for greenImageReference (to be transformed) to greenImageSession (anchor image)
% load images, rescale
movingImg=(double(imread(structMovingImg.greenImg)));%imread(convertCharsToStrings(filenameStructRef.greenImg));
fixedImg=(double(imread(structFixedImg.greenImg)));%imread(convertCharsToStrings(filenameStructSession.greenImg));

imgSizeX=[1 size(fixedImg,1)];
imgSizeY=[1 size(fixedImg,2)];

nexttile
imagesc(imfuse(movingImg,fixedImg,'ColorChannels','red-cyan')); axis square
addPix2MM(min(imgSizeX),max(imgSizeX),min(imgSizeY),max(imgSizeY),1,2,2);
title('Raw')

%% Process vasculature
switch method
    case {'none'}
    case {'fermifilter'}
            movingImgPost = FilterFermi2D(movingImg,0.01,inf,1);
            fixedImgPost = FilterFermi2D(fixedImg,0.01,inf,1);
    case {'threshold'}
        %movingImg>0
    case {'hessian'}
        % ----- Mask positive vasculature -----
        %       .FrangiScaleRange : The range of sigmas used, default [1 8]
        %       .FrangiScaleRatio : Step size between sigmas, default 2
        %       .FrangiBetaOne : Frangi correction constant, default 0.5
        %       .FrangiBetaTwo : Frangi correction constant, default 15
        %       .BlackWhite : Detect black ridges (default) set to true, for
        %                       white ridges set to false.
        options.FrangiScaleRange=[6 9];%[1.5 8];
        options.FrangiScaleRatio=.5;
        options.FrangiBetaOne=2;
        options.FrangiBetaTwo=7;
        options.BlackWhite=true;
        
        
        movingImgPost = (FrangiFilter2D(movingImg,options));
        fixedImgPost = (FrangiFilter2D(fixedImg,options));

        %mask=double(movingImg<0.01);
        %mask(mask==0)=NaN;
        
end

%display
nexttile
imagesc(imfuse(movingImgPost,fixedImgPost,'ColorChannels','red-cyan')); axis square
addPix2MM(min(imgSizeX),max(imgSizeX),min(imgSizeY),max(imgSizeY),1,2,2);
title('Post-processing')

%transform type: auto, similarity
nexttile
method='Auto'; %or 'Manual'
[movingImgPostCoreg,coregStats,~]=coregisterGreenImages(movingImgPost,fixedImgPost,method,structFixedImg);
imagesc(imfuse(movingImgPostCoreg,fixedImgPost,'ColorChannels','red-cyan')); axis square
title({['Post-' method], ...
    sprintf(', r=%.2g (T-x=%.2g, T-y=%.2g, Rot=%.2g, Scale=%.2g)',...
    coregStats.Similarity,coregStats.Translate(1),coregStats.Translate(2),coregStats.Rotation,coregStats.Rotation)});
addPix2MM(min(imgSizeX),max(imgSizeX),min(imgSizeY),max(imgSizeY),3,2,2);
similarityScore=coregStats.Similarity;


%% Plotting
if plotOrNot==0
    set(groot,'defaultFigureVisible','on')
end

end