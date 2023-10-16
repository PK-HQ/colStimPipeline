function [columnarBitmaps,columnarPCAs]=coregisterBitmap2GreenImgV2(dataStructReference,dataStructSession,bitmapRef,VERpca,plotFlag,saveFlag)
%% Plotting
if plotFlag==0
    set(groot,'defaultFigureVisible','off')
else
    set(groot,'defaultFigureVisible','on')
end

%% Get filenames
filenameStructSession=generateFilenames(dataStructSession);
filenameStructRef=generateFilenames(dataStructReference);
pdfFilename=filenameStructSession.neurometricPDF; %pdfFilename='Y:/Chip/Chip20221102/coregGreen.pdf'%[filenameStructSession.plotPath 'M' dataStructSession.monkeyNo 'R' dataStructReference.date 'S' dataStructSession.date '.pdf'];

%% Load
if ispc
  mainPath='Y:/';
elseif contains(getenv('HOSTNAME'),'psy.utexas.edu')
  mainPath='/eslab/data/';
end
%{
% load green images
filenameReferenceY=[mainPath  dataStructReference.monkey '/' dataStructReference.monkey dataStructReference.date '/green_binned.bmp'];
filenameReferenceOI=['D:/' dataStructReference.monkey dataStructReference.date '/green_binned.bmp'];
filenameTargetY=[mainPath  dataStructSession.monkey '/' dataStructSession.monkey dataStructSession.date '/green2_binned.bmp'];
filenameTargetOI=['D:/'  dataStructSession.monkey dataStructSession.date '/green2_binned.bmp'];

if isfile(filenameReferenceY)
    imgReference=imread(filenameReferenceY);
else
    imgReference=imread(filenameReferenceOI);
end

if isfile(filenameTargetY)
    imgTarget=imread(filenameTargetY);
elseif isfile(filenameTargetOI)
    imgTarget=imread(filenameTargetOI);
end
%}
%imgReference=imread(['Y:/'  dataStructReference.monkey '/' dataStructReference.monkey dataStructReference.date '/green_binned.bmp']);%imread(convertCharsToStrings(filenameStructRef.greenImg));
%imgTarget=imread(['Y:/'  dataStructSession.monkey '/' dataStructSession.monkey dataStructSession.date '/green2_binned.bmp']);%imread(convertCharsToStrings(filenameStructSession.greenImg));
imgSizeX=[1 size(bitmapRef,1)];
imgSizeY=[1 size(bitmapRef,2)];
% ROI mask
load(filenameStructRef.Orientation,'Mask'); %contains mask, RespCondPCA
Mask=double(Mask);
Mask(Mask==0)=NaN;

%% Get coregistration params for greenImageReference (to be transformed) to greenImageSession (anchor image)
%transform type: auto, similarity
sameSession=isequal(dataStructSession.date,dataStructReference.date);
switch sameSession
    case {1}
        method='auto'; %manual/auto
    case {0}
        method='manual';
end

[imgReferenceCoreg,coregStats,transformParams,imgReference,imgTarget]=coregisterGreenImages(dataStructReference,dataStructSession,Mask);
[~,h]=suplabel('Identifying transformation matrix from reference to current session','t',[0.08 0.08 .84 .80]);
set(h,'FontSize',16)

switch saveFlag
  case {1}
    export_fig(pdfFilename,'-pdf','-nocrop','-append');
end

%% Use transform to align bitmap with greenImageSession, and threshold
VERpca(isnan(VERpca)) = 0;

for i=1:size(bitmapRef,3)
  columnarBitmaps(:,:,i)=double(imwarp(bitmapRef(:,:,i),transformParams,'OutputView',imref2d(size(imgTarget))));
  columnarPCAs(:,:,i)=double(imwarp(VERpca(:,:,i),transformParams,'OutputView',imref2d(size(imgTarget))));
end

%transform and plot aligned images
bitmapRef0 =bitmapRef(:,:,1);
bitmapRef90 =bitmapRef(:,:,2);%floor(size(bitmapRef,3)/2) + 1);

%transform the ref VSD column map
bitmapRef0transformed =  columnarBitmaps(:,:,1);
bitmapRef90transformed = columnarBitmaps(:,:,2);%floor(size(bitmapRef,3)/2) + 1);

%0
figure('name','Coregistered columnar map')
subplot(2,3,1)
imagesc(imfuse(imgTarget,bitmapRef0,'ColorChannels','red-cyan')); axis square
title([sprintf('Pre-coregistration %0g',0) char(176) ],'FontWeight','normal');
addPix2MM(min(imgSizeX),max(imgSizeX),min(imgSizeY),max(imgSizeY),1,2,2);
addPix2MM(min(imgSizeX),max(imgSizeX),min(imgSizeY),max(imgSizeY),1,2,3);

subplot(2,3,2)
imagesc(imfuse(imgTarget,bitmapRef0transformed,'ColorChannels','red-cyan')); axis square
title([sprintf('Post-coregistration %0g',0) char(176)],'FontWeight','normal'); 
addPix2MM(min(imgSizeX),max(imgSizeX),min(imgSizeY),max(imgSizeY),1,2,2);
addPix2MM(min(imgSizeX),max(imgSizeX),min(imgSizeY),max(imgSizeY),1,2,3);

subplot(2,3,3)
imagesc(bitmapRef0transformed); axis square; colormap(gray)
title([sprintf('Bitmap %0g',0) char(176)],'FontWeight','normal');
addPix2MM(min(imgSizeX),max(imgSizeX),min(imgSizeY),max(imgSizeY),1,2,2);

%90
subplot(2,3,4)
imagesc(imfuse(imgTarget,bitmapRef90,'ColorChannels','red-cyan')); axis square
title([sprintf('Pre-coregistration %0g',90) char(176)],'FontWeight','normal');
addPix2MM(min(imgSizeX),max(imgSizeX),min(imgSizeY),max(imgSizeY),1,2,2);
addPix2MM(min(imgSizeX),max(imgSizeX),min(imgSizeY),max(imgSizeY),1,2,3);

subplot(2,3,5)
imagesc(imfuse(imgTarget,bitmapRef90transformed,'ColorChannels','red-cyan')); axis square
title([sprintf('Post-coregistration %0g',90) char(176)],'FontWeight','normal');
addPix2MM(min(imgSizeX),max(imgSizeX),min(imgSizeY),max(imgSizeY),1,2,2);
addPix2MM(min(imgSizeX),max(imgSizeX),min(imgSizeY),max(imgSizeY),1,2,3);

subplot(2,3,6)
imagesc(bitmapRef90transformed); axis square; colormap(gray)
title([sprintf('Bitmap %0g',90) char(176)],'FontWeight','normal');
addPix2MM(min(imgSizeX),max(imgSizeX),min(imgSizeY),max(imgSizeY),1,2,2);
upFontSize(14,.015)

[~,h]=suplabel('Co-registration of columnar map to current session vasculature','t',[.08 .08 .84 .90]);
set(h,'FontSize',16)
switch saveFlag
  case {1}
    export_fig(pdfFilename,'-pdf','-nocrop','-append');
end
%% Plotting
if plotFlag==0
    set(groot,'defaultFigureVisible','on')
end

end