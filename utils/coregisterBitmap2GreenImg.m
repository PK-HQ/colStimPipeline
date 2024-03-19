function [columnarBitmaps, similarityScore]=coregisterBitmap2GreenImg(dataStructReference,dataStructSession,bitmapRef, plotOrNot)
projParam=[-3 +10 6 .95];
thresholdVal=90;

%% Plotting
if plotOrNot==0
    set(groot,'defaultFigureVisible','off')
else
    set(groot,'defaultFigureVisible','on')
end

%% Get filenames
filenameStructSession=generateFilenames(dataStructSession);
filenameStructRef=generateFilenames(dataStructReference);
pdfFilename='Y:\Chip\Chip20221102\coregGreen.pdf'%[filenameStructSession.plotPath 'M' dataStructSession.monkeyNo 'R' dataStructReference.date 'S' dataStructSession.date '.pdf'];

%% Get coregistration params for greenImageReference (to be transformed) to greenImageSession (anchor image)
% load images
imgReference=imread(['Y:\'  dataStructReference.monkey '\' dataStructReference.monkey dataStructReference.date '\green_binned.bmp']);%imread(convertCharsToStrings(filenameStructRef.greenImg));
imgTarget=imread(['Y:\'  dataStructSession.monkey '\' dataStructSession.monkey dataStructSession.date '\green2_binned.bmp']);%imread(convertCharsToStrings(filenameStructSession.greenImg));

imgSizeX=[1 size(bitmapRef,1)];
imgSizeY=[1 size(bitmapRef,2)];
%{
figure()
colormap(gray)
subplot(2,4,1)
imagesc(imgTarget); axis square
addPix2MM(min(imgSizeX),max(imgSizeX),min(imgSizeY),max(imgSizeY),1,2,2);
title('Current vasculature')
subplot(2,4,2)
imagesc(imgReference); axis square
addPix2MM(min(imgSizeX),max(imgSizeX),min(imgSizeY),max(imgSizeY),2,2,2);
title('Orientation block vasculature')
% highpass filter to emphasize vasculature
subplot(2,4,5)
imgTarget = FilterFermi2D(imgTarget,0.01,inf,1);
imagesc(imgTarget); axis square
title('Current vasculature (high-pass)')
addPix2MM(min(imgSizeX),max(imgSizeX),min(imgSizeY),max(imgSizeY),3,2,2);
subplot(2,4,6)
imgReference = FilterFermi2D(imgReference,0.01,inf,1);
imagesc(imgReference); axis square
title('Orientation block vasculature (high-pass)')
addPix2MM(min(imgSizeX),max(imgSizeX),min(imgSizeY),max(imgSizeY),4,2,2);
%export_fig(pdfFilename,'-pdf','-nocrop');

%display
suplabel('Coregistration of reference and session green images','t');
subplot(2,4,[3 4 7 8])
imagesc(imfuse(imgTarget,imgReference,'ColorChannels','red-cyan')); axis square
addPix2MM(min(imgSizeX),max(imgSizeX),min(imgSizeY),max(imgSizeY),1,2,2);
title('Original')
%}

%transform type: auto, similarity
method='Manual'; %or 'Manual'
[imgReferenceCoreg,coregStats,transformParams]=coregisterGreenImages(imgReference,imgTarget,method,filenameStructSession);
figure
imagesc(imfuse(imgTarget,imgReferenceCoreg,'ColorChannels','red-cyan')); axis square
title(sprintf([method ', r=%.2g (T-x=%.2g, T-y=%.2g, Rot=%.2g, Scale=%.2g)'],...
    coregStats.Similarity,coregStats.Translate(1),coregStats.Translate(2),coregStats.Rotation,coregStats.Rotation));
addPix2MM(min(imgSizeX),max(imgSizeX),min(imgSizeY),max(imgSizeY),3,2,2);
1;
similarityScore=coregStats.Similarity;
%transform type: manual, landmark
%{
figure
method='manual';
[imgReferenceCoreg,coregStats,transformParams]=coregisterGreenImages(imgReference,imgTarget,method,[],filenameStructRef,filenameStructSession);
imagesc(imfuse(imgTarget,imgReferenceCoreg,'ColorChannels','red-cyan')); axis square
title(sprintf('Manual, r=%.2g (T-x=%.2g, T-y=%.2g, Rot=%.2g, Scale=%.2g)',...
    coregStats.Similarity,coregStats.Translate(1),coregStats.Translate(2),coregStats.Rotation));
addPix2MM(min(imgSizeX),max(imgSizeX),min(imgSizeY),max(imgSizeY),3,2,2);
%export_fig(pdfFilename,'-pdf','-nocrop','-append');
%}


%transform type: affine
%{
transformType = 'affine';
[imgReferenceCoreg,coregStats,~]=coregisterGreenImages(imgReference,imgTarget,transformType,transformParams,maxIters);
subplot(2,2,4)
imagesc(imfuse(imgTarget,imgReferenceCoreg,'ColorChannels','red-cyan')); axis square
title(sprintf('Affine-transform (R=%.2g, T-x=%.2g, T-y=%.2g, Rot=%.2g, Scaling=%.2g)',...
    coregStats.Similarity,coregStats.Translate(1),coregStats.Translate(2),coregStats.Rotation,coregStats.Scale));
addPix2MM(min(imgSizeX),max(imgSizeX),min(imgSizeY),max(imgSizeY),4,2,2);
%}
%% Use transform to align bitmap with greenImageSession, and threshold
for i=1:size(bitmapRef,3)
	columnarBitmaps(:,:,i)=double(imwarp(bitmapRef(:,:,i),transformParams,'OutputView',imref2d(size(imgTarget))));
    tmp=columnarBitmaps(:,:,i);
    x=columnarBitmaps(:,:,i)>prctile(tmp(:),thresholdVal);
    columnarBitmaps(:,:,i) =double(x);
end

%transform and plot aligned images
bitmapRef0 =bitmapRef(:,:,1);
bitmapRef90 =bitmapRef(:,:,floor(size(bitmapRef,3)/2) + 1);

%transform the ref VSD column map
bitmapRef0transformed =  columnarBitmaps(:,:,1);
bitmapRef90transformed = columnarBitmaps(:,:,floor(size(bitmapRef,3)/2) + 1);

%0
figure()
subplot(2,3,1)
imagesc(imfuse(imgTarget,bitmapRef0,'ColorChannels','red-cyan')); axis square
title([sprintf('Pre-coregistration (%0g',0) '\circ)'],'FontWeight','normal');
addPix2MM(min(imgSizeX),max(imgSizeX),min(imgSizeY),max(imgSizeY),1,2,2);
addPix2MM(min(imgSizeX),max(imgSizeX),min(imgSizeY),max(imgSizeY),1,2,3);

subplot(2,3,2)
imagesc(imfuse(imgTarget,bitmapRef0transformed,'ColorChannels','red-cyan')); axis square
title([sprintf('Post-coregistration (%0g',0) '\circ)'],'FontWeight','normal');
addPix2MM(min(imgSizeX),max(imgSizeX),min(imgSizeY),max(imgSizeY),1,2,2);
addPix2MM(min(imgSizeX),max(imgSizeX),min(imgSizeY),max(imgSizeY),1,2,3);

subplot(2,3,3)
imagesc(bitmapRef0transformed); axis square; colormap(gray)
title([sprintf('Bitmap %0g',0) '\circ)'],'FontWeight','normal');
addPix2MM(min(imgSizeX),max(imgSizeX),min(imgSizeY),max(imgSizeY),1,2,2);

%90
subplot(2,3,4)
imagesc(imfuse(imgTarget,bitmapRef90,'ColorChannels','red-cyan')); axis square
title([sprintf('Pre-coregistration (%0g',90) '\circ)'],'FontWeight','normal');
addPix2MM(min(imgSizeX),max(imgSizeX),min(imgSizeY),max(imgSizeY),1,2,2);
addPix2MM(min(imgSizeX),max(imgSizeX),min(imgSizeY),max(imgSizeY),1,2,3);

subplot(2,3,5)
imagesc(imfuse(imgTarget,bitmapRef90transformed,'ColorChannels','red-cyan')); axis square
title([sprintf('Post-coregistration (%0g',90) '\circ)'],'FontWeight','normal');
addPix2MM(min(imgSizeX),max(imgSizeX),min(imgSizeY),max(imgSizeY),1,2,2);
addPix2MM(min(imgSizeX),max(imgSizeX),min(imgSizeY),max(imgSizeY),1,2,3);

subplot(2,3,6)
imagesc(bitmapRef90transformed); axis square; colormap(gray)
title([sprintf('Bitmap %0g',90) '\circ)'],'FontWeight','normal');
addPix2MM(min(imgSizeX),max(imgSizeX),min(imgSizeY),max(imgSizeY),1,2,2);

%% Plotting
if plotOrNot==0
    set(groot,'defaultFigureVisible','on')
end

end