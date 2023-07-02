function bitmapSession=coregisterBitmap(dataStructReference,dataStructSession,bitmapTar,bitmapRef)
useROI=0;
projParam=[-3 +10 6 .95];
thresholdVal=90;
%% Get filenames
filenameStructSession=generateFilenames(dataStructSession);
filenameStructRef=generateFilenames(dataStructReference);
pdfFilename=[filenameStructSession.plotPath 'M' dataStructSession.monkeyNo 'R' dataStructReference.date 'S' dataStructSession.date '.pdf'];

%% Get coregistration params for greenImageReference (to be transformed) to greenImageSession (anchor image)
% load images
%load(filenameStructRef.initTransformParams)
imgTarget=imread('Y:\Chip\Chip20221017\greenR_binned.bmp');%imread(convertCharsToStrings(filenameStructSession.greenImg));
imgReference=imread('Y:\Chip\Chip20221017\green_binned.bmp');%imread(convertCharsToStrings(filenameStructRef.greenImg));
%imgTarget=rgb2gray(imgTarget);
%imgReference=rgb2gray(imgReference);

%ROI
switch useROI
    case {1}
        frameExpansion=75;
        offset=-7;
        ROIx=150-frameExpansion:385+frameExpansion;
        ROIy=200-frameExpansion-offset:430+frameExpansion-offset;
        imgTarget=imgTarget(ROIy,ROIx);
        bitmapTar=bitmapTar(ROIy,ROIx,:);
        ROIx=150-frameExpansion:385+frameExpansion;
        ROIy=200-frameExpansion:430+frameExpansion;
        imgReference=imgReference(ROIy,ROIx);
        bitmapRef=bitmapRef(ROIy,ROIx,:);
    case {0}
        ROIx=1:512;%150:385;
        ROIy=1:512;%200:430;
        %imgTarget=imgTarget(ROIy,ROIx);
        %bitmapTar=bitmapTar(ROIy,ROIx,:);
        ROIx=1:512;%150-50:385+50;
        ROIy=1:512;%200-80:430;
        %imgReference=imgReference(ROIy,ROIx);
        %bitmapRef=bitmapRef(ROIy,ROIx,:);
end



figure()
colormap(gray)
subplot(2,2,1)
imagesc(imgTarget); axis square
addPix2MM(min(ROIx),max(ROIx),min(ROIy),max(ROIy),1,2,2);
title('GCaMP target image')
subplot(2,2,2)
imagesc(imgReference); axis square
addPix2MM(min(ROIx),max(ROIx),min(ROIy),max(ROIy),2,2,2);
title('VSD reference image')
% highpass filter to emphasize vasculature
subplot(2,2,3)
imgTarget = FilterFermi2D(imgTarget,0.01,inf,1);
imagesc(imgTarget); axis square
title('GCaMP target image, high-pass')
addPix2MM(min(ROIx),max(ROIx),min(ROIy),max(ROIy),3,2,2);
subplot(2,2,4)
imgReference = FilterFermi2D(imgReference,0.01,inf,1);
imagesc(imgReference); axis square
title('VSD reference image, high-pass')
addPix2MM(min(ROIx),max(ROIx),min(ROIy),max(ROIy),4,2,2);
%export_fig(pdfFilename,'-pdf','-nocrop');

maxIters=300;
%display
figure()
suplabel('Coregistration of reference and session green images','t');
subplot(1,2,1)
imagesc(imfuse(imgTarget,imgReference,'ColorChannels','red-cyan')); axis square
addPix2MM(min(ROIx),max(ROIx),min(ROIy),max(ROIy),1,2,2);
title('Original')

%transform type: manual, landmark
figure
method='manual';
[imgReferenceCoreg,coregStats,transformParams]=coregisterGreenImages(imgReference,imgTarget,method,[],filenameStructRef,filenameStructSession);
imagesc(imfuse(imgTarget,imgReferenceCoreg,'ColorChannels','red-cyan')); axis square
title(sprintf('Manual, r=%.2g (T-x=%.2g, T-y=%.2g, Rot=%.2g, Scale=%.2g)',...
    coregStats.Similarity,coregStats.Translate(1),coregStats.Translate(2),coregStats.Rotation));
addPix2MM(min(ROIx),max(ROIx),min(ROIy),max(ROIy),3,2,2);
%export_fig(pdfFilename,'-pdf','-nocrop','-append');

%transform type: auto, similarity
%{
method='auto';
[imgReferenceCoreg,coregStats,transformParams]=coregisterGreenImages(imgReference,imgTarget,method,[]);
subplot(1,3,3)
imagesc(imfuse(imgTarget,imgReferenceCoreg,'ColorChannels','red-cyan')); axis square
title(sprintf('Manual, r=%.2g (T-x=%.2g, T-y=%.2g, Rot=%.2g, Scale=%.2g)',...
    coregStats.Similarity,coregStats.Translate(1),coregStats.Translate(2),coregStats.Rotation,coregStats.Scale));
addPix2MM(min(ROIx),max(ROIx),min(ROIy),max(ROIy),2,2,2);
1;
%}
%transform type: affine
%{
transformType = 'affine';
[imgReferenceCoreg,coregStats,~]=coregisterGreenImages(imgReference,imgTarget,transformType,transformParams,maxIters);
subplot(2,2,4)
imagesc(imfuse(imgTarget,imgReferenceCoreg,'ColorChannels','red-cyan')); axis square
title(sprintf('Affine-transform (R=%.2g, T-x=%.2g, T-y=%.2g, Rot=%.2g, Scaling=%.2g)',...
    coregStats.Similarity,coregStats.Translate(1),coregStats.Translate(2),coregStats.Rotation,coregStats.Scale));
addPix2MM(min(ROIx),max(ROIx),min(ROIy),max(ROIy),4,2,2);
%}
%% Use transform to align bitmap with greenImageSession
%transform and plot aligned images
bitmapTar0nan =bitmapTar(:,:,1);
bitmapTar90nan =bitmapTar(:,:,7);
bitmapRef0nan =bitmapRef(:,:,1);
bitmapRef90nan =bitmapRef(:,:,7);

bitmapTar0 = double(bitmapTar(:,:,1));%imwarp(,transformParams,'OutputView',imref2d(size(imgTarget)));
bitmapTar90 = double(bitmapTar(:,:,7));%imwarp(,transformParams,'OutputView',imref2d(size(imgTarget)));
bitmapTar0(isnan(bitmapTar0)) = 0;
bitmapTar90(isnan(bitmapTar90)) = 0;

bitmapRef0 = double(bitmapRef(:,:,1));%imwarp(,transformParams,'OutputView',imref2d(size(imgTarget)));
bitmapRef90 = double(bitmapRef(:,:,7));%imwarp(,transformParams,'OutputView',imref2d(size(imgTarget)));
bitmapRef0(isnan(bitmapRef0)) = 0;
bitmapRef90(isnan(bitmapRef90)) = 0;

%transform the ref VSD column map
bitmapRef0 =  double(imwarp(bitmapRef0,transformParams,'OutputView',imref2d(size(imgTarget))));
bitmapRef90 = double(imwarp(bitmapRef90,transformParams,'OutputView',imref2d(size(imgTarget))));

%0.1DC = 15x15px, 0.2DC=30x30px

%0
figure()
bitmapRef0nan=double(bitmapRef0nan);
bitmapRef90nan=double(bitmapRef90nan);
tMask = bitmapTar0>0;
subplot(2,4,1)
imagesc(imfuse(imgTarget,bitmapRef0.*tMask,'ColorChannels','red-cyan')); axis square
title([sprintf('Session + Reference columns (%0g',0) '\circ)'],'FontWeight','normal');
addPix2MM(min(ROIx),max(ROIx),min(ROIy),max(ROIy),1,2,2);
addPix2MM(min(ROIx),max(ROIx),min(ROIy),max(ROIy),1,2,3);

subplot(2,4,2)
imagesc(imfuse(imgTarget,bitmapTar0.*tMask,'ColorChannels','red-cyan')); axis square
title([sprintf('%0g',0) '\circ'],'FontWeight','normal');
title([sprintf('Session + Session columns (%0g',0) '\circ)'],'FontWeight','normal');
addPix2MM(min(ROIx),max(ROIx),min(ROIy),max(ROIy),2,2,3);

subplot(2,4,3)
bitmapTar0nan.*tMask;
x=ProjectImage(bitmapRef0nan,projParam,size(bitmapRef0nan));
imagesc(imfuse(mat2gray(bitmapTar0nan.*tMask),mat2gray(x),'ColorChannels','red-cyan','Scaling','independent')); axis square

imagesc(imfuse(mat2gray(bitmapTar0nan.*tMask),mat2gray(bitmapRef0nan.*tMask),'ColorChannels','red-cyan','Scaling','independent')); axis square
tCorr = corrcoef(bitmapTar0nan.*tMask,bitmapRef0nan.*tMask);
title([sprintf('Reference columns + Session columns (R=%.1g, %0g',tCorr(1,2),0) '\circ)'],'FontWeight','normal');
addPix2MM(min(ROIx),max(ROIx),min(ROIy),max(ROIy),3,2,3);

subplot(2,4,4)
greyNaN((bitmapTar0nan>prctile(bitmapTar0nan(:),thresholdVal))); axis square
title([sprintf('Bitmap (%0g',0) '\circ)'],'FontWeight','normal');
addPix2MM(min(ROIx),max(ROIx),min(ROIy),max(ROIy),3,2,3);


%90
tMask = double(bitmapTar90>0);
subplot(2,4,5)
imagesc(imfuse(imgTarget,bitmapRef90.*tMask,'ColorChannels','red-cyan')); axis square
title([sprintf('Session + Reference columns (%0g',90) '\circ)'],'FontWeight','normal');
addPix2MM(min(ROIx),max(ROIx),min(ROIy),max(ROIy),4,2,3);

subplot(2,4,6)
imagesc(imfuse(imgTarget,bitmapTar90.*tMask,'ColorChannels','red-cyan')); axis square
title([sprintf('%0g',0) '\circ'],'FontWeight','normal');
title([sprintf('Session + Session columns (%0g',90) '\circ)'],'FontWeight','normal');
addPix2MM(min(ROIx),max(ROIx),min(ROIy),max(ROIy),5,2,3);



subplot(2,4,7)
x=ProjectImage(bitmapRef90nan,projParam,size(bitmapRef90nan));
imagesc(imfuse(mat2gray(bitmapTar90nan.*tMask),mat2gray(x),'ColorChannels','red-cyan','Scaling','independent')); axis square

subplot(2,4,7)
imagesc(imfuse(mat2gray(bitmapTar90nan.*tMask),mat2gray(bitmapRef90nan.*tMask),'ColorChannels','red-cyan','Scaling','independent')); axis square
tCorr = corrcoef(bitmapTar90nan.*tMask,bitmapRef90nan.*tMask);
title([sprintf('Reference columns + Session columns (R=%.1g, %0g',tCorr(1,2),90) '\circ)'],'FontWeight','normal');
addPix2MM(min(ROIx),max(ROIx),min(ROIy),max(ROIy),3,2,3);

subplot(2,4,8)
greyNaN((bitmapTar90nan>prctile(bitmapTar90nan(:),thresholdVal))); axis square
title([sprintf('Bitmap %0g',90) '\circ)'],'FontWeight','normal');
addPix2MM(min(ROIx),max(ROIx),min(ROIy),max(ROIy),3,2,3);

%bitmapSession(:,:,1)=(bitmapTar0nan>prctile(bitmapTar0nan(:),thresholdVal));
%bitmapSession(:,:,2)=(bitmapTar90nan>prctile(bitmapTar90nan(:),thresholdVal));

for i=1:size(bitmapTar,3)
    collapsed=double(bitmapTar(:,:,i));
    x=bitmapTar(:,:,i)>prctile(collapsed(:),thresholdVal);
    bitmapSession(:,:,i) =double(x);
end
%export_fig(pdfFilename,'-pdf','-nocrop','-append');

%% Account for PRF

%% Adapt to projector center (translate + rotate + black bars)

end