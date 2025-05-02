function bitmapData=coregisterBitmap2GreenImgV2(currentBlockStruct,referenceBlockStruct,imagingData,bitmapData,blockID,...
    analysisMode, pdfFilename, plotFlag,saveFlag)
disp('Getting transformation matrix to map reference map to current green image...')
%% Plotting
if plotFlag==0
    set(groot,'defaultFigureVisible','off')
else
    set(groot,'defaultFigureVisible','on')
end

%% Get coregistration params for greenImageReference (to be transformed) to greenImageSession (anchor image)
%transform type: auto, similarity
[imgReferenceCoreg,coregStats,transformParams,imgReference,imgTarget]=coregisterGreenImages(currentBlockStruct,referenceBlockStruct,imagingData.nanmask, analysisMode);
[~,h]=suplabel('Identifying transformation matrix from reference to current session','t',[0.08 0.08 .84 .80]);
set(h,'FontSize',16)
switch saveFlag
  case {1}
    export_fig(pdfFilename,'-pdf','-nocrop','-append');
end

%% Use transform to align bitmap with greenImageSession, and threshold
bitmapData.transformParams(blockID)=transformParams;
for i=1:size(bitmapData.columnarbitmap(:,:,:,blockID),3)
  bitmapData.columnarbitmapCoreg(:,:,i,blockID)=double(imwarp(bitmapData.columnarbitmap(:,:,i,blockID),bitmapData.transformParams(blockID),'OutputView',imref2d(size(imgTarget))));
end

%transform and plot aligned images
bitmapRef0 =bitmapData.columnarbitmap(:,:,1,blockID);
bitmapRef90 =bitmapData.columnarbitmap(:,:,2,blockID);%floor(size(bitmapData.columnarbitmap,3)/2) + 1);

%transform the ref VSD column map
bitmapRef0transformed =  bitmapData.columnarbitmapCoreg(:,:,1,blockID);
bitmapRef90transformed = bitmapData.columnarbitmapCoreg(:,:,2,blockID);%floor(size(bitmapData.columnarbitmap,3)/2) + 1);

%0
figure('name','Coregistered columnar map')
subplot(2,3,1)
imagesc(imfuse(imgTarget,bitmapRef0,'ColorChannels','red-cyan')); axis square
title([sprintf('Pre-coregistration %0g',0) char(176) ],'FontWeight','normal');
addPix2MM(bitmapData.columnarbitmap(:,:,1,blockID),1,2,2);
addPix2MM(bitmapData.columnarbitmap(:,:,1,blockID),1,2,3);

subplot(2,3,2)
imagesc(imfuse(imgTarget,bitmapRef0transformed,'ColorChannels','red-cyan')); axis square
title([sprintf('Post-coregistration %0g',0) char(176)],'FontWeight','normal'); 
addPix2MM(bitmapData.columnarbitmapCoreg(:,:,1,blockID),1,2,2);
addPix2MM(bitmapData.columnarbitmapCoreg(:,:,1,blockID),1,2,3);

subplot(2,3,3)
imagesc(bitmapRef0transformed); axis square; colormap(gray)
title([sprintf('Columns %0g',0) char(176)],'FontWeight','normal');
addPix2MM(bitmapData.columnarbitmapCoreg(:,:,1,blockID),1,2,2);

%90
subplot(2,3,4)
imagesc(imfuse(imgTarget,bitmapRef90,'ColorChannels','red-cyan')); axis square
title([sprintf('Pre-coregistration %0g',90) char(176)],'FontWeight','normal');
addPix2MM(bitmapData.columnarbitmap(:,:,1,blockID),1,2,2);
addPix2MM(bitmapData.columnarbitmap(:,:,1,blockID),1,2,3);

subplot(2,3,5)
imagesc(imfuse(imgTarget,bitmapRef90transformed,'ColorChannels','red-cyan')); axis square
title([sprintf('Post-coregistration %0g',90) char(176)],'FontWeight','normal');
addPix2MM(bitmapData.columnarbitmapCoreg(:,:,1,blockID),1,2,2);
addPix2MM(bitmapData.columnarbitmapCoreg(:,:,1,blockID),1,2,3);

subplot(2,3,6)
imagesc(bitmapRef90transformed); axis square; colormap(gray)
title([sprintf('Columns %0g',90) char(176)],'FontWeight','normal');
addPix2MM(bitmapData.columnarbitmapCoreg(:,:,1,blockID),1,2,2);
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

if plotFlag==0
    close
end
end