function [diffBitmap,VERpca,PCAExplTotal,pixelDensities,bitmapEnergies]=getColumnarBitmapV4(...
  mainPath,dataStructRef, dataStructSession,bitmapParams, plotFlag)
%% Change log
% 2. F0 footprint crop by using defining ROI per orientation map that is > 70% amplitude, then crop all to that ROI

%% What it does
% Loads raw and PCA-ed responses and mask. 
% Collects SNR and ROI masked blank subtracted VER and PCA-ed VER (VER,
% VERpca)
% 1. SNR and ROI masking 
% 2. Rescale + histogram equalization, i.e. histEqCol
% 3. Gamma correction, i.e. gammaCol
% 4. Thresholding, i.e. columnarBitmap

%% Plotting
disablePlots(plotFlag)

%% Define column selection method
orangeLightPD=56.1;

%% Get filenames
filenameStructSession=generateFilenames(dataStructSession);
load(filenameStructSession.TS)

%% Load raw and PCA-ed activity for all orientations, and SNR mask
filenameStructRef=generateFilenames(dataStructRef);
load(filenameStructRef.FFTAmp,'DataCond');
load(filenameStructRef.Orientation,'Ort','RespCondPCA','Mask','MapAmpOrt','nPCAComp','PCAExpl'); %contains mask, RespCondPCA

%conds
nOrt=numel(Ort); %excl blanks
selectedOrt=find(ismember(Ort,bitmapParams.desiredOrts)==1);

%ROI
imgDims=size(RespCondPCA(:,:,1));
imgX=[0 imgDims(1)];
imgY=[0 imgDims(2)];
ROIx=[1 imgDims(1)];%[round(.016/.016) imgDims(1)];
ROIy=[1 imgDims(2)];%[round(3.2/.016) imgDims(2)];
%ROISquareMask=NaN(imgX(2),imgY(2));
%ROISquareMask(ROIy(1):ROIy(2),ROIx(1):ROIx(2))=1;

%n
nBlanks=size(DataCond,3)-nOrt;
DataCondBlankSubt=DataCond(:,:,nBlanks+1:nOrt+nBlanks)-mean(DataCond(:,:,1:nBlanks),3);

%% Extract masked VER and pca-ed VER
%ROI mask
ROImask=double(Mask); % converts binary columnar mask to double
ROItuningMap=ROImask.*MapAmpOrt; % converts binary columnar mask into graded amplitude mask

%% Figure params
nRows=2;
nCols=6;

%% Processing
[VERpca,normVERpca,histEqCol,gammaCol,columnarBitmap,columnarBitmapStats,optostimBitmap]=generateBitmap(selectedOrt,RespCondPCA,ROImask,bitmapParams);















%% Setup for plot
figure('name','PCA-ed responses')
colormap(gray);
[hAx,~]=tight_subplot(nRows,nCols,[.025 .025]);
plotIdx=[1, ((nRows*nCols)/2 +1)];plotCounter=1;
PCAExplTotal=sum(PCAExpl(1:nPCAComp));

plotCounter=plotter(VERpca,ROItuningMap,selectedOrts,plotIdx,plotCounter,['PCA, nComp=' num2str(nPCAComp,'%.0f') ', Expl. var=' num2str(PCAExplTotal,'%.0f') '%']);
plotCounter=plotter(normVERpca,ROItuningMap,selectedOrts,plotIdx,plotCounter,'Normalized');
plotCounter=plotter(histEqCol,ROItuningMap,selectedOrts,plotIdx,plotCounter,'Equalization');
plotCounter=plotter(gammaCol,ROItuningMap,selectedOrts,plotIdx,plotCounter,'Gamma corr.');
for i=1:2
  plotCounter=plotter(optostimBitmap(:,:,i),ROItuningMap,selectedOrts(i),plotIdx,plotCounter,['Adapt. threshold', ' (DC=' sprintf('%.2f',estDutycycle(i)) ', ' sprintf('%.1f%%',pixelDensity(i)) ')']);
end
1;


%% 1 Visualize PCA-ed data


%% Plotting
for ortNo=selectedOrt
  axes(hAx(plotIdx(plotCounter)));
  img=VERpca(:,:,ortNo);
  greyNaN(img)
  
  titleL1=[sprintf('%0g',Ort(ortNo)) '\circ'];
 
  if ortNo==1
      title({['PCA-ed (nComp=' num2str(nPCAComp,'%.0f') ', ' 'Expl. var=' num2str(PCAExplTotal,'%.0f') '%'],titleL1},'FontWeight','normal');
  else
      xlabel('X (mm)');
      ylabel('Y (mm)');
  end
end
  





%%
for ortNo=selectedOrt
    axes(hAx(plotIdx(plotCounter)));
    greyNaN(VERpca(:,:,ortNo));
    
    titleL1=[sprintf('%0g',Ort(ortNo)) '\circ'];
    if ortNo==1
        title({['PCA-ed (nComp=' num2str(nPCAComp,'%.0f') ')'],['Expl. var=' num2str(PCAExplTotal,'%.0f') '%'],titleL1},'FontWeight','normal');
    else
        xlabel('X (mm)');
        ylabel('Y (mm)');

    end
    axis square;
    addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),ortNo,nRows,nCols);
    colormap(gray);colorbar;caxis([-2 2]*10^-3)
    plotCounter=plotCounter+1;
end
upFontSize(14,.015)
[~,h]=suplabel('Columnar bitmap generation','t',[.08 .08 .84 .80]);
set(h,'FontSize',16)

%% 2 Normalize it
colormap(gray);
plotIdx=[1, ((nRows*nCols)/2 +1)]+1;plotCounter=1;
for ortNo=selectedOrt
    axes(hAx(plotIdx(plotCounter)));
    %rescale to norm
    map=VERpca(:,:,ortNo);
    normMap=rescale(map,0,1); 
    %clips image, remap to -1 to 1
    normMap = normMap .* ROImaskNaN;
    normVERpca(:,:,ortNo)=normMap;

    %plot
    greyNaN(normVERpca(:,:,ortNo).*ROItuningMap);
    if plotCounter==1
        title({'Normalized'},'FontWeight','normal');
    else
        title({titleL1},'FontWeight','normal');
    end
    axis square;
    addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),ortNo,nRows,nCols);
    colormap(gray);colorbar;
    plotCounter=plotCounter+1;
end

%% 3 Histogram equalization to maximize dynamic range across image
plotIdx=[1, ((nRows*nCols)/2 +1)]+2;plotCounter=1;
for ortNo=selectedOrt
    axes(hAx(plotIdx(plotCounter)));
    histEqCol(:,:,ortNo)=adapthisteq(normVERpca(:,:,ortNo),'NumTiles',imgDims./bitmapParams.gridSize,'Range','full');
    greyNaN(histEqCol(:,:,ortNo).*ROItuningMap);colorbar;colormap(fireice)
    if plotCounter==1
        title({'Equalization'},'FontWeight','normal');
    else
        title({[sprintf('%0g',Ort(ortNo)) '\circ']},'FontWeight','normal');
    end
    axis square;
    addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),ortNo,nRows,nCols);
    colormap(gray);colorbar;caxis('auto')
    plotCounter=plotCounter+1;
    offwarning
end
upFontSize(14,.015)
set(h,'FontSize',16)

%% 4 Gamma correction
plotIdx=[1, ((nRows*nCols)/2 +1)]+3;plotCounter=1;
for ortNo=selectedOrt
    axes(hAx(plotIdx(plotCounter)));
    gammaCol(:,:,ortNo)=imadjust(histEqCol(:,:,ortNo),[],[],bitmapParams.gammaCorr);
    greyNaN(gammaCol(:,:,ortNo).*ROItuningMap);
    if plotCounter==1
        title({'Gamma corr.'},'FontWeight','normal');
    else
        title({[sprintf('%0g',Ort(ortNo)) '\circ']},'FontWeight','normal');
    end
    axis square;
    addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),ortNo,nRows,nCols);
    colormap(gray);colorbar;caxis('auto')   
    plotCounter=plotCounter+1;
end
upFontSize(14,.015)
set(h,'FontSize',16)

%% 5 Adaptive (local) threshold of the gamma map to obtain bitmap
plotIdx=[1, ((nRows*nCols)/2 +1)]+4;plotCounter=1;
for ortNo=selectedOrt
    axes(hAx(plotIdx(plotCounter)));
    % do adaptive thresholding
    threshold = adaptthresh(double(gammaCol(:,:,ortNo)), bitmapParams.sensitivity,...
        'NeighborhoodSize',2*floor(size(gammaCol(:,:,ortNo))/16)+1);
    columnarBitmap(:,:,ortNo) = double(imbinarize(double(gammaCol(:,:,ortNo)),threshold)).*ROImask; %zero mask

    %calculate pixel density
    bitmapNaN=columnarBitmap(:,:,ortNo).*ROImaskNaN;
    pixelDensity=nansum(bitmapNaN(:))*100/nansum(ROImaskNaN(:));
    equivalentDutyCycle=sqrt(pixelDensity*2/100)*5/10;
    %plot
    greyNaN(columnarBitmap(:,:,ortNo).*ROImask.*ROItuningMap,.7);
    titleL1=['DC=' sprintf('%.2f',equivalentDutyCycle) ' (' sprintf('%.1f%%',pixelDensity) ')'];
    if plotCounter==1
        title({'Adapt. threshold',titleL1},'FontWeight','normal');
    else
        title({titleL1},'FontWeight','normal');
    end
    axis square;
    addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),ortNo,nRows,nCols);
    colormap(gray);colorbar;caxis([0 1]);
    plotCounter=plotCounter+1;
end
upFontSize(14,.015)
set(h,'FontSize',16)

%% 6 Difference bitmaps
plotIdx=[1, ((nRows*nCols)/2 +1)]+5;plotCounter=1;
for ortNo=selectedOrt(plotCounter)
    axes(hAx(plotIdx(plotCounter)));
    % do adaptive thresholding
    diffBitmap(:,:,ortNo)=double(columnarBitmap(:,:,1)>columnarBitmap(:,:,7));
    %calculate pixel density
    diffbitmapNaN=diffBitmap(:,:,ortNo).*ROImaskNaN;
    pixelDensity=nansum(diffbitmapNaN(:))*100/nansum(ROImaskNaN(:));
    equivalentDutyCycle=sqrt(pixelDensity*2/100)*5/10;

    if isSummary
      %calculate bitmap energy
      bitmapEnergy=calculateBitmapEnergy(TS, orangeLightPD, pixelDensity);
      pixelDensities(ortNo)=pixelDensity;
      bitmapEnergies(ortNo)=bitmapEnergy;
    else
      pixelDensities(ortNo)=pixelDensity;
      bitmapEnergies(ortNo)=NaN;
    end

    %plot
    greyNaN(diffBitmap(:,:,ortNo).*ROImask.*ROItuningMap,.7);
    titleL1=['DC=' sprintf('%.2f',equivalentDutyCycle) ' (' sprintf('%.1f%%',pixelDensity) ')'];
    if plotCounter==1
        title({'Diff. bitmap',titleL1},'FontWeight','normal');
    else
        title({titleL1},'FontWeight','normal');
    end
    axis square;
    addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),ortNo,nRows,nCols);
    colormap(gray);colorbar;caxis([0 1]);
    plotCounter=plotCounter+1;
end

for ortNo=selectedOrt(plotCounter)
    axes(hAx(plotIdx(plotCounter)));
    % do adaptive thresholding
    diffBitmap(:,:,ortNo)=double(columnarBitmap(:,:,7)>columnarBitmap(:,:,1));
    %calculate pixel density
    diffbitmapNaN=diffBitmap(:,:,ortNo).*ROImaskNaN;
    pixelDensity=nansum(diffbitmapNaN(:))*100/nansum(ROImaskNaN(:));
    equivalentDutyCycle=sqrt(pixelDensity*2/100)*5/10;

    if isSummary
      %calculate bitmap energy
      bitmapEnergy=calculateBitmapEnergy(TS, orangeLightPD, pixelDensity);
      pixelDensities(ortNo)=pixelDensity;
      bitmapEnergies(ortNo)=bitmapEnergy;
    else
      pixelDensities(ortNo)=pixelDensity;
      bitmapEnergies(ortNo)=NaN;
    end

    %plot
    greyNaN(diffBitmap(:,:,ortNo).*ROItuningMap,.7);
    titleL1=['DC=' sprintf('%.2f',equivalentDutyCycle) ' (' sprintf('%.1f%%',pixelDensity) ')'];
    if plotCounter==1
        title({'Diff. bitmap',titleL1},'FontWeight','normal');
    else
        title({titleL1},'FontWeight','normal');
    end
    axis square;
    addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),ortNo,nRows,nCols);
    colormap(gray);colorbar;caxis([0 1]);
    plotCounter=plotCounter+1;
end
upFontSize(14,.015)
set(h,'FontSize',16)

export_fig(filenameStructSession.neurometricPDF,'-pdf','-nocrop');

diffBitmap(:,:,selectedOrt);
bitmapEnergies=bitmapEnergies(selectedOrt);


%% Save maps
save(filenameStructSession.colmapMat,'diffBitmap','columnarBitmap','gammaCol','histEqCol','VERpca','normVERpca');

% re-enable plots if turned off
enablePlots

end