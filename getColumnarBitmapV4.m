function bitmapData=getColumnarBitmapV4(currentBlockStruct, imagingData, bitmapData, blockID, ...
  pdfFilename, plotFlag, saveFlag)
disp('Getting column positions...')
%% Change log
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

%% Load raw and PCA-ed activity for all orientations, and SNR mask
%conds
nOrt=numel(imagingData.orts); %excl blanks
selectedOrt=find(ismember(imagingData.orts(:,:,blockID),bitmapData.orts(:,:,blockID))==1);

%ROI
imgDims=size(imagingData.ortpca(:,:,1));
imgX=[0 imgDims(1)];
imgY=[0 imgDims(2)];
ROIx=[1 imgDims(1)];%[round(.016/.016) imgDims(1)];
ROIy=[1 imgDims(2)];%[round(3.2/.016) imgDims(2)];
%ROISquareMask=NaN(imgX(2),imgY(2));
%ROISquareMask(ROIy(1):ROIy(2),ROIx(1):ROIx(2))=1;


%% Extract masked VER and pca-ed VER
%Mask(:,400:512)=0; %(420:512)=0;

%% Figure params
nRows=2;
nCols=6;

%% Processing
[VERpca,normVERpca,histEqCol,gammaCol,columnarBitmap,bitmapData]=generateBitmap(selectedOrt,imagingData,bitmapData,blockID);

%% Setup for plot
figure('name','Columnar map')
colormap(gray);
gap=.025;marginV=.01;marginH=.01;
[hAx,~]=tight_subplot(nRows,nCols,[gap gap], [marginV marginV+.15], [marginH+.05 marginH+.05]);
plotIdx=[1 7 2 8 3 9 4 10 5 11 6 12];plotCounter=1;
PCAExplTotal=sum(imagingData.pcaexpl(1:imagingData.npca,blockID));
selectedOrts=[1 2];
plotCounter=plotter(hAx,VERpca,imagingData.mask(:,:,blockID),selectedOrts,plotIdx,plotCounter,...
  ['PCA, nComp=' num2str(imagingData.npca(:,blockID),'%.0f') ', Expl. var=' num2str(PCAExplTotal,'%.0f') '%'],2,6); addPix2MM(VERpca(:,:,1),7,nRows,nCols);
plotCounter=plotter(hAx,normVERpca,imagingData.ortampmap(:,:,:,blockID),selectedOrts,plotIdx,plotCounter,'Normalized',2,6);
plotCounter=plotter(hAx,histEqCol,imagingData.ortampmap(:,:,:,blockID),selectedOrts,plotIdx,plotCounter,'Equalization',2,6);
plotCounter=plotter(hAx,gammaCol,imagingData.ortampmap(:,:,:,blockID),selectedOrts,plotIdx,plotCounter,'Gamma corr.',2,6);
for i=1:2
  plotCounter=plotter(hAx,bitmapData.columnarbitmap(:,:,i,blockID),imagingData.ortampmap(:,:,:,blockID),selectedOrts(i),plotIdx,plotCounter,...
    'Adaptive thresholding',2,6);
end
for i=1:2
  plotCounter=plotter(hAx,bitmapData.columnarbitmap(:,:,i,blockID),ones(512,512),selectedOrts(i),plotIdx,plotCounter,...
    ['Column positions (' sprintf('%.1f%%',bitmapData.pixDensity(i,blockID)) ')'],2,6);
end
[~,h]=suplabel(['Identifying column positions, ' currentBlockStruct.date 'R'  currentBlockStruct.run],'t',[0.08 0.08 .84 .80]);
set(h,'FontSize',16)
offwarning;

%% Save maps
save(currentBlockStruct.colmapMat,'columnarBitmap','gammaCol','histEqCol','VERpca','normVERpca');
switch saveFlag
  case {1}
    export_fig(pdfFilename,'-pdf','-nocrop');
end
% re-enable plots if turned off
enablePlots

if plotFlag==0
    close
end
end