function [diffBitmap,VERpca,columnarBitmapStats,PCAExplTotal]=getColumnarBitmapV4(...
  mainPath,dataStructRef, dataStructSession,bitmapParams, plotFlag, saveFlag)
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

%% Load raw and PCA-ed activity for all orientations, and SNR mask
filenameStructRef=generateFilenames(dataStructRef);
%load(filenameStructSession.TS)
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
Mask(:,400:512)=0; %(420:512)=0;
ROImask=double(Mask); % converts binary columnar mask to double
ROItuningMap=ROImask.*MapAmpOrt; % converts binary columnar mask into graded amplitude mask

%% Figure params
nRows=2;
nCols=6;

%% Processing
[VERpca,normVERpca,histEqCol,gammaCol,columnarBitmap,columnarBitmapStats,optostimBitmap]=generateBitmap(selectedOrt,RespCondPCA,ROImask,bitmapParams);


%% Setup for plot
figure('name','Columnar map')
colormap(gray);
[hAx,~]=tight_subplot(nRows,nCols,[.025 .025]);
plotIdx=[1 7 2 8 3 9 4 10 5 11 6 12];plotCounter=1;
PCAExplTotal=sum(PCAExpl(1:nPCAComp));
selectedOrts=[1 2];
plotCounter=plotter(hAx,VERpca,ROImask,selectedOrts,plotIdx,plotCounter,...
  ['PCA, nComp=' num2str(nPCAComp,'%.0f') ', Expl. var=' num2str(PCAExplTotal,'%.0f') '%'],2,6);
plotCounter=plotter(hAx,normVERpca,ROItuningMap,selectedOrts,plotIdx,plotCounter,'Normalized',2,6);
plotCounter=plotter(hAx,histEqCol,ROItuningMap,selectedOrts,plotIdx,plotCounter,'Equalization',2,6);
plotCounter=plotter(hAx,gammaCol,ROItuningMap,selectedOrts,plotIdx,plotCounter,'Gamma corr.',2,6);
for i=1:2
  plotCounter=plotter(hAx,optostimBitmap(:,:,i),ROItuningMap,selectedOrts(i),plotIdx,plotCounter,...
    ['Adapt. threshold', ' (DC=' sprintf('%.2f',columnarBitmapStats.DC(i)) ', ' sprintf('%.1f%%',columnarBitmapStats.pixDensity(i)) ')'],2,6);
end
for i=1:2
  plotCounter=plotter(hAx,optostimBitmap(:,:,i),ones(512,512),selectedOrts(i),plotIdx,plotCounter,...
    ['Bitmap', ' (DC=' sprintf('%.2f',columnarBitmapStats.DC(i)) ', ' sprintf('%.1f%%',columnarBitmapStats.pixDensity(i)) ')'],2,6);
end
[~,h]=suplabel(['Extracting columnar map (S' dataStructSession.date 'R'  dataStructSession.run ')'],'t',[0.08 0.08 .84 .80]);
set(h,'FontSize',16)
offwarning;

%% Save maps
diffBitmap=optostimBitmap;
columnarBitmap=optostimBitmap;
save(filenameStructSession.colmapMat,'diffBitmap','columnarBitmap','gammaCol','histEqCol','VERpca','normVERpca');
switch saveFlag
  case {1}
    export_fig(filenameStructSession.neurometricPDF,'-pdf','-nocrop');
end
% re-enable plots if turned off
enablePlots

end