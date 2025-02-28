function analyzeColumnarSignal(dataStruct, uniqueStr, plotOrNot)
%% Change log
% 1. Suspect that there are pixels with abnormally high value in 90-165deg patterns, try solve by clipping outliers 
% outside +-3sd at [-1 ,1]
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
if plotOrNot==0
    set(groot,'defaultFigureVisible','off')
else
    set(groot,'defaultFigureVisible','on')
end

%% Define column selection method
columnSelectionMethod='MapOrt';

%% Get filenames
filenameStruct=generateFilenames(dataStruct);
pdfFilename=filenameStruct.psychneuroPDF;


%% Load raw and PCA-ed activity for all orientations, and SNR mask
load(filenameStruct.Intg,'DataCond');
load(filenameStruct.Orientation,'Ort','RespCondPCA','Mask','MapOrt','MapAmpOrt','ColorMap'); %contains mask, RespCondPCA
load(filenameStruct.colmapMat,'')
load(filenameStruct.TS);

% ROI mask
ROImask=double(Mask);
ROImaskNaN=ROImask;
ROImaskNaN(ROImaskNaN==0)=NaN;
ROItuningMap=ROImask.*MapAmpOrt;

%% Subtract blank
DataCond=DataCond;
DataCondBlankSubt=DataCond(:,:,2:end)-DataCond(:,:,1);
for condID=1:size(DataCondBlankSubt,3)
    DataCondColumn(:,:,condID)=FilterFermi2D(DataCondBlankSubt(:,:,condID), 0.8, 1.2, 0.016);
end

%% Indexes for blank, congruent, incongruent
blank0=[2 5 8 11];
blank90=[14 17 20 23];
vis0_BCI=[blank0; blank0+1; blank0+2]-1; %-1 remove blank ID
vis90_BCI=[blank90; blank90+2; blank90+1]-1; %-1 remove blank ID
vis090_BCI=[vis0_BCI;vis90_BCI]; 
vis090_BCI=vis090_BCI';
%% Raw
figure
nRows=size(vis090_BCI,1);
nCols=size(vis090_BCI,2);
[hAx,~]=tight_subplot(nRows,nCols,[.05 .05]);
for condID=1:size(DataCondBlankSubt,3)
    condImgID=vis090_BCI(condID);
    axes(hAx(condID))
    imgsc(DataCondBlankSubt(:,:,condImgID).*ROImaskNaN)
    addPix2MM(DataCondBlankSubt(:,:,condImgID).*ROImaskNaN,condID,nRows,nCols);
    caxis([0 3]*10^-1)
end
upFontSize(14,.015)
suplabel('Raw signals','t',[.08 .08 .84 .9]);
export_fig(pdfFilename,'-pdf','-nocrop');

%% Bandpassed columnar (.8-1.2cpmm)
figure
[hAx,~]=tight_subplot(nRows,nCols,[.05 .05]);
for condID=1:size(DataCondColumn,3)
    condImgID=vis090_BCI(condID);
    axes(hAx(condID))
    imgsc(DataCondColumn(:,:,condImgID).*ROImaskNaN)
    addPix2MM(DataCondColumn(:,:,condImgID).*ROImaskNaN,condID,nRows,nCols);
    caxis([-3.5 3.5]*10^-2)
end
upFontSize(14,.015)
suplabel('Columnar signals (0.8-1.2cpmm)','t',[.08 .08 .84 .9]);
export_fig(pdfFilename,'-pdf','-nocrop','-append');

%% Bandpassed columnar (.8-1.2cpmm)
figure
[hAx,~]=tight_subplot(nRows,nCols,[.05 .05]);
for condID=1:size(DataCondColumn,3)
    condImgID=vis090_BCI(condID);
    axes(hAx(condID))
    imgsc(DataCondColumn(:,:,condImgID).*colMask(:,:,1))
    addPix2MM(1,512,1,512,condID,nRows,nCols);
    caxis([-3.5 3.5]*10^-3)
end
suplabel('Masked columnar signals (0.8-1.2cpmm)','t',[.08 .08 .84 .9]);


ROIamp=squeeze(max(DataCondColumn.*ROImaskNaN,[],[1 2]));
normROIamp=rescale(ROIamp);

figure
[hAx,~]=tight_subplot(2,1,[.05 .05]);
axes(hAx(1))
x=1:numel(vis0_BCI);
y=ROIamp(reshape(vis0_BCI',[1 numel(vis0_BCI)]))';
g=[repmat(1,numel(x)/3)]
gscatter(x,y,g,symb)
axes(hAx(2))
scatter(1:numel(vis90_BCI),ROIamp(reshape(vis90_BCI',[1 numel(vis90_BCI)]))')



if plotOrNot==0
    set(groot,'defaultFigureVisible','on')
end
end