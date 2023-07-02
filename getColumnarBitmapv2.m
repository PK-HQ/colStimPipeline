function [columnarBitmap,columnarMaps, colVERraw,colVERpca]=getColumnarBitmap(dataStruct,gridSize,gammaCorr,threshPrctile,uniqueStr, plotOrNot)
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
pdfFilename=[filenameStruct.plotPath 'M' dataStruct(1).monkeyNo uniqueStr dataStruct(1).date '.pdf'];

%% Load raw and PCA-ed activity for all orientations, and SNR mask
load(filenameStruct.FFTAmp,'DataCond');
load(filenameStruct.Orientation,'Ort','RespCondPCA','Mask','MapOrt'); %contains mask, RespCondPCA
%conds
nOrt=numel(Ort); %excl blanks
%ROI
imgDims=size(RespCondPCA(:,:,1));
imgX=[0 imgDims(1)];
imgY=[0 imgDims(2)];
ROIx=[1 imgDims(1)];%[round(.016/.016) imgDims(1)];
ROIy=[1 imgDims(2)];%[round(3.2/.016) imgDims(2)];
ROIMask=NaN(imgX(2),imgY(2));
ROIMask(ROIy(1):ROIy(2),ROIx(1):ROIx(2))=1;
%n
nBlanks=size(DataCond,3)-nOrt;
DataCondBlankSubt=DataCond(:,:,nBlanks+1:nOrt+nBlanks)-mean(DataCond(:,:,1:nBlanks),3);

%% Extract masked VER and pca-ed VER
%figure('units','normalized','outerposition',[0 0 1 1])

nRows=2;
nCols=6;
%ROI mask
Mask=double(Mask);
Mask(Mask==0)=NaN;
for ortNo=1:nOrt
    %VER for orientation, masked with SNR mask
    VERraw(:,:,ortNo)=DataCondBlankSubt(:,:,ortNo).*Mask;
    VERpca(:,:,ortNo)=RespCondPCA(:,:,ortNo).*Mask;
    
    %VER for orientation, masked with ROI mask (to isolate site)
    VERraw(:,:,ortNo)=VERraw(:,:,ortNo).*ROIMask;
    VERpca(:,:,ortNo)=VERpca(:,:,ortNo).*ROIMask;

    colVERraw(:,:,ortNo)=VERraw(:,:,ortNo).*ROIMask;
    colVERpca(:,:,ortNo)=VERpca(:,:,ortNo).*ROIMask;
    colPositions(:,:,ortNo)=colVERpca(:,:,ortNo)>0;
    

end

%% PCA-ed data
figure('name','PCA-ed responses')
colormap(gray);
for ortNo=1:nOrt   
    subplot(nRows,nCols,ortNo);
    greyNaN(VERpca(:,:,ortNo));
    title([sprintf('%0g',Ort(ortNo)) '\circ'],'FontWeight','normal');
    axis square;
    addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),ortNo,nRows,nCols);
    colormap(gray);colorbar;caxis([-2 2]*10^-3)
end
suplabel('ROI-masked visual-evoked responses (PCA-ed)','t');
%ColorMap = [hsv(18)]*.9;
%export_fig(pdfFilename,'-pdf','-nocrop');

%% Diff PCA-ed data
figure('name','Subtr. PCA-ed responses')
colormap(gray);
for ortNo=1:nOrt/2  
    subplot(nRows,nCols,ortNo);
    greyNaN(VERpca(:,:,ortNo)-VERpca(:,:,ortNo+6));
    title([sprintf('%0g - %0g',Ort(ortNo),Ort(ortNo+6)) '\circ'],'FontWeight','normal');
    axis square;
    addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),ortNo,nRows,nCols);
    colormap(gray);colorbar;caxis([-2 2]*10^-3)
end
for ortNo=nOrt/2+1:nOrt 
    subplot(nRows,nCols,ortNo);
    greyNaN(VERpca(:,:,ortNo)-VERpca(:,:,ortNo-6));
    title([sprintf('%0g - %0g',Ort(ortNo),Ort(ortNo-6)) '\circ'],'FontWeight','normal');
    axis square;
    addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),ortNo,nRows,nCols);
    colormap(gray);colorbar;caxis([-2 2]*10^-3)
end
suplabel('ROI-masked visual-evoked responses (PCA-ed)','t');
%ColorMap = [hsv(18)]*.9;
%export_fig(pdfFilename,'-pdf','-nocrop');

%% Rescale
figure('name','Rescale')
for ortNo=1:nOrt
    ax=subplot(nRows,nCols,ortNo);
    rescaledCol(:,:,ortNo)=rescale(colVERpca(:,:,ortNo),-1,1);
    greyNaN(rescaledCol(:,:,ortNo));colorbar;
    
    title([sprintf('%0g',Ort(ortNo)) '\circ'],'FontWeight','normal');
    axis square;
    addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),ortNo,nRows,nCols);
    colormap(gray);colorbar;caxis('auto')%colormap(ax,[[0:.1:1]'.*ColorMap(ortNo,:)]);
end
suplabel(['Visual-evoked response (rescaled-PCA)'],'t');
%export_fig(pdfFilename,'-pdf','-nocrop','-append');

%% Rescale and histogram equalization to maximize dynamic range across image
figure('name','Contrast Eq.')
for ortNo=1:nOrt
    ax=subplot(nRows,nCols,ortNo);
    histEqCol(:,:,ortNo)=adapthisteq(colVERpca(:,:,ortNo),'NumTiles',imgDims./gridSize,'Range','full');

    %histEqCol(:,:,ortNo)=adapthisteq(rescaledCol(:,:,ortNo),'NumTiles',imgDims./gridSize,'Range','full');
    greyNaN(histEqCol(:,:,ortNo));colorbar;
    title([sprintf('%0g',Ort(ortNo)) '\circ'],'FontWeight','normal');
    axis square;
    addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),ortNo,nRows,nCols);
    colormap(gray);colorbar;caxis('auto')%colormap(ax,[[0:.1:1]'.*ColorMap(ortNo,:)]);
end
suplabel(['Visual-evoked response (rescaled-PCA + histogram equalization)'],'t');
%export_fig(pdfFilename,'-pdf','-nocrop','-append');

%% Gamma correction
figure('name','Gamma corr.')
for ortNo=1:nOrt
    ax=subplot(nRows,nCols,ortNo);
    gammaCol(:,:,ortNo)=imadjust(histEqCol(:,:,ortNo),[],[],gammaCorr);
    greyNaN(gammaCol(:,:,ortNo));
    title([sprintf('%0g',Ort(ortNo)) '\circ'],'FontWeight','normal');
    axis square;
    addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),ortNo,nRows,nCols);
    colormap(gray);colorbar;caxis('auto')%colormap(ax,[[0:.1:1]'.*ColorMap(ortNo,:)]);   
end
suplabel(['Visual-evoked response (rescaled-PCA + HE + gamma-corrected, gamma=' num2str(gammaCorr) ')'],'t');
%export_fig(pdfFilename,'-pdf','-nocrop','-append');


%% Threshold for bitmap
figure('name','Threshold')
for ortNo=1:nOrt
    ax=subplot(nRows,nCols,ortNo);

    imgDist=gammaCol(:,:,ortNo);
    threshold=prctile(imgDist(:),threshPrctile);

    thresholdMask=double(gammaCol(:,:,ortNo)>threshold);
    
    columnarMaps(:,:,ortNo)=(gammaCol(:,:,ortNo).*thresholdMask);
    
    thresholdMask(thresholdMask==0)=NaN;
    columnarBitmap(:,:,ortNo)=(gammaCol(:,:,ortNo).*thresholdMask)>0;
    
    greyNaN(columnarBitmap(:,:,ortNo),.7);
    title([sprintf('%0g',Ort(ortNo)) '\circ'],'FontWeight','normal');
    axis square;
    addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),ortNo,nRows,nCols);
    colormap(gray);colorbar;caxis([0 1])%colormap(ax,[[0:.1:1]'.*ColorMap(ortNo,:)]);
end
suplabel(['Visual-evoked response (rescaled-PCA + HE + gamma-corrected + threshold, Threshold=' num2str(threshPrctile) '%)'],'t');
%export_fig(pdfFilename,'-pdf','-nocrop','-append');

if plotOrNot==0
    set(groot,'defaultFigureVisible','on')
end
end