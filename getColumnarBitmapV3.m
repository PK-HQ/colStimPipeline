function [columnarBitmap,VERpca]=getColumnarBitmapV3(dataStruct,gridSize,gammaCorr, sensitivity, uniqueStr, plotOrNot)
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
pdfFilename=filenameStruct.colmapPDF;
mapsFilename=filenameStruct.colmapMat;
%% Load raw and PCA-ed activity for all orientations, and SNR mask
load(filenameStruct.FFTAmp,'DataCond');
load(filenameStruct.Orientation,'Ort','RespCondPCA','Mask','MapOrt','MapAmpOrt','ColorMap'); %contains mask, RespCondPCA
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

%% Figure params
nRows=2;
nCols=6;

%% Extract masked VER and pca-ed VER
%ROI mask
ROImask=double(Mask);
ROImaskNaN=ROImask;
ROImaskNaN(ROImaskNaN==0)=NaN;
ROItuningMap=ROImask.*MapAmpOrt;

for ortNo=1:nOrt
    %VER for orientation, masked with SNR mask
    VERraw(:,:,ortNo)=DataCondBlankSubt(:,:,ortNo).*ROImaskNaN;
    VERpca(:,:,ortNo)=RespCondPCA(:,:,ortNo).*ROImaskNaN;
    
    %VER for orientation, masked with ROI mask (to isolate site)
    VERraw(:,:,ortNo)=VERraw(:,:,ortNo).*ROIMask;
    VERpca(:,:,ortNo)=VERpca(:,:,ortNo).*ROIMask;
end

%% PCA-ed data
figure('name','PCA-ed responses')
colormap(gray);
[hAx,~]=tight_subplot(2,6);
for ortNo=1:nOrt
    axes(hAx(ortNo));
    greyNaN(VERpca(:,:,ortNo));
    map=VERpca(:,:,ortNo);
    mapSD=nanstd(map,[],'all');
    mapMin=min(map(:));
    mapMax=max(map(:));
    titleL1=[sprintf('%0g',Ort(ortNo)) '\circ'];
    titleL2=['\sigma=' sprintf('%.1d',mapSD)];
    titleL3=['Range=[' sprintf('%.1d',mapMin) ', ' sprintf('%.1d',mapMax) ']'];
    title({titleL1,titleL2,titleL3},'FontWeight','normal');
    axis square;
    addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),ortNo,nRows,nCols);
    colormap(gray);colorbar;caxis([-2 2]*10^-3)
end
upFontSize(14,.015)
[ax,h]=suplabel('PCA visual-evoked responses (ROI masked, no-clip)','t');
set(h,'FontSize',16)
export_fig(pdfFilename,'-pdf','-nocrop');


%% PCA-ed data clipped
%figure('name','PCA-ed responses (clipped)')
%colormap(gray);
%[hAx,~]=tight_subplot(2,6);
for ortNo=1:nOrt
    %axes(hAx(ortNo));
    
    %rescale to norm
    map=VERpca(:,:,ortNo);
    normMap=rescale(map,0,1); 
    % stats
    mapMean=nanmean(normMap,'all');
    mapSD=nanstd(normMap,[],'all');
    mapMin=min(normMap(:));
    mapMax=max(normMap(:));   
    pos3SD=mapMean+3*mapSD;
    neg3SD=mapMean-3*mapSD;
   
    %clips image, remap to -1 to 1
    normMap = normMap .* ROImaskNaN;%imadjust(normMap,[neg3SD pos3SD],[0 1]) .* ROImaskNaN;
    normVERpca(:,:,ortNo)=normMap;
    mapSD=nanstd(normMap,[],'all');
    mapMin=min(normMap(:));
    mapMax=max(normMap(:));  
    
    %plot
    %{
    greyNaN(clippedMap);
    titleL1=[sprintf('%0g',Ort(ortNo)) '\circ'];
    titleL2=['\sigma=' sprintf('%.1d',mapSD)];
    titleL3=['Range=[' sprintf('%.1d',mapMin) ', ' sprintf('%.1d',mapMax) ']'];
    title({titleL1,titleL2,titleL3},'FontWeight','normal');
    axis square;
    addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),ortNo,nRows,nCols);
    colormap(gray);colorbar;%caxis([-2 2]*10^-3)
    %}
end
%suplabel('ROI-masked visual-evoked responses (PCA-ed, noclip)','t');
%export_fig(pdfFilename,'-pdf','-nocrop');


%% Histogram equalization to maximize dynamic range across image
figure('name','Contrast Equalization')
[hAx,~]=tight_subplot(2,6);
for ortNo=1:nOrt
    axes(hAx(ortNo));
    if ortNo==1 || ortNo==7
        histEqCol(:,:,ortNo)=adapthisteq(normVERpca(:,:,ortNo),'NumTiles',imgDims./gridSize,'Range','full');     %histEqImg=histeq(rescaled,1000);
    else
        histEqCol(:,:,ortNo)=zeros(size(normVERpca(:,:,ortNo),1),size(normVERpca(:,:,ortNo),1));
    end
    greyNaN(histEqCol(:,:,ortNo));colorbar;colormap(fireice)
    
    title([sprintf('%0g',Ort(ortNo)) '\circ'],'FontWeight','normal');
    axis square;
    addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),ortNo,nRows,nCols);
    colormap(gray);colorbar;caxis('auto')%colormap(ax,[[0:.1:1]'.*ColorMap(ortNo,:)]);
end
upFontSize(14,.015)
[ax,h]=suplabel('PCA visual-evoked responses (histogram equalization)','t');
set(h,'FontSize',16)
export_fig(pdfFilename,'-pdf','-nocrop','-append');

%% Gamma correction
figure('name','Gamma correction')
[hAx,~]=tight_subplot(2,6);
for ortNo=1:nOrt
    axes(hAx(ortNo));
    gammCorrImg=imadjust(histEqCol(:,:,ortNo),[],[],gammaCorr);
    gammaCol(:,:,ortNo)=gammCorrImg;
    greyNaN(gammaCol(:,:,ortNo));
    
    title([sprintf('%0g',Ort(ortNo)) '\circ'],'FontWeight','normal');
    axis square;
    addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),ortNo,nRows,nCols);
    colormap(gray);colorbar;caxis('auto')%colormap(ax,[[0:.1:1]'.*ColorMap(ortNo,:)]);   
end
upFontSize(14,.015)
[ax,h]=suplabel(['PCA visual-evoked responses (gamma correction=' num2str(gammaCorr) ')'],'t');
set(h,'FontSize',16)
export_fig(pdfFilename,'-pdf','-nocrop','-append');


%% Adaptive (local) threshold of the gamma map to obtain bitmap
figure('name','Adaptive threshold')
[hAx,~]=tight_subplot(2,6);
for ortNo=1:nOrt
    axes(hAx(ortNo));        
    % do adaptive thresholding
    threshold = adaptthresh(double(gammaCol(:,:,ortNo)), sensitivity);
    columnarBitmap(:,:,ortNo) = double(imbinarize(double(gammaCol(:,:,ortNo)),threshold)).*ROImask; %zero mask

    %calculate pixel density
    BWnan=columnarBitmap(:,:,ortNo).*ROImaskNaN;
    pixelDensity=nansum(BWnan(:))/numel(ROImaskNaN);
    
    %plot
    greyNaN(columnarBitmap(:,:,ortNo).*ROImask,.7);
    titleL1=[sprintf('%0g',Ort(ortNo)) '\circ'];
    titleL2=['PD=' sprintf('%.1d',pixelDensity)];
    title({titleL1,titleL2},'FontWeight','normal');
    axis square;
    addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),ortNo,nRows,nCols);
    colormap(gray);colorbar;caxis([0 1]);
end
upFontSize(14,.015)
[ax,h]=suplabel(['PCA visual-evoked responses (threshold sensitivity=' num2str(sensitivity) '%)'],'t');
set(h,'FontSize',16)
upFontSize(14,.015)
export_fig(pdfFilename,'-pdf','-nocrop','-append');


%% Overlay with columnar tuning map
figure('name','Overlay tuning map')
[hAx,~]=tight_subplot(2,6);
for ortNo=1:nOrt
    axes(hAx(ortNo));
    %plot
    greyNaN(columnarBitmap(:,:,ortNo).*ROItuningMap,.7);
    titleL1=[sprintf('%0g',Ort(ortNo)) '\circ'];
    title(titleL1,'FontWeight','normal');
    axis square;
    addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),ortNo,nRows,nCols);
    colormap(gray);colorbar;caxis('auto');
end
axes('Position',[.96,0.065,0.01,.85]);
imagesc(permute(flipdim(ColorMap,1),[1,3,2]));
set(gca,'YAxisLocation','Right', ...
        'XTick',[], ...
        'YTick',[1,9.5,18], ...
        'YTickLabel',{'180','90','0'});
ylabel('Orientation ^o', ...
       'FontSize',10, ...
       'FontWeight','Bold');
upFontSize(14,.015)
[ax,h]=suplabel('PCA visual-evoked responses (orientation map overlay)','t');
set(h,'FontSize',16)
export_fig(pdfFilename,'-pdf','-nocrop','-append');

%% Save maps
save(mapsFilename,'columnarBitmap','gammaCol','histEqCol','VERpca','normVERpca');

%% Show plot and goldilocks plot (for iterating sensitivity)
%{
figure('name','Threshold')


[~,OrtSubset]=find(Ort==[0 90]');


for ortNo=1:OrtSubset
    subplot(2,2,ortNo);
    
    imgDist=gammaCol(:,:,ortNo);
    
    % do adaptive thresholding
    sensitivity=0.01;
    threshold = adaptthresh(double(gammaCol(:,:,ortNo)), sensitivity);
    columnarBitmap(:,:,ortNo) = double(imbinarize(double(gammaCol(:,:,ortNo)),threshold));

    BWnan=columnarBitmap(:,:,ortNo).*ROImaskNaN;
    greyNaN(columnarBitmap(:,:,ortNo).*ROImask,.7);
    
    pixelDensity=nansum(BWnan(:))/numel(ROImaskNaN);
    
    titleL1=[sprintf('%0g',Ort(ortNo)) '\circ'];
    titleL2=['PD=' sprintf('%.1d',pixelDensity)];
    title({titleL1,titleL2},'FontWeight','normal');
    axis square;
    addPix2MM(min(imgX),max(imgX),min(imgY),max(imgY),ortNo,nRows,nCols);
    colormap(gray);colorbar;caxis([0 1])%colormap(ax,[[0:.1:1]'.*ColorMap(ortNo,:)]);
end
suplabel(['Visual-evoked response (rescaled-PCA + HE + gamma-corrected + threshold, Threshold=' num2str(threshPrctile) '%)'],'t');
export_fig(pdfFilename,'-pdf','-nocrop','-append');

%}



if plotOrNot==0
    set(groot,'defaultFigureVisible','on')
end
end