function assessOptoEffect(dataStructReference,dataStructCurrent)
%% Get filenames
filenameStructReference=generateFilenames(dataStructReference);
filenameStructCurrent=generateFilenames(dataStructCurrent);


%% Load raw and PCA-ed activity for all orientations, and SNR mask
%load(filenameStruct.Orientation,'Ort','RespCondPCA','Mask','MapOrt','MapAmpOrt','ColorMap'); %contains mask, RespCondPCA
%load(filenameStruct.colmapMat,'')
%load(filenameStruct.TS);

%% Short circuit
load(filenameStructCurrent.Intg, 'DataTrial') %DataTrial, DataCond
load(filenameStructCurrent.TS,'TS') %TS.Trial.Outcome, TS.Trial.CurrCond
load(filenameStructReference.Orientation,'Mask','MapAmpOrt','RespCondPCA')
load(filenameStructReference.colmapMat, 'columnarBitmap')
pdfFilename='test.pdf';

% ROI mask
ROImask=double(Mask);
ROImaskNaN=ROImask;
ROImaskNaN(ROImaskNaN==0)=NaN;
ROItuningMap=ROImask.*MapAmpOrt;

%% ---- Condition ID per trial, including OI discarded
%getTrialOutcomes, change to numel(find(outcomeTrials == 10 | outcomeTrials == 11 | outcomeTrials == -10 | outcomeTrials == -11))
% Get IDs of usable trials
outcomeTrials=extractfield(TS.Trial, 'Outcome');
condTrial=extractfield(TS.Trial, 'CurrCond');
uniqueConds=unique(condTrial);
nConds=numel(uniqueConds)-1;
successCode=10;
failedCode=11;
successSpecialCode=-10;
failedSpecialCode=-11;
completedTrialIdx=find(outcomeTrials == successCode | outcomeTrials == failedCode |...
    outcomeTrials == successSpecialCode | outcomeTrials == failedSpecialCode)'; 
usableTrialIdxTmp=find(outcomeTrials == successCode | outcomeTrials == failedCode); % indexes it to be within 1:nTotalTrials, e.g. 1:363
% get usable trials, with trial index based out of completed trials
[usableTrialIdx,~]=find(completedTrialIdx==usableTrialIdxTmp); % re-indexes index of 1:nTotalTrials to be within 1:nCompletedTrials, e.g. 1:266

%% Extract data for plotting
% get condition of each trial
usableTrialCond=condTrial(usableTrialIdxTmp)'; % grabs trials using index of 1:nTotalTrials
% get integrated images of usable trials
usableImages=DataTrial(:,:,usableTrialIdx);

% define condition indices
blank=1;
b0v0=[2 5 8 11];
b0v90=[14 17 20 23];
o0v0=b0v0+1;
o0v90=b0v90+1;
o90v0=b0v0+2;
o90v90=b0v90+2;
plotOrder=[b0v0' o0v0' o90v0' b0v90' o0v90' o90v90']-1;
plotOrder=reshape(plotOrder',[1 nConds]); %to plot correctly with tight_subplot
[condTrialIdx,~]=find(usableTrialCond==blank);
blankTrialsAvg=mean(usableImages(:,:,condTrialIdx),3);

% Spatial filtering for columns
for condID=2:nConds+1
    [condTrialIdx,~]=find(usableTrialCond==condID);
    condTrialAvg(:,:,condID-1)=mean(usableImages(:,:,condTrialIdx),3)-blankTrialsAvg;
    condTrialAvgColumnar(:,:,condID-1)=FilterFermi2D(condTrialAvg(:,:,condID-1), 0.8, 3, TS.Header.Imaging.SizePxl);
end

%% Plot raw integrated images
figure
nCols=6;
nRows=nConds/nCols;
[hAx,~]=tight_subplot(nRows,nCols,[.02 .02]);
for plotNo=1:nConds
    plotID=plotOrder(plotNo);
    img=condTrialAvg(:,:,plotID);
    
    axes(hAx(plotNo))
    imgsc(img.*ROImaskNaN)
    addPix2MM(1,512,1,512,plotNo,nRows,nCols);
    caxis([0 .1])
end
upFontSize(14,.015)
suplabel({'Raw signals',' ','Visual 0                                                                                                                                                                                            Visual 90'},'t',[.08 .08 .84 .76]);
export_fig(pdfFilename,'-pdf','-nocrop');

%% Plot spatially filtered integrated images
figure
nCols=6;
nRows=nConds/nCols;
[hAx,~]=tight_subplot(nRows,nCols,[.02 .02]);
colormap(fireice)
for plotNo=1:nConds
    colormap(fireice)
    plotID=plotOrder(plotNo);
    img=condTrialAvgColumnar(:,:,plotID);
    
    axes(hAx(plotNo))
    imgsc(img.*ROImaskNaN)
    addPix2MM(1,512,1,512,plotNo,nRows,nCols);
    caxis([-.007 .007])
end
colormap(fireice)
upFontSize(14,.015)
suplabel({'Columnar signals (bandpassed 0.8-3.0 cpmm)',' ','Visual 0                                                                                                                                                      Visual 90'},'t',[.08 .08 .84 .76]);
export_fig(pdfFilename,'-pdf','-nocrop','-append');


%% Plot correlation of each opto+vis condition with baseline+vis
% Calculate corr of a condSet of images () to a reference img (0)
% opto0 and opto90 versus b0v0
condSet=[o0v0,o0v90,o90v0,o90v90];
for i=1:length(condSet)
    img=condTrialAvgColumnar(:,:,condSet(i)-1).*ROImask;
    imgRef=condTrialAvgColumnar(:,:,b0v0(end)-1).*ROImask;
    b0v0results(i)=nansum(img.*imgRef);%corr2(img,imgRef);
end
for i=1:length(condSet)
    img=condTrialAvgColumnar(:,:,condSet(i)-1).*ROImask;
    imgRef=condTrialAvgColumnar(:,:,b0v90(end)-1).*ROImask;
    b0v90results(i)=nansum(img.*imgRef);%corr2(img,imgRef);
end

% split into lines for each comparison condSet
nComparisons=numel(b0v0results);
o0v0_b0v0=b0v0results(1:nComparisons/4);
o0v90_b0v0=b0v0results(nComparisons/4 + 1 : nComparisons*2/4);
o90v0_b0v0=b0v0results(nComparisons*2/4 + 1 : nComparisons*3/4);
o90v90_b0v0=b0v0results(nComparisons*3/4 + 1 : nComparisons);

% split into lines for each comparison condSet
o0v0_b0v90=b0v90results(1:nComparisons/4);
o0v90_b0v90=b0v90results(nComparisons/4 + 1 : nComparisons*2/4);
o90v0_b0v90=b0v90results(nComparisons*2/4 + 1 : nComparisons*3/4);
o90v90_b0v90=b0v90results(nComparisons*3/4 + 1 : nComparisons);

stimCon=unique(TS.Header.Conditions.StimCon);
stimCon=stimCon(stimCon>0);

% Plotting
figure
% Plot versus visual 0°
nCols=3;
nRows=1;
[hAx,~]=tight_subplot(nRows,nCols,[.1 .1]);
axes(hAx(1))
plot(stimCon,o0v0_b0v0,'b-','LineWidth',2,'DisplayName','O0 + V0'); hold on;
plot(stimCon,o90v0_b0v0,'r:.','LineWidth',2,'DisplayName','O90 + V0'); hold on;
plot(stimCon,o0v90_b0v0,'b:','LineWidth',2,'DisplayName','O0 + V90'); hold on;
plot(stimCon,o90v90_b0v0,'r-','LineWidth',2,'DisplayName','O90 + V90'); hold on;
title(sprintf('Similarity to 0° gabor (contrast=%.0f%%)',max(stimCon)),'FontWeight', 'Normal')
xlabel('Contrast (%)');xlim([min(stimCon) max(stimCon)])
ylabel('Correlation (r)');ylim([.0 1])
legend('location','southeast')
upFontSize(14,.01)

% Plot versus visual 90°
axes(hAx(2))
plot(stimCon,o0v0_b0v90,'b-','LineWidth',2,'DisplayName','O0 + V0'); hold on;
plot(stimCon,o90v0_b0v90,'r:','LineWidth',2,'DisplayName','O90 + V0'); hold on;
plot(stimCon,o0v90_b0v90,'b:','LineWidth',2,'DisplayName','O0 + V90'); hold on;
plot(stimCon,o90v90_b0v90,'r-','LineWidth',2,'DisplayName','O90 + V90'); hold on;
title(sprintf('Similarity to 90° gabor (contrast=%.0f%%)',max(stimCon)),'FontWeight', 'Normal')
xlabel('Contrast (%)');xlim([min(stimCon) max(stimCon)])
ylabel('Correlation (r)');ylim([.0 1])
legend('location','southeast')
upFontSize(14,.01)

% Plot versus visual 90°
axes(hAx(3))
plot(stimCon,o0v0_b0v0-o0v0_b0v90,'b-','LineWidth',2,'DisplayName','O0 + V0'); hold on;
plot(stimCon,o90v0_b0v0-o90v0_b0v90,'r:','LineWidth',2,'DisplayName','O90 + V0'); hold on;
plot(stimCon,o0v90_b0v0-o0v90_b0v90,'b:','LineWidth',2,'DisplayName','O0 + V90'); hold on;
plot(stimCon,o90v90_b0v0-o90v90_b0v90,'r-','LineWidth',2,'DisplayName','O90 + V90'); hold on;
% zero line
line([0 50],[0 0],'Color',[0 0 0],'LineStyle','--','LineWidth',.25,'HandleVisibility','off')
title(sprintf('Similarity index (similarity to 0° - 90°'),'FontWeight', 'Normal')
xlabel('Contrast (%)');xlim([min(stimCon) max(stimCon)])
ylabel('Similarity index (r)');ylim([-1 1])
legend('location','southeast')
upFontSize(14,.01)
suplabel({'Comparison of opto+visual with black bmp+visual'},'t',[.08 .08 .84 .76]);
upFontSize(14,.01)

          % + similarity to original PCA map and bitmap?
          
          
%% Plot correlation of each opto+vis condition with bitmap/PCA map
% Calculate corr of a condSet of images () to a reference img (0)
% opto0 and opto90 versus b0v0
condSet=[o0v0,o0v90,o90v0,o90v90];
for i=1:length(condSet)
    img=condTrialAvgColumnar(:,:,condSet(i)-1).*ROImask;
    imgRef=columnarBitmap(:,:,1).*ROImask;
    pca0results(i)=nansum(img.*imgRef);%corr2(img,imgRef);
end
for i=1:length(condSet)
    img=condTrialAvgColumnar(:,:,condSet(i)-1).*ROImask;
    imgRef=columnarBitmap(:,:,7).*ROImask;
    pca90results(i)=nansum(img.*imgRef);%corr2(img,imgRef);
end

% split into lines for each comparison condSet
nComparisons=numel(pca0results);
o0v0_pca0=pca0results(1:nComparisons/4);
o0v90_pca0=pca0results(nComparisons/4 + 1 : nComparisons*2/4);
o90v0_pca0=pca0results(nComparisons*2/4 + 1 : nComparisons*3/4);
o90v90_pca0=pca0results(nComparisons*3/4 + 1 : nComparisons);

% split into lines for each comparison condSet
o0v0_pca90=pca90results(1:nComparisons/4);
o0v90_pca90=pca90results(nComparisons/4 + 1 : nComparisons*2/4);
o90v0_pca90=pca90results(nComparisons*2/4 + 1 : nComparisons*3/4);
o90v90_pca90=pca90results(nComparisons*3/4 + 1 : nComparisons);

stimCon=unique(TS.Header.Conditions.StimCon);
stimCon=stimCon(stimCon>0);

% Plotting
figure
% Plot versus visual 0Â°
nCols=3;
nRows=1;
[hAx,~]=tight_subplot(nRows,nCols,[.1 .1]);
axes(hAx(1))
plot(stimCon,o0v0_pca0,'b-','LineWidth',2,'DisplayName','O0 + V0'); hold on;
plot(stimCon,o90v0_pca0,'r:.','LineWidth',2,'DisplayName','O90 + V0'); hold on;
plot(stimCon,o0v90_pca0,'b:','LineWidth',2,'DisplayName','O0 + V90'); hold on;
plot(stimCon,o90v90_pca0,'r-','LineWidth',2,'DisplayName','O90 + V90'); hold on;
title(sprintf('Similarity to 0° PCA (contrast=%.0f%%)',max(stimCon)),'FontWeight', 'Normal')
xlabel('Contrast (%)');xlim([min(stimCon) max(stimCon)])
ylabel('Correlation (r)');ylim([-1 1])
legend('location','southeast')
upFontSize(14,.01)

% Plot versus visual 90Â°
axes(hAx(2))
plot(stimCon,o0v0__pca90,'b-','LineWidth',2,'DisplayName','O0 + V0'); hold on;
plot(stimCon,o90v0__pca90,'r:','LineWidth',2,'DisplayName','O90 + V0'); hold on;
plot(stimCon,o0v90__pca90,'b:','LineWidth',2,'DisplayName','O0 + V90'); hold on;
plot(stimCon,o90v90__pca90,'r-','LineWidth',2,'DisplayName','O90 + V90'); hold on;
title(sprintf('Similarity to 90° PCA (contrast=%.0f%%)',max(stimCon)),'FontWeight', 'Normal')
xlabel('Contrast (%)');xlim([min(stimCon) max(stimCon)])
ylabel('Correlation (r)');ylim([-1 1])
legend('location','southeast')
upFontSize(14,.01)

% Plot versus visual 90Â°
axes(hAx(3))
plot(stimCon,o0v0_pca0-o0v0_pca90,'b-','LineWidth',2,'DisplayName','O0 + V0'); hold on;
plot(stimCon,o90v0_pca0-o90v0_pca90,'r:','LineWidth',2,'DisplayName','O90 + V0'); hold on;
plot(stimCon,o0v90_pca0-o0v90_pca90,'b:','LineWidth',2,'DisplayName','O0 + V90'); hold on;
plot(stimCon,o90v90_pca0-o90v90_pca90,'r-','LineWidth',2,'DisplayName','O90 + V90'); hold on;
% zero line
line([0 50],[0 0],'Color',[0 0 0],'LineStyle','--','LineWidth',.25,'HandleVisibility','off')
title(sprintf('Similarity index (similarity to 0° - 90°'),'FontWeight', 'Normal')
xlabel('Contrast (%)');xlim([min(stimCon) max(stimCon)])
ylabel('Similarity to 0 (r)');ylim([-1 1])
legend('location','southeast')
upFontSize(14,.01)
suplabel({'Comparison of opto+visual with PCA'},'t',[.08 .08 .84 .76]);
upFontSize(14,.01)
export_fig(pdfFilename,'-pdf','-nocrop','-append');

end
          
          
