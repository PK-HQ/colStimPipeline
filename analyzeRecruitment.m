clear all;clc

%% From mainColstimPipeline
if ispc
  mainPath='Y:/';
elseif contains(getenv('HOSTNAME'),'psy.utexas.edu')
  mainPath='/eslab/data/';
end
run([mainPath 'users/PK/colStimPipeline/exptListBiasingFull.m'])
currentSessID=64;%analysisSessID(sessionID);      
dsCurrentSess=datastruct(currentSessID); %fixed, to get ort map projected onto
referenceSessID=dsCurrentSess.referenceSessionEntryNo;
dsReferenceSess=datastruct(referenceSessID); %fixed, to get ort map projected onto

filenameStructCurrent=generateFilenames(dsCurrentSess);
filenameStructReference=generateFilenames(dsReferenceSess);

%% 1. Load integrated cortical response, TS, gaussian ROI mask
if isfile(filenameStructCurrent.Intg) % Load imaging data if it exists 
    %% Load data
    % TS
    load(filenameStructCurrent.TS,'TS') %TS.Trial.Outcome, TS.Trial.CurrCond
    % .intg file
    load(filenameStructCurrent.Intg, 'DataTrial') %DataTrial, DataCond
    % .gaussian file
    load(filenameStructCurrent.gaussianFit)
    % d' ROI mask, RespCondPCA
    load(filenameStructReference.Orientation,'Mask','RespCondPCA')
    ortmap=RespCondPCA;
    ortmap=transformImage(RespCondPCA,filenameStructCurrent);
    
    ROImask=double(Mask);
    %ROImask=transformImage(ROImask,filenameStructCurrent);
    ROImaskNaN=ROImask;
    ROImaskNaN(ROImaskNaN==0)=NaN;
end

%% 2. For each intg image, split image into stimulated and nonstimulated area with gaussian ROI mask
%% Process integrated image to get the average image per condition (successfully completed trials only)
[trials,images,condIDs]=getUsableTrials(TS,DataTrial);
[images.averageColumn]=columnarfilter(TS,images.average);

%% Split the image into optostim and recruitment ROI
optostimVisual=[];optostimOrientation=[];
recruitmentVisual=[];recruitmentOrientation=[];
optostimROI=[];recruitROI=[];
nOptoConds=numel(find(TS.Header.Conditions.TypeCond==3));
blank=condIDs.blankConds;
nblank=numel(blank);
gaussLevels=dsCurrentSess.gaussianContourLevel; % grab the gaussian levels used
for gaussNo=1:numel(gaussLevels)
    %% Define the optostim+visual and visual-only conditions for the optostim ROI mask (0 or 90), based on which optostim was applied during experiment
    if gaussNo==1 % if opto 0 was applied
        % Q: does opto 0 drive recruitment similar to visual 0? 
        condsVisual=condIDs.V0;
        condsOptostim=[condIDs.V0O0 condIDs.V90O0];
    else % if opto 90 was applied
        % Q: does opto 90 drive recruitment similar to visual 90? 
        condsVisual=condIDs.V90;
        condsOptostim=[condIDs.V0O90 condIDs.V90O90];
    end
    
    %% Create binary mask of each image, with the mask = gaussian or inverse of the gaussian
    % grab the gaussian level previously used to create optostim ROI (out of gaussian max level)
    gaussLevel=max(gaussLevels); %gaussLevels(gaussNo); 
    % define the stimulated ROI (gaussian) and nonstimulated ROI (inverse)
    optostimMask=double(ROIMaskgaussian(gaussLevel).area);recruitmentMask=double(~optostimMask);
    optostimMask(optostimMask==0)=NaN;
    recruitmentMask(recruitmentMask==0)=NaN;
   
    %% Create difference map from the visual-only images (reference map)
    condIDV0max=max(condIDs.V0);condIDV90max=max(condIDs.V90);
    visualSubtMap=images.averageColumn(:,:,condIDV90max)-images.averageColumn(:,:,condIDV0max);
    ortSubtMap=ortmap(:,:,7)-ortmap(:,:,1);
    
    % mask the optostim images
    for cond=1:numel(condsOptostim)
        condID=condsOptostim(cond);
        % unmasked for sanity check
        fullROI(:,:,condID-nblank)=(images.averageColumn(:,:,condID)) .* ROImaskNaN;

        % mask it twice, once for d' mask, then split into stimulated or nonstimulated area
        optostimROI(:,:,condID-nblank)=(images.averageColumn(:,:,condID)) .* optostimMask .* ROImaskNaN;
        recruitROI(:,:,condID-nblank)=(images.averageColumn(:,:,condID)) .* recruitmentMask .* ROImaskNaN;
        
        %% Correlate responses in the recruitment area to a difference map (high contrast visual, subtracted)
        [~,~,corrValue,~]=calculateSimilarity(recruitROI(:,:,condID-nblank), visualSubtMap .* recruitmentMask .* ROImaskNaN);
        recruitmentVisual(condID-nblank)=corrValue;
        
        [~,~,corrValue,~]=calculateSimilarity(recruitROI(:,:,condID-nblank), ortSubtMap .* recruitmentMask .* ROImaskNaN);
        recruitmentOrientation(condID-nblank)=corrValue;
        
        %% Correlate responses in the stimulated area to a difference map (high contrast visual, subtracted)
        [~,~,corrValue,~]=calculateSimilarity(optostimROI(:,:,condID-nblank), visualSubtMap .* optostimMask .* ROImaskNaN);
        optostimVisual(condID-nblank)=corrValue;
        
        [~,~,corrValue,~]=calculateSimilarity(optostimROI(:,:,condID-nblank), ortSubtMap .* optostimMask .* ROImaskNaN);
        optostimOrientation(condID-nblank)=corrValue;
        
    end
    
    % mask the visual images
    for cond=1:numel(condsVisual)
        condID=condsVisual(cond);
        % unmasked for sanity check
        fullROI(:,:,condID-nblank)=(images.averageColumn(:,:,condID)) .* ROImaskNaN;

        % mask it twice, once for d' mask, then split into stimulated or nonstimulated area
        optostimROI(:,:,condID-nblank)=(images.averageColumn(:,:,condID)) .* optostimMask .* ROImaskNaN;
        recruitROI(:,:,condID-nblank)=(images.averageColumn(:,:,condID)) .* recruitmentMask .* ROImaskNaN;
        
        %% Correlate responses in the recruitment area to a difference map (high contrast visual, subtracted)
        [~,~,corrValue,~]=calculateSimilarity(recruitROI(:,:,condID-nblank), visualSubtMap .* recruitmentMask .* ROImaskNaN);
        recruitmentVisual(condID-nblank)=corrValue;
        
        [~,~,corrValue,~]=calculateSimilarity(recruitROI(:,:,condID-nblank), ortSubtMap .* recruitmentMask .* ROImaskNaN);
        recruitmentOrientation(condID-nblank)=corrValue;
        
        %% Correlate responses in the stimulated area to a difference map (high contrast visual, subtracted)
        [~,~,corrValue,~]=calculateSimilarity(optostimROI(:,:,condID-nblank), visualSubtMap .* optostimMask .* ROImaskNaN);
        optostimVisual(condID-nblank)=corrValue;
        
        [~,~,corrValue,~]=calculateSimilarity(optostimROI(:,:,condID-nblank), ortSubtMap .* optostimMask .* ROImaskNaN);
        optostimOrientation(condID-nblank)=corrValue;
        
    end
end

%% Expectation
%%-incon +con



%% Plot images
%% 1-A. Plot spatially filtered images
figure('name','Columnar responses')
nCols=6;

nRows=nOptoConds/nCols;
% Define plot order
plotOrder=[condIDs.V0' condIDs.V0O0' condIDs.V0O90' condIDs.V90' condIDs.V90O0' condIDs.V90O90']-nblank;
plotOrder=reshape(plotOrder',[1 nOptoConds]); %to plot correctly with tight_subplot
[hAx,~]=tight_subplot(nRows,nCols,[.02 .02]);
colormap(fireice)

% --- Plot ---
for plotNo=1:nOptoConds
    colormap(fireice)
    plotID=plotOrder(plotNo);
    img=recruitROI(:,:,plotID);

    axes(hAx(plotNo))
    imgsc(img)
    addPix2MM(1,512,1,512,plotNo,nRows,nCols);
    caxis([-1 1]*10^-2)

    if plotID==1
      title('Baseline','FontWeight','Normal')
    elseif plotID==2
      title('Opto 0\circ','FontWeight','Normal')
    elseif plotID==3
      title('Opto 90\circ','FontWeight','Normal')
    elseif plotID==nRows*nCols*0.5 + 1
      title('Baseline','FontWeight','Normal')
    elseif plotID==nRows*nCols*0.5 + 2
      title('Opto 0\circ','FontWeight','Normal')
    elseif plotID==nRows*nCols*0.5 + 3
      title('Opto 90\circ','FontWeight','Normal')
    end
    colorbar
end
%colormap(gray)
upFontSize(14,.015)
[~,h]=suplabel({'Columnar signals (bandpassed 0.8-3.0 cpmm)',' ','Visual 0\circ                                                                                                                             Visual 90\circ'},'t',[.08 .08 .84 .81]);

colormap(fireice)




% Save the figures
set(h,'FontSize',16)
if saveFlag
    export_fig(pdfFilename,'-pdf','-nocrop','-append');
end


%% X. Plotting correlation of responses to visual and PCA difference reference map  
% --- Visual reference map ---
figure('name','Correlation: Responses vs visual & PCA reference ')
% Plot versus visual 0-90°
nCols=2;
nRows=2;
[hAx,~]=tight_subplot(nRows,nCols,[.1 .15],[.05 .15],[.04 .04]);
axes(hAx(2))

%% Define plot x-y
% stimulus conditions
stimCon=unique(TS.Header.Conditions.StimCon(find(TS.Header.Conditions.TypeCond==3)));
% visual-only
V0xVisualr=recruitmentVisual([condIDs.V0]-2);
V90xVisualr=recruitmentVisual([condIDs.V90]-2);
% congruent optostim+visual
O0V0xVisualr=recruitmentVisual([condIDs.V0O0]-2);
O90V90xVisualr=recruitmentVisual([condIDs.V90O90]-2);
% incongruent optostim+visual
O90V0xVisualr=recruitmentVisual([condIDs.V0O90]-2);
O0V90xVisualr=recruitmentVisual([condIDs.V90O0]-2);


% visual-only
V0xVisualr=recruitmentOrientation([condIDs.V0]-2);
V90xVisualr=recruitmentOrientation([condIDs.V90]-2);
% congruent optostim+visual
O0V0xVisualr=recruitmentOrientation([condIDs.V0O0]-2);
O90V90xVisualr=recruitmentOrientation([condIDs.V90O90]-2);
% incongruent optostim+visual
O90V0xVisualr=recruitmentOrientation([condIDs.V0O90]-2);
O0V90xVisualr=recruitmentOrientation([condIDs.V90O0]-2);

% Plot baselines: visual 0 and 90, without opto
plot(stimCon,V0xVisualr,'k-','LineWidth',2,'Marker','o','MarkerSize',10,'MarkerFaceColor','b','MarkerEdgeColor','k','DisplayName','Vis-0'); hold on;
plot(stimCon,V90xVisualr,'k-','LineWidth',2,'Marker','o','MarkerSize',10,'MarkerFaceColor','r','MarkerEdgeColor','k','DisplayName','Vis-90'); hold on;

% Plot congruent: visual 0 and 90, with congruent opto
plot(stimCon,O0V0xVisualr,'b-','LineWidth',2,'Marker','v','MarkerSize',10,'MarkerFaceColor','b','MarkerEdgeColor','k','DisplayName','Opto-0 + Vis-0'); hold on;
plot(stimCon,O90V90xVisualr,'r-','LineWidth',2,'Marker','^','MarkerSize',10,'MarkerFaceColor','r','MarkerEdgeColor','k','DisplayName','Opto-90 + Vis-90'); hold on;

% Plot incongruent: visual 0 and 90, with incongruent opto
plot(stimCon,O90V0xVisualr,'b--','LineWidth',2,'Marker','^','MarkerSize',10,'MarkerFaceColor','r','MarkerEdgeColor','k','DisplayName','Opto-90 + Vis-0'); hold on;
plot(stimCon,O0V90xVisualr,'r--','LineWidth',2,'Marker','v','MarkerSize',10,'MarkerFaceColor','b','MarkerEdgeColor','k','DisplayName','Opto-0 + Vis-90'); hold on;

% title and labels
title({'Subtracted visual response', sprintf('Single gabor (S.F.: 6cpd, Con.: %.0f%%)',max(stimCon))},'FontWeight', 'Normal')
xlabel('Absolute gabor contrast (%)');xlim([min(stimCon) max(stimCon)])
ylabel('Correlation (r)');
hLeg=legend('location','southeast','NumColumns',3,'FontSize',10);

% cosmetic
upFontSize(14,.01)
ylim([-1 1])


% --- PCA reference map ---
axes(hAx(2))
% Plot baselines: visual 0 and 90, without opto
plot(stimCon,V0xPCAr,'b-','LineWidth',2,'Marker','o','MarkerSize',10,'MarkerFaceColor','w','MarkerEdgeColor','k','DisplayName','Vis-0'); hold on;
plot(stimCon,V90xPCAr,'r-','LineWidth',2,'Marker','o','MarkerSize',10,'MarkerFaceColor','w','MarkerEdgeColor','k','DisplayName','Vis-90'); hold on;

% Plot congruent: visual 0 and 90, with congruent opto
plot(stimCon,O0V0xPCAr,'b--','LineWidth',2,'Marker','v','MarkerSize',10,'MarkerFaceColor','b','MarkerEdgeColor','k','DisplayName','Opto-0 + Vis-0'); hold on;
plot(stimCon,O90V90xPCAr,'r--','LineWidth',2,'Marker','^','MarkerSize',10,'MarkerFaceColor','r','MarkerEdgeColor','k','DisplayName','Opto-90 + Vis-90'); hold on;

% Plot incongruent: visual 0 and 90, with incongruent opto
plot(stimCon,O90V0xPCAr,'b:','LineWidth',2,'Marker','^','MarkerSize',10,'MarkerFaceColor','r','MarkerEdgeColor','k','DisplayName','Opto-90 + Vis-0'); hold on;
plot(stimCon,O0V90xPCAr,'r:','LineWidth',2,'Marker','v','MarkerSize',10,'MarkerFaceColor','b','MarkerEdgeColor','k','DisplayName','Opto-0 + Vis-90'); hold on;

% title and labels
title({'Subtracted PCA response', sprintf('Flashed 4Hz grating (S.F.: 2cpd, Con.: %.0f%%)',100)},'FontWeight', 'Normal')
xlim([min(stimCon) max(stimCon)])
ylim([-1 1]);
upFontSize(14,.01)

% cosmetic
hLeg.FontSize=11;
[~,hSupLabel]=suplabel({'Correlation between opto-evoked and desired response'},'t',[.08 .08 .84 .84]);
set(hSupLabel,'FontSize',16)




