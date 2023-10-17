clear all;clc

%% From mainColstimPipeline
if ispc
  mainPath='Y:/';
elseif contains(getenv('HOSTNAME'),'psy.utexas.edu')
  mainPath='/eslab/data/';
end
saveFlag=0;
run([mainPath 'users/PK/colStimPipeline/exptListBiasingFull.m'])
currentSessID=62;%analysisSessID(sessionID);      
dsCurrentSess=datastruct(currentSessID); %fixed, to get ort map projected onto
referenceSessID=dsCurrentSess.referenceSessionEntryNo;
dsReferenceSess=datastruct(referenceSessID); %fixed, to get ort map projected onto

filenameStructCurrent=generateFilenames(dsCurrentSess);
filenameStructReference=generateFilenames(dsReferenceSess);

%% Stage 1. Load integrated cortical response, TS, gaussian ROI mask
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
    %make ROI mask  
    ROImask=double(Mask);
    ROImaskNaN=ROImask;
    ROImaskNaN(ROImaskNaN==0)=NaN; 
    %mask ort map, apply vasculature transform (morphs ort map to the imaging window view from experiment)
    PCAmap=RespCondPCA .* ROImaskNaN;
    PCAmap=transformImage(RespCondPCA,filenameStructCurrent);
    

end

%% Stage 2. For each intg image, split image into stimulated and nonstimulated area with gaussian ROI mask
%% Process integrated image to get the average image per condition (successfully completed trials only)
[trials,images,condIDs]=getUsableTrials(TS,DataTrial);
[images.averageColumn]=columnarfilter(TS,images.average);

%% Stage 3. Split the image into optostim and recruitment ROI
% Init
optostimVisual=[];optostimOrientation=[];
recruitmentVisMap=[];recruitmentPCAMap=[];
optostimROI=[];recruitROI=[];
% Define n-optostim conds and n-blanks
nOptoConds=numel(find(TS.Header.Conditions.TypeCond==3));
blank=condIDs.blankConds;
nblank=numel(blank);
% Extract gaussian levels used in experiment
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
    optostimMask=double(ROIMaskgaussian(gaussLevel).area);
    recruitmentMask=double(~optostimMask);
    optostimMask(optostimMask==0)=NaN;
    recruitmentMask(recruitmentMask==0)=NaN;
    
    %% Mask the optostim images
    for cond=1:numel(condsOptostim)
        condID=condsOptostim(cond);
        % unmasked for sanity check
        fullROI(:,:,condID-nblank)=(images.averageColumn(:,:,condID)) .* ROImaskNaN;

        % mask it twice, once for d' mask, then split into stimulated or nonstimulated area
        optostimROI(:,:,condID-nblank)=(images.averageColumn(:,:,condID)) .* optostimMask .* ROImaskNaN;
        recruitROI(:,:,condID-nblank)=(images.averageColumn(:,:,condID)) .* recruitmentMask .* ROImaskNaN;
    end
    
    %% Mask the visual images
    for cond=1:numel(condsVisual)
        condID=condsVisual(cond);
        % unmasked for sanity check
        fullROI(:,:,condID-nblank)=(images.averageColumn(:,:,condID)) .* ROImaskNaN;

        % mask it twice, once for d' mask, then split into stimulated or nonstimulated area
        optostimROI(:,:,condID-nblank)=(images.averageColumn(:,:,condID)) .* optostimMask .* ROImaskNaN;
        recruitROI(:,:,condID-nblank)=(images.averageColumn(:,:,condID)) .* recruitmentMask .* ROImaskNaN;
    end
end

%% Stage 4. Correlate to get difference map
%% Create difference map from the visual-only images (reference map)
condIDV0max=max(condIDs.V0);condIDV90max=max(condIDs.V90);
visualSubtMap=(images.averageColumn(:,:,condIDV90max)-images.averageColumn(:,:,condIDV0max)) .* ROImaskNaN;
pcaSubtMap=(PCAmap(:,:,7)-PCAmap(:,:,1)); %ROI mask already applied

for visualStim=1:2
    %% Define the optostim+visual and visual-only conditions for the correlation (congruent or incongruent), based on visual stimulus used
    if visualStim==1 % if opto 0 was applied
        % Q: for given visual 0, does opto 0 drive recruitment that resembles 0 more than opto 90? 
        condsVisual=condIDs.V0;
        condsOptostimCon=[condIDs.V0O0];
        condsOptostimIncon=[condIDs.V0O90];
    else % if opto 90 was applied
        % Q: for given visual 90, does opto 90 drive recruitment that resembles 90 more than opto 0? 
        condsVisual=condIDs.V90;
        condsOptostimCon=[condIDs.V90O90];
        condsOptostimIncon=[condIDs.V90O0];
    end
    
    %% Correlate the optostim images
    for cond=1:numel(condsOptostimCon)
        condCon=condsOptostimCon(cond);
        condIncon=condsOptostimIncon(cond);
        
        %compute correlation between difference image and difference map
        [~,~,corrValue,~]=calculateSimilarity(recruitROI(:,:,condCon-nblank) - recruitROI(:,:,condIncon-nblank), visualSubtMap);
        recruitmentVisMap(visualStim,cond)=corrValue;
        [~,~,corrValue,~]=calculateSimilarity(recruitROI(:,:,condCon-nblank) - recruitROI(:,:,condIncon-nblank), pcaSubtMap);
        recruitmentPCAMap(visualStim,cond)=corrValue;
        

        %compute correlation between difference image and difference map
        [~,~,corrValue,~]=calculateSimilarity(fullROI(:,:,condCon-nblank) - fullROI(:,:,condIncon-nblank), visualSubtMap);
        fullVisMap(visualStim,cond)=corrValue;
        [~,~,corrValue,~]=calculateSimilarity(fullROI(:,:,condCon-nblank) - fullROI(:,:,condIncon-nblank), pcaSubtMap);
        fullPCAMap(visualStim,cond)=corrValue;
    end
end

%% Stage 4A. Plot columnar responses
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
    img=fullROI(:,:,plotID);

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

%% Stage 4B. Plotting correlation of responses to visual and PCA difference reference map  
figure('name','Correlation: Responses vs visual & PCA reference ')
% Plot versus visual 0-90°
nCols=2;
nRows=1;
[hAx,~]=tight_subplot(nRows,nCols,[.2 .1],[.2 .2],[.1 .1]);

%% Full, versus PCA map
axes(hAx(1))
% stimulus conditions
stimCon=unique(TS.Header.Conditions.StimCon(find(TS.Header.Conditions.TypeCond==3)));
% congruent optostim+visual
opto0effect=recruitmentPCAMap(1,:);
opto90effect=recruitmentPCAMap(2,:);
% Plot congruent: visual 0 and 90, with congruent opto
plot(stimCon,opto0effect,'b-','LineWidth',2,'Marker','v','MarkerSize',10,'MarkerFaceColor','b','MarkerEdgeColor','k','DisplayName','Opto-0 (Vis 0 - 90)'); hold on;
plot(stimCon,opto90effect,'r-','LineWidth',2,'Marker','^','MarkerSize',10,'MarkerFaceColor','r','MarkerEdgeColor','k','DisplayName','Opto-90 (Vis 90 - 0)'); hold on;
% title and labels
title({'Recruitment', sprintf('4Hz flashed grating (SF: 2cpd, 100%%)')},'FontWeight', 'Normal')
xlim([min(stimCon) max(stimCon)]); xlabel('Absolute gabor contrast (%)');
line(xlim,[0 0],'color',.3*[1 1 1],'lineWidth',1.5,'LineStyle','--','HandleVisibility','off')
ylabel('Correlation with PCA map (r)');
hLeg=legend('location','southeast','NumColumns',1,'FontSize',14);
% cosmetic
upFontSize(14,.01)
ylim([-1 1])
% Add inset of subtracted image
axes('Position',[.33 .65 .15 .15])
imgsc(fullROI(:,:,max(condIDs.V0O0)-nblank) - fullROI(:,:,max(condIDs.V0O90)-nblank));
axis off; colormap fireice; caxis([-7.5 7.5]*10^-3)

%% Recruitment, versus PCA map
axes(hAx(2))
% stimulus conditions
stimCon=unique(TS.Header.Conditions.StimCon(find(TS.Header.Conditions.TypeCond==3)));
% congruent optostim+visual
opto0effect=-recruitmentPCAMap(1,:)*100./fullPCAMap(1,:);
opto90effect=recruitmentPCAMap(2,:)*100./fullPCAMap(2,:);
% Plot congruent: visual 0 and 90, with congruent opto
plot(stimCon,opto0effect,'b-','LineWidth',2,'Marker','v','MarkerSize',10,'MarkerFaceColor','b','MarkerEdgeColor','k','DisplayName','Opto-0 (Vis 0 - 90)'); hold on;
plot(stimCon,opto90effect,'r-','LineWidth',2,'Marker','^','MarkerSize',10,'MarkerFaceColor','r','MarkerEdgeColor','k','DisplayName','Opto-90 (Vis 90 - 0)'); hold on;
% title and labels
title({'Ratio of recruitment to full response', sprintf('4Hz flashed grating (SF: 2cpd, 100%%)')},'FontWeight', 'Normal')
xlim([min(stimCon) max(stimCon)]);ylabel('Ratio of recruitment to full response (%)');
line(xlim,[0 0],'color',.3*[1 1 1],'lineWidth',1.5,'LineStyle','--','HandleVisibility','off')
hLeg=legend('location','southeast','NumColumns',1,'FontSize',14);
% cosmetic
upFontSize(14,.01)
ylim([-100 100])
% Add inset of subtracted image
axes('Position',[.78 .65 .15 .15])
imgsc(recruitROI(:,:,max(condIDs.V0O0)-nblank) - recruitROI(:,:,max(condIDs.V0O90)-nblank));
axis off; colormap fireice; caxis([-7.5 7.5]*10^-3)
