function [visualr,visualrRatios,PCAr,PCArRatios,optoIntensity]=assessOptoEffectChoice(dataStructReference,dataStructCurrent,columnarBitmapCoregistered,columnarPCAsCoregistered,saveFlag)
%% Get filenames
filenameStructReference=generateFilenames(dataStructReference);
filenameStructCurrent=generateFilenames(dataStructCurrent);
trialsOutcomeLabels={'all','correct','error'};
%% Load integrated cortical response, TS, ROI mask
if isfile(filenameStructCurrent.Intg) % Load imaging data if it exists 
    %% Load data
    load(filenameStructCurrent.Intg, 'DataTrial') %DataTrial, DataCond
    load(filenameStructCurrent.TS,'TS') %TS.Trial.Outcome, TS.Trial.CurrCond
    load(filenameStructReference.Orientation,'Mask','RespCondPCA')
    pdfFilename=filenameStructCurrent.neurometricPDF;

    %% Coregister ROI mask
    % Coregister
    [~,~,transformParams]=coregisterGreenImages(dataStructReference,dataStructCurrent,Mask);
    
    % set ROI mask
    ROImask=double(Mask);
    ROImask=double(imwarp(ROImask,transformParams,'OutputView',imref2d(size(imgTarget))));
    ROImaskNaN=ROImask;
    ROImaskNaN(ROImaskNaN==0)=NaN;

    %% Get usable trials, indexed appropriately to match cortical response images (find trials in TS with success or failure code only)
    [trialIdx,images,condIDs]=getUsableTrials(TS,DataTrial);

    %% Define conditions, plot order
    nOptoConds=numel(find(TS.Header.Conditions.TypeCond==3));
    blank=condIDs.blankConds;
    nblank=numel(blank);
    nConds=6;
    
    % Define condition indices
    [V0, V90, V0O0, V90O0, V0O90, V90O90]=getCondIDs(condIDs);

    % Define plot order
    plotOrder=[V0' V0O0' V0O90' V90' V90O0' V90O90']-nblank;
    plotOrder=reshape(plotOrder',[1 nOptoConds]); %to plot correctly with tight_subplot
    [condTrialIdx,~]=find(trialIdx.All==blank);
    blankTrialsAvg=mean(images.Average(:,:,condTrialIdx),3);

    
    for trialOutcome=1:3
      % Get the trial outcome label (all, correct, error)
      trialsOutcomeLabel=trialsOutcomeLabels{trialOutcome};
      
      % Define type of trials used to analyze
      switch trialsOutcomeLabel
        case {'all'}
          trialsWanted=trialIdx.All;
        case {'correct'}
          trialsWanted=trialIdx.Correct;
        case {'error'}
          trialsWanted=trialIdx.Error;
      end
      
    
      %% 1. Spatial bandpass to extract columnar-scale responses
      for optoCondID=1:nOptoConds%2:nConds+1
          [condTrialIdx,~]=find(trialsWanted==optoCondID+nblank); %to get opto cond ID, use 1:nOptoConds + number of blank
          condTrialAvg(:,:,optoCondID)=mean(images.Average(:,:,condTrialIdx),3)-blankTrialsAvg;
          columnarResponse(:,:,optoCondID)=FilterFermi2D(condTrialAvg(:,:,optoCondID), 0.8, 3, TS.Header.Imaging.SizePxl);
      end

      %% 1-A. Plot spatially filtered images
      figure('name',['Columnar responses ' trialsOutcomeLabel])
      nCols=6;
      nRows=nOptoConds/nCols;
      [hAx,~]=tight_subplot(nRows,nCols,[.02 .02]);
      colormap(fireice)

      % --- Plot ---
      for plotNo=1:nOptoConds
          colormap(fireice)
          plotID=plotOrder(plotNo);
          img=columnarResponse(:,:,plotID);

          axes(hAx(plotNo))
          imgsc(img.*ROImaskNaN)
          addPix2MM(1,512,1,512,plotNo,nRows,nCols);
          caxis([-.007 .007])

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
      colormap(gray)
      upFontSize(14,.015)
      [~,h]=suplabel({'Columnar signals (bandpassed 0.8-3.0 cpmm)',' ','Visual 0\circ                                                                                                                             Visual 90\circ'},'t',[.08 .08 .84 .81]);
      % Save the figures
      set(h,'FontSize',16)
      if saveFlag
          export_fig(pdfFilename,'-pdf','-nocrop','-append');
      end

      %% 1-B. Plot subtracted spatially filtered images
      figure('name',['Subtracted columnar responses ' trialsOutcomeLabel])
      nCols=4;
      nRows=1;
      [hAx,~]=tight_subplot(nRows,nCols,[.02 .02]);
      colormap(fireice)

      % --- Get subtracted images ---
      % Difference between max contrast visual stimuli (visual 0 - visual 90
      visual_difference=columnarResponse(:,:,V90(end)-nblank)-columnarResponse(:,:,V0(end)-nblank);

      % Difference between PCA images
      pca_difference=columnarPCAsCoregistered(:,:,7)-columnarPCAsCoregistered(:,:,1);

      % Difference between opto 0 - opto 90
      %opto_difference_visual0=(columnarResponse(:,:,V0O90(1)-nblank)-columnarResponse(:,:,V0O0(1)-nblank));
      %opto_difference_visual90=(columnarResponse(:,:,V90O90(1)-nblank)-columnarResponse(:,:,V90O0(1)-nblank));
      %opto_difference=(opto_difference_visual0+opto_difference_visual90)/2;

      opto_0=((columnarResponse(:,:,V0O0(1))+columnarResponse(:,:,V90O0(1)))-nblank)./2;
      opto_90=((columnarResponse(:,:,V0O90(1))+columnarResponse(:,:,V90O90(1)))-nblank)./2;
      opto_difference = opto_0-opto_90;

      % bmp difference
      bitmap_difference=columnarBitmapCoregistered(:,:,7)-columnarBitmapCoregistered(:,:,1);

      % compile
      condTrialAvgDiffColumnar(:,:,1)=visual_difference;
      condTrialAvgDiffColumnar(:,:,2)=pca_difference;
      condTrialAvgDiffColumnar(:,:,3)=bitmap_difference;
      condTrialAvgDiffColumnar(:,:,4)=opto_difference;


      % --- Plot ---
      for plotNo=1:4
          img=condTrialAvgDiffColumnar(:,:,plotNo);
          axes(hAx(plotNo))
          imgsc(img.*ROImaskNaN)
          addPix2MM(1,512,1,512,plotNo,nRows,nCols);
          %caxis([-.007 .007])
          colorbar
          if plotNo==1
            title('Visual-only','FontWeight','Normal')
            caxis([-8.5 8.5]*10^-3)
          elseif plotNo==2
            title('Denoised visual-only','FontWeight','Normal')
            caxis([-4 4]*10^-3)
          elseif plotNo==3
            title('Bitmap','FontWeight','Normal')
            caxis([-1 1])
          elseif plotNo==4
            title('Optostimulation','FontWeight','Normal')
            caxis([-5 5]*10^-3)
          end
          ax = gca;
          ax.XLim = [1.5 6.5] * 512/8.22;
          ax.YLim = [1.5 6] * 512/8.22;
      end
      colormap(fireice)
      upFontSize(14,.0125)
      [~,h]=suplabel('Differences in columnar signals (H - V)','t',[.08 .08 .84 .81]);
      %% Save the figures
      set(h,'FontSize',16)
      if saveFlag
          export_fig(pdfFilename,'-pdf','-nocrop','-append');
      end
















      %% Calculate correlation and intensity of each opto+vis condition with baseline+vis
      % --- Calculate correlation & intensity between each condition's average image to a reference visual baseline
      % (max contrast) ---  
      conditionSet=[V0O0,V90O0,V0O90,V90O90,V0,V90];
      conditionSetversusV0=[V0O90];
      conditionSetversusV90=[V90O0];
      for i=1:length(conditionSet)
          img=columnarResponse(:,:,conditionSet(i)-nblank).*ROImask;
          imgRef=visual_difference.*ROImask;%columnarBitmap(:,:,1).*ROImask;
          [~,~,corrValue,~]=calculateSimilarity(img,imgRef);
          visual_correlation(i)=corrValue;
      end

      % for correlations, split into lines for each condition
      nComparisons=numel(visual_correlation);
      sectionLength=nComparisons/nConds;
      O0V0xVisualr=visual_correlation(1:sectionLength);
      O0V90xVisualr=visual_correlation(sectionLength + 1 : sectionLength*2);
      O90V0xVisualr=visual_correlation(sectionLength*2 + 1 : sectionLength*3);
      O90V90xVisualr=visual_correlation(sectionLength*3 + 1 : sectionLength*4);
      V0xVisualr=visual_correlation(sectionLength*4 + 1 : sectionLength*5);
      V90xVisualr=visual_correlation(sectionLength*5 + 1 : sectionLength*6);

      stimCon=unique(TS.Header.Conditions.StimCon);
      stimCon=stimCon(stimCon>0);


      % --- Calculate correlation & intensity between each condition's average image to a reference PCA baseline (max
      % contrast)  --- 
      for i=1:length(conditionSet)
          img=columnarResponse(:,:,conditionSet(i)-nblank).*ROImask;
          imgRef=pca_difference.*ROImask;%columnarBitmap(:,:,1).*ROImask;
          [~,~,corrValue,~]=calculateSimilarity(img,imgRef);
          pca_correlation(i)=corrValue;
      end

      % split into lines for each condition
      O0V0xPCAr=pca_correlation(1:sectionLength);
      O0V90xPCAr=pca_correlation(sectionLength + 1 : sectionLength*2);
      O90V0xPCAr=pca_correlation(sectionLength*2 + 1 : sectionLength*3);
      O90V90xPCAr=pca_correlation(sectionLength*3 + 1 : sectionLength*4);
      V0xPCAr=pca_correlation(sectionLength*4 + 1 : sectionLength*5);
      V90xPCAr=pca_correlation(sectionLength*5 + 1 : sectionLength*6);


      %% Plotting correlation of responses to visual and PCA difference reference map
      % --- Visual reference map ---
      figure('name',['Correlation: Responses vs visual & PCA reference ' trialsOutcomeLabel])
      % Plot versus visual 0-90°
      nCols=2;
      nRows=1;
      [hAx,~]=tight_subplot(nRows,nCols,[.1 .1]);
      axes(hAx(1))

      % Plot baselines: visual 0 and 90, without opto
      plot(stimCon,V0xVisualr,'b-','LineWidth',2,'Marker','o','MarkerSize',10,'MarkerFaceColor','w','MarkerEdgeColor','k','DisplayName','Vis-0'); hold on;
      plot(stimCon,V90xVisualr,'r-','LineWidth',2,'Marker','o','MarkerSize',10,'MarkerFaceColor','w','MarkerEdgeColor','k','DisplayName','Vis-90'); hold on;

      % Plot congruent: visual 0 and 90, with congruent opto
      plot(stimCon,O0V0xVisualr,'b--','LineWidth',2,'Marker','v','MarkerSize',10,'MarkerFaceColor','b','MarkerEdgeColor','k','DisplayName','Opto-0 + Vis-0'); hold on;
      plot(stimCon,O90V90xVisualr,'r--','LineWidth',2,'Marker','^','MarkerSize',10,'MarkerFaceColor','r','MarkerEdgeColor','k','DisplayName','Opto-90 + Vis-90'); hold on;

      % Plot incongruent: visual 0 and 90, with incongruent opto
      plot(stimCon,O90V0xVisualr,'b:','LineWidth',2,'Marker','^','MarkerSize',10,'MarkerFaceColor','r','MarkerEdgeColor','k','DisplayName','Opto-90 + Vis-0'); hold on;
      plot(stimCon,O0V90xVisualr,'r:','LineWidth',2,'Marker','v','MarkerSize',10,'MarkerFaceColor','b','MarkerEdgeColor','k','DisplayName','Opto-0 + Vis-90'); hold on;

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
      if saveFlag
          export_fig(pdfFilename,'-pdf','-nocrop','-append');
      end

      %% Extract fitted value at low contrast
      % --- Versus max visual response ---
      % Fit line, get 10% correlation value per condition
      fitExponent=1; %linear
      desired_xvalue=[10:1:15]; % average over 10-15% contrast
      [O0V0corr]=mean(getYfit(stimCon,O0V0xVisualr,fitExponent,desired_xvalue)); % opto 0 con
      [O0V90corr]=mean(getYfit(stimCon,O0V90xVisualr,fitExponent,desired_xvalue)); % opto 0 incon
      [O90V0corr]=mean(getYfit(stimCon,O90V0xVisualr,fitExponent,desired_xvalue)); % opto 90 incon
      [O90V90corr]=mean(getYfit(stimCon,O90V90xVisualr,fitExponent,desired_xvalue)); % opto 90 con
      visualr=[O0V90corr,O90V0corr]; % incon opto (H, V)
      visualrRatios=[O0V90corr./O0V0corr,O90V0corr./O90V90corr]; % incon opto / con opto (max possible during biasing)    

      % --- Versus PCA ---
      % Fit line, get 10% correlation value per condition
      [O0V0corrPCA]=mean(getYfit(stimCon,O0V0xPCAr,fitExponent,desired_xvalue)); % opto 0 con
      [O0V90corrPCA]=mean(getYfit(stimCon,O0V90xPCAr,fitExponent,desired_xvalue)); % opto 0 incon
      [O90V0corrPCA]=mean(getYfit(stimCon,O90V0xPCAr,fitExponent,desired_xvalue)); % opto 90 incon
      [O90V90corrPCA]=mean(getYfit(stimCon,O90V90xPCAr,fitExponent,desired_xvalue)); % opto 90 con
      PCAr=[O0V90corrPCA,O90V0corrPCA]; % incon opto (H, V)
      PCArRatios=[O0V90corrPCA./O0V0corrPCA, O90V0corrPCA./O90V90corrPCA]; % incon opto / con opto (max possible during biasing)

      %to reconstitute con opto O0V0corrPCA O90V90corrPCA
      %PCArFull=1/(PCArRatios/PCAr);


      %% Plotting intensity of responses to visual and PCA difference reference map

      % Intensity of opto 90 vs visual 0
      for i=1:length(V0O90)
          img=columnarResponse(:,:,V0O90(i)-nblank).*ROImask;
          imgRef=columnarResponse(:,:,V0(i)-nblank).*ROImask;
          [~,~,~,intensityPrct]=calculateSimilarity(img,imgRef);
          opto90IntensityFull(i)=intensityPrct;
      end

      % Intensity of opto 0 vs visual 90
      for i=1:length(V90O0)
          img=columnarResponse(:,:,V90O0(1)-nblank).*ROImask;
          imgRef=columnarResponse(:,:,V90(i)-nblank).*ROImask;
          [~,~,~,intensityPrct]=calculateSimilarity(img,imgRef);
          opto0IntensityFull(i)=intensityPrct;
      end



      %% Plot intensity difference for full, peaks and dips for O90 effect versus V0
      figure('name',['Opto 90 intensity ' trialsOutcomeLabel])
      % Plot versus visual 0-90°
      nCols=5;
      nRows=numel(V0O90);    
      [hAx,~]=tight_subplot(nRows,nCols,[.045 .045]);

      % get cortical response, and reference
      imgOpto=bitmap_difference;%columnarResponse(:,:,V0O90(1)-nblank).*ROImaskNaN;

      % --- Matched visual, incongruent opto vs visual ---
      % Intensity of opto 90 vs visual 0
      for i=1:length(V0O90)
          % get cortical response, and reference
          img=columnarResponse(:,:,V0O90(i)-nblank).*ROImaskNaN;
          imgRef=columnarResponse(:,:,V0(end)-nblank).*ROImaskNaN;

          % get image as percentage of reference
          imgPrct=100*(img - imgRef) ./ imgRef;        


          % get peaks and dips as percentage of reference
          imgPeaks=img; % Peaks
          imgRefPeaks=imgRef;
          imgPeaksMask=imgOpto>0;
          imgPeaks(~imgPeaksMask)=NaN;
          imgRefPeaks(~imgPeaksMask)=NaN;
          imgPeaksPrct=100*(imgPeaks - imgRefPeaks) ./ imgRefPeaks;        

          imgDips=img; % Dips
          imgRefDips=imgRef;
          imgDipsMask=imgOpto<0;
          imgDips(~imgDipsMask)=NaN;
          imgRefDips(~imgDipsMask)=NaN;
          imgDipsPrct=100*(imgDips - imgRefDips) ./ imgRefDips;

          axes(hAx((i-1)*5 + 1))
          imgsc(img,'Image');colorbar;caxis([-.01 .01])

          axes(hAx((i-1)*5 + 2))
          imgsc(imgRef, 'Reference');colorbar;caxis([-.01 .01])

          axes(hAx((i-1)*5 + 3))
          imgsc(rescale(imgPrct), ['\Delta_{Full}: ', sprintf('%.0f%%',nanmean(imgPrct(:)))]);colorbar

          axes(hAx((i-1)*5 + 4))
          imgsc(rescale(imgPeaksPrct), ['\Delta_{Peaks}: ', sprintf('%.0f%%',nanmean(imgPeaksPrct(:)))]);colorbar

          axes(hAx((i-1)*5 + 5))
          imgsc(rescale(imgDipsPrct), ['\Delta_{Dips}: ', sprintf('%.0f%%',nanmean(imgDipsPrct(:)))]);colorbar

          % save value
          %[~,~,~,intensityPrct]=calculateSimilarity(img,imgRef);
          opto90IntensityFull(i)=nanmean(imgPrct(:));
          opto90IntensityPeaks(i)=nanmean(imgPeaksPrct(:));
      end

     1; 

      %% Plot intensity difference for full, peaks and dips for O0 effect versus V90
      figure('name',['Opto 0 intensity ' trialsOutcomeLabel])
      % Plot versus visual 0-90°
      nCols=5;
      nRows=numel(V0O90); 
      [hAx,~]=tight_subplot(nRows,nCols,[.045 .045]);

      % get cortical response, and reference
      imgOpto=bitmap_difference;

      % --- Matched visual, incongruent opto vs visual ---
      % Intensity of opto 90 vs visual 0
      for i=1:length(V90O0)
          % get cortical response, and reference
          img=columnarResponse(:,:,V90O0(i)-nblank).*ROImaskNaN;
          imgRef=columnarResponse(:,:,V90(end)-nblank).*ROImaskNaN;

          % get image as percentage of reference
          imgPrct=100*(img - imgRef) ./ imgRef;        


          % get peaks and dips as percentage of reference
          imgPeaks=img; % Peaks
          imgRefPeaks=imgRef;
          imgPeaksMask=imgOpto<0;
          imgPeaks(~imgPeaksMask)=NaN;
          imgRefPeaks(~imgPeaksMask)=NaN;
          imgPeaksPrct=100*(imgPeaks - imgRefPeaks) ./ imgRefPeaks;        

          imgDips=img; % Dips
          imgRefDips=imgRef;
          imgDipsMask=imgOpto>0;
          imgDips(~imgDipsMask)=NaN;
          imgRefDips(~imgDipsMask)=NaN;
          imgDipsPrct=100*(imgDips - imgRefDips) ./ imgRefDips;

          axes(hAx((i-1)*5 + 1))
          imgsc(img,'Image');colorbar;caxis([-.01 .01])

          axes(hAx((i-1)*5 + 2))
          imgsc(imgRef, 'Reference');colorbar;caxis([-.01 .01])

          axes(hAx((i-1)*5 + 3))
          imgsc(rescale(imgPrct), ['\Delta_{Full}: ', sprintf('%.0f%%',nanmean(imgPrct(:)))]);colorbar

          axes(hAx((i-1)*5 + 4))
          imgsc(rescale(imgPeaksPrct), ['\Delta_{Peaks}: ', sprintf('%.0f%%',nanmean(imgPeaksPrct(:)))]);colorbar

          axes(hAx((i-1)*5 + 5))
          imgsc(rescale(imgDipsPrct), ['\Delta_{Dips}: ', sprintf('%.0f%%',nanmean(imgDipsPrct(:)))]);colorbar

          % save value
          %[~,~,~,intensityPrct]=calculateSimilarity(img,imgRef);
          opto0IntensityFull(i)=nanmean(imgPrct(:));
          opto0IntensityPeaks(i)=nanmean(imgPeaksPrct(:));
      end

     1; 
      %{
      % Strength of opto 0 vs visual 90
      for i=1:length(V90O0)
          img=columnarResponse(:,:,V90O0(1)-nblank).*ROImask;
          imgRef=columnarResponse(:,:,V90(i)-nblank).*ROImask;
          [~,~,~,intensityPrct]=calculateSimilarity(img,imgRef);
          opto0StrengthFull(i)=intensityPrct;
      end
     %}



      opto0IntensityFull=rmoutliers(opto0IntensityFull);
      opto0IntensityPeaks=rmoutliers(opto0IntensityPeaks);
      opto90IntensityFull=rmoutliers(opto90IntensityFull);
      opto90IntensityPeaks=rmoutliers(opto90IntensityPeaks);


      %% Opto intensity ratios
      optoIntensity=[mean(opto0IntensityFull),mean(opto90IntensityFull), mean(opto0IntensityPeaks),mean(opto90IntensityPeaks)];%[opto0BiasingStrength,opto90BiasingStrength];
    end
else
    visualr=nan(1,2);
    visualrRatios=nan(1,2);
    PCAr=nan(1,2);
    PCArRatios=nan(1,2);
    optoIntensity=nan(1,4);
end
end