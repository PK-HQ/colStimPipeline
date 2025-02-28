%% Define reference and current session
pipelineMode='';% beta/stable
analysisType='psychfit'; %summary/psychfit/beta/neurometric/stability/neurometricsfix

%for biasing expt
currentSessID=112; %81;
%for analysis (demo:[79 78 81])
%analysisBlockID=[120:122 116:118 112:114 108:109 89 84 83 68 58];%[106 102 98]; %sort([109 108* 105* 104 101 100 93 92 90 89 88 87 84 83 81 79 78 74 73 72 67 66 65 64 62 59 76 77]);%sort([93 92 90 89 88 87 84 83 81 79 78 74 73 72 67 66 65 64 62 59 76 77]); % 51 50 %[12 14 15 17 19 20 21 23 24 26 29 31 32 34 37 40]; %10 excluded because baseline different

%% Load dataStruct
if ispc
  mainPath='Y:/';
elseif contains(getenv('HOSTNAME'),'psy.utexas.edu')
  mainPath='/eslab/data/';
end
run([mainPath 'users/PK/colStimPipeline/exptListBiasingFull.m'])
nColumnsWanted=[];
analysisBlockID = organizeBlocks(datastruct, 'L', nColumnsWanted);
datastructChamberLeft=datastruct(analysisBlockID);
%analysisBlockID=setdiff(datastructChamberRight,10)%setdiff(datastructChamberLeft,[44 45 45 52:55 95]);
set(0,'DefaultFigureWindowStyle','docked')

%% Run the desired bitmap generation pipeline variant
switch pipelineMode
    case {'beta'}
        saveFlag=0;
        saveFlagBMP=1;
        plotFlag=1;
        
        %% Fast pipeline (coregisters within session ort map to optostim block green image)
        % Get correction for camera-projector alignment
        % imregtform of orignal and recovered bitmap, apply to correct for
        % camera-projector alignment
        open projectorCameraCalibration
        
        % Get projector-camera alignment transformation matrix
        load([alignmentSessionFolder 'alignmentTransform.mat'],'alignmentTransform')

        % Define sessions
        dsReferenceSess=datastruct(referenceSessID); %moving, e.g. ort map
        dsCurrentSess=datastruct(currentSessID); %fixed, to get ort map projected onto
        filenameStructCurrent=generateFilenames(dsCurrentSess);
        pdfFilename=filenameStructCurrent.neurometricPDF;

        % Get within session reference orientation map
        [columnarBitmap,VERpca,columnarmapStats]=getColumnarBitmapV4(mainPath,dsReferenceSess,dsCurrentSess,bitmapParams, ...
          plotFlag, saveFlag, pdfFilename);

        % Coregister to most current green image
        [columnarBitmapCoregistered, columnarPCAsCoregistered]=coregisterBitmap2GreenImgV2(dsReferenceSess,dsCurrentSess, ...
          columnarBitmap,VERpca,plotFlag,saveFlag, pdfFilename);

        % Account for PRF
        %?
        
        % Correct for camera-projector alignment
        orts=[0 90];%0:15:165;
        HE=1000;
        
        % Define sessions
        dsReferenceSess=datastruct(referenceSessID); %moving, e.g. ort map
        dsCurrentSess=datastruct(currentSessID); %fixed, to get ort map projected onto
        
        [projBitmapTRBB,bitmapsCamSpace,nBlobs,medianBlobAreas]=convertForProjector(dsReferenceSess,dsCurrentSess,columnarBitmapCoregistered,orts,...
            bitmapParams.gridSize,bitmapParams.gammaCorrFactor,bitmapParams.sensitivity,'cam2proj',alignmentTransform,plotFlag,saveFlagBMP,saveFlag, pdfFilename);
    1;
  case {'stable'} %within session
        saveFlag=0;
        saveFlagBMP=1;
        plotFlag=1;
        
        %% Fast pipeline (coregisters within session ort map to optostim block green image)
        % Get correction for camera-projector alignment
        % imregtform of orignal and recovered bitmap, apply to correct for
        % camera-projector alignment
        open projectorCameraCalibration
        
        % Get projector-camera alignment transformation matrix
        load([alignmentSessionFolder 'alignmentTransform.mat'],'alignmentTransform')

        % Define sessions
        dsReferenceSess=datastruct(referenceSessID); %moving, e.g. ort map
        dsCurrentSess=datastruct(currentSessID); %fixed, to get ort map projected onto
        
        % Get within session reference orientation map
        [columnarBitmap,VERpca,columnarmapStats]=getColumnarBitmapV4(mainPath,dsReferenceSess,dsCurrentSess,bitmapParams, ...
          plotFlag, saveFlag, pdfFilename);

        % Coregister to most current green image
        [columnarBitmapCoregistered, columnarPCAsCoregistered]=coregisterBitmap2GreenImgV2(dsReferenceSess,dsCurrentSess, ...
          columnarBitmap,VERpca,plotFlag,saveFlag, pdfFilename);

        % Correct for camera-projector alignment
        orts=[0 90];%0:15:165;
        HE=1000;
        [projBitmapTRBB,bitmapsCamSpace,nBlobs,medianBlobAreas]=convertForProjector(dsReferenceSess,dsCurrentSess,columnarBitmapCoregistered,orts,...
            bitmapParams.gridSize,bitmapParams.gammaCorrFactor,bitmapParams.sensitivity,'cam2proj',alignmentTransform,plotFlag,saveFlagBMP,saveFlag, pdfFilename);
end

%% Run the desired analysis pipeline variant
switch analysisType
    case {'summary2024'}
        if ~exist('behavioralData','var')
            behavioralData=[];
            imagingData=[];
            bitmapData=[];
        end
        
        for blockID=1:numel(analysisBlockID)
            disp(['=== Block ' num2str(blockID)  '/' num2str(numel(analysisBlockID)) '==='])
            tic
            saveFlag=1;
            saveFlagBMP=0; 
            plotFlag=1;
            if isfield(behavioralData,'auc') && size(behavioralData.auc,3)>=blockID
                 disp('(Skipping completed block)')
                 continue
            end
            % === Load block data === (CHANGE SO THAT THIS WILL SAVE PAST DATA)
             [currentBlockStruct,referenceBlockStruct,...
                behavioralData, imagingData, bitmapData, successFlag]=loadBlockData(datastruct, analysisBlockID, behavioralData, imagingData, bitmapData, blockID);
             if ~successFlag
                 continue
             end
            % === Select pdf save file name ===
             pdfFilename=currentBlockStruct.psychneuroPDF;
    
              % ===  Get orientation map ===
              bitmapData=getColumnarBitmapV4(currentBlockStruct, imagingData, bitmapData, blockID, ...
                pdfFilename, plotFlag, saveFlag);
    
              % === Transform columnar positions to current cortical view ===
              bitmapData=coregisterBitmap2GreenImgV2(currentBlockStruct,referenceBlockStruct, ...
                imagingData,bitmapData,blockID, ...
                pdfFilename, plotFlag,saveFlag);
    
              % === Generate bitmap, correct for projector properties and camera-projector alignment===
              orts=[0 90];%0:15:165;
              HE=1000;

              %{
               [bitmapData]=convertForProjectorGPT(behavioralData, imagingData, bitmapData,...
                    currentBlockStruct,'cam2proj', blockID, ...
                    pdfFilename, plotFlag,saveFlagBMP,saveFlag)
              %}
               [bitmapData]=convertForProjectorGPT(behavioralData, imagingData, bitmapData,...
                    currentBlockStruct,'proj2cam', blockID, ...
                    pdfFilename, plotFlag,saveFlagBMP,saveFlag);
    
                % === Plot behavioral biasing results ===
                reportType='summary';
                behavioralData=analyzeBlockPsychometrics(currentBlockStruct, behavioralData, bitmapData, blockID,...
                    pdfFilename, reportType, saveFlag);
                toc
                CF;
                imagingData.intg=[];behavioralData.TS=[];
                
                %save('Y:/Chip/Meta/summary/statisticsR.mat','bitmapData','behavioralData','imagingData','analysisBlockID', '-v7.3')
                
                %pause(5)
        end
        1;
        %save('Y:/Chip/Meta/summary/statistics.mat','bitmapData','behavioralData','imagingData','analysisBlockID')
        
        %save('Y:/Chip/Meta/summary/statistics73.mat','bitmapData','behavioralData','imagingData','analysisBlockID', '-v7.3')
        
       % load('Y:\Chip\Meta\summary\statisticsHPC.mat')

        %% Psychometrics
        % === Session summary ===
        analyzeSessionPsychometrics(behavioralData, bitmapData, datastruct, analysisBlockID, saveFlag)
        %export_fig('Y:\users\PK\Eyal\meetings\summary\summarypsychometrics1to40.pdf','-pdf','-nocrop');

        % === Power series ===
        nColumns=20;
        plotPowerSeries(bitmapData, behavioralData, nColumns);
        %export_fig('Y:/Chip/Meta/powerSeries/powerSeries.pdf','-pdf','-nocrop');

        % === Min columns ===
        plotMinColumns(bitmapData, behavioralData, analysisBlockID)
        %export_fig('Y:/Chip/Meta/minColumnSeries/minColumns.pdf','-pdf','-nocrop');

        %% Modelling
        open npmodelParameterEval.m

    
    case {'psychometrics'} %% Psychfit
        psychometrics=[];
        for blockID=1:numel(analysisBlockID)%1:numel(analysisBlockID)%1:numel(analysisBlockID)
            currentSessID=analysisBlockID(blockID);
            dsCurrentSess=datastruct(currentSessID); %fixed, to get ort map projected onto
            dsReferenceSess=datastruct(dsCurrentSess.referenceBlockNo);
            alignmentSession=datastruct(dsCurrentSess.alignmentBlockNo).date;

            saveFlag=1;
            plotFlag=1;
            isSummary=1;

            % --- Behavioral biasing ---
            hAx=gca;
            filenameStructCurrent=generateFilenames(dsCurrentSess);
            reportType='optostim';
            [betaSorted, muSorted, sigmaSorted, contrastsSorted, percentVerticalSorted]=analyzeBlockPsychometrics(dsCurrentSess, filenameStructCurrent, reportType, saveFlag);


            if saveFlag
                if blockID==1
                    export_fig(filenameStructCurrent.psychfitPDF,'-pdf','-nocrop');
                elseif blockID>1
                    export_fig(filenameStructCurrent.psychfitPDF,'-pdf','-nocrop','-append');
                end
            end
            %[betaSorted,muSorted,sigmaSorted,contrastsSorted,percentVerticalSorted]=analyzePsychometrics(dsCurrentSess,hAx,saveFlag,'psychometrics');
            
            %{
            psychometrics(blockID).betaSorted=betaSorted;
            psychometrics(blockID).muSorted=muSorted;
            psychometrics(blockID).sigmaSorted=sigmaSorted;
            psychometrics(blockID).contrastsSorted=contrastsSorted;
            psychometrics(blockID).percentVerticalSorted=percentVerticalSorted;
            psychometrics(blockID).percentHorizontalSorted=100-percentVerticalSorted;
            
            % calculate correct
            horzCond=find(contrastsSorted<0);
            percentVerticalSorted(horzCond)=100-percentVerticalSorted(horzCond);
            psychometrics(blockID).percentCorrectSorted=percentVerticalSorted;
            %nColumnsCont(blockID,:)=dsCurrentSess.nColumns;
            %}
        end
        % session summary
        sessionSEM
        
        % stats
        ranksum(betaSortedCont(:,3),betaSortedCont(:,1),'tail','both')
        
        % plot minimum columns
        hBeta=betaSortedCont(:,1);
        vBeta=betaSortedCont(:,3);
        deltaBeta=round(vBeta-hBeta,2);
        figure('name','Minimum columns for biasing')
        scatter(nColumnsCont,deltaBeta,100,'ok','linewidth',2)
        xticks([1 5:5:40])
        xlim([0 40])
        upFontSize(13,0.005)
        1;
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    
    
    case {'summary'}        
        if ~exist('behavioralData','var')
            behavioralData=[];
            imagingData=[];
            bitmapData=[];
        end
        for blockID=find(analysisBlockID==76)%1:numel(analysisBlockID)%1:numel(analysisBlockID)%1:numel(analysisBlockID)%1:numel(analysisBlockID)
            saveFlag=0;
            saveFlagBMP=0;
            plotFlag=1;

        % === Load block data ===
         [currentBlockStruct,referenceBlockStruct,...
            behavioralData, imagingData, bitmapData]=loadBlockData(datastruct, analysisBlockID, behavioralData, imagingData, bitmapData, blockID);
       
         pdfFilename=currentBlockStruct.psychneuroPDF;

          % ===  Get orientation map ===
          bitmapData=getColumnarBitmapV4(currentBlockStruct, imagingData, bitmapData, blockID, ...
            pdfFilename, plotFlag, saveFlag);

          % === Transform columnar positions to current cortical view ===
          bitmapData=coregisterBitmap2GreenImgV2(currentBlockStruct,referenceBlockStruct, ...
            imagingData,bitmapData,blockID,pdfFilename,plotFlag,saveFlag);

          % === Generate bitmap, correct for projector properties and camera-projector alignment===
          orts=[0 90];%0:15:165;
          HE=1000;
          [bitmapData]=convertForProjector(behavioralData,imagingData,bitmapData,...
              currentBlockStruct,'cam2proj',blockID,pdfFilename,plotFlag,saveFlagBMP,saveFlag);

            % === Plot behavioral biasing results ===
            hAx=gca;
            reportType='optostim';
            [behavioralData]=...
                analyzeBlockPsychometrics(currentBlockStruct, behavioralData, bitmapData, blockID, pdfFilename, reportType, saveFlag);
            if saveFlag
                if blockID==1
                    export_fig(filenameStructCurrent.psychfitPDF,'-pdf','-nocrop');
                elseif blockID>1
                    export_fig(filenameStructCurrent.psychfitPDF,'-pdf','-nocrop','-append');
                end
            end

            % === Storing results ===
            bitmapData.nColumns
            bitmapData.adjustedSPD_uW
            behavioralData.auc(2:3,3)' %behavioralData.auc(2:3,3)-behavioralData.auc(1,3)'
            
            % --- Power series ---
            plot(bitmapData.adjustedSPD_uW,behavioralData.auc(2:3,3)')


            %% Old
            % Static images
            %[blurEstimate,gfpStaticMax,mcherryStaticMax]=QCstatic(dsCurrentSess,plotFlag,saveFlag);
            
            %{
            % slope and intercept
            figure('name','Psychoneurometric plot')
            nRows=1;nCols=3;
            gap=.05;marginV=.01;marginH=.03;
            [hAx,~]=tight_subplot(nRows,nCols,[gap gap], [marginV+.1 marginV+.2], [marginH+.05 marginH]);
            % --- Behavioral biasing ---
            [betaSorted,muSorted,sigmaSorted]=analyzeBiasingBlock(dsCurrentSess,hAx(1),0,'full'); upFontSize(32,.01);
            %}
            % --- Neurometric biasing (integrated responses) --- 
                if isfile(filenameStructCurrent.gaussianFit) % load newly fit one
                    load(filenameStructCurrent.gaussianFit) % .gaussian file
                else
                    ROIMaskgaussian=[];
                end
            phiFull=analyzeProjection1D(dsCurrentSess, TSstruct, DataTrialStruct, Mask, ROIMaskgaussian, bitmapsCamSpace, RespCondPCA, hAx, 2, 'full', saveFlag); upFontSize(32,.01);
            colormap(inferno)
            % --- Neurometric biasing (integrated responses) --- 
            phiRecruit=analyzeProjection1D(dsCurrentSess, TSstruct, DataTrialStruct, Mask, ROIMaskgaussian, bitmapsCamSpace, RespCondPCA, hAx, 3, 'nonstim', saveFlag); upFontSize(32,.01);
            colormap(inferno);caxis([-.0025 0.0025]);caxis([-.0015 0.0015])
            1;
            %[visualr,visualrRatio,PCAr,PCArRatio,optoIntensity]=assessOptoEffect2(dsReferenceSess,dsCurrentSess,columnarBitmapCoregistered,columnarPCAsCoregistered,saveFlag);
            % --- Neurometric biasing (dynamic frame-wise responses) --- 
            %analyzeTrajectory(dsReferenceSess,dsCurrentSess)
            

           %% Measures (one per measure, per session)

           % --- Imaging quality ---
           %{
            % Green image blur
            %blurEstimateCont(blockID,1)=blurEstimate;
            % GFP static intensity
            gfpStaticMaxCont(blockID,1)=gfpStaticMax;
            % mCherry static intensity
            mcherryStaticMaxCont(blockID,1)=mcherryStaticMax;

            % --- Stimulation power and positioning ---
            % Orientation map quality
            %PCAExplTotalCont(blockID,1)=PCAExplTotal;
            % Pixel density = Percentage pixels on within ROI
            %pixelDensitiesCont(blockID,1:2)=nonzeros(columnarmapStats.pixDensity);
            % Power density = Percentage pixels on within ROI * power * ON-OFF
            %powerDensitiesCont(blockID,1:2)=powerDensities;
            % No. of blobs
            %nBlobsCont(blockID,1:2)=nonzeros(bitmapStats.nBlobs);
           %}
            % Size of blobs (median)
            %contourDC(blockID,1:2)=nonzeros(bitmapStats.medianBlobAreas);
            contourColumns(blockID,1:2)=nonzeros(bitmapStats.nBlobs);
            contourPixels(blockID,1:2)=bitmapStats.contourPix;
            contourPixelDensities(blockID,1:2)=bitmapStats.contourPixDensity;
            
            %pixDensityContour
            %estDutycycleContour
            %pixDensityTotal
            %estDutycycleTotal
            
            % Correlation with bitmap
            
            % --- Session ID ---
            blockIDCont(blockID)=currentSessID;
            
            % --- Behavioral biasing ---
            betaSortedCont(blockID,:)=betaSorted; % H-BL-V
            muSortedCont(blockID,:)=muSorted;
            sigmaSortedCont(blockID,:)=sigmaSorted;

            % --- Neurometric biasing ---
            % Opto similarity vs visual
            phiSortedCont(blockID,:)=phiFull; % H-BL-V
            phiRecruitSortedCont(blockID,:)=phiRecruit; % H-BL-V

            % Opto intensity vs visual
            %optoIntensityCont(blockID,1:4)=optoIntensity; %2x overall mean, 2x peak mean         

            CF
            toc
        end

        getTemporalStimParams

        save([mainPath 'Chip/Meta/summary/neurometricsMinimalColumns' num2str(numel(analysisBlockID)) 'X08112023.mat'],...
            'analysisBlockID',...
            'betaSortedCont','muSortedCont','sigmaSortedCont',...
            'phiSortedCont','phiRecruitSortedCont',...
            'temporalTotal','temporalDC','temporalON')

        save([mainPath 'Chip/Meta/summary/neurometricsMinimalColumns' num2str(numel(analysisBlockID)) 'X08112023.mat'],...
            'analysisBlockID',...
            'contourColumns','contourPixels','contourPixelDensities',...
            'betaSortedCont','muSortedCont','sigmaSortedCont',...
            'phiSortedCont','phiRecruitSortedCont',...
            'temporalTotal','temporalDC','temporalON','tbl')

        tbl=array2table([analysisBlockID' contourColumns (contourPixels.*temporalON') betaSortedCont],...
            'VariableNames',{'EntryID' 'Column-H' 'Column-V', 'Energy-H', 'Energy-V','Beta-H','Beta-Baseline','Beta-V'});
        
        load('Y:\Chip\Meta\summary\neurometricsMinimalColumns25tbl.mat')
        %% P-values
        % Define x y z, averaged across 2 patterns       
        zPower=contourPixels.*temporalON'; 
        xColumns=contourColumns;
        yBeta=betaSortedCont;
        
        outlierBlocks=isoutlier(zPower); % 
        zPower(outlierBlocks)=NaN;
        xColumns(outlierBlocks)=NaN;
        yBeta(outlierBlocks)=NaN;
        [p_hbl,~]=ranksum(yBeta(:,1),yBeta(:,2),'tail','both')
        [p_vbl,~]=ranksum(yBeta(:,3),yBeta(:,2),'tail','both')
        [p_hv,~]=ranksum(yBeta(:,1),yBeta(:,3),'tail','both')
        
        
        
        %% Bias x time
        figure('name','Psychoneurometric plot')
        nRows=1;nCols=3;
        gap=.05;marginV=.01;marginH=.03;
        [hAx,~]=tight_subplot(nRows,nCols,[gap gap], [marginV+.1 marginV+.2], [marginH+.05 marginH]);
        axes(hAx(1))
        yline(0,'--','Color',[.3 .3 .3],'LineWidth',1.5,'HandleVisibility','off'); hold on;
        x=1:numel(analysisBlockID);
        y=betaSortedCont(:,1)-betaSortedCont(:,2);
        scatter(x,y,350,'b','filled','v','MarkerFaceColor','b','MarkerEdgeColor','k','LineWidth',2); hold on;
        x=1:numel(analysisBlockID);
        y=betaSortedCont(:,3)-betaSortedCont(:,2);
        scatter(x,y,350,'r','filled','^','MarkerFaceColor','r','MarkerEdgeColor','k','LineWidth',2); hold on;
        legend({'H-optostim','V-optostim'})
        xlabel('Session');addSkippedTicks(0,roundup(numel(analysisBlockID),4),2,'x'); xlim([0,roundup(numel(analysisBlockID),4)])
        ylabel('\beta_{Optostim-Baseline}');addSkippedTicks(-.7,.7,.1,'y'); ylim([-.7 .7])
        legend({'H-optostim','V-optostim'})
        title('Psychometrics across sessions','FontWeight','Normal')
        upFontSize(28,.008); axis square
        
        %% Phi x time
        axes(hAx(2))
        yline(0,'--','Color',[.3 .3 .3],'LineWidth',1.5,'HandleVisibility','off'); hold on;
        x=1:numel(analysisBlockID);
        y=phiSortedCont(:,2)-phiSortedCont(:,1);
        scatter(x,y,250,'b','filled','v','MarkerFaceColor','b','MarkerEdgeColor','k','LineWidth',2); hold on;
        x=1:numel(analysisBlockID);
        y=phiSortedCont(:,3)-phiSortedCont(:,1);
        scatter(x,y,250,'r','filled','^','MarkerFaceColor','r','MarkerEdgeColor','k','LineWidth',2); hold on;
        legend({'H-optostim','V-optostim'})
        xlabel('Session');addSkippedTicks(0,roundup(numel(analysisBlockID),4),2,'x'); xlim([0,roundup(numel(analysisBlockID),4)])
        ylabel('\phi_{Optostim-Baseline}');addSkippedTicks(-3,3,.5,'y');ylim([-3 3])
        legend({'H-optostim','V-optostim'})
        title('Neurometrics across sessions','FontWeight','Normal')
        upFontSize(28,.008); axis square

       %% Minimal columns
        % Define x y z, averaged across 2 patterns       
        zPower=contourPixels.*temporalON'; 
        xColumns=contourColumns;
        yBeta=betaSortedCont(:,[1 3]);
        
        % 2D line fit (ncolumns X beta)
        metricStr='beta';
        yBeta=betaSortedCont(:,[1 3]);
        plotMinimumColumns(xColumns,yBeta,zPower,analysisBlockID,metricStr);upFontSize(32,.01);
        title('No. of columns versus \beta','FontWeight','Normal');ylabel('\beta_{optostim}')
         
        metricStr='deltabeta';
        yBeta=betaSortedCont(:,[1 3]);
        plotMinimumColumns(xColumns,yBeta,zPower,analysisBlockID,metricStr);upFontSize(32,.01);
        title('No. of columns versus \Delta\beta','FontWeight','Normal');
         
        metricStr='deltaphi';
        yBeta=phiSortedCont(:,[2 3]);
        plotMinimumColumns(xColumns,yBeta,zPower,analysisBlockID,metricStr);upFontSize(32,.01);
        title('No. of columns versus ϕ (full ROI)','FontWeight','Normal'); ylabel('\Delta\phi'); ylim
        
        metricStr='deltaphi';
        yBeta=phiRecruitSortedCont(:,[2 3]);
        plotMinimumColumns(xColumns,yBeta,zPower,analysisBlockID,metricStr);upFontSize(32,.01);
        title('No. of columns versus ϕ (recruitment ROI)','FontWeight','Normal'); ylabel('\Delta\phi'); ylim([-.1 .5]); addSkippedTicks(-.2,.5,.1,'y')

        metricStr='deltabetaphi';
        yBeta=betaSortedCont(:,[1 3]);
        xPhiRecruit=phiRecruitSortedCont(:,[2 3]);
        plotMinimumColumns(xPhiRecruit,yBeta,zPower,analysisBlockID,metricStr);upFontSize(32,.01);
        title('ɸ versus β (recruitment ROI)','FontWeight','Normal'); xlabel('\Delta\phi'), ylabel('\Delta\beta'); %ylim([-.1 .5]); addSkippedTicks(-.2,.5,.1,'y')
        xlim([-.05 .5])
        ylim([0 .7])
        addSkippedTicks(-.2,.6,.1,'x')
        addSkippedTicks(0,.8,.1,'y')

        metricStr='deltabetaphi';
        yBeta=betaSortedCont(:,[1 3]);
        xPhi=phiSortedCont(:,[2 3]);
        plotMinimumColumns(xPhi,yBeta,zPower,analysisBlockID,metricStr);upFontSize(32,.01);
        title('ɸ versus β (full ROI)','FontWeight','Normal'); xlabel('\Delta\phi'), ylabel('\Delta\beta'); %ylim([-.1 .5]); addSkippedTicks(-.2,.5,.1,'y')

        

        
        
        % 3D surface fit (power x ncolumns X beta)
        plotMinimumColumns3d(xColumns,yBeta,zPower,analysisBlockID)

        
        zPower(controlBlocks)=NaN; %control
        outlierBlocks=zPower<70000 | zPower>100000; zPower(controlBlocks)=NaN; zPower(outlierBlocks)=NaN;
        zPower(isoutlier(zPower))=NaN; exclBlocks=isnan(zPower); zPower(exclBlocks)=[]; zPower(controlBlocks)=0;
        xColumns(exclBlocks)=NaN;
        yDeltaBeta(exclBlocks)=NaN;
        xColumns(isnan(xColumns))=[];yDeltaBeta(isnan(yDeltaBeta))=[];
        
        
        yline(0,'--','HandleVisibility','off'); hold on % y=0 line
        h=scatter(xColumns,yDeltaBeta,100,zPower,'filled','MarkerEdgeColor','k','LineWidth',2)%,'w','o','filled','MarkerEdgeColor','k','HandleVisibility','on','LineWidth',2);
        hold off
        yticks([-.5:.1:1]); ylim([-.2 1]);  
        xlim([0 40])
        upFontSize(14,.0025)
        title('Columns x biasing effect')
        xlabel('Average of columns for H- and V-optostim bitmap','FontName','Arial')
        ylabel('\Delta\beta_{V-H optostim}','FontName','Arial','FontSize',18)
        legend({'\Delta\beta_{V-H optostim}'})
        
        
        
        %% other params
        PCArConOpto=1./(PCArRatiosCont./PCArCont);       
        PCArInconOpto=PCArCont;
        PCArFullCont=[PCArConOpto(:,1) PCArInconOpto(:,:) PCArConOpto(:,2)];
        %save([mainPath 'Chip/Meta/summary/neurometricsMinimalColumns.mat'],'visualrCont','visualrRatiosCont','PCArCont','PCArFullCont','PCArRatiosCont','optoIntensityCont','blurEstimateCont','gfpStaticMaxCont','mcherryStaticMaxCont','PCAExplTotalCont','pixelDensitiesCont','powerDensitiesCont','nBlobsCont','medianBlobAreasCont','betaSortedCont','muSortedCont','sigmaSortedCont')
        %load('Y:\Chip\Meta\summary\neurometricsMinimalColumns25tbl.mat')
        1;
        
        %% Analyze parameters across sessions
        analyzeParameters
        
        %% Plot histogram
        plotBiasingHistogram(mainPath,psychometric,saveFlag)
        
        %% Beta
        x_baseline=betaSortedCont(:,2);
        y_hopto=betaSortedCont(:,1);
        y_vopto=betaSortedCont(:,3);
        
        figure
        h1=scatter(x_baseline,y_hopto,160,'b','v','filled','MarkerEdgeColor','k'); hold on
        h2=scatter(x_baseline,y_vopto,160,'r','^','filled','MarkerEdgeColor','k');
        hold off
        xticks([0:.25:1])
        yticks([0:.25:1])
        xlabel('Baseline (\beta)')
        ylabel('Optostim (\beta)')
        title('Optostim biasing wrt. baseline')
        xlim([0 1]);ylim([0 1])
        
        % diagonal reference
        hline=refline(1,0);
        hline.LineStyle='--';hline.Color=ones(1,3)*.6;hline.LineWidth=1;
        
        %stats
        [p_hv,~]=ranksum(y_hopto,y_vopto,'tail','both')
        [p_hbl,~]=ranksum(y_hopto,x_baseline,'tail','both')
        [p_vbl,~]=ranksum(y_vopto,x_baseline,'tail','both')
        upFontSize(24,0.012)
        1;

        
       %% Beta subtracted
        x_baseline=betaSortedCont(:,2);
        y_vhopto=betaSortedCont(:,3)-betaSortedCont(:,1);
        
        figure
        h1=scatter(x_baseline,y_vhopto,160,[138,43,226]/255,'diamond','filled','MarkerEdgeColor','k'); hold on
        hold off
        xticks([0:.25:1])
        yticks([0:.25:1])
        xlabel('Baseline (\beta)')
        ylabel('Optostim V-H (\beta)')
        title('Diff. in optostim biasing wrt. baseline')
        xlim([0 1]);ylim([-1 1])
        
        % diagonal reference
        hline=refline(0,0);vline=xline(0.5);
        hline.LineStyle='--';hline.Color=ones(1,3)*.6;hline.LineWidth=1;
        vline.LineStyle='--';vline.Color=ones(1,3)*.6;hline.LineWidth=1;
        
        %stats
        [p_vh,~]=ranksum(y_vhopto,zeros(size(y_vhopto)),'tail','both')
        upFontSize(24,0.012)
        1;
        
       %% mu
        x_baseline=muSortedCont(:,2);
        y_hopto=muSortedCont(:,1);
        y_vopto=muSortedCont(:,3);
        
        figure
        h1=scatter(x_baseline,y_hopto,160,'b','o','filled','MarkerEdgeColor','k'); hold on
        h2=scatter(x_baseline,y_vopto,160,'r','o','filled','MarkerEdgeColor','k');
        hold off
        xticks([0:10:40])
        yticks([0:10:40])
        xlabel('Baseline (\mu)')
        ylabel('Optostim (\mu)')
        title('Optostim masking wrt. baseline')
        xlim([0 40]);ylim([0 40])
       
        % diagonal reference
        hline=refline(1,0);
        hline.LineStyle='--';hline.Color=ones(1,3)*.6;hline.LineWidth=1;
        upFontSize(24,0.012)
        
        %stats
        [p_hv,~]=ranksum(y_hopto,y_vopto,'tail','both')
        [p_hbl,~]=ranksum(y_hopto,x_baseline,'tail','both')
        [p_vbl,~]=ranksum(y_vopto,x_baseline,'tail','both')
        1;        
        
        %% sigma
        x_baseline=sigmaSortedCont(:,2);
        y_hopto=sigmaSortedCont(:,1);
        y_vopto=sigmaSortedCont(:,3);
        
        figure
        h1=scatter(x_baseline,y_hopto,160,'b','.','filled','MarkerEdgeColor','b'); hold on
        h2=scatter(x_baseline,y_vopto,160,'r','.','filled','MarkerEdgeColor','r');
        hold off
        xticks([0:1:4])
        yticks([0:1:4])
        xlabel('Baseline (\sigma)')
        ylabel('Optostim (\sigma)')
        title('Optostim masking wrt. baseline')
        xlim([0 4]);ylim([0 4])
       
        % diagonal reference
        hline=refline(1,0);
        hline.LineStyle='--';hline.Color=ones(1,3)*.6;hline.LineWidth=1;
        
        %stats
        [p_hv,~]=ranksum(y_hopto,y_vopto,'tail','both')
        [p_hbl,~]=ranksum(y_hopto,x_baseline,'tail','both')
        [p_vbl,~]=ranksum(y_vopto,x_baseline,'tail','both')
        1;        
        
        upFontSize(24,0.012)
        1;        
        
        %Robust fit
        plotRobustFit(betaSortedCont(:,1)-betaSortedCont(:,3),optoIntensityCont(:,1)-optoIntensityCont(:,2))
        
        %% Do all cross-session statistics here
        %plotBiasingHistogram(betaSortedCont)  
    case {'neurometricsfix'}
        for blockID=1:numel(analysisBlockID)%1:numel(analysisBlockID)%1:numel(analysisBlockID)%1:numel(analysisBlockID)
            tic
            % define session IDs, session meta info
            currentSessID=analysisBlockID(blockID);
            dsCurrentSess=datastruct(currentSessID); %fixed, to get ort map projected onto
            dsReferenceSess=datastruct(dsCurrentSess.referenceBlockNo);
            alignmentSession=datastruct(dsCurrentSess.alignmentBlockNo).date;
            
            bitmapParams.gridSize=datastruct(currentSessID).gridSize; %smaller = more rubbish; splits image into NumTiles=512x512/gridSize for adaptive histeq
            bitmapParams.gammaCorrFactor=datastruct(currentSessID).gammaCorrFactor; %smaller = more rubbish
            bitmapParams.sensitivity=datastruct(currentSessID).sensitivity; % higher = less stringent
            bitmapParams.desiredOrts=[0 90];

            
            % define alignment and current session data folders
            if ispc
              switch getenv('COMPUTERNAME')
                case {'LA-CPSD077020WD'} %ARC RM4
                    alignmentSessionFolder=['D:/Chip' alignmentSession '/'];
                case {'LA-CPSA07019WD'} %ARC RM5
                    alignmentSessionFolder=[mainPath 'Chip/Chip' alignmentSession '/'];
                case {'CVIS-A64882', 'PSYC-A77304'} % CPS L4
                    alignmentSessionFolder=[mainPath 'Chip/Chip' alignmentSession '/'];
                    currentSessionFolder=[mainPath 'Chip/Chip' currentSession '/'];
              end
            elseif contains(getenv('HOSTNAME'),'psy.utexas.edu')
              alignmentSessionFolder=[mainPath 'Chip/Chip' alignmentSession '/'];
              currentSessionFolder=[mainPath 'Chip/Chip' currentSession '/'];
            end
                        
            % Get projector-camera alignment transformation matrix
            load([alignmentSessionFolder 'alignmentTransform.mat'],'alignmentTransform')

            % Get all runs, iterate
            runIDs=dsCurrentSess.run;
            clear TSstruct;
            
            %% Get bitmap for mask
            saveFlag=1;
            saveFlagBMP=0;
            plotFlag=1;

            %% Plot
            phiRecruitment=[];
            for runNo=1:numel(runIDs)
                
                dsCurrentSess.run=num2str(runIDs(runNo));
                filenameStructCurrent=generateFilenames(dsCurrentSess);
                filenameStructReference=generateFilenames(dsReferenceSess);

                if runNo==1
                        pdfFilename=filenameStructCurrent.recruitmentPDF;
                        % Get within session reference orientation map
                        fastSwitch=1;
                        desiredOrts=[0 90];

                        % Get within session reference orientation map
                        [columnarBitmap,VERpca,columnarmapStats]=getColumnarBitmapV4(mainPath,dsReferenceSess,dsCurrentSess,bitmapParams, ...
                        plotFlag, saveFlag, pdfFilename);

                        % Coregister to most current green image
                        [columnarBitmapCoregistered, columnarPCAsCoregistered]=coregisterBitmap2GreenImgV2(dsReferenceSess,dsCurrentSess, ...
                        columnarBitmap,VERpca,plotFlag,saveFlag, pdfFilename);

                        % Correct for camera-projector alignment
                        orts=[0 90];%0:15:165;
                        HE=1000;
                        [projBitmapTRBB,bitmapsCamSpace,bitmapStats,bitmapStats2]=convertForProjector(dsReferenceSess,dsCurrentSess,columnarBitmapCoregistered,orts,...
                          bitmapParams.gridSize,bitmapParams.gammaCorrFactor,bitmapParams.sensitivity,'cam2proj',alignmentTransform,plotFlag,saveFlagBMP,saveFlag, pdfFilename);

                         figure('name','Neurometric plot (Recruitment dynamics)')

                end

                % Load imaging data if it exists
                if ~isempty(filenameStructCurrent.Intg) & isempty(dsCurrentSess.baselineTS) % if no separate baseline run (combined)
                    TSstruct=load(filenameStructCurrent.TS,'TS'); %TS
                    DataTrialStruct=load(filenameStructCurrent.Intg,'DataTrial'); %integrated data
                    load(filenameStructReference.Orientation,'Mask','RespCondPCA','Ort') % d' ROI mask, RespCondPCA
                    if isfile(filenameStructCurrent.gaussianFit)
                        load(filenameStructCurrent.gaussianFit) % .gaussian file
                    else
                        ROIMaskgaussian=[];
                    end
                elseif ~isempty(filenameStructCurrent.Intg) & ~isempty(dsCurrentSess.baselineTS) % if has separate baseline run
                    TSstruct(1)=load(filenameStructCurrent.TS,'TS'); %TS
                    DataTrialStruct(1)=load(filenameStructCurrent.Intg,'DataTrial'); %integrated data
                    load(filenameStructCurrent.gaussianFit); % .gaussian file
                    load(filenameStructReference.Orientation,'Mask','RespCondPCA','Ort'); % d' ROI mask, RespCondPCA
                    if ~isempty(filenameStructCurrent.baselineIntg)
                        %TSstruct(2)=load(filenameStructCurrent.baselineTS,'TS');
                        %DataTrialStruct(2)=load(filenameStructCurrent.baselineIntg,'DataTrial'); %integrated data
                    end
                    if isfile(filenameStructCurrent.gaussianFit)
                        load(filenameStructCurrent.gaussianFit) % .gaussian file
                    else
                        ROIMaskgaussian=[];
                    end
                end

                % --- Neurometric biasing (integrated responses) --- 
                disp([unique(TSstruct.TS.Header.Conditions.ProjTTLPulseOn) unique(TSstruct.TS.Header.Conditions.TTL2PulseOn)])
                powerRatio=round(TSstruct.TS.Header.Conditions.ProjTTLPulseOn / (TSstruct.TS.Header.Conditions.ProjTTLPulseOn+TSstruct.TS.Header.Conditions.ProjTTLPulseOff),2) * 100;

                [phiRecruit, images]=analyzeProjection1DFix(dsCurrentSess, TSstruct, DataTrialStruct, Mask, ROIMaskgaussian, bitmapsCamSpace, RespCondPCA, 'nonstim', powerRatio, saveFlag); upFontSize(32,.01);
                 colormap(inferno);caxis([-.0025 0.0025]);caxis([-.0015 0.0015])
                
                 % Save for plotting
                 phiRecruitment=[phiRecruitment;phiRecruit];

                 plotImg.fullROI(:,:,:,runNo)=images.fullROI;
                 plotImg.recruitROI(:,:,:,runNo)=images.recruitROI;

                 plotImg.fullROIcolumnar(:,:,:,runNo)=images.fullROIcolumn;
                 plotImg.recruitROIcolumnar(:,:,:,runNo)=images.recruitROIcolumn;

                 plotImg.temporalON(runNo)=unique(TSstruct.TS.Header.Conditions.ProjTTLPulseOn * 100 ./ (TSstruct.TS.Header.Conditions.ProjTTLPulseOn + TSstruct.TS.Header.Conditions.ProjTTLPulseOff));
            end
            if saveFlag
                export_fig(pdfFilename,'-pdf','-nocrop','-append');
            end

            %% Plot
            % sort according to powerr
            [~,sortOrder]=sort(plotImg.temporalON);
            plotImg.temporalON=plotImg.temporalON(sortOrder);
            phiRecruitment=phiRecruitment(sortOrder);
            plotImg.fullROI=plotImg.fullROI(:,:,:,sortOrder);
            plotImg.recruitROI=plotImg.recruitROI(:,:,:,sortOrder);
            plotImg.fullROIcolumnar=plotImg.fullROIcolumnar(:,:,:,sortOrder);
            plotImg.recruitROIcolumnar=plotImg.recruitROIcolumnar(:,:,:,sortOrder);
            
            % Usage example for 'Raw full ROI'
            titleStr='Raw (Recruitment ROI)';
            figure('name', titleStr);
            nCols = size(plotImg.recruitROI, 4); nRows = size(plotImg.recruitROI, 3) + 1;
            plotRecruitment(plotImg, 'recruitROI', nRows, nCols, .8, titleStr);
            if saveFlag
                export_fig(pdfFilename,'-pdf','-nocrop','-append');
            end
            titleStr='Columnar (Recruitment ROI)';
            figure('name', titleStr);
            nCols = size(plotImg.recruitROIcolumnar, 4); nRows = size(plotImg.recruitROIcolumnar, 3) + 1;
            plotRecruitment(plotImg, 'recruitROIcolumnar', nRows, nCols, .3, titleStr);
            if saveFlag
                export_fig(pdfFilename,'-pdf','-nocrop','-append');
            end
        end

    case {'psychfit'} %% Psychfit
        psychometrics=[];
        contrastsSortedCont=nan(6*2,3,numel(analysisBlockID));
        clusterIdx=[5,5,5,5,3,4,3,4,3,3,3,4,3,4,3,3,3,3,3,4,4,4,3,3,3,3,4,4,3,3,3,2,3,2,3,1,2,3,3,1,2,1,1];
        for cluster=[3]%sort(unique(clusterIdx))
            disp(['=== Cluster ' num2str(cluster) '===='])
            clusterBlocks=find(clusterIdx==cluster);
            
            for blockID=18%clusterBlocks(find(bitmapData.nColumns(2,clusterBlocks) ==1))%clusterBlocks%1:numel(analysisBlockID)%1:numel(analysisBlockID)
                currentSessID=analysisBlockID(blockID);
                dsCurrentSess=datastruct(currentSessID); %fixed, to get ort map projected onto
                dsReferenceSess=datastruct(dsCurrentSess.referenceBlockNo);
                alignmentSession=datastruct(dsCurrentSess.alignmentBlockNo).date;

                saveFlag=0;
                plotFlag=1;
                isSummary=1;

                % slope and intercept
                % --- Behavioral biasing ---
                figure()
                hAx=gca;
                [betaSorted,muSorted,sigmaSorted,contrastsSorted,percentVerticalSorted]=analyzeBiasingBlock(dsCurrentSess, bitmapData, blockID,hAx,saveFlag,'psychfit');

                % Behavioral biasing
                betaSortedCont(blockID,:)=betaSorted;
                muSortedCont(blockID,:)=muSorted;
                sigmaSortedCont(blockID,:)=sigmaSorted;
                contrastsSortedCont(1:size(contrastsSorted,1),:,blockID)=contrastsSorted;
                percentVerticalSortedCont(1:size(contrastsSorted,1),:,blockID)=percentVerticalSorted;

                psychometrics(blockID).betaSorted=betaSorted;
                psychometrics(blockID).muSorted=muSorted;
                psychometrics(blockID).sigmaSorted=sigmaSorted;
                psychometrics(blockID).contrastsSorted=contrastsSorted;
                psychometrics(blockID).percentVerticalSorted=percentVerticalSorted;
                psychometrics(blockID).percentHorizontalSorted=100-percentVerticalSorted;

                % calculate correct
                horzCond=find(contrastsSorted<0);
                percentVerticalSorted(horzCond)=100-percentVerticalSorted(horzCond);
                psychometrics(blockID).percentCorrectSorted=percentVerticalSorted;
                %nColumnsCont(blockID,:)=dsCurrentSess.nColumns;
                xticks(-100:25:100)
            end
        end
        
        %% Per cluster: Minimum columns x beta\
        nClusters=sort(unique(clusterIdx));
        mkrColors=slanCM('bold',20);mkrColors=mkrColors([2 6 4 8 10],:);
        figure('Name','Min. columns')
        yline(0,'LineStyle','--', 'LineWidth', 2.5); hold on
        for cluster=nClusters
            disp(['=== Cluster ' num2str(cluster) '===='])
            clusterBlocks=find(clusterIdx==cluster);
            nClusterBlocks=numel(clusterBlocks);
            % Get values
            nColumnsCluster=mean(bitmapData.nColumns(:,clusterBlocks)',2);
            bitmapEnergyCluster=mean(squeeze(bitmapData.energy(:,:,clusterBlocks))',2);
            betaCluster = vertcat(psychometrics(clusterBlocks).betaSorted) .* 100;
            deltaBeta=(betaCluster(:,3)-betaCluster(:,1)) ./ 2;
            % plot
            plot(nColumnsCluster,deltaBeta, 'square', 'MarkerSize', 15, 'MarkerFaceColor', mkrColors(cluster,:), 'MarkerEdgeColor', 'k', 'linewidth', 2.5); hold on
            ylabel('$\bar{\Delta\beta}$', 'FontName', 'Arial','Interpreter','latex'); % Label the x-axis as 'AUC'
            xlabel('No. of columns (block average)'); % Label the y-axis as 'Mean nColumns'
            title('Effect of power and no. of columns stimulated','FontWeight','normal', 'Interpreter', 'none') % Title for the plot
            legend({'Chance','C1','C2','C3','C4','C5'},'Location','eastoutside', 'NumColumns', 1)
            xticks([0:5:40])
            yticks([-50:10:50])
            xlim([0 40])
            ylim([-50 50])
            upFontSize(24,.01)
        end
        hold off

        
        
        % Session summary
        sessionSEM
        
        % stats
        ranksum(betaSortedCont(:,3),betaSortedCont(:,1),'tail','both')
        
        % plot minimum columns
        hBeta=betaSortedCont(:,1);
        vBeta=betaSortedCont(:,3);
        deltaBeta=round(vBeta-hBeta,2);
        figure('name','Minimum columns for biasing')
        scatter(nColumnsCont,deltaBeta,100,'ok','linewidth',2)
        xticks([1 5:5:40])
        xlim([0 40])
        upFontSize(13,0.005)
        1;
    
    case {'stability'}
        referenceSessID=12; % 1-VSD or 4-GCaMP 20220920
        dsReferenceSess=datastruct(referenceSessID);

        nOrtStandard=12;
        ortGap12=165/(nOrtStandard-1);
        ortAllStandard=0:(ortGap12*nOrtStandard/nOrtStandard):165;

        sessionsWanted=[15 14 12 11 9 8];
        nSessions=numel(sessionsWanted);
        sessionMaps=cell(1,nSessions);
        columnarMaps=cell(1,nOrtStandard);
        similarityScores=cell(1,nOrtStandard);
        for sessionNo=1:nSessions
            currentSessID=sessionsWanted(sessionNo);
            % Define sessions
            dsCurrentSess=datastruct(currentSessID);
            % Get within session ort map
            gammaCorrFactor=1;
            thresholdPrctile=0;
            
            [~, columnarMapHE,colVERRef,colVERpcaRef]=getColumnarBitmap(dsCurrentSess,gridSize,gammaCorrFactor,thresholdPrctile,'R',0); %R for reference
            % Coregister to most current green image
            [columnarMapHE]=coregisterBitmap2GreenImgV2(dsCurrentSess,dsReferenceSess,columnarMapHE,1);
            % Make sure its formatted for 12 ort
            nOrtSession=size(columnarMapHE,3);
            ortAll=0:(ortGap12*nOrtStandard/nOrtSession):165;
            [~,ortIdx]=intersect(ortAll,ortAllStandard);
            for ort=ortIdx'
                % each cell in columnarMaps contains maps of all sessions for that orientation (care 8 ort sessions)
                columnarMaps{ort}=cat(3,columnarMaps{ort},columnarMapHE(:,:,ort));
            end
        end
        
        
        %% Plot stack overlay per orientation
        [~,ortIdx]=intersect([0 45 90 135],ortAllStandard);
        imgDimX=size(columnarMaps{1},1);
        imgDimY=size(columnarMaps{1},2);

        for ort=ortIdx
            figure('name','PCA-ed map across sessions')
            [hA,~]=tight_subplot(1,numel(sessionsWanted)+1);
            columnarMapSessions=columnarMaps{ort};
            for sess=1:size(columnarMapSessions,3)
                axes(hA(sess))
                % display
                
                imgDims=size(columnarMapSessions(:,:,sess));
                
                % adaptive contrast
                rescaled=rescale(columnarMapSessions(:,:,sess),-1,1);
                histEqImg(:,:,sess)=adapthisteq(rescaled,'NumTiles',imgDims./gridSize,'Range','full');     %histEqImg=histeq(rescaled,1000);
                
                %gamma
                gammaCorr=40;
                gammaCol(:,:,sess)=imadjust(histEqImg(:,:,sess),[],[],gammaCorr);
                
                % adaptive threshold
                threshold = adaptthresh(gammaCol(:,:,sess), sensitivity);
                columnarBitmap(:,:,sess) = imbinarize(gammaCol(:,:,sess),threshold); %zero mask
                
                % plot
                imagesc(columnarBitmap(:,:,sess));axis square;colormap(fire); hold on; colorbar;
                title(datastruct(sessionsWanted(sess)).date) % add similarity of each session to reference (output from coregisterBitmap2GreenImg)
                addPix2MM(1,imgDimX,1,imgDimY,2,2,2);
                upFontSize(13, 0.02)
            end
            axes(hA(sess+1))
            imagesc(mean(columnarBitmap(:,:,:),3));axis square;colormap(fire); colorbar;
            title('Average across all sessions') % add similarity of each session to reference (output from coregisterBitmap2GreenImg)
            addPix2MM(1,imgDimX,1,imgDimY,2,2,2);
            upFontSize(13, 0.02)
            suplabel(sprintf('%.0f columns(averaged, %d sessions)',ortAllStandard(ort),size(columnarMapSessions,3)),'t',[.08 .08 .84 .86]);
            export_fig(pdfFilename,'-pdf','-nocrop');
        end
        
    case {'neurometric'}
        open analyzeNeurometrics
        saveFlag=1;
        plotFlag=1;

        %% Fast pipeline (coregisters within session ort map to optostim block green image)
        % Get correction for camera-projector alignment
        % imregtform of orignal and recovered bitmap, apply to correct for
        % camera-projector alignment

        % Get projector-camera alignment transformation matrix
        load([alignmentSessionFolder 'alignmentTransform.mat'],'alignmentTransform')

        % Define sessions
        dsAnalysisSess=datastruct(analysisBlockID);

        % Get within session reference orientation map
        analyzeColumnarSignal(dsAnalysisSess,'R',plotFlag); %R for reference
end
