%% Dev notes
% Pre 12/15/2022
% To add DC squares over centroid
% apply PRF map

% 12/15/2022
% Smaller blobs, higher ON-OFF ratio
% Added analyzeNeurometrics

%% Saving figures
%{
%% Save the figures
[~,h]=suplabel(['Converting bitmap for projector dimensions (' num2str(ort) '\circ)'],'t',[.08 .08 .84 .87]);
set(h,'FontSize',16)
export_fig(pdfFilename,'-pdf','-nocrop','-append');
%}

%% Calibration
%open projectorCameraCalibration

%% Load dataStruct
% def mainPath
if ispc
  mainPath='Y:/';
elseif contains(getenv('HOSTNAME'),'psy.utexas.edu')
  mainPath='/eslab/data/';
end
%open([mainPath 'users/PK/colStimPipeline/exptListBiasingFull.m'])
run([mainPath 'users/PK/colStimPipeline/exptListBiasingFull.m'])

%% Define reference and current session
pipelineMode='';% beta/stable
analysisType='summary'; %summary/psychfit/beta/neurometric/stability

%for biasing expt
currentSessID=92;
%for analysis
analysisSessID=81%fliplr([93 92 90 89 88 87 84 83 81 79 78 74 73 72 67 66 65 64 62 59 51 50]);%sort([83 82 81 79 78 77 76 74 73 72 67 66 65 64 62 59 58 51 50],'ascend'); %63 %[50 51]; %fliplr([43 44 45 47 48, 52 53 54 57 58 59 62 63 64 65 66]); %[12 14 15 17 19 20 21 23 24 26 29 31 32 34 37 40]; %17; %10 excluded because baseline different
%[79 78 77 76 74 73 72 68 67 66 65 64 62 59 58 51 50]
% Columnar bitmap processing params
referenceSessID=datastruct(currentSessID).referenceSessionEntryNo; %choose ID with a session and run containing orientation map
alignmentSessID=datastruct(currentSessID).alignmentSessionEntryNo;
currentSession=datastruct(currentSessID).date;
alignmentSession=datastruct(alignmentSessID).date;
bitmapParams.gridSize=datastruct(currentSessID).gridSize; %smaller = more rubbish; splits image into NumTiles=512x512/gridSize for adaptive histeq
bitmapParams.gammaCorrFactor=datastruct(currentSessID).gammaCorrFactor; %smaller = more rubbish
bitmapParams.sensitivity=datastruct(currentSessID).sensitivity; % higher = less stringent
bitmapParams.desiredOrts=[0 90];

% Datapath
if ispc
  switch getenv('COMPUTERNAME')
    case {'LA-CPSD077020WD'} %ARC RM4
        alignmentSessionFolder=['D:/Chip' alignmentSession '/'];
    case {'LA-CPSA07019WD'} %ARC RM5
        alignmentSessionFolder=[mainPath 'Chip/Chip' alignmentSession '/'];
    case {'SEIDEMANN2PANAL'} % CPS L4
        alignmentSessionFolder=[mainPath 'Chip/Chip' alignmentSession '/'];
        currentSessionFolder=[mainPath 'Chip/Chip' currentSession '/'];
  end
elseif contains(getenv('HOSTNAME'),'psy.utexas.edu')
  alignmentSessionFolder=[mainPath 'Chip/Chip' alignmentSession '/'];
  currentSessionFolder=[mainPath 'Chip/Chip' currentSession '/'];
end

%plot
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
        
        % Get within session reference orientation map
        [columnarBitmap,VERpca,columnarmapStats]=getColumnarBitmapV4(mainPath,dsReferenceSess,dsCurrentSess,bitmapParams, ...
          plotFlag, saveFlag);

        % Coregister to most current green image
        [columnarBitmapCoregistered, columnarPCAsCoregistered]=coregisterBitmap2GreenImgV2(dsReferenceSess,dsCurrentSess, ...
          columnarBitmap,VERpca,plotFlag,saveFlag);

        % Account for PRF
        %?
        
        % Correct for camera-projector alignment
        orts=[0 90];%0:15:165;
        HE=1000;
        
        [projBitmapTRBB,nBlobs,medianBlobAreas]=convertForProjector(dsReferenceSess,dsCurrentSess,columnarBitmapCoregistered,orts,...
            bitmapParams.gridSize,bitmapParams.gammaCorrFactor,bitmapParams.sensitivity,'cam2proj',alignmentTransform,saveFlagBMP,saveFlag);
    1;
  case {'stable'} %within session
        saveOrNot=1;
        plotFlag=1;
        
        %% Fast pipeline (coregisters within session ort map to optostim block green image)
        % Get correction for camera-projector alignment
        % imregtform of orignal and recovered bitmap, apply to correct for
        % camera-projector alignment
        % open projectorCameraCalibration
        
        
        % Get projector-camera alignment transformation matrix
        load([alignmentSessionFolder 'alignmentTransform.mat'],'alignmentTransform')

        % Define sessions
        dsReferenceSess=datastruct(referenceSessID); %moving, e.g. ort map
        dsCurrentSess=datastruct(currentSessID); %fixed, to get ort map projected onto
        
        % Get within session reference orientation map
        [columnarBitmapReference,~]=getColumnarBitmapV3(dsReferenceSess,gridSize,gammaCorrFactor,sensitivity,'R',plotFlag); %R for reference
        
        % Coregister to most current green image
        [columnarBitmapCoregistered]=coregisterBitmap2GreenImgV2(dsReferenceSess,dsCurrentSess,columnarBitmapReference,plotFlag);

        % Account for PRF
        %?
        
        % Correct for camera-projector alignment
        orts=[0:15:165];
        HE=1000;

        %make save folder
        bmpPath=sprintf('X:/PK/ColSeries/%s/',datestr(now,'yyyymmdd'));
        if ~exist(bmpPath, 'dir')
            mkdir(bmpPath)
        end
        saveOrNot=1; % save bitmap
        [projBitmapTRBB]=convertForProjector(dsCurrentSess,columnarBitmapCoregistered,orts,...
            HE,gammaCorrFactor,thresholdPrctile,'cam2proj',alignmentTransform,saveOrNot);

end

%% Run the desired analysis pipeline variant
switch analysisType
    case {'summary'}        
        for sessionID=1:numel(analysisSessID)%1:numel(analysisSessID)%1:numel(analysisSessID)%1:numel(analysisSessID)
            tic
            saveFlag=0;
            saveFlagBMP=0;
            plotFlag=1;
            % define session IDs, session meta info
            currentSessID=analysisSessID(sessionID);
            dsCurrentSess=datastruct(currentSessID); %fixed, to get ort map projected onto
            dsReferenceSess=datastruct(dsCurrentSess.referenceSessionEntryNo);
            alignmentSession=datastruct(dsCurrentSess.alignmentSessionEntryNo).date;
            
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
                case {'SEIDEMANN2PANAL'} % CPS L4
                    alignmentSessionFolder=[mainPath 'Chip/Chip' alignmentSession '/'];
                    currentSessionFolder=[mainPath 'Chip/Chip' currentSession '/'];
              end
            elseif contains(getenv('HOSTNAME'),'psy.utexas.edu')
              alignmentSessionFolder=[mainPath 'Chip/Chip' alignmentSession '/'];
              currentSessionFolder=[mainPath 'Chip/Chip' currentSession '/'];
            end
            

            
           %% Fast pipeline (coregisters within session ort map to optostim block green image)
            % Get correction for camera-projector alignment
            % imregtform of orignal and recovered bitmap, apply to correct for
            % camera-projector alignment
            % open projectorCameraCalibration

            % Get projector-camera alignment transformation matrix
            load([alignmentSessionFolder 'alignmentTransform.mat'],'alignmentTransform')


            % Get within session reference orientation map
            fastSwitch=1;
            desiredOrts=[0 90];
                
          % Get within session reference orientation map
          [columnarBitmap,VERpca,columnarmapStats]=getColumnarBitmapV4(mainPath,dsReferenceSess,dsCurrentSess,bitmapParams, ...
            plotFlag, saveFlag);

          % Coregister to most current green image
          [columnarBitmapCoregistered, columnarPCAsCoregistered]=coregisterBitmap2GreenImgV2(dsReferenceSess,dsCurrentSess, ...
            columnarBitmap,VERpca,plotFlag,saveFlag);

          % Account for PRF
          %?

          % Correct for camera-projector alignment
          orts=[0 90];%0:15:165;
          HE=1000;

          [projBitmapTRBB,bitmapStats,bitmapStats2]=convertForProjector(dsReferenceSess,dsCurrentSess,columnarBitmapCoregistered,orts,...
              bitmapParams.gridSize,bitmapParams.gammaCorrFactor,bitmapParams.sensitivity,'cam2proj',alignmentTransform,saveFlagBMP,saveFlag);

            % Static images
            [blurEstimate,gfpStaticMax,mcherryStaticMax]=QCstatic(dsCurrentSess,plotFlag,saveFlag);
            
            % slope and intercept
            % --- Behavioral biasing ---
            [betaSorted,muSorted,sigmaSorted]=analyzeBiasingBlock(dsCurrentSess,saveFlag,'full');
            % --- Neurometric biasing (integrated responses) --- 
            
            %[visualr,visualrRatio,PCAr,PCArRatio,optoIntensity]=assessOptoEffect2(dsReferenceSess,dsCurrentSess,columnarBitmapCoregistered,columnarPCAsCoregistered,saveFlag);
            % --- Neurometric biasing (dynamic frame-wise responses) --- 
            %analyzeTrajectory(dsReferenceSess,dsCurrentSess)
            

           %% Measures (one per measure, per session)

           % --- Imaging quality ---
            % Green image blur
            %blurEstimateCont(sessionID,1)=blurEstimate;
            % GFP static intensity
            gfpStaticMaxCont(sessionID,1)=gfpStaticMax;
            % mCherry static intensity
            mcherryStaticMaxCont(sessionID,1)=mcherryStaticMax;

            % --- Stimulation power and positioning ---
            % Orientation map quality
            %PCAExplTotalCont(sessionID,1)=PCAExplTotal;
            % Pixel density = Percentage pixels on within ROI
            %pixelDensitiesCont(sessionID,1:2)=nonzeros(columnarmapStats.pixDensity);
            % Power density = Percentage pixels on within ROI * power * ON-OFF
            %powerDensitiesCont(sessionID,1:2)=powerDensities;
            % No. of blobs
            %nBlobsCont(sessionID,1:2)=nonzeros(bitmapStats.nBlobs);
            % Size of blobs (median)
            %contourDC(sessionID,1:2)=nonzeros(bitmapStats.medianBlobAreas);
            contourColumns(sessionID,1:2)=nonzeros(bitmapStats.nBlobs);
            contourPixels(sessionID,1:2)=bitmapStats.contourPix;
            contourPixelDensities(sessionID,1:2)=bitmapStats.contourPixDensity;

            %pixDensityContour
            %estDutycycleContour
            %pixDensityTotal
            %estDutycycleTotal
            
            % Correlation with bitmap
            
            % --- Behavioral biasing ---
            betaSortedCont(sessionID,:)=betaSorted;
            muSortedCont(sessionID,:)=muSorted;
            sigmaSortedCont(sessionID,:)=sigmaSorted;


            % --- Neurometric biasing ---
            % Opto similarity vs visual
            %visualrCont(sessionID,1:2)=visualr;
            %visualrRatiosCont(sessionID,1:2)=visualrRatio;
            
            %PCArCont(sessionID,1:2)=PCAr;
            %PCArRatiosCont(sessionID,1:2)=PCArRatio;
           
            % Opto intensity vs visual
            %optoIntensityCont(sessionID,1:4)=optoIntensity; %2x overall mean, 2x peak mean         

            CF
            toc
        end

        getTemporalStimParams


        save([mainPath 'Chip/Meta/summary/neurometricsMinimalColumns' num2str(numel(analysisSessID)) '.mat'],...
            'analysisSessID',...
            'contourColumns','contourPixels','contourPixelDensities',...
            'betaSortedCont','muSortedCont','sigmaSortedCont',...
            'temporalTotal','temporalDC','temporalON','tbl')

        tbl=array2table([analysisSessID' contourColumns (contourPixels.*temporalON') betaSortedCont],...
            'VariableNames',{'EntryID' 'Column-H' 'Column-V', 'Energy-H', 'Energy-V','Beta-H','Beta-Baseline','Beta-V'});
       %% Minimal columns
        % Define x y z, averaged across 2 patterns       
        zPower=contourPixels.*temporalON'; 
        xColumns=contourColumns;
        yBeta=betaSortedCont(:,[1 3]);
        
        % 2D line fit (ncolumns X beta)
        metricStr='beta';
        plotMinimumColumns(xColumns,yBeta,zPower,analysisSessID,metricStr)
        metricStr='deltabeta';
        plotMinimumColumns(xColumns,yBeta,zPower,analysisSessID,metricStr)
        
        % 3D surface fit (power x ncolumns X beta)
        plotMinimumColumns3d(xColumns,yBeta,zPower,analysisSessID)

        
        
        
        
        
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
    
    
        
    case {'psychfit'}               
        for sessionID=1:numel(analysisSessID)%1:numel(analysisSessID)%1:numel(analysisSessID)
            currentSessID=analysisSessID(sessionID);
            dsCurrentSess=datastruct(currentSessID); %fixed, to get ort map projected onto
            dsReferenceSess=datastruct(dsCurrentSess.referenceSessionEntryNo);
            alignmentSession=datastruct(dsCurrentSess.alignmentSessionEntryNo).date;

            saveFlag=1;
            plotFlag=1;
            isSummary=1;

            % slope and intercept
            % --- Behavioral biasing ---
            [betaSorted,muSorted,sigmaSorted]=analyzeBiasingBlock(dsCurrentSess,saveFlag,'psychfit');
            % Behavioral biasing
            betaSortedCont(sessionID,:)=betaSorted;
            muSortedCont(sessionID,:)=muSorted;
            sigmaSortedCont(sessionID,:)=sigmaSorted;
            %nColumnsCont(sessionID,:)=dsCurrentSess.nColumns;
        end
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
        dsAnalysisSess=datastruct(analysisSessID);

        % Get within session reference orientation map
        analyzeColumnarSignal(dsAnalysisSess,'R',plotFlag); %R for reference
    
    case {'green'} % show and coregister green images
       
        referenceSessID=5; % 1-VSD or 4-GCaMP 20220920
        dsMovingImg=datastruct(referenceSessID);

        nOrtStandard=12;
        ortGap12=165/(nOrtStandard-1);
        ortAllStandard=0:(ortGap12*nOrtStandard/nOrtStandard):165;

        sessionsWanted=3:9;
        nSessions=numel(sessionsWanted);
        sessionMaps=cell(1,nSessions);
        columnarMaps=cell(1,nOrtStandard);
        similarityScores=cell(1,nOrtStandard);
        
        % figure with n-session rows x 2 columns
        figure
        tiledlayout(nSessions,3, 'Padding', 'loose', 'TileSpacing', 'compact');
        for sessionNo=1:nSessions
            currentSessID=sessionsWanted(sessionNo);
            % Define sessions
            dsFixedImg=datastruct(currentSessID);
            [similarityScore]=showGreenImage(dsMovingImg,dsFixedImg, 1);
        end
        suplabel('Coregistration of reference and session green images','t');
      
    case {'visopto'}
        %
end


%% Storage

%[columnarBitmapReference,~]=getColumnarBitmapV3(dsReferenceSess,gridSize,gammaCorrFactor,sensitivity,'R',plotFlag); %R for reference

%{
switch testStr
    case {1}
        for nOrt=1:size(bitmapSession,3)
            ds=dataStruct(4);
            imagesc(cdata>(prctile(cdata(:),99)));
            img=cdata>(prctile(cdata(:),99));
            ort=90;
            [projBitmapTRBB]=convertForProjector(ds,img,ort,'cam2proj');
            bmpFoldername = sprintf('X:/PK/ColSeries/%s/',datestr(now,'yyyymmdd'));

            bmpFilename = sprintf('X:/PK/ColSeries/%s/Col%04gHE%04gG%03gT%05g_imregtform.bmp',datestr(now,'yyyymmdd'),90,1,1,1*100);
            %exportgraphics(projBitmapTRBB,bmpFilename);%,'Resolution',300)
            if ~exist(bmpFoldername, 'dir')
               mkdir(bmpFoldername)
            end
            imwrite(projBitmapTRBB,bmpFilename);
        end
end


 case {'slow'}
    % Reference: Get column VER, column positions
    [columnarBitmapReference,~,colVERRef,colVERpcaRef]=getColumnarBitmap(dataStructReference,gammaCorrFactor,thresholdPrctile,'R'); %R for reference

    % Session: Get column VER, column positions
    [columnarBitmapSession,colVER,colVERpca]=getColVER(dataStructSession,gammaCorrFactor,thresholdPrctile,'S'); %S for session

    % Obtain bitmap by coregistering to reference bitmap and green image to the session's green image
    [bitmapSession,~]=coregisterBitmap(dataStructReference,dataStructSession,columnarBitmapSession,columnarBitmapReference,0);
    orts=[0:22.5:180];%[0:15:165];
    HE=1000;
    %make save folder
    bmpPath=sprintf('X:/PK/ColSeries/%s/',datestr(now,'yyyymmdd'));
    if ~exist(bmpPath, 'dir')
        mkdir(bmpPath)
    end
    for nOrtStandard=1:size(bitmapSession,3)
        ds=ds(4);
        img=bitmapSession(:,:,nOrtStandard);
        ort=orts(nOrtStandard);
        [projBitmapTRBB]=convertForProjector(ds,img,ort,'cam2proj');
        bmpFilename = sprintf('X:/PK/ColSeries/%s/Col%04gHE%04gG%03gT%05g_%s.bmp',datestr(now,'yyyymmdd'),orts(nOrtStandard),HE,gammaCorrFactor,thresholdPrctile*100,datestr(now,'yyyymmdd'));
        %exportgraphics(projBitmapTRBB,bmpFilename);%,'Resolution',300)
        imwrite(projBitmapTRBB,bmpFilename);
    end


%}





%% PRF/Projector gaussian light correction
%get PRF look-up table

% assign PD per pixel

% assign pixel-ON per pixel, if nPixel >=1 (convert from PD, where nPixel= desiredPD/maxPD * maxPixels)

% gamma correct and threshold

% check coregistration with green image

%% Projector check
%projectorCheck
