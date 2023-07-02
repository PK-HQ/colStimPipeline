function analyzeTrajectory(dataStructReference,dataStructCurrent)
%% README
% Visualizes 1P calcium trajectory of population. Defines PC on trialwise data, then projects trialwise data
% (grouped into conditions) into the top-3 PCs

%% LOAD FILES AND SET PARAMETERS
%%% Generate filenames
if ispc
  mainPath='Y:/';
elseif contains(getenv('HOSTNAME'),'psy.utexas.edu')
  mainPath='/eslab/data/';
end
filenameStructReference=generateFilenames(dataStructReference);
filenameStructCurrent=generateFilenames(dataStructCurrent);
filenameReferenceY=[mainPath  dataStructReference.monkey '/' dataStructReference.monkey dataStructReference.date '/green_binned.bmp'];
filenameReferenceOI=['D:/' dataStructReference.monkey dataStructReference.date '/green_binned.bmp'];
filenameTargetY=[mainPath  dataStructCurrent.monkey '/' dataStructCurrent.monkey dataStructCurrent.date '/green2_binned.bmp'];
filenameTargetOI=['D:/'  dataStructCurrent.monkey dataStructCurrent.date '/green2_binned.bmp'];
pdfFilename=filenameStructCurrent.neurometricPDF; % PDF save file name

%% ANALYZE SESSION DATA
if isfile(filenameStructCurrent.Intg) % Check if intg response file exist
    
    
    
    %% LOAD TS, RAW AND INTEGRATED CORTICAL RESPONSE, ROI MASK, GREEN IMAGES FOR COREGISTRATION
    load(filenameStructCurrent.TS,'TS') % TS file, for extracting usable trials (no error codes). Excluded TS.Trial.Outcome, TS.Trial.CurrCond 
    intgDataTrial=load(filenameStructCurrent.Intg, 'DataTrial'); intgDataTrial=intgDataTrial.DataTrial; % Integrated data, for defining PCs
    rawDataTrial=load(filenameStructCurrent.Raw, 'DataTrial'); rawDataTrial=rawDataTrial.DataTrial; % Raw data
    load(filenameStructReference.Orientation,'Mask') % Load ROI mask, to be used with coregistration tranformation matrix below
    
    if isfile(filenameReferenceY) % Load reference green images for coregistration
        imgReference=imread(filenameReferenceY);    
    else
        imgReference=imread(filenameReferenceOI);    
    end
    if isfile(filenameTargetY) % Load target green images for coregistration
        imgTarget=imread(filenameTargetY);
    elseif isfile(filenameTargetOI)
        imgTarget=imread(filenameTargetOI);
    end
    
        
    %% COREGISTER ROI MASK
    % Get coregistration params for greenImageReference (to be transformed) to greenImageSession (anchor image)
    sameSession=isequal(dataStructCurrent.date,dataStructReference.date);
    switch sameSession
        case {1}
            method='auto'; %manual/auto
        case {0}
            method='manual';
    end
    % Coregister
    [imgReferenceCoreg,coregStats,transformParams]=coregisterGreenImages(imgReference,imgTarget,method,filenameStructCurrent,Mask);
    % transform and define ROI mask
    ROImask=double(Mask);
    ROImask=double(imwarp(ROImask,transformParams,'OutputView',imref2d(size(imgTarget))));
    ROImaskNaN=ROImask;
    ROImaskNaN(ROImaskNaN==0)=NaN;
    
    %% GET USABLE TRIALS, INDEXING MATCH CORTICAL RESPONSE IMAGES (FIND TRIALS IN TS WITH SUCCESS OR FAILURE CODE ONLY)
    [trialIdx,RespCondTrial,condIDs]=getUsableTrials(TS,rawDataTrial);  % 64 pixel x 64 pixel x 24 bins x 335 trials
    RespCondTrialIntg=intgDataTrial(:,:,trialIdx.All);
    
    % Order conditions
    nOptoConds=numel(find(TS.Header.Conditions.TypeCond==3));
    blank=condIDs.blankConds;
    nblank=numel(blank);
    iCondStim = find(TS.Header.Conditions.TypeCond>0);
    nCondStim = length(iCondStim);
    
    % define condition indices
    V0=condIDs.baselineConds(1:numel(condIDs.baselineConds)/2);%[1:3:nOptoConds/2]+nblank; 
    V90=condIDs.baselineConds(numel(condIDs.baselineConds)/2 + 1 : numel(condIDs.baselineConds));%[nOptoConds/2+1:3:nOptoConds]+nblank;
    V0O0=condIDs.opto0Conds(1:numel(condIDs.opto0Conds)/2);
    V90O0=condIDs.opto0Conds(numel(condIDs.opto0Conds)/2 + 1 : numel(condIDs.opto0Conds));
    V0O90=condIDs.opto90Conds(1:numel(condIDs.opto90Conds)/2);
    V90O90=condIDs.opto90Conds(numel(condIDs.opto90Conds)/2 + 1 : numel(condIDs.opto90Conds));
    
    % plot stuff
    plotOrder=[V0' V0O0' V0O90' V90' V90O0' V90O90']-nblank;
    %plotOrder=reshape(plotOrder',[1 nOptoConds]); %to plot correctly with tight_subplot
    
    
    %% SPATIAL FFT FOR COLUMNS
    % Define average blank response, to be subtracted
    [condTrialIdx,~]=find(trialIdx.All==blank);
    blankTrialsAvg=mean(RespCondTrial.All(:,:,:,condTrialIdx),4);
    
    % Spatial filtering for columns
    stimFrames=4:4+(6-1);
    pixX=size(RespCondTrial.All,1);
    pixY=size(RespCondTrial.All,2);
    nFrames=size(RespCondTrial.All,3);
    nTrials=size(RespCondTrial.All,4);
    colRespAvgPCA=[];
    
    for optoCondID=1:nOptoConds
        [condTrialIdx,~]=find(trialIdx.All==optoCondID+nblank); %to get opto cond ID, use 1:nOptoConds + number of blank
        
        % get condition response, blank subtracted (pixX, pixY, nFrames, nTrials)
        condResp=RespCondTrial.All(:,:,:,condTrialIdx)-blankTrialsAvg;
        
        % filter for columnar response
        colRespTrial=FilterFermi2D(condResp, 0.8, 3, TS.Header.Imaging.SizePxl);
        
        % average across trials
        colRespAvg=squeeze(mean(colRespTrial(:,:,stimFrames,:),4)); % get stimulus frames, averaged across all trials
        
        colRespAvgPCA=[colRespAvgPCA reshape(colRespAvg,pixX*pixY,size(colRespAvg,3))]; %makes a pixX*pixY x frames*optoCond matrix for PCA
        
    end
    
    %% Remove the mean response map
    %columnarResponse = columnarResponse-repmat(mean(columnarResponse,4),[1,1,,1,nCondStim]);
    
    %% PCA for orientation data
    % --- Get top-3 PCs from intg data ---
    binSize=8;
    RespCondTrialIntgBin = BinND(RespCondTrialIntg,[binSize binSize]);
    nPCAComp = 3;
    %reshape into trials x neurons
    reshapedIntgData=reshape(RespCondTrialIntgBin,[size(RespCondTrialIntgBin,1)*size(RespCondTrialIntgBin,2),numel(trialIdx.All)])';
    % Get PCs
    [PCACoef,PCAScore,~,~,PCAExpl] = ...
      pca(reshapedIntgData); % compress 2D image into 1D, 2nd dim = cond
    % Extract top-3 PCs
    top3PCs=PCACoef(:,1:nPCAComp);
    
    % project trial data into PC space
    RespCondTrialProjAll=[];
    blankFrames=[];
    stimFrames=4:10;
    % --- Project trialwise data into above PC-space, one PC at a time ---
    RespCondTrialProjAll = projectPCspace(RespCondTrial.All,trialIdx.All,blankTrialsAvg,nOptoConds,nblank,[blankFrames stimFrames],PCACoef,nPCAComp);
    RespCondTrialProjCorrect = projectPCspace(RespCondTrial.All,trialIdx.Correct,blankTrialsAvg,nOptoConds,nblank,[blankFrames stimFrames],PCACoef,nPCAComp);
    RespCondTrialProjError = projectPCspace(RespCondTrial.All,trialIdx.Error,blankTrialsAvg,nOptoConds,nblank,[blankFrames stimFrames],PCACoef,nPCAComp);
    
    %{
    for optoCondID=1:nOptoConds
        [condTrialIdx,~]=find(trialIdx.All==optoCondID+nblank); %to get opto cond ID, use 1:nOptoConds + number of blank
        % get condition response averaged across all trials, blank subtracted (pixX, pixY, nFrames, nTrials)
        condResp=mean(RespCondTrial.All(:,:,:,condTrialIdx)-blankTrialsAvg,4);
        counter=1;
        for frame=[blankFrames stimFrames]
            % reshape 64 px x 64 px x 1 frame x 317 trials into 4096 x 317 trials
            reshapedData=reshape(condResp(:,:,frame),...
            [size(condResp(:,:,frame),1)*size(condResp(:,:,frame),2),size(condResp(:,:,frame),3)])';
            for nPC=1:nPCAComp
              % project data into 3-dimensional PC space (top-3 components)
              projectedData = reshapedData*top3PCs(:,nPC);
              % store as frames x opto cond x PC-coordinates 1-3
              RespCondTrialProj(counter,optoCondID,nPC)=projectedData;
            end
            counter=counter+1;
        end
    end
    %}
    % --- Get condition average (per frame, get the mean across trials) ---
    
    %% Plot trajectory in top 3 PCs
    plotOrder=[V0' V0O0' V0O90' V90' V90O0' V90O90']-nblank;
    plotOrder=plotOrder';
    nContrasts=numel(V0);
    blues=[0 0 1];%[0.2:.8/(nContrasts-1):1]'.*[0 0 1];
    reds=[1 0 0];%[0.2:.8/(nContrasts-1):1]'.*[1 0 0];
    colors=[repmat(blues,3,1); repmat(reds,3,1)];
    
    
    %% Plot lowest contrast
    figure('name','Pop. Trajectories')
    colormap(gray);
    nRows=3;
    nCols=3;
    [hAx,~]=tight_subplot(nRows,nCols,[.1 .025]);
    
    
    % All contrasts
    axes(hAx(1));
    for optoCondTypeID=1:size(plotOrder,1)
        % def conditions belonging to same type
        conds=plotOrder(optoCondTypeID,:);
        % get projected data
        data3D=nanmean(RespCondTrialProjAll(:,conds,:),2);
        
        % set linecolor
        if intersect(optoCondTypeID,1:3) % if V0
            linecolor=blues;
        elseif intersect(optoCondTypeID,4:6) % if V0
            linecolor=reds;
        end
        
        % set linestyle
        if intersect(optoCondTypeID,[1 4]) % if V0
            linestyle='-';
        elseif intersect(optoCondTypeID,[2 6]) % if V0
            linestyle='--';
        elseif intersect(optoCondTypeID,[3 5]) % if V0
            linestyle=':';
        end
        
        % set marker
        if intersect(optoCondTypeID,[1 4]) % if V0
            marker='o';
            markercolor=[0 0 0];
        elseif intersect(optoCondTypeID,[3 6]) % if V0
            marker='^';
            markercolor=reds;
        elseif intersect(optoCondTypeID,[2 5]) % if V0
            marker='v';
            markercolor=blues;
        end
        
        % Plot
        plot3(data3D(:,:,1),data3D(:,:,2),data3D(:,:,3),strcat(linestyle,marker),...
          'Color',linecolor,'LineWidth',2,...
            'MarkerEdgeColor',markercolor,'MarkerFaceColor',[1 1 1]); 
        % to show top 2 PCs with PC2 as the X-axis
        view([-90 90]);
        hold on
    end
    xlabel(['PC1 (' num2str(PCAExpl(1),'%.1f') '%)'])
    ylabel(['PC2 (' num2str(PCAExpl(2),'%.1f') '%)'])
    zlabel(['PC3 (' num2str(PCAExpl(3),'%.1f') '%)'])
    title('All contrast optostim condition (all)')
    legend({'V0','V0+O0','V0+O90','V90','V90+O0','V90+O90'})
    upFontSize(14,0.012)
    
    axes(hAx(2));
    for optoCondTypeID=1:size(plotOrder,1)
        % def conditions belonging to same type
        conds=plotOrder(optoCondTypeID,:);
        % get projected data
        data3D=nanmean(RespCondTrialProjCorrect(:,conds,:),2);
        
        % set linecolor
        if intersect(optoCondTypeID,1:3) % if V0
            linecolor=blues;
        elseif intersect(optoCondTypeID,4:6) % if V0
            linecolor=reds;
        end
        
        % set linestyle
        if intersect(optoCondTypeID,[1 4]) % if V0
            linestyle='-';
        elseif intersect(optoCondTypeID,[2 6]) % if V0
            linestyle='--';
        elseif intersect(optoCondTypeID,[3 5]) % if V0
            linestyle=':';
        end
        
        % set marker
        if intersect(optoCondTypeID,[1 4]) % if V0
            marker='o';
            markercolor=[0 0 0];
        elseif intersect(optoCondTypeID,[3 6]) % if V0
            marker='^';
            markercolor=reds;
        elseif intersect(optoCondTypeID,[2 5]) % if V0
            marker='v';
            markercolor=blues;
        end
        
        % Plot
        plot3(data3D(:,:,1),data3D(:,:,2),data3D(:,:,3),strcat(linestyle,marker),...
          'Color',linecolor,'LineWidth',2,...
            'MarkerEdgeColor',markercolor,'MarkerFaceColor',[1 1 1]);
        % to show top 2 PCs with PC2 as the X-axis
        view([-90 90]);
        hold on
    end
    xlabel([' '])
    ylabel([' '])
    zlabel([' '])
    title('All contrast optostim condition (correct)')
    legend({'V0','V0+O0','V0+O90','V90','V90+O0','V90+O90'})
    upFontSize(14,0.012)
    
    axes(hAx(3));
    for optoCondTypeID=1:size(plotOrder,1)
        % def conditions belonging to same type
        conds=plotOrder(optoCondTypeID,:);
        % get projected data
        data3D=nanmean(RespCondTrialProjError(:,conds(:),:),2);
        
        % set linecolor
        if intersect(optoCondTypeID,1:3) % if V0
            linecolor=blues;
        elseif intersect(optoCondTypeID,4:6) % if V0
            linecolor=reds;
        end
        
        % set linestyle
        if intersect(optoCondTypeID,[1 4]) % if V0
            linestyle='-';
        elseif intersect(optoCondTypeID,[2 6]) % if V0
            linestyle='--';
        elseif intersect(optoCondTypeID,[3 5]) % if V0
            linestyle=':';
        end
        
        % set marker
        if intersect(optoCondTypeID,[1 4]) % if V0
            marker='o';
            markercolor=[0 0 0];
        elseif intersect(optoCondTypeID,[3 6]) % if V0
            marker='^';
            markercolor=reds;
        elseif intersect(optoCondTypeID,[2 5]) % if V0
            marker='v';
            markercolor=blues;
        end
        
        % Plot
        plot3(data3D(:,:,1),data3D(:,:,2),data3D(:,:,3),strcat(linestyle,marker),...
          'Color',linecolor,'LineWidth',2,...
            'MarkerEdgeColor',markercolor,'MarkerFaceColor',[1 1 1]);     
        % to show top 2 PCs with PC2 as the X-axis
        view([-90 90]);
        hold on
    end
    xlabel([' '])
    ylabel([' '])
    zlabel([' '])
    title('All contrast optostim condition (error)')
    legend({'V0','V0+O0','V0+O90','V90','V90+O0','V90+O90'})
    upFontSize(14,0.012)
    
    
    
    
    
    
    %low contrasts
    axes(hAx(4));
    for optoCondTypeID=1:size(plotOrder,1)
        % def conditions belonging to same type
        conds=plotOrder(optoCondTypeID,:);
        % get projected data
        data3D=nanmean(RespCondTrialProjAll(:,conds(1),:),2);
        
        % set linecolor
        if intersect(optoCondTypeID,1:3) % if V0
            linecolor=blues;
        elseif intersect(optoCondTypeID,4:6) % if V0
            linecolor=reds;
        end
        
        % set linestyle
        if intersect(optoCondTypeID,[1 4]) % if V0
            linestyle='-';
        elseif intersect(optoCondTypeID,[2 6]) % if V0
            linestyle='--';
        elseif intersect(optoCondTypeID,[3 5]) % if V0
            linestyle=':';
        end
        
        % set marker
        if intersect(optoCondTypeID,[1 4]) % if V0
            marker='o';
            markercolor=[0 0 0];
        elseif intersect(optoCondTypeID,[3 6]) % if V0
            marker='^';
            markercolor=reds;
        elseif intersect(optoCondTypeID,[2 5]) % if V0
            marker='v';
            markercolor=blues;
        end
        
        % Plot
        plot3(data3D(:,:,1),data3D(:,:,2),data3D(:,:,3),strcat(linestyle,marker),...
          'Color',linecolor,'LineWidth',2,...
            'MarkerEdgeColor',markercolor,'MarkerFaceColor',[1 1 1]);       
        % to show top 2 PCs with PC2 as the X-axis
        view([-90 90]);
        hold on
    end
    xlabel([' '])
    ylabel([' '])
    zlabel([' '])
    title('Low contrast optostim condition (all)')
    legend({'V0','V0+O0','V0+O90','V90','V90+O0','V90+O90'})
    upFontSize(14,0.012)
    
    %% Plot lowest contrast
    axes(hAx(5));
    for optoCondTypeID=1:size(plotOrder,1)
        % def conditions belonging to same type
        conds=plotOrder(optoCondTypeID,:);
        % get projected data
        data3D=nanmean(RespCondTrialProjCorrect(:,conds(1),:),2);
        
        % set linecolor
        if intersect(optoCondTypeID,1:3) % if V0
            linecolor=blues;
        elseif intersect(optoCondTypeID,4:6) % if V0
            linecolor=reds;
        end
        
        % set linestyle
        if intersect(optoCondTypeID,[1 4]) % if V0
            linestyle='-';
        elseif intersect(optoCondTypeID,[2 6]) % if V0
            linestyle='--';
        elseif intersect(optoCondTypeID,[3 5]) % if V0
            linestyle=':';
        end
        
        % set marker
        if intersect(optoCondTypeID,[1 4]) % if V0
            marker='o';
            markercolor=[0 0 0];
        elseif intersect(optoCondTypeID,[3 6]) % if V0
            marker='^';
            markercolor=reds;
        elseif intersect(optoCondTypeID,[2 5]) % if V0
            marker='v';
            markercolor=blues;
        end
        
        % Plot
        plot3(data3D(:,:,1),data3D(:,:,2),data3D(:,:,3),strcat(linestyle,marker),...
          'Color',linecolor,'LineWidth',2,...
            'MarkerEdgeColor',markercolor,'MarkerFaceColor',[1 1 1]);      
        % to show top 2 PCs with PC2 as the X-axis
        view([-90 90]);
        hold on
    end
    xlabel([' '])
    ylabel([' '])
    zlabel([' '])
    title('Low contrast optostim condition (correct)')
    legend({'V0','V0+O0','V0+O90','V90','V90+O0','V90+O90'})
    upFontSize(14,0.012)
    
    %% Plot lowest contrast
    axes(hAx(6));
    for optoCondTypeID=1:size(plotOrder,1)
        % def conditions belonging to same type
        conds=plotOrder(optoCondTypeID,:);
        % get projected data
        data3D=nanmean(RespCondTrialProjError(:,conds(1),:),2);
        
        % set linecolor
        if intersect(optoCondTypeID,1:3) % if V0
            linecolor=blues;
        elseif intersect(optoCondTypeID,4:6) % if V0
            linecolor=reds;
        end
        
        % set linestyle
        if intersect(optoCondTypeID,[1 4]) % if V0
            linestyle='-';
        elseif intersect(optoCondTypeID,[2 6]) % if V0
            linestyle='--';
        elseif intersect(optoCondTypeID,[3 5]) % if V0
            linestyle=':';
        end
        
        % set marker
        if intersect(optoCondTypeID,[1 4]) % if V0
            marker='o';
            markercolor=[0 0 0];
        elseif intersect(optoCondTypeID,[3 6]) % if V0
            marker='^';
            markercolor=reds;
        elseif intersect(optoCondTypeID,[2 5]) % if V0
            marker='v';
            markercolor=blues;
        end
        
        % Plot
        plot3(data3D(:,:,1),data3D(:,:,2),data3D(:,:,3),strcat(linestyle,marker),...
          'Color',linecolor,'LineWidth',2,...
            'MarkerEdgeColor',markercolor,'MarkerFaceColor',[1 1 1]);      
        % to show top 2 PCs with PC2 as the X-axis
        view([-90 90]);
        hold on
    end
    xlabel([' '])
    ylabel([' '])
    zlabel([' '])
    title('Low contrast optostim condition (error)')
    legend({'V0','V0+O0','V0+O90','V90','V90+O0','V90+O90'})
    upFontSize(14,0.012)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    %% Plot lowest contrast
    axes(hAx(7));
    for optoCondTypeID=1:size(plotOrder,1)
        % def conditions belonging to same type
        conds=plotOrder(optoCondTypeID,:);
        % get projected data
        data3D=nanmean(RespCondTrialProjAll(:,conds(end-1),:),2);
        
        % set linecolor
        if intersect(optoCondTypeID,1:3) % if V0
            linecolor=blues;
        elseif intersect(optoCondTypeID,4:6) % if V0
            linecolor=reds;
        end
        
        % set linestyle
        if intersect(optoCondTypeID,[1 4]) % if V0
            linestyle='-';
        elseif intersect(optoCondTypeID,[2 6]) % if V0
            linestyle='--';
        elseif intersect(optoCondTypeID,[3 5]) % if V0
            linestyle=':';
        end
        
        % set marker
        if intersect(optoCondTypeID,[1 4]) % if V0
            marker='o';
            markercolor=[0 0 0];
        elseif intersect(optoCondTypeID,[3 6]) % if V0
            marker='^';
            markercolor=reds;
        elseif intersect(optoCondTypeID,[2 5]) % if V0
            marker='v';
            markercolor=blues;
        end
        
        % Plot
        plot3(data3D(:,:,1),data3D(:,:,2),data3D(:,:,3),strcat(linestyle,marker),...
          'Color',linecolor,'LineWidth',2,...
            'MarkerEdgeColor',markercolor,'MarkerFaceColor',[1 1 1]); 
        % to show top 2 PCs with PC2 as the X-axis
        view([-90 90]);
        hold on
    end
    xlabel([' '])
    ylabel([' '])
    zlabel([' '])
    title('High contrast optostim condition (all)')
    legend({'V0','V0+O0','V0+O90','V90','V90+O0','V90+O90'})
    upFontSize(14,0.012)
    
    %% Plot lowest contrast
    axes(hAx(8));
    for optoCondTypeID=1:size(plotOrder,1)
        % def conditions belonging to same type
        conds=plotOrder(optoCondTypeID,:);
        % get projected data
        data3D=nanmean(RespCondTrialProjCorrect(:,conds(end-1),:),2);
        
        % set linecolor
        if intersect(optoCondTypeID,1:3) % if V0
            linecolor=blues;
        elseif intersect(optoCondTypeID,4:6) % if V0
            linecolor=reds;
        end
        
        % set linestyle
        if intersect(optoCondTypeID,[1 4]) % if V0
            linestyle='-';
        elseif intersect(optoCondTypeID,[2 6]) % if V0
            linestyle='--';
        elseif intersect(optoCondTypeID,[3 5]) % if V0
            linestyle=':';
        end
        
        % set marker
        if intersect(optoCondTypeID,[1 4]) % if V0
            marker='o';
            markercolor=[0 0 0];
        elseif intersect(optoCondTypeID,[3 6]) % if V0
            marker='^';
            markercolor=reds;
        elseif intersect(optoCondTypeID,[2 5]) % if V0
            marker='v';
            markercolor=blues;
        end
        
        % Plot
        plot3(data3D(:,:,1),data3D(:,:,2),data3D(:,:,3),strcat(linestyle,marker),...
          'Color',linecolor,'LineWidth',2,...
            'MarkerEdgeColor',markercolor,'MarkerFaceColor',[1 1 1]);
        % to show top 2 PCs with PC2 as the X-axis
        view([-90 90]);
        hold on
    end
    xlabel([' '])
    ylabel([' '])
    zlabel([' '])
    title('High contrast optostim condition (correct)')
    legend({'V0','V0+O0','V0+O90','V90','V90+O0','V90+O90'})
    upFontSize(14,0.012)
    
    %% Plot lowest contrast
    axes(hAx(9));
    for optoCondTypeID=1:size(plotOrder,1)
        % def conditions belonging to same type
        conds=plotOrder(optoCondTypeID,:);
        % get projected data
        data3D=nanmean(RespCondTrialProjError(:,conds(end-1),:),2);
        
        % set linecolor
        if intersect(optoCondTypeID,1:3) % if V0
            linecolor=blues;
        elseif intersect(optoCondTypeID,4:6) % if V0
            linecolor=reds;
        end
        
        % set linestyle
        if intersect(optoCondTypeID,[1 4]) % if V0
            linestyle='-';
        elseif intersect(optoCondTypeID,[2 6]) % if V0
            linestyle='--';
        elseif intersect(optoCondTypeID,[3 5]) % if V0
            linestyle=':';
        end
        
        % set marker
        if intersect(optoCondTypeID,[1 4]) % if V0
            marker='o';
            markercolor=[0 0 0];
        elseif intersect(optoCondTypeID,[3 6]) % if V0
            marker='^';
            markercolor=reds;
        elseif intersect(optoCondTypeID,[2 5]) % if V0
            marker='v';
            markercolor=blues;
        end
        
        % Plot
        plot3(data3D(:,:,1),data3D(:,:,2),data3D(:,:,3),strcat(linestyle,marker),...
          'Color',linecolor,'LineWidth',2,...
            'MarkerEdgeColor',markercolor,'MarkerFaceColor',[1 1 1]); 
        % to show top 2 PCs with PC2 as the X-axis
        view([-90 90]);
        hold on
    end
    xlabel([' '])
    ylabel([' '])
    zlabel([' '])
    title('High contrast optostim condition (error)')
    legend({'V0','V0+O0','V0+O90','V90','V90+O0','V90+O90'})
    upFontSize(14,0.012)
    1;
    
    
    [~,h]=suplabel(['Population dynamics during stimulus period (300ms)'],'t',[.08 .08 .84 .87]);
    export_fig(pdfFilename,'-pdf','-nocrop','-append');
else
    %
end