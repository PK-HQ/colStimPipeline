function phi=analyzeProjection1D(dsCurrentSess, TSstruct, DataTrialStruct, Mask, ROIMaskgaussian,  bitmapsCamSpace, RespCondPCA, hAx, hAxPanel, optostimAreaFlag, saveFlag)
filenameStructCurrent=generateFilenames(dsCurrentSess);
pdfFilename=filenameStructCurrent.neurometricPDF;

%% STAGE 2. For each intg image, split image into stimulated and nonstimulated area with gaussian ROI mask
%% Get behavioral and neural responses from valid completed trials
% Process integrated image to get the average image per condition (successfully completed trials only)
nTS=max(size(TSstruct));
switch nTS
    case {1}
        % Extract TS and Datatrial
        TS=TSstruct(1).TS;
        DataTrial=DataTrialStruct(1).DataTrial;
        % Split by cond, blanks subtracted and removed
        [~,images,condIDs]=getUsableTrials(TS,DataTrial);
    case {2}
        % Here the baseline and opto are in separate blocks, so extract TS and Datatrial
        TSopto=TSstruct(1).TS;
        TSbaseline=TSstruct(2).TS;
        DataTrialOpto=DataTrialStruct(1).DataTrial;
        DataTrialBaseline=DataTrialStruct(2).DataTrial;
        
        % Split by cond, blanks subtracted and removed
        [~,imagesbl,condIDsbl]=getUsableTrials(TSbaseline,DataTrialBaseline);
        [~,imagesOpto,condIDsOpto]=getUsableTrials(TSopto,DataTrialOpto);

        fields=fieldnames(condIDsOpto);
        for fn=1:numel(fields)
           name=fields{fn};
          % reassign
           if ~isempty(condIDsOpto.(name))
               condsOpto(fn)=max(condIDsOpto.(name));
           else
               condsOpto(fn)=NaN;
           end
        end
        fields=fieldnames(condIDsbl);
        for fn=1:numel(fields)
           name=fields{fn};
          % reassign
           if ~isempty(condIDsbl.(name))
               condsBL(fn)=max(condIDsbl.(name));
           else
               condsBL(fn)=NaN;
           end
        end
        
        % Append baseline condIDs and images to the back of opto trials
        maxOptoConds=max(condsOpto);
        maxBLConds=max(condsBL);
        condIDsOpto.baselineConds=condIDsbl.baselineConds + maxOptoConds;
        condIDsOpto.V0=condIDsbl.V0 + maxOptoConds;
        condIDsOpto.V90=condIDsbl.V90 + maxOptoConds;
        imagesOpto.trials(:,:,maxOptoConds+1:maxOptoConds+maxBLConds,1:size(imagesbl.trials,4))=imagesbl.trials;
        imagesOpto.trialsCorrect(:,:,maxOptoConds+1:maxOptoConds+maxBLConds,1:size(imagesbl.trials,4))=imagesbl.trialsCorrect;
        imagesOpto.trialsError(:,:,maxOptoConds+1:maxOptoConds+maxBLConds,1:size(imagesbl.trials,4))=imagesbl.trialsError;
        imagesOpto.average(:,:,maxOptoConds+1:maxOptoConds+maxBLConds,:)=imagesbl.average;
        imagesOpto.averageCorrect(:,:,maxOptoConds+1:maxOptoConds+maxBLConds,:)=imagesbl.averageCorrect;
        imagesOpto.averageError(:,:,maxOptoConds+1:maxOptoConds+maxBLConds,:)=imagesbl.averageError;
        
        % rename TSopto and imagesOpto to TS and images for rest of the script
        images=imagesOpto;
        TS=TSopto;
        condIDs=condIDsOpto;
end
  
[images.averageColumn]=columnarFilter(TS,images.average);
%% Get reference PCA map for comparison
% Process make ROI mask
ROImask=double(Mask); snrMask=ROImask; snrMask(snrMask==0)=NaN;
gaussMask=ROIMaskgaussian; %rename
% mask ort map, apply vasculature transform (morphs ort map to the imaging window view from experiment)
PCAmap=RespCondPCA .* snrMask;
PCAmap=transformImage(PCAmap,filenameStructCurrent);

%% STAGE 3. Split the image into optostim and recruitment ROI
[fullROI, optostimROI, recruitROI]=parcellation(condIDs,dsCurrentSess,images,snrMask,gaussMask, bitmapsCamSpace);

%% STAGE 4. Project into PCA space, plot vectors
%% Create difference map from the visual-only images (reference map)
axes(hAx(hAxPanel))
mkrSize=15;lineWidth=3;
plotLimitX=[];plotLimitY=[];
visualOrientations=unique(TS.Header.Conditions.GaborOrt(find(TS.Header.Conditions.TypeCond==3)));
colorCodes={'-RdPu','-Blues','-BuGn','Purples'};colorSchemeA=brewermap(43,colorCodes{4});colorSchemeB=brewermap(43,colorCodes{3});
for condSet=1:3
    contrasts=[];
    amplitudes=[];
    % Define conditions to plot
    switch condSet
        case {1} %BL V-only
            conds=[condIDs.V0;...
                condIDs.V90]; % store as 2 rows, one for each set of conds
            condsetStr='Visual-only';
            condMarker='o';
            condMarkerColor='w';
            condLineColor='k';
            contrast=unique(TSstruct(end).TS.Header.Conditions.StimCon(find(TSstruct(end).TS.Header.Conditions.TypeCond==3)));
        case {2} %Con V+O
            conds=[condIDs.V0O0;...
                condIDs.V90O0];
            condsetStr='Visual + optostim-0°';
            condMarker='v';
            colorRow=[156, 14, 254]/255;
            condMarkerColor=colorRow;
            condLineColor=colorRow;
            contrast=unique(TSstruct(1).TS.Header.Conditions.StimCon(find(TSstruct(1).TS.Header.Conditions.TypeCond==3)));
        case {3} %Incon V+O
            conds=[condIDs.V0O90;...
                condIDs.V90O90]; % flipped to match visual 0 and 90
            condsetStr='Visual + optostim-90°';
            condMarker='^';
            colorRow=[0, 225, 80]/255;
            condMarkerColor=colorRow;
            condLineColor=colorRow;
            contrast=unique(TSstruct(1).TS.Header.Conditions.StimCon(find(TSstruct(1).TS.Header.Conditions.TypeCond==3)));
    end

    %% Plot lineplot
    for visualOrt=1:size(conds,1)
        %code 0° as -ve
        if visualOrt==1
            contrastSigned=contrast.*-1; 
        else
            contrastSigned=contrast*1;
        end
        condsWanted=conds(visualOrt,:); % select the row of conds wanted
        nConds=numel(condsWanted);
        % Extract images
        visualOrientation=visualOrientations(visualOrt);
        referenceImage=rescale(PCAmap(:,:,7),-1,1)-rescale(PCAmap(:,:,1),-1,1);%rescale(fullROI(:,:,condIDs.V90(end)),-1,1)-rescale(fullROI(:,:,condIDs.V0(end)),-1,1); % this assumes that visual dominates and all responses resemble visual (may not be true for incon conditions)
        
        % Choose the area analyzed
        switch optostimAreaFlag
            case {'full'}
                columnarResp=fullROI(:,:,condsWanted);
                titleStr='Neurometric (full ROI)';
            case {'stim'}
                columnarResp=optostimROI(:,:,condsWanted);
                titleStr='Neurometric (Stimulated ROI)';
            case {'nonstim'}
                columnarResp=recruitROI(:,:,condsWanted);
                titleStr='Neurometric (Recruitment ROI)';
        end
        
       %compute projection vector between recruitment difference image and difference map
        [~, amplitude, ~, ~]=calculateProjectionVector(columnarResp, referenceImage, visualOrientation);
        contrasts=[contrasts,contrastSigned];
        amplitudes=[amplitudes,amplitude];
        plotLimitX=[plotLimitX,contrastSigned];
        plotLimitY=[plotLimitY,amplitude];
    end
    
    % Check if 0 exists, if yes merge contrasts and amplitudes
    zeroExist=ismember(0,contrasts);
    if zeroExist
        zeroIdx=find(contrasts==0);
        contrasts(zeroIdx(2))=[];
        amplitudes(zeroIdx(1))=mean(amplitudes(zeroIdx));
        amplitudes(zeroIdx(2))=[];
    else
        %
    end
    
    %% Get ϕ
    [contrastsSort,idx]=sort(contrasts);
    amplitudesSort=amplitudes(idx);
    if zeroExist
        % Take single point in the middle
        phi(condSet)=amplitudesSort(find(contrastsSort==0));
    elseif ~zeroExist && mod(numel(contrastsSort),2)==0
        % Take average of two points in the middle, if no. contrasts is even
        phi(condSet)=mean([amplitudesSort(ceil(end/2)) amplitudesSort(ceil(end/2)+1)]);
    end
    
    if condSet==2 || condSet==3
        %amplitudesSort=rescale(amplitudesSort,rescaleVisLimits(1),rescaleVisLimits(2));
    end
    
    plot(contrastsSort,amplitudesSort,'-','LineWidth',lineWidth,...
    'Marker',condMarker,'MarkerSize',mkrSize,'MarkerFaceColor',condMarkerColor,'MarkerEdgeColor','k',...
    'Color',condLineColor,...
    'DisplayName',condsetStr); hold on;
    yline(0,'--','Color',[.4 .4 .4],'LineWidth',1.5,'HandleVisibility','off')           
    xline(0,'--','Color',[.4 .4 .4],'LineWidth',1.5,'HandleVisibility','off')
    upFontSize(30,.005)
    
    if condSet==1
        rescaleVisLimits=[min(amplitudesSort), max(amplitudesSort)];
    end
end
maxY=roundup(max(abs(plotLimitY)),1);
addSkippedTicks(-maxY, maxY,1,'y')
maxX=roundup(max(abs(plotLimitX)),20);
addSkippedTicks(-maxX,maxX,20,'x')
ylim([-maxY maxY])
xlim([-maxX maxX])

xlabel('Contrast (%)')
ylabel('Amplitude')

h=title({titleStr,...
         ['${\phi}$: ' num2str(phi(2),'%.2f') ', '  num2str(phi(1),'%.2f') ', ' num2str(phi(3),'%.2f')]},...
         'FontWeight','normal','FontSize',14,'FontName','Arial','interpreter','latex');

yLims=ylim;
x=-.9*max(plotLimitX);y=((yLims(2)-yLims(1)) * .95)+yLims(1);text(x,y,'Vertical','Color',[0, 225, 80]/255,'FontWeight','Bold')
x=-.9*max(plotLimitX);y=((yLims(2)-yLims(1)) * .05)+yLims(1);text(x,y,'Horizontal','Color',[156, 14, 254]/255,'FontWeight','Bold')
upFontSize(14,.01)
axis square

%% Add inset of subtracted image
switch hAxPanel
    case {2}
        ax2=axes('Position',[.36 .55 .15 .15]);
        %caxisLimits=[-.025 .025];
    case {3}
        ax2=axes('Position',[.665 .55 .15 .15]);
        %caxisLimits=[-.005 .005];
end
% Choose the area analyzed
switch optostimAreaFlag
    case {'full'}
        columnarResp=fullROI;
    case {'stim'}
        columnarResp=optostimROI;
    case {'nonstim'}
        columnarResp=recruitROI;
end
imgsc(columnarResp(:,:,condIDs.V90O90(1)) - columnarResp(:,:,condIDs.V0O0(1)));
cax()
axis off; colormap fireice;
%% Stage 4B. Plotting correlation of responses to visual and PCA difference reference map  
if saveFlag
    export_fig(pdfFilename,'-pdf','-nocrop','-append');
end
end



