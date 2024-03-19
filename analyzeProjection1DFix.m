function [phi,images]=analyzeProjection1DFix(dsCurrentSess, TSstruct, DataTrialStruct, Mask, ROIMaskgaussian, bitmapMask, RespCondPCA, optostimAreaFlag, powerRatio, saveFlag)
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
[images.fullROIcolumn, images.optostimROIcolumn, images.recruitROIcolumn]=parcellation(condIDs,dsCurrentSess,images.averageColumn,snrMask,gaussMask,bitmapMask);
[images.fullROI, images.optostimROI, images.recruitROI]=parcellation(condIDs,dsCurrentSess,images.average,snrMask,gaussMask,bitmapMask);
% Display for sanity check
%axes(hAx(1)); title('Subtracted map (full ROI)')
%imgsc(images.fullROIcolumn(:,:,2)-images.fullROIcolumn(:,:,1)); hold on; addPix2MM(images.fullROIcolumn(:,:,2)-images.fullROIcolumn(:,:,1),1,1,3); cax;
%axes(hAx(2)); title('Subtracted map (recruitment ROI)')
%imgsc(images.recruitROIcolumn(:,:,2)-images.recruitROIcolumn(:,:,1)); hold on; addPix2MM(images.recruitROIcolumn(:,:,2)-images.recruitROIcolumn(:,:,1),2,1,3); cax(.4);
%% STAGE 4. Project into PCA space, plot vectors
%% Create difference map from the visual-only images (reference map)
%axes(hAx(hAxPanel))
mkrSize=15;lineWidth=3;
plotLimitX=[];plotLimitY=[];
optostimOrientations=extractOptostimOrt(TS);
%colorCodes={'-RdPu','-Blues','-BuGn','Purples'};colorSchemeA=brewermap(43,colorCodes{4});colorSchemeB=brewermap(43,colorCodes{3});
for condSet=1:2
    contrasts=[];
    amplitudes=[];
    % Define conditions to plot
    switch condSet
        case {1} %Con V+O
            conds=[condIDs.opto0Conds];
            condsetStr='Visual + optostim-0°';
            condMarker='v';
            colorRow=[156, 14, 254]/255;
            condMarkerColor=colorRow;
            condLineColor=colorRow;
            contrast=unique(TSstruct(1).TS.Header.Conditions.GGContrastCond(find(TSstruct(1).TS.Header.Conditions.TypeCond==2)));
        case {2} %Incon V+O
            conds=[condIDs.opto90Conds]; % flipped to match visual 0 and 90
            condsetStr='Visual + optostim-90°';
            condMarker='^';
            colorRow=[0, 225, 80]/255;
            condMarkerColor=colorRow;
            condLineColor=colorRow;
            contrast=unique(TSstruct(1).TS.Header.Conditions.GGContrastCond(find(TSstruct(1).TS.Header.Conditions.TypeCond==2)));
    end

    %% Plot lineplot
    for cond=conds
        %code 0° as -ve
        if cond==1
            contrastSigned=contrast.*-1; 
        else
            contrastSigned=contrast*1;
        end
        condsWanted=cond; % select the row of conds wanted
        nConds=numel(condsWanted);
        % Extract images
        visualOrientation=optostimOrientations(cond);
        referenceImage=rescale(PCAmap(:,:,7),-1,1)-rescale(PCAmap(:,:,1),-1,1);
        
        % Choose the area analyzed
        switch optostimAreaFlag
            case {'full'}
                columnarResp=images.fullROIcolumn(:,:,condsWanted);
                titleStr='Neurometric (full ROI)';
            case {'stim'}
                columnarResp=images.optostimROIcolumn(:,:,condsWanted);
                titleStr='Neurometric (Stimulated ROI)';
            case {'nonstim'}
                columnarResp=images.recruitROIcolumn(:,:,condsWanted);
                titleStr='Neurometric (Recruitment ROI)';
        end
        
       %compute projection vector between recruitment difference image and difference map
        [~, amplitude, ~, ~]=calculateProjectionVector(columnarResp, referenceImage, cond);
        contrasts=[contrasts,contrastSigned];
        amplitudes=[amplitudes,amplitude];
        plotLimitX=[plotLimitX,contrastSigned];
        plotLimitY=[plotLimitY,amplitude];
    end
    
    % Check if 0 exists, if yes merge contrasts and amplitudes
    zeroExist=ismember(0,contrasts);
    if zeroExist
        zeroIdx=find(contrasts==0);
        if numel(zeroIdx)>1
            contrasts(zeroIdx(2))=[];
            amplitudes(zeroIdx(1))=mean(amplitudes(zeroIdx));
            amplitudes(zeroIdx(2))=[];
        else
            amplitudes(zeroIdx(1))=mean(amplitudes(zeroIdx));
        end
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
    
    plot(powerRatio,amplitudesSort,'-','LineWidth',lineWidth,...
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
maxX=roundup(max(abs(plotLimitX)),20);
if maxX==0
    maxX=1;
    xtickSpacing=1;
else
    xtickSpacing=20;
end
ylim([-maxY maxY])

xlabel('Power (LED time-on %)')
ylabel('\Phi')

if numel(phi)==3
    h=title({titleStr,...
             ['\Phi: ' num2str(phi(2),'%.2f') ', '  num2str(phi(1),'%.2f') ', ' num2str(phi(3),'%.2f')]},...
             'FontWeight','normal','FontSize',14,'FontName','Arial');%,'interpreter','latex');
elseif numel(phi)==2
    h=title({titleStr,...
         ['\Phi: ' num2str(phi(2),'%.2f') ', '  num2str(phi(1),'%.2f')]},...
         'FontWeight','normal','FontSize',14,'FontName','Arial');%,'interpreter','latex');
end

yLims=ylim;
x=-.9*max(plotLimitX);y=((yLims(2)-yLims(1)) * .95)+yLims(1);%text(x,y,'Vertical','Color',[0, 225, 80]/255,'FontWeight','Bold')
x=-.9*max(plotLimitX);y=((yLims(2)-yLims(1)) * .05)+yLims(1);%text(x,y,'Horizontal','Color',[156, 14, 254]/255,'FontWeight','Bold')
upFontSize(14,.01)
axis square; hold on;

%% Add inset of subtracted image
%{
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
        columnarResp=images.fullROIcolumn;
    case {'stim'}
        columnarResp=images.optostimROIcolumn;
    case {'nonstim'}
        columnarResp=images.recruitROIcolumn;
end
imgsc(columnarResp(:,:,condIDs.opto90Conds(1)) - columnarResp(:,:,condIDs.opto0Conds(1))); 
cax()
axis off; colormap fireice; hold on;
%}
hold on;
%% Stage 4B. Plotting correlation of responses to visual and PCA difference reference map  
if saveFlag
    export_fig(pdfFilename,'-pdf','-nocrop','-append');
end
hold on;
end



