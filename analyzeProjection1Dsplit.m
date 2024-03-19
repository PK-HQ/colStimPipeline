clear all; clc
saveFlag=0;

%% From mainColstimPipeline
if ispc
  mainPath='Y:/';
elseif contains(getenv('HOSTNAME'),'psy.utexas.edu')
  mainPath='/eslab/data/';
end
run([mainPath 'users/PK/colStimPipeline/exptListBiasingFull.m'])

analysisSessID=fliplr([72 67 66 65 64 62]); %fliplr([93 92 90 89 88 87 84 83 81 79 78 74 73 72 67 66 65 64 62 59]);

for sessionID=1:numel(analysisSessID)
    %% STAGE 1. Load integrated cortical response, TS, gaussian ROI mask
    % Get filenames for loading data
    currentSessID=analysisSessID(sessionID);
    dsCurrentSess=datastruct(currentSessID); %fixed, to get ort map projected onto
    dsReferenceSess=datastruct(dsCurrentSess.referenceSessionEntryNo); %fixed, to get ort map projected onto
    filenameStructCurrent=generateFilenames(dsCurrentSess);
    filenameStructReference=generateFilenames(dsReferenceSess);
    pdfFilename=filenameStructCurrent.recruitmentPDF;
    
    % Load imaging data if it exists 
    if isfile(filenameStructCurrent.Intg) 
        load(filenameStructCurrent.TS,'TS') %TS 
        load(filenameStructCurrent.Intg,'DataTrial') %integrated data
        load(filenameStructCurrent.gaussianFit) % .gaussian file
        load(filenameStructReference.Orientation,'Mask','RespCondPCA','Ort') % d' ROI mask, RespCondPCA
    end

    %% STAGE 2. For each intg image, split image into stimulated and nonstimulated area with gaussian ROI mask
    %% Get behavioral and neural responses from valid completed trials
    % Process integrated image to get the average image per condition (successfully completed trials only)
    [trials,images,condIDs]=getUsableTrials(TS,DataTrial);
    [images.averageColumn]=columnarFilter(TS,images.average);
    %% Get reference PCA map for comparison
    % Process make ROI mask
    ROImask=double(Mask); snrMask=ROImask; snrMask(snrMask==0)=NaN; 
    gaussMask=ROIMaskgaussian;
    % mask ort map, apply vasculature transform (morphs ort map to the imaging window view from experiment)
    PCAmap=RespCondPCA .* snrMask;
    PCAmap=transformImage(PCAmap,filenameStructCurrent);
    
    %% STAGE 3. Split the image into optostim and recruitment ROI
    [fullROI, optostimROI, recruitROI]=parcellation(condIDs,dsCurrentSess,images,snrMask,gaussMask);

    %% STAGE 4. Project into PCA space, plot vectors
    %% Create difference map from the visual-only images (reference map)
    figure('name','Vector plot')
    nRows=3;nCols=2;
    gap=.05;marginV=.04;marginH=.04;
    [hAx,~]=tight_subplot(nRows,nCols,[gap+.05 gap], [marginV marginV+.1], [marginH+.02 marginH]);
    contrasts=unique(TS.Header.Conditions.StimCon(find(TS.Header.Conditions.TypeCond==3)));
    mkrSize=rescale(contrasts, 100, 250);
    plotOrder=[1 2; 3 4; 5 6];
    plotLimit=[];
    for condSet=1:3
        axes(hAx(plotOrder(condSet,1)));
        % Define conditions to plot
        switch condSet
            case {1} %BL V-only
                conds=[condIDs.V0;...
                    condIDs.V90]-nblank; % store as 2 rows, one for each set of conds
                condsetStr='Baseline (visual-only)';
                condMarker={'o','o'};
                condMarkerColor={'b','r'};
                condLineColor={'b','r'};
            case {2} %Con V+O
                conds=[condIDs.V0O0;...
                    condIDs.V90O90]-nblank;
                condsetStr='Visual + optostim (congruent)';
                condMarker={'v','^'};
                condMarkerColor={'b','r'};
                condLineColor={'b','r'};
            case {3} %Incon V+O
                conds=[condIDs.V0O90;...
                    condIDs.V90O0]-nblank; % flipped to match visual 0 and 90
                condsetStr='Visual + optostim (incongruent)';
                condMarker={'^','v'};
                condMarkerColor={'r','b'};
                condLineColor={'b','r'};
        end

        % Plot a rose plot per condset
        visualOrientations=[0+270 90]; %hacky solution to show 90 as 'north' and 0 as 'south' (270)
        recruitmentProjections=[];
        for visualOrt=1:size(conds,1)
            switch visualOrt
                case {1}
                    pcaOrt=1;
                    multiplier=-1; % flip to -ve if 0
                case {2}
                    pcaOrt=7;
                    multiplier=1;
            end
            condsWanted=conds(visualOrt,:); % select the row of conds wanted
            nConds=numel(condsWanted);
            % Extract images
            visualOrientation=visualOrientations(visualOrt);
            referenceImage=PCAmap(:,:,pcaOrt); % this assumes that visual dominates and all responses resemble visual (may not be true for incon conditions)
            fullImages=fullROI(:,:,condsWanted);
            recruitImages=recruitROI(:,:,condsWanted);

            %compute projection vector between recruitment difference image and difference map
             [angle, amplitude, angleAverage, amplitudeAverage]=calculateProjectionVector(recruitImages, referenceImage, visualOrientation);

            recruitmentProjections(visualOrt,1:nConds,1)=angle; %recruitmentPCAMap(i,j,1) contains angle
            recruitmentProjections(visualOrt,1:nConds,2)=amplitude.*multiplier; %recruitmentPCAMap(i,j,2) contains amplitude, flip if 0
            %{
            polarplot(deg2rad(angle),amplitude.*multiplier,...
                condMarker{visualOrt},'MarkerFaceColor',condMarkerColor{visualOrt},'MarkerEdgeColor','k',...
                'Color',condLineColor{visualOrt},'LineWidth',2); hold on;
            %}
            polarscatter(deg2rad(angle),amplitude.*multiplier,mkrSize,condMarkerColor{visualOrt},...
                'MarkerEdgeColor','k','MarkerFaceColor',condMarkerColor{visualOrt},'LineWidth',1.5,'Marker',condMarker{visualOrt}); hold on
            upFontSize(14,.01)
            plotLimit=[plotLimit,amplitude]; 
        end
        hold off;
        pax=gca;
        pax.RLim=[0 .6];%[0 roundup(max(plotLimit,[],'all'),.2)];
        pax.ThetaDir='counterclockwise';
        ThetaZeroLocation='bottom';
        pax.ThetaLim=[90 270];
        pax.ThetaTick=[90 180 270];
        pax.ThetaTickLabels={'90°',condsetStr,'0°'};
        if condSet==1
            title('Projection vector','FontWeight','normal');
        end
        pax.GridAlpha=1;

        
        
        %% Scatter
        for visualOrt=1:size(conds,1)
            switch visualOrt
                case {1}
                    pcaOrt=1;
                    multiplier=-1;
                case {2}
                    pcaOrt=7;
                    multiplier=1;
            end
            condsWanted=conds(visualOrt,:); % select the row of conds wanted
            nConds=numel(condsWanted);
            % Extract images
            visualOrientation=visualOrientations(visualOrt);
            referenceImage=PCAmap(:,:,pcaOrt); % this assumes that visual dominates and all responses resemble visual (may not be true for incon conditions)
            fullImages=fullROI(:,:,condsWanted);
            recruitImages=recruitROI(:,:,condsWanted);

            %compute projection vector between recruitment difference image and difference map
             [angle, amplitude, angleAverage, amplitudeAverage]=calculateProjectionVector(recruitImages, referenceImage, visualOrientation);
            
            axes(hAx(plotOrder(condSet,2)));
            
            legLineStr={'Visual 0°','Visual 90°'};
            legMkrStr={'Opto 0°','Opto 90°'};

            plot(contrasts.*multiplier,amplitude,'-','LineWidth',1.5,...
                'Marker',condMarker{visualOrt},'MarkerSize',1,'MarkerFaceColor',condMarkerColor{visualOrt},'MarkerEdgeColor','k',...
                'Color',condLineColor{visualOrt},...
                'DisplayName',legLineStr{visualOrt}); hold on;
            scatter(contrasts.*multiplier,amplitude,mkrSize,'MarkerEdgeColor','k','MarkerFaceColor',condMarkerColor{visualOrt},...
                'LineWidth',1.5,'Marker',condMarker{visualOrt},...
                'DisplayName',legMkrStr{visualOrt}); hold on;
            xlim([-max(contrasts),max(contrasts)])
            yLims=[-1 1] .* max(pax.RLim);
            ylim(yLims)
            addSkippedTicks(yLims(1),yLims(2),.1,'y')
            addSkippedTicks(-max(contrasts),max(contrasts),5,'x')
            yline(0,'--','Color',[.4 .4 .4],'LineWidth',1.5,'HandleVisibility','off')           
            upFontSize(14,.01)
            if condSet==1
                xlabel('Contrast (%)')
                ylabel('Amplitude')
                title('Vector x contrast','FontWeight','Normal')
                x=(max(contrasts)-min(contrasts))/10;y=((yLims(2)-yLims(1)) * .9)+yLims(1);text(x,y,'Vertical','Color','r')
                x=(max(contrasts)-min(contrasts))/10;y=((yLims(2)-yLims(1)) * .1)+yLims(1);text(x,y,'Horizontal','Color','b')
                upFontSize(14,.01)
            end
        end
    end
    [~,h]=suplabel('Vector projections (reference map space)','t',[.08 .08 .84 .86]); set(h,'FontWeight','normal');
        
    %% Stage 4B. Plotting correlation of responses to visual and PCA difference reference map  
    if saveFlag
        export_fig(pdfFilename,'-pdf','-nocrop','-append');
    end
end



