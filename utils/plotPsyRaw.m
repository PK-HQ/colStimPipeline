function plotPsyRaw(datastruct, mainPath,chambers, chamberIDs, filterColumns, saveFlag)
    % Function to plot psy-phi correlation for clusters
    %
    % Args:
    %   chambers: List of chamber IDs
    %   datastruct: Data structure with analysis block information
    %   bitmapData: Bitmap data structure
    %   mdlStruct: Model structure containing psy data
    %   neuroStruct: Neuro data structure containing phi data
    for chamberID = chamberIDs % Change this to loop over chambers if needed
        chamberWanted = chambers{chamberID};
        
        loadFlag=exist('dataTag');
        if ~loadFlag
            load([mainPath 'Chip/Meta/summary/statistics' chamberWanted '.mat'], 'bitmapData');
            load([mainPath 'Chip/Meta/psychometrics/' chamberWanted '-chamber/weibullfreeAll/mdlStruct' chamberWanted '.mat'], 'mdlStruct');
            load([mainPath 'Chip/Meta/neurometric/neuroStruct' chamberWanted '.mat'], 'neuroStruct');
        elseif loadFlag
            if ~strcmp(dataTag,chamberWanted)
                load([mainPath 'Chip/Meta/summary/statistics' chamberWanted '.mat'], 'bitmapData');
                load([mainPath 'Chip/Meta/psychometrics/' chamberWanted '-chamber/weibullfreeAll/mdlStruct' chamberWanted '.mat'], 'mdlStruct');
                load([mainPath 'Chip/Meta/neurometric/neuroStruct' chamberWanted '.mat'], 'neuroStruct');
            end
        end
        
        % Init predictor struct
        data.(chamberWanted).psy=[];
        data.(chamberWanted).phi=[];
        data.(chamberWanted).columns=[];
        data.(chamberWanted).energy=[];
        
        % Get datastruct block IDs
        nColumnsWanted = []; %all
        analysisBlockID = organizeBlocks(datastruct, chamberWanted, nColumnsWanted);

        % Get cluster information
        [~, ~, clusterIdx] = clusterEnergy(squeeze(bitmapData.energy), 'bin', 2);
        nClusters = numel(unique(clusterIdx));

        % Define colors and markers
        conColor = [0, 225, 80] / 255; % Green
        inconColor = [156, 14, 254] / 255; % Purple
        conInconColor = 'k'; % Black for differences
        conMarker = '^'; % Upward triangle
        inconMarker = 'v'; % Downward triangle
        conInconMarker = 's'; % Square for differences
        edgeColor = 'k'; % Black edge

        % Create figure
        figure('Name', ['Psy-Phi Correlation']);
        make_it_tight = true;
        hmarg = 0.01;
        subplot = @(m, n, p) subtightplot(m, n, p, [0.05 0.05], [hmarg hmarg], [0.1 0.1]);
        if ~make_it_tight, clear subplot; end

        % Iterate through clusters
        clusterData=[];
        for cluster = 1:nClusters
            % Select blocks for the cluster
            clusterBlocksAll = find(clusterIdx == cluster);
            
            if filterColumns
                columnsDesired = 20;
                columnSpread = 4;
                [clusterBlocks,blockColumns,blockEnergy] = selectNColumnBlocks(bitmapData, clusterBlocksAll, columnsDesired, columnSpread);
                clusterBlocksIdx = find(clusterBlocks == clusterBlocks); % make idx from itself
                1;
            else
                columnsDesired = 30;
                columnSpread = 30;
                [clusterBlocks,blockColumns,blockEnergy] = selectNColumnBlocks(bitmapData, clusterBlocksAll, columnsDesired, columnSpread);
                clusterBlocksIdx = find(clusterBlocksAll == clusterBlocksAll); % make idx from itself
            end
            clusterData=[clusterData,repmat(cluster,1,numel(clusterBlocksIdx))]; % for plotting later
            % Get psy data
            clusterPsyCon = nanmean(squeeze(mdlStruct.(['LweibullfreeAllC' num2str(cluster)]).yBlock(2, :, clusterBlocksIdx)));
            clusterPsyIncon = nanmean(squeeze(mdlStruct.(['LweibullfreeAllC' num2str(cluster)]).yBlock(3, :, clusterBlocksIdx)));
            clusterPsy = [clusterPsyCon; clusterPsyIncon]';

            % Second row: Scatterplot of differences
            clusterPsyDiff = clusterPsy(:, 1) - clusterPsy(:, 2);
            
            % Store predictor data
            data.(chamberWanted).psy=[data.(chamberWanted).psy;clusterPsy];
            data.(chamberWanted).columns=[data.(chamberWanted).columns;blockColumns];
            data.(chamberWanted).energy=[data.(chamberWanted).energy;blockEnergy];
        end
    end
    
    %% PLOT RAW DATA
    condStrs={'Con-opto','Incon-opto'};
    condSaveStrs={'Con','Incon'};
    nConds=size(data.(chamberWanted).psy,2);
    nBlocks=size(data.(chamberWanted).psy,1);
    % Setup predictors
    figure('Name', 'Psy');
    for cluster=3%1:max(clusterData)

        clusterDataIdx=find(clusterData==cluster);
        blockNo=1:numel(clusterDataIdx);
        dataPsyCon=data.(chamberWanted).psy(clusterDataIdx,1);
        dataPsyIncon=data.(chamberWanted).psy(clusterDataIdx,2);
        dataEnergyAvg=mean(data.(chamberWanted).energy(clusterDataIdx,1:2),2);
    
        % === Psy ===
        % Grouping colors
        conColors=slanCM('Reds',7);conColors=conColors(2:6,:);
        inconColors=slanCM('Blues',7);inconColors=inconColors(2:6,:);
        fprintf('=============== %s ===============\n','Psy')
        scatter(dataEnergyAvg, dataPsyCon, 250, 'k^', 'LineWidth', 2.5,'markerFaceColor',conColors(cluster,:),'markerFaceAlpha',.95); hold on;%[0.9294, 0.1098, 0.1373] * 1.05
        scatter(dataEnergyAvg, dataPsyIncon, 250, 'kv', 'LineWidth', 2.5,'markerFaceColor',inconColors(cluster,:),'markerFaceAlpha',.95); hold on;%[0.9294, 0.1098, 0.1373] * 1.05
    end
    yline(50,'--','LineWidth',2,'Color',[.65 .65 .65])
    xlabel('Power (mW mm^{-2})');
    ylabel('Performance (%)');
    xlim([0 .6]); addSkippedTicks(0,.6,.1,'x')
    ylim([0 100]); addSkippedTicks(0,100,12.5,'y')
    title('Behavior','FontSize',24,'FontWeight','normal')
    upFontSize(24, .01); axis square
    legend({'Con','Incon'},'Location','northeast','FontSize', 18)
    axis square
    % Save
    switch saveFlag
        case 1
            set(gcf, 'Renderer', 'painters'); % Use painters for vector graphics
            print(gcf, [mainPath 'Chip/Meta/psychometrics/averagePsy.png'], '-dpng', '-r600'); % High-res PNG
            savefig(gcf, [mainPath 'Chip/Meta/psychometrics/averagePsy.fig']);           % FIG
            print(gcf, [mainPath 'Chip/Meta/psychometrics/averagePsy.svg'], '-dsvg');        % SVG
    end

    %% PLOT RAW DATA
    condStrs={'Con-opto','Incon-opto'};
    condSaveStrs={'Con','Incon'};
    nConds=size(data.(chamberWanted).psy,2);
    nBlocks=size(data.(chamberWanted).psy,1);
    % Setup predictors
    figure('Name', 'Psy');        
    diffColors=slanCM('sunburst',5);diffColors=diffColors(1:5,:);
    for cluster=1:max(clusterData)
        clusterDataIdx=find(clusterData==cluster);
        blockNo=1:numel(clusterDataIdx);
        dataPsyCon=data.(chamberWanted).psy(clusterDataIdx,1);
        dataPsyIncon=data.(chamberWanted).psy(clusterDataIdx,2);
        dataEnergyAvg=mean(data.(chamberWanted).energy(clusterDataIdx,1:2),2);
        %scatter(dataEnergyAvg, dataPsyCon-dataPsyIncon, 250, 'ksq', 'LineWidth', 2.5,'markerFaceColor',[156, 14, 254] / 255,'markerFaceAlpha',.8); hold on;%[0.9294, 0.1098, 0.1373] * 1.05
        % Grouping colors
        scatter(dataEnergyAvg, dataPsyCon-dataPsyIncon, 250, 'ksq', 'LineWidth', 2.5,'markerFaceColor',diffColors(cluster,:),'markerFaceAlpha',1); hold on;%[0.9294, 0.1098, 0.1373] * 1.05
    end
    yline(0,'--','LineWidth',2,'Color',[.65 .65 .65])
    xlabel('Power (mW mm^{-2})');
    ylabel('\DeltaPerformance (%)');
    xlim([0 .6]); addSkippedTicks(0,.6,.1,'x')
    ylim([-10 50]); addSkippedTicks(-10,50,5,'y')
    title('Behavior','FontSize',24,'FontWeight','normal')
    upFontSize(24, .01); axis square
    legend({'Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5',},'Location','southeast','FontSize', 16)
    axis square
    % Save
    switch saveFlag
        case 1
            set(gcf, 'Renderer', 'painters'); % Use painters for vector graphics
            print(gcf, [mainPath 'Chip/Meta/psychometrics/averagePsyDiffAll.png'], '-dpng', '-r600'); % High-res PNG
            savefig(gcf, [mainPath 'Chip/Meta/psychometrics/averagePsyDiffAll.fig']);           % FIG
            print(gcf, [mainPath 'Chip/Meta/psychometrics/averagePsyDiffAll.svg'], '-dsvg');        % SVG
    end

    %% PLOT RAW DATA
    condStrs={'Con-opto','Incon-opto'};
    condSaveStrs={'Con','Incon'};
    nConds=size(data.(chamberWanted).psy,2);
    nBlocks=size(data.(chamberWanted).psy,1);
    % Setup predictors
    figure('Name', 'Psy');
    for cluster=2:max(clusterData)
        clusterDataIdx=find(clusterData==cluster);
        blockNo=1:numel(clusterDataIdx);
        dataPsyCon=data.(chamberWanted).psy(clusterDataIdx,1);
        dataPsyIncon=data.(chamberWanted).psy(clusterDataIdx,2);
        dataColumnsAvg=mean(data.(chamberWanted).columns(clusterDataIdx,1:2),2);

        % Grouping colors
        scatter(dataColumnsAvg, dataPsyCon-dataPsyIncon, 250, 'ksq', 'LineWidth', 2.5,'markerFaceColor',diffColors(cluster,:),'markerFaceAlpha',1); hold on;%[0.9294, 0.1098, 0.1373] * 1.05
        %scatter(dataColumnsAvg, dataPsyCon-dataPsyIncon, 250, 'ksq', 'LineWidth', 2.5,'markerFaceColor',[156, 14, 254] / 255,'markerFaceAlpha',.8); hold on;%[0.9294, 0.1098, 0.1373] * 1.05
    end
    yline(0,'--','LineWidth',2,'Color',[.65 .65 .65])
    xlabel('No. of columns');
    ylabel('\DeltaPerformance (%)');
    xlim([0 40]); addSkippedTicks(0,40,5,'x')
    ylim([-10 50]); addSkippedTicks(-10,50,5,'y')
    title('Behavior','FontSize',24,'FontWeight','normal')
    upFontSize(24, .01); axis square
    legend({'Cluster 2','Cluster 3','Cluster 4','Cluster 5',},'Location','southeast','FontSize', 16)
    axis square
    % Save
    switch saveFlag
        case 1
            set(gcf, 'Renderer', 'painters'); % Use painters for vector graphics
            print(gcf, [mainPath 'Chip/Meta/psychometrics/averagePsyDiff25.png'], '-dpng', '-r600'); % High-res PNG
            savefig(gcf, [mainPath 'Chip/Meta/psychometrics/averagePsyDiff25.fig']);           % FIG
            print(gcf, [mainPath 'Chip/Meta/psychometrics/averagePsyDiff25.svg'], '-dsvg');        % SVG
    end
end

%fitNonlinearModel(data.(chamberWanted).psy(:,cond), predictors, predictorNames,'elasticnet', 0, 1, 'false') %fits incon well
%fitNonlinearModel2(data.(chamberWanted).psy(:,cond), predictors, predictorNames,'elasticnet', 0, 0, 'false', 0, 0, 0,1) %fits con incon R2=.7
%fitNonlinearModel3(data.(chamberWanted).psy(:,cond), predictors, predictorNames,'elasticnet', includeInteractions, transformType, predictorSelection, flagAnalyticalTools) %fits con incon R2>.7

function updateMarkers(condStr, metricStr)
    % Subfunction to update scatterplot markers and colors
    
    % Define colors and markers
    switch metricStr
        case 'psy'
            conColor =  [0.9294, 0.1098, 0.1373] * 1.05; %[0, 225, 80] / 255; % Green
            inconColor =  [0, 0.0941, 0.6627] * 1.25; %[156, 14, 254] / 255; % Purple
        case 'phi'
            conColor = [0, 225, 80] / 255; % Green
            inconColor = [156, 14, 254] / 255; % Purple
    end
    conInconColor = 'k'; % Black for differences
    conMarker = '^'; % Upward triangle
    inconMarker = 'v'; % Downward triangle
    conInconMarker = 's'; % Square for differences
    markerSize=250;
    markerAlpha=.7;
    % Retrieve the current figure and axes
    figureHandle = gcf;
    axesHandle = gca;

    % Check the condition and update scatter points
    if contains(condStr, 'Con', 'IgnoreCase', false)
        % Update to con settings
        scatterData = findobj(axesHandle, 'Type', 'Scatter'); % Find scatter plots
        for i = 1:numel(scatterData)
            scatterData(i).MarkerEdgeColor = 'k'; % Update color
            scatterData(i).MarkerFaceColor = conColor; % Update face color
            scatterData(i).Marker = conMarker;        % Update marker
            scatterData(i).SizeData = markerSize;        % Update marker
            scatterData(i).MarkerFaceAlpha = markerAlpha;        % Update marker
        end
    elseif contains(condStr, 'Incon', 'IgnoreCase', false)
        % Update to incon settings
        scatterData = findobj(axesHandle, 'Type', 'Scatter');
        for i = 1:numel(scatterData)
            scatterData(i).MarkerEdgeColor = 'k';
            scatterData(i).MarkerFaceColor = inconColor;
            scatterData(i).Marker = inconMarker;
            scatterData(i).SizeData = markerSize;        % Update marker
            scatterData(i).MarkerFaceAlpha = markerAlpha;        % Update marker
        end
    else
        % Default settings for differences
        scatterData = findobj(axesHandle, 'Type', 'Scatter');
        for i = 1:numel(scatterData)
            scatterData(i).MarkerEdgeColor = 'k';
            scatterData(i).MarkerFaceColor = conInconColor;
            scatterData(i).Marker = conInconMarker;
            scatterData(i).SizeData = markerSize;        % Update marker
            scatterData(i).MarkerFaceAlpha = markerAlpha;        % Update marker
        end
    end
end
