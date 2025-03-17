function plotPsyPhiDist(datastruct, mainPath,chambers, chamberIDs, filterColumns, saveFlag)
    % Function to plot psy-phi correlation for clusters
    %
    % Args:
    %   chambers: List of chamber IDs
    %   datastruct: Data structure with analysis block information
    %   bitmapData: Bitmap data structure
    %   mdlStruct: Model structure containing psy data
    %   neuroStruct: Neuro data structure containing phi data

    % Init predictor struct
    data.deltaPsy=[];
    data.deltaPhi=[];
    data.meanColumns=[];
    data.meanEnergy=[];

    % Check if chambers contains only L, only R, or both
    hasL = any(strcmp(chambers, 'L'));
    hasR = any(strcmp(chambers, 'R'));
    if hasL && hasR
        chamberStr='LR';
    elseif hasL && ~hasR
        chamberStr='L';
    elseif hasR && ~hasL
        chamberStr='R';
    end
    for chamberID = chamberIDs % Change this to loop over chambers if needed
        chamberWanted = chambers{chamberID};
        if strcmp(chamberWanted,'L')
            clustersDesired=2:4;
        elseif strcmp(chamberWanted,'R')
            clustersDesired=1;
        end
        loadFlag=exist('dataTag');
        if ~loadFlag
            load([mainPath 'Chip/Meta/summary/statistics' chamberWanted '.mat'], 'bitmapData');
            load([mainPath 'Chip/Meta/psychometrics/' chamberWanted '-chamber/weibullfreeAll/mdlStruct' chamberWanted '.mat'], 'mdlStruct');
            %load([mainPath 'Chip/Meta/neurometric/neuroStruct' chamberWanted '.mat'], 'neuroStruct');
        elseif loadFlag
            if ~strcmp(dataTag,chamberWanted)
                load([mainPath 'Chip/Meta/summary/statistics' chamberWanted '.mat'], 'bitmapData');
                load([mainPath 'Chip/Meta/psychometrics/' chamberWanted '-chamber/weibullfreeAll/mdlStruct' chamberWanted '.mat'], 'mdlStruct');
                load([mainPath 'Chip/Meta/neurometric/neuroStruct' chamberWanted '.mat'], 'neuroStruct');
            end
        end
               
        % Get datastruct block IDs
        nColumnsWanted = []; %all
        analysisBlockID = organizeBlocks(datastruct, chamberWanted, nColumnsWanted);

        % Get cluster information
        analysisParams=[];
        [~, ~, clusterIdx, ~] = clusterEnergy(squeeze(bitmapData.energy), squeeze(bitmapData.nColumns), 'bin', 2, analysisParams);

        nClusters = numel(unique(clusterIdx));

        % Define colors and markers
        conColor = [0, 225, 80] / 255; % Green
        inconColor = [156, 14, 254] / 255; % Purple
        conInconColor = 'k'; % Black for differences
        conMarker = '^'; % Upward triangle
        inconMarker = 'v'; % Downward triangle
        conInconMarker = 's'; % Square for differences
        edgeColor = 'k'; % Black edge

        % Iterate through clusters
        clusterData=[];
        for cluster = clustersDesired
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
            clusterPsyCon = nanmean(squeeze(mdlStruct.([chamberWanted 'weibullfreeAllC' num2str(cluster)]).yBlock(2, :, clusterBlocksIdx)));
            clusterPsyIncon = nanmean(squeeze(mdlStruct.([chamberWanted 'weibullfreeAllC' num2str(cluster)]).yBlock(3, :, clusterBlocksIdx)));
            clusterPsy = [clusterPsyCon; clusterPsyIncon]';

            % Get phi data
            %clusterPhiCon = neuroStruct.(['C' num2str(cluster)]).refmap.averageProjection(:, 2)';
            %clusterPhiIncon = neuroStruct.(['C' num2str(cluster)]).refmap.averageProjection(:, 3)';
            %clusterPhi = [clusterPhiCon; clusterPhiIncon]';

            % Second row: Scatterplot of differences
            clusterPsyDiff = clusterPsy(:, 1) - clusterPsy(:, 2);
            %clusterPhiDiff = clusterPhi(:, 1) - clusterPhi(:, 2);

            % Store predictor data
            data.deltaPsy=[data.deltaPsy;clusterPsyDiff];
            data.meanColumns=[data.meanColumns;mean(blockColumns,2)];
            data.meanEnergy=[data.meanEnergy;mean(blockEnergy,2)];
        end
    end
   
    %% PLOT DELTA OPTO SCATTER
    %{
    condStrs={'Con-opto','Incon-opto'};
    condSaveStrs={'Con','Incon'};
    nConds=size(data.(chamberWanted).psy,2);
    nBlocks=size(data.(chamberWanted).psy,1);
    % Setup predictors
    figure('Name', 'Psy');
    for cluster=clustersDesired
        clusterDataIdx=find(clusterData==cluster);
        blockNo=1:numel(clusterDataIdx);
        dataPsyCon=data.(chamberWanted).psy(clusterDataIdx,1);
        dataPsyIncon=data.(chamberWanted).psy(clusterDataIdx,2);
        dataPhiCon=data.(chamberWanted).phi(clusterDataIdx,1);
        dataPhiIncon=data.(chamberWanted).phi(clusterDataIdx,2);
        dataEnergyAvg=mean(data.(chamberWanted).energy(clusterDataIdx,1:2),2);
        scatter(dataEnergyAvg, dataPsyCon-dataPsyIncon, 250, 'ksq', 'LineWidth', 2.5,'markerFaceColor',[156, 14, 254] / 255,'markerFaceAlpha',.8); hold on;%[0.9294, 0.1098, 0.1373] * 1.05
    end
    yline(0,'--','LineWidth',2,'Color',[.65 .65 .65])
    xlabel('Power (mW mm-2)');
    ylabel('\DeltaPerformance (%)');
    xlim([0 .6]); addSkippedTicks(0,.6,.1,'x')
    ylim([-10 50]); addSkippedTicks(-10,50,5,'y')
    title('Behavior','FontSize',24,'FontWeight','normal')
    upFontSize(24, .01); axis square
    legend({'Con-Incon'},'Location','northeast','FontSize', 18)
    axis square
    % Save
    switch saveFlag
        case 1
            set(gcf, 'Renderer', 'painters'); % Use painters for vector graphics
            print(gcf, [mainPath 'Chip/Meta/psychometrics/averagePsyDiff3.png'], '-dpng', '-r600'); % High-res PNG
            savefig(gcf, [mainPath 'Chip/Meta/psychometrics/averagePsyDiff3.fig']);           % FIG
            print(gcf, [mainPath 'Chip/Meta/psychometrics/averagePsyDiff3.svg'], '-dsvg');        % SVG
    end
    %}
    %% PLOT DELTA OPTO HIST
    % Setup predictors
    figure('Name', 'Psy');

    % Define bin edges to center bins at [0, 5, 10, ..., 50]
    binEdges = -2.5:5:52.5;  % This creates edges at [-2.5, 2.5, 7.5, 12.5, ..., 52.5]
    
    % Create the histogram with specified bin edges
    histogram(data.deltaPsy, binEdges, 'EdgeColor', 'k', 'FaceColor', [156, 14, 254]/255, 'FaceAlpha', 0.8, 'LineWidth',2);
    
    % Add formatting
    xline(0,'--','LineWidth',2,'Color',[.65 .65 .65])
    xlabel('\DeltaPerformance (%)');
    ylabel('Count');
    title('Distribution of \DeltaCorrect', 'FontSize', 24, 'FontWeight', 'normal');
    % cosmetics
    addSkippedTicks(-60,60,10,'x'); xlim([-50 50]);
    addSkippedTicks(0,8,1,'y'); ylim([0 8]);
    title('Behavior','FontSize',24,'FontWeight','normal')
    upFontSize(24, .01); axis square
    legend({'Con-Incon'},'Location','northwest','FontSize', 18)
    axis square

    % === Significance test ===
    % Perform two-tailed Wilcoxon signed-rank test against 0
    [p, h] = signrank(data.deltaPsy, 0, 'alpha', 0.05, 'tail', 'both');
    
    % Format p-value string for display in figure
    if p < 0.0001
        % Get the exponent and mantissa for scientific notation
        exponent = floor(log10(p));
        mantissa = p / 10^exponent;
        
        % Create formatted string with italicized p
        pvalueStr = ['{\it p} = ' num2str(mantissa, '%.2f') '×10^{' num2str(exponent) '}'];
    else
        % For p ≥ 0.0001, use regular formatting with 4 decimal places
        pvalueStr = ['{\it p} = ' num2str(p, '%.4f')];
    end
    nStr=['n = ' num2str(numel(data.deltaPsy))];
    % Add text to bottom left of the current figure
    % Get the current axis limits
    xLims = xlim;
    yLims = ylim;
    
    % Calculate position (5% from left and bottom edges)
    xPos1 = xLims(1) + 0.06 * (xLims(2) - xLims(1));
    yPos1 = yLims(1) + 0.07 * (yLims(2) - yLims(1));
    xPos2 = xPos1;
    yPos2 = yLims(1) + 0.12 * (yLims(2) - yLims(1));
    % Add the text
    text(xPos1, yPos1, pvalueStr, 'FontSize', 14, 'Interpreter', 'tex');
    text(xPos2, yPos2, nStr, 'FontSize', 14, 'Interpreter', 'tex');

    % Save
    switch saveFlag
        case 1
            set(gcf, 'Renderer', 'painters'); % Use painters for vector graphics
            print(gcf, [mainPath 'Chip/Meta/psychometrics/averageDeltaPsy' chamberStr '.png'], '-dpng', '-r600'); % High-res PNG
            savefig(gcf, [mainPath 'Chip/Meta/psychometrics/averageDeltaPsy' chamberStr '.fig']);           % FIG
            print(gcf, [mainPath 'Chip/Meta/psychometrics/averageDeltaPsy' chamberStr '.svg'], '-dsvg');        % SVG
    end
end

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
