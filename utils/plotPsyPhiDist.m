function plotPsyPhiDist(datastruct, mainPath, monkeyName, chambers, chamberIDs, nBlockStr, filterColumns, saveFlag)
    % Function to plot psy-phi correlation for clusters
    %
    % Args:
    %   chambers: List of chamber IDs
    %   datastruct: Data structure with analysis block information
    %   bitmapData: Bitmap data structure
    %   mdlStruct: Model structure containing psy data
    %   neuroStruct: Neuro data structure containing phi data
    
    xlims=[-50 50];

    % Init predictor struct
    data.deltaPsy=[];
    data.deltaPhi=[];
    data.meanColumns=[];
    data.meanEnergy=[];
    % colors
     purple=[156, 14, 254]/255;

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
            clustersDesired=1;
        elseif strcmp(chamberWanted,'R')
            clustersDesired=1;
        end
        loadFlag=exist('dataTag');
        if ~loadFlag
            %load([mainPath monkeyName '/Meta/summary/statistics' chamberWanted 'tag.mat'], 'bitmapData');
            %load([mainPath monkeyName '/Meta/psychometrics/' chamberWanted '-chamber/weibullfreeAll/mdlStruct' chamberWanted '20.mat'], 'mdlStruct');
            load([mainPath '/' monkeyName '/Meta/summary/statistics' chamberWanted '-final' nBlockStr '.mat'],'blockData','bitmapData','behavioralData','analysisBlockID','datastruct','dataTag','mdlStruct')

            %load([mainPath 'Chip/Meta/neurometric/neuroStruct' chamberWanted '.mat'], 'neuroStruct');
        elseif loadFlag
            if ~strcmp(dataTag,chamberWanted)
                %load([mainPath monkeyName '/Meta/summary/statistics' chamberWanted 'tag.mat'], 'bitmapData');
                %load([mainPath monkeyName '/Meta/psychometrics/' chamberWanted '-chamber/weibullfreeAll/mdlStruct' chamberWanted '20.mat'], 'mdlStruct');
                load([mainPath '/' monkeyName '/Meta/summary/statistics' chamberWanted '-final' nBlockStr '.mat'],'blockData','bitmapData','behavioralData','analysisBlockID','datastruct','dataTag','mdlStruct')

                %load([mainPath monkeyName '/Meta/neurometric/neuroStruct' chamberWanted '.mat'], 'neuroStruct');
            end
        end
        
        if strcmp(datastruct(1).monkeyNo,'28')
            monkeyNo='1';
        elseif strcmp(datastruct(1).monkeyNo,'32')
            monkeyNo='2';
        end
       chamberStr=chambers{chamberIDs};

        % Get datastruct block IDs
        nColumnsWanted = []; %all
        analysisBlockID = organizeBlocks(datastruct, chamberWanted, nColumnsWanted);

        % Get cluster information
        analysisParams=[];
        [~, ~, clusterIdx, ~] = clusterEnergy(squeeze(bitmapData.totalPowerToOnPixelsWithinROI_mW), squeeze(bitmapData.nColumns), 'kmeans', 1, analysisParams);

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
    %% 20 column power x biasing
    figure
    deltaBiasY = mdlStruct.([chamberWanted, 'weibullfreeAll' , 'C1']).deltaBias;
    deltaMaskY = mdlStruct.([chamberWanted, 'weibullfreeAll' , 'C1']).deltaMask;
    bitmapEnergy=mean(squeeze(bitmapData.totalPowerToOnPixelsWithinROI_mW(:,:,:))',2);
    columns=mean(bitmapData.nColumns(:,:))';
    blocksActual = find(columns >= columnsDesired-columnSpread &...
        columns <= columnsDesired+columnSpread &...
        mean(squeeze(bitmapData.orts)==[0;90])');
    
    blocksControl = find(columns >= columnsDesired-columnSpread &...
        columns <= columnsDesired+columnSpread &...
        mean(squeeze(bitmapData.orts)==[45;135])');
    
    yline(0,'LineStyle','--','color',[.5 .5 .5],'LineWidth',1.5,'HandleVisibility','off'); hold on
    % Plot actual data
    mdlScatter=fitSaturatingCurve(bitmapEnergy(blocksActual),deltaBiasY(blocksActual),  'k', 1); hold on
    mdlScatter=fitSaturatingCurve(bitmapEnergy(blocksActual),deltaMaskY(blocksActual),  'k', 1); hold on
    scatter(bitmapEnergy(blocksActual),deltaBiasY(blocksActual),200,'Marker', 'square', ...
        'MarkerFaceColor', purple, 'MarkerEdgeColor', 'k', 'LineWidth',2,'DisplayName','Con-Incon (0\circ, 90\circ)'); hold on
    scatter(bitmapEnergy(blocksActual),deltaMaskY(blocksActual),200,'Marker', 'square', ...
        'MarkerFaceColor', [.5 .5 .5], 'MarkerEdgeColor', 'k', 'LineWidth',2,'DisplayName','Base-Opto'); hold on
    % Plot control data
    if ~isempty(blocksControl)
        scatter(bitmapEnergy(blocksControl)+rand(1,5)'/10,deltaBiasY(blocksControl),200,'Marker', 'square', ...
            'MarkerFaceColor', [1 1 1]*1, 'MarkerEdgeColor', 'k', 'LineWidth',2,'DisplayName','Con-Incon (45\circ, 135\circ)'); hold on
    end
    legend('Location','best')
    title({['Power x behavior (M' monkeyNo '-' chamberStr ')'], sprintf('%.0f ± %.0f columns',columnsDesired,round(std(columns(blocksActual))))})
    axis square

    addSkippedTicks(0,4.5,.25,'x')
    addSkippedTicks(-10,50,5,'y')
    xlim([0 4.5])
    ylim([-10 50])
    %[p,h,stats] = ranksum(actualY(actualX>=3.5), [controlY(controlX>=3.5)' 2 -1 -3], 'tail','both')
    xlabel('Total power delivered (mW)')
    ylabel('Δ correct %')
    upFontSize(20,.01)
    
    % Save
    switch saveFlag
        case 1
            figureName=[mainPath monkeyName '/Meta/psychometrics/M' monkeyNo '-powercurve-' chamberStr];
            set(gcf, 'Renderer', 'painters'); % Use painters for vector graphics
            print(gcf, [figureName '.png'], '-dpng', '-r600'); % High-res PNG
            savefig(gcf, [figureName '.fig']);           % FIG
            print(gcf, [figureName '.svg'], '-dsvg');        % SVG
    end

    %% PLOT DELTA OPTO HIST
    % Setup predictors
    figure('Name', 'Psy');
   
    % Define bin edges to center bins at [0, 5, 10, ..., 50]
    binEdges = -2.5:5:52.5;  % This creates edges at [-2.5, 2.5, 7.5, 12.5, ..., 52.5]
    
    % Create the histogram with specified bin edges
    histogram(data.deltaPsy, binEdges, 'EdgeColor', 'k', 'FaceColor', purple, 'FaceAlpha', 0.8, 'LineWidth',2); hold on
    %histogram([11 -1 2 8 7], binEdges, 'EdgeColor', 'k', 'FaceColor',[255, 166, 0]/255, 'FaceAlpha',1, 'LineWidth',2); hold on
    % Add formatting
    xline(0,'--','LineWidth',2,'Color',[.65 .65 .65])
    xlabel('\DeltaCorrect (%)');
    ylabel('Count');
    title('Distribution of \DeltaCorrect', 'FontSize', 24, 'FontWeight', 'normal');
    % cosmetics
    xlim(xlims);addSkippedTicks(xlims(1),xlims(2),xlims(2)/4,'x');
    ylim([0 16]); addSkippedTicks(0,16,2,'y'); 
    title(['Behavior (M' monkeyNo '-' chamberStr ')'],'FontSize',24,'FontWeight','normal')
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
    xPos1 = xLims(1) + 0.045 * (xLims(2) - xLims(1));
    yPos1 = yLims(1) + 0.07 * (yLims(2) - yLims(1));
    xPos2 = xLims(1) + 0.06 * (xLims(2) - xLims(1));
    yPos2 = yLims(1) + 0.12 * (yLims(2) - yLims(1));
    % Add the text
    text(xPos1, yPos1, pvalueStr, 'FontSize', 14, 'Interpreter', 'tex');
    text(xPos2, yPos2, nStr, 'FontSize', 14, 'Interpreter', 'tex');

    % Save
    switch saveFlag
        case 1
            figureName=[mainPath monkeyName '/Meta/psychometrics/M' monkeyNo '-histogram-' chamberStr];
            set(gcf, 'Renderer', 'painters'); % Use painters for vector graphics
            print(gcf, [figureName '.png'], '-dpng', '-r600'); % High-res PNG
            savefig(gcf, [figureName '.fig']);           % FIG
            print(gcf, [figureName '.svg'], '-dsvg');        % SVG
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
