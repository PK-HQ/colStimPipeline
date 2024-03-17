function plotPowerSeries(bitmapData, behavioralData, nColumns)
    % Ensure dependent functions or colormaps like fireice are available

    % Set parameters based on the 'combined' condition
    setStr = 'combined';
    switch setStr
        case 'combined'
            setID = 3;
            offset = 1;
    end

    % Determine usable blocks excluding specific conditions
    usableBlocks = setdiff(1:size(bitmapData.nColumns,2),...
        unique(round(find(abs(bitmapData.nColumns-nColumns) <= 1 == 0)/2)));

    % Define colors and markers
    cmap = fireice(); % Ensure this colormap function is defined or available
    lineColor = [cmap([33, 11], :); .52, 0, .84];
    markerColor = [cmap([33, 11], :); .52, 0, .84];
    markerType = {'^', 'v', 's'};

    % Plotting setup
    figure('Name', 'Power series (20 column)');
    nRows = 1; nCols = 2;
    gap = .05; marginV = .02; marginH = .03;
    [hAx, ~] = tight_subplot(nRows, nCols, [gap gap], [marginV + .1 marginV + .15], [marginH + .05 marginH]);
    
    % Plot for Congruent and Incongruent Optostim
    axes(hAx(1));
    for optoConIncon = 1:2
        x = squeeze(bitmapData.adjustedSPD_uW(1, optoConIncon, usableBlocks));
        y = squeeze(behavioralData.auc(optoConIncon + offset, setID, usableBlocks));

        hold on;
        %plot(x, y, markerType{optoConIncon}, 'MarkerSize', 15, 'Color', lineColor(optoConIncon, :), 'MarkerFaceColor', markerColor(optoConIncon, :), 'MarkerEdgeColor', 'k', 'linewidth', 2.5);
        %x=nColumns(set,:);
        %y=squeeze(behavioralData.auc(set+1,3,usableBlocks))';
        plotFitCurve(x, y, markerColor(optoConIncon,:), markerType{optoConIncon}, hAx, 1);

    end
    title('Congruent versus incongruent','FontWeight','normal')
    ylabel('AUC (%)'); % Label the x-axis as 'AUC'
    xlabel(['Surface power density (' char(0181) 'W mm^{-2})'])
    xlim([0 1]); ylim([50 100])
    line([-100 100], [50 50], 'Color', .5*[1 1 1], 'LineStyle', '--', 'LineWidth', 2, 'HandleVisibility', 'on');upFontSize(24,.01)
    legend({'Congruent','Congruent-Fit','Incongruent','Incongruent-Fit','Chance'},'Location','southwest')
    upFontSize(24,.01)

    % Plot for Congruent minus Incongruent Optostim
    axes(hAx(2));
    for optoConIncon = 3
        x = mean(squeeze(bitmapData.adjustedSPD_uW(1, 1:2, usableBlocks)));
        y = squeeze(behavioralData.auc(optoConIncon + offset, setID, usableBlocks));

        hold on;
        %plot(x, y, markerType{optoConIncon}, 'MarkerSize', 15, 'Color', lineColor(optoConIncon, :), 'MarkerFaceColor', markerColor(optoConIncon, :), 'MarkerEdgeColor', 'k', 'linewidth', 2.5);
        plotFitCurve(x, y, markerColor(optoConIncon,:), markerType{optoConIncon}, hAx, 2);
    end
    title('Difference','FontWeight','normal')
    ylabel('\Delta AUC (%)'); % Label the x-axis as 'AUC'
    xlabel(['Surface power density (' char(0181) 'W mm^{-2})'])
    xlim([0 1]); ylim([-40 40])
    line([-100 100], [0 0], 'Color', .5*[1 1 1], 'LineStyle', '--', 'LineWidth', 2, 'HandleVisibility', 'on');
    legend({'Difference','Difference-Fit','Chance'},'Location','southwest')
    upFontSize(24,.01)

    suplabel({'Power series',...
        [num2str(nColumns) ' ' char(0177) ' 1 columns, n_{blocks}=' num2str(numel(usableBlocks))]}, 't', [0.1 0.1 .8 .84])
    
    upFontSize(24,.01)

end

function plotFitCurve(x, y, markerColor, markerType, hAx, set)
% Function to plot fitted curve with scatter plot overlay
% Implement plotting based on provided example

% Select axis for plotting
axes(hAx(set));
hold on;

% Curve plotting
if markerColor==[1 1 1]
    lineColor='k';
else
    lineColor=markerColor;
end
plot(x, y, markerType, 'MarkerSize', 15, 'Color', lineColor, 'MarkerFaceColor', markerColor, 'MarkerEdgeColor', 'k', 'linewidth', 2.5);
fitSaturatingCurve(x, y, markerColor, hAx, set)
end
