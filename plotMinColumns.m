function plotMinColumns(bitmapData, behavioralData, analysisBlockID)

    % Filter by power
    selectedBlocks=powerfilter(bitmapData, analysisBlockID, 'fixedpower');

    % Extract the relevant data for plotting
    aucData = squeeze(behavioralData.auc(4,3,usableBlocks));
    nColumns = squeeze(bitmapData.nColumns(:,usableBlocks));
    nColumnsMean = squeeze(mean(bitmapData.nColumns(:,usableBlocks), 1));

    % Scatter plot
    cmap=fireice; markerColor = [cmap([43-10 1+10],:); 1 1 1; .52 0 .84]; % colors:  Black, blue, red
    markerType = {'^','v','o','s'}; % marker type

    figure('Name','Minimum columns'); % Create a new figure
    [hAx, ~] = tight_subplot(1, 2); % Adjust subplot parameters as needed

    for set=1:2
        % Select axis for plotting
        hold on;
        
        % Curve plotting
        if markerColor==[1 1 1]
            lineColor='k';
        else
            lineColor=markerColor;
        end
        x=nColumns(set,:);
        y=squeeze(behavioralData.auc(set+1,3,usableBlocks))';
        plotFitCurve(x, y, markerColor(set,:), markerType{set}, hAx, 1);
        upFontSize(24,0.01)
        xlim([0 60])
        ylim([0 100])
    end
    line([-100 100],[50 50],'Color',.6*[1 1 1], 'LineStyle','--','LineWidth',.15,'HandleVisibility','on')
    ylabel('AUC (%)'); % Label the x-axis as 'AUC'
    xlabel('No. of columns'); % Label the y-axis as 'Mean nColumns'
    legend({'Con','Con-Fit','Incon','Incon-Fit','Chance'},'Location','southwest','NumColumns',1)
    title('Congruent vs incongruent','FontWeight','normal') % Title for the plot
    xlim([0 40])
    ylim([0 100])

    for set=4
        % Select axis for plotting
        hold on;
        
        % Curve plotting
        if markerColor==[1 1 1]
            lineColor='k';
        else
            lineColor=markerColor;
        end
        x=mean(nColumns(:,:),1);
        y=squeeze(behavioralData.auc(4,3,usableBlocks))';
        plotFitCurve(x, y, markerColor(set,:), markerType{set}, hAx, 2);
        upFontSize(24,0.01)
        xlim([0 60])
        ylim([0 100])
    end
    line([-100 100],[0 0],'Color',.6*[1 1 1], 'LineStyle','--','LineWidth',.15,'HandleVisibility','on')
    ylabel('\Delta AUC (%)'); % Label the x-axis as 'AUC'
    xlabel('No. of columns (average)'); % Label the y-axis as 'Mean nColumns'
    title('Difference','FontWeight','normal') % Title for the plot
    legend({'Difference','Fit','Chance'},'Location','southwest')
    xlim([0 40])
    ylim([-50 50])
    hold off; % Release the plot hold


    suplabel({'Minimum column series',...
        ['n_{blocks}=' num2str(numel(usableBlocks))]}, 't', [0.1 0.1 .8 .80])
    upFontSize(24,0.01)
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