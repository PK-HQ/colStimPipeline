function compareAndPlotDistributions(paramsStruct, chamberWanted, savefilename, saveFlag)
    monkeyName='Chip';

    % Naming and colors
    fieldNames = {'beta', 'rmax', 'exponent', 'c50', 'xDelta'};
    titleStrs = {'\beta', 'R_{max}', 'Exponent', 'C_{50}', '\Deltax'};
    lineColors = {[0.2, 0.2, 0.2], [0.9294, 0.1098, 0.1373] * 1.05, [0, 0.0941, 0.6627] * 1.25};
    legendNames = {'Baseline', 'Con-opto', 'Incon-opto'}; % Legend names for param1, param2, param3

    for fieldIndex = 1:numel(fieldNames)
        fieldName = fieldNames{fieldIndex};
        % Extract the full field data
        fieldData = vertcat(paramsStruct.(fieldName));
        
        % Prepare the data subsets
        param1 = fieldData(:,1);
        param2 = fieldData(:,2);
        param3 = fieldData(:,3);

        % Ensure vectors are the same length for pairwise comparison
        minLength = min([length(param1), length(param2), length(param3)]);
        param1 = param1(1:minLength);
        param2 = param2(1:minLength);
        param3 = param3(1:minLength);

        % Statistical tests
        [pVal12, ~, ~] = signrank(param1, param2);
        [pVal13, ~, ~] = signrank(param1, param3);
        [pVal23, ~, ~] = signrank(param2, param3);

        % Plot distributions
        figure('Name', sprintf('%s Comparisons', fieldName));
        pairs = {{param1, param2}, {param1, param3}, {param2, param3}};
        pVals = [pVal12, pVal13, pVal23];
        colorPairs = {[1, 2], [1, 3], [2, 3]};
        pairLegendNames = {{[legendNames{1} ' (' num2str(median(param1),2) ')'],...
                                            [legendNames{2} ' (' num2str(median(param2),2) ')']},...
                                            {[legendNames{1} ' (' num2str(median(param1),2) ')'],...
                                            [legendNames{3} ' (' num2str(median(param3),2) ')']},...
                                            {[legendNames{2} ' (' num2str(median(param2),2) ')'],...
                                            [legendNames{3} ' (' num2str(median(param3),2) ')']}};

        for i = 1:length(pairs)
            subplot(1, 3, i);
            pair = pairs{i};
            combinedParams = [pair{1}; pair{2}];
            currentMax = max(combinedParams);
            currentMin = min(combinedParams);

            if currentMax < 10
                binWidth = 1;
                xlimits = [0, currentMax+1];
            elseif currentMax >= 10 && currentMax < 50
                binWidth = 2.5;
                xlimits = [0, currentMax+5];
            else
                binWidth = 10;
                xlimits = [0, 100];
            end

            colorIndex1 = colorPairs{i}(1);
            colorIndex2 = colorPairs{i}(2);

            plotCenteredHistogram(pair{1}, binWidth, xlimits, lineColors{colorIndex1});
            xline(median(pair{1}),'--','LineWidth',2,'HandleVisibility','off');
            hold on;
            plotCenteredHistogram(pair{2}, binWidth, xlimits, lineColors{colorIndex2})
            xline(median(pair{2}),'--','LineWidth',2,'HandleVisibility','off'); hold on;

            xlim(xlimits);
            ylim([0 1]);
            axis square;
            title(sprintf('p = %.3f', pVals(i)));
            xlabel([upper(fieldName(1)) fieldName(2:end)]);
            ylabel('Probability');
            upFontSize()
            legend(pairLegendNames{i}, 'Location', 'northwest'); offwarning()
        end
        suplabel(titleStrs{fieldIndex},'t',[.1 .1 .835 .78]);upFontSize();

        % Save as PDF
        savePDF(savefilename, monkeyName, saveFlag, fieldIndex)
    end
end

%% === Subfunctions ===
function plotCenteredHistogram(data, binWidth, xlimits, faceColor)
    % plotCustomHistogram Plots a histogram with specified bin width and color
    % 
    % Inputs:
    %   data        - Numeric array of data to be plotted.
    %   binWidth    - Width of each bin.
    %   colorIndex1 - Index for selecting the color from lineColors.
    %   xl          - 2-element array [xmin xmax] specifying x-limits.
    %   lineColors  - Cell array of colors.

    % Extend the limits to fit the bin width exactly
    if mod(diff(xlimits), binWidth) ~= 0
        %xlimits(2) = xlimits(1) + ceil(diff(xlimits)/binWidth) * binWidth;
    end

    % Create bin edges based on the x-limits
    binEdges = (xlimits(1)-binWidth/2):binWidth:(xlimits(2)+binWidth/2); %xlimits(1):binWidth:xlimits(2);

    % Plot the histogram
    histogram(data, 'BinEdges', binEdges, 'Normalization', 'probability', ...
              'FaceColor', faceColor, 'FaceAlpha', 0.7);

    % Set the x-axis limits to ensure it uses your predefined limits
    xlim(xlimits);
end
