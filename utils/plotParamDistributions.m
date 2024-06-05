function plotParamDistributions(data1,data2,fieldName,titleStr, modelType, monkeyName, savefilename, saveFlag)
    
    nPlots=size(fieldName,1);
    
    lineColors = {.3*[1 1 1], [0.9294, 0.1098, 0.1373] * 1.07, [0, 0.0941, 0.6627] * 1.25};
    legendNames = { [(modelType{1}) '_{median}'], [(modelType{2}) '_{median}' ]}; % Legend names for param1, param2, param3
    
    % Prepare the data subsets
    param1 = data1.fittedParams(:,end);
    param2 = data2.fittedParams(:,end);
    fieldData=[param1;param2];
    
    % Statistical tests
    [pVal12, ~, ~] = signrank(param1, param2, 'tail', 'right');

    % Plot distributions
    figure('Name', sprintf('%s Comparisons', fieldName));
    pairs = {{param1, param2}};
    pVals = [pVal12];
    colorPairs = {[1, 2]};
    pairLegendNames = {[legendNames{1} ' = ' num2str(round(median(param1)))], [legendNames{2} ' = ' num2str(round(median(param2)))]};

    for subplt = 1%length(pairs)
        subplot(1, 1, subplt);
        pair = pairs{subplt};
        currentMax = max(abs(fieldData(:)));

        if currentMax < 5
            binWidth = 1;
            xMax=roundup((currentMax+binWidth),binWidth);
            xlimits = [findNearestMin(fieldData(:), -xMax), xMax];
        elseif currentMax >= 5 && currentMax < 30
            binWidth = 3;
            xMax=roundup((currentMax+binWidth),binWidth);
            xlimits = [findNearestMin(fieldData(:), -xMax), xMax];
        else
            binWidth = 5;
            xMax=roundup((currentMax+binWidth),binWidth);
            xlimits = [findNearestMin(fieldData(:), -xMax), xMax];
        end

        colorIndex1 = colorPairs{subplt}(1);
        colorIndex2 = colorPairs{subplt}(2);

        plotCenteredHistogram(pair{1}, binWidth, xlimits, lineColors{colorIndex1});
        xline(median(pair{1}),'--','LineWidth',2,'Color',lineColors{colorIndex1},'HandleVisibility','off');
        hold on;
        plotCenteredHistogram(pair{2}, binWidth, xlimits, lineColors{colorIndex2})
        xline(median(pair{2}),'--','LineWidth',2,'Color',lineColors{colorIndex2},'HandleVisibility','off'); hold on;

        xlim(xlimits);
        ylim([0 1]);
        axis square;
        title({titleStr,sprintf('p_{%s < %s} = %0.1f x 10^{%i}', modelType{2}, modelType{1}, 10^mod(log10(pVals(subplt)),1),floor(log10(pVals(subplt))))});

        if subplt==1
            xlabel('AICc');
            ylabel('Probability');
        end
        upFontSize()
        legend(pairLegendNames{[1:2]}, 'Location', 'northwest'); offwarning()
    end

    % Save as PDF
    savePDF(savefilename, monkeyName, saveFlag, 1, nPlots)
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
              'FaceColor', faceColor, 'FaceAlpha', 0.5);

    % Set the x-axis limits to ensure it uses your predefined limits
    xlim(xlimits);
end

function closestValue = findNearestMin(fieldData, limit)
    % Find the minimum value in the array
    minValue = min(fieldData);

    % Calculate the distance to 0 and -500
    distanceToZero = abs(minValue - 0);
    distanceToNegativeLimit = abs(minValue + limit);

    % Determine which is closer and return the closest value
    if distanceToZero < distanceToNegativeLimit
        closestValue = 0;
    else
        closestValue = distanceToNegativeLimit;
    end
end


