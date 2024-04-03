function [binnedX, meanY, binnedY] = plotSEMv2(x, y, lineColor, markerType)
    bins = -2.5:5:102.5; % Define bins

    

    % Initialize binData structure to store raw Y values
    binData = arrayfun(@(binMin, binMax) initBin(x, y, binMin, binMax), bins(1:end-1), bins(2:end), 'UniformOutput', false);
    binData = [binData{:}];
    
    if size(x,2) > 1
        % Merging logic
        binData = mergeBinData(binData);
    end
    % After merging, process bins to calculate mean and sem
    binData = arrayfun(@(bin) processBinAfterMerging(bin), binData, 'UniformOutput', false);
    binData = [binData{:}];
    
    % Extract final binned data
    binnedX = [binData.binCenter];
    meanY = [binData.meanY];
    binnedY = {binData.binnedY};
    semY = [binData.sems];
    counts = [binData.counts];
    
    % Filter out entries with NaN meanY or empty binnedY
    validIndices = ~isnan(meanY) & ~cellfun(@isempty, binnedY);
    binnedX = binnedX(validIndices);
    meanY = meanY(validIndices);
    binnedY = binnedY(validIndices);
    semY = semY(validIndices);
    if size(x,2) > 1
        counts = counts(validIndices);
    else
        counts = counts(validIndices)*10;
    end
    % Plotting remains unchanged...
    plotData(binnedX, meanY, semY, y, lineColor, markerType, counts);
end

function bin = initBin(x, y, binMin, binMax)
    binIndices = x >= binMin & x < binMax;
    bin.binnedY = y(binIndices); % Store raw Y values directly
    bin.binCenter = (binMin + binMax) / 2;
end

function binData = mergeBinData(binData)
    
    i = 1;
    while i < length(binData)
        if binData(i).binCenter == 0 || numel((binData(i).binnedY)) >= 10
            i = i + 1;
            continue; % Skip bin centered at 0 and bins with sufficient counts
        end

        % Check if it's not the last bin and merge if the current bin has fewer than 10 counts
        binData(i) = mergeBins(binData(i), binData(i+1));
        binData(i+1) = []; % Remove the merged bin
        %{
        if i < length(binData) - 1
            binData(i) = mergeBins(binData(i), binData(i+1));
            binData(i+1) = []; % Remove the merged bin
        else
            break; % If it's the last bin, stop the process
        end
        %}
        % Move to the next bin if we've reached or exceeded the target counts
        if numel((binData(i).binnedY)) >= 10
            i = i + 1;
        end
    end
end

function mergedBin = mergeBins(bin1, bin2)
    mergedBin.binCenter = mean([bin1.binCenter, bin2.binCenter]);
    mergedBin.binnedY = [bin1.binnedY; bin2.binnedY]; % Combine raw Y values from both bins
end

function bin = processBinAfterMerging(bin)
    if ~isempty(bin.binnedY)
        binValues = bin.binnedY(~isnan(bin.binnedY));
        bin.meanY = mean(binValues, 'omitnan');
        bin.sems = std(binValues, 'omitnan') / sqrt(length(binValues));
        bin.counts = numel(binValues);
    else
        bin.meanY = NaN;
        bin.sems = NaN;
        bin.counts = 0;
    end
end

% Define the plotting subfunction
function plotData(binnedX, meanY, semY, y, lineColor, markerType, counts)
    hold on; % Ensure the plot stays for adding text
    
    if strcmp(lineColor,'k') || isequal(lineColor,[0 0 0])
        markerFaceColor = [1 1 1];
    else
        markerFaceColor = lineColor;
    end
    
    if size(y,2) > 1 % if no sem, don't shade
        patchSaturationVal = 0.2;
    else
        patchSaturationVal = 0;
    end
    
    markerSize = 20;
    shadedErrorBar(binnedX, meanY, semY, 'patchSaturation', patchSaturationVal, 'lineprops', ...
                   {'Color', lineColor, 'LineStyle', 'none', 'LineWidth', 2, 'Marker', markerType, ...
                    'MarkerFaceColor', markerFaceColor, 'MarkerEdgeColor', 'k', 'MarkerSize', markerSize});
    
    % Assuming markerFaceColor is defined, e.g., markerFaceColor = [0.2, 0.5, 0.7];
    % Convert markerFaceColor to grayscale using luminance to decide on text color
    annotateDataPoints(binnedX, meanY, counts, markerFaceColor);

    hold off; % Release the plot hold
    
    % Cosmetic adjustments
    upFontSize(24, 0.01);
end
