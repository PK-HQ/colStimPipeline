function [binnedX, meanY, binnedY] = plotSEM(x, y, lineColor, markerType)
    bins = -2.5:5:102.5; % Define bins
    % Adjust binData structure to include raw Y values
    binData = arrayfun(@(binMin, binMax) processBin(x, y, binMin, binMax), bins(1:end-1), bins(2:end), 'UniformOutput', false);
    binData = [binData{:}];
    
    % Skip merging for the bin centered at 0 and ensure merging left to right with target counts
    i = 1;
    while i < length(binData)
        if binData(i).binCenter == 0 || binData(i).counts >= 10
            i = i + 1;
            continue; % Skip to the next bin
        end

        % Find the best neighbor to merge based on counts and proximity
        if i < length(binData) - 1 % Check if not the last bin
            % Merge with the next bin if the current bin has fewer than 10 counts
            binData(i) = mergeBins(binData(i), binData(i+1));
            binData(i+1) = []; % Remove the merged bin
        else
            break; % Stop if we reached the last bin
        end
        
        % If after merging we've reached or exceeded the target counts, move to the next bin
        if binData(i).counts >= 10
            i = i + 1;
        end
    end
    
    % Extract final binned data
    binnedX = [binData.binCenter];
    meanY = [binData.meanY];
    binnedY = {binData.binnedY}; % Adjusted for the cell array of raw Y values per bin
    semY = [binData.sems];
    counts = [binData.counts];
    
    % Filtering entries where meanY is NaN or binnedY is empty
    validIndices = ~isnan(meanY) & ~cellfun(@isempty, binnedY);
    binnedX = binnedX(validIndices);
    meanY = meanY(validIndices);
    binnedY = binnedY(validIndices);
    semY = semY(validIndices);
    counts = counts(validIndices);
    
    % Plotting
    hold on; % Ensure the plot stays for adding text
    % Plotting
    if strcmp(lineColor,'k') || isequal(lineColor,[0 0 0])
        markerFaceColor=[1 1 1];
    else
        markerFaceColor=lineColor;
    end
    if size(y,2)>1 % if no sem, don't shade
        patchSaturationVal=0.2;
    else
        patchSaturationVal=0;
    end
    markerSize=20;
    shadedErrorBar(binnedX, meanY, semY, 'patchSaturation',  patchSaturationVal, 'lineprops',{'Color', lineColor,'LineStyle','none','LineWidth',2,...
        'Marker',markerType,'MarkerFaceColor',markerFaceColor,'MarkerEdgeColor','k','MarkerSize',markerSize});
    
    %% Add count of datapoints depending on luminance of datapoint color
    % Assuming markerFaceColor is defined, e.g., markerFaceColor = [0.2, 0.5, 0.7];
    % Convert markerFaceColor to grayscale using luminance to decide on text color
    annotateDataPoints(binnedX, meanY, counts, markerFaceColor)

    hold off; % Release the plot hold

    % Cosmetic
    upFontSize(24,0.01)    
end

function bin = processBin(x, y, binMin, binMax)
    binIndices = x >= binMin & x < binMax;
    binValues = y(binIndices);
    binnedY = binValues; % Capture raw Y values before filtering NaNs
    binValues = binValues(~isnan(binValues));
    bin.binCenter = (binMin + binMax) / 2;
    bin.counts = numel(binValues);
    bin.binnedY = binnedY; % Store raw Y values in the bin structure
    if ~isempty(binValues)
        bin.meanY = mean(binValues,'omitnan');
        bin.sems = std(binValues,'omitnan') / sqrt(numel(binValues));
        binValues'
        bin.sems'
    else
        bin.meanY = NaN;
        bin.sems = NaN;
        binValues'
        bin.sems'
    end
end

function mergedBin = mergeBins(bin1, bin2)
    % Adjusting the mergeBins function to correctly handle the merging of bins
    % including the 'binnedY' field to ensure structure consistency.

    % Calculate the mean center
    mergedBin.binCenter = mean([bin1.binCenter, bin2.binCenter]);
    
    % Combine the meanY values while excluding NaNs
    combinedMeanY = [bin1.meanY * ones(1, bin1.counts), bin2.meanY * ones(1, bin2.counts)];
    combinedMeanY = combinedMeanY(~isnan(combinedMeanY));
    
    % Calculate the new meanY and SEM
    mergedBin.meanY = mean(combinedMeanY, 'omitnan');
    mergedBin.counts = numel(combinedMeanY); 
    
    % Combine the raw Y values (binnedY) from both bins
    mergedBin.binnedY = [bin1.binnedY; bin2.binnedY];
    
    if ~isempty(combinedMeanY)
        % Calculate SEM based on combined meanY values
        mergedBin.sems = std(combinedMeanY, 'omitnan') / sqrt(numel(combinedMeanY));
    else
        mergedBin.meanY = NaN;
        mergedBin.sems = NaN;
    end
end
