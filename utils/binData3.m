function [xOptoAvg, yOptoAvg] = binData3(xOptoAll, yOptoAll)
    % Define the bins
    binCenters = -100:5:100;  % Bin centers from -100 to 100
    binEdges = (-102.5):5:102.5;  % Bin edges extending 2.5 units on each side of bin centers

    % Initialize output arrays
    xOptoAvg = binCenters;  % Use bin centers for the x-axis values
    yOptoAvg = NaN(1, numel(binCenters));  % Initialize yOptoAvg with NaNs

    % Process each bin
    for i = 1:length(binEdges)-1
        % Find indices of xOptoAll within the current bin range
        binIndices = xOptoAll >= binEdges(i) & xOptoAll < binEdges(i+1);

        % Average yOptoAll values corresponding to the binIndices
        % Calculate mean only for non-NaN entries
        meanValues = arrayfun(@(row) mean(yOptoAll(row, binIndices(row, :)), 'omitnan'), 1:size(yOptoAll, 1));

        % Check if there are any values to average
        validMeans = meanValues(~isnan(meanValues));
        if ~isempty(validMeans)
            yOptoAvg(i) = mean(validMeans);
        else
            yOptoAvg(i) = NaN;  % Maintain NaN if no valid entries are found
        end
    end

    % Remove NaN values and their corresponding x-axis values
    validIdx = ~isnan(yOptoAvg);
    xOptoAvg = xOptoAvg(validIdx);
    yOptoAvg = yOptoAvg(validIdx);
    
    % Special processing for zero
    if sum(xOptoAll == 0) > 0
        % Additional processing for data exactly at zero can be added here
    end
end
