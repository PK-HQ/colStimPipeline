function [xOptoAvg, yOptoAvg] = binData2(xOptoAll, yOptoAll)
    % Define the bins
    %binEdges = -100:5:100;  % Bin edges from -100 to 100 with a width of 5
    %binCenters = -97.5:5:97.5;  % Bin centers from -97.5 to 97.5
    % Define the bins
    binEdges = [-100, -95:5:-5, 0, 0, 5:5:95, 100];  % Include special bin edge 0 to 0
    binCenters = [-97.5:5:-2.5, 0, 2.5:5:97.5];  % Include special bin center at 0
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
    validIdx=~isnan(yOptoAvg);
    xOptoAvg=xOptoAvg(validIdx);
    yOptoAvg=yOptoAvg(validIdx);
    
    if sum(xOptoAll==0)>0
        1;
    end
    
end
