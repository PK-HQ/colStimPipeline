function [powerFiltered, powerFilteredMask] = powerFiltering(power)
    % Check the size of the input
    [nBlocks, cols] = size(power);
    if cols ~= 2
        error('The input matrix power must be of size nBlocks x 2.');
    end
    
    % Initialize the filtered power array
    powerFiltered = power;
    
    % Loop through each column to detect and replace outliers
    for j = 1:cols
        colData = power(:, j);
        medianVal = median(colData);
        lowerBound = 0.75 * medianVal; % 75% of the median
        upperBound = 1.25 * medianVal; % 125% of the median
        
        % Find outliers
        outlierIdx = colData < lowerBound | colData > upperBound;
        
        % Replace outliers with NaN
        powerFiltered(outlierIdx, j) = NaN;
    end
    
    % Create the powerFilteredMask
    powerFilteredMask = ones(nBlocks, 1); % Initialize mask
    for i = 1:nBlocks
        if any(isnan(powerFiltered(i, :)))
            powerFilteredMask(i) = NaN;
        end
    end
end
