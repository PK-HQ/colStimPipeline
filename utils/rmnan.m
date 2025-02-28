function [array, idx] = rmnan(array)
    % Initialize a cell array to temporarily hold cleaned rows
    cleanedRows = cell(size(array, 1), 1);
    % Initialize a cell array to store the indices of non-NaN elements per row
    indices = cell(size(array, 1), 1);
    
    % Loop through each row to remove NaNs and store indices
    for i = 1:size(array, 1)
        nonNanIndices = ~isnan(array(i, :));
        cleanedRows{i} = array(i, nonNanIndices);
        indices{i} = find(nonNanIndices);
    end
    
    % Determine the maximum length of the cleaned rows
    maxLength = max(cellfun(@length, cleanedRows));
    
    % Initialize the final array with NaNs
    cleanedArray = NaN(size(array, 1), maxLength);
    % Initialize the indices array with zeros
    idx = zeros(size(array, 1), maxLength);
    
    % Fill the final array with cleaned data and indices
    for i = 1:size(array, 1)
        numElements = length(cleanedRows{i});
        cleanedArray(i, 1:numElements) = cleanedRows{i};
        idx(i, 1:numElements) = indices{i};
    end
    
    % Return the cleaned and padded array and the indices
    array = cleanedArray;
end
