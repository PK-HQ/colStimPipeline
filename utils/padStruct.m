function paddedStruct = padStruct(originalStruct, targetSize)
    % Determine the current size of the struct
    currentSize = numel(originalStruct);
    
    % Calculate how many new elements need to be added
    elementsToAdd = targetSize - currentSize;
    
    % Check if we need to add elements
    if elementsToAdd <= 0
        % If the original struct is already the target size or larger, return it directly
            paddedStruct = originalStruct;
        return;
    end
    
    % Create a template for new elements based on the existing struct fields
    % Assuming the struct has at least one element and all elements have the same fields
    fieldNames ={'threshPercentile','area'}; %fieldnames(originalStruct);
    newElementTemplate = cell2struct(cell(length(fieldNames), 1), fieldNames);
    
    % Initialize an array of new elements to add
    newElements(elementsToAdd) = newElementTemplate; % Preallocate for efficiency
    
    % Fill each new element with default values (e.g., empty arrays, NaNs, or zeros)
    for i = 1:elementsToAdd
        for f = 1:length(fieldNames)
            fieldName = fieldNames{f};
            % Assuming numeric values for simplicity; adjust based on actual field types
            newElements(i).(fieldName) = NaN; % Or use NaN, 0, or '' for strings, etc.
        end
    end
    
    % Append the new elements to the original struct
    paddedStruct = [originalStruct, newElements];
end
