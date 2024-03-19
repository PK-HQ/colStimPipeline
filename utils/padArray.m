function paddedArray = padArray(inputArray, targetLength, dim, padValue)
    if isempty(inputArray)
        paddedArray=inputArray;
    else
        % inputArray: The N-dimensional array to pad
        % targetLength: The desired length of the array along dimension 'dim'
        % dim: The dimension along which to pad
        % padValue: The value to use for padding
        
        % Get the size of the input array along the specified dimension
        currentLength = size(inputArray, dim);
        
        % Calculate the total padding required
        totalPadding = targetLength - currentLength;
        
        if totalPadding <= 0
            % No padding needed, return the original array
            paddedArray = inputArray;
            return;
        end
        
        % Calculate padding to add before and after
        % For simplicity, this example adds all padding at the end
        paddingBefore = 0;
        paddingAfter = totalPadding;
        
        % Create padding arrays
        paddingSize = size(inputArray);
        paddingSize(dim) = paddingBefore;
        paddingBeforeArray = repmat(padValue, paddingSize);
        
        paddingSize(dim) = paddingAfter;
        paddingAfterArray = repmat(padValue, paddingSize);
        
        % Concatenate the original array with the padding arrays
        paddedArray = cat(dim, paddingBeforeArray, inputArray, paddingAfterArray);

    end
end
