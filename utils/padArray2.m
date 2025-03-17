function outputArray = padArray2(inputArray, targetLength, dim, padValue, flag)
    if isempty(inputArray)
        outputArray = inputArray;
    else
        switch flag
            case 'pad'
                % Padding logic
                currentLength = size(inputArray, dim);
                totalPadding = targetLength - currentLength;

                if totalPadding <= 0
                    % No padding needed, return the original array
                    outputArray = inputArray;
                else
                    % Add padding at the end for simplicity
                    paddingSize = size(inputArray);
                    paddingSize(dim) = totalPadding;
                    paddingArray = repmat(padValue, paddingSize);
                    outputArray = cat(dim, inputArray, paddingArray);
                end

            case 'rmpad'
                % Remove padding logic
                currentLength = size(inputArray, dim);
                indices = repmat({':'}, 1, ndims(inputArray));
                isPadSlice = false(1, currentLength);

                for i = 1:currentLength
                    indices{dim} = i;
                    % Check if entire slice consists only of NaNs (for NaN padding)
                    % Use `isnan` for NaN checks, as `== NaN` is not valid
                    if padValue == NaN
                        isPadSlice(i) = all(all(isnan(inputArray(indices{:}))));
                    else
                        isPadSlice(i) = all(all(inputArray(indices{:}) == padValue));
                    end
                end

                % Remove padding from the beginning and end of the dimension
                firstNonPad = find(~isPadSlice, 1, 'first');
                lastNonPad = find(~isPadSlice, 1, 'last');
                
                if isempty(firstNonPad) || isempty(lastNonPad)
                    outputArray = []; % All slices are padding
                else
                    indices{dim} = firstNonPad:lastNonPad;
                    outputArray = inputArray(indices{:});
                end

            otherwise
                error('Unsupported flag. Use ''pad'' or ''rmpad''.');
        end
    end
end
