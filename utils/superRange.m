function [minVal, maxVal] = superRange(image)
%GETIMAGERANGE Gets the minimum, maximum, and range of values in an image.
%
%   [minVal, maxVal, imageRange] = getImageRange(image)
%
%   Inputs:
%       image - The input image (can be any numeric type: uint8, int16, double, etc.).
%
%   Outputs:
%       minVal    - The minimum pixel value in the image.
%       maxVal    - The maximum pixel value in the image.
%       imageRange - The range of pixel values (maxVal - minVal).

    % Input validation
    if ~isnumeric(image)
        error('Input image must be a numeric array.');
    end

    % Handle NaNs and Infs
    image = image(:);  % Convert to a column vector for easier handling
    image = image(~isinf(image) & ~isnan(image));  % Remove NaNs and Infs
    
    if isempty(image)
        minVal = [];
        maxVal = [];
        imageRange = [];
        warning('Input image contains only NaN or Inf values. Returning empty outputs.');
        return;  % Important: Return early if the image is all NaN/Inf
    end
   

    % Get the minimum and maximum values
    minVal = min(image);
    maxVal = max(image);

    % Calculate the range
    imageRange = maxVal - minVal;

end