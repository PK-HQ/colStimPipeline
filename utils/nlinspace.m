function spacedValues = nlinspace(startValue, endValue, numPoints, spacingType)
    % nlinspace: Generates numbers between startValue and endValue with
    % linearly or nonlinearly increasing spacing.
    %
    % Input:
    %   startValue  - Starting value of the range
    %   endValue    - Ending value of the range
    %   numPoints   - Number of points in the range
    %   spacingType - Type of spacing ('linear' or 'nonlinear')
    %
    % Output:
    %   spacedValues - Array of spaced numbers
    
    if strcmp(spacingType, 'linear')
        % Linearly increasing spacing
        linSpace = linspace(1, 10, numPoints); % Linearly increasing vector
        linSpace = linSpace / sum(linSpace); % Normalize to sum to 1
        spacedValues = cumsum(linSpace); % Cumulative sum to get linearly increasing spacing
    elseif strcmp(spacingType, 'nonlinear')
        % Nonlinearly increasing spacing (cubic)
        expSpace = linspace(1, 10, numPoints).^3; % Cubic increase
        expSpace = expSpace / sum(expSpace); % Normalize to sum to 1
        spacedValues = cumsum(expSpace); % Cumulative sum to get nonlinearly increasing spacing
    else
        error('Invalid spacing type. Choose ''linear'' or ''nonlinear''.');
    end
    
    % Scale the values to the specified range
    spacedValues = startValue + (endValue - startValue) * spacedValues / spacedValues(end);
end