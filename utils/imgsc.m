function imgsc(varargin)
    img = double(varargin{1});
    if size(varargin, 2) > 1
        titleStr = varargin{2};
    else
        titleStr = '';
    end

    % Define colormap and NaN color
    nanColor = 0.75;  % Grey color for NaNs
    [colormapRB, colormapR, colormapB] = fireice;
    
    % Determine min and max of the non-NaN data
    cmin = min(img(:), [], 'omitnan');
    cmax = max(img(:), [], 'omitnan');

    if cmax<=cmin
        cmin=-1;
        cmax=1;
    end

    % Replace NaNs with a value below the data range (e.g., cmin - 1)
    imgWithNaNs = img;
    imgWithNaNs(isnan(img)) = cmin - 1;

    % Plot the image
    imagesc(imgWithNaNs);
    caxis([cmin, cmax]);  % Set color axis to include the NaN placeholder value

    % Define a custom colormap that includes grey for NaN values
    nanColorRGB = [1 1 1];  % RGB values for nan
    customColormap = [nanColorRGB; colormapRB];  % Add grey as the first color

    % Apply the custom colormap
    colormap(customColormap);

    % Adjust the axis properties
    if size(img, 1) == size(img, 2)
        axis square
    else
        axis image
    end
    
    % Add title with correct interpreter
    if contains(titleStr, '\') || contains(titleStr, '{')
        title(titleStr, 'FontWeight', 'Normal', 'Interpreter', 'tex');
    else
        title(titleStr, 'FontWeight', 'Normal');
    end

    % Set the figure background color
    set(gcf, 'color', 'w');
end
