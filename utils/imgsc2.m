function imgsc2(varargin)
    img = double(varargin{1});
    if size(varargin, 2) > 1
        titleStr = varargin{2};
    else
        titleStr = '';
    end

    % Define colormap and NaN color
    nanColor = [0.75 0.75 0.75];  % Grey color for NaNs, as an RGB triplet
    [colormapRB, colormapR, colormapB] = fireice; % Assuming fireice is defined

    % Determine min and max of the non-NaN data
    cmin = min(img(:), [], 'omitnan');
    cmax = max(img(:), [], 'omitnan');

    if cmax <= cmin
        cmin = -1;
        cmax = 1;
    end

    % Replace NaNs with a value *outside* the main data range, but included in caxis
    imgWithNaNs = img;
    nanValue = cmin - (cmax - cmin) * 0.1; % 10% below cmin
    imgWithNaNs(isnan(img)) = nanValue;

    % Plot the image
    imagesc(imgWithNaNs);

    % Set the color axis to include the NaN value AND the data range
    caxis([nanValue, cmax]);

    % Define a custom colormap. Add nanColor as the FIRST entry.
    customColormap = [nanColor; colormapRB];

    % Apply the custom colormap
    colormap(customColormap);

    % Add a colorbar (important for interpretation)
    colorbar;

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