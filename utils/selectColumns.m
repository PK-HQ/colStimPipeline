function [bitmapCamSpace, gaussianMask] = selectColumns(bitmapData, imagingData, bitmapNo, blockID, camspaceBitmap)
    % Check if nColumnsWanted is specified and not empty
    if ~isempty(bitmapData.nColumnsWanted)
        nColumnsWanted = bitmapData.nColumnsWanted(1, bitmapNo, blockID);
        contourLevel = bitmapData.gaussianContourLevel(1, bitmapNo, blockID);
        disp(['Columns wanted: ' num2str(nColumnsWanted)]);

        if nColumnsWanted > 0
            % For single or multiple columns, but not zero columns
            if contourLevel == 200 && nColumnsWanted == 1
                % Specific case for single column with contour level 200
                disp('Single column with specific contour level');
                [bitmapCamSpace, gaussianMask] = processSingleColumnSpecialCase(camspaceBitmap, imagingData.centerCoords);
            else
                % General case for single or multiple columns
                disp(['Processing for ' num2str(nColumnsWanted) ' columns']);
                gaussianMask = getGaussianMask(imagingData, blockID, contourLevel, nColumnsWanted);
                bitmapCamSpace = gaussianMask .* camspaceBitmap;
            end
        elseif nColumnsWanted == 0 || isnan(nColumnsWanted)
            % Handle cases with zero columns or unspecified number of columns
            disp('No specific columns wanted or number of columns is NaN');
            gaussianMask = getDefaultMask(imagingData, nColumnsWanted);
            bitmapCamSpace = gaussianMask .* camspaceBitmap;
        end
    else
        disp('nColumnsWanted is not specified.');
        bitmapCamSpace = []; % Return empty if conditions are not met
    end
end

function [bitmapCamSpace, gaussianMask] = processSingleColumnSpecialCase(camspaceBitmap, centerCoords)
    L = bwlabel(camspaceBitmap);
    S = regionprops(L, 'centroid');
    C = vertcat(S.Centroid); % Centers of all labeled regions
    c0 = centerCoords; % Center of image
    DS = sum((C - c0(1, :)).^2, 2); % Squared Euclidean distance
    [~, idx] = min(DS); % Find the closest centroid
    bitmapCamSpace = double(L == idx);
    gaussianMask = bitmapCamSpace; % Use the selected region as the mask
end

function gaussianMask = getGaussianMask(imagingData, blockID, contourLevel, nColumnsWanted)
    % Placeholder for generating a gaussian mask based on contour level and number of columns wanted
    % Adjust the logic as needed to match your specific requirements
    if contourLevel < 200 || nColumnsWanted > 1
        gaussianMask = imagingData.gaussfit(1, contourLevel, blockID).area;
    else
        gaussianMask = ones(imagingData.pixels(1), imagingData.pixels(1)); % Example default mask
    end
end

function defaultMask = getDefaultMask(imagingData, nColumnsWanted)
    % Generate a default mask based on whether columns are wanted or not
    if nColumnsWanted == 0 || isnan(nColumnsWanted)
        defaultMask = zeros(imagingData.pixels(1), imagingData.pixels(1));
    else
        defaultMask = ones(imagingData.pixels(1), imagingData.pixels(1));
    end
end
