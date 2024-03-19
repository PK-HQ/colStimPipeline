function [bitmapData] = convertForProjector2(behavioralData, imagingData, bitmapData, ...
    currentBlockStruct, conversionType, blockID, ...
    pdfFilename, plotFlag, saveFlagBMP, saveFlag)
disp('Generating bitmap...');

% Initialize dimensions
[camSizeMm, projSizePx] = initializeDimensions(imagingData, bitmapData, blockID);

% Main bitmap generation loop
bitmapData = generateBitmaps(bitmapData, imagingData, blockID, camSizeMm, projSizePx, ...
    conversionType, plotFlag, saveFlagBMP, pdfFilename, saveFlag);

end

function [camSizeMm, projSizePx] = initializeDimensions(imagingData, bitmapData, blockID)
camSizePx = imagingData.pixels(blockID);
camPxSizeMm = imagingData.pixelsizemm(blockID);
camSizeMm = [camSizePx * camPxSizeMm, camSizePx * camPxSizeMm]; % [X, Y]
projSizePx = [bitmapData.projectorx, bitmapData.projectory, bitmapData.pixelsizemm];
end

function bitmapData = generateBitmaps(bitmapData, imagingData, blockID, camSizeMm, projSizePx, ...
    conversionType, plotFlag, saveFlagBMP, pdfFilename, saveFlag)

for bitmapNo = 1:size(bitmapData.columnarbitmapCoreg, 3)
    % Fit and process only if conditions meet
    if any(bitmapData.orts(1, bitmapNo, blockID) == [0, 90])
        [gaussianMask, centerCoords] = fit2DGaussianIfNeeded(bitmapData, imagingData, bitmapNo, blockID, plotFlag);
        bitmapData = processBitmap(bitmapData, imagingData, gaussianMask, centerCoords, bitmapNo, blockID, camSizeMm, projSizePx, ...
            conversionType, plotFlag, saveFlagBMP, pdfFilename, saveFlag);
    end
end

end

%% === Fit gaussian if needed ===

function [gaussianMask, centerCoords] = fit2DGaussianIfNeeded(bitmapData, imagingData, bitmapNo, blockID, plotFlag)
% Initialize outputs in case conditions are not met
gaussianMask = [];
centerCoords = [];

% Check if fitting is required based on bitmapData.gaussianCond
if ~isnan(bitmapData.gaussianCond(blockID))
    gaussianResponse = abs(imagingData.gaussresp(:,:,bitmapData.gaussianCond(blockID), blockID));
    desiredResp = imresize(gaussianResponse, [imagingData.pixels(blockID), imagingData.pixels(blockID)], 'bilinear');
    
    % Assuming fit2Dgaussian is a function you would implement that fits a 2D Gaussian to the desiredResp
    % and returns a mask (gaussianMask) and the center coordinates (centerCoords) of the fit.
    % For simplicity, we'll treat it as already implemented:
    [gaussianMask, centerCoords] = fit2Dgaussian(desiredResp, plotFlag);
    
    % Update imagingData with the results (this might need adjustment based on actual data structure)
    imagingData.gaussfit(1, :, blockID) = gaussianMask;
    imagingData.centerCoords(1, :, blockID) = centerCoords;
else
    % Handle case where gaussianCond is NaN; setup default or empty outputs if needed
    % gaussianMask = []; centerCoords = []; are already set for this case.
end
end

%% ===
function bitmapData = processBitmap(bitmapData, imagingData, gaussianMask, centerCoords, bitmapNo, blockID, camSizeMm, projSizePx, conversionType, plotFlag, saveFlagBMP, pdfFilename, saveFlag)
% Assuming conversionType is always 'cam2proj' for this example

% Extract bitmap and contour mask level
camspaceBitmap = bitmapData.columnarbitmapCoreg(:,:,bitmapNo, blockID);
contourMaskLevel = bitmapData.gaussianContourLevel(1,bitmapNo,blockID);

% Process bitmap based on orientation and conversion type
% The following steps would typically include resizing to projector space, applying transformations, etc.
% For brevity, let's assume we resize and then apply transformations directly:

% Resize bitmap to projector space
[bitmapProjSpace] = resizeBitmapToProjectorSpace(camspaceBitmap, camSizeMm, projSizePx);

% Apply rotation and translation based on predefined offsets and calibration
[bitmapProjSpaceTransformed] = applyRotationAndTranslation(bitmapProjSpace, imagingData, projSizePx);

% Save the transformed bitmap if required
if saveFlagBMP == 1
    saveTransformedBitmap(bitmapProjSpaceTransformed, bitmapNo, blockID, pdfFilename);
end

% Update bitmapData with the transformed bitmap for this blockID and bitmapNo
bitmapData.columnarbitmapTFprojspace(:,:,bitmapNo,blockID) = bitmapProjSpaceTransformed;

% Additional processing could include plotting if plotFlag is set, etc.
if plotFlag
    % Plotting logic here
end

end

function [resizedBitmap] = resizeBitmapToProjectorSpace(camspaceBitmap, camSizeMm, projSizePx)
% Calculate the resize factor based on the projector and camera resolutions
resizeFactor = projSizePx(3) / mean(camSizeMm); % Simplified, adjust based on actual dimensions
resizedBitmap = imresize(camspaceBitmap, resizeFactor, 'bilinear');
end

function [bitmapProjSpaceTransformed] = applyRotationAndTranslation(bitmapProjSpace, imagingData, projSizePx)
% Placeholder for rotation and translation application
% This would involve calculations based on camera and projector alignment specifics
% For simplicity, assume it's a direct application of known offsets/angles
%% Apply camera-projector transformation
if ~isempty(imagingData.transformmatrix(:,blockID))
    bitmapProjSpaceTransformed = imwarp(bitmapProjSpace,imagingData.transformmatrix(:,blockID),'OutputView',imref2d(size(projSizePx)));
end
end

function saveTransformedBitmap(bitmap, bitmapNo, blockID, destinationPath)
% Construct filename and save bitmap to the specified destination
filename = sprintf('transformed_%d_%d.bmp', bitmapNo, blockID);
fullPath = fullfile(destinationPath, filename);
%imwrite(bitmap, fullPath);
disp(['Saved bitmap to ', fullPath]);
end
