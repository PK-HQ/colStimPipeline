function [camSizeMm, projSizePx] = initializeCamProjDims(imagingData, bitmapData, blockID)
camSizePx = imagingData.pixels(blockID);
camPxSizeMm = imagingData.pixelsizemm(blockID);
camSizeMm = [camSizePx * camPxSizeMm, camSizePx * camPxSizeMm]; % [X, Y]
projSizePx = [bitmapData.projectorx, bitmapData.projectory, bitmapData.pixelsizemm];
end