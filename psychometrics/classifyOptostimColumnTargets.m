function [experimentalMask, controlMask] = classifyOptostimColumnTargets(orts, nBlocks)
% Classify blocks by the pair of stimulated orientation columns.

    orientationPairs = squeeze(orts);
    if nBlocks == 1 && isvector(orientationPairs)
        orientationPairs = orientationPairs(:);
    elseif size(orientationPairs, 2) == nBlocks
        % Expected layout: orientations by block.
    elseif size(orientationPairs, 1) == nBlocks
        orientationPairs = orientationPairs';
    else
        error('bitmapData.orts cannot be aligned to %d blocks.', nBlocks);
    end

    if size(orientationPairs, 1) ~= 2
        error('Each block in bitmapData.orts must contain two orientations.');
    end

    orientationPairs = sort(orientationPairs, 1);
    tolerance = 1e-6;
    experimentalMask = all(abs(orientationPairs - [0; 90]) <= tolerance, 1)';
    controlMask = all(abs(orientationPairs - [45; 135]) <= tolerance, 1)';
end
