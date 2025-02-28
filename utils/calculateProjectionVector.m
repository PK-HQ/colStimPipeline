function [angle, amplitude, similarityRef, similarityData, angleAverage, amplitudeAverage]=...
    calculateProjectionVector(bitmapData, imagingData, images, referenceImage, blockID,...
    condIDs, orientation, snrMask, analysisFlag)
%% Load ort map, but separated into each color
% init
nContrasts=size(images,3);
orientations=repmat(orientation,1,nContrasts);
angle=nan(1, nContrasts);
amplitude=nan(1, nContrasts);
x=nan(1, nContrasts);
y=nan(1, nContrasts);
similarityRef=[];
similarityData=[];
angleAverage=[];
amplitudeAverage=[];
switch analysisFlag
    case 'similarity'
        %% Calculate angle and amplitude of image projected into reference image space
        referenceImage=rescale(referenceImage(:,:,2),-1,1)-rescale(referenceImage(:,:,1),-1,1);
        % Initialize the rescaled image ([-1 1])
        referenceImageRescaled = zeros(size(referenceImage));
        
        % Separate and rescale the positive values
        positivePart = referenceImage > 0; % Logical mask for positive values
        if any(positivePart, 'all')
            maxPositive = max(referenceImage(positivePart), [], 'all'); % Maximum positive value
            referenceImageRescaled(positivePart) = referenceImage(positivePart) / maxPositive; % Rescale to [0, 1]
        end
        
        % Separate and rescale the negative values
        negativePart = referenceImage < 0; % Logical mask for negative values
        if any(negativePart, 'all')
            minNegative = min(referenceImage(negativePart), [], 'all'); % Minimum negative value
            referenceImageRescaled(negativePart) = referenceImage(negativePart) / abs(minNegative); % Rescale to [-1, 0]
        end
        referenceImageRescaledUnmasked=referenceImageRescaled;
        referenceImageRescaledUnmasked(isnan(snrMask))=0;
        
        %% Coregister by linear translation with crosscorr
        imagesHighCon=images(:,:,condIDs.V90O90(end))-images(:,:,condIDs.V0O0(end));
        plotFlag=0;
        [imagesCoreg, tformParams] = coregisterImages(referenceImageRescaled, snrMask, imagesHighCon, images, plotFlag);
        %coregisterImages(images, referenceImageRescaledUnmasked, snrMask, condIDs, 'crosscorr');
        
        
        %% Calculate projection values per image
        for contrast=1:nContrasts
            % Intensity weighted Pearson
            [rWeighted,correlationImage,correlationImageRef] = weightedCorrelation(referenceImageRescaled, imagesCoreg(:,:,contrast));
            
            % Get angle and amplitude of vector
            angle(contrast)=orientations(contrast); amplitude(contrast)=rWeighted*100; 
            
            % Get slice image
            similarityData.Img(:,:,contrast)=correlationImage*100;
            [similarityData.sliceImg(:,:,contrast), similarityData.sliceAmp(:,:,contrast)] = getImageSlice(similarityData.Img(:,:,contrast), 250, 10);
        
            % Convert to 2D vector coordinates
            x(contrast)=cos(angle(contrast)) * amplitude(contrast); %x
            y(contrast)=sin(angle(contrast)) * amplitude(contrast); %y
        
            %plotSubplots(imagesCoreg, projectedMap, referenceImageRescaled, contrast, nContrasts);
        
        end
        [similarityRef.sliceImg, similarityRef.sliceAmp] = getImageSlice(correlationImageRef, 250, 20);
        similarityRef.Img=correlationImageRef;
        % Get vector average (angle and amplitude)
        angleAverage=atan(mean(y)/mean(x));
        amplitudeAverage=mean(y)/sin(angleAverage);
    case 'similarityColumn'
        % Bitmap mask for removing stimulated pixels
        roiMask=double(~bitmapData.bitmapMask(:,:,blockID));
        roiMask(roiMask==0)=nan;            
        
        % Create reference map of 0- and 90-tuned columns, with snr mask and bitmap mask
        tileGridSize = [32,32];
        gammaVal = 5;
        sensitivity = 1e-3;  % Example sensitivity value
    
        [~, referenceImageColumn0, ~, ~] = processActivityMap(...
            rescale(imagingData.refmapCoreg(:,:,1,blockID), -1, 1), tileGridSize, gammaVal, sensitivity);
        [~, referenceImageColumn90, ~, ~] = processActivityMap(...
            rescale(imagingData.refmapCoreg(:,:,2,blockID), -1, 1), tileGridSize, gammaVal, sensitivity);
        
        referenceImageColumn0 = -referenceImageColumn0 .* snrMask .* roiMask;
        referenceImageColumn90 = referenceImageColumn90 .* snrMask .* roiMask;

        % Initialize result vectors
        numPairs = numel(condIDs.V0O0);  % Assuming condIDs.V0O0 and condIDs.V90O90 have the same length
        dataO0C0 = nan(numPairs,1);
        dataO0C90 = nan(numPairs,1);
        dataO90C0 = nan(numPairs,1);
        dataO90C90 = nan(numPairs,1);
        
        % Loop through each contrast, 1st row = 0 column and 2nd row = 90 column
        for contrast=1:nContrasts
            % Define data for current iteration
            resp = images(:,:,contrast) .* snrMask;
            respROI = resp .* roiMask;
        
            % Compute responses for each image, separately for 0 and 90 columns
            resp0 = respROI .* referenceImageColumn0 .* roiMask; resp0(resp0==0) = nan;
            resp90 = respROI .* referenceImageColumn90 .* roiMask; resp90(resp90==0) = nan;
        
        % Save response, 1st-row = 0 columns and 2nd-row = 90 columns
            amplitude(1, contrast) = mean(removenan(resp0));
            amplitude(2,contrast)= mean(removenan(resp90));
        end

end
end

function [imagesCoreg, tformParams] = coregisterImages(...
    refMap, snrMask, diffImage, dataStack, plotFlag)
%COREGBYTRANSLATIONIMREGTFORM
%   Estimate a pure translation using imregtform (multimodal, Mattes MI) 
%   to align 'diffImage' (e.g., 90-0 deg difference) to 'refMap', then
%   apply that translation to 'snrMask' and each slice of 'dataStack'.
%
%   USAGE:
%       [imagesCoreg, tformParams] = coregByTranslationImregtform(...
%           refMap, snrMask, diffImage, dataStack, plotFlag)
%
%   INPUTS:
%       refMap    : [NxM], reference map (normalized [-1, +1]) 
%                   (NaNs possible outside valid SNR region)
%       snrMask   : [NxM], logical or double mask (1 = inside ROI)
%       diffImage : [NxM], difference image (e.g. 90 deg minus 0 deg)
%       dataStack : [NxMxZ], stack of images to be co-registered
%       plotFlag  : 0 or 1; if 1, display pre- and post-registration overlays
%
%   OUTPUTS:
%       imagesCoreg : [NxMxZ], the coregistered data stack (shifted + masked)
%       tformParams : struct with fields:
%                       .shiftX, .shiftY : (in pixels, approximate)
%                       .affine2d        : the final MATLAB affine2d object
%
%   NOTES:
%       - Uses 'multimodal' metric -> Mattes Mutual Information by default.
%       - Only a translation transform is allowed (no rotation/scaling).
%       - We mask out NaNs in the refMap by setting them to zero, but 
%         'imregtform' itself doesn't exclude them automatically. 
%         You can also force them to be outside the FOV or use a 
%         masked registration approach if needed.
%       - The final shiftX, shiftY is estimated from tform.T(3,1:2), 
%         but keep in mind that imregtform can do sub-pixel alignment.

    if nargin < 5
        plotFlag = 0;
    end

    %% 1) Basic checks
    if ndims(refMap) ~= 2 || ndims(diffImage) ~= 2
        error('refMap and diffImage must be 2D.');
    end
    [nRows, nCols] = size(refMap);

    if any(size(diffImage) ~= [nRows, nCols])
        error('diffImage must match refMap size.');
    end
    if ndims(dataStack) ~= 3
        error('dataStack must be NxMxZ.');
    end
    if any(size(dataStack,1:2) ~= [nRows, nCols])
        error('dataStack must match refMap size in first two dims.');
    end
    if any(size(snrMask) ~= [nRows, nCols])
        error('snrMask must match refMap size.');
    end

    %% 2) Replace NaNs in refMap with 0 or some neutral value
    refMapNoNaN = refMap;
    nanMask     = isnan(refMapNoNaN);
    refMapNoNaN(nanMask) = 0;  

    %% 3) Setup imref2d for the fixed image
    Rfixed = imref2d([nRows, nCols]);

    %% 4) Prepare the 'optimizer' and 'metric' for 'multimodal' (Mattes MI)
    %     This returns default settings for a multimodal scenario.
    [optimizer, metric] = imregconfig('multimodal');
    % Sometimes you may want to tweak optimizer settings:
    optimizer.MaximumIterations = 1000; 
    optimizer.InitialRadius = 0.009; 
    optimizer.Epsilon = 1e-6;


    %% 5) imregtform: solve for translation
    %   fixed = refMapNoNaN, moving = diffImage
    diffImage(isnan(diffImage))=0;
    tformEstimate = imregtform(diffImage, Rfixed, refMapNoNaN, Rfixed, ...
                               'translation', optimizer, metric);

    % Extract the shift from the affine matrix 
    % tformEstimate.T = [ 1  0  0
    %                     0  1  0
    %                    Tx  Ty 1 ]
    shiftX = tformEstimate.T(3,1);
    shiftY = tformEstimate.T(3,2);

    fprintf('imregtform translation => shiftX=%.3f, shiftY=%.3f\n', ...
            shiftX, shiftY);

    %% 6) Optional: Show pre- and post-registration overlays
    if plotFlag
        %{
        figure('Name','Translation - imregtform','Color','w');
        
        % (A) Pre-registration
        subplot(1,2,1);
        preFused = imfuse(refMapNoNaN, diffImage, 'blend', 'Scaling','joint');
        imshow(preFused);
        title('Pre-registration Overlay (blend)');

        % (B) Post-registration
        subplot(1,2,2);
        % Warp diffImage with tformEstimate
        diffReg = imwarp(diffImage, tformEstimate, 'OutputView', Rfixed);
        postFused = imfuse(refMapNoNaN, diffReg, 'blend', 'Scaling','joint');
        imshow(postFused);
        title(sprintf(['Post-registration Overlay (blend)\n' ...
                       'shiftX=%.3f, shiftY=%.3f'], shiftX, shiftY));
        %}
    end

    %% 7) Warp the snrMask
    % Use nearest-neighbor so mask remains binary-like
    snrMaskReg = imwarp(snrMask, tformEstimate, 'OutputView', Rfixed, ...
                        'Interpolation','nearest');

    %% 8) Warp each slice of dataStack, apply the mask
    nSlices = size(dataStack,3);
    imagesCoreg = zeros(nRows, nCols, nSlices, 'like', dataStack);

    for z = 1:nSlices
        thisSlice = dataStack(:,:,z);
        thisSliceReg = imwarp(thisSlice, tformEstimate, ...
                              'OutputView', Rfixed, 'Interpolation','linear');
        % Apply the warped mask
        imagesCoreg(:,:,z) = thisSliceReg .* snrMaskReg;
    end

    %% 9) Collect transform results
    tformParams = struct('shiftX', shiftX, ...
                         'shiftY', shiftY, ...
                         'affine2d', tformEstimate);

end

function [rWeighted,correlationImage,correlationImageRef] = weightedCorrelation(refMap, data2D)
%WEIGHTEDCORRELATIONIGNORENAN2D Compute a weighted correlation between a 2D reference map
% (normalized to [-1, 1]) and a single 2D image of the same size, ignoring NaNs.
%
%   rVal = weightedCorrelationIgnoreNaN2D(refMap, data2D)
%
%   INPUTS:
%       refMap : [NxM] reference map, normalized to [-1, 1], may contain NaNs
%       data2D : [NxM] image, may contain NaNs
%
%   OUTPUT:
%       rVal   : Weighted correlation coefficient in the range [-1, 1].
%                NaN if insufficient valid data or zero variance.
%
%   Weights are W = abs(refMap). Pixels where refMap or data2D is NaN are excluded.

    % Basic input checks
    if ndims(refMap) ~= 2 || ndims(data2D) ~= 2
        error('Both refMap and data2D must be 2D.');
    end
    [nRows, nCols] = size(refMap);
    if any(size(data2D) ~= [nRows, nCols])
        error('data2D must have the same dimensions as refMap.');
    end

    % Flatten
    refMapVec = refMap(:);
    dataVec   = data2D(:);

    % Identify valid pixels (non-NaN in both)
    validIdx  = ~isnan(refMapVec) & ~isnan(dataVec);

    % If no valid pixels, return NaN
    if ~any(validIdx)
        rVal = NaN;
        return;
    end

    % Extract valid pixels
    refMapValid = refMapVec(validIdx);
    dataValid   = dataVec(validIdx);

    % Weights from absolute reference
    wVec = abs(refMapValid);
    wSum = sum(wVec);

    % If all weights are zero (e.g. reference is zero or NaN everywhere), return NaN
    if wSum == 0
        rVal = NaN;
        return;
    end

    % Weighted means
    wMeanRef   = sum(wVec .* refMapValid) / wSum;
    wMeanData  = sum(wVec .* dataValid  ) / wSum;

    % Weighted deviations
    refDev     = refMapValid - wMeanRef;
    dataDev    = dataValid   - wMeanData;

    % Numerator = sum of weighted product of deviations
    numerator  = sum(wVec .* refDev .* dataDev);

    % Denominator = sqrt of (weighted sum of squares of deviations)
    denomRef   = sum(wVec .* refDev.^2);
    denomData  = sum(wVec .* dataDev.^2);

    denominator = sqrt(denomRef * denomData);

    % If denominator is zero (e.g. uniform data or reference), return NaN
    if denominator == 0
        rWeighted = NaN;
    else
        rWeighted = numerator / denominator;  % Weighted correlation in [-1, +1]
    end

    % Unflatten the correlation values
    correlationImageRef = nan(size(refMap)); % Initialize with NaNs
    correlationValuesRef = wVec .* refDev;
    correlationImageRef(validIdx) = correlationValuesRef;

    % Unflatten the correlation values
    correlationImage = nan(size(refMap)); % Initialize with NaNs
    correlationValues = wVec .* refDev .* dataDev;
    correlationImage(validIdx) = correlationValues;
end

function [slicedImage, sliceAmplitude] = getImageSlice(image, sliceCenter, sliceWidth)
% Extracts a vertical slice from an image, masks the rest with NaNs, and
% returns the slice values.
%
% Args:
%   image: The input 2D image.
%   sliceCenter: The column index representing the center of the vertical slice.
%   sliceWidth: The total width of the slice in pixels.
%
% Returns:
%   slicedImage: An image of the same size as the input, where only the
%                pixels within the vertical slice are preserved, and the
%                rest are set to NaN.
%   sliceAmplitude: A 2D matrix containing the pixel values within the slice.
%                   The dimensions will be [number of rows in the image, sliceWidth].

    % Input checks
    if ~ismatrix(image)
        error('Input image must be a 2D matrix.');
    end
    if sliceCenter < 1 || sliceCenter > size(image, 2)
        error('sliceCenter must be within the image bounds.');
    end
    if sliceWidth < 1
        error('sliceWidth must be a positive integer.');
    end

    % Calculate slice boundaries
    halfWidth = floor(sliceWidth / 2);
    startCol = sliceCenter - halfWidth;
    endCol = sliceCenter + halfWidth;

    % Ensure slice boundaries are within image bounds
    startCol = max(1, startCol);
    endCol = min(size(image, 2), endCol);

    % Extract the slice values
    sliceAmplitude = nanmean(image(:, startCol:endCol),2);

    % Create a mask
    mask = false(size(image));
    mask(:, startCol:endCol) = true;

    % Apply mask and set out-of-slice pixels to NaN
    slicedImage = image;
    slicedImage(~mask) = NaN;
end

function [procMapGamma, procMapThresh, activityMapProc, Rfixed] = processActivityMap(activityMap, tileGridSize, gammaVal, sensitivity)
%PROCESSACTIVITYMAP Processes an activity map using CLAHE, gamma correction, and adaptive thresholding.
%
%   [procMapGamma, procMapThresh, activityMapProc, Rfixed] = processActivityMap(activityMap, tileGridSize, gammaVal, sensitivity)
%
%   Inputs:
%       activityMap  - The input activity map (2D matrix, may contain NaNs).
%       tileGridSize - [rows, cols] for the CLAHE algorithm (e.g., [16, 16]).
%       gammaVal     - Gamma correction factor (e.g., 0.5).
%       sensitivity  - Sensitivity for adaptive thresholding (0 to 1).
%
%   Outputs:
%       procMapGamma   - The activity map after CLAHE and gamma correction.
%       procMapThresh  - The thresholded activity map (binary, with NaNs).
%       activityMapProc- The processed map for registration (NaNs filled with 0).
%       Rfixed         - Spatial referencing object for the processed map.

    %% 2) HIST EQ & GAMMA CORRECTION ON ACTIVITY MAP (IGNORING NANS)
    validMask = ~isnan(activityMap);
    fillVal   = mean(activityMap(validMask), 'omitnan');  % neutral value
    % 2a) Fill NaNs, run CLAHE
    activityFilled = activityMap;
    activityFilled(~validMask) = fillVal;
    procMap = adapthisteq(activityFilled, ...
                          'NumTiles', tileGridSize, ...
                          'Range',   'full');  % or 'original'
    procMap(~validMask) = NaN;  % restore NaNs

    % 2b) Gamma correction
    temp2 = procMap;
    temp2(~validMask) = fillVal;
    procMapGamma = imadjust(temp2, [], [], gammaVal);
    procMapGamma(~validMask) = NaN;  % restore NaNs again

    %% 3) ADAPTTHRESH ON THE PROCESSED MAP
    %    This is for visualization/analysis, showing thresholded columns, etc.
    temp3 = procMapGamma;
    temp3(~validMask) = fillVal;

    nbhdSize = 2*floor(size(temp3)/16) + 1; % e.g. [9,9] if image is ~128
    threshVal = adaptthresh(double(temp3), sensitivity, ...
                            'NeighborhoodSize', nbhdSize);
    procMapThresh = double(imbinarize(double(temp3), threshVal));
    procMapThresh(procMapThresh==0) = NaN;  % keep invalid regions as NaN

    %% 4) PREPARE THE PROCESSED MAP FOR REGISTRATION
    %    We fill its NaNs with 0 so imregtform doesn't crash.
    activityMapProc = procMapGamma;
    nanMask = isnan(activityMapProc);
    activityMapProc(nanMask) = 0;
    Rfixed = imref2d(size(activityMapProc));
    %activityMapProc=rescale(activityMapProc);

end