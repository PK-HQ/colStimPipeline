function [fullROI, optostimROI, recruitROI]=parcellation(conditionIDs,dsCurrentSess,bitmapData,...
    imagingData,imageSeries,trialOutcomeType,gaussMask,bitmapMask,nColumns,...
    optostimMaskType,blockID)
%% Splits images into directly stimulated or non-stimulated parcels, per condition (visual + opto combination)
% Inputs:
%         conditionIDs = struct of condition IDs sorted into visual + opto combination
%         dsCurrentSess = datastruct containing info about experiment, here we extract info on the contour map applied to the gaussian
%         imageSeries = any raw/bandpassed GCaMP images with 512 * 512 * nCondition
%         imagingData.snrMaskCoreg(:,:,blockID) = mask defined by the signal-noise ratio of reference map (flashed 12 orientations), defined by d'>~6
%         gaussMask = mask defined by gaussian-evoked activity (flashed gaussian), for identifying the center of activity
imageSeries=imageSeries.(trialOutcomeType);
% Initialize
optostimROI=[];recruitROI=[];fullROI=[];
% Define n-optostim conds and n-blanks
blank=conditionIDs.blankConds;
nblank=numel(blank);
optostimROI=[];
recruitROI=[];
% Resize if necessary
if size(imagingData.snrMaskCoreg(:,:,blockID),1) > size(imageSeries,1) % if binned image
    for i=1:size(imageSeries,3)
        imageSeries(1:size(imagingData.snrMaskCoreg(:,:,blockID),1),1:size(imagingData.snrMaskCoreg(:,:,blockID),1),i)=imresize(imageSeries(:,:,i),size(imagingData.snrMaskCoreg(:,:,blockID),[1 2]) ,'bilinear');
    end
end
%% If gaussian mask used, use it else no
gaussStructExist=0;%~isempty(gaussMask);
distancePrct=20/200;%15/200;
switch gaussStructExist
    %{
    case {1}
        disp('Gaussian parcellation')
        % Extract gaussian levels used in experiment
        orts=size(bitmapData.gaussianContourLevel(:,:,blockID),2); % grab the gaussian levels used


        for ortNo=1:orts
            gaussLevel=bitmapData.gaussianContourLevel(:,ortNo,blockID);
            maxGaussLevel=bitmapData.gaussianContourLevelMax(ortNo,blockID);
            bufferDist=0;%round(maxGaussLevel * distancePrct);
            %% Define the optostim+visual and visual-only conditions for the optostim ROI mask (0 or 90), based on which optostim was applied during experiment
            if ortNo==1 % if opto 0 was applied
                % Q: does opto 0 drive recruitment similar to visual 0? 
                condsVisual=conditionIDs.V0;
                condsOptostim=[conditionIDs.V0O0 conditionIDs.V90O0];
                if isempty(condsOptostim) & ~isempty(conditionIDs.opto0Conds) % for fixtation optostim blocks
                    condsOptostim=[conditionIDs.opto0Conds];
                end
            else % if opto 90 was applied
                % Q: does opto 90 drive recruitment similar to visual 90? 
                condsVisual=conditionIDs.V90;
                condsOptostim=[conditionIDs.V0O90 conditionIDs.V90O90];
                if isempty(condsOptostim) & ~isempty(conditionIDs.opto0Conds) % for fixtation optostim blocks
                    condsOptostim=[conditionIDs.opto90Conds];
                end
            end
            %% Create binary mask of each image, with the mask = gaussian or inverse of the gaussian
            % grab the gaussian level previously used to create optostim ROI (out of gaussian max level)
            gaussLevel=gaussLevel - bufferDist; % take the larger one of both, then add a buffer distance
            % select the type of mask used to define the stimulated ROI (gaussian) and nonstimulated ROI (inverse)
            switch optostimMaskType
                case {'gaussian'}
                    optostimAreaMask=double(gaussMask(gaussLevel).area);
                case {'bitmap'}
                     optostimAreaMask=bitmapMask(:,:,ortNo,blockID);
            end
            recruitmentMask=double(~optostimAreaMask);
            optostimAreaMask(optostimAreaMask==0)=NaN;
            recruitmentMask(recruitmentMask==0)=NaN;

            %% Mask the baseline images, split into stimulated and non-stimulated areas
            for condsWanted=1:numel(condsVisual)
                condID=condsVisual(condsWanted);
                % unmasked for sanity check
                fullROI(:,:,condID)=(imageSeries(:,:,condID));
                % mask it twice, once for d' mask, then split into stimulated or nonstimulated area
                optostimROI(:,:,condID)=(imageSeries(:,:,condID)) .* optostimAreaMask;
                recruitROI(:,:,condID)=(imageSeries(:,:,condID)) .* recruitmentMask;
            end

            %% Mask the optostim images, split into stimulated and non-stimulated areas
            for condsWanted=1:numel(condsOptostim)
                condID=condsOptostim(condsWanted);
                % unmasked for sanity check
                fullROI(:,:,condID)=(imageSeries(:,:,condID));
                % mask it twice, once for d' mask, then split into stimulated or nonstimulated area
                optostimROI(:,:,condID)=(imageSeries(:,:,condID)) .* optostimAreaMask;
                recruitROI(:,:,condID)=(imageSeries(:,:,condID)) .* recruitmentMask;
            end
        end
        %}
    case {0}       % Bitmap masking
        disp('Bitmap parcellation')       
        %{
        % Define the affine transformation matrix
        T = [1 0 0; 0 1 0; imagingData.Stab(blockID,1) imagingData.Stab(blockID,2) 1]; % Translation matrix in affine form
        
        % Create the affine2d object
        tform = affine2d(T);
        sizeImage=size(bitmapMask,[1 2]);
        % Apply the translation to the image
        for i=1:size(bitmapMask,3)
            bitmapMaskCoreg(:,:,i,blockID) = imwarp(bitmapMask(:,:,i,blockID), tform, 'OutputView', imref2d(sizeImage));
        end
        %}
        for visualOrt=1:2
            %% Define the optostim+visual and visual-only conditions for the optostim ROI mask (0 or 90), based on which optostim was applied during experiment
            if visualOrt==1 % if opto 0 was applied
                % Q: does opto 0 drive recruitment similar to visual 0? 
                condsVisual=conditionIDs.V0;
                condsOptostim=[conditionIDs.V0O0 conditionIDs.V90O0];
            else % if opto 90 was applied
                % Q: does opto 90 drive recruitment similar to visual 90? 
                condsVisual=conditionIDs.V90;
                condsOptostim=[conditionIDs.V0O90 conditionIDs.V90O90];
            end
            
            activityMap=(imageSeries(:,:,conditionIDs.V90O90(1))-imageSeries(:,:,conditionIDs.V0O0(1))).*imagingData.snrMaskCoreg(:,:,blockID);

           %% Coreg bitmap to account for some slight offset/movement
            bitmapDiff=bitmapMask(:,:,2,blockID)-bitmapMask(:,:,1,blockID);
            if nColumns<8
                tileGridSize = [64,64];  % for CLAHE
                gammaVal = 1;     % gamma correction factor
                sensitivity= 1;     % threshold sensitivity
            else
                tileGridSize = [16,16];  % for CLAHE
                gammaVal = .5;     % gamma correction factor
                sensitivity= 1;     % threshold sensitivity
            end
            plotFlag=1;

            %{
            [bitmapMaskCoreg(:,:,2,blockID), tformParams, procMapThresh] = coregisterImages(...
               activityMap, bitmapDiff, bitmapMask(:,:,2,blockID), tileGridSize, gammaVal, sensitivity, plotFlag);
            [bitmapMaskCoreg(:,:,1,blockID), tformParams, procMapThresh] = coregisterImages(...
               activityMap, bitmapDiff, bitmapMask(:,:,1,blockID), tileGridSize, gammaVal, sensitivity, plotFlag);

            % Compute union of both bitmap masks
            bitmapMaskCoregUnion = double( ...
                (bitmapMaskCoreg(:,:,2,blockID) > 0) | (bitmapMaskCoreg(:,:,1,blockID) > 0));
            %}

            
            %% Define mask using bitmap
            bufferPercent=1;
            [optostimAreaMask, ~] = fitBestEllipse(bitmapData.bitmapMask(:,:,blockID), optostimMaskType, bufferPercent);
            recruitmentMask=double(~optostimAreaMask);
            optostimAreaMask=double(optostimAreaMask);
            optostimAreaMask(optostimAreaMask==0)=NaN;
            optostimAreaMask=optostimAreaMask.*imagingData.snrMaskCoreg(:,:,blockID);
            recruitmentMask(recruitmentMask==0)=NaN;
            recruitmentMask=recruitmentMask.*imagingData.snrMaskCoreg(:,:,blockID);
           
            %% Mask the baseline images, split into stimulated and non-stimulated areas
            for condsWanted=1:numel(condsVisual)
                condID=condsVisual(condsWanted);
                % unmasked for sanity check
                fullROI(:,:,condID)=(imageSeries(:,:,condID));
                % mask it twice, once for d' mask, then split into stimulated or nonstimulated area
                optostimROI(:,:,condID)=(imageSeries(:,:,condID)) .* optostimAreaMask;
                recruitROI(:,:,condID)=(imageSeries(:,:,condID)) .* recruitmentMask;
            end

            %% Mask the optostim images, split into stimulated and non-stimulated areas
            for condsWanted=1:numel(condsOptostim)
                condID=condsOptostim(condsWanted);
                % unmasked for sanity check
                fullROI(:,:,condID)=(imageSeries(:,:,condID));
                % mask it twice, once for d' mask, then split into stimulated or nonstimulated area
                optostimROI(:,:,condID)=(imageSeries(:,:,condID)) .* optostimAreaMask;
                recruitROI(:,:,condID)=(imageSeries(:,:,condID)) .* recruitmentMask;
            end
        end
                    %{
        figure('Name','Recruitment mask');
        subplot(2,4,1)
        imagesc(fullROI(:,:,conditionIDs.V0O0(end))); axis square; cax()
        subplot(2,4,2)
        imagesc(recruitROI(:,:,conditionIDs.V0O0(end))); axis square; cax()
        subplot(2,4,3)
        imagesc(bitmapMask(:,:,1,blockID)); axis square
        subplot(2,4,4)
        imagesc(bitmapMaskCoreg(:,:,1,blockID)); axis square

        subplot(2,4,5)
        imagesc(fullROI(:,:,conditionIDs.V90O90(end))); axis square; cax()
        subplot(2,4,6)
        imagesc(recruitROI(:,:,conditionIDs.V90O90(end))); axis square; cax()
        subplot(2,4,7)
        imagesc(bitmapMask(:,:,2,blockID)); axis square
        subplot(2,4,8)
        imagesc(bitmapMaskCoreg(:,:,2,blockID)); axis square
                    %}
        
end
end

function [outputMask, ellipseParams] = fitBestEllipse(optostimAreaMask, mode, bufferPercent)
    % FITBESTELLIPSE Fits a best-fit ellipse around non-zero blobs in a 2D mask
    % and generates an output mask based on the specified mode: 'gaussian' or 'bitmap'.
    %
    % Inputs:
    %   optostimAreaMask - 2D binary or non-binary image with blobs (non-zero regions)
    %   mode - 'gaussian' or 'bitmap'
    %   bufferPercent - Percentage to expand the area (e.g., 10 for 10%)
    %
    % Outputs:
    %   outputMask - Binary mask with the ellipse or expanded blobs
    %   ellipseParams - Struct containing ellipse parameters:
    %                   - Center: [xc, yc]
    %                   - Axes: [a, b] (semi-major and semi-minor axes)
    %                   - Angle: Rotation angle in degrees

    % Validate inputs
    if nargin < 3
        bufferPercent = 0; % Default buffer is 0%
    end
    if ~ismember(mode, {'gaussian', 'bitmap'})
        error('Mode must be ''gaussian'' or ''bitmap''.');
    end

    % Find non-zero pixel coordinates
    [rowIdx, colIdx] = find(optostimAreaMask);
    
    % Perform Principal Component Analysis (PCA) on the pixel coordinates
    coords = [colIdx, rowIdx]; % x = colIdx, y = rowIdx
    [coeff, ~, latent] = pca(coords); % coeff contains eigenvectors, latent eigenvalues
    
    % Calculate the center of the ellipse
    center = mean(coords, 1);
    
    % Semi-major and semi-minor axes (sqrt of eigenvalues scaled by std)
    a = sqrt(latent(1)) * 2; % Semi-major axis length
    b = sqrt(latent(2)) * 2; % Semi-minor axis length
    
    % Angle of rotation (in degrees)
    angle = atan2d(coeff(2, 1), coeff(1, 1)); % Angle of the first principal component
    
    % Ensure angle aligns bottom-left to top-right
    if angle < 0
        angle = angle + 180; % Adjust angle to ensure positive rotation
    end

    % Apply bufferPercent to the axes
    bufferFactor = 1 + bufferPercent / 100;
    a = a * bufferFactor;
    b = b * bufferFactor;

    % Package ellipse parameters
    ellipseParams = struct('Center', center, ...
                           'Axes', [a, b], ...
                           'Angle', angle);
    
    % Generate the output mask based on the mode
    switch mode
        case 'gaussian'
            % Create an elliptical mask with the buffered ellipse
            outputMask = generateEllipseMask(size(optostimAreaMask), center, a, b, angle);

        case 'bitmap'
            % Expand each blob individually
            outputMask = expandBlobs(optostimAreaMask, bufferPercent);
    end
end

function ellipseMask = generateEllipseMask(imageSize, center, a, b, angle)
    % GENERATEELLIPSEMASK Generates a binary mask with a fitted ellipse
    %
    % Input:
    %   imageSize - [rows, cols] of the output mask
    %   center - [xc, yc], center of the ellipse
    %   a, b - semi-major and semi-minor axes
    %   angle - rotation angle in degrees
    %
    % Output:
    %   ellipseMask - Binary mask with the ellipse (1 inside, 0 outside)
    
    [cols, rows] = meshgrid(1:imageSize(2), 1:imageSize(1)); % Coordinate grids
    
    % Shift coordinates to center the ellipse
    xc = center(1);
    yc = center(2);
    xShifted = cols - xc;
    yShifted = rows - yc;
    
    % Rotation matrix
    cosAngle = cosd(angle);
    sinAngle = sind(angle);
    xRotated = cosAngle * xShifted + sinAngle * yShifted;
    yRotated = -sinAngle * xShifted + cosAngle * yShifted;
    
    % Equation of the ellipse
    ellipseEq = (xRotated.^2) / (a^2) + (yRotated.^2) / (b^2);
    
    % Generate mask (1 inside ellipse, 0 outside)
    ellipseMask = ellipseEq <= 1;
end

function expandedMask = expandBlobs(inputMask, bufferPercent)
    % EXPANDBLOBS Expands each blob in the input mask by bufferPercent
    %
    % Input:
    %   inputMask - Binary mask with blobs
    %   bufferPercent - Percentage to expand each blob
    %
    % Output:
    %   expandedMask - Binary mask with expanded blobs

    % Compute the structural element for dilation
    seSize = ceil(bufferPercent / 100 * min(size(inputMask)));
    se = strel('disk', seSize);

    % Apply morphological dilation to expand blobs
    expandedMask = imdilate(inputMask > 0, se);
end
function [bitmapCoreg, tformParams, procMapThresh] = coregisterImages(...
    activityMap, bitmapDiff, bitmapOrt, tileGridSize, gammaVal, sensitivity, plotFlag)
%COREGISTERIMAGES 
%   1) Applies adaptive histogram equalization (CLAHE) + gamma correction
%      to an activityMap that may contain NaNs.
%   2) Also uses adaptthresh on the *processed* activity map (for visualization).
%   3) Estimates a Rigid-register transform from (translate + rotate) 
%       binary reference bitmap (bitmap 90-0) to reference activity (activity 90-0, min contrast)
%       processed map using 'imregtform' (Mattes Mutual Info).
%   4) Applies transform to bitmap (0 or 90) for masking the activity (0 or 90) downstream
%   5) Optionally shows a 2×2 figure:
%         (A) Original activity map
%         (B) Processed map (CLAHE + gamma)
%         (C) Thresholded map via adaptthresh
%         (D) Post-registration overlay
%
%   USAGE:
%   [bitmapCoreg, tformParams, procMapThresh] = coregisterImages(...
%       activityMap, bitmap, tileGridSize, gammaVal, plotFlag)
%
%   INPUTS:
%       activityMap : [NxM], real 2D image (NaNs where invalid).
%       bitmap      : [NxM], 0/1 binary image (moving), must match activityMap size.
%       tileGridSize: 2-element vector for 'NumTiles' in adapthisteq, e.g. [8 8].
%       gammaVal    : scalar gamma for imadjust (e.g. 1.0 => no change).
%       plotFlag    : 0 or 1; if 1, plots the 2×2 figure described above.
%
%   OUTPUTS:
%       bitmapCoreg : [NxM], warped bitmap in registration with the processed map
%       tformParams : struct with fields:
%                       .tform    -> the final affine2d
%                       .shiftX   -> X translation (pixels)
%                       .shiftY   -> Y translation (pixels)
%                       .rotation -> rotation in degrees
%       procMapThresh: [NxM], the thresholded map (NaNs left in place)
%
%   NOTES:
%       - The code uses 'rigid' => translation + rotation only (no scaling).
%       - We do hist eq + gamma on valid (non-NaN) pixels by temporarily
%         filling NaNs with a neutral value.
%       - 'adaptthresh' is applied to the final processed map for
%         visualization. We restore NaNs in that thresholded map as well.
%       - The final registration is between the processed map (fixed)
%         and the input 'bitmap' (moving).
%
%       => So all hist eq, gamma corr, and adaptthresh are on the activityMap
%          prior to doing the coreg.

    %% 1) Basic checks & default plotFlag
    if nargin < 5, plotFlag = 0; end
    if ndims(activityMap) ~= 2 || ndims(bitmapDiff) ~= 2
        error('Both activityMap and bitmap must be 2D arrays.');
    end
    [nRows, nCols] = size(activityMap);
    if any(size(bitmapDiff) ~= [nRows, nCols])
        error('bitmap must match activityMap in size.');
    end
    % Ensure bitmap is double so imregtform can handle it
    if islogical(bitmapDiff)
        bitmapDiff = double(bitmapDiff);
    end
    % Convert bitmap to 1s
    bitmapDiff=abs(bitmapDiff);
    activityMap=abs(activityMap);

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
    procMapThresh(~validMask) = NaN;  % keep invalid regions as NaN

    %% 4) PREPARE THE PROCESSED MAP FOR REGISTRATION
    %    We fill its NaNs with 0 so imregtform doesn't crash.
    activityMapProc = procMapGamma;
    nanMask = isnan(activityMapProc);
    activityMapProc(nanMask) = 0;
    Rfixed = imref2d(size(activityMapProc));
    activityMapProc=rescale(activityMapProc);
    %% 5) CONFIGURE REGISTRATION (MULTIMODAL => MATTES MUTUAL INFO), RIGID
    %optimizer
    [optimizer,metric] = imregconfig('multimodal');
    optimizer.InitialRadius = 6.25*10^-4;
    optimizer.GrowthFactor = 1.05;

    tformEstimate = imregtform(bitmapDiff, Rfixed, ...        % moving
                               activityMapProc, Rfixed, ... % fixed
                               'translation', optimizer, metric,...
                               'PyramidLevels',3);

    % Extract shift/rotation
    shiftX = tformEstimate.T(3,1);
    shiftY = tformEstimate.T(3,2);
    %{
    rotRad = atan2(tformEstimate.T(2,1), tformEstimate.T(1,1));
    rotationDeg = rotRad * (180/pi);

    fprintf('Rigid transform => shiftX=%.3f, shiftY=%.3f, rotation=%.2f deg\n',...
            shiftX, shiftY, rotationDeg);
    %}
    %% 6) WARP THE BITMAP WITH NEAREST-NEIGHBOR
    bitmapCoreg = imwarp(bitmapOrt, tformEstimate, ...
                         'OutputView', Rfixed, 'Interpolation','nearest');

    %% 7) PLOT (OPTIONAL) - 4 SUBPLOTS
    if plotFlag
        %{
        figure('Name','coregisterImages: HistEq + Gamma + AdaptThresh + Rigid Reg','Color','w');
        
        % (A) Top-left: original
        subplot(2,2,1);
        imagesc(activityMap); axis image; colormap gray;
        title('Original Activity Map');
        colorbar;

        % (B) Top-right: processed map
        subplot(2,2,2);
        imagesc(procMapGamma); axis image; colormap gray;
        title('Processed Map');
        colorbar;

        % (C) Bottom-left: thresholded map
        preFuse = imfuse(activityMapProc, bitmapDiff, ...
                                            'falsecolor', 'Scaling','independent', ...
                                            'ColorChannels','green-magenta');
        subplot(2,2,3);
        imagesc(preFuse); axis image; colormap(gray);
        title('Pre-registration');
        colorbar;

        % (D) Bottom-right: post-registration overlay
        postFuse = imfuse(activityMapProc, bitmapCoreg, ...
                                            'falsecolor', 'Scaling','independent', ...
                                            'ColorChannels','green-magenta');
        subplot(2,2,4);
        imshow(postFuse); offwarning()
        title(sprintf('Post-registration (Tx=%.2f, Ty=%.2f)',...
                      shiftX, shiftY));
        %}
    end

    %% 8) RETURN TRANSFORM & OUTPUT
    tformParams = struct('tform',    tformEstimate, ...
                         'shiftX',   shiftX, ...
                         'shiftY',   shiftY);

end
