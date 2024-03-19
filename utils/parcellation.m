function [fullROI, optostimROI, recruitROI]=parcellation(conditionIDs,dsCurrentSess,imageSeries,snrMask,gaussMask,bitmapMask)
%% Splits images into directly stimulated or non-stimulated parcels, per condition (visual + opto combination)
% Inputs:
%         conditionIDs = struct of condition IDs sorted into visual + opto combination
%         dsCurrentSess = datastruct containing info about experiment, here we extract info on the contour map applied to the gaussian
%         imageSeries = any raw/bandpassed GCaMP images with 512 * 512 * nCondition
%         snrMask = mask defined by the signal-noise ratio of reference map (flashed 12 orientations), defined by d'>~6
%         gaussMask = mask defined by gaussian-evoked activity (flashed gaussian), for identifying the center of activity

optostimMaskType='bitmap';
% Initialize
optostimROI=[];recruitROI=[];fullROI=[];
% Define n-optostim conds and n-blanks
blank=conditionIDs.blankConds;
nblank=numel(blank);
% Resize if necessary
if size(snrMask,1) > size(imageSeries,1) % if binned image
    for i=1:size(imageSeries,3)
        imageSeries(1:size(snrMask,1),1:size(snrMask,1),i)=imresize(imageSeries(:,:,i),size(snrMask,[1 2]) ,'bilinear');
    end
end
%% If gaussian mask used, use it else no
gaussStructExist=~isempty(gaussMask);
distancePrct=20/200;%15/200;
switch gaussStructExist
    case {1}
        % Extract gaussian levels used in experiment
        orts=dsCurrentSess.gaussianContourLevel; % grab the gaussian levels used
        maxGaussLevel=dsCurrentSess.gaussianContourLevelMax;
        bufferDist=round(maxGaussLevel * distancePrct);

        for ortNo=1:numel(orts)
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
            gaussLevel=max(orts) - bufferDist; % take the larger one of both, then add a buffer distance
            % select the type of mask used to define the stimulated ROI (gaussian) and nonstimulated ROI (inverse)
            switch optostimMaskType
                case {'gaussianfit'}
                    optostimAreaMask=double(gaussMask(gaussLevel).area);
                case {'bitmap'}
                     optostimAreaMask=bitmapMask(:,:,ortNo);
            end
            recruitmentMask=double(~optostimAreaMask);
            optostimAreaMask(optostimAreaMask==0)=NaN;
            recruitmentMask(recruitmentMask==0)=NaN;
            %% Mask the baseline images, split into stimulated and non-stimulated areas
            for condsWanted=1:numel(condsVisual)
                condID=condsVisual(condsWanted);
                % unmasked for sanity check
                fullROI(:,:,condID)=(imageSeries(:,:,condID)) .* snrMask;
                % mask it twice, once for d' mask, then split into stimulated or nonstimulated area
                optostimROI(:,:,condID)=(imageSeries(:,:,condID)) .* optostimAreaMask .* snrMask;
                recruitROI(:,:,condID)=(imageSeries(:,:,condID)) .* recruitmentMask .* snrMask;
            end

            %% Mask the optostim images, split into stimulated and non-stimulated areas
            for condsWanted=1:numel(condsOptostim)
                condID=condsOptostim(condsWanted);
                % unmasked for sanity check
                fullROI(:,:,condID)=(imageSeries(:,:,condID)) .* snrMask;
                % mask it twice, once for d' mask, then split into stimulated or nonstimulated area
                optostimROI(:,:,condID)=(imageSeries(:,:,condID)) .* optostimAreaMask .* snrMask;
                recruitROI(:,:,condID)=(imageSeries(:,:,condID)) .* recruitmentMask .* snrMask;
            end
        end
        
    case {0}       % No gaussian masking
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

            %% Mask the baseline images, split into stimulated and non-stimulated areas
            for condsWanted=1:numel(condsVisual)
                condID=condsVisual(condsWanted);
                % unmasked for sanity check
                fullROI(:,:,condID)=(imageSeries(:,:,condID)) .* snrMask;
            end

            %% Mask the optostim images, split into stimulated and non-stimulated areas
            for condsWanted=1:numel(condsOptostim)
                condID=condsOptostim(condsWanted);
                % unmasked for sanity check
                fullROI(:,:,condID)=(imageSeries(:,:,condID)) .* snrMask;
            end
        end
end
end
    