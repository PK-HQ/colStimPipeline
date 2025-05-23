function neuroStruct=analyzeProjection1D(mainPath, datastruct, behavioralData, imagingData, bitmapData, analysisBlockID, trialOutcomeType, ...
    hAx, hAxPanel, optostimAreaFlag, clusterIdx, filterColumns,optostimMaskType, saveFlag)
% EXEMPLAR: #81: cluster=4 block=5

% Colors
colorCon = [0, 225, 80] / 255;colorIncon = [156, 14, 254] / 255;
errorIntgBlocks={};
nClusters=numel(unique(clusterIdx));
for cluster=3%
    clusterBlocksAll=find(clusterIdx==cluster);

    % filter for n columns or not
    switch filterColumns
        case 1
            columnsDesired = 20;
            columnSpread = 4;
            %energyThreshold = 3.5;
            clusterBlocks = selectNColumnBlocks(bitmapData, clusterBlocksAll, columnsDesired, columnSpread);
        case 0
            clusterBlocks=clusterBlocksAll;
    end
    neuroStruct.(['C' num2str(cluster)]).analyzeBlockIDs=clusterBlocks;
    %Skip to next cluster if none match
    if ~isempty(clusterBlocks)
        for blockNo=19%14%1:numel(clusterBlocks)%[3 5 10 11 14]
            tic
            blockID=clusterBlocks(blockNo);
            imagingData.optoIntg=[];
            imagingData.baselineIntg=[];
    
            % Load block data individually as imaging data too large
             [currentBlockStruct,~,...
                behavioralData, imagingData, ~, ~]=loadBlockData(datastruct, analysisBlockID, behavioralData, imagingData, bitmapData, blockID, '');

                filenameStructCurrent=generateFilenames(currentBlockStruct);

            monkey=datastruct(analysisBlockID(blockID)).monkey;
            date=datastruct(analysisBlockID(blockID)).date;
            run=datastruct(analysisBlockID(blockID)).run;
            pdfFilename=filenameStructCurrent.recruitmentPDF;
            blockColumns=bitmapData.nColumns(:,blockID);
           %% STAGE 2. For each intg image, split image into stimulated and nonstimulated area with gaussian ROI imagingData.mask
            % Determine if opto+baseline are combined (nTS==1) or in different blocks (nTS==2)
            nTS = 1;
            if isfield(behavioralData, 'baselineTS') && ~isempty(behavioralData.baselineTS(blockID).Header)
                nTS = 2;
            end
            % Get behavioral and neural responses from valid completed trials
            [images, TSbaseline, TSopto, condIDs] = processImagingData(nTS, behavioralData, imagingData, blockID);
            
           %% Define SNR mask, reference map, gaussian footprint mask, spatially filter for columns
            % Define snrMask
            ROImask=double(imagingData.mask(:,:,blockID)); snrMask=ROImask; snrMask(snrMask==0)=NaN;
            imagingData.snrMask(:,:,blockID)=snrMask;
            % Define refMap(pcaMap)
             imagingData.refmap(:,:,1:2,blockID)=imagingData.ortpca(:,:,[1 7],blockID);
            % Define gaussianMask
            gaussMask=imagingData.gaussfit(:,:,blockID); %rename
            % Option for bandpass spatial filtering for columns
            images=columnarFilter(behavioralData.optoTS(blockID),images,trialOutcomeType);

           %% Coregister all images to block responses
            % Coregister reference map (pca90-0) to activity
             plotFlag=0;
             [imagingData] = coregisterRefMap2Responses(...
                images, imagingData, blockColumns, condIDs, blockID, trialOutcomeType, plotFlag);
            % Coregister bitmap (pca90-0) to activity
            plotFlag=0;
             [bitmapData] = coregisterBitmap2Responses(...
                images, imagingData, bitmapData, blockColumns, condIDs, blockID, trialOutcomeType, plotFlag);
            
           %% STAGE 3. Split the image into optostim and recruitment ROI
            columnarResp = splitImageIntoROIs(blockColumns, condIDs, currentBlockStruct, bitmapData, imagingData, ...
                images, trialOutcomeType, gaussMask, optostimMaskType, blockID, optostimAreaFlag);
            
           %% STAGE 4. Plot processed images for ROI
            %plotProcessedImages(imagingData, columnarResp, condIDs, blockID);
            plotTypeFlag=2;
            plotNeurometricData(bitmapData, imagingData, cluster, blockNo, columnarResp,  blockID, condIDs,...
                TSbaseline, TSopto,...
                blockColumns, neuroStruct, colorCon, colorIncon, optostimAreaFlag,...
                plotTypeFlag); 
            %{
            plotTypeFlag=3;
            plotNeurometricDataSlice(imagingData, bitmapData, cluster, blockNo, blockID, columnarResp, condIDs, TSbaseline, TSopto,...
                blockColumns, neuroStruct, colorCon, colorIncon, optostimAreaFlag, plotTypeFlag); 
            %}
           %% Saving
            fileSuffix='proj';
            saveFigures(saveFlag, mainPath, monkey, cluster, date, run, optostimAreaFlag, fileSuffix)
            toc
        end   
    else
        disp('Cluster empty')
    end
end
end



%% Subfunctions

% --- Subfunction to process imaging data ---
function [images, TSbaseline, TSopto, condIDs] = processImagingData(nTS, behavioralData, imagingData, blockID)

    switch nTS
        case 1
            % Single time series case
            if isempty(imagingData.optoIntg)
                images = []; TSbaseline = []; TSopto = []; condIDs = [];
                return; % Early exit if opto data is empty
            else
                [~, images, condIDs] = getUsableTrials(behavioralData.optoTS(blockID), imagingData.optoIntg(:,:,:, blockID));
            end

            TSbaseline = behavioralData.optoTS(blockID);
            TSopto = behavioralData.optoTS(blockID);

        case 2
            % Two time series case (baseline and opto)
            if isempty(imagingData.baselineIntg) || isempty(imagingData.optoIntg)
                images = []; TSbaseline = []; TSopto = []; condIDs = [];
                return; % Early exit if baseline or opto data is empty
            else
                [~, imagesbl, condIDsbl] = getUsableTrials(behavioralData.baselineTS(blockID), imagingData.baselineIntg(:,:,:, blockID));
                [~, imagesOpto, condIDsOpto] = getUsableTrials(behavioralData.optoTS(blockID), imagingData.optoIntg(:,:,:, blockID));
            end

            [imagesOpto, condIDsOpto] = appendBaselineData(imagesbl, condIDsbl, imagesOpto, condIDsOpto);

            % Rename variables for consistency
            images = imagesOpto;
            TSbaseline = behavioralData.baselineTS(blockID);
            TSopto = behavioralData.optoTS(blockID);
            condIDs = condIDsOpto;
    end
end

% --- Subfunction to append baseline data to opto data ---
function [imagesOpto, condIDsOpto] = appendBaselineData(imagesbl, condIDsbl, imagesOpto, condIDsOpto)

    condsOpto = processCondIDs(condIDsOpto);
    condsBL = processCondIDs(condIDsbl);

    % Append baseline conditions and images
    maxOptoConds = max(condsOpto);
    BLCondIdx = maxOptoConds + (1:max(condsBL));

    condIDsOpto.V0 = [condIDsOpto.V0, condIDsbl.V0 + maxOptoConds];
    condIDsOpto.V90 = [condIDsOpto.V90, condIDsbl.V90 + maxOptoConds];
    imagesOpto.average(:,:,BLCondIdx) = imagesbl.average;
    imagesOpto.averageCorrect(:,:,BLCondIdx) = imagesbl.averageCorrect;
    imagesOpto.averageError(:,:,BLCondIdx) = imagesbl.averageError;

end

% --- Helper function to process condition IDs ---
function conds = processCondIDs(condIDs)
    % Processes condition IDs from getUsableTrials output.
    fields = fieldnames(condIDs);
    conds = NaN(1, numel(fields));
    for fn = 1:numel(fields)
        name = fields{fn};
        if ~isempty(condIDs.(name))
            conds(fn) = max(condIDs.(name));
        end
    end
end

function formattedTitle = createColoredTitle(mapName, meanProjection, phi, colorCon, colorIncon)
    % createColoredTitle generates a formatted LaTeX title with colored text for mean values.
    %
    % Inputs:
    %   mapName       - String representing the name of the map.
    %   meanProjection - Vector containing projection mean values, with congruent at index 2 and incongruent at index 3.
    %   phi           - Scalar value of phi to display in the title.
    %   colorCon      - 1x3 RGB vector for the congruent color.
    %   colorIncon    - 1x3 RGB vector for the incongruent color.
    %
    % Outputs:
    %   formattedTitle - String formatted for use as a LaTeX title.

    % Format meanProjection values
    meanCon = sprintf('%.1f', meanProjection(2));
    meanIncon = sprintf('%.1f', meanProjection(3));
    phiStr = sprintf('%.2f', phi);

    % Convert RGB colors to LaTeX color strings
    colorConStr = sprintf('{\\color[rgb]{%.3f,%.3f,%.3f}', colorCon(1), colorCon(2), colorCon(3));
    colorInconStr = sprintf('{\\color[rgb]{%.3f,%.3f,%.3f}', colorIncon(1), colorIncon(2), colorIncon(3));

    % Create LaTeX-formatted title
    formattedTitle = sprintf('%s\nSimilarity = %s%s}, %s%s}', ...
                              mapName, colorConStr, meanCon, colorInconStr, meanIncon);
    %formattedTitle = sprintf(['%s\nMean proj. = %s%s}, %s%s}, \\Phi = %s'], ...
    %                          mapName, colorConStr, meanCon, colorInconStr, meanIncon, phiStr);
end

function yAtZero = handleZeroData(xData, yData)
    % Initialize the output
    yAtZero = [];
    
    % Find indices where xData is exactly 0
    zeroIdx = xData == 0;
    
    if any(zeroIdx)
        % Case 1: xData contains 0
        yAtZero = mean(yData(zeroIdx));
    else
        % Case 2: No exact 0 in xData, find flanking points
        lowerIdx = find(xData < 0, 1, 'last'); % Last value less than 0
        upperIdx = find(xData > 0, 1, 'first'); % First value greater than 0
        
        if ~isempty(lowerIdx) && ~isempty(upperIdx)
            % Perform linear interpolation to estimate yData at xData == 0
            xLower = xData(lowerIdx);
            xUpper = xData(upperIdx);
            yLower = yData(lowerIdx);
            yUpper = yData(upperIdx);
            
            % Interpolating yData at xData == 0
            yAtZero = yLower + (yUpper - yLower) * (0 - xLower) / (xUpper - xLower);
        elseif ~isempty(zeroIdx)
            % Handle edge cases where no flanking points exist
            disp('Unable to find flanking points around 0 in xData.');
        end
    end
end

%% Plotting subfuncs
% Helper function to generate grid indices for plotting
function plotIndices = generateGridIdx(nRows, nCols)
    plotIndices = reshape(1:nRows*nCols, nCols, nRows)';
end

% --- Function to split the image into ROIs (Stage 3) ---
function columnarResp = splitImageIntoROIs(blockColumns, condIDs, currentBlockStruct, bitmapData, imagingData, ...
    images, trialOutcomeType, gaussMask, optostimMaskType, blockID, optostimAreaFlag)

    meanColumns = mean(blockColumns);
    [fullROI, optostimROI, recruitROI] = parcellation(condIDs, currentBlockStruct, bitmapData, imagingData, ...
        images, trialOutcomeType, gaussMask, bitmapData.columnarbitmapTFcamspace, meanColumns, ...
        optostimMaskType, blockID);

    % Choose the area analyzed
    switch optostimAreaFlag
        case {'full'}
            columnarResp = fullROI;
            titleStr = 'Neurometric (full ROI)';
        case {'stim'}
            columnarResp = optostimROI;
            titleStr = 'Neurometric (Stimulated ROI)';
        case {'nonstim'}
            columnarResp = recruitROI;
            titleStr = 'Neurometric (Recruitment ROI)';
    end
end

% --- Function to plot processed images for ROI (Stage 4) ---
function plotProcessedImages(imagingData, columnarResp, condIDs, blockID);

    % SNR mask
    snrMaskCoreg=imagingData.snrMaskCoreg(:,:,blockID);
    snrMaskCoreg(snrMaskCoreg==0)=NaN;
    figure
    nCols = 6;
    nRows = size(columnarResp, 3) / nCols;
    plotIndices = generateGridIdx(nRows, nCols);
    colormap(fireice);
    gap = .0075;
    marginV = .01;
    marginH = .01;
    [hAx, ~] = tight_subplot(nRows, nCols, [gap + 0.025 gap], [marginV marginV + .06], [marginH + .06 marginH]);

    for condSet = 1:3
        % Define conditions to plot
        switch condSet
            case {1} % BL V-only
                conds = [condIDs.V0; condIDs.V90];
                yData = columnarResp(:,:,conds);
                plotIdx = plotIndices(:,[1 4])';
                superYData = yData.*snrMaskCoreg;
            case {2} % V0
                conds = [condIDs.V0O0; condIDs.V0O90];
                yData = columnarResp(:,:,conds);
                plotIdx = plotIndices(:,[2 3])';
                superYData = columnarResp(:,:,[condIDs.V0O0; condIDs.V0O90; condIDs.V90O90; condIDs.V90O0]).*snrMaskCoreg;
            case {3} % V90
                conds = [condIDs.V90O0; condIDs.V90O90];
                yData = columnarResp(:,:,conds);
                plotIdx = plotIndices(:,[5 6])';
                superYData = columnarResp(:,:,[condIDs.V0O0; condIDs.V0O90; condIDs.V90O90; condIDs.V90O0]).*snrMaskCoreg;
        end

        for cond = 1:numel(plotIdx)
            axes(hAx(plotIdx(cond)))
            imgsc2(yData(:,:,cond).*snrMaskCoreg); caxis([-10 10]*10^-3);
            %cax2(.7, superYData.*snrMaskCoreg)
            upFontSize(10, .005)
            axis square
            colorbar
            addPix2MM(yData(:,:,cond), nCols, nRows, plotIdx(cond))
        end
        upFontSize(12, .005)

        % Store min-max
        minValue=[];maxValue=[];
        numImages = size(superYData, 3);  % Get the number of images (N)

        for i = 1:numImages
            currentImage = superYData(:,:,i);  % Extract the i-th image
            minValue(i) = min(currentImage(:));   % Find the minimum value
            maxValue(i) = max(currentImage(:));   % Find the maximum value
            %fprintf('Image %d: Min = %f, Max = %f\n', i, minValue, maxValue);
        end
        imgdata.(['Cond' num2str(condSet)])=[minValue;maxValue];

    end

    % Headers for visual contrast and conds
    h = suplabel('Visual contrast', 'y', [.1 .1 .8 .8]);
    h.FontSize = 14;
    h.FontWeight = 'normal';

    h = suplabel({'Visual and optostim conditions', ['     V0                                             V0O0                                          V0O90' ...
        '                                          V90                                            V90O0                                      V90O90']}, ...
        't', [.1 .1 .83 .87]);
    h.FontSize = 14;
    h.FontWeight = 'normal';
end

function plotNeurometricData(bitmapData, imagingData, cluster, blockNo, columnarResp, blockID, condIDs, TSbaseline,...
    TSopto, blockColumns, neuroStruct, colorCon, colorIncon,...
    optostimAreaFlag, plotType)
    snrMask=imagingData.snrMaskCoreg(:,:,blockID);
    if plotType == 1 % Sanity Check Plot (Stage 5a)
        figure('Name', 'Neurometrics')
        yLabelStr = 'Projection amplitude';
        yLimit = [-3 3];
        yTicks = -3:3;
        xLimit = [-100 100];
        xTicks = -100:25:100;
        
    elseif plotType == 2 % Projection into Visual/Reference Space (Stage 5b)
        figure('Name', ['C' num2str(cluster) 'N' num2str(blockNo) ': Projection into Visual/Reference Space'])
    else
        error('Invalid plotType. Choose 1 for Sanity Check or 2 for Projection into Visual/Reference Space.');
    end

    nCols = 3;
    nRows = 1;
    colormap(fireice);
    gap = 0.05;
    marginV = .15;
    marginH = .1;
    [hAx, ~] = tight_subplot(nRows, nCols, [gap + 0.025 gap], [marginV marginV + .06], [marginH + .06 marginH]);
    mkrSize = 15;
    lineWidth = 3;
    visualOrientation = NaN;

    % Calculate similarities and amplitudes
    referenceImage = imagingData.refmapCoreg(:,:,1:2,blockID).*snrMask;
    [~, amplitude, similarityRef, similarityData, ~, ~] = calculateProjectionVector(...
        bitmapData, imagingData, columnarResp, referenceImage,  blockID, condIDs, visualOrientation, ...
        snrMask, 'similarity');
    [~, amplitudeColumn, similarityRef, similarityData, ~, ~] = calculateProjectionVector(...
        bitmapData, imagingData, columnarResp, referenceImage,  blockID, condIDs, visualOrientation, ...
        snrMask, 'similarityColumn');

    % Loop over each map
    for plotID = 1:3
        yZeroData = [];
        
        if plotID == 1
            % Similarity over all pixels with weighted projection
            axes(hAx(plotID))
            mapName = '{Reference map}';
            similarityValues=amplitude;
            yLabelStr = 'Similarity (%)';
            % Plot axis limits
            yLimit = [-100 100];
            yTicks = -100:25:100;
            xLimit = [-100 100];
            xTicks = -100:12.5:100;
            modifier=[-1;1];
        elseif plotID == 2
            % Assumes opto dominant
            axes(hAx(plotID))
            mapName = '{Iso-tuned columns wrt. opto}';
            yLabelStr = '';

            % Recruitment            
            similarityValues=[amplitudeColumn(1,condIDs.V0),amplitudeColumn(2,condIDs.V90),...
                amplitudeColumn(1,condIDs.V0O0),amplitudeColumn(2,condIDs.V0O90),...
                amplitudeColumn(1,condIDs.V90O0),amplitudeColumn(2,condIDs.V90O90)];
            
            % Generate sorting indices based on condIDs fields
            sortedIndices = [condIDs.V0, condIDs.V90, condIDs.V0O0, condIDs.V0O90, condIDs.V90O0, condIDs.V90O90];
            % Sort similarityValues to preserve the order in amplitude/images
            [~, sortOrder] = sort(sortedIndices);  % Get sorting order
            similarityValues = similarityValues(sortOrder);
            
            % Plot axis limits
            yLimit = [-1 1];
            yTicks = [-1:.2:1];
            xLimit = [0 100];
            xTicks = 0:12.5:100;
            modifier=[-1;1];
        elseif plotID == 3
            % Assumes opto dominant
            axes(hAx(plotID))
            mapName = '{Ortho-tuned columns wrt. opto}';
            yLabelStr = '';

            % Suppression
            similarityValues=[amplitudeColumn(2,condIDs.V0),amplitudeColumn(1,condIDs.V90),...
                amplitudeColumn(2,condIDs.V0O0),amplitudeColumn(1,condIDs.V0O90),...
                amplitudeColumn(2,condIDs.V90O0),amplitudeColumn(1,condIDs.V90O90)];
           
            % Generate sorting indices based on condIDs fields
            sortedIndices = [condIDs.V0, condIDs.V90, condIDs.V0O0, condIDs.V0O90, condIDs.V90O0, condIDs.V90O90];
            % Sort similarityValues to preserve the order in amplitude/images
            [~, sortOrder] = sort(sortedIndices);  % Get sorting order
            similarityValues = similarityValues(sortOrder);
            % Plot axis limits
            yLimit = [-1 1];
            yTicks = [-1:.2:1];
            xLimit = [0 100];
            xTicks = 0:12.5:100;
            modifier=[-1;1];
        end

        % Store y-values for x=0 by marker type
        yCircleValues = [];
        yNonCircleValues = [];

        % Loop over each condition set
        for condSet = 1:3
            % Define conditions to plot
            switch condSet
                case {1} % BL V-only
                    conds = [condIDs.V0; condIDs.V90];
                    condMarker = 'o';
                    condMarkerColor = 'w';
                    condLineColor = 'k';
                    offsetBL = min(conds(:)) - 1;
                    condsOffset = conds - offsetBL;
                    gaborOrientations = TSbaseline(1).Header.Conditions.GaborOrt(TSbaseline(1).Header.Conditions.TypeCond == 3);
                    stimContrasts = TSbaseline(1).Header.Conditions.StimCon(TSbaseline(1).Header.Conditions.TypeCond == 3);
                    stimContrasts(gaborOrientations == 0) = -abs(stimContrasts(gaborOrientations == 0));
                    stimContrasts(gaborOrientations == 90) = abs(stimContrasts(gaborOrientations == 90));
                    
                    if plotType == 1
                        xData = stimContrasts(condsOffset);
                        yData = similarityValues(conds);
                        xData = [fliplr(xData(1, :)), xData(2, :)];
                        yData = [fliplr(yData(1, :)), yData(2, :)];
                    elseif plotType == 2
                        xData = stimContrasts(condsOffset(2,:));
                        yData = similarityValues(conds);
                        yData = mean(yData .* modifier);
                    end
                    condsetStr = 'Visual-only';

                    % mergeZeros
                    zeroIdx = find(xData == 0);
                    if sum(zeroIdx) > 0
                        yData(zeroIdx) = mean(yData(zeroIdx));
                        xData(zeroIdx(1)) = [];
                        yData(zeroIdx(1)) = [];
                    end
                case {2} % Con V+O / O0 V+O
                    if plotType == 1
                        conds = [condIDs.V0O0; condIDs.V90O0];
                        condMarker = 'v';
                        colorRow = [156, 14, 254] / 255;
                    elseif plotType == 2
                        conds = [condIDs.V0O0; condIDs.V90O90];
                        condMarker = '^';
                        colorRow = [0, 225, 80] / 255;
                    end
                    condMarkerColor = colorRow;
                    condLineColor = colorRow;
                    gaborOrientations = TSopto(1).Header.Conditions.GaborOrt(TSopto(1).Header.Conditions.TypeCond == 3);
                    stimContrasts = TSopto(1).Header.Conditions.StimCon(TSopto(1).Header.Conditions.TypeCond == 3);
                    stimContrasts(gaborOrientations == 0) = -abs(stimContrasts(gaborOrientations == 0));
                    stimContrasts(gaborOrientations == 90) = abs(stimContrasts(gaborOrientations == 90));
                    
                    if plotType == 1
                        xData = stimContrasts(conds);
                        yData = similarityValues(conds);
                        xData = [fliplr(xData(1, :)), xData(2, :)];
                        yData = [fliplr(yData(1, :)), yData(2, :)];
                    elseif plotType == 2
                        xData = stimContrasts(conds(2,:));
                        yData = similarityValues(conds);
                        yData = mean(yData .* modifier);
                    end
                    
                    if plotType == 1
                        condsetStr = 'Visual + opto-0';
                    elseif plotType == 2
                        condsetStr = 'Visual + con-opto';
                        zeroConstrast=find(xData==0);
                        if ~isempty(zeroConstrast)
                            yZeroData=[yZeroData, yData(zeroConstrast)];
                        else
                            yZeroData=[yZeroData, NaN];
                        end
                    end

                case {3} % Incon V+O / O90 V+O
                    if plotType == 1
                        conds = [condIDs.V0O90; condIDs.V90O90];
                        condMarker = '^';
                        colorRow = [0, 225, 80] / 255;
                    elseif plotType == 2
                        conds = [condIDs.V0O90; condIDs.V90O0];
                        condMarker = 'v';
                        colorRow = [156, 14, 254] / 255;
                    end

                    condMarkerColor = colorRow;
                    condLineColor = colorRow;
                    gaborOrientations = TSopto(1).Header.Conditions.GaborOrt(TSopto(1).Header.Conditions.TypeCond == 3);
                    stimContrasts = TSopto(1).Header.Conditions.StimCon(TSopto(1).Header.Conditions.TypeCond == 3);
                    stimContrasts(gaborOrientations == 0) = -abs(stimContrasts(gaborOrientations == 0));
                    stimContrasts(gaborOrientations == 90) = abs(stimContrasts(gaborOrientations == 90));
                    
                    if plotType == 1
                        xData = stimContrasts(conds);
                        yData = similarityValues(conds);
                        xData = [fliplr(xData(1, :)), xData(2, :)];
                        yData = [fliplr(yData(1, :)), yData(2, :)];
                    elseif plotType == 2
                        xData = -stimContrasts(conds(1,:));
                        yData = similarityValues(conds);
                        yData = mean(yData .* modifier);
                    end

                    if plotType == 1
                        condsetStr = 'Visual + opto-90';
                    elseif plotType == 2
                        condsetStr = 'Visual + incon-opto';
                        zeroConstrast=find(xData==0);
                        if ~isempty(zeroConstrast)
                            yZeroData=[yZeroData, yData(zeroConstrast)];
                        else
                            yZeroData=[yZeroData, NaN];
                        end
                    end

            end

            % Plot the data
            plot(xData, yData, '-', 'LineWidth', lineWidth, ...
                'Marker', condMarker, 'MarkerSize', mkrSize, 'MarkerFaceColor', condMarkerColor, ...
                'MarkerEdgeColor', 'k', 'Color', condLineColor, 'DisplayName', condsetStr);
            hold on;
            yline(0, '--', 'Color', [.4 .4 .4], 'LineWidth', 1.5, 'HandleVisibility', 'off');
            
            if plotType == 1
                xline(0, '--', 'Color', [.4 .4 .4], 'LineWidth', 1.5, 'HandleVisibility', 'off');
            end

            % Collect y-values at x=0 for calculating neuroStruct.phi
            if plotType == 1
                yAtZero = handleZeroData(xData, yData);
                % Append the calculated value based on condMarker
                if ~isempty(yAtZero)
                    if strcmp(condMarker, 'o')
                        yCircleValues = [yCircleValues, yAtZero];
                    else
                        yNonCircleValues = [yNonCircleValues, yAtZero];
                    end
                else
                    warning('Could not calculate yData at xData == 0 for this block.');
                end
            end

            % Average
            if plotType == 1
                neuroStruct.(['C' num2str(cluster)]).visualmap.averageProjection(blockNo, condSet) = mean(yData(:));
            elseif plotType == 2
                neuroStruct.(['C' num2str(cluster)]).refmap.averageProjection(blockNo, condSet) = mean(yData(:));
                % Plot mean
                meanLen = 7; % The half-width of the line segment
                % Plot the line from (xVal - halfW, yVal) to (xVal + halfW, yVal)
                plot([90, 90 + meanLen], [mean(yData(:)), mean(yData(:))], ...
                    'LineWidth', 4, 'Color', condMarkerColor, 'HandleVisibility', 'off');
            end
        end

        % Compute neuroStruct.phi (difference between means of y-values at x=0)
        if plotType == 1
            if ~isempty(yCircleValues) && ~isempty(yNonCircleValues)
                neuroStruct.(['C' num2str(cluster)]).visualmap.phi(blockNo) = mean(yNonCircleValues - yCircleValues);
            else
                neuroStruct.(['C' num2str(cluster)]).visualmap.phi(blockNo) = nan;
            end
        elseif plotType == 2
             neuroStruct.(['C' num2str(cluster)]).refmap.phi(blockNo) = yZeroData(1)-yZeroData(2);
        end

        % Set plot limits and labels
        
        if plotType == 1
            phi = round(neuroStruct.(['C' num2str(cluster)]).visualmap.phi(blockNo), 2, 'significant');
            meanProjection = round(neuroStruct.(['C' num2str(cluster)]).visualmap.averageProjection(blockNo, :), 2, 'significant');
        elseif plotType == 2
            phi = round(neuroStruct.(['C' num2str(cluster)]).refmap.phi(blockNo), 2, 'significant');
            meanProjection = round(neuroStruct.(['C' num2str(cluster)]).refmap.averageProjection(blockNo, :), 2, 'significant');
        end
        
        titleStr = createColoredTitle(mapName, meanProjection, phi, colorCon, colorIncon);
        title({titleStr, ...
            ['Columns: ' num2str(blockColumns(1)) ', ' num2str(blockColumns(2))]}, 'FontWeight', 'normal', 'Interpreter', 'tex')
        xlabel('Visual contrast (%)');
        ylabel(yLabelStr);
        axis square;
        upFontSize(20, .01);
        
        ylim(yLimit);
        yticks(yTicks)
        xlim(xLimit);
        xticks(xTicks);
        
        if plotType == 2
            %addSkippedTicks(-100,100,25, 'y');
            addSkippedTicks(0, 100, 12.5, 'x');
        elseif plotType == 1
            addSkippedTicks(-100, 100,25, 'x');
            addSkippedTicks(-100,100,25, 'y');
        end
        if plotID==1
            legend('location','southeast')
        end

        %% Add inset of subtracted image
        % Choose the area analyzed
        
        if plotType == 2 && plotID==1
            gcf;
            if plotID == 1
                axInset = axes('Position', [.1 + .043 .257 .18 .18]);
            elseif plotID == 2
                axInset = axes('Position', [.35 + .043 .257 .18 .18]);
            end
            imgsc(columnarResp(:,:,condIDs.V90O90(1)).*snrMask - columnarResp(:,:,condIDs.V0O0(1)).*snrMask);
            hold on;
            if strcmp(optostimAreaFlag, 'nonstim')
                if mean(blockColumns) >= 15
                    cax(.8);
                elseif mean(blockColumns) >= 3 && mean(blockColumns) < 15
                    cax(.6);
                else
                    cax(.3);
                end
            else
                cax();
            end
            box off;
            set(axInset, 'XTick', [], 'YTick', [], 'Box', 'off', 'XColor', 'none', 'YColor', 'none', 'Color', 'none'); % No ticks, no box
        end
        
        
    end

    % Colorbar struff
    if plotType == 2
        hFig = gcf();
        allKids = hFig.Children;
        for k = 1:numel(allKids)
            if isa(allKids(k), 'matlab.graphics.illustration.ColorBar')
                allKids(k).FontSize = 10;
            end
        end
    end
end

%% Saving code
function saveFigures(saveFlag, mainPath, monkey, cluster, date, run, optostimAreaFlag, fileSuffix)
    if saveFlag
        % Define paths for different file types
        pngPath = [mainPath monkey '\Meta\neurometric\projections\png\'];
        svgPath = [mainPath monkey '\Meta\neurometric\projections\svg\'];
        figPath = [mainPath monkey '\Meta\neurometric\projections\fig\'];

        % Create directories if they don't exist
        if ~exist(pngPath, 'dir')
            mkdir(pngPath);
        end
        if ~exist(svgPath, 'dir')
            mkdir(svgPath);
        end
        if ~exist(figPath, 'dir')
            mkdir(figPath);
        end

        % Set font and vector renderer
        set(findall(gcf, '-property', 'FontName'), 'FontName', 'SansSerif');
        set(gcf, 'Renderer', 'painters'); % Use painters for vector graphics

        % Save the figure in different formats
        fullFigFilename=['C' num2str(cluster) 'M28D' date 'R' run optostimAreaFlag fileSuffix];
        exportgraphics(gcf, [pngPath fullFigFilename '.png'], 'Resolution', 600); % High-res PNG
        savefig(gcf, [figPath fullFigFilename '.fig']); % MATLAB figure file
        print(gcf, [svgPath fullFigFilename '.svg'], '-dsvg'); % Vector graphics SVG
        %exportgraphics(gcf, [svgPath fullFigFilename '.svg'], 'ContentType', 'vector'); % Vector graphics SVG
        disp(['Figure saved as: ', fullFigFilename]);
    else
        disp('saveFlag is 0, figure not saved.');
    end
end



function plotNeurometricDataSlice(imagingData, bitmapData, cluster, blockNo, blockID, columnarResp, condIDs, TSbaseline,...
    TSopto, blockColumns, neuroStruct, colorCon, colorIncon,...
    optostimAreaFlag, plotType)

    figure('Name','Similarity slice')
    nCols = 4;
    nRows = 2;
    colormap(fireice);
    gap = 0.05;
    marginV = .1;
    marginH = .05;
    [hAx, ~] = tight_subplot(nRows, nCols, [gap + 0.025 gap], [marginV marginV], [marginH marginH]);
    mkrSize = 15;
    lineWidth = 3;

    % Loop over each map
    for map = 2%1:2
        visualOrientation = NaN;
        yZeroData = [];

        if map == 1
            axes(hAx(map))
            mapName = '{Visual map}';
            referenceImage = rescale(columnarResp(:,:,condIDs.V90(end)), -1, 1) - ...
                rescale(columnarResp(:,:,condIDs.V0(end)), -1, 1);
            snrMask= imagingData.snrMaskCoreg(:,:,blockID);
            [~, amplitude, similarityRef, similarityData, ~, ~] = calculateProjectionVector(columnarResp, referenceImage, condIDs, visualOrientation, snrMask);
        elseif map == 2
            axes(hAx(map))
            mapName = '{Reference map}';
            snrMask= imagingData.snrMaskCoreg(:,:,blockID);
            snrMask(snrMask==0)=nan;
            referenceImage = imagingData.refmapCoreg(:,:,1:2,blockID).*snrMask;
            referenceImage0=nan(size(referenceImage)); referenceImage90=nan(size(referenceImage));
            referenceImage0(referenceImage0<0)=1;referenceImage0=referenceImage0.*snrMask;
            referenceImage90(referenceImage90>0)=1;referenceImage90=referenceImage90.*snrMask;

            roiMask=double(~bitmapData.bitmapMask(:,:,blockID));
            roiMask(roiMask==0)=nan;
            [~, amplitude, similarityRef, similarityData, ~, ~] = calculateProjectionVector(columnarResp, referenceImage, condIDs, visualOrientation, snrMask);
        end
        
        % Initialize result vectors
        numPairs = numel(condIDs.V0O0);  % Assuming condIDs.V0O0 and condIDs.V90O90 have the same length
        dataO0C0 = nan(numPairs,1);
        dataO0C90 = nan(numPairs,1);
        dataO90C0 = nan(numPairs,1);
        dataO90C90 = nan(numPairs,1);
        
        % Loop through each condition pair
        for i = 1:numPairs
            % Define data for current iteration
            respO0 = columnarResp(:,:,condIDs.V0O0(i)) .* snrMask;
            respO90 = columnarResp(:,:,condIDs.V90O90(i)) .* snrMask;
            respO0ROI = respO0 .* roiMask;
            respO90ROI = respO90 .* roiMask;
        
            % Reference map processing
            tileGridSize = [32,32];
            gammaVal = 5;
            sensitivity = 1e-3;  % Example sensitivity value
        
            [~, referenceImageColumn0, ~, ~] = processActivityMap(...
                rescale(imagingData.refmapCoreg(:,:,1,blockID), 0, 1), tileGridSize, gammaVal, sensitivity);
            [~, referenceImageColumn90, ~, ~] = processActivityMap(...
                rescale(imagingData.refmapCoreg(:,:,2,blockID), 0, 1), tileGridSize, gammaVal, sensitivity);
        
            referenceImageColumn0 = referenceImageColumn0 .* snrMask;
            referenceImageColumn90 = referenceImageColumn90 .* snrMask;
        
            % Column masks
            columns0 = ones(size(referenceImage));
            columns0(referenceImage >= 0) = NaN;
            columns0 = columns0 .* snrMask .* roiMask;
        
            columns90 = ones(size(referenceImage));
            columns90(referenceImage <= 0) = NaN;
            columns90 = columns90 .* snrMask .* roiMask;
        
            % Compute responses
            respO0C0 = respO0ROI .* referenceImageColumn0; respO0C0(respO0C0==0) = nan;
            respO0C90 = respO0ROI .* referenceImageColumn90; respO0C90(respO0C90==0) = nan;
            respO90C0 = respO90ROI .* referenceImageColumn0; respO90C0(respO90C0==0) = nan;
            respO90C90 = respO90ROI .* referenceImageColumn90; respO90C90(respO90C90==0) = nan;
        
            % Store results as vectors
            dataO0C0(i) = mean(removenan(respO0C0));
            dataO0C90(i) = mean(removenan(respO0C90));
            dataO90C0(i) = mean(removenan(respO90C0));
            dataO90C90(i) = mean(removenan(respO90C90));
        end
        
        % Store all responses in a 4D matrix
        superRespO090C090 = cat(4, dataO0C0, dataO0C90, dataO90C0, dataO90C90);

        % Figure
        figure('Name','Similarity slice')
        nCols = 3;
        nRows = 1;
        colormap(fireice);
        gap = 0.05;
        marginV = .1;
        marginH = .05;
        [hAx, ~] = tight_subplot(nRows, nCols, [gap + 0.025 gap], [marginV marginV], [marginH marginH]);
        mkrSize = 15;
        lineWidth = 3;

        axes(hAx(1))
        imgsc(referenceImage)
        axes(hAx(2))
        imgsc(referenceImageColumn0)
        axes(hAx(3))
        imgsc(referenceImageColumn90)
            
 % Figure
        figure('Name','Similarity slice')
        nCols = 4;
        nRows = 1;
        colormap(fireice);
        gap = 0.05;
        marginV = .1;
        marginH = .05;
        [hAx, ~] = tight_subplot(nRows, nCols, [gap + 0.025 gap], [marginV marginV], [marginH marginH]);
        mkrSize = 15;
        lineWidth = 3;
        validValues = superRespO090C090(~isnan(superRespO090C090));
        climits = max([max(validValues, [], 'all'), abs(min(validValues, [], 'all'))]);
        axes(hAx(1))
        imgsc(respO0C0); title(['Con: O0C0 (mean=' num2str(nanmean(respO0C0(:)),2) ')']); cax(1,[-climits climits])
        axes(hAx(2))
        imgsc(respO0C90); title(['Incon: O0C90 (mean=' num2str(nanmean(respO0C90(:)),2) ')']);  cax(1,[-climits climits])
        axes(hAx(3))
        imgsc(respO90C0); title(['Incon: O90C0 (mean=' num2str(nanmean(respO90C0(:)),2) ')']);  cax(1,[-climits climits])
        axes(hAx(4))
        imgsc(respO90C90); title(['Con: O90C90 (mean=' num2str(nanmean(respO90C90(:)),2) ')']);  cax(1,[-climits climits])

        figure('Name','Similarity slice')
        nCols = 4;
        nRows = 2;
        colormap(fireice);
        gap = 0.05;
        marginV = .1;
        marginH = .05;
        [hAx, ~] = tight_subplot(nRows, nCols, [gap + 0.05 gap], [marginV marginV], [marginH marginH]);
        mkrSize = 15;
        lineWidth = 3;
        % Raw
        axes(hAx(1))
        imgsc(respO0);cax; title('Opto 0')
        upFontSize; addPix2MM(respO0,1,nRows,nCols)

        axes(hAx(5))
        imgsc(respO90);cax; title('Opto 90')
        upFontSize;addPix2MM(respO0,5,nRows,nCols)
        
        % Similarity Img
        axes(hAx(2))
        imgsc(respO0C0);cax; title('Opto 0, Column 0')
        upFontSize;addPix2MM(respO0,2,nRows,nCols)

        axes(hAx(6))
        imgsc(respO90C0);cax; title('Opto 90, Column 0')
        upFontSize;addPix2MM(respO0,6,nRows,nCols)
        
        % Similarity Img
        axes(hAx(3))
        imgsc(respO0C90);cax; title('Opto 0, Column 90')
        upFontSize;addPix2MM(respO0,3,nRows,nCols)

        axes(hAx(7))
        imgsc(respO90C90);cax; title('Opto 90, Column 90')
        upFontSize;addPix2MM(respO0,7,nRows,nCols)

        axes(hAx(4))
        % Combine data
        conData=[dataO0C0; dataO0C90];
        conLabel=[ones(size(dataO0C0));2*ones(size(dataO0C90))];
        colors = [[0, 225, 80] / 255; [156, 14, 254] / 255];  % Green con, Purple incon\
        plot(conData)
        %{
        h = daboxplot(conData, 'groups', conLabel, 'colors', colors, ...
            'mean' ,1, 'outliers', 0, 'whiskers', 1,'scatter',0, 'jitter',0, 'xtlabels', {'0-col','90-col'});
        title('Distributions (Opto 0)','FontWeight','normal')
        %}
        upFontSize; 

        axes(hAx(8))
        % Combine data
        inconData=[dataO90C0; dataO90C90];
        inconLabel=[ones(size(dataO90C0));2*ones(size(dataO90C90))];
        colors = [[156, 14, 254] / 255; [0, 225, 80] / 255];  % Green con, Purple incon
        h = daboxplot(inconData, 'groups', inconLabel, 'colors', colors, ...
            'mean' ,1, 'outliers', 0, 'whiskers', 1,'scatter',0, 'jitter',0, 'xtlabels', {'0-col','90-col'});
        title('Distributions (Opto 90)','FontWeight','normal')
        upFontSize;

        % Similarity slice  amp
        sliceAmpRef=similarityRef.sliceAmp;
        sliceAmpV0O0=similarityData.sliceAmp(:,:,condIDs.V0O0(1));
        sliceAmpV90O90=similarityData.sliceAmp(:,:,condIDs.V90O90(1));






        
        xIdx = -numel(sliceAmpV90O90)/2+1 : 1 : numel(sliceAmpV90O90)/2;
        refSliceV0 = double(nanmean(-sliceAmpRef,2)>0);
        
        % Shaded bars for refSliceV0
        for i = 1:length(refSliceV0)
            if refSliceV0(i) > 0
                shadePlot([xIdx(i)-0.5, xIdx(i)+0.5], [0, 1], 'b', 0.2); % Adjust y-bounds and alpha as needed
            end
        end
        
        plot(xIdx,-sliceAmpV0O0,'k'); axis square; cax; xlim([-numel(sliceAmpV90O90)/2+1, numel(sliceAmpV90O90)/2]); ylim([-.3 1]); hold on;
        
        axes(hAx(8))
        refSliceV90 = double(nanmean(sliceAmpRef,2)>0);
        
        % Shaded bars for refSliceV90
        for i = 1:length(refSliceV90)
            if refSliceV90(i) > 0
                shadePlot([xIdx(i)-0.5, xIdx(i)+0.5], [0, 2], 'r', 0.2); % Adjust y-bounds and alpha as needed
            end
        end
        
        plot(xIdx,sliceAmpV90O90,'k'); axis square; cax; xlim([-numel(sliceAmpV90O90)/2+1, numel(sliceAmpV90O90)/2]); ylim([-.3 2]); hold on;
    end
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
    procMapThresh(~validMask) = NaN;  % keep invalid regions as NaN

    %% 4) PREPARE THE PROCESSED MAP FOR REGISTRATION
    %    We fill its NaNs with 0 so imregtform doesn't crash.
    activityMapProc = procMapGamma;
    nanMask = isnan(activityMapProc);
    activityMapProc(nanMask) = 0;
    Rfixed = imref2d(size(activityMapProc));
    %activityMapProc=rescale(activityMapProc);

end