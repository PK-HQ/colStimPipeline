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
        for blockNo=18%1:numel(clusterBlocks)%[3 5 10 11 14]
            tic
            blockID=clusterBlocks(blockNo);
            imagingData.optoIntg=[];
            imagingData.baselineIntg=[];
    
            % Load block data individually as imaging data too large
             [currentBlockStruct,~,...
                behavioralData, imagingData, ~, ~]=loadBlockData(datastruct, analysisBlockID, behavioralData, imagingData, bitmapData, blockID,'');
                filenameStructCurrent=generateFilenames(currentBlockStruct);
            %{
            if size(imagingData.optoIntg,3)>=blockID
                if isempty(imagingData.optoIntg(:,:,:,blockID))
                    disp(['Check INTG: D' datastruct(analysisBlockID(blockID)).date 'R' datastruct(analysisBlockID(blockID)).run])
                    errorIntgBlocks{end+1}= ['D' datastruct(analysisBlockID(blockID)).date 'R' datastruct(analysisBlockID(blockID)).run];
                    continue
                end
            elseif size(imagingData.optoIntg,3)<blockID
                disp(['Check INTG: D' datastruct(analysisBlockID(blockID)).date 'R' datastruct(analysisBlockID(blockID)).run])
                errorIntgBlocks{end+1}= ['D' datastruct(analysisBlockID(blockID)).date 'R' datastruct(analysisBlockID(blockID)).run];
                continue
            end
            %}
            monkey=datastruct(analysisBlockID(blockID)).monkey;
            date=datastruct(analysisBlockID(blockID)).date;
            run=datastruct(analysisBlockID(blockID)).run;
            pdfFilename=filenameStructCurrent.recruitmentPDF;
            blockColumns=bitmapData.nColumns(:,blockID);
            %% STAGE 2. For each intg image, split image into stimulated and nonstimulated area with gaussian ROI imagingData.mask
            %% Get behavioral and neural responses from valid completed trials
            % Process integrated image to get the average image per condition (successfully completed trials only)
            nTS=1;

            if isfield(behavioralData,'baselineTS')
                if ~isempty(behavioralData.baselineTS(blockID).Header)
                    nTS=2;
                end
            else
                nTS=1;
            end
            switch nTS
                case {1}
                    % Split by cond, blanks subtracted and removed
                    if isempty(imagingData.optoIntg)
                        continue
                    else
                        [~,images,condIDs]=getUsableTrials(behavioralData.optoTS(blockID), imagingData.optoIntg(:,:,:,blockID));
                    end
    
                   % rename TSopto and imagesOpto to TSstruct and images for rest of the script
                   %images=imagingData.optoIntg(:,:,:,blockID);
                    TSbaseline=behavioralData.optoTS(blockID);
                    TSopto=behavioralData.optoTS(blockID);

                    condIDs=condIDs;
                case {2}
                    % Here the baseline and opto are in separate blocks, so extract TSstruct and Datatrial
                    % Split by cond, blanks subtracted and removed
                    if isempty(imagingData.baselineIntg) || isempty(imagingData.optoIntg)
                        continue
                    else
                        [~,imagesbl,condIDsbl]=getUsableTrials(behavioralData.baselineTS(blockID), imagingData.baselineIntg(:,:,:,blockID));
                        [~,imagesOpto,condIDsOpto]=getUsableTrials(behavioralData.optoTS(blockID),imagingData.optoIntg(:,:,:,blockID));
                    end
                    fields=fieldnames(condIDsOpto);
                    for fn=1:numel(fields)
                       name=fields{fn};
                      % reassign
                       if ~isempty(condIDsOpto.(name))
                           condsOpto(fn)=max(condIDsOpto.(name));
                       else
                           condsOpto(fn)=NaN;
                       end
                    end
                    fields=fieldnames(condIDsbl);
                    for fn=1:numel(fields)
                       name=fields{fn};
                      % reassign
                       if ~isempty(condIDsbl.(name))
                           condsBL(fn)=max(condIDsbl.(name));
                       else
                           condsBL(fn)=NaN;
                       end
                    end
    
                    % Append baseline condIDs and images to the back of opto trials
                    maxOptoConds=max(condsOpto);
                    maxBLConds=max(condsBL);
                    BLCondIdx=maxOptoConds+1:maxOptoConds+maxBLConds;
                    
                    %condIDsOpto.baselineConds=condIDsbl.baselineConds + maxOptoConds;
                    
                    condIDsOpto.V0=condIDsbl.V0 + maxOptoConds;
                    condIDsOpto.V90=condIDsbl.V90 + maxOptoConds;
                    %imagesOpto.trials(:,:,BLCondIdx,1:size(imagesbl.trials,4))=imagesbl.trials;
                    %imagesOpto.trialsCorrect(:,:,BLCondIdx,1:size(imagesbl.trialsCorrect,4))=imagesbl.trialsCorrect;
                    %imagesOpto.trialsError(:,:,BLCondIdx,1:size(imagesbl.trialsError,4))=imagesbl.trialsError;
                    imagesOpto.average(:,:,BLCondIdx)=imagesbl.average;
                    imagesOpto.averageCorrect(:,:,BLCondIdx)=imagesbl.averageCorrect;
                    imagesOpto.averageError(:,:,BLCondIdx)=imagesbl.averageError;
    
                    % rename TSopto and imagesOpto to TSstruct and images for rest of the script
                    images=imagesOpto;
                    TSbaseline=behavioralData.baselineTS(blockID);
                    TSopto=behavioralData.optoTS(blockID); 
                    condIDs=condIDsOpto;
            end
    
            [images.averageColumn]=columnarFilter(behavioralData.optoTS(blockID),images,trialOutcomeType);
            %% Get reference PCA map for comparison
            % Process make ROI imagingData.mask
            ROImask=double(imagingData.mask(:,:,blockID)); snrMask=ROImask; snrMask(snrMask==0)=NaN;
            gaussCond=datastruct(analysisBlockID(blockID)).gaussianCond;
            gaussMask=imagingData.gaussfit(:,:,blockID); %rename
            % imagingData.mask ort map, apply vasculature transform (morphs ort map to the imaging window view from experiment)
            PCAmap=imagingData.ortpca(:,:,:,blockID) .* snrMask;
            %PCAmap=transformImage(PCAmap,filenameStructCurrent);
    
            %% STAGE 3. Split the image into optostim and recruitment ROI
            meanColumns=mean(blockColumns);
            [fullROI, optostimROI, recruitROI]=parcellation(condIDs,currentBlockStruct,bitmapData,imagingData,...
                images,trialOutcomeType,snrMask,gaussMask, bitmapData.columnarbitmapTFcamspace,meanColumns,...
                optostimMaskType,blockID);

            [fullROI, optostimROI, recruitROI] = parcellation(condIDs, currentBlockStruct, bitmapData, imagingData, ...
            images, trialOutcomeType, gaussMask, bitmapData.columnarbitmapTFcamspace, meanColumns, ...
            optostimMaskType, blockID);
            
            % Choose the area analyzed
            switch optostimAreaFlag
                case {'full'}
                    columnarResp=fullROI;
                    titleStr='Neurometric (full ROI)';
                case {'stim'}
                    columnarResp=optostimROI;
                    titleStr='Neurometric (Stimulated ROI)';
                case {'nonstim'}
                    columnarResp=recruitROI;
                    titleStr='Neurometric (Recruitment ROI)';
            end
            
            %% STAGE 4. Plot processed images for ROI
            %{
            figure
            nCols=6;nRows=size(columnarResp,3)/nCols;
            plotIndices = generateGridIdx(nRows,nCols);
            %gap=.01;marginV=.0;marginH=.0;
            %[hAx,~]=tight_subplot(nRows,nCols,[gap gap], [marginV marginV], [marginH marginH]);
            colormap(fireice); 
            gap=.0075;marginV=.01;marginH=.01;
            [hAx,~]=tight_subplot(nRows,nCols,[gap+0.025 gap], [marginV marginV+.06], [marginH+.06 marginH]);
            
            for condSet=1:3
                contrasts=[];
                amplitudes=[];
                % Define conditions to plot
                switch condSet
                    case {1} %BL V-only
                        conds=[condIDs.V0; condIDs.V90]; % store as 2 rows, one for each set of conds
                        yData=columnarResp(:,:,conds);
                        plotIdx=plotIndices(:,[1 4])';
                        condsetStr='Visual-only';
                        superYData=yData;
                    case {2} %V0
                        conds=[condIDs.V0O0; condIDs.V0O90];
                        yData=columnarResp(:,:,conds);
                        plotIdx=plotIndices(:,[2 3])';
                        condsetStr='Visual + con/incon opto';
                        superYData=columnarResp(:,:,[condIDs.V0O0; condIDs.V0O90; condIDs.V90O90; condIDs.V90O0]);
                    case {3} %V90
                        conds=[condIDs.V90O0; condIDs.V90O90];
                        yData=columnarResp(:,:,conds);
                        plotIdx=plotIndices(:,[5 6])';
                        condsetStr='Visual + con/incon opto';
                        superYData=columnarResp(:,:,[condIDs.V0O0; condIDs.V0O90; condIDs.V90O90; condIDs.V90O0]);
                end
    
                for cond=1:numel(plotIdx)
                    axes(hAx(plotIdx(cond)))
                    imgsc(yData(:,:,cond));                
                    cax(.7,superYData)
                    upFontSize(10,.005)
                    axis square
                    colorbar
                    addPix2MM(yData(:,:,cond),nCols,nRows,plotIdx(cond))
                end
                upFontSize(12,.005)
            end
            % Headers for visual contrast and conds
            h=suplabel('Visual contrast','y',[.1 .1 .8 .8]);
            h.FontSize=14;
            h.FontWeight='normal';
    
            h=suplabel({'Visual and optostim conditions',['     V0                                             V0O0                                          V0O90' ...
                '                                          V90                                            V90O0                                      V90O90']}...
                ,'t',[.1 .1 .83 .87]);
            h.FontSize=14;
            h.FontWeight='normal';
            if saveFlag
                monkey=datastruct(analysisBlockID(blockID)).monkey;
                date=datastruct(analysisBlockID(blockID)).date;
                run=datastruct(analysisBlockID(blockID)).run;
                figFilename=[mainPath monkey '\Meta\neurometric\projections\M28D' date 'R' run 'averageconds'];
                set(gcf, 'Renderer', 'painters'); % Use painters for vector graphics
                set(findall(gcf, '-property', 'FontName'), 'FontName', 'SansSerif');
                print(gcf, [figFilename '.png'], '-dpng', '-r600'); % High-res PNG
                savefig(gcf, [figFilename '.fig']);           % FIG
                print(gcf, [figFilename '.svg'], '-dsvg');        % SVG
            end
            %}
            %% Stage 5a. Sanity check raw plot
            %{
            figure('Name','Neurometrics')
            nCols=2; nRows=1;
            colormap(fireice); 
            gap=0; marginV=.25; marginH=.2;
            [hAx,~]=tight_subplot(nRows,nCols,[gap+0.025 gap], [marginV marginV+.06], [marginH+.06 marginH]);
            mkrSize=15; lineWidth=3;
    
            % Loop over each map
            for map=1:2
                visualOrientation=NaN; 
                if map==1
                    visualOrientation=NaN; 
                    axes(hAx(map))
                    mapName='{Visual map}';
                    referenceImage=rescale(columnarResp(:,:,condIDs.V90(end)),-1,1) - ...
                                  rescale(columnarResp(:,:,condIDs.V0(end)),-1,1); 
                    [~, amplitude, ~, ~] = calculateProjectionVector(columnarResp, referenceImage, condIDs, visualOrientation, snrMask);
                    maxY = 3;
                elseif map==2
                    axes(hAx(map))
                    mapName='{Reference map}';
                    referenceImage=rescale(PCAmap(:,:,7),0,1) - rescale(PCAmap(:,:,1),0,1);
                    [~, amplitude, ~, ~] = calculateProjectionVector(columnarResp, referenceImage, condIDs, visualOrientation, snrMask);
                    maxY = 5;
                end
                
                % Store y-values for x=0 by marker type
                yCircleValues = [];
                yNonCircleValues = [];
                
                % Loop over each condition set
                for condSet=1:3
                    contrasts=[];
                    amplitudes=[];
                    
                    % Define conditions to plot
                    switch condSet
                        case {1} % BL V-only
                            conds=[condIDs.V0; condIDs.V90];
                            condMarker = 'o';
                            condMarkerColor = 'w';
                            condLineColor = 'k';
                            offsetBL = min(conds(:)) - 1;
                            condsOffset = conds - offsetBL;
                            gaborOrientations = TSbaseline(1).Header.Conditions.GaborOrt(TSbaseline(1).Header.Conditions.TypeCond == 3);
                            stimContrasts = TSbaseline(1).Header.Conditions.StimCon(TSbaseline(1).Header.Conditions.TypeCond == 3);
                            stimContrasts(gaborOrientations == 0) = -abs(stimContrasts(gaborOrientations == 0));
                            stimContrasts(gaborOrientations == 90) = abs(stimContrasts(gaborOrientations == 90));
                            xData = stimContrasts(condsOffset);
                            yData = amplitude(conds);
                            xData = [fliplr(xData(1, :)), xData(2, :)];
                            yData = [fliplr(yData(1, :)), yData(2, :)];
                            condsetStr = 'Visual-only';
                            % mergeZeros
                            zeroIdx=find(xData==0);
                            if sum(zeroIdx)>0
                                yData(zeroIdx)=mean(yData(zeroIdx));
                                xData(zeroIdx(1))=[];yData(zeroIdx(1))=[];
                            end
                        case {2} % O0 V+O
                            conds=[condIDs.V0O0; condIDs.V90O0];
                            condMarker = 'v';
                            colorRow = [156, 14, 254] / 255;
    
                            condMarkerColor = colorRow;
                            condLineColor = colorRow;
                            gaborOrientations = TSopto(1).Header.Conditions.GaborOrt(TSopto(1).Header.Conditions.TypeCond == 3);
                            stimContrasts = TSopto(1).Header.Conditions.StimCon(TSopto(1).Header.Conditions.TypeCond == 3);
                            stimContrasts(gaborOrientations == 0) = -abs(stimContrasts(gaborOrientations == 0));
                            stimContrasts(gaborOrientations == 90) = abs(stimContrasts(gaborOrientations == 90));
                            xData = stimContrasts(conds);
                            yData = amplitude(conds);
                            xData = [fliplr(xData(1, :)), xData(2, :)];
                            yData = [fliplr(yData(1, :)), yData(2, :)];
                            condsetStr = 'Visual + opto-0';
                        case {3} % O90 V+O
                            conds = [condIDs.V0O90; condIDs.V90O90];
                            condMarker = '^';
                            colorRow = [0, 225, 80] / 255;
                            condMarkerColor = colorRow;
                            condLineColor = colorRow;
                            gaborOrientations = TSopto(1).Header.Conditions.GaborOrt(TSopto(1).Header.Conditions.TypeCond == 3);
                            stimContrasts = TSopto(1).Header.Conditions.StimCon(TSopto(1).Header.Conditions.TypeCond == 3);
                            stimContrasts(gaborOrientations == 0) = -abs(stimContrasts(gaborOrientations == 0));
                            stimContrasts(gaborOrientations == 90) = abs(stimContrasts(gaborOrientations == 90));
                            xData = stimContrasts(conds);
                            yData = amplitude(conds);
                            xData = [fliplr(xData(1, :)), xData(2, :)];
                            yData = [fliplr(yData(1, :)), yData(2, :)];
                            condsetStr = 'Visual + opto-90';
                    end
                    
                    % Plot the data
                    plot(xData, yData, '-', 'LineWidth', lineWidth, ...
                        'Marker', condMarker, 'MarkerSize', mkrSize, 'MarkerFaceColor', condMarkerColor, ...
                        'MarkerEdgeColor', 'k', 'Color', condLineColor, 'DisplayName', condsetStr); 
                    hold on;
                    yline(0, '--', 'Color', [.4 .4 .4], 'LineWidth', 1.5, 'HandleVisibility', 'off');
                    xline(0, '--', 'Color', [.4 .4 .4], 'LineWidth', 1.5, 'HandleVisibility', 'off');
                    
                    % Collect y-values at x=0 for calculating neuroStruct.phi
                    % Call the subfunction to handle the logic
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
                    %{
                    zeroIdx = xData == 0;
                    if any(zeroIdx)
                        if strcmp(condMarker, 'o')
                            yCircleValues = [yCircleValues, mean(yData(zeroIdx))];
                        else
                            yNonCircleValues = [yNonCircleValues, mean(yData(zeroIdx))];
                        end
                    end
                    %}
                    % Average
                    neuroStruct.(['C' num2str(cluster)]).visualmap.averageProjection(blockNo,condSet)=mean(yData(:));
                end
                
                % Compute neuroStruct.phi (difference between means of y-values at x=0)
                if ~isempty(yCircleValues) && ~isempty(yNonCircleValues)
                    neuroStruct.(['C' num2str(cluster)]).visualmap.phi(blockNo) = mean(yNonCircleValues - yCircleValues);
                else
                    neuroStruct.(['C' num2str(cluster)]).visualmap.phi(blockNo)=nan;
                end
    
                % Set plot limits and labels
                ylim([-maxY maxY]);
                xlim([-100 100]);
                phi=round(neuroStruct.(['C' num2str(cluster)]).visualmap.phi(blockNo),2,'significant');
                meanProjection=round(neuroStruct.(['C' num2str(cluster)]).visualmap.averageProjection(blockNo,:),2,'significant');
                % Call subfunction
                titleStr = createColoredTitle(mapName, meanProjection, phi, colorCon, colorIncon);
                % Set title with LaTeX interpreter
                title({titleStr,...
                    ['Columns: ' num2str(blockColumns(1)) ',' num2str(blockColumns(2))]}, 'FontWeight', 'normal', 'Interpreter', 'tex')
                xlabel('Visual contrast (%)');
                ylabel('Projection amplitude');
                axis square;
                upFontSize(20, .01);
            end
            %}
            
    
    
            %% STAGE 5b. Projection into visual-only space and reference map space
            figure('Name',['C' num2str(cluster) 'N' num2str(blockID) ': D' date 'R' run])
            nCols=2; nRows=1;
            colormap(fireice); 
            gap=0; marginV=.25; marginH=.2;
            [hAx,~]=tight_subplot(nRows,nCols,[gap+0.025 gap], [marginV marginV+.06], [marginH+.06 marginH]);
            mkrSize=15; lineWidth=3;
            % Loop over each map
            for map=2%1:2
                visualOrientation=NaN; 
                yZeroData=[];
                if map==1
                    visualOrientation=NaN; 
                    axes(hAx(map))
                    mapName='{Visual map}';
                    referenceImage=rescale(columnarResp(:,:,condIDs.V90(end)),0,1) - ...
                                  rescale(columnarResp(:,:,condIDs.V0(end)),0,1); 
                    [~, amplitude, ~, ~] = calculateProjectionVector(columnarResp, referenceImage, condIDs, visualOrientation, snrMask);
                    yMax = 50;yTickInterval=yMax/4;
                elseif map==2
                    axes(hAx(map))
                    mapName='{Reference map}';
                    referenceImage=rescale(PCAmap(:,:,7),0,1) - rescale(PCAmap(:,:,1),0,1);
                    [~, amplitude, ~, ~] = calculateProjectionVector(columnarResp, referenceImage, condIDs, visualOrientation, snrMask);
                    yMax = 50;yTickInterval=yMax/4;
                end
                
                % Store y-values for x=0 by marker type
                yCircleValues = [];
                yNonCircleValues = [];
    
                % Loop over each condition set
                for condSet=1:3                    
                    % Define conditions to plot
                    switch condSet
                        case {1} % BL V-only
                            conds=[condIDs.V0; condIDs.V90];
                            condMarker = 'o';
                            condMarkerColor = 'w';
                            condLineColor = 'k';
                            offsetBL = min(conds(:)) - 1;
                            condsOffset = conds - offsetBL;
                            gaborOrientations = TSbaseline(1).Header.Conditions.GaborOrt(TSbaseline(1).Header.Conditions.TypeCond == 3);
                            stimContrasts = TSbaseline(1).Header.Conditions.StimCon(TSbaseline(1).Header.Conditions.TypeCond == 3);
                            stimContrasts(gaborOrientations == 0) = -abs(stimContrasts(gaborOrientations == 0));
                            stimContrasts(gaborOrientations == 90) = abs(stimContrasts(gaborOrientations == 90));
                            xData = stimContrasts(condsOffset(2,:));
                            yData = amplitude(conds);
                            yData = mean(yData .* [-1;1]);
                            condsetStr = 'Visual-only';
    
                            % mergeZeros
                            zeroIdx=find(xData==0);
                            if sum(zeroIdx)>0
                                yData(zeroIdx)=mean(yData(zeroIdx));
                                xData(zeroIdx(1))=[];yData(zeroIdx(1))=[];
                            end
                        case {2} % Con V+O
                            conds=[condIDs.V0O0; condIDs.V90O90];
                            condMarker = '^';
                            colorRow = [0, 225, 80] / 255;
                            condMarkerColor = colorRow;
                            condLineColor = colorRow;
                            gaborOrientations = TSopto(1).Header.Conditions.GaborOrt(TSopto(1).Header.Conditions.TypeCond == 3);
                            stimContrasts = TSopto(1).Header.Conditions.StimCon(TSopto(1).Header.Conditions.TypeCond == 3);
                            stimContrasts(gaborOrientations == 0) = -abs(stimContrasts(gaborOrientations == 0));
                            stimContrasts(gaborOrientations == 90) = abs(stimContrasts(gaborOrientations == 90));
                            xData = stimContrasts(conds(2,:));
                            yData = amplitude(conds);
                            yData = mean(yData .* [-1;1]);
                            condsetStr = 'Visual + con-opto';

                            zeroConstrast=find(xData==0);
                            if ~isempty(zeroConstrast)
                                yZeroData=[yZeroData, yData(zeroConstrast)];
                            else
                                yZeroData=[yZeroData, NaN];
                            end
                        case {3} % Incon V+O
                            conds = [condIDs.V0O90; condIDs.V90O0];
                            condMarker = 'v';
                            colorRow = [156, 14, 254] / 255;
                            condMarkerColor = colorRow;
                            condLineColor = colorRow;
                            gaborOrientations = TSopto(1).Header.Conditions.GaborOrt(TSopto(1).Header.Conditions.TypeCond == 3);
                            stimContrasts = TSopto(1).Header.Conditions.StimCon(TSopto(1).Header.Conditions.TypeCond == 3);
                            stimContrasts(gaborOrientations == 0) = -abs(stimContrasts(gaborOrientations == 0));
                            stimContrasts(gaborOrientations == 90) = abs(stimContrasts(gaborOrientations == 90));
                            xData = -stimContrasts(conds(1,:));
                            yData = amplitude(conds);
                            yData = mean(yData .* [-1;1]);
                            condsetStr = 'Visual + incon-opto';

                            zeroConstrast=find(xData==0);
                            if ~isempty(zeroConstrast)
                                yZeroData=[yZeroData, yData(zeroConstrast)];
                            else
                                yZeroData=[yZeroData, NaN];
                            end
                    end
                    
                    % Plot the data
                    xlim([-100 100])
                    yline(0, '--', 'Color', [.4 .4 .4], 'LineWidth', 1.5, 'HandleVisibility', 'off');
                    plot(xData, yData, '-', 'LineWidth', lineWidth, ...
                        'Marker', condMarker, 'MarkerSize', mkrSize, 'MarkerFaceColor', condMarkerColor, ...
                        'MarkerEdgeColor', 'k', 'Color', condLineColor, 'DisplayName', condsetStr); 
                    hold on;

                    % Average
                    neuroStruct.(['C' num2str(cluster)]).refmap.averageProjection(blockNo,condSet)=mean(yData(:));
                    % Plot mean
                    meanLen  = 7;          % The half-width of the line segment
                    % Plot the line from (xVal - halfW, yVal) to (xVal + halfW, yVal)
                    plot([90, 90 + meanLen], [mean(yData(:)), mean(yData(:))], ...
                         'LineWidth', 4, 'Color', condMarkerColor, 'HandleVisibility','off');
                end
                
                % Compute neuroStruct.phi (difference between means of y-values at x=0)
                
                neuroStruct.(['C' num2str(cluster)]).refmap.phi(blockNo) = yZeroData(1)-yZeroData(2);
    
                % Set plot limits and labels
                %{
                meanCon=sprintf('%.1f', meanProjection(2));
                meanIncon=sprintf('%.1f', meanProjection(3));
                phi=sprintf('%.1f', round(neuroStruct.phi(map),2,'significant'));
                title({[mapName],['Mean = ' meanCon ', ' meanIncon,...
                    ', ɸ = ' phi]}, 'FontWeight', 'normal');
                %}
                phi=round(neuroStruct.(['C' num2str(cluster)]).refmap.phi(blockNo),2,'significant');
                meanProjection=round(neuroStruct.(['C' num2str(cluster)]).refmap.averageProjection(blockNo,:),2,'significant');

                % Call subfunction
                titleStr = createColoredTitle(mapName, meanProjection, phi, colorCon, colorIncon);
                
                % Set title with LaTeX interpreter
                title({titleStr,...
                    ['Columns: ' num2str(blockColumns(1)) ', ' num2str(blockColumns(2))]}, 'FontWeight', 'normal', 'Interpreter', 'tex')
                ylim([-100 100]); addSkippedTicks(-100,100,25, 'y');
                xlim([0 100]); addSkippedTicks(0, 100, 12.5, 'x');
                xlabel('Visual contrast (%)');
                ylabel('Similarity (%)');
                axis square;
                upFontSize(20, .01);


               %% Add inset of subtracted image
                % Choose the area analyzed
                gcf;
                if map==1
                    axInset = axes('Position',[.29+.043 .257 .18 .18]);

                    %axInset = axes('Position',[.3 .52 .13 .13]);
                elseif map==2
                    axInset = axes('Position',[.56+.043 .257 .18 .18]);

                    %axInset = axes('Position',[.565 .52 .13 .13]);
                end
                imgsc(columnarResp(:,:,condIDs.V90O90(1)).*snrMask - columnarResp(:,:,condIDs.V0O0(1)).*snrMask); hold on; 
                if strcmp(optostimAreaFlag,'nonstim')
                    if mean(blockColumns)>=15
                        cax(.8);
                    elseif mean(blockColumns)>=3 & mean(blockColumns)<15
                        cax(.6);
                    else
                        cax(.3);
                    end
                else
                    cax();
                end
                box off; set(axInset, 'XTick', [], 'YTick', [], 'Box', 'off', 'XColor', 'none', 'YColor', 'none','Color','none'); % No ticks, no box

            end

            % Colorbar struff
            hFig = gcf();
            allKids = hFig.Children;
            for k = 1:numel(allKids)
                % A ColorBar in MATLAB is a 'matlab.graphics.illustration.ColorBar' object
                if isa(allKids(k), 'matlab.graphics.illustration.ColorBar')
                    allKids(k).FontSize = 10;
                    %allKids(k).Color='none';
                end
            end

            %% Saving
            if saveFlag
                % Filenames
                pngPath=[mainPath monkey '\Meta\neurometric\projections\png\'];                
                svgPath=[mainPath monkey '\Meta\neurometric\projections\svg\'];
                figPath=[mainPath monkey '\Meta\neurometric\projections\fig\'];
                figFilename=['C' num2str(cluster) 'M28D' date 'R' run optostimAreaFlag 'slides'];

                % Set font and vector renderer
                set(findall(gcf, '-property', 'FontName'), 'FontName', 'SansSerif');                
                set(gcf, 'Renderer', 'painters'); % Use painters for vector graphics

                % Adjust figure to remove white space
                exportgraphics(gcf, [pngPath figFilename '.png'], 'Resolution', 600); % High-res PNG %print(gcf, [figFilename '.png'], '-dpng', '-r600'); 
                %savefig(gcf, [figPath figFilename '.fig']);           % FIG
                %savefig(gcf, [svgPath figFilename '.svg']);           % FIG
                %exportgraphics(gcf, [svgPath figFilename '.svg'],'ContentType','vector');        % SVG
            end
            %{
            h=title({titleStr,...
                             ['${\neuroStruct.phi}$: ' num2str(neuroStruct.phi(2),'%.2f') ', '  num2str(neuroStruct.phi(1),'%.2f') ', ' num2str(neuroStruct.phi(3),'%.2f')]},...
                             'FontWeight','normal','FontSize',14,'FontName','Arial','interpreter','latex');
            %}
            
            toc
        end   
    else
        disp('Cluster empty')
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
