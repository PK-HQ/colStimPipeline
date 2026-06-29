function [mdl, reportState]=plotNakaRushtonFit5(behavioralData, bitmapData, datastruct, analysisBlockID,...
    mdl, fitParams, x, monkeyName, clusterBlocks, plotAverageFlag, plotLine,...
    saveFlag, cluster, modelTypeStr, savefilename, clusterLabel, reportState, plotOpts)
    if nargin < 16
        clusterLabel = '';
    end
    if nargin < 17
        reportState = [];
    end
    if nargin < 18 || isempty(plotOpts)
        plotOpts = struct();
    end
    plotOpts = applyDeltaPermutationPlotDefaults(plotOpts);

    % Reuse the settled whole-experiment delta-bias permutation test for
    % individual weibullFreeAll pages. The test uses the empirical merged
    % Con-Opto/Incon-Opto points and is independent of fitted model parameters.
    if strcmpi(modelTypeStr, 'weibullFreeAll') && ~plotAverageFlag
        plotOpts.showDeltaPermutationStats = true;
    end
    if plotOpts.showDeltaPermutationStats
        fprintf('plotNakaRushtonFit5 received showDeltaPermutationStats=true (%d permutations).\n', ...
            plotOpts.nDeltaPermutations);
    end
    endIdx=size(mdl.headers,2);
    % Number of blocks and conditions
    nConditions = size(behavioralData.gaborContrasts, 1);
    nBlocks = numel(clusterBlocks);
    if plotAverageFlag==1
        nBlocks=1;
    end
    % Normalize cluster labels for plotting/saving.  A plotting call already
    % represents one numeric cluster, so when the caller does not supply
    % labels (or supplies one shared label), repeat that label for every
    % rendered block.  This keeps saved pages annotated without forcing the
    % caller to construct a redundant per-block label array.
    clusterLabel = normalizeClusterLabels(clusterLabel, cluster, nBlocks);
    fprintf('plotNakaRushtonFit5 render setup: nRenderBlocks=%d | numel(clusterLabel)=%d | plotAverageFlag=%d\n', ...
        nBlocks, numel(clusterLabel), logical(plotAverageFlag));
    mdl.cluster=cluster;
    mdl.clusterBlocksIdx=clusterBlocks;
    
    % Dynamic x-axis limits for this animal/chamber/cluster set.
    % These are computed from all currently available mdl rows, so every block
    % in this plotting call uses the same x-limits and same mean-bar geometry.
    [xLimPre, xLimMerged, meanBar, tickCfg] = getPsychometricAxisLimits(mdl);
    xLimMerged = [0 100];
    tickCfg.mergedMajor = 100 / 8;
    tickCfg.mergedSkip = tickCfg.mergedMajor;
    meanBar.rightEdgeRange = [95 100];

    for block = 1:nBlocks
        blockInfo = getPlotBlockInfo(datastruct, analysisBlockID, clusterBlocks, block, plotAverageFlag);
        pageClusterLabel = getClusterLabelForRenderBlock(clusterLabel, block);
        appendPage = block > 1;
        fprintf('plotNakaRushtonFit5 render page: renderIdx=%d | globalBlockIdx=%d | clusterLabel=%s | appendPage=%d\n', ...
            block, blockInfo.blockIdx, char(pageClusterLabel), appendPage);
        [baselineModeThis, baselineSeparateThis, baselineSourceThis] = ...
            detectBaselineModeByBlock(datastruct, analysisBlockID, blockInfo.blockIdx);
        baselineModeThis = baselineModeThis(1);
        baselineSourceThis = baselineSourceThis(1);
        mdl.baselineMode(block, 1) = baselineModeThis;
        mdl.baselineModeSourceField = 'datastruct(analysisBlockID(blockIdx)).baselineTS';
        mdl.baselineTSValue(block, 1) = baselineSourceThis;
        mdl.combinedBL(block) = ~baselineSeparateThis(1);
        % Init figure
        dat=[];
        make_it_tight = true;
        hmarg = .18;
        wmarg = [0.12 0.2];
        panelGap = [0.09 0.055];
        subplot = @(m,n,p) subtightplot(m, n, p, panelGap, [hmarg hmarg], wmarg);
        if ~make_it_tight,  clear subplot;  end
       
        % Create a real standalone figure. Saving overrides any upstream
        % figureVisible='off' setting because hidden/docked figures in R2018b
        % remain at MATLAB's default 560x420 canvas and do not export WYSIWYG.
        renderVisible = plotOpts.figureVisible;
        if saveFlag
            renderVisible = 'on';
        end

        oldRootUnits = get(groot, 'Units');
        set(groot, 'Units', 'pixels');
        screenRect = get(groot, 'ScreenSize');
        set(groot, 'Units', oldRootUnits);

        targetWidth = min(1800, max(1200, screenRect(3) - 100));
        targetHeight = min(1050, max(750, screenRect(4) - 150));
        targetLeft = max(1, round(screenRect(1) + ...
            (screenRect(3) - targetWidth) / 2));
        targetBottom = max(1, round(screenRect(2) + ...
            (screenRect(4) - targetHeight) / 2));
        targetFigurePosition = [targetLeft targetBottom targetWidth targetHeight];

        fig = figure('Name', ['Block ', blockInfo.label], ...
            'WindowStyle', 'normal', ...
            'Visible', renderVisible, ...
            'Color', 'w', ...
            'Units', 'pixels', ...
            'Position', targetFigurePosition, ...
            'PaperPositionMode', 'auto');

        % Reapply after creation because a docked default can otherwise win.
        set(fig, 'WindowStyle', 'normal', ...
            'Units', 'pixels', ...
            'Position', targetFigurePosition, ...
            'Visible', renderVisible);
        drawnow;
        
        sideData = getPreMergedSideData(mdl, block);
        signedBX0DisplayParams = [];
        if strcmp(modelTypeStr, 'weibullSignedBX0')
            signedBX0DisplayParams = getSignedBX0DisplayParamsForBlock(...
                fitParams(block, :));
        end

        % Row 1, columns 1-2: split pre-merged data by visual stimulus side
        axRow1Col1 = subplot(2,3,1);
        [mdl, sideData.horizontal] = plotSidePsychometricPanel(mdl, block, sideData.horizontal, ...
            xLimMerged, meanBar, tickCfg, 'Horizontal visual stimulus');

        axRow1Col2 = subplot(2,3,2);
        [mdl, sideData.vertical] = plotSidePsychometricPanel(mdl, block, sideData.vertical, ...
            xLimMerged, meanBar, tickCfg, 'Vertical visual stimulus');

        signedBX0PanelFits = [];
        bx0HorizontalDeltaFit = [];
        bx0VerticalDeltaFit = [];
        if strcmp(modelTypeStr, 'weibullSignedBX0')
            globalDeltaX0 = getSignedBX0GlobalDeltaX0(mdl, block, ...
                signedBX0DisplayParams);
            [mdl, signedBX0PanelFits] = fitSignedBX0PanelFitsForBlock( ...
                mdl, block, sideData, globalDeltaX0, ...
                signedBX0DisplayParams, blockInfo);

            bx0XGrid = linspace(xLimMerged(1), xLimMerged(2), 400);
            displayCurves.horizontal = predictSignedBX0PanelFitCurves( ...
                bx0XGrid, signedBX0PanelFits.horizontal.fitParams, ...
                signedBX0PanelFits.globalDeltaX0, 'horizontal');
            displayCurves.vertical = predictSignedBX0PanelFitCurves( ...
                bx0XGrid, signedBX0PanelFits.vertical.fitParams, ...
                signedBX0PanelFits.globalDeltaX0, 'vertical');

            assertSignedBX0PanelCurveSource(displayCurves.horizontal, 'horizontal');
            assertSignedBX0PanelCurveSource(displayCurves.vertical, 'vertical');
            overlaySignedBX0SideCurves(axRow1Col1, bx0XGrid, ...
                displayCurves.horizontal);
            overlaySignedBX0SideCurves(axRow1Col2, bx0XGrid, ...
                displayCurves.vertical);

            bx0HorizontalDeltaFit = struct('x', bx0XGrid, ...
                'viewCurves', displayCurves.horizontal, ...
                'sourcePanel', 'horizontal');
            bx0VerticalDeltaFit = struct('x', bx0XGrid, ...
                'viewCurves', displayCurves.vertical, ...
                'sourcePanel', 'vertical');
            validateSignedBX0PanelFitFields(mdl.signedBX0.panelFits.horizontal, ...
                mdl.signedBX0.panelFits.vertical, ...
                mdl.signedBX0.panelFits.merged, block);
        end

        % Row 2, columns 1-2: deltas from the side-specific split data
        subplot(2,3,4)
        mdl = plotSideDeltaPanel(mdl, block, sideData.horizontal, ...
            xLimMerged, meanBar, tickCfg, 'Horizontal visual stimulus', ...
            'Horizontal', blockInfo, baselineModeThis, bx0HorizontalDeltaFit);

        subplot(2,3,5)
        mdl = plotSideDeltaPanel(mdl, block, sideData.vertical, ...
            xLimMerged, meanBar, tickCfg, 'Vertical visual stimulus', ...
            'Vertical', blockInfo, baselineModeThis, bx0VerticalDeltaFit);
        
        % Row 1, column 3: merged fitted data
        axMergedPanel = subplot(2,3,3);
        axMergedDeltaPanel = gobjects(1);
        optoStatsTextHandle = gobjects(1);
        hold on;
        yline(50,'--','LineWidth',1.5,'Color',.4*[1 1 1],'HandleVisibility','off'); hold on;

        for cond = 1:nConditions+2
            % Extract fitted parameters for current condition and block
            beta = fitParams(block, 1);
            n = fitParams(block, 2);
            C50 = fitParams(block, 3);
            if cond==1
                deltax = 0;
            elseif cond==2
                deltax = fitParams(block, end-2);
            elseif cond==3
                deltax = fitParams(block, end-1);
            end
            aicc = fitParams(block, end);
            
            % Define the Naka-Rushton functions with the fitted parameters
            betaBL=50;
            switch cond
                case 1 % Baseline
                    xPlot = x;
                    if strcmp(modelTypeStr, 'weibullSignedBX0')
                        displayCurves = predictSignedBX0PanelFitCurves(...
                            xPlot, signedBX0PanelFits.merged.fitParams, ...
                            signedBX0PanelFits.globalDeltaX0, 'merged');
                        assertSignedBX0PanelCurveSource(displayCurves, 'merged');
                        predictedCurve = displayCurves.baseline;
                    elseif strcmp(modelTypeStr, 'weibullSignedX0')
                        displayCurves = predictSignedX0DisplayCurves(...
                            xPlot, fitParams(block, :));
                        predictedCurve = displayCurves.merged.baseline;
                    else
                        predictedCurve = mdl.mdlBaseline(xPlot, fitParams(block, 1:end-1)).pcntrl;
                    end
                    lineColor = [0 0 0]; % Black for baseline
                    markerFaceColor = [1 1 1]; 
                    edgeColor = 'k';
                    markerType = 'o';
                    xBlock = rmnan(mdl.xBaseline(block, :));
                    yBlock = rmnan(mdl.yBaseline(block, :));
                    [~, idx] = find(xBlock>=0);
                    xBlock = xBlock(idx);
                    yBlock = yBlock(idx);
                    predictedCurve = replaceNanSections(predictedCurve);
                    mdl.xFitted(cond,:,block)=padArray(xPlot,400,2,nan);
                    mdl.yFitted(cond,:,block)=padArray(predictedCurve,400,2,nan);
                    mdl.xBlock(cond,:,block)=padArray(xBlock, 120, 2, nan);
                    mdl.yBlock(cond,:,block)=padArray(yBlock, 120, 2, nan);
                    %get AUC
                    idxx=mdl.xFitted(cond,:,block)>=0 & ~isnan(mdl.yFitted(cond,:,block));
                    contr=mdl.xFitted(cond,idxx,block);
                    curve=mdl.yFitted(cond,idxx,block);
                    %get contrast
                    mdl.thresholdContrast(block, 1) = getThreshold(mdl.xFitted(cond,:,block), mdl.yFitted(cond,:,block), 70);

                    %mdl.fittedParams(block,endIdx) = trapz(contr, curve) / (max(contr) - min(contr));
                case 2 % Con-Opto
                    xPlot = x; % Positive contrasts for congruent condition
                    if strcmp(modelTypeStr, 'weibullSignedBX0')
                        displayCurves = predictSignedBX0PanelFitCurves(...
                            xPlot, signedBX0PanelFits.merged.fitParams, ...
                            signedBX0PanelFits.globalDeltaX0, 'merged');
                        assertSignedBX0PanelCurveSource(displayCurves, 'merged');
                        predictedCurve = displayCurves.con;
                    elseif strcmp(modelTypeStr, 'weibullSignedX0')
                        displayCurves = predictSignedX0DisplayCurves(...
                            xPlot, fitParams(block, :));
                        predictedCurve = displayCurves.merged.con;
                    else
                        predictedCurve = mdl.mdlOpto(xPlot, fitParams(block, 1:end-1)).pc;
                    end
                    xPlot = x(x >= 0); % Positive contrasts for congruent condition
                    predictedCurve=predictedCurve(x >= 0);
                    lineColor = [0.9294, 0.1098, 0.1373] * 1.05; % Red for con-opto
                    markerFaceColor = lineColor;
                    edgeColor = 'k';
                    markerType = '^';
                    xBlock = rmnan(mdl.xConOpto(block, :));
                    yBlock = rmnan(mdl.yConOpto(block, :));
                    [~, idx] = find(xBlock >= 0);
                    xBlock = xBlock(idx);
                    yBlock = yBlock(idx);
                    predictedCurve = replaceNanSections(predictedCurve);
                    mdl.xFitted(cond,:,block)=padArray(xPlot,400,2,nan);
                    mdl.yFitted(cond,:,block)=padArray(predictedCurve,400,2,nan);
                    mdl.xBlock(cond,:,block)=padArray(xBlock, 120, 2, nan);
                    mdl.yBlock(cond,:,block)=padArray(yBlock, 120, 2, nan);
                    %get AUC
                    idxx=mdl.xFitted(cond,:,block)>=0 & ~isnan(mdl.yFitted(cond,:,block));
                    contr=mdl.xFitted(cond,idxx,block);
                    curve=mdl.yFitted(cond,idxx,block);
                    mdl.thresholdContrast(block, 2) = getThreshold(mdl.xFitted(cond,:,block), mdl.yFitted(cond,:,block), 70);

                    %mdl.fittedParams(block,endIdx+1) = trapz(contr, curve) / (max(contr) - min(contr));
                case 3 % Incon-Opto
                    xPlot = x; % Positive contrasts for congruent condition
                    if strcmp(modelTypeStr, 'weibullSignedBX0')
                        displayCurves = predictSignedBX0PanelFitCurves(...
                            xPlot, signedBX0PanelFits.merged.fitParams, ...
                            signedBX0PanelFits.globalDeltaX0, 'merged');
                        assertSignedBX0PanelCurveSource(displayCurves, 'merged');
                        predictedCurve = displayCurves.incon;
                    elseif strcmp(modelTypeStr, 'weibullSignedX0')
                        displayCurves = predictSignedX0DisplayCurves(...
                            xPlot, fitParams(block, :));
                        predictedCurve = displayCurves.merged.incon;
                    else
                        predictedCurve = mdl.mdlOpto(xPlot, fitParams(block, 1:end-1)).pic;
                    end

                    lineColor = [0, 0.0941, 0.6627] * 1.25; % Blue for incon-opto
                    markerFaceColor = lineColor;
                    edgeColor = 'k';
                    markerType = 'v';
                    xBlock = rmnan(mdl.xInconOpto(block, :));
                    yBlock = rmnan(mdl.yInconOpto(block, :));
                    xBlock = 1.*xBlock(idx);
                    yBlock = yBlock(idx);
                    predictedCurve = replaceNanSections(predictedCurve);
                    mdl.xFitted(cond,:,block)=padArray(xPlot,400,2,nan);
                    mdl.yFitted(cond,:,block)=padArray(predictedCurve,400,2,nan);
                    mdl.xBlock(cond,:,block)=padArray(xBlock, 120, 2, nan);
                    mdl.yBlock(cond,:,block)=padArray(yBlock, 120, 2, nan);
                    mdl.thresholdContrast(block, 3) = getThreshold(mdl.xFitted(cond,:,block), mdl.yFitted(cond,:,block), 70);

                case 4
                    idxPos=mdl.xFitted(2,:, block)>=0;
                    idxNeg=mdl.xFitted(3,:, block)>=0;

                    % Biasing datapoints: Con - incon datapoints
                    xBlock = mdl.xBlock(2,:, block);
                    yBlock = rmnan(mdl.yBlock(2,:, block))-rmnan((mdl.yBlock(3,:, block)));
                    % Biasing curve: Con - incon curve
                    xPlot = mdl.xFitted(2,idxPos, block);
                    predictedCurve = mdl.yFitted(2,idxPos, block)-mdl.yFitted(3,idxNeg, block);

                    markerType = 'square';
                    ylimMax=50;%roundup(max(abs(yBlock(:))),10);

                    lineColor =[127, 0, 255]/255; % Gray for combined case
                    markerFaceColor =lineColor;
                    edgeColor = 'k';
                    %TEMPORARY FIX
                    predictedCurve = replaceNanSections(predictedCurve);
                    mdl.xFitted(cond,:,block)=padArray([nan(size(xPlot,1),size(xPlot,2)), xPlot],400,2,nan);
                    mdl.yFitted(cond,:,block)=padArray([nan(size(predictedCurve,1),size(predictedCurve,2)), predictedCurve],400,2,nan);
                    mdl.xBlock(cond,:,block)=padArray(xBlock, 120, 2, nan);
                    mdl.yBlock(cond,:,block)=padArray(yBlock, 120, 2, nan);
                    %get AUC
                    idxx=mdl.xFitted(cond,:,block)>=0 & ~isnan(mdl.yFitted(cond,:,block));
                    contr=mdl.xFitted(cond,idxx,block);
                    curve=mdl.yFitted(cond,idxx,block);
                    %mdl.fittedParams(block,endIdx+3) = trapz(contr, curve) / (max(contr) - min(contr));
                    
                    %mdl.fittedParams(block,endIdx+2) = mdl.fittedParams(block,end-3) - mdl.fittedParams(block,end);

                case 5 % masking line
                    idxPos=mdl.xFitted(2,:, block)>=0;
                    idxNeg=mdl.xFitted(3,:, block)>=0;

                    % Biasing datapoints: Con - incon datapoints
                    if length(rmnan(mdl.yBlock(1,:, block))) == length(rmnan(mdl.yBlock(2,:, block)))
                        xBlock = mdl.xBlock(2,:, block);
                        yBlock = rmnan(mdl.yBlock(1,:, block))-(rmnan(mdl.yBlock(2,:, block))+rmnan((mdl.yBlock(3,:, block))))/2;
                        % Biasing curve: Con - incon curve
                        xPlot = mdl.xFitted(2,idxPos, block);
                        predictedCurve = mdl.yFitted(1,idxPos, block)-(mdl.yFitted(2,idxPos, block)+mdl.yFitted(3,idxNeg, block))/2;
                    else
                        optoX=rmnan(mdl.xBlock(2,:, block));
                        baseX=rmnan(mdl.xBlock(1,:, block));
                        baseY=rmnan(mdl.yBlock(1,:, block));
                        [baseXpad, baseYpad, missingIdxInOpto] = padMissingX(baseX, baseY, optoX);
                        xBlock = baseXpad;
                        yBlock = baseYpad-(rmnan(mdl.yBlock(2,:, block))+rmnan((mdl.yBlock(3,:, block))))/2;
                        % Biasing curve: Con - incon curve
                        xPlot = mdl.xFitted(2,idxPos, block);
                        predictedCurve = mdl.yFitted(1,idxPos, block)-(mdl.yFitted(2,idxPos, block)+mdl.yFitted(3,idxNeg, block))/2;
                        mdl.xBlock(1,:, block);
                    end

                    markerType = 'square';
                    ylimMax=50;%roundup(max(abs(yBlock(:))),10);

                    lineColor =[125, 125, 125]/255; % Gray for combined case
                    markerFaceColor =lineColor;
                    edgeColor = 'k';
                    %TEMPORARY FIX
                    predictedCurve = replaceNanSections(predictedCurve);
                    mdl.xFitted(cond,:,block)=padArray([nan(size(xPlot,1),size(xPlot,2)), xPlot],400,2,nan);
                    mdl.yFitted(cond,:,block)=padArray([nan(size(predictedCurve,1),size(predictedCurve,2)), predictedCurve],400,2,nan);
                    mdl.xBlock(cond,:,block)=padArray(xBlock, 120, 2, nan);
                    mdl.yBlock(cond,:,block)=padArray(yBlock, 120, 2, nan);
                    %get AUC
                    idxx=mdl.xFitted(cond,:,block)>=0 & ~isnan(mdl.yFitted(cond,:,block));
                    contr=mdl.xFitted(cond,idxx,block);
                    curve=mdl.yFitted(cond,idxx,block);
                    %mdl.fittedParams(block,endIdx+3) = trapz(contr, curve) / (max(contr) - min(contr));
                    
                    %mdl.fittedParams(block,endIdx+2) = mdl.fittedParams(block,end-3) - mdl.fittedParams(block,end);

                otherwise
                    lineColor = 'g'; % Fallback lineColor
                    markerFaceColor = lineColor;
                    edgeColor = 'g';
                    markerType='.';
            end


            %% Plots
            if cond<=3
                % Plot line fit
                if plotLine==1
                    plot(mdl.xFitted(cond,:,block), mdl.yFitted(cond,:,block), 'Color', lineColor, 'LineWidth', 3, 'HandleVisibility', 'on'); hold on;
                end
                upFontSize(32, 0.01); legend; axis square

                % Plot scatterplot
                semY=std(mdl.yBlock(cond,:,block), 'omitnan') / sqrt(length(yBlock));
                counts=size(mdl.yBlock(cond,:,block),2);
                if size(mdl.yBlock(cond,:,block),1) > 1 % if no sem, don't shade
                    patchSaturationVal = 0.2;
                else
                    semY=zeros(counts,1);
                    patchSaturationVal = 0;
                end

                % Add data points and shaded error bar
                markerSize = 12;
                patchSaturationVal=1;
                % Data points
                shadedErrorBar(mdl.xBlock(cond,:,block), mdl.yBlock(cond,:,block), semY', 'patchSaturation', patchSaturationVal, 'lineprops', ...
                               {'Color', lineColor, 'LineStyle', 'none', 'LineWidth', 3, 'Marker', markerType, ...
                                'MarkerFaceColor', markerFaceColor, 'MarkerEdgeColor', edgeColor, 'MarkerSize', markerSize}); hold on;
                % Average
                barLength = meanBar.rightEdgeRange;
                meanMergedVal = weightedMeanForPlot(mdl.xBlock(cond,:,block), mdl.yBlock(cond,:,block), 'merged');
                plot(barLength, repmat(meanMergedVal,1,numel(barLength)), '-', 'Color', lineColor, 'LineWidth', 3, 'HandleVisibility', 'off')

                switch cond
                    case 1
                        mdl.meanBaselineMerged(block) = meanMergedVal;
                    case 2
                        mdl.meanConOptoMerged(block) = meanMergedVal;
                    case 3
                        mdl.meanInconOptoMerged(block) = meanMergedVal;
                        mdl.meanPsychometricMerged(block,:) = [ ...
                            mdl.meanBaselineMerged(block), ...
                            mdl.meanConOptoMerged(block), ...
                            mdl.meanInconOptoMerged(block)];
                        mdl.meanPsychometricHeaders = {'Baseline', 'ConOpto', 'InconOpto'};
                end

                % Add datapoints's count annotation
                zeroConstrastPoint=xBlock==0;
                if sum(zeroConstrastPoint)>0 % data contains 0-point
                    nTrials=[20*ones(1,numel(zeroConstrastPoint)) 20*ones(1,numel(~zeroConstrastPoint))];
                    if cond>1
                            nTrials=[40*zeroConstrastPoint + 20*~zeroConstrastPoint];
                    end
                elseif sum(zeroConstrastPoint)>0 && plotAverageFlag==1 % data contains 0-point & average
                    nTrials=[nBlocks*ones(1,numel(zeroConstrastPoint)) nBlocks*ones(1,numel(~zeroConstrastPoint))];
                    if cond>1
                            nTrials=[40*zeroConstrastPoint + 20*~zeroConstrastPoint];
                    end
                elseif sum(zeroConstrastPoint)==0
                     nTrials=[20*ones(1,numel(~zeroConstrastPoint))];
                elseif sum(zeroConstrastPoint)==0 && plotAverageFlag==1
                     nTrials=[nBlocks*ones(1,numel(~zeroConstrastPoint))];
                end
                %annotateDataPoints(mdl.xBlock(cond,:,block), mdl.yBlock(cond,:,block), nTrials, markerFaceColor); hold on;
            
                
                if cond==3
                    combinedBLStr = char(baselineModeThis);

                    if isfield(bitmapData, 'meanPowerDensityWithinROI_mWmm2') & ~isempty(bitmapData.meanPowerDensityWithinROI_mWmm2(:,:,blockInfo.blockIdx))
                        if plotAverageFlag % for plotting the average across all blocks
                            bitmapSPD=squeeze(bitmapData.meanPowerDensityWithinROI_mWmm2(:,:,clusterBlocks));
                            bitmapColumns=bitmapData.nColumns(:,clusterBlocks)';
                        else
                            bitmapSPD=squeeze(bitmapData.meanPowerDensityWithinROI_mWmm2(:,:,blockInfo.blockIdx));
                            bitmapColumns=bitmapData.nColumns(:,blockInfo.blockIdx)';
                        end
                    else
                        bitmapSPD=[0 0];
                    end
                    bitmapColumnhv=nanmean(bitmapColumns,1);
                    bitmapSPDhv=nanmean(bitmapSPD,1);
                    bitmapSPDmean=nanmean(bitmapSPD,'all');
                    bitmapSPDstd=nanstd(bitmapSPD,[],'all');
                    if plotAverageFlag==1
                        bitmapSPDhv1 = formatPowerMetricsForDisplay( ...
                            bitmapSPDhv(1), NaN);
                        bitmapSPDhv2 = formatPowerMetricsForDisplay( ...
                            bitmapSPDhv(2), NaN);
                        title({[modelTypeStr ', cluster ' num2str(cluster) ' average'],...
                            ['meanPowerDensityWithinROI_mWmm2: ' ...
                            bitmapSPDhv1.PDROI ' & ' bitmapSPDhv2.PDROI ...
                            ' mW mm^{-2} (' sprintf('%.4f', bitmapSPDmean) ...
                            ' \pm ' sprintf('%.4f', bitmapSPDstd) ...
                            ' mW mm^{-2})',...
                            ', Columns: ' num2str(bitmapColumnhv(1),2) ' & ' num2str(bitmapColumnhv(2),2)]});
                    else
                        %{
                        title({[blockInfo.label ' (' combinedBLStr ')'],...
                            ['meanPowerDensityWithinROI_mWmm2: ' num2str(bitmapSPDhv(1),2) ' & ' num2str(bitmapSPDhv(2),2) ' mW, ',...
                            'Columns: ' num2str(bitmapColumnhv(1),2) ' & ' num2str(bitmapColumnhv(2),2)]});
                        %}
                    end
                    % Labels etc
                    %axis square
                    xlim(xLimMerged); ylim([0 100]);
                    xticks(xLimMerged(1):tickCfg.mergedMajor:xLimMerged(2));
                    addSkippedTicks(xLimMerged(1), xLimMerged(2), tickCfg.mergedSkip, 'x');
                    addSkippedTicks(0, 100, 10, 'y');
                    axis square
                    % Adding legend after plotting to ensure it covers all conditions
                    moveLines()
                    h2 = get(gca,'Children');
                    legend({'Baseline', 'Con-Opto', 'Incon-Opto'}, 'Location', 'southeast',...
                       'NumColumns',1,'FontSize',32);
                    %{
                    legend(h2([end-2:end]), {'Baseline', 'Con-Opto', 'Incon-Opto'}, 'Location', 'east',...
                       'NumColumns',1,'FontSize',32);
                    %}                    
                    upFontSize(32, 0.01)
            
                    
                    % Add text for biasing
                    xOffset=.53;
                    yOffset=.05;
                    %text('Units', 'Normalized', 'Position', [1 1]-[xOffset yOffset], 'string', 'More biasing', 'color', 'k','FontWeight','bold', 'Fontsize',14)
                    ax = gca;
                    ylabel('Correct (%)'); set(gca,'ycolor','k') 
                    xlabel('Gabor contrast (%)');
                    % --- Optostim metadata annotation ---
                    clusterNo = cluster;
                    
                    blockIdx = blockInfo.blockIdx;
                    
                    powerMetrics = computePowerMetricsFromSource(bitmapData, blockIdx);
                    powerSummary = powerMetrics.summary;
                    if ~powerMetrics.pass
                        warning('plotNakaRushtonFit5:PowerAuditMismatch', ...
                            ['%s: stored power fields differ from canonical ' ...
                            'PDDMD*area*duty-cycle recomputation. Using ' ...
                            'recomputed values for plot text.'], blockInfo.label);
                    end
                    meanColumns = powerSummary.columns;
                    meanProjectorPD = powerSummary.projectorPowerDensity;
                    meanAreaROI = powerSummary.areaROI;
                    meanAreaON = powerSummary.areaON;
                    meanSpatialDC = powerSummary.spatialDutyCycleFraction * 100;
                    meanTemporalDC = powerSummary.temporalDutyCycleFraction * 100;
                    meanROIPD = powerSummary.roiPowerDensityRecomputed;
                    meanTotalPower = powerSummary.totalPowerRecomputed;
                    powerDisplay = formatPowerMetricsForDisplay( ...
                        meanROIPD, meanTotalPower);
                    
                    title('Merged fitted', 'Interpreter', 'none');

                    clusterLabelLine = formatClusterLabelLine(pageClusterLabel);
                    optoText = sprintf([ ...
                        '%s' ...
                        'BL: %s\n' ...
                        '%.0f cols\n' ...
                        'PD_{DMD} %.2f mW mm^{-2}\n' ...
                        'Area_{ROI} %.2f mm^2\n' ...
                        'Area_{ON} %.2f mm^2\n' ...
                        'sDC %.1f%% | tDC %.1f%%\n' ...
                        'PD_{ROI} %s mW mm^{-2}\n' ...
                        'P_{total} %s mW'], ...
                        clusterLabelLine, ...
                        char(baselineModeThis), ...
                        meanColumns, ...
                        meanProjectorPD, ...
                        meanAreaROI, ...
                        meanAreaON, ...
                        meanSpatialDC, meanTemporalDC, ...
                        powerDisplay.PDROI, ...
                        powerDisplay.Ptotal);

                    if strcmp(modelTypeStr, 'weibullSignedX0')
                        optoText = sprintf('%s\n%s', optoText, ...
                            signedX0AnnotationText(mdl, block));
                    end
                    optoStatsTextHandle = addOptoStatsText(axRow1Col1, optoText);

                    if strcmp(modelTypeStr, 'weibullSignedBX0')
                        addSignedBX0FitParameterTable(axMergedPanel, ...
                            signedBX0PanelFits.merged.fitParams, ...
                            signedBX0PanelFits.globalDeltaX0, ...
                            mdl.signedBX0.deltaAICcX0(block));
                    else
                        addFitParameterTable(axMergedPanel, mdl.headers, fitParams(block,:), modelTypeStr);
                    end
                    
                    1;
                end
            elseif cond==4
                subplot(2,3,6)
                axMergedDeltaPanel = gca;
                % Plot line fit
                if plotLine==1
                    plot(mdl.xFitted(cond,:,block), mdl.yFitted(cond,:,block), 'Color', lineColor, 'LineWidth', 3, 'HandleVisibility', 'on'); hold on;
                end
                upFontSize(32, 0.01); legend; axis square

                % Plot scatterplot
                semY=std(mdl.yBlock(cond,:,block), 'omitnan') / sqrt(length(yBlock));
                counts=size(mdl.yBlock(cond,:,block),2);
                if size(mdl.yBlock(cond,:,block),1) > 1 % if no sem, don't shade
                    patchSaturationVal = 0.2;
                else
                    semY=zeros(counts,1);
                    patchSaturationVal = 0;
                end

                % Add data points and shaded error bar
                markerSize = 15;
                patchSaturationVal=0.14;
                % Data points
                shadedErrorBar(mdl.xBlock(cond,:,block), mdl.yBlock(cond,:,block), semY', 'patchSaturation', patchSaturationVal, 'lineprops', ...
                               {'Color', lineColor, 'LineStyle', 'none', 'LineWidth', 3, 'Marker', markerType, ...
                                'MarkerFaceColor', markerFaceColor, 'MarkerEdgeColor', edgeColor, 'MarkerSize', markerSize}); hold on;
                % Average
                [baseMeanMerged, conMeanMerged, inconMeanMerged] = computeMergedConditionMeans(mdl, block);
                [deltaBiasMergedForPlot, ~] = computeDeltaFromConditionMeans(baseMeanMerged, conMeanMerged, inconMeanMerged);
                plot(barLength, repmat(deltaBiasMergedForPlot,1,numel(barLength)), '-', 'Color', lineColor, 'LineWidth', 3, 'HandleVisibility', 'off')

                % Add datapoints's count annotation
                zeroConstrastPoint=xBlock==0;
                if sum(zeroConstrastPoint)>0
                    nTrials=[20*ones(1,numel(zeroConstrastPoint)) 20*ones(1,numel(~zeroConstrastPoint))];
                    if cond>1
                            nTrials=[40*zeroConstrastPoint + 20*~zeroConstrastPoint];
                    end
                elseif sum(zeroConstrastPoint)>0 && plotAverageFlag==1
                    nTrials=[nBlocks*ones(1,numel(zeroConstrastPoint)) nBlocks*ones(1,numel(~zeroConstrastPoint))];
                    if cond>1
                            nTrials=[40*zeroConstrastPoint + 20*~zeroConstrastPoint];
                    end
                elseif sum(zeroConstrastPoint)==0
                     nTrials=[20*ones(1,numel(~zeroConstrastPoint))];
                elseif sum(zeroConstrastPoint)==0 && plotAverageFlag==1
                     nTrials=[nBlocks*ones(1,numel(~zeroConstrastPoint))];
                end
                upFontSize(32, 0.01);
                %annotateDataPoints(mdl.xBlock(cond,:,block), mdl.yBlock(cond,:,block), nTrials, markerFaceColor); hold on;
                xlim(xLimMerged)
                ylim([-75 75])
                xticks(xLimMerged(1):tickCfg.mergedMajor:xLimMerged(2));
                addSkippedTicks(xLimMerged(1), xLimMerged(2), tickCfg.mergedSkip, 'x');
                addSkippedTicks(-75, 75, 15, 'y');
                yline(0,'--','LineWidth',1.5,'Color',.4*[1 1 1],'HandleVisibility','off'); hold on;
                ylabel('\DeltaCorrect (%)')
                axis square
                
            elseif cond==5
                % Plot line fit
                if plotLine==1
                    plot(mdl.xFitted(cond,:,block), mdl.yFitted(cond,:,block), 'Color', lineColor, 'LineWidth', 3, 'HandleVisibility', 'on'); hold on;
                end
                upFontSize(32, 0.01); legend; axis square

                % Plot scatterplot
                semY=std(mdl.yBlock(cond,:,block), 'omitnan') / sqrt(length(yBlock));
                counts=size(mdl.yBlock(cond,:,block),2);
                if size(mdl.yBlock(cond,:,block),1) > 1 % if no sem, don't shade
                    patchSaturationVal = 0.2;
                else
                    semY=zeros(counts,1);
                    patchSaturationVal = 0;
                end

                % Add data points and shaded error bar
                markerSize = 15;
                patchSaturationVal=0.14;
                % Data points
                shadedErrorBar(mdl.xBlock(cond,:,block), mdl.yBlock(cond,:,block), semY', 'patchSaturation', patchSaturationVal, 'lineprops', ...
                               {'Color', lineColor, 'LineStyle', 'none', 'LineWidth', 3, 'Marker', markerType, ...
                                'MarkerFaceColor', markerFaceColor, 'MarkerEdgeColor', edgeColor, 'MarkerSize', markerSize}); hold on;
                % Average
                [baseMeanMerged, conMeanMerged, inconMeanMerged] = computeMergedConditionMeans(mdl, block);
                [~, deltaMaskMergedForPlot] = computeDeltaFromConditionMeans(baseMeanMerged, conMeanMerged, inconMeanMerged);
                plot(barLength, repmat(deltaMaskMergedForPlot,1,numel(barLength)), '-', 'Color', lineColor, 'LineWidth', 3, 'HandleVisibility', 'off')

                % Add datapoints's count annotation
                zeroConstrastPoint=xBlock==0;
                if sum(zeroConstrastPoint)>0
                    nTrials=[20*ones(1,numel(zeroConstrastPoint)) 20*ones(1,numel(~zeroConstrastPoint))];
                    if cond>1
                            nTrials=[40*zeroConstrastPoint + 20*~zeroConstrastPoint];
                    end
                elseif sum(zeroConstrastPoint)>0 && plotAverageFlag==1
                    nTrials=[nBlocks*ones(1,numel(zeroConstrastPoint)) nBlocks*ones(1,numel(~zeroConstrastPoint))];
                    if cond>1
                            nTrials=[40*zeroConstrastPoint + 20*~zeroConstrastPoint];
                    end
                elseif sum(zeroConstrastPoint)==0
                     nTrials=[20*ones(1,numel(~zeroConstrastPoint))];
                elseif sum(zeroConstrastPoint)==0 && plotAverageFlag==1
                     nTrials=[nBlocks*ones(1,numel(~zeroConstrastPoint))];
                end
                upFontSize(32, 0.01);
                %annotateDataPoints(mdl.xBlock(cond,:,block), mdl.yBlock(cond,:,block), nTrials, markerFaceColor); hold on;
                xlim(xLimMerged)
                ylim([-75 75])
                xticks(xLimMerged(1):tickCfg.mergedMajor:xLimMerged(2));
                addSkippedTicks(xLimMerged(1), xLimMerged(2), tickCfg.mergedSkip, 'x');
                addSkippedTicks(-75, 75, 15, 'y');
                yline(0,'--','LineWidth',1.5,'Color',.4*[1 1 1],'HandleVisibility','off'); hold on;
                xlabel('Gabor contrast (%)')
                
                % Explicit legend handles for third panel
                % Purple = biasing, gray = masking
                hBiasLegend = plot(nan, nan, 's-', ...
                    'Color', [127, 0, 255]/255, ...
                    'MarkerFaceColor', [127, 0, 255]/255, ...
                    'MarkerEdgeColor', 'k', ...
                    'LineWidth', 3, ...
                    'MarkerSize', 15);
                
                hMaskLegend = plot(nan, nan, 's-', ...
                    'Color', [125, 125, 125]/255, ...
                    'MarkerFaceColor', [125, 125, 125]/255, ...
                    'MarkerEdgeColor', 'k', ...
                    'LineWidth', 3, ...
                    'MarkerSize', 15);
                
                legend([hBiasLegend, hMaskLegend], ...
                    {'Biasing', 'Masking'}, ...
                    'Location', 'southeast', ...
                    'NumColumns', 1, ...
                    'FontSize', 32);
                
                axis square


                [baseMeanMerged, conMeanMerged, inconMeanMerged] = computeMergedConditionMeans(mdl, block);
                [deltaBias, deltaMask] = computeDeltaFromConditionMeans(baseMeanMerged, conMeanMerged, inconMeanMerged);
                mdl.deltaBias(block)=deltaBias;
                mdl.deltaMask(block)=deltaMask;
                mdl.deltaBiasMerged(block)=deltaBias;
                mdl.deltaMaskMerged(block)=deltaMask;
                mdl.meanDeltaMerged(block,:) = [deltaBias, deltaMask];
                mdl.meanDeltaHeaders = {'Biasing', 'Masking'};
                title('');
                addDeltaSummaryText(deltaBias, deltaMask);
            end
        end
        
        % Customize the starting position and spacing
        startPos = [30, 30]; % Starting position in data coordinates
        xSpacing =5; % Horizontal spacing between columns
        ySpacing = 3; % Vertical spacing between rows

        % Call the function to create the table
        subplot(2, 3, 3); ax1=gca;
        %createCustomTable2(ax1, modelTypeStr, mdl.headers, mdl.fittedParams(block,:), startPos, xSpacing, ySpacing);
        
        if plotAverageFlag==1
            [nConditions, ~, nBlocks] = size(behavioralData.gaborContrasts(:, :, clusterBlocks));
            nBlocks=1;
            block=1;
        end
        upFontSize(21, .01);
        addBlockSuplabel(blockInfo.label);
        setPlotAnnotationFontSizes();
        if plotOpts.showDeltaPermutationStats && isgraphics(axMergedDeltaPanel)
            axes(axMergedDeltaPanel);
            permSeed = stableDeltaPermutationSeed(plotOpts.deltaPermutationBaseSeed, blockInfo.blockIdx, cluster);
            permResult = computeDeltaBiasPermutationFromEmpiricalPoints(...
                mdl.xBlock(2,:,block), mdl.yBlock(2,:,block), ...
                mdl.xBlock(3,:,block), mdl.yBlock(3,:,block), ...
                plotOpts.nDeltaPermutations, permSeed);
            if isempty(permResult.contrast)
                error('plotNakaRushtonFit5:NoDeltaPermutationContrasts', ...
                    'No exact con/incon contrasts for block %d.', blockInfo.blockIdx);
            end
            assertPermutationMatchesDisplayedBias(permResult, ...
                mdl.xBlock(4,:,block), mdl.yBlock(4,:,block));
            context = struct('experimentID', string(blockInfo.label), ...
                'modelRow', block, 'powerClusterID', cluster, ...
                'type', 'experiment');
            fprintf('Individual permutation ON | block %d | contrasts %d\n', ...
                blockInfo.blockIdx, numel(permResult.contrast));
            permAudit = addDeltaBiasPermutationVisualization(axMergedDeltaPanel, permResult, context);
            assertDeltaPermutationVisualizationAudit(permAudit, permResult, ...
                sprintf('block %d', blockInfo.blockIdx));
            validateDeltaPermutationLegend(axMergedDeltaPanel, sprintf('block %d', blockInfo.blockIdx));
            saveDeltaPermutationExampleFigure(gcf, plotOpts, 'individual');
            mdl = appendDeltaPermutationAudit(mdl, permAudit);
            printDeltaPermutationSummary('experiment', blockInfo.label, ...
                numel(permResult.contrast), permResult.observedMeanDeltaBias, ...
                permResult.rawOverallTwoSidedP, permResult.overallSignificant, ...
                permResult.rawOverallPositiveOneSidedP);
            printDeltaPermutationContrastDiagnostics('experiment', blockInfo.label, permResult);
        end
        if saveFlag
            assertClusterAnnotationPresent(optoStatsTextHandle, blockInfo, pageClusterLabel, block);
        end
        % Saving
        if saveFlag
            forceFigureSansSerif(fig);

            % Save one SVG for this experiment/block FIRST, before any PDF
            % helper can alter PaperPosition/PaperSize or the current figure.
            monkey = valueToChar(datastruct(blockInfo.datastructIdx).monkey);
            dateStr = valueToChar(blockInfo.date);
            runStr = valueToChar(blockInfo.run);
            chamber = valueToChar(datastruct(blockInfo.datastructIdx).chamber);
            monkeyNo = valueToChar(datastruct(blockInfo.datastructIdx).monkeyNo);

            if ispc
                svgRoot = 'Y:\';
            elseif contains(getenv('HOSTNAME'), 'psy.utexas.edu')
                svgRoot = '/eslab/data/';
            else
                svgRoot = pwd;
            end

            svgPath = fullfile(svgRoot, monkey, 'Meta', 'psychometrics', ...
                [chamber '-chamber'], modelTypeStr, 'svg');
            if ~exist(svgPath, 'dir')
                mkdir(svgPath);
            end

            svgName = sprintf('C%dM%sD%sR%s.svg', ...
                cluster, monkeyNo, dateStr, runStr);
            svgFile = fullfile(svgPath, svgName);

            % Force the same standalone visible geometry immediately before
            % export. This prevents caller options or a docked default from
            % silently reverting the figure to 560x420.
            set(fig, 'WindowStyle', 'normal', ...
                'Visible', 'on', ...
                'Units', 'pixels', ...
                'Position', targetFigurePosition, ...
                'Renderer', 'painters');
            set(findall(fig, '-property', 'FontName'), 'FontName', 'Arial');
            drawnow;

            figPixelsBeforeExport = getpixelposition(fig, true);
            if figPixelsBeforeExport(3) < 1000 || figPixelsBeforeExport(4) < 700
                error('plotNakaRushtonFit5:UnexpectedFigureCanvas', ...
                    ['Figure canvas remained %.0f x %.0f px before SVG export. ' ...
                    'Expected at least 1000 x 700 px.'], ...
                    figPixelsBeforeExport(3), figPixelsBeforeExport(4));
            end

            % Force the SVG paper canvas to match the actual on-screen figure.
            oldFigUnits = get(fig, 'Units');
            oldPaperUnits = get(fig, 'PaperUnits');
            oldPaperPosition = get(fig, 'PaperPosition');
            oldPaperSize = get(fig, 'PaperSize');
            oldPaperPositionMode = get(fig, 'PaperPositionMode');

            set(fig, 'Units', 'inches');
            figPositionInches = get(fig, 'Position');
            set(fig, ...
                'PaperUnits', 'inches', ...
                'PaperPositionMode', 'manual', ...
                'PaperPosition', [0 0 figPositionInches(3) figPositionInches(4)], ...
                'PaperSize', figPositionInches(3:4));

            figPixels = getpixelposition(fig, true);
            fprintf(['SVG WYSIWYG export | visible=%s | screen=%.0f x %.0f px ' ...
                '| paper=%.2f x %.2f in\n'], ...
                get(fig, 'Visible'), figPixels(3), figPixels(4), ...
                figPositionInches(3), figPositionInches(4));

            print(fig, svgFile, '-dsvg', '-painters');
            fprintf('Saved SVG: %s\n', svgFile);

            % Restore figure properties after SVG export.
            set(fig, ...
                'PaperUnits', oldPaperUnits, ...
                'PaperPosition', oldPaperPosition, ...
                'PaperSize', oldPaperSize, ...
                'PaperPositionMode', oldPaperPositionMode, ...
                'Units', oldFigUnits);

            % Existing report/PDF output: one page per rendered block.
            if isempty(reportState)
                saveCompressedPDFPage(savefilename, monkeyName, fig, appendPage);
            else
                reportState = stageReportPDFPage(reportState, fig);
            end
        end

        % Do not close merely because the figure was saved.  This preserves
        % visible WYSIWYG behavior.  The caller can still explicitly request
        % automatic closure through plotOpts.closeAfterRender.
        if ~saveFlag && plotOpts.closeAfterRender && isgraphics(fig)
            close(fig);
        end
    end
end


function assertClusterAnnotationPresent(annotationHandle, blockInfo, expectedClusterLabel, renderIdx)
    if nargin < 3 || isempty(expectedClusterLabel)
        error('plotNakaRushtonFit5:MissingExpectedClusterLabel', ...
            ['No expected cluster label was supplied for saved page %d ' ...
            '(%s, source block %d).'], ...
            renderIdx, char(string(blockInfo.label)), blockInfo.blockIdx);
    end

    expectedClusterText = char(string(expectedClusterLabel));
    expectedNormalized = normalizeClusterAnnotationText(expectedClusterText);
    actualString = '<unavailable>';
    actualNormalized = '<unavailable>';
    actualClass = class(annotationHandle);
    actualVisible = '<unavailable>';

    if isempty(annotationHandle) || ~isgraphics(annotationHandle)
        error('plotNakaRushtonFit5:MissingClusterAnnotationHandle', ...
            ['Cluster annotation handle is not a valid graphics object for ' ...
            'page %d (%s, source block %d). Expected normalized text: %s. ' ...
            'Actual class: %s. Actual string: %s'], ...
            renderIdx, char(string(blockInfo.label)), blockInfo.blockIdx, ...
            expectedNormalized, actualClass, actualString);
    end

    actualClass = class(annotationHandle);
    if isprop(annotationHandle, 'String')
        actualString = annotationStringToChar(get(annotationHandle, 'String'));
        actualNormalized = normalizeClusterAnnotationText(actualString);
    end
    if isprop(annotationHandle, 'Visible')
        actualVisible = char(string(get(annotationHandle, 'Visible')));
    end
    if ~strcmpi(actualVisible, 'on')
        error('plotNakaRushtonFit5:ClusterAnnotationNotVisible', ...
            ['Cluster annotation handle is not visible for page %d ' ...
            '(%s, source block %d). Expected normalized text: %s. ' ...
            'Visible: %s. Actual class: %s. Actual string: %s'], ...
            renderIdx, char(string(blockInfo.label)), blockInfo.blockIdx, ...
            expectedNormalized, actualVisible, actualClass, actualString);
    end

    if isempty(strfind(actualNormalized, expectedNormalized))
        error('plotNakaRushtonFit5:MissingClusterAnnotation', ...
            ['Missing expected visible cluster annotation for page %d ' ...
            '(%s, source block %d). Expected normalized text: %s. ' ...
            'Actual normalized string: %s. Actual class: %s. ' ...
            'Actual string: %s'], ...
            renderIdx, char(string(blockInfo.label)), blockInfo.blockIdx, ...
            expectedNormalized, actualNormalized, actualClass, actualString);
    end

    if renderIdx == 1
        clusterValue = regexprep(expectedNormalized, '^Cluster:\s*', '');
        fprintf('Cluster annotation validated | experiment %s | cluster %s | class %s\n', ...
            char(string(blockInfo.label)), clusterValue, actualClass);
    end
end

function textString = annotationStringToChar(textString)
    lineBreak = sprintf('\n');
    if iscell(textString)
        parts = cell(size(textString));
        for ii = 1:numel(textString)
            parts{ii} = char(string(textString{ii}));
        end
        textString = strjoin(parts(:)', lineBreak);
    else
        textString = char(string(textString));
    end
end

function textString = normalizeClusterAnnotationText(textString)
    lineBreak = sprintf('\n');
    textString = annotationStringToChar(textString);
    textString = strrep(textString, sprintf('\r\n'), lineBreak);
    textString = strrep(textString, sprintf('\r'), lineBreak);
    textString = regexprep(textString, '\s+', ' ');
    textString = strtrim(textString);
end
function moveLines()
% Get the current handles of all children
h = get(gca, 'Children');

% Initialize an array to hold the indices of lines to move
indicesToMove = [];

% Loop through each child handle
for idx = 1:length(h)
    % Check if the handle is a line object and isline
    if isa(h(idx), 'matlab.graphics.chart.primitive.Line') && ~strcmp(get(h(idx), 'LineStyle'),'none')
        % Add the index to the list
        indicesToMove = [indicesToMove, idx];
    end
end

% Extract the elements you want to move
elementsToMove = h(indicesToMove);%fliplr(indicesToMove));

% Remove these elements from the original array
h(indicesToMove) = [];

% Concatenate the removed elements at the end of the array
h = [h; elementsToMove];

% Set the new order of children
set(gca, 'Children', h);


%{
if billplot==1
    clf
    plot(mdl.xFitted(cond,:,block), mdl.yFitted(cond,:,block), 'Color', lineColor, 'LineWidth', 3, 'HandleVisibility', 'on'); hold on;
end
%}
end

function thresholdContrast=getThreshold(xData,yData,target);
% Find threshold value

% Find the index of the closest value
[~, index] = min(abs(yData - target));

% Get the closest value
thresholdContrast = xData(index);
end

function createCustomTable(ax, headers, fitParams, startPos, xSpacing, ySpacing)
    % Create a custom table on axes using text objects with subscripts.
    %
    % Parameters:
    %   ax - The axes handle where the table will be drawn.
    %   headers - A cell array of header names with subscripts.
    %   fitParams - A 1xN numeric array containing the values.
    %   startPos - Starting position of the table ([x, y]) in data coordinates.
    %   xSpacing - Horizontal spacing between columns.
    %   ySpacing - Vertical spacing between rows.

    % Ensure fitParams and headers are the same length
    if length(headers) ~= length(fitParams)
        error('Headers and fitParams must have the same length.');
    end
    fitParams=round(fitParams,1);
    % Plot headers
    for i = 1:length(headers)
        text(startPos(1) + (i - 1) * xSpacing, startPos(2), ...
             headers{i}, ...
             'Interpreter', 'tex', ...
             'HorizontalAlignment', 'center', ...
             'FontWeight', 'bold', ...
             'Parent', ax); % Use the specified axes
    end

    % Plot values
    for i = 1:length(fitParams)
        text(startPos(1) + (i - 1) * xSpacing, startPos(2) - ySpacing, ...
             sprintf('%.1f', fitParams(i)), ...
             'HorizontalAlignment', 'center', ...
             'Parent', ax); % Use the specified axes
    end
end

function createCustomTable2(ax, modelTypeStr, headers, fitParams, startPos, xSpacing, ySpacing)
    % Create a custom table on axes using text objects with subscripts.
    %
    % Parameters:
    %   ax - The axes handle where the table will be drawn.
    %   headers - A cell array of header names with subscripts.
    %   fitParams - A 1xN numeric array containing the values.
    %   startPos - Starting position of the table ([x, y]) in data coordinates.
    %   xSpacing - Horizontal spacing between columns.
    %   ySpacing - Vertical spacing between rows.
    
    % Exclude AICc for now
    % Strings to remove
    stringsToRemove = {'AICc^{total}', 'AUC^{con-incon}'};

    % Find the indices of the strings to remove
    idxToRemove = find(ismember(headers, stringsToRemove));
    headers(idxToRemove)=[];
    fitParams(idxToRemove)=[];

    % Resort for weibullfreeAll
    if strcmp(modelTypeStr,'weibullfreeAll')
        % Re-sort values to have rows = bl; con; incon.
        sortIdx=[1:4 13, 9:12 14, 5:8 15];
        % Sort
        fitParams=fitParams(sortIdx);
        headers=headers(sortIdx);
        
        % Indices for parameters
        idxA=1:5:numel(headers);
        idxB=idxA+1;
        idxAlpha=idxA+2;
        idxBeta=idxA+3;
        idxAUC=idxA+4;
        % Convert A & B to 100%
        fitParams([idxA,idxB]) = fitParams([idxA,idxB]) * 100;

        % A: Convert A to 100-A, calculate delta from baseline
        fitParams(idxA) = 100-fitParams(idxA);
        fitParams(idxA(2:3)) = fitParams(idxA(2:3)) - fitParams(idxA(1));
        % B: Flip, calculate delta from baseline
        fitParams(idxB(3))=fitParams(idxB(3))-50;
        fitParams(idxB(2))=-fitParams(idxB(3));
        % Alpha: calculate delta from baseline
        fitParams(idxAlpha(2:3)) = fitParams(idxAlpha(2:3)) - fitParams(idxAlpha(1));
        % Beta: calculate delta from baseline
        fitParams(idxBeta(2:3)) = fitParams(idxBeta(2:3)) - fitParams(idxBeta(1));
        % AUC: calculate delta from baseline
        fitParams(idxAUC(2:3)) = fitParams(idxAUC(2:3)) - fitParams(idxAUC(1));
        1;
    elseif strcmp(modelTypeStr,'weibull-beta')
        sortIdx=[1:4, 9:12, 5:8];
        % C50
        fitParams(10)=1-fitParams(6)-.5;
        fitParams(6)=fitParams(6)-.5;
        % Sort
        
        fitParams([2 6 10])=fitParams([2 6 10])*100;
        fitParams([1 5 9])=fitParams([1 5 9])*100;
        
        %A
        fitParams([5 9]) = fitParams([5 9]) - fitParams(1);
        %B
        fitParams([7 11]) = fitParams([7 11]) - fitParams(3);
        %Beta
        fitParams([8 12]) = fitParams([12 12]) - fitParams(4);
    end

    %fitParams = round(fitParams, 3, 'significant');

    % Determine the number of items per row
    numRows = 3;
    itemsPerRow = floor(length(headers) / numRows);

    % Plot headers and values
    for row = 1:numRows
        % Extract headers and fitParams for the current row
        startIdx = (row - 1) * itemsPerRow + 1;
        endIdx = min(row * itemsPerRow, length(headers));
        currentHeaders = headers(startIdx:endIdx);
        currentParams = fitParams(startIdx:endIdx);

        % Plot headers for the current row
        for i = 1:length(currentHeaders)
            text(startPos(1) + (i - 1) * xSpacing, ...
                 startPos(2) - (row - 1) * ySpacing-(5 * (row-1)), ...
                 currentHeaders{i}, ...
                 'Interpreter', 'tex', ...
                 'HorizontalAlignment', 'center', ...
                 'FontWeight', 'bold', ...
                 'Parent', ax); % Use the specified axes
        end
        
        for i = 1:length(currentParams)
            value = currentParams(i);
            if 1==1% Check if the value is an integer
                %formattedValue = sprintf('%.0f', value); % Format integers with .0

                % Force 2 significant figures while ensuring consistent formatting
                formattedValue = sprintf('%.2f', round(value, 2)); 
            else
            end
            text(startPos(1) + (i - 1) * xSpacing, ...
                 startPos(2) - row * ySpacing - (5 * (row - 1)), ...
                 formattedValue, ...
                 'HorizontalAlignment', 'center', ...
                 'Parent', ax); % Use the specified axes
        end


        
        %{
        for i = 1:length(currentParams)
            text(startPos(1) + (i - 1) * xSpacing, ...
                 startPos(2) - row * ySpacing-(5 * (row-1)), ...
                 sprintf('%.2f', currentParams(i)), ...
                 'HorizontalAlignment', 'center', ...
                 'Parent', ax); % Use the specified axes
        end
        %}
    end
end
function y = replaceNanSections(y)
    % Replace NaN sections in an array for a Weibull curve.
    % If NaN sections are flanked by equal values, replace them with those values.
    % For NaN sections at the end of the curve, use the last valid plateau value
    % if the last two valid values are similar within 0.01.
    y=real(y);
    % Ensure input is a row vector for easier indexing
    if iscolumn(y)
        y = y';
    end

    % Identify NaN sections
    isnanIndices = isnan(y);
    nanStarts = find(diff([0 isnanIndices]) == 1); % Start indices of NaN sections
    nanEnds = find(diff([isnanIndices 0]) == -1);  % End indices of NaN sections

    % Iterate through each NaN section
    for i = 1:length(nanStarts)
        startIdx = nanStarts(i);
        endIdx = nanEnds(i);

        % Case 1: Regular NaN sections with left and right flanks
        if startIdx > 1 && endIdx < length(y)
            leftFlank = y(startIdx - 1);
            rightFlank = y(endIdx + 1);

            % Replace NaNs if left and right flanks are equal (or nearly equal)
            if abs(leftFlank - rightFlank) <= 0.01
                y(startIdx:endIdx) = leftFlank;
            end
        % Case 2: NaN sections at the end of the curve
        elseif startIdx > 1 && endIdx == length(y)
            % Check the last two valid values to the left of the NaN section
            leftFlank1 = y(startIdx - 1);
            leftFlank2 = y(startIdx - 2);

            % Replace NaNs if the last two values are similar
            if abs(leftFlank1 - leftFlank2) <= 0.01
                y(startIdx:endIdx) = leftFlank1;
            end
        end
    end
end

function [shortXpad, shortYpad, missingIdxInLong] = padMissingX(shortX, shortY, longX)
% padMissingX
%
% Finds where shortX is missing one x-position relative to longX,
% then pads shortX and shortY with NaN at that position.
%
% Example:
%   shortX = [13 18 22 30 60];
%   longX  = [0 14 20 24 35 80];
%
%   returns missingIdxInLong = 1
%   shortXpad = [NaN 13 18 22 30 60]

    shortX = shortX(:)';   % force row
    shortY = shortY(:)';   % force row
    longX  = longX(:)';    % force row

    nShort = numel(shortX);
    nLong  = numel(longX);

    if nLong == nShort
        shortXpad = shortX;
        shortYpad = shortY;
        missingIdxInLong = [];
        return;
    end

    if nLong ~= nShort + 1
        error('This simple function expects longX to have exactly one more point than shortX.');
    end

    % Try removing each point from longX and see which removal best matches shortX
    err = nan(1, nLong);

    for candidateMissingIdx = 1:nLong

        longX_withoutCandidate = longX;
        longX_withoutCandidate(candidateMissingIdx) = [];

        % Error between shortX and longX with this candidate point removed
        err(candidateMissingIdx) = nansum(abs(shortX - longX_withoutCandidate));

    end

    % Best missing index
    [~, missingIdxInLong] = min(err);

    % Pad shortX and shortY with NaN at that position
    shortXpad = nan(1, nLong);
    shortYpad = nan(1, nLong);

    keepIdx = true(1, nLong);
    keepIdx(missingIdxInLong) = false;

    shortXpad(keepIdx) = shortX;
    shortYpad(keepIdx) = shortY;
end

function sideData = getPreMergedSideData(mdl, block)
    xBaseline = rmnan(mdl.xBaselinePreMerge(block, :));
    yBaseline = rmnan(mdl.yBaselinePreMerge(block, :));

    xH = rmnan(mdl.xHorizontalOptoPreMerge(block, :));
    yH = rmnan(mdl.yHorizontalOptoPreMerge(block, :));
    tagH = rmnan(mdl.visualTagHorizontalOptoPreMerge(block, :));
    congrH = rmnan(mdl.congruencyHorizontalOptoPreMerge(block, :));

    xV = rmnan(mdl.xVerticalOptoPreMerge(block, :));
    yV = rmnan(mdl.yVerticalOptoPreMerge(block, :));
    tagV = rmnan(mdl.visualTagVerticalOptoPreMerge(block, :));
    congrV = rmnan(mdl.congruencyVerticalOptoPreMerge(block, :));

    xOpto = [xH, xV];
    yOpto = [yH, yV];
    tagOpto = [tagH, tagV];
    congrOpto = [congrH, congrV];

    idxBaseHorizontal = xBaseline <= 0;
    idxBaseVertical = xBaseline >= 0;

    % Split by visual stimulus tag, not contrast sign. This is especially
    % important at x=0, where the contrast sign does not identify the side.
    idxConHorizontal = tagOpto == 0 & congrOpto == 1;
    idxInconHorizontal = tagOpto == 0 & congrOpto == -1;

    idxConVertical = tagOpto == 90 & congrOpto == 1;
    idxInconVertical = tagOpto == 90 & congrOpto == -1;

    sideData.horizontal = makeVisualSideData('Horizontal', ...
        abs(xBaseline(idxBaseHorizontal)), yBaseline(idxBaseHorizontal), ...
        abs(xOpto(idxConHorizontal)), yOpto(idxConHorizontal), ...
        abs(xOpto(idxInconHorizontal)), yOpto(idxInconHorizontal));

    sideData.vertical = makeVisualSideData('Vertical', ...
        xBaseline(idxBaseVertical), yBaseline(idxBaseVertical), ...
        xOpto(idxConVertical), yOpto(idxConVertical), ...
        xOpto(idxInconVertical), yOpto(idxInconVertical));
end

function sideData = makeVisualSideData(sideName, xBaseline, yBaseline, xConOpto, yConOpto, xInconOpto, yInconOpto)
    [xBaseline, yBaseline] = cleanAndSortXY(xBaseline, yBaseline);
    [xConOpto, yConOpto] = cleanAndSortXY(xConOpto, yConOpto);
    [xInconOpto, yInconOpto] = cleanAndSortXY(xInconOpto, yInconOpto);

    sideData = struct();
    sideData.sideName = sideName;
    sideData.xBaseline = xBaseline;
    sideData.yBaseline = yBaseline;
    sideData.wBaseline = getPlotMeanWeights(xBaseline, 'sideBaseline');
    sideData.xConOpto = xConOpto;
    sideData.yConOpto = yConOpto;
    sideData.wConOpto = getPlotMeanWeights(xConOpto, 'sideOpto');
    sideData.xInconOpto = xInconOpto;
    sideData.yInconOpto = yInconOpto;
    sideData.wInconOpto = getPlotMeanWeights(xInconOpto, 'sideOpto');
end

function [mdl, sideData] = plotSidePsychometricPanel(mdl, block, sideData, xLimMerged, meanBar, tickCfg, titleStr)
    hold on;

    style = getPreMergedPlotStyle();
    markerSize = 12;
    lineWidth = 3;

    yline(50, '--', ...
        'LineWidth', 1.5, ...
        'Color', .4 * [1 1 1], ...
        'HandleVisibility', 'off');

    hBaseline = plotConditionSeries(sideData.xBaseline, sideData.yBaseline, ...
        style.baselineColor, 'o', style.baselineFaceColor, markerSize, lineWidth, 'Baseline', 'on');

    hCon = plotConditionSeries(sideData.xConOpto, sideData.yConOpto, ...
        style.conColor, '^', style.conColor, markerSize, lineWidth, 'Con-Opto', 'on');

    hIncon = plotConditionSeries(sideData.xInconOpto, sideData.yInconOpto, ...
        style.inconColor, 'v', style.inconColor, markerSize, lineWidth, 'Incon-Opto', 'on');

    sideData.meanBaseline = weightedMeanForPlot(sideData.xBaseline, sideData.yBaseline, 'sideBaseline');
    sideData.meanConOpto = weightedMeanForPlot(sideData.xConOpto, sideData.yConOpto, 'sideOpto');
    sideData.meanInconOpto = weightedMeanForPlot(sideData.xInconOpto, sideData.yInconOpto, 'sideOpto');
    meanVals = [sideData.meanBaseline, sideData.meanConOpto, sideData.meanInconOpto];

    if strcmp(sideData.sideName, 'Horizontal')
        mdl.meanBaselineHorizontal(block) = sideData.meanBaseline;
        mdl.meanConOptoHorizontal(block) = sideData.meanConOpto;
        mdl.meanInconOptoHorizontal(block) = sideData.meanInconOpto;
        mdl.meanPanel1Horizontal(block,:) = meanVals;
        mdl.meanPsychometricHorizontal(block,:) = meanVals;
    else
        mdl.meanBaselineVertical(block) = sideData.meanBaseline;
        mdl.meanConOptoVertical(block) = sideData.meanConOpto;
        mdl.meanInconOptoVertical(block) = sideData.meanInconOpto;
        mdl.meanPanel1Vertical(block,:) = meanVals;
        mdl.meanPsychometricVertical(block,:) = meanVals;
    end

    mdl.meanPanel1Headers = {'Baseline', 'ConOpto', 'InconOpto'};
    mdl.meanPanel1SideHeaders = {'Horizontal', 'Vertical'};
    mdl.meanPsychometricHeaders = {'Baseline', 'ConOpto', 'InconOpto'};

    meanValsPlot = jitterOverlappingMeans(meanVals, .5);
    plotMeanBarAtY(meanValsPlot(1), meanBar.rightEdgeRange, style.baselineColor, lineWidth);
    plotMeanBarAtY(meanValsPlot(2), meanBar.rightEdgeRange, style.conColor, lineWidth);
    plotMeanBarAtY(meanValsPlot(3), meanBar.rightEdgeRange, style.inconColor, lineWidth);

    xlim(xLimMerged);
    ylim([0 100]);
    xticks(xLimMerged(1):tickCfg.mergedMajor:xLimMerged(2));
    addSkippedTicks(xLimMerged(1), xLimMerged(2), tickCfg.mergedSkip, 'x');
    addSkippedTicks(0, 100, 10, 'y');

    xlabel('Gabor contrast (%)');
    ylabel('Correct (%)');
    title(titleStr);

    legend([hBaseline, hCon, hIncon], ...
        {'Baseline', 'Con-Opto', 'Incon-Opto'}, ...
        'Location', 'southeast', ...
        'NumColumns', 1, ...
        'FontSize', 24);

    axis square;
    upFontSize(32, 0.01);
end

function mdl = plotSideDeltaPanel(mdl, block, sideData, xLimMerged, meanBar, tickCfg, titleStr, sideFieldName, blockInfo, baselineMode, fittedDeltaInput)
    hold on;

    lineWidth = 3;
    markerSize = 15;
    biasColor = [127, 0, 255] / 255;
    maskColor = [125, 125, 125] / 255;

    if nargin < 11
        fittedDeltaInput = [];
    end

    [xBias, yBias, xMask, yMask, deltaAudit] = computeSideDeltaData(sideData);
    deltaAudit = addDeltaPointAuditMetadata( ...
        deltaAudit, block, blockInfo, baselineMode);
    if ~isfield(mdl, 'deltaPointAudit') || isempty(mdl.deltaPointAudit)
        mdl.deltaPointAudit = emptyDeltaPointAuditTable();
    elseif ~isequal(mdl.deltaPointAudit.Properties.VariableNames, ...
            deltaAudit.Properties.VariableNames)
        warning('plotNakaRushtonFit5:DeltaAuditSchemaChanged', ...
            ['Existing mdl.deltaPointAudit schema does not match the ' ...
            'current strict audit schema; resetting the table.']);
        mdl.deltaPointAudit = emptyDeltaPointAuditTable();
    end
    mdl.deltaPointAudit = [mdl.deltaPointAudit; deltaAudit];
    printDeltaPointAuditSummary(deltaAudit, blockInfo.label);

    plotFittedSideDeltaCurves(fittedDeltaInput, biasColor, maskColor, lineWidth);

    hBias = plotConditionSeries(xBias, yBias, ...
        biasColor, 's', biasColor, markerSize, lineWidth, 'Biasing', 'on');

    hMask = plotConditionSeries(xMask, yMask, ...
        maskColor, 's', maskColor, markerSize, lineWidth, 'Masking', 'on');

    [deltaBias, deltaMask] = computeDeltaFromConditionMeans( ...
        sideData.meanBaseline, sideData.meanConOpto, sideData.meanInconOpto);

    mdl.(['deltaBias' sideFieldName])(block) = deltaBias;
    mdl.(['deltaMask' sideFieldName])(block) = deltaMask;
    mdl.(['meanDelta' sideFieldName])(block,:) = [deltaBias, deltaMask];
    mdl.meanDeltaHeaders = {'Biasing', 'Masking'};

    deltaMeansPlot = jitterOverlappingMeans([deltaBias, deltaMask], .5);
    plotMeanBarAtY(deltaMeansPlot(1), meanBar.rightEdgeRange, biasColor, lineWidth);
    plotMeanBarAtY(deltaMeansPlot(2), meanBar.rightEdgeRange, maskColor, lineWidth);

    xlim(xLimMerged);
    ylim([-75 75]);
    xticks(xLimMerged(1):tickCfg.mergedMajor:xLimMerged(2));
    addSkippedTicks(xLimMerged(1), xLimMerged(2), tickCfg.mergedSkip, 'x');
    addSkippedTicks(-75, 75, 15, 'y');
    yline(0, '--', 'LineWidth', 1.5, 'Color', .4 * [1 1 1], 'HandleVisibility', 'off');

    xlabel('Gabor contrast (%)');
    ylabel('\DeltaCorrect (%)');
    title('');

    legend([hBias, hMask], ...
        {'Biasing', 'Masking'}, ...
        'Location', 'southeast', ...
        'NumColumns', 1, ...
        'FontSize', 24);

    axis square;
    upFontSize(32, 0.01);
    addDeltaSummaryText(deltaBias, deltaMask);
end


function plotFittedSideDeltaCurves(fittedDeltaInput, biasColor, maskColor, lineWidth)
    if isempty(fittedDeltaInput) || ~isstruct(fittedDeltaInput) || ...
            ~isfield(fittedDeltaInput, 'x') || ...
            ~isfield(fittedDeltaInput, 'viewCurves')
        return;
    end
    viewCurves = fittedDeltaInput.viewCurves;
    requiredFields = {'baseline', 'con', 'incon'};
    if any(~isfield(viewCurves, requiredFields))
        return;
    end

    sourcePanel = '';
    if isfield(fittedDeltaInput, 'sourcePanel')
        sourcePanel = char(lower(fittedDeltaInput.sourcePanel));
    elseif isfield(viewCurves, 'sourcePanel')
        sourcePanel = char(lower(viewCurves.sourcePanel));
    end

    switch sourcePanel
        case 'horizontal'
            biasHorizontal = viewCurves.con - viewCurves.incon;
            maskHorizontal = viewCurves.baseline - 0.5 .* (viewCurves.con + viewCurves.incon);
            biasCurve = biasHorizontal;
            maskCurve = maskHorizontal;
        case 'vertical'
            biasVertical = viewCurves.con - viewCurves.incon;
            maskVertical = viewCurves.baseline - 0.5 .* (viewCurves.con + viewCurves.incon);
            biasCurve = biasVertical;
            maskCurve = maskVertical;
        case 'merged'
            biasMerged = viewCurves.con - viewCurves.incon;
            maskMerged = viewCurves.baseline - 0.5 .* (viewCurves.con + viewCurves.incon);
            biasCurve = biasMerged;
            maskCurve = maskMerged;
        otherwise
            error('plotNakaRushtonFit5:SignedBX0DeltaSourceMissing', ...
                'Signed-BX0 fitted delta curves require a horizontal, vertical, or merged sourcePanel.');
    end

    xGrid = fittedDeltaInput.x;
    plot(xGrid, biasCurve, '-', ...
        'Color', biasColor, ...
        'LineWidth', lineWidth, ...
        'Tag', 'SignedBX0FittedBiasCurve', ...
        'HandleVisibility', 'off');
    plot(xGrid, maskCurve, '-', ...
        'Color', maskColor, ...
        'LineWidth', lineWidth, ...
        'Tag', 'SignedBX0FittedMaskCurve', ...
        'HandleVisibility', 'off');
end
function audit = addDeltaPointAuditMetadata(audit, block, blockInfo, baselineMode)
    if isempty(audit) || height(audit) == 0
        return;
    end
    audit.sessionRowIndex(:) = block;
    audit.blockIndex(:) = blockInfo.blockIdx;
    audit.experimentID(:) = string(blockInfo.label);
    audit.baselineMode(:) = string(baselineMode);
end

function printDeltaPointAuditSummary(audit, experimentID)
    experimentLabel = char(string(experimentID));
    if isempty(audit) || height(audit) == 0
        warning('plotNakaRushtonFit5:NoDeltaPairs', ...
            '%s: no side-specific delta pairs were plotted.', experimentLabel);
        return;
    end
    metrics = unique(audit.deltaMetric, 'stable');
    for metricIdx = 1:numel(metrics)
        metricRows = audit.deltaMetric == metrics(metricIdx);
        metricAudit = audit(metricRows, :);
        plottedRows = metricAudit.isPlotted;
        expectedCount = firstFiniteValue(metricAudit.expectedPointCount);
        plottedCount = sum(plottedRows);
        skippedCount = sum(~plottedRows);
        sideLabel = char(metricAudit.visualSide(1));
        metricLabel = char(metrics(metricIdx));
        methodLabel = char(strjoin(unique(metricAudit.pairingMethod), ','));
        fprintf(['Delta point audit: %s | side=%s | metric=%s | ' ...
            'plotted=%d/%g | skipped=%d | method=%s | ' ...
            'max contrast diff=%0.3g | x=%s | delta=%s\n'], ...
            experimentLabel, sideLabel, metricLabel, ...
            plottedCount, expectedCount, skippedCount, methodLabel, ...
            maxFiniteValue(metricAudit.maxAbsContrastDiff(plottedRows)), ...
            formatAuditNumberList(metricAudit.xDelta(plottedRows)), ...
            formatAuditNumberList(metricAudit.deltaValue(plottedRows)));
        if skippedCount > 0 || plottedCount < expectedCount
            skippedReasons = unique(metricAudit.reasonExcluded(~plottedRows));
            skippedReasons = skippedReasons(strlength(skippedReasons) > 0);
            warning('plotNakaRushtonFit5:DeltaPointMismatch', ...
                ['%s %s %s plotted %d of expected %g delta points. ' ...
                'Skipped reasons: %s'], experimentLabel, sideLabel, ...
                metricLabel, plottedCount, expectedCount, ...
                char(strjoin(skippedReasons, '; ')));
        end
    end
end

function value = maxFiniteValue(values)
    values = values(isfinite(values));
    if isempty(values)
        value = NaN;
    else
        value = max(values);
    end
end

function value = firstFiniteValue(values)
    idx = find(isfinite(values), 1, 'first');
    if isempty(idx)
        value = NaN;
    else
        value = values(idx);
    end
end

function label = formatAuditNumberList(values)
    values = values(:)';
    values = values(isfinite(values));
    if isempty(values)
        label = '[]';
    else
        label = ['[', strtrim(sprintf('%0.3g ', values)), ']'];
    end
end

function addDeltaSummaryText(deltaBias, deltaMask)
    ax = gca;
    axPos = get(ax, 'Position');
    textWidth = 0.58 * axPos(3);
    textX = axPos(1) + 0.5 * axPos(3) - 0.5 * textWidth;
    textY = max(0.001, axPos(2) - 0.23);
    textHeight = 0.085;

    annotation(gcf, 'textbox', [textX, textY, textWidth, textHeight], ...
        'String', sprintf('biasing: %.1f%%\nmasking: %.1f%%', deltaBias, deltaMask), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', ...
        'FontSize', 13.5, ...
        'Interpreter', 'tex', ...
        'EdgeColor', 'none', ...
        'BackgroundColor', 'none', ...
        'FitBoxToText', 'off', ...
        'Tag', 'PlotNakaDeltaText');
end

function annotationHandle = addOptoStatsText(ax, optoText)
    axPos = get(ax, 'Position');
    statsRight = axPos(1) - 0.012;
    statsX = 0.002;
    statsWidth = max(0.02, statsRight - statsX);
    statsHeight = 0.27;
    statsY = axPos(2) + 0.5 * axPos(4) - 0.5 * statsHeight;

    annotationHandle = annotation(gcf, 'textbox', [statsX, statsY, statsWidth, statsHeight], ...
        'String', optoText, ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', 10.5, ...
        'Interpreter', 'tex', ...
        'EdgeColor', 'none', ...
        'BackgroundColor', 'none', ...
        'FitBoxToText', 'off', ...
        'Tag', 'PlotNakaStatsText');
end

function displayParams = getSignedBX0DisplayParamsForBlock(fitParamRow)
    displayParams = squeeze(fitParamRow);
    displayParams = displayParams(:)';
    if numel(displayParams) < 11 || any(~isfinite(displayParams(1:11)))
        error('plotNakaRushtonFit5:SignedBX0ParamCount', ...
            'weibullSignedBX0 plotting requires 11 full-model parameters.');
    end
    displayParams = displayParams(1:11);
end

function globalDeltaX0 = getSignedBX0GlobalDeltaX0(mdl, block, signedBX0DisplayParams)
    globalDeltaX0 = NaN;
    if isfield(mdl, 'signedBX0') && isfield(mdl.signedBX0, 'deltaX0') && ...
            numel(mdl.signedBX0.deltaX0) >= block
        globalDeltaX0 = mdl.signedBX0.deltaX0(block);
    end
    if ~isfinite(globalDeltaX0)
        globalDeltaX0 = signedBX0DisplayParams(11);
    end
    if numel(globalDeltaX0) ~= 1 || ~isfinite(globalDeltaX0)
        error('plotNakaRushtonFit5:InvalidSignedBX0GlobalX0', ...
            'Expected exactly one finite global deltaX0 for block %d.', block);
    end
end

function [mdl, panelFits] = fitSignedBX0PanelFitsForBlock(mdl, block, sideData, globalDeltaX0, primaryParams, blockInfo)
    if numel(globalDeltaX0) ~= 1 || ~isfinite(globalDeltaX0)
        error('plotNakaRushtonFit5:InvalidSignedBX0GlobalX0', ...
            'Panel fitting requires one fixed global deltaX0.');
    end

    nRows = size(mdl.fittedParams, 1);
    mdl = ensureSignedBX0PanelFitStorage(mdl, nRows);
    panelFits = struct();
    panelFits.globalDeltaX0 = globalDeltaX0;

    panelNames = {'horizontal', 'vertical', 'merged'};
    for panelIdx = 1:numel(panelNames)
        panelName = panelNames{panelIdx};
        switch panelName
            case 'horizontal'
                panelData = buildSignedBX0PanelFitData(panelName, ...
                    sideData.horizontal.xBaseline, sideData.horizontal.yBaseline, 'sideBaseline', ...
                    sideData.horizontal.xConOpto, sideData.horizontal.yConOpto, 'sideOpto', ...
                    sideData.horizontal.xInconOpto, sideData.horizontal.yInconOpto, 'sideOpto');
            case 'vertical'
                panelData = buildSignedBX0PanelFitData(panelName, ...
                    sideData.vertical.xBaseline, sideData.vertical.yBaseline, 'sideBaseline', ...
                    sideData.vertical.xConOpto, sideData.vertical.yConOpto, 'sideOpto', ...
                    sideData.vertical.xInconOpto, sideData.vertical.yInconOpto, 'sideOpto');
            case 'merged'
                panelData = buildSignedBX0MergedPanelFitData(mdl, block);
        end

        fitResult = fitSignedBX0SinglePanel(panelData, globalDeltaX0, primaryParams, blockInfo);
        fitResult.sourcePanel = panelName;
        assertSignedBX0PanelFit(fitResult, panelName, globalDeltaX0);
        panelFits.(panelName) = fitResult;
        mdl = storeSignedBX0PanelFit(mdl, block, panelName, fitResult);
        fprintf(['weibullSignedBX0 panel fit audit | experiment %s | panel %s | ' ...
            'sourceChecksum %.12g | fitChecksum %.12g\n'], ...
            char(blockInfo.label), panelName, fitResult.sourceChecksum, fitResult.fitChecksum);
    end
end

function mdl = ensureSignedBX0PanelFitStorage(mdl, nRows)
    if ~isfield(mdl, 'signedBX0')
        mdl.signedBX0 = struct();
    end
    if ~isfield(mdl.signedBX0, 'panelFits') || isempty(mdl.signedBX0.panelFits)
        mdl.signedBX0.panelFits = struct();
    end

    panelNames = {'horizontal', 'vertical', 'merged'};
    for panelIdx = 1:numel(panelNames)
        panelName = panelNames{panelIdx};
        if ~isfield(mdl.signedBX0.panelFits, panelName) || ...
                isempty(mdl.signedBX0.panelFits.(panelName))
            mdl.signedBX0.panelFits.(panelName) = struct();
        end
        panel = mdl.signedBX0.panelFits.(panelName);
        panel.fitParams = initializePanelField(panel, 'fitParams', nRows, 10, NaN);
        panel.nLL = initializePanelField(panel, 'nLL', nRows, 1, NaN);
        panel.fitStatus = initializePanelCellField(panel, 'fitStatus', nRows, '');
        panel.bestStartIndex = initializePanelField(panel, 'bestStartIndex', nRows, 1, NaN);
        panel.exitFlag = initializePanelField(panel, 'exitFlag', nRows, 1, NaN);
        panel.boundHit = initializePanelField(panel, 'boundHit', nRows, 1, false);
        panel.boundHitFields = initializePanelCellField(panel, 'boundHitFields', nRows, {});
        panel.maxSlopeOverall = initializePanelField(panel, 'maxSlopeOverall', nRows, 1, NaN);
        panel.slopeCapActive = initializePanelField(panel, 'slopeCapActive', nRows, 1, false);
        panel.sourceChecksum = initializePanelField(panel, 'sourceChecksum', nRows, 1, NaN);
        panel.fitChecksum = initializePanelField(panel, 'fitChecksum', nRows, 1, NaN);
        panel.globalDeltaX0 = initializePanelField(panel, 'globalDeltaX0', nRows, 1, NaN);
        panel.parameterNames = signedBX0PanelParameterNames();
        mdl.signedBX0.panelFits.(panelName) = panel;
    end
end

function value = initializePanelField(panel, fieldName, nRows, nCols, fillValue)
    if isfield(panel, fieldName) && ~isempty(panel.(fieldName)) && ...
            size(panel.(fieldName), 1) >= nRows && size(panel.(fieldName), 2) >= nCols
        value = panel.(fieldName);
    else
        value = repmat(fillValue, nRows, nCols);
    end
end

function value = initializePanelCellField(panel, fieldName, nRows, fillValue)
    if isfield(panel, fieldName) && ~isempty(panel.(fieldName)) && ...
            numel(panel.(fieldName)) >= nRows
        value = panel.(fieldName);
    else
        value = repmat({fillValue}, nRows, 1);
    end
end

function names = signedBX0PanelParameterNames()
    names = {'A_baseline', 'alpha_baseline', 'beta_baseline', ...
        'A_con', 'alpha_con', 'beta_con', ...
        'A_incon', 'alpha_incon', 'beta_incon', 'deltaB_panel'};
end

function panelData = buildSignedBX0MergedPanelFitData(mdl, block)
    xBaseline = rmnan(mdl.xBaseline(block, :));
    yBaseline = rmnan(mdl.yBaseline(block, :));
    xCon = rmnan(mdl.xConOpto(block, :));
    yCon = rmnan(mdl.yConOpto(block, :));
    xIncon = rmnan(mdl.xInconOpto(block, :));
    yIncon = rmnan(mdl.yInconOpto(block, :));

    idxBaseline = xBaseline >= 0;
    idxCon = xCon >= 0;
    idxIncon = xIncon >= 0;
    panelData = buildSignedBX0PanelFitData('merged', ...
        abs(xBaseline(idxBaseline)), yBaseline(idxBaseline), 'merged', ...
        abs(xCon(idxCon)), yCon(idxCon), 'merged', ...
        abs(xIncon(idxIncon)), yIncon(idxIncon), 'merged');
end

function panelData = buildSignedBX0PanelFitData(panelName, xBaseline, yBaseline, baselineWeightMode, xCon, yCon, conWeightMode, xIncon, yIncon, inconWeightMode)
    [xBaseline, yBaseline] = cleanAndSortXY(abs(xBaseline), yBaseline);
    [xCon, yCon] = cleanAndSortXY(abs(xCon), yCon);
    [xIncon, yIncon] = cleanAndSortXY(abs(xIncon), yIncon);

    nBaseline = getPlotMeanWeights(xBaseline, baselineWeightMode);
    nCon = getPlotMeanWeights(xCon, conWeightMode);
    nIncon = getPlotMeanWeights(xIncon, inconWeightMode);

    panelData = struct();
    panelData.panelName = panelName;
    panelData.xBaseline = xBaseline;
    panelData.yBaseline = yBaseline;
    panelData.nBaseline = nBaseline;
    panelData.successBaseline = round((yBaseline ./ 100) .* nBaseline);
    panelData.xCon = xCon;
    panelData.yCon = yCon;
    panelData.nCon = nCon;
    panelData.successCon = round((yCon ./ 100) .* nCon);
    panelData.xIncon = xIncon;
    panelData.yIncon = yIncon;
    panelData.nIncon = nIncon;
    panelData.successIncon = round((yIncon ./ 100) .* nIncon);
    panelData.sourceChecksum = computeSignedBX0PanelSourceChecksum(panelData);

    if isempty(xBaseline) || isempty(xCon) || isempty(xIncon)
        error('plotNakaRushtonFit5:SignedBX0PanelSourceMissing', ...
            'weibullSignedBX0 panel %s has missing baseline/con/incon source data.', panelName);
    end
end

function checksum = computeSignedBX0PanelSourceChecksum(panelData)
    values = [panelData.xBaseline, panelData.yBaseline, panelData.nBaseline, panelData.successBaseline, ...
        panelData.xCon, panelData.yCon, panelData.nCon, panelData.successCon, ...
        panelData.xIncon, panelData.yIncon, panelData.nIncon, panelData.successIncon];
    values = values(isfinite(values));
    checksum = sum(values .* (1:numel(values)));
end

function fitResult = fitSignedBX0SinglePanel(panelData, globalDeltaX0, primaryParams, blockInfo)
    [initialParams] = getWeibullSignedBX0InitParams();
    [lb, ub] = getWeibullSignedBX0PanelBounds();
    params0 = initialParams(1:10);
    if numel(primaryParams) >= 10 && all(isfinite(primaryParams(1:10)))
        params0 = primaryParams(1:10);
    end
    params0 = min(max(params0, lb), ub);

    starts = makeSignedBX0PanelDeterministicStarts(params0, lb, ub);
    bestNLL = Inf;
    bestParams = nan(1, 10);
    bestStartIndex = NaN;
    bestExitFlag = NaN;
    exitFlags = nan(size(starts, 1), 1);
    attemptedMaxSlopes = nan(size(starts, 1), 1);

    for startIdx = 1:size(starts, 1)
        [candidateParams, candidateNLL, exitFlag] = fitSignedBX0PanelWithExit( ...
            starts(startIdx, :), lb, ub, panelData, globalDeltaX0);
        exitFlags(startIdx) = exitFlag;
        slopeDiagnostics = getSignedBX0PanelSlopeDiagnostics(candidateParams, globalDeltaX0);
        attemptedMaxSlopes(startIdx) = slopeDiagnostics.maxSlopeOverall;
        if numel(candidateParams) == 10 && all(isfinite(candidateParams)) && ...
                isfinite(candidateNLL) && slopeDiagnostics.isValid && candidateNLL < bestNLL
            bestNLL = candidateNLL;
            bestParams = candidateParams;
            bestStartIndex = startIdx;
            bestExitFlag = exitFlag;
        end
    end

    if ~isfinite(bestNLL)
        error('plotNakaRushtonFit5:SignedBX0PanelFitFailed', ...
            ['All conditional weibullSignedBX0 panel starts failed for %s panel %s. ' ...
            'Exit flags: %s. Lowest attempted max slope: %.4g.'], ...
            blockInfo.label, panelData.panelName, mat2str(exitFlags'), ...
            min(attemptedMaxSlopes, [], 'omitnan'));
    end

    diagnostics = makeSignedBX0PanelFitDiagnostics(bestParams, bestStartIndex, ...
        bestExitFlag, lb, ub, globalDeltaX0);
    fitResult = struct( ...
        'fitParams', bestParams, ...
        'nLL', bestNLL, ...
        'fitStatus', 'ok', ...
        'bestStartIndex', bestStartIndex, ...
        'exitFlag', bestExitFlag, ...
        'boundHit', diagnostics.boundHit, ...
        'boundHitFields', {diagnostics.boundHitFields}, ...
        'maxSlopeOverall', diagnostics.maxSlopeOverall, ...
        'slopeCapActive', diagnostics.slopeCapActive, ...
        'sourceChecksum', panelData.sourceChecksum, ...
        'fitChecksum', sum(bestParams .* (1:10)) + 11 .* globalDeltaX0, ...
        'globalDeltaX0', globalDeltaX0, ...
        'parameterNames', {signedBX0PanelParameterNames()});
end

function starts = makeSignedBX0PanelDeterministicStarts(params0, lb, ub)
    starts = repmat(params0(:)', 7, 1);
    starts(2, [3 6 9]) = 1.5;
    starts(3, [3 6 9]) = 5;
    starts(4, 10) = 10;
    starts(5, 10) = -10;
    starts(6, [3 6 9 10]) = [2.5 2.5 2.5 5];
    starts(7, [3 6 9 10]) = [6 6 6 -5];
    starts = min(max(starts, lb), ub);
    starts = unique(starts, 'rows', 'stable');
end

function [fittedParams, nLL, exitFlag] = fitSignedBX0PanelWithExit(params0, lb, ub, panelData, globalDeltaX0)
    toParams = @(u) lb + u(:)' .* (ub - lb);
    toUnit = @(p) (p(:)' - lb) ./ (ub - lb);
    u0 = min(max(toUnit(params0), 0), 1);
    obj = @(u) signedBX0PanelObjective(toParams(min(max(u, 0), 1)), ...
        panelData, globalDeltaX0);

    haveFMC = exist('fmincon','file') == 2;
    if haveFMC
        fopts = optimoptions('fmincon', ...
            'Algorithm','interior-point', ...
            'Display','off', ...
            'MaxFunctionEvaluations', 5000, ...
            'FiniteDifferenceType','central');
        [uFit, ~, exitFlag] = fmincon(obj, u0, [], [], [], [], ...
            zeros(size(u0)), ones(size(u0)), [], fopts);
    else
        opts = optimset('Display', 'off', 'MaxFunEvals', 5000, ...
            'MaxIter', 5000);
        [uFit, ~, exitFlag] = fminsearchbnd(obj, u0, zeros(size(u0)), ...
            ones(size(u0)), opts);
    end

    fittedParams = toParams(min(max(uFit, 0), 1));
    nLL = signedBX0PanelObjective(fittedParams, panelData, globalDeltaX0);
end

function nLL = signedBX0PanelObjective(panelParams, panelData, globalDeltaX0)
    epsilon = 1e-10;
    if numel(panelParams) ~= 10 || any(~isfinite(panelParams(:))) || ...
            ~isfinite(globalDeltaX0)
        nLL = 1e12;
        return;
    end

    slopeDiagnostics = getSignedBX0PanelSlopeDiagnostics(panelParams, globalDeltaX0);
    if ~slopeDiagnostics.isValid
        excess = slopeDiagnostics.slopeExcess;
        excess(~isfinite(excess)) = 5.0;
        nLL = 1e12 + 1e6 .* sum(excess .^ 2);
        return;
    end

    curves = predictSignedBX0PanelFitCurvesForData(panelParams, globalDeltaX0, ...
        panelData.xBaseline, panelData.xCon, panelData.xIncon);
    pBaseline = min(max(curves.baseline ./ 100, epsilon), 1 - epsilon);
    pCon = min(max(curves.con ./ 100, epsilon), 1 - epsilon);
    pIncon = min(max(curves.incon ./ 100, epsilon), 1 - epsilon);

    nLL = signedBX0ConditionNLL(panelData.successBaseline, panelData.nBaseline, pBaseline) + ...
        signedBX0ConditionNLL(panelData.successCon, panelData.nCon, pCon) + ...
        signedBX0ConditionNLL(panelData.successIncon, panelData.nIncon, pIncon);
end

function nLL = signedBX0ConditionNLL(successes, nTrials, probability)
    failures = nTrials - successes;
    nLL = -sum(successes .* log(probability) + failures .* log(1 - probability));
end

function curves = predictSignedBX0PanelFitCurvesForData(panelParams, globalDeltaX0, xBaseline, xCon, xIncon)
    validateSignedBX0BranchIndependence(panelParams, globalDeltaX0);
    baselineCurve = predictSignedBX0PanelBaseline(abs(xBaseline), panelParams);
    [conCurve, inconCurve] = predictSignedBX0PanelOpto(abs(xCon), abs(xIncon), ...
        panelParams, globalDeltaX0);
    curves = struct('baseline', baselineCurve, 'con', conCurve, 'incon', inconCurve);
end

function displayCurves = predictSignedBX0PanelFitCurves(xMagnitude, panelParams, globalDeltaX0, sourcePanel)
    validateSignedBX0BranchIndependence(panelParams, globalDeltaX0);
    xMagnitude = abs(xMagnitude);
    displayCurves.baseline = predictSignedBX0PanelBaseline(xMagnitude, panelParams);
    [displayCurves.con, displayCurves.incon] = predictSignedBX0PanelOpto( ...
        xMagnitude, xMagnitude, panelParams, globalDeltaX0);
    displayCurves.sourcePanel = sourcePanel;
    validateSignedBX0PanelCurves(displayCurves, size(xMagnitude), sourcePanel);
end

function yBaseline = predictSignedBX0PanelBaseline(c, panelParams)
    yBaseline = predictShiftedWeibullBranch(c, ...
        panelParams(1), 50, panelParams(2), panelParams(3), 0);
end

function [yCon, yIncon] = predictSignedBX0PanelOpto(cCon, cIncon, panelParams, globalDeltaX0)
    deltaB_panel = panelParams(10);
    B_con = 50 + deltaB_panel;
    B_incon = 50 - deltaB_panel;
    X0_con = -globalDeltaX0;
    X0_incon = +globalDeltaX0;
    yCon = predictShiftedWeibullBranch(abs(cCon), ...
        panelParams(4), B_con, panelParams(5), panelParams(6), X0_con);
    yIncon = predictShiftedWeibullBranch(abs(cIncon), ...
        panelParams(7), B_incon, panelParams(8), panelParams(9), X0_incon);
end

function diagnostics = getSignedBX0PanelSlopeDiagnostics(panelParams, globalDeltaX0)
    maxSlopePctPerContrast = 5.0;
    deltaB_panel = panelParams(10);
    B_con = 50 + deltaB_panel;
    B_incon = 50 - deltaB_panel;
    amplitudes = [(100 - panelParams(1)) - 50, ...
        (100 - panelParams(4)) - B_con, ...
        (100 - panelParams(7)) - B_incon];
    alphas = [panelParams(2), panelParams(5), panelParams(8)];
    betas = [panelParams(3), panelParams(6), panelParams(9)];
    slopeValues = nan(1, 3);
    if all(isfinite([amplitudes, alphas, betas, globalDeltaX0])) && ...
            all(amplitudes > 0) && all(alphas > 0) && all(betas > 1)
        slopeValues = getWeibullHalfMaxSlope(amplitudes, alphas, betas);
    end
    diagnostics.maxSlopeBaseline = slopeValues(1);
    diagnostics.maxSlopeCon = slopeValues(2);
    diagnostics.maxSlopeIncon = slopeValues(3);
    diagnostics.maxSlopeValues = slopeValues;
    diagnostics.maxSlopeOverall = max(slopeValues, [], 'omitnan');
    diagnostics.maxAllowedSlopePctPerContrast = maxSlopePctPerContrast;
    diagnostics.slopeConstraintActive = true;
    diagnostics.slopeCapActive = isfinite(diagnostics.maxSlopeOverall) && ...
        diagnostics.maxSlopeOverall > maxSlopePctPerContrast;
    diagnostics.isValid = all(isfinite(slopeValues)) && ...
        all(slopeValues <= maxSlopePctPerContrast);
    diagnostics.slopeExcess = max(0, slopeValues - maxSlopePctPerContrast);
end

function diagnostics = makeSignedBX0PanelFitDiagnostics(params, bestStartIndex, exitFlag, lb, ub, globalDeltaX0)
    parameterNames = signedBX0PanelParameterNames();
    tol = max(1e-6, 1e-3 .* (ub - lb));
    hit = abs(params - lb) <= tol | abs(params - ub) <= tol;
    slopeDiagnostics = getSignedBX0PanelSlopeDiagnostics(params, globalDeltaX0);
    diagnostics = struct( ...
        'bestStartIndex', bestStartIndex, ...
        'exitFlag', exitFlag, ...
        'boundHit', any(hit), ...
        'boundHitFields', {parameterNames(hit)}, ...
        'maxSlopeOverall', slopeDiagnostics.maxSlopeOverall, ...
        'slopeCapActive', slopeDiagnostics.slopeCapActive);
end

function assertSignedBX0PanelFit(fitResult, expectedPanel, globalDeltaX0)
    actualPanel = '<missing>';
    if isfield(fitResult, 'sourcePanel')
        actualPanel = fitResult.sourcePanel;
    end
    if ~strcmp(actualPanel, expectedPanel)
        error('plotNakaRushtonFit5:SignedBX0PanelMismatch', ...
            'Expected %s panel fit, got %s.', expectedPanel, actualPanel);
    end
    if numel(fitResult.fitParams) ~= 10 || any(~isfinite(fitResult.fitParams))
        error('plotNakaRushtonFit5:SignedBX0PanelParamCount', ...
            'Panel %s must have exactly 10 finite fitted parameters.', expectedPanel);
    end
    if fitResult.globalDeltaX0 ~= globalDeltaX0
        error('plotNakaRushtonFit5:SignedBX0PanelX0Mismatch', ...
            'Panel %s did not hold the global deltaX0 fixed.', expectedPanel);
    end
end

function assertSignedBX0PanelCurveSource(displayCurves, expectedPanel)
    actualPanel = '<missing>';
    if isfield(displayCurves, 'sourcePanel')
        actualPanel = displayCurves.sourcePanel;
    end
    if ~strcmp(actualPanel, expectedPanel)
        error('plotNakaRushtonFit5:SignedBX0CurveSourceMismatch', ...
            'Expected %s curves, got %s.', expectedPanel, actualPanel);
    end
end

function validateSignedBX0BranchIndependence(panelParams, globalDeltaX0)
    c = [0 10 25 50 75 100];
    [con0, incon0] = predictSignedBX0PanelOpto(c, c, panelParams, globalDeltaX0);
    conPerturbed = panelParams;
    conPerturbed(4:6) = min(conPerturbed(4:6) + [0.5 1 0.1], [24.5 59 7.9]);
    [~, inconAfterConChange] = predictSignedBX0PanelOpto(c, c, conPerturbed, globalDeltaX0);
    inconPerturbed = panelParams;
    inconPerturbed(7:9) = min(inconPerturbed(7:9) + [0.5 1 0.1], [24.5 59 7.9]);
    [conAfterInconChange, ~] = predictSignedBX0PanelOpto(c, c, inconPerturbed, globalDeltaX0);
    if any(abs(incon0 - inconAfterConChange) > 1e-10) || ...
            any(abs(con0 - conAfterInconChange) > 1e-10)
        error('plotNakaRushtonFit5:SignedBX0BranchSwitching', ...
            'Displayed con/incon panel predictions are not branch-independent.');
    end
end
function validateSignedBX0PanelCurves(displayCurves, expectedSize, sourcePanel)
    conditionNames = {'baseline', 'con', 'incon'};
    for conditionIdx = 1:numel(conditionNames)
        conditionName = conditionNames{conditionIdx};
        y = displayCurves.(conditionName);
        if ~isequal(size(y), expectedSize) || ~isreal(y) || ...
                any(~isfinite(y(:))) || any(y(:) < 0) || any(y(:) > 100)
            error('plotNakaRushtonFit5:InvalidSignedBX0PanelCurve', ...
                ['Invalid weibullSignedBX0 %s panel curve for %s. ' ...
                'Expected size %s, got %s.'], ...
                sourcePanel, conditionName, mat2str(expectedSize), mat2str(size(y)));
        end
    end
end

function mdl = storeSignedBX0PanelFit(mdl, block, panelName, fitResult)
    panel = mdl.signedBX0.panelFits.(panelName);
    panel.fitParams(block,:) = fitResult.fitParams;
    panel.nLL(block) = fitResult.nLL;
    panel.fitStatus{block} = fitResult.fitStatus;
    panel.bestStartIndex(block) = fitResult.bestStartIndex;
    panel.exitFlag(block) = fitResult.exitFlag;
    panel.boundHit(block) = logical(fitResult.boundHit);
    panel.boundHitFields{block} = fitResult.boundHitFields;
    panel.maxSlopeOverall(block) = fitResult.maxSlopeOverall;
    panel.slopeCapActive(block) = logical(fitResult.slopeCapActive);
    panel.sourceChecksum(block) = fitResult.sourceChecksum;
    panel.fitChecksum(block) = fitResult.fitChecksum;
    panel.globalDeltaX0(block) = fitResult.globalDeltaX0;
    panel.parameterNames = fitResult.parameterNames;
    mdl.signedBX0.panelFits.(panelName) = panel;
end
function validateSignedBX0PanelFitFields(panelFitsHorizontal, panelFitsVertical, panelFitsMerged, block)
    panelFits.horizontal = panelFitsHorizontal;
    panelFits.vertical = panelFitsVertical;
    panelFits.merged = panelFitsMerged;
    names = {'horizontal', 'vertical', 'merged'};
    for ii = 1:numel(names)
        panelName = names{ii};
        params = panelFits.(panelName).fitParams(block, :);
        if numel(params) ~= 10 || any(~isfinite(params))
            error('plotNakaRushtonFit5:SignedBX0PanelStorageInvalid', ...
                'Stored %s panel fit for block %d is not a finite 10-parameter vector.', ...
                panelName, block);
        end
    end
end
function displayCurves = predictSignedBX0DisplayCurves(xMagnitude, fitParamRow)
    params = squeeze(fitParamRow);
    params = params(:)';
    if numel(params) < 11
        error('plotNakaRushtonFit5:SignedBX0ParamCount', ...
            'weibullSignedBX0 plotting requires 11 fitted parameters.');
    end
    params = params(1:11);
    c = abs(xMagnitude);

    ABase = params(1);
    alphaBase = params(2);
    betaBase = params(3);
    ACon = params(4);
    alphaCon = params(5);
    betaCon = params(6);
    AIncon = params(7);
    alphaIncon = params(8);
    betaIncon = params(9);
    deltaB = params(10);
    deltaX0 = params(11);

    BBase = 50;
    X0Base = 0;
    BHorizontal = BBase - deltaB;
    BVertical = BBase + deltaB;
    X0Horizontal = deltaX0;
    X0Vertical = -deltaX0;

    pBLneg = weibullSignedBX0Mdl(-c, ABase, alphaBase, betaBase, ...
        ABase, alphaBase, betaBase, BBase, X0Base);
    pBLpos = weibullSignedBX0Mdl(+c, ABase, alphaBase, betaBase, ...
        ABase, alphaBase, betaBase, BBase, X0Base);

    pHneg = weibullSignedBX0Mdl(-c, ACon, alphaCon, betaCon, ...
        AIncon, alphaIncon, betaIncon, BHorizontal, X0Horizontal);
    pHpos = weibullSignedBX0Mdl(+c, ACon, alphaCon, betaCon, ...
        AIncon, alphaIncon, betaIncon, BHorizontal, X0Horizontal);

    pVneg = weibullSignedBX0Mdl(-c, AIncon, alphaIncon, betaIncon, ...
        ACon, alphaCon, betaCon, BVertical, X0Vertical);
    pVpos = weibullSignedBX0Mdl(+c, AIncon, alphaIncon, betaIncon, ...
        ACon, alphaCon, betaCon, BVertical, X0Vertical);

    displayCurves.horizontal.baseline = 100 - pBLneg;
    displayCurves.horizontal.con = 100 - pHneg;
    displayCurves.horizontal.incon = 100 - pVneg;

    displayCurves.vertical.baseline = pBLpos;
    displayCurves.vertical.con = pVpos;
    displayCurves.vertical.incon = pHpos;

    displayCurves.merged.baseline = 0.5 .* (...
        displayCurves.horizontal.baseline + ...
        displayCurves.vertical.baseline);
    displayCurves.merged.con = 0.5 .* (...
        displayCurves.horizontal.con + ...
        displayCurves.vertical.con);
    displayCurves.merged.incon = 0.5 .* (...
        displayCurves.horizontal.incon + ...
        displayCurves.vertical.incon);

    validateSignedX0DisplayCurves(displayCurves, size(xMagnitude));
end

function overlaySignedBX0SideCurves(ax, xPlot, viewCurves)
    axes(ax);
    hold(ax, 'on');
    plot(ax, xPlot, viewCurves.baseline, 'Color', [0 0 0], ...
        'LineWidth', 3, 'HandleVisibility', 'off');
    plot(ax, xPlot, viewCurves.con, ...
        'Color', [0.9294, 0.1098, 0.1373] * 1.05, ...
        'LineWidth', 3, 'HandleVisibility', 'off');
    plot(ax, xPlot, viewCurves.incon, ...
        'Color', [0, 0.0941, 0.6627] * 1.25, ...
        'LineWidth', 3, 'HandleVisibility', 'off');
end
function displayCurves = predictSignedX0DisplayCurves(xPlot, fitParamRow)
    params = squeeze(fitParamRow);
    params = params(:)';
    if numel(params) < 10
        error('plotNakaRushtonFit5:SignedX0ParamCount', ...
            'weibullSignedX0 plotting requires 10 fitted parameters.');
    end
    params = params(1:10);
    c = abs(xPlot);

    pBLneg = weibullSignedX0Mdl(-c, params(1), params(2), params(3), 0);
    pBLpos = weibullSignedX0Mdl(+c, params(1), params(2), params(3), 0);

    pHneg = weibullSignedX0Mdl(-c, params(4), params(5), params(6), params(10));
    pHpos = weibullSignedX0Mdl(+c, params(4), params(5), params(6), params(10));

    pVneg = weibullSignedX0Mdl(-c, params(7), params(8), params(9), -params(10));
    pVpos = weibullSignedX0Mdl(+c, params(7), params(8), params(9), -params(10));

    displayCurves.horizontal.baseline = 100 - pBLneg;
    displayCurves.horizontal.con = 100 - pHneg;
    displayCurves.horizontal.incon = 100 - pVneg;

    displayCurves.vertical.baseline = pBLpos;
    displayCurves.vertical.con = pVpos;
    displayCurves.vertical.incon = pHpos;

    displayCurves.merged.baseline = 0.5 .* (...
        displayCurves.horizontal.baseline + ...
        displayCurves.vertical.baseline);
    displayCurves.merged.con = 0.5 .* (...
        displayCurves.horizontal.con + ...
        displayCurves.vertical.con);
    displayCurves.merged.incon = 0.5 .* (...
        displayCurves.horizontal.incon + ...
        displayCurves.vertical.incon);

    validateSignedX0DisplayCurves(displayCurves, size(xPlot));
end

function validateSignedX0DisplayCurves(displayCurves, expectedSize)
    viewNames = {'horizontal', 'vertical', 'merged'};
    conditionNames = {'baseline', 'con', 'incon'};
    for viewIdx = 1:numel(viewNames)
        viewName = viewNames{viewIdx};
        for conditionIdx = 1:numel(conditionNames)
            conditionName = conditionNames{conditionIdx};
            y = displayCurves.(viewName).(conditionName);
            if ~isequal(size(y), expectedSize) || ~isreal(y) || ...
                    any(~isfinite(y(:))) || any(y(:) < 0) || any(y(:) > 100)
                error('plotNakaRushtonFit5:InvalidSignedX0DisplayCurve', ...
                    ['Invalid weibullSignedX0 display curve for %s %s. ' ...
                    'Expected size %s, got %s.'], ...
                    viewName, conditionName, mat2str(expectedSize), ...
                    mat2str(size(y)));
            end
        end
    end
end
function addFitParameterTable(ax, headers, fitParamRow, modelTypeStr)
    if nargin < 4 || isempty(modelTypeStr)
        modelTypeStr = 'weibullfreeAll';
    end
    if isempty(headers) || isempty(fitParamRow)
        return;
    end
    if strcmp(modelTypeStr, 'weibullSignedBX0')
        addSignedBX0FitParameterTable(ax, fitParamRow);
        return;
    elseif strcmp(modelTypeStr, 'weibullSignedX0')
        addSignedX0FitParameterTable(ax, fitParamRow);
        return;
    end

    headers = headers(:)';
    fitParamRow = squeeze(fitParamRow);
    fitParamRow = fitParamRow(:)';

    nParams = min(numel(headers), numel(fitParamRow));
    headers = headers(1:nParams);
    fitParamRow = fitParamRow(1:nParams);

    paramHeaders = {'A', 'B', '\alpha', '\beta'};
    paramDisplayHeaders = {'A', 'B', '\alpha', '\beta'};
    rowLabels = {'Baseline', 'Con-Opto', 'Incon-Opto'};
    rowColors = [0 0 0; 0.55 0 0; 0 0.05 0.45];
    tableValues = nan(numel(rowLabels), numel(paramHeaders));

    for ii = 1:nParams
        header = headers{ii};
        if startsWith(header, 'AUC') || contains(header, 'AICc')
            continue;
        end

        [rowIdx, paramIdx] = parseFitParameterHeader(header, paramHeaders);
        if ~isnan(rowIdx) && ~isnan(paramIdx)
            tableValues(rowIdx, paramIdx) = fitParamRow(ii);
        end
    end

    if all(isnan(tableValues(:)))
        return;
    end

    % Opto parameters are fitted as deltas from baseline. B is symmetric:
    % B_con = 0.5 + deltaB_con and B_incon = 0.5 - deltaB_con.
    baselineValues = tableValues(1,:);
    conDeltas = tableValues(2,:);
    inconDeltas = tableValues(3,:);
    bIdx = strcmp(paramHeaders, 'B');
    nonBIdx = ~bIdx;

    tableValues(1, bIdx) = 0.5;
    tableValues(2, nonBIdx) = baselineValues(nonBIdx) + conDeltas(nonBIdx);
    tableValues(3, nonBIdx) = baselineValues(nonBIdx) + inconDeltas(nonBIdx);
    tableValues(2, bIdx) = 0.5 + conDeltas(bIdx);
    tableValues(3, bIdx) = 0.5 - conDeltas(bIdx);

    axPos = get(ax, 'Position');
    tableGap = 0.010;
    maxTableRight = 0.992;
    tableX = axPos(1) + axPos(3) + tableGap;
    availableWidth = maxTableRight - tableX;
    tableWidth = min(0.22, availableWidth);
    if tableWidth < 0.18
        tableWidth = 0.18;
        tableX = max(0.01, maxTableRight - tableWidth);
    end

    tableHeight = 0.44 * axPos(4);
    tableY = axPos(2) + 0.5 * axPos(4) - 0.5 * tableHeight;
    rowHeight = tableHeight / 4;
    labelWidth = 0.070;
    labelGap = 0.0015;
    paramGap = 0.0100;
    valueWidth = (tableWidth - labelWidth - labelGap - (numel(paramHeaders) - 1) * paramGap) / numel(paramHeaders);
    fontSize = 10.5;

    addFitTableCell(tableX, tableY + 3 * rowHeight, labelWidth, rowHeight, '', [0 0 0], fontSize, 'bold', 'left');
    for col = 1:numel(paramHeaders)
        addFitTableCell(tableX + labelWidth + labelGap + (col - 1) * (valueWidth + paramGap), ...
            tableY + 3 * rowHeight, valueWidth, rowHeight, paramDisplayHeaders{col}, ...
            [0 0 0], fontSize, 'bold', 'center');
    end

    for row = 1:numel(rowLabels)
        yPos = tableY + (3 - row) * rowHeight;
        addFitTableCell(tableX, yPos, labelWidth, rowHeight, rowLabels{row}, ...
            rowColors(row,:), fontSize, 'bold', 'left');
        for col = 1:numel(paramHeaders)
            valueStr = formatFitParameterValue(tableValues(row, col), paramHeaders{col});
            addFitTableCell(tableX + labelWidth + labelGap + (col - 1) * (valueWidth + paramGap), ...
                yPos, valueWidth, rowHeight, valueStr, rowColors(row,:), ...
                fontSize, 'normal', 'center');
        end
    end
end

function addSignedBX0FitParameterTable(ax, panelFitParams, globalDeltaX0, deltaAICcX0)
    if nargin < 4 || isempty(deltaAICcX0)
        deltaAICcX0 = NaN;
    end
    if nargin < 3 || isempty(globalDeltaX0) || ~isfinite(globalDeltaX0)
        return;
    end
    panelFitParams = squeeze(panelFitParams);
    panelFitParams = panelFitParams(:)';
    if numel(panelFitParams) < 10 || any(~isfinite(panelFitParams(1:10)))
        return;
    end

    params = panelFitParams(1:10);
    deltaB_panel = params(10);
    tableValues = [ ...
        params(1), 50, params(2), params(3), 0; ...
        params(4), 50 + deltaB_panel, params(5), params(6), -globalDeltaX0; ...
        params(7), 50 - deltaB_panel, params(8), params(9), +globalDeltaX0];
    paramHeaders = {'A', 'B', '\alpha', '\beta', 'X0'};
    rowLabels = {'Baseline', 'Con-Opto', 'Incon-Opto'};
    rowColors = [0 0 0; 0.55 0 0; 0 0.05 0.45];
    strongColor = [0 0 0];
    weakColor = [0.45 0.45 0.45];

    axPos = get(ax, 'Position');
    tableGap = 0.004;
    maxTableRight = 0.988;
    tableX = axPos(1) + axPos(3) + tableGap;
    availableWidth = maxTableRight - tableX;
    tableWidth = min(0.205, availableWidth);
    if tableWidth < 0.185
        tableWidth = 0.185;
        tableX = max(0.01, maxTableRight - tableWidth);
    end

    tableHeight = 0.44 * axPos(4);
    tableY = axPos(2) + 0.50 * axPos(4) - 0.5 * tableHeight;
    rowHeight = tableHeight / 5;
    fontSize = 10.5;
    columnX = [0.00, 0.43, 0.56, 0.69, 0.82, 0.94];
    columnWidth = [0.41, 0.105, 0.105, 0.105, 0.105, 0.055];
    columnX = tableX + tableWidth .* columnX;
    columnWidth = tableWidth .* columnWidth;
    if columnX(end) + columnWidth(end) > 0.99
        error('plotNakaRushtonFit5:SignedBX0TableClipped', ...
            'Signed-BX0 X0 table column exceeds the normalized figure boundary.');
    end

    addFitTableCell(columnX(1), tableY + 4 * rowHeight, columnWidth(1), ...
        rowHeight, '', [0 0 0], fontSize, 'bold', 'left');
    for col = 1:numel(paramHeaders)
        addFitTableCell(columnX(col + 1), tableY + 4 * rowHeight, ...
            columnWidth(col + 1), rowHeight, paramHeaders{col}, ...
            [0 0 0], fontSize, 'bold', 'center');
    end

    for row = 1:numel(rowLabels)
        yPos = tableY + (4 - row) * rowHeight;
        addFitTableCell(columnX(1), yPos, columnWidth(1), rowHeight, ...
            rowLabels{row}, rowColors(row,:), fontSize, 'bold', 'left');
        for col = 1:numel(paramHeaders)
            valueStr = formatSignedBX0ParameterValue( ...
                tableValues(row, col), paramHeaders{col});
            addFitTableCell(columnX(col + 1), yPos, columnWidth(col + 1), ...
                rowHeight, valueStr, rowColors(row,:), fontSize, 'normal', 'center');
        end
    end

    if isfinite(deltaAICcX0) && deltaAICcX0 >= 8
        aicColor = strongColor;
    else
        aicColor = weakColor;
    end
    aicText = sprintf('DeltaAICc_X0    %.1f', deltaAICcX0);
    addFitTableCell(tableX, tableY, tableWidth, rowHeight, ...
        aicText, aicColor, fontSize, 'normal', 'left');
end
function valueStr = formatSignedBX0ParameterValue(value, paramName)
    if isnan(value)
        valueStr = '';
        return;
    end
    if strcmp(paramName, 'X0')
        if abs(value) < 0.05
            valueStr = sprintf('%.1f', 0);
        else
            valueStr = sprintf('%+.1f', value);
        end
    else
        valueStr = sprintf('%.1f', value);
    end
end
function addSignedX0FitParameterTable(ax, fitParamRow)
    fitParamRow = squeeze(fitParamRow);
    fitParamRow = fitParamRow(:)';
    if numel(fitParamRow) < 10
        return;
    end

    deltaX0 = fitParamRow(10);
    tableValues = [ ...
        fitParamRow(1), fitParamRow(2), fitParamRow(3), 0; ...
        fitParamRow(4), fitParamRow(5), fitParamRow(6), deltaX0; ...
        fitParamRow(7), fitParamRow(8), fitParamRow(9), -deltaX0];
    paramHeaders = {'A', '\alpha', '\beta', 'X0'};
    rowLabels = {'Baseline', 'H-Opto', 'V-Opto'};
    rowColors = [0 0 0; 0.55 0 0; 0 0.05 0.45];

    axPos = get(ax, 'Position');
    tableGap = 0.010;
    maxTableRight = 0.992;
    tableX = axPos(1) + axPos(3) + tableGap;
    availableWidth = maxTableRight - tableX;
    tableWidth = min(0.22, availableWidth);
    if tableWidth < 0.18
        tableWidth = 0.18;
        tableX = max(0.01, maxTableRight - tableWidth);
    end

    tableHeight = 0.44 * axPos(4);
    tableY = axPos(2) + 0.5 * axPos(4) - 0.5 * tableHeight;
    rowHeight = tableHeight / 4;
    labelWidth = 0.060;
    labelGap = 0.002;
    paramGap = 0.010;
    valueWidth = (tableWidth - labelWidth - labelGap - ...
        (numel(paramHeaders) - 1) * paramGap) / numel(paramHeaders);
    fontSize = 10.5;

    addFitTableCell(tableX, tableY + 3 * rowHeight, labelWidth, ...
        rowHeight, '', [0 0 0], fontSize, 'bold', 'left');
    for col = 1:numel(paramHeaders)
        addFitTableCell(tableX + labelWidth + labelGap + ...
            (col - 1) * (valueWidth + paramGap), ...
            tableY + 3 * rowHeight, valueWidth, rowHeight, ...
            paramHeaders{col}, [0 0 0], fontSize, 'bold', 'center');
    end

    for row = 1:numel(rowLabels)
        yPos = tableY + (3 - row) * rowHeight;
        addFitTableCell(tableX, yPos, labelWidth, rowHeight, ...
            rowLabels{row}, rowColors(row,:), fontSize, 'bold', 'left');
        for col = 1:numel(paramHeaders)
            valueStr = formatSignedX0ParameterValue(...
                tableValues(row, col), paramHeaders{col});
            addFitTableCell(tableX + labelWidth + labelGap + ...
                (col - 1) * (valueWidth + paramGap), yPos, ...
                valueWidth, rowHeight, valueStr, rowColors(row,:), ...
                fontSize, 'normal', 'center');
        end
    end
end

function valueStr = formatSignedX0ParameterValue(value, paramName)
    if isnan(value)
        valueStr = '';
        return;
    end
    if strcmp(paramName, 'A')
        value = 100 .* value;
        valueStr = sprintf('%.1f', value);
    elseif strcmp(paramName, 'X0')
        valueStr = sprintf('%+.1f', value);
    else
        valueStr = sprintf('%.1f', value);
    end
end
function [rowIdx, paramIdx] = parseFitParameterHeader(header, paramHeaders)
    rowIdx = nan;
    paramIdx = nan;

    if contains(header, 'incon-bl')
        rowIdx = 3;
    elseif contains(header, 'con-bl')
        rowIdx = 2;
    elseif contains(header, '^{bl}')
        rowIdx = 1;
    end

    paramName = regexprep(header, '^\\Delta', '');
    paramName = regexprep(paramName, '\^\{[^}]+\}', '');
    paramName = strtrim(paramName);

    for ii = 1:numel(paramHeaders)
        if strcmp(paramName, paramHeaders{ii})
            paramIdx = ii;
            return;
        end
    end
end

function addFitTableCell(xPos, yPos, widthVal, heightVal, textVal, colorVal, fontSize, fontWeight, horizontalAlignment)
    xPos = max(0, min(0.999, xPos));
    yPos = max(0, min(0.999, yPos));
    widthVal = max(0.001, min(widthVal, 1 - xPos));
    heightVal = max(0.001, min(heightVal, 1 - yPos));

    annotation(gcf, 'textbox', [xPos, yPos, widthVal, heightVal], ...
        'String', textVal, ...
        'HorizontalAlignment', horizontalAlignment, ...
        'VerticalAlignment', 'middle', ...
        'FontSize', fontSize, ...
        'FontWeight', fontWeight, ...
        'Interpreter', 'tex', ...
        'Color', colorVal, ...
        'EdgeColor', 'none', ...
        'BackgroundColor', 'none', ...
        'FitBoxToText', 'off', ...
        'Tag', 'PlotNakaFitParamText');
end

function valueStr = formatFitParameterValue(value, paramName)
    if isnan(value)
        valueStr = '';
        return;
    end

    isPercentParam = any(strcmp(paramName, {'A', 'B'}));
    if isPercentParam
        value = value * 100;
    end

    valueStr = sprintf('%.1f', value);
end

function textOut = signedX0AnnotationText(mdl, block)
    textOut = '';
    if ~isfield(mdl, 'signedX0')
        return;
    end
    requiredFields = {'X0Horizontal', 'X0Vertical', ...
        'deltaAICcX0', 'akaikeWeightX0'};
    for idx = 1:numel(requiredFields)
        if ~isfield(mdl.signedX0, requiredFields{idx}) || ...
                numel(mdl.signedX0.(requiredFields{idx})) < block
            return;
        end
    end
    textOut = sprintf(['X0H: %+.1f%%\n' ...
        'X0V: %+.1f%%\n' ...
        'DeltaAICc_X0: %.2g\n' ...
        'w_X0: %.2f'], ...
        mdl.signedX0.X0Horizontal(block), ...
        mdl.signedX0.X0Vertical(block), ...
        mdl.signedX0.deltaAICcX0(block), ...
        mdl.signedX0.akaikeWeightX0(block));
end
function addBlockSuplabel(blockLabel)
    [~, hLabel] = suplabel(blockLabel, 't', [.1 .1 .82 .88]);
    set(hLabel, ...
        'FontSize', 16, ...
        'FontWeight', 'normal', ...
        'Interpreter', 'none', ...
        'Tag', 'PlotNakaBlockLabel');
end

function setPlotAnnotationFontSizes()
    set(findall(gcf, 'Tag', 'PlotNakaStatsText'), 'FontSize', 10.5);
    set(findall(gcf, 'Tag', 'PlotNakaDeltaText'), 'FontSize', 13.5);
    set(findall(gcf, 'Tag', 'PlotNakaFitParamText'), 'FontSize', 10.5);
    set(findall(gcf, 'Tag', 'PlotNakaBlockLabel'), 'FontSize', 16);
end

function blockInfo = getPlotBlockInfo(datastruct, analysisBlockID, clusterBlocks, block, plotAverageFlag)
    if plotAverageFlag == 1
        blockInfo.blockIdx = clusterBlocks(1);
    else
        blockInfo.blockIdx = clusterBlocks(block);
    end

    if numel(analysisBlockID) >= blockInfo.blockIdx
        blockInfo.datastructIdx = analysisBlockID(blockInfo.blockIdx);
    else
        blockInfo.datastructIdx = blockInfo.blockIdx;
    end

    blockInfo.date = valueToChar(datastruct(blockInfo.datastructIdx).date);
    blockInfo.run = valueToChar(datastruct(blockInfo.datastructIdx).run);
    blockInfo.label = [blockInfo.date 'R' blockInfo.run];
end

function labels = normalizeClusterLabels(labels, cluster, nBlocks)
    if nargin < 3 || isempty(nBlocks)
        nBlocks = 1;
    end

    defaultLabel = sprintf('Cluster: %s', valueToChar(cluster));

    if nargin < 1 || isempty(labels)
        labels = repmat({defaultLabel}, 1, nBlocks);
        return;
    end

    % Convert supported scalar/vector inputs into a row cell array.
    if ischar(labels)
        labels = {labels};
    elseif isstring(labels)
        labels = cellstr(labels(:)');
    elseif isnumeric(labels) || islogical(labels)
        labels = arrayfun(@(x) valueToChar(x), labels(:)', ...
            'UniformOutput', false);
    elseif iscell(labels)
        labels = labels(:)';
        for ii = 1:numel(labels)
            labels{ii} = valueToChar(labels{ii});
        end
    else
        labels = {valueToChar(labels)};
    end

    % One shared cluster label is valid for every page in this cluster.
    if numel(labels) == 1 && nBlocks > 1
        labels = repmat(labels, 1, nBlocks);
    end

    if numel(labels) ~= nBlocks
        error('plotNakaRushtonFit5:ClusterLabelCountMismatch', ...
            ['clusterLabel has %d entries, but %d render blocks were ' ...
            'requested. Supply either one shared label or one label per block.'], ...
            numel(labels), nBlocks);
    end
end

function label = getClusterLabelForRenderBlock(clusterLabel, block)
    if isempty(clusterLabel)
        label = '';
    elseif iscell(clusterLabel)
        label = clusterLabel{block};
    elseif numel(clusterLabel) >= block
        label = clusterLabel(block);
    else
        label = clusterLabel(1);
    end
end

function labelLine = formatClusterLabelLine(clusterLabel)
    if nargin < 1 || isempty(clusterLabel)
        labelLine = '';
        return;
    end

    labelLine = sprintf('%s\n', char(clusterLabel));
end

function str = valueToChar(value)
    if ischar(value)
        str = value;
    elseif isstring(value)
        str = char(value);
    elseif isnumeric(value)
        str = num2str(value);
    else
        str = char(string(value));
    end
end


function opts = applyDeltaPermutationPlotDefaults(opts)
    if ~isfield(opts, 'showDeltaPermutationStats') || isempty(opts.showDeltaPermutationStats)
        opts.showDeltaPermutationStats = false;
    end
    if ~isfield(opts, 'nDeltaPermutations') || isempty(opts.nDeltaPermutations)
        opts.nDeltaPermutations = 500;
    end
    if ~isfield(opts, 'deltaPermutationBaseSeed') || isempty(opts.deltaPermutationBaseSeed)
        opts.deltaPermutationBaseSeed = 99173;
    end
    if ~isfield(opts, 'figureVisible') || isempty(opts.figureVisible)
        opts.figureVisible = 'on';
    end
    if ~isfield(opts, 'closeAfterRender') || isempty(opts.closeAfterRender)
        opts.closeAfterRender = false;
    end
    opts.showDeltaPermutationStats = logical(opts.showDeltaPermutationStats);
    opts.closeAfterRender = logical(opts.closeAfterRender);
end

function seed = stableDeltaPermutationSeed(baseSeed, primaryID, secondaryID)
    values = double(char(string(primaryID)));
    if isnumeric(primaryID)
        values = [values, double(primaryID(:)')];
    end
    if nargin > 2 && ~isempty(secondaryID)
        if isnumeric(secondaryID)
            values = [values, double(secondaryID(:)')];
        else
            values = [values, double(char(string(secondaryID)))];
        end
    end
    seed = mod(double(baseSeed) + sum((1:numel(values)) .* values), 2^31 - 1);
    if seed <= 0
        seed = double(baseSeed);
    end
end

function assertPermutationMatchesDisplayedBias(permResult, displayedX, displayedY)
    displayedX = displayedX(:);
    displayedY = displayedY(:);
    keep = isfinite(displayedX) & isfinite(displayedY);
    displayedX = displayedX(keep);
    displayedY = displayedY(keep);
    [displayedX, order] = sort(displayedX);
    displayedY = displayedY(order);
    if numel(displayedX) ~= numel(permResult.contrast) || ...
            any(abs(displayedX - permResult.contrast(:)) > 1e-10) || ...
            any(abs(displayedY - permResult.observedDeltaBias(:)) > 1e-8)
        error('plotNakaRushtonFit5:DeltaPermutationMismatch', ...
            ['Permutation inputs do not reproduce the displayed merged ' ...
            'purple deltaBias points.']);
    end
end

function assertDeltaPermutationVisualizationAudit(audit, permResult, contextLabel)
    if ~isfield(audit, 'visualization') || isempty(audit.visualization)
        error('plotNakaRushtonFit5:MissingDeltaPermutationVisualizationAudit', ...
            'Missing permutation visualization audit for %s.', contextLabel);
    end

    if isfield(audit.visualization, 'nOverallNullBands')
        if audit.visualization.nNullIntervals ~= 0 || ...
                audit.visualization.nNullMedians ~= 0 || ...
                audit.visualization.nLabels ~= 0 || ...
                audit.visualization.nOverallNullBands ~= 1 || ...
                audit.visualization.nOverallNullMedians ~= 1 || ...
                audit.visualization.nOverall < 1
            error('plotNakaRushtonFit5:DeltaPermutationVisualizationCountMismatch', ...
                ['Permutation visualization count mismatch for %s: expected ' ...
                '0 contrast intervals/medians/labels and 1 overall band/median; ' ...
                'got intervals %d, medians %d, labels %d, bands %d, overall medians %d, overall %d.'], ...
                contextLabel, audit.visualization.nNullIntervals, ...
                audit.visualization.nNullMedians, audit.visualization.nLabels, ...
                audit.visualization.nOverallNullBands, ...
                audit.visualization.nOverallNullMedians, audit.visualization.nOverall);
        end
    else
        nExpected = numel(permResult.contrast);
        if audit.visualization.nNullIntervals ~= nExpected || ...
                audit.visualization.nNullMedians ~= nExpected || ...
                audit.visualization.nLabels ~= nExpected || ...
                audit.visualization.nOverall < 1
            error('plotNakaRushtonFit5:DeltaPermutationVisualizationCountMismatch', ...
                ['Permutation visualization count mismatch for %s: expected %d, ' ...
                'intervals %d, medians %d, labels %d, overall %d.'], ...
                contextLabel, nExpected, audit.visualization.nNullIntervals, ...
                audit.visualization.nNullMedians, audit.visualization.nLabels, ...
                audit.visualization.nOverall);
        end
        if ~audit.visualization.labelsShareY
            error('plotNakaRushtonFit5:DeltaPermutationLabelYMismatch', ...
                'Permutation labels do not share one y-coordinate for %s.', contextLabel);
        end
    end

    if isfield(audit.visualization, 'allHandlesHidden') && ~audit.visualization.allHandlesHidden
        error('DeltaPermutation:HandleVisibility', ...
            'Permutation visualization handles are not hidden from legend for %s.', char(string(contextLabel)));
    end
    fprintf('Added permutation visuals | contrast boxes %d | contrast labels %d | overall bands %d | overall annotations %d\n', ...
        audit.visualization.nNullIntervals, audit.visualization.nLabels, ...
        getAuditScalarField(audit.visualization, 'nOverallNullBands', 0), ...
        audit.visualization.nOverall);
end

function value = getAuditScalarField(s, fieldName, defaultValue)
    if isfield(s, fieldName)
        value = s.(fieldName);
    else
        value = defaultValue;
    end
end

function saveDeltaPermutationExampleFigure(fig, plotOpts, exampleType)
    if ~isfield(plotOpts, 'saveDeltaPermutationExamples') || ...
            ~plotOpts.saveDeltaPermutationExamples || ...
            ~isfield(plotOpts, 'deltaPermutationExampleDir') || ...
            isempty(plotOpts.deltaPermutationExampleDir)
        return;
    end
    if strcmp(exampleType, 'individual')
        fileName = 'PepperR_regularIndividual_permStats.png';
    else
        return;
    end
    outputPath = fullfile(plotOpts.deltaPermutationExampleDir, fileName);
    if ~exist(fileparts(outputPath), 'dir')
        mkdir(fileparts(outputPath));
    end
    if isfile(outputPath)
        return;
    end
    exportgraphics(fig, outputPath, 'Resolution', 200);
    info = dir(outputPath);
    assert(isfile(outputPath) && info.bytes > 0, ...
        'Failed to write delta permutation example PNG: %s', outputPath);
end

function validateDeltaPermutationLegend(ax, contextLabel)
    legends = findobj(ancestor(ax, 'figure'), 'Type', 'Legend');
    expected = {'Biasing', 'Masking'};
    legendStrings = {};
    for ii = 1:numel(legends)
        candidateStrings = cellstr(string(legends(ii).String));
        if numel(candidateStrings) == 2 && isequal(candidateStrings(:)', expected)
            legendStrings = candidateStrings;
            break;
        end
    end
    if isempty(legendStrings)
        allStrings = cell(size(legends));
        for ii = 1:numel(legends)
            allStrings{ii} = strjoin(cellstr(string(legends(ii).String)), ', ');
        end
        error('plotNakaRushtonFit5:DeltaPermutationLegendMismatch', ...
            'Legend mismatch for %s. Expected Biasing/Masking. Found legends: %s', ...
            contextLabel, strjoin(allStrings, ' | '));
    end
    fprintf('Legend validation passed for %s: {%s, %s}\n', ...
        contextLabel, legendStrings{1}, legendStrings{2});
end
function mdl = appendDeltaPermutationAudit(mdl, audit)
    if ~isfield(mdl, 'deltaBiasPermutationExperimentContrasts') || ...
            isempty(mdl.deltaBiasPermutationExperimentContrasts)
        mdl.deltaBiasPermutationExperimentContrasts = audit.experimentContrasts;
    else
        mdl.deltaBiasPermutationExperimentContrasts = [ ...
            mdl.deltaBiasPermutationExperimentContrasts; audit.experimentContrasts];
    end
    if ~isfield(mdl, 'deltaBiasPermutationExperimentSummary') || ...
            isempty(mdl.deltaBiasPermutationExperimentSummary)
        mdl.deltaBiasPermutationExperimentSummary = audit.experimentSummary;
    else
        mdl.deltaBiasPermutationExperimentSummary = [ ...
            mdl.deltaBiasPermutationExperimentSummary; audit.experimentSummary];
    end
end

function printDeltaPermutationSummary(labelType, label, nContrasts, meanDelta, pValue, significant, positiveOneSidedP)
    if significant
        sigLabel = '*';
    else
        sigLabel = 'n.s.';
    end
    fprintf('%s %s | n contrasts %d | mean DeltaBias %.3f | raw two-sided p %.4g | raw positive one-sided p %.4g | %s\n', ...
        labelType, char(string(label)), nContrasts, meanDelta, pValue, positiveOneSidedP, sigLabel);
end


function printDeltaPermutationContrastDiagnostics(labelType, label, permResult)
    fprintf('%s %s permutation contrast diagnostics:\n', ...
        labelType, char(string(label)));
    for contrastIdx = 1:numel(permResult.contrast)
        conErr = permResult.conReconstructionErrorPct(contrastIdx);
        inconErr = permResult.inconReconstructionErrorPct(contrastIdx);
        if abs(conErr) > 0.25 || abs(inconErr) > 0.25
            flag = ' | RECONSTRUCTION ERROR > 0.25 pct-pt';
        else
            flag = '';
        end
        fprintf(['  x %.4g | con %.3f%% incon %.3f%% | nCon %d nIncon %d | ' ...
            'conCorrect %d inconCorrect %d | conErr %.4f inconErr %.4f | ' ...
            'raw two-sided p %.4g | raw positive one-sided p %.4g%s\n'], ...
            permResult.contrast(contrastIdx), ...
            permResult.conPct(contrastIdx), permResult.inconPct(contrastIdx), ...
            permResult.nCon(contrastIdx), permResult.nIncon(contrastIdx), ...
            permResult.conCorrect(contrastIdx), permResult.inconCorrect(contrastIdx), ...
            conErr, inconErr, permResult.rawTwoSidedP(contrastIdx), ...
            permResult.rawPositiveOneSidedP(contrastIdx), flag);
    end
end
function [baseMean, conMean, inconMean] = computeMergedConditionMeans(mdl, block)
    baseMean = weightedMeanForPlot(mdl.xBlock(1,:,block), mdl.yBlock(1,:,block), 'merged');
    conMean = weightedMeanForPlot(mdl.xBlock(2,:,block), mdl.yBlock(2,:,block), 'merged');
    inconMean = weightedMeanForPlot(mdl.xBlock(3,:,block), mdl.yBlock(3,:,block), 'merged');
end

function [deltaBias, deltaMask] = computeDeltaFromConditionMeans(baseMean, conMean, inconMean)
    deltaBias = conMean - inconMean;
    optoMean = mean([conMean, inconMean], 'omitnan');
    deltaMask = baseMean - optoMean;
end

function yMean = weightedMeanForPlot(x, y, weightMode)
    x = x(:)';
    y = y(:)';

    validIdx = ~isnan(x) & ~isnan(y);
    if ~any(validIdx)
        yMean = NaN;
        return;
    end

    x = x(validIdx);
    y = y(validIdx);
    weights = getPlotMeanWeights(x, weightMode);

    yMean = sum(y .* weights, 'omitnan') / sum(weights, 'omitnan');
end

function weights = getPlotMeanWeights(x, weightMode)
    x = x(:)';
    weights = ones(size(x));

    switch weightMode
        case 'sideBaseline'
            % Split baseline panels: 0 contrast combines the two visual tags.
            weights(:) = 10;
            weights(abs(x) == 0) = 20;
        case 'sideOpto'
            % Split opto panels keep con/incon separated by visual tag.
            weights(:) = 10;
        case 'merged'
            % Merged panel: 0 contrast combines both visual tags/opto sides.
            weights(:) = 20;
            weights(abs(x) == 0) = 40;
        otherwise
            error('Unknown plot mean weight mode: %s', weightMode);
    end
end

function [xBias, yBias, xMask, yMask, audit] = computeSideDeltaData(sideData)
    [xBias, yBias, xMask, yMask, audit] = ...
        matchPsychometricDeltaPoints(sideData, ...
        struct('contrastPairTolerance', 5));
end

function [xOut, yOut] = collapseConditionData(xIn, yIn)
    [xIn, yIn] = cleanAndSortXY(xIn, yIn);

    if isempty(xIn)
        xOut = [];
        yOut = [];
        return;
    end

    [xOut, ~, groupIdx] = unique(xIn, 'stable');
    yOut = nan(size(xOut));

    for ii = 1:numel(xOut)
        yOut(ii) = nanmean(yIn(groupIdx == ii));
    end
end

function [xOut, yOut] = cleanAndSortXY(xIn, yIn)
    xIn = xIn(:)';
    yIn = yIn(:)';

    validIdx = ~isnan(xIn) & ~isnan(yIn);
    xOut = xIn(validIdx);
    yOut = yIn(validIdx);

    [xOut, sortIdx] = sort(xOut);
    yOut = yOut(sortIdx);
end

function h = plotConditionSeries(xPlot, yPlot, colorVal, markerType, markerFaceColor, markerSize, lineWidth, displayName, handleVisibility)
    [xPlot, yPlot] = cleanAndSortXY(xPlot, yPlot);

    if isempty(xPlot)
        h = plot(nan, nan, ...
            'LineStyle', 'none', ...
            'Color', colorVal, ...
            'LineWidth', lineWidth, ...
            'Marker', markerType, ...
            'MarkerFaceColor', markerFaceColor, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerSize', markerSize, ...
            'DisplayName', displayName, ...
            'HandleVisibility', handleVisibility);
        return;
    end

    h = plot(xPlot, yPlot, ...
        'LineStyle', 'none', ...
        'Color', colorVal, ...
        'LineWidth', lineWidth, ...
        'Marker', markerType, ...
        'MarkerFaceColor', markerFaceColor, ...
        'MarkerEdgeColor', 'k', ...
        'MarkerSize', markerSize, ...
        'DisplayName', displayName, ...
        'HandleVisibility', handleVisibility);
end

function style = getPreMergedPlotStyle()
    baselineColor = [0 0 0];
    conColor = min([0.9294, 0.1098, 0.1373] * 1.05, 1);
    inconColor = min([0, 0.0941, 0.6627] * 1.25, 1);

    desatAmount = 0.5;
    lighten = @(c) c * (1 - desatAmount) + [1 1 1] * desatAmount;

    style.baselineColor = lighten(baselineColor);
    style.baselineFaceColor = style.baselineColor;
    style.conColor = lighten(conColor);
    style.inconColor = lighten(inconColor);
end

function mdl = plotPreMergedPanel(mdl, block, xLimPre, meanBar, tickCfg)
    % Plot pre-merged/split data:
    %   x = signed Gabor contrast
    %   y = percent correct
    %
    % Also save panel-1 side means:
    %   baseline/con/incon x horizontal/vertical side.

    hold on;

    yline(50, '--', ...
        'LineWidth', 1.5, ...
        'Color', .4 * [1 1 1], ...
        'HandleVisibility', 'off');

    xline(0, '--', ...
        'LineWidth', 1.5, ...
        'Color', .4 * [1 1 1], ...
        'HandleVisibility', 'off');

    % Original-ish colors
    baselineColor = [0 0 0];
    conColor = [0.9294, 0.1098, 0.1373] * 1.05;
    inconColor = [0, 0.0941, 0.6627] * 1.25;

    conColor = min(conColor, 1);
    inconColor = min(inconColor, 1);

    % Lighten/desaturate for pre-merged panel
    desatAmount = 0.5;
    lighten = @(c) c * (1 - desatAmount) + [1 1 1] * desatAmount;

    baselineColorLight = lighten(baselineColor);
    conColorLight = lighten(conColor);
    inconColorLight = lighten(inconColor);

    markerSize = 12;
    lineWidth = 3;

    %% Baseline
    xBaseline = rmnan(mdl.xBaselinePreMerge(block, :));
    yBaseline = rmnan(mdl.yBaselinePreMerge(block, :));

    hBaseline = plot(xBaseline, yBaseline, ...
        'LineStyle', 'none', ...
        'Color', baselineColorLight, ...
        'LineWidth', lineWidth, ...
        'Marker', 'o', ...
        'MarkerFaceColor', baselineColorLight, ...
        'MarkerEdgeColor', 'k', ...
        'MarkerSize', markerSize, ...
        'DisplayName', 'Baseline');

    %% Pool horizontal + vertical opto by congruency
    xH = rmnan(mdl.xHorizontalOptoPreMerge(block, :));
    yH = rmnan(mdl.yHorizontalOptoPreMerge(block, :));
    cH = rmnan(mdl.congruencyHorizontalOptoPreMerge(block, :));

    xV = rmnan(mdl.xVerticalOptoPreMerge(block, :));
    yV = rmnan(mdl.yVerticalOptoPreMerge(block, :));
    cV = rmnan(mdl.congruencyVerticalOptoPreMerge(block, :));

    xOpto = [xH, xV];
    yOpto = [yH, yV];
    cOpto = [cH, cV];

    hCon = plotByCongruencyPooled( ...
        xOpto, yOpto, cOpto, 1, ...
        conColorLight, '^', markerSize, lineWidth, 'Con-Opto');

    hIncon = plotByCongruencyPooled( ...
        xOpto, yOpto, cOpto, -1, ...
        inconColorLight, 'v', markerSize, lineWidth, 'Incon-Opto');

    %% Compute and save panel-1 side means
    idxBaseH = xBaseline <= 0;
    idxBaseV = xBaseline >= 0;

    idxConH = cOpto == 1  & xOpto <= 0;
    idxConV = cOpto == 1  & xOpto >= 0;

    idxInconH = cOpto == -1 & xOpto <= 0;
    idxInconV = cOpto == -1 & xOpto >= 0;

    baseHMean  = computeMeanForPlot(xBaseline, yBaseline, idxBaseH);
    conHMean   = computeMeanForPlot(xOpto, yOpto, idxConH);
    inconHMean = computeMeanForPlot(xOpto, yOpto, idxInconH);

    baseVMean  = computeMeanForPlot(xBaseline, yBaseline, idxBaseV);
    conVMean   = computeMeanForPlot(xOpto, yOpto, idxConV);
    inconVMean = computeMeanForPlot(xOpto, yOpto, idxInconV);

    mdl.meanBaselineHorizontal(block) = baseHMean;
    mdl.meanConOptoHorizontal(block) = conHMean;
    mdl.meanInconOptoHorizontal(block) = inconHMean;

    mdl.meanBaselineVertical(block) = baseVMean;
    mdl.meanConOptoVertical(block) = conVMean;
    mdl.meanInconOptoVertical(block) = inconVMean;

    % Also save compact vector forms for convenience
    mdl.meanPanel1Horizontal(block,:) = [baseHMean, conHMean, inconHMean];
    mdl.meanPanel1Vertical(block,:)   = [baseVMean, conVMean, inconVMean];
    mdl.meanPanel1Headers = {'Baseline', 'ConOpto', 'InconOpto'};
    mdl.meanPanel1SideHeaders = {'Horizontal', 'Vertical'};

    %% Jitter overlapping mean bars for visibility
    % This only jitters the plotted bars. Saved mdl means remain unjittered.
    jitterStep = .5;

    horizontalMeans = [baseHMean, conHMean, inconHMean];
    verticalMeans   = [baseVMean, conVMean, inconVMean];

    horizontalMeansPlot = jitterOverlappingMeans(horizontalMeans, jitterStep);
    verticalMeansPlot   = jitterOverlappingMeans(verticalMeans, jitterStep);

    %% Mean bars for left and right visual-stimulus sides
    meanLineWidth = 3;

    % Left-side / horizontal-stimulus mean bars
    plotMeanBarAtY(horizontalMeansPlot(1), meanBar.leftEdgeRange, baselineColorLight, meanLineWidth);
    plotMeanBarAtY(horizontalMeansPlot(2), meanBar.leftEdgeRange, conColorLight, meanLineWidth);
    plotMeanBarAtY(horizontalMeansPlot(3), meanBar.leftEdgeRange, inconColorLight, meanLineWidth);

    % Right-side / vertical-stimulus mean bars
    plotMeanBarAtY(verticalMeansPlot(1), meanBar.rightEdgeRangePre, baselineColorLight, meanLineWidth);
    plotMeanBarAtY(verticalMeansPlot(2), meanBar.rightEdgeRangePre, conColorLight, meanLineWidth);
    plotMeanBarAtY(verticalMeansPlot(3), meanBar.rightEdgeRangePre, inconColorLight, meanLineWidth);

    %% Axes/labels
    xlim(xLimPre);
    ylim([0 100]);

    xticks(xLimPre(1):tickCfg.preMajor:xLimPre(2));
    addSkippedTicks(xLimPre(1), xLimPre(2), tickCfg.preSkip, 'x');
    addSkippedTicks(0, 100, 10, 'y');

    xlabel('Signed Gabor contrast (%)');
    ylabel('Correct (%)');

    title('Split data before con/incon merge');

    legend([hBaseline, hCon, hIncon], ...
        {'Baseline', 'Con-Opto', 'Incon-Opto'}, ...
        'Location', 'southeast', ...
        'NumColumns', 1, ...
        'FontSize', 24);

    axis square;
    upFontSize(32, 0.01);
end

function h = plotByCongruencyPooled(x, y, congr, congrValue, colorVal, markerType, markerSize, lineWidth, displayName)
    idx = congr == congrValue & ~isnan(x) & ~isnan(y);

    if ~any(idx)
        h = plot(nan, nan, ...
            'LineStyle', 'none', ...
            'Color', colorVal, ...
            'LineWidth', lineWidth, ...
            'Marker', markerType, ...
            'MarkerFaceColor', colorVal, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerSize', markerSize, ...
            'DisplayName', displayName);
        return;
    end

    xUse = x(idx);
    yUse = y(idx);

    % Split by visual side using signed x.
    % This prevents the line from connecting the two x=0 points
    % belonging to horizontal-tagged and vertical-tagged visual stimuli.
    %
    % Negative-side branch includes the first zero encountered.
    % Positive-side branch includes the second zero encountered.
    zeroIdx = find(xUse == 0);

    if numel(zeroIdx) >= 2
        idxNegBranch = xUse < 0;
        idxPosBranch = xUse > 0;

        % Assign duplicate zero points separately.
        idxNegBranch(zeroIdx(1)) = true;
        idxPosBranch(zeroIdx(2)) = true;

    elseif numel(zeroIdx) == 1
        idxNegBranch = xUse < 0;
        idxPosBranch = xUse > 0;

        % Single zero: attach to both branches only if both sides exist;
        % otherwise attach to the existing side.
        if any(xUse < 0) && any(xUse > 0)
            idxNegBranch(zeroIdx(1)) = true;
            idxPosBranch(zeroIdx(1)) = true;
        elseif any(xUse < 0)
            idxNegBranch(zeroIdx(1)) = true;
        else
            idxPosBranch(zeroIdx(1)) = true;
        end

    else
        idxNegBranch = xUse < 0;
        idxPosBranch = xUse > 0;
    end

    % Plot negative branch. This handle is used for the legend.
    h = plotOneBranch(xUse(idxNegBranch), yUse(idxNegBranch), ...
        colorVal, markerType, markerSize, lineWidth, displayName, 'on');

    % Plot positive branch, but hide from legend.
    plotOneBranch(xUse(idxPosBranch), yUse(idxPosBranch), ...
        colorVal, markerType, markerSize, lineWidth, displayName, 'off');
end


function h = plotOneBranch(xPlot, yPlot, colorVal, markerType, markerSize, lineWidth, displayName, handleVisibility)
    if isempty(xPlot)
        h = plot(nan, nan, ...
            'LineStyle', 'none', ...
            'Color', colorVal, ...
            'LineWidth', lineWidth, ...
            'Marker', markerType, ...
            'MarkerFaceColor', colorVal, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerSize', markerSize, ...
            'DisplayName', displayName, ...
            'HandleVisibility', handleVisibility);
        return;
    end

    [xPlot, sortIdx] = sort(xPlot);
    yPlot = yPlot(sortIdx);

    h = plot(xPlot, yPlot, ...
        'LineStyle', 'none', ...
        'Color', colorVal, ...
        'LineWidth', lineWidth, ...
        'Marker', markerType, ...
        'MarkerFaceColor', colorVal, ...
        'MarkerEdgeColor', 'k', ...
        'MarkerSize', markerSize, ...
        'DisplayName', displayName, ...
        'HandleVisibility', handleVisibility);
end

function [xLimPre, xLimMerged, meanBar, tickCfg] = getPsychometricAxisLimits(mdl)
    % Compute dynamic x-limits, ticks, and mean-bar positions.
    %
    % Main rule:
    %   If max contrast <= 50, use 50.
    %   If max contrast > 50 and <= 100, use 100.
    %   If max contrast > 100, round upward to nearest 50.
    %
    % Panel 1 uses [-xMax xMax].
    % Panels 2/3 use [0 xMax].

    %% Collect all x values
    xPreAll = [];

    preFields = { ...
        'xBaselinePreMerge', ...
        'xHorizontalOptoPreMerge', ...
        'xVerticalOptoPreMerge'};

    for ii = 1:numel(preFields)
        if isfield(mdl, preFields{ii})
            xPreAll = [xPreAll, rmnan(mdl.(preFields{ii})(:))'];
        end
    end

    xMergedAll = [];

    mergedFields = { ...
        'xBaseline', ...
        'xConOpto', ...
        'xInconOpto'};

    for ii = 1:numel(mergedFields)
        if isfield(mdl, mergedFields{ii})
            xMergedAll = [xMergedAll, rmnan(mdl.(mergedFields{ii})(:))'];
        end
    end

    maxAbsPre = max(abs(xPreAll), [], 'omitnan');
    maxMerged = max(xMergedAll, [], 'omitnan');

    maxData = max([maxAbsPre, maxMerged], [], 'omitnan');

    if isempty(maxData) || isnan(maxData) || maxData == 0
        maxData = 50;
    end

    %% Choose clean xMax
    if maxData <= 50
        xMax = 50;
    elseif maxData <= 100
        xMax = 100;
    else
        xMax = ceil(maxData / 50) * 50;
    end

    xLimPre = [-xMax, xMax];
    xLimMerged = [0, xMax];

    %% Tick density
    % 4x denser than previous version.
    %
    % For xMax = 100:
    %   Panel 1 ticks every 25: -100 -75 -50 -25 0 25 50 75 100
    %   Panels 2/3 ticks every 12.5: 0 12.5 25 ... 100
    %
    % For xMax = 50:
    %   Panel 1 ticks every 12.5
    %   Panels 2/3 ticks every 6.25

    tickCfg.preMajor = xMax / 4;
    tickCfg.mergedMajor = xMax / 8;

    % For addSkippedTicks, use same intervals as xticks.
    tickCfg.preSkip = tickCfg.preMajor;
    tickCfg.mergedSkip = tickCfg.mergedMajor;

    %% Mean-bar geometry
    % Half the previous visual length.
    %
    % Old merged bar width was 10% of panel 2/3 range.
    % New merged bar width is 5% of panel 2/3 range.
    barWidthMerged = 0.05 * range(xLimMerged);

    % Panel 1 is twice as wide, so scale bar width accordingly
    % to preserve visual length across panels.
    barWidthPre = barWidthMerged * (range(xLimPre) / range(xLimMerged));

    %% Panel 1 edge bars
    meanBar.leftEdgeRange = [xLimPre(1), xLimPre(1) + barWidthPre];
    meanBar.rightEdgeRangePre = [xLimPre(2) - barWidthPre, xLimPre(2)];

    %% Panels 2/3 right-edge bars
    meanBar.rightEdgeRange = [xLimMerged(2) - barWidthMerged, xLimMerged(2)];

    %% Save metadata
    meanBar.xMax = xMax;
    meanBar.barWidthMerged = barWidthMerged;
    meanBar.barWidthPre = barWidthPre;
end

function h = plotSideMeanBar(x, y, idx, xRange, colorVal, lineWidth)
    % Plot a short horizontal mean bar for a subset of points.
    % xRange = [xStart xEnd].
    % Hidden from legend.

    idx = idx & ~isnan(x) & ~isnan(y);

    if ~any(idx)
        h = plot(nan, nan, ...
            'LineStyle', 'none', ...
            'Color', colorVal, ...
            'LineWidth', lineWidth, ...
            'HandleVisibility', 'off');
        return;
    end

    yMean = mean(y(idx), 'omitnan');

    h = plot(xRange, [yMean yMean], '-', ...
        'Color', colorVal, ...
        'LineWidth', lineWidth, ...
        'HandleVisibility', 'off');
end
function forceFigureSansSerif(figHandle)
    % Force a sans-serif font before export.
    % This helps prevent PDF/export functions from falling back to serif fonts.

    if nargin < 1 || isempty(figHandle)
        figHandle = gcf;
    end

    fontName = 'Arial';

    set(findall(figHandle, '-property', 'FontName'), 'FontName', fontName);
    set(findall(figHandle, 'Type', 'axes'), 'FontName', fontName);
    set(findall(figHandle, 'Type', 'text'), 'FontName', fontName);
    set(findall(figHandle, 'Type', 'legend'), 'FontName', fontName);

    set(figHandle, 'Renderer', 'painters');
end

function yMean = computeMeanForPlot(x, y, idx)
    idx = idx & ~isnan(x) & ~isnan(y);

    if ~any(idx)
        yMean = NaN;
        return;
    end

    yMean = mean(y(idx), 'omitnan');
end
function yPlot = jitterOverlappingMeans(y, jitterStep)
    % Apply small vertical jitter only when mean bars are identical.
    % Saved values remain unjittered.
    %
    % Example:
    %   [75 75 75] -> [74.65 75 75.35] if jitterStep = 0.35
    %
    % Non-identical values are unchanged:
    %   [74.8 75.0 75.2] stays [74.8 75.0 75.2]

    if nargin < 2
        jitterStep = 0.35;
    end

    yPlot = y;
    validIdx = find(~isnan(y));

    if numel(validIdx) <= 1
        return;
    end

    % Tiny tolerance: only truly identical/effectively identical values jitter.
    tol = 1e-9;

    used = false(size(validIdx));

    for ii = 1:numel(validIdx)
        if used(ii)
            continue;
        end

        thisOriginalIdx = validIdx(ii);

        sameGroupLocal = abs(y(validIdx) - y(thisOriginalIdx)) <= tol;
        sameGroupLocal = sameGroupLocal & ~used;

        groupOriginalIdx = validIdx(sameGroupLocal);
        used(sameGroupLocal) = true;

        nGroup = numel(groupOriginalIdx);

        if nGroup > 1
            offsets = ((1:nGroup) - (nGroup + 1)/2) * jitterStep;
            yPlot(groupOriginalIdx) = y(groupOriginalIdx) + offsets;
        end
    end
end
function h = plotMeanBarAtY(yMean, xRange, colorVal, lineWidth)
    % Plot a short horizontal mean bar at yMean.
    % Hidden from legend.

    if isnan(yMean)
        h = plot(nan, nan, ...
            'LineStyle', 'none', ...
            'Color', colorVal, ...
            'LineWidth', lineWidth, ...
            'HandleVisibility', 'off');
        return;
    end

    h = plot(xRange, [yMean yMean], '-', ...
        'Color', colorVal, ...
        'LineWidth', lineWidth, ...
        'HandleVisibility', 'off');
end


