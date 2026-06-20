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
    if saveFlag && ~plotAverageFlag && isempty(clusterLabel)
        error('plotNakaRushtonFit5:MissingClusterLabels', ...
            ['Saved individual-page render requested without cluster labels. ' ...
            'Expected one label per rendered block.']);
    end
    if ~plotAverageFlag && ~isempty(clusterLabel) && ...
            numel(clusterLabel) ~= nBlocks
        error('plotNakaRushtonFit5:ClusterLabelCountMismatch', ...
            ['clusterLabel has %d entries, but %d render blocks were ' ...
            'requested.'], numel(clusterLabel), nBlocks);
    end
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
       
        figure('Name', ['Block ', blockInfo.label]);
        
        sideData = getPreMergedSideData(mdl, block);

        % Row 1, columns 1-2: split pre-merged data by visual stimulus side
        axRow1Col1 = subplot(2,3,1);
        [mdl, sideData.horizontal] = plotSidePsychometricPanel(mdl, block, sideData.horizontal, ...
            xLimMerged, meanBar, tickCfg, 'Horizontal visual stimulus');

        subplot(2,3,2)
        [mdl, sideData.vertical] = plotSidePsychometricPanel(mdl, block, sideData.vertical, ...
            xLimMerged, meanBar, tickCfg, 'Vertical visual stimulus');

        % Row 2, columns 1-2: deltas from the side-specific split data
        subplot(2,3,4)
        mdl = plotSideDeltaPanel(mdl, block, sideData.horizontal, ...
            xLimMerged, meanBar, tickCfg, 'Horizontal visual stimulus', ...
            'Horizontal', blockInfo, baselineModeThis);

        subplot(2,3,5)
        mdl = plotSideDeltaPanel(mdl, block, sideData.vertical, ...
            xLimMerged, meanBar, tickCfg, 'Vertical visual stimulus', ...
            'Vertical', blockInfo, baselineModeThis);
        
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
                    predictedCurve = mdl.mdlBaseline(xPlot, fitParams(block, 1:end-1)).pcntrl; 
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
                    predictedCurve = mdl.mdlOpto(xPlot, fitParams(block, 1:end-1)).pc; 
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
                    predictedCurve = mdl.mdlOpto(xPlot, fitParams(block, 1:end-1)).pic; 

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
                        'BL: %s\n' ...
                        '%s' ...
                        '%.0f cols\n' ...
                        'PD_{DMD} %.2f mW mm^{-2}\n' ...
                        'Area_{ROI} %.2f mm^2\n' ...
                        'Area_{ON} %.2f mm^2\n' ...
                        'sDC %.1f%% | tDC %.1f%%\n' ...
                        'PD_{ROI} %s mW mm^{-2}\n' ...
                        'P_{total} %s mW'], ...
                        char(baselineModeThis), ...
                        clusterLabelLine, ...
                        meanColumns, ...
                        meanProjectorPD, ...
                        meanAreaROI, ...
                        meanAreaON, ...
                        meanSpatialDC, meanTemporalDC, ...
                        powerDisplay.PDROI, ...
                        powerDisplay.Ptotal);

                    optoStatsTextHandle = addOptoStatsText(axRow1Col1, optoText);

                    addFitParameterTable(axMergedPanel, mdl.headers, fitParams(block,:));
                    
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
        %Saving
        if saveFlag
            forceFigureSansSerif(gcf);
            if isempty(reportState)
                saveCompressedPDFPage(savefilename, monkeyName, gcf, appendPage);
            else
                reportState = stageReportPDFPage(reportState, gcf);
            end
            close(gcf);
            % Png/SVG
            %{
            monkey=datastruct(blockInfo.datastructIdx).monkey;
            date= blockInfo.date;
            run=blockInfo.run;
            if ispc
              mainPath='Y:/';
            elseif contains(getenv('HOSTNAME'),'psy.utexas.edu')
              mainPath='/eslab/data/';
            end
            figPath=[mainPath monkey '\Meta\psychometrics\' datastruct(blockInfo.datastructIdx).chamber '-chamber\' modelTypeStr];
            figName=['\C' num2str(cluster) 'M' datastruct(blockInfo.datastructIdx).monkeyNo 'D' date 'R' run];
            set(findall(gcf, '-property', 'FontName'), 'FontName', 'SansSerif');                
            set(gcf, 'Renderer', 'painters'); % Use painters for vector graphics
            %print(gcf, [figPath '\png' figName '.png'], '-dpng', '-r600'); % High-res PNG
            %savefig(gcf, [figPath '\fig' figName '.fig']);           % FIG
            %print(gcf, [figPath '\svg' figName '.svg'], '-dsvg');        % SVG
            %}
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

function mdl = plotSideDeltaPanel(mdl, block, sideData, xLimMerged, meanBar, tickCfg, titleStr, sideFieldName, blockInfo, baselineMode)
    hold on;

    lineWidth = 3;
    markerSize = 15;
    biasColor = [127, 0, 255] / 255;
    maskColor = [125, 125, 125] / 255;

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

function addFitParameterTable(ax, headers, fitParamRow)
    if isempty(headers) || isempty(fitParamRow)
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
    opts.showDeltaPermutationStats = logical(opts.showDeltaPermutationStats);
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
    if isfield(audit.visualization, 'allHandlesHidden') && ~audit.visualization.allHandlesHidden
        error('DeltaPermutation:HandleVisibility', ...
            'Permutation visualization handles are not hidden from legend for %s.', char(string(contextLabel)));
    end
    fprintf('Added null intervals %d | labels %d | overall annotations %d\n', ...
        audit.visualization.nNullIntervals, audit.visualization.nLabels, ...
        audit.visualization.nOverall);
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
        h = plot(nan, nan, '-', ...
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

    h = plot(xPlot, yPlot, '-', ...
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

    hBaseline = plot(xBaseline, yBaseline, '-', ...
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
        h = plot(nan, nan, '-', ...
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
        h = plot(nan, nan, '-', ...
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

    h = plot(xPlot, yPlot, '-', ...
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
        h = plot(nan, nan, '-', ...
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
        h = plot(nan, nan, '-', ...
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


