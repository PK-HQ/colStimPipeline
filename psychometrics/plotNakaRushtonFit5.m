function mdl=plotNakaRushtonFit5(behavioralData, bitmapData, datastruct, analysisBlockID,...
    mdl, fitParams, x, monkeyName, clusterBlocks, plotAverageFlag, plotLine,...
    saveFlag, cluster, modelTypeStr, savefilename)
    endIdx=size(mdl.headers,2);
    % Number of blocks and conditions
    [nConditions, ~, nBlocks] = size(behavioralData.gaborContrasts(:, :, clusterBlocks));
    if plotAverageFlag==1
        nBlocks=1;
    end
    mdl.cluster=cluster;
    mdl.clusterBlocksIdx=clusterBlocks;
    
    % Dynamic x-axis limits for this animal/chamber/cluster set.
    % These are computed from all currently available mdl rows, so every block
    % in this plotting call uses the same x-limits and same mean-bar geometry.
    [xLimPre, xLimMerged, meanBar, tickCfg] = getPsychometricAxisLimits(mdl);
    
    for block = 1:nBlocks
        % Init figure
        dat=[];
        make_it_tight = true;
        hmarg = .18;
        subplot = @(m,n,p) subtightplot(m, n, p, [0.1 0.1], [hmarg hmarg], [0.1 0.1]);
        if ~make_it_tight,  clear subplot;  end
       
        figure('Name', ['Block #', datastruct(clusterBlocks(block)).date]);
        
        % Panel 1: pre-merged raw/split data
        subplot(1,3,1)
        mdl = plotPreMergedPanel(mdl, block, xLimPre, meanBar, tickCfg);
        
        % Panel 2: merged fitted data
        subplot(1,3,2)
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
                markerSize = 16;
                patchSaturationVal=1;
                % Data points
                shadedErrorBar(mdl.xBlock(cond,:,block), mdl.yBlock(cond,:,block), semY', 'patchSaturation', patchSaturationVal, 'lineprops', ...
                               {'Color', lineColor, 'LineStyle', 'none', 'LineWidth', 3, 'Marker', markerType, ...
                                'MarkerFaceColor', markerFaceColor, 'MarkerEdgeColor', edgeColor, 'MarkerSize', markerSize}); hold on;
                % Average
                barLength = meanBar.rightEdgeRange;
                plot(barLength, repmat(nanmean(mdl.yBlock(cond,:,block)),1,numel(barLength)), '-', 'Color', lineColor, 'LineWidth', 3, 'HandleVisibility', 'off')

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
                    [combinedBLStr,mdl.combinedBL(block)]=checkMergedBaseline(datastruct,analysisBlockID,block);

                    if isfield(bitmapData, 'meanPowerDensityWithinROI_mWmm2') & ~isempty(bitmapData.meanPowerDensityWithinROI_mWmm2(:,:,clusterBlocks(block)))
                        if plotAverageFlag % for plotting the average across all blocks
                            bitmapSPD=squeeze(bitmapData.meanPowerDensityWithinROI_mWmm2(:,:,clusterBlocks));
                            bitmapColumns=bitmapData.nColumns(:,clusterBlocks)';
                        else
                            bitmapSPD=squeeze(bitmapData.meanPowerDensityWithinROI_mWmm2(:,:,clusterBlocks(block)));
                            bitmapColumns=bitmapData.nColumns(:,block)';
                        end
                    else
                        bitmapSPD=[0 0];
                    end
                    bitmapColumnhv=nanmean(bitmapColumns,1);
                    bitmapSPDhv=nanmean(bitmapSPD,1);
                    bitmapSPDmean=nanmean(bitmapSPD,'all');
                    bitmapSPDstd=nanstd(bitmapSPD,[],'all');
                    if plotAverageFlag==1
                        title({[modelTypeStr ', cluster ' num2str(cluster) ' average'],...
                            ['meanPowerDensityWithinROI_mWmm2: ' num2str(bitmapSPDhv(1),2) ' & ' num2str(bitmapSPDhv(2),2) ' mW (' num2str(bitmapSPDmean,2) ' \pm ' num2str(bitmapSPDstd,1) ' mW)',...
                            ', Columns: ' num2str(bitmapColumnhv(1),2) ' & ' num2str(bitmapColumnhv(2),2)]});
                    else
                        %{
                        title({[datastruct(clusterBlocks(block)).date 'R' datastruct(analysisBlockID(block)).run ' (' combinedBLStr ')'],...
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
                    
                    blockIdx = clusterBlocks(block);
                    
                    meanColumns = mean(bitmapData.nColumns(:, blockIdx), 'omitnan');
                    
                    meanProjectorPD = mean(bitmapData.projectorPowerDensity_mWmm2(:, blockIdx), 'omitnan');
                    meanAreaROI     = mean(bitmapData.areaFinalROI(:, blockIdx), 'omitnan');
                    meanAreaON      = mean(bitmapData.areaPixelsONWithinROI(:, blockIdx), 'omitnan');
                    meanSpatialDC   = mean(bitmapData.spatialDutyCycleWithinROI(:, blockIdx), 'omitnan') * 100;
                    meanTemporalDC  = mean(bitmapData.temporalDutyCycle(:, blockIdx), 'omitnan') * 100;
                    
                    meanROIPD       = mean(bitmapData.meanPowerDensityWithinROI_mWmm2(:, blockIdx), 'omitnan');
                    meanTotalPower  = mean(bitmapData.totalPowerToOnPixelsWithinROI_mW(:, blockIdx), 'omitnan');
                    
                    % Keep title short
                    title(sprintf('%sR%s', ...
                          datastruct(analysisBlockID(block)).date, ...
                          datastruct(analysisBlockID(block)).run), ...
                          'Interpreter', 'tex');
                    
                    % Bottom-centered annotation, normalized to axes
                    optoText = sprintf([ ...
                        '%.0f cols\n' ...
                        'PD_{DMD} %.2f mW/mm^2 | Area_{ROI} %.2f mm^2 | Area_{ON} %.2f mm^2\n' ...
                        'sDC %.1f%% | tDC %.1f%%\n' ...
                        'PD_{ROI} %.2f mW/mm^2 | P_{total} %.2f mW'], ...
                        meanColumns, ...
                        meanProjectorPD, meanAreaROI, meanAreaON, ...
                        meanSpatialDC, meanTemporalDC, ...
                        meanROIPD, meanTotalPower);
                    
                        text(0.0, -0.28, optoText, ...
                            'Units', 'normalized', ...
                            'HorizontalAlignment', 'left', ...
                            'VerticalAlignment', 'top', ...
                            'FontSize', 10, ...
                            'Interpreter', 'tex', ...
                            'Clipping', 'off');
                    1;
                end
            elseif cond==4
                subplot(1,3,3)
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
                markerSize = 20;
                patchSaturationVal=1;
                % Data points
                shadedErrorBar(mdl.xBlock(cond,:,block), mdl.yBlock(cond,:,block), semY', 'patchSaturation', patchSaturationVal, 'lineprops', ...
                               {'Color', lineColor, 'LineStyle', 'none', 'LineWidth', 3, 'Marker', markerType, ...
                                'MarkerFaceColor', markerFaceColor, 'MarkerEdgeColor', edgeColor, 'MarkerSize', markerSize}); hold on;
                % Average
                plot(barLength, repmat(nanmean(mdl.yBlock(cond,:,block)),1,numel(barLength)), '-', 'Color', lineColor, 'LineWidth', 3, 'HandleVisibility', 'off')

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
                ylim([-50 50])
                xticks(xLimMerged(1):tickCfg.mergedMajor:xLimMerged(2));
                addSkippedTicks(xLimMerged(1), xLimMerged(2), tickCfg.mergedSkip, 'x');
                addSkippedTicks(-50, 50, 10, 'y');
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
                markerSize = 20;
                patchSaturationVal=1;
                % Data points
                shadedErrorBar(mdl.xBlock(cond,:,block), mdl.yBlock(cond,:,block), semY', 'patchSaturation', patchSaturationVal, 'lineprops', ...
                               {'Color', lineColor, 'LineStyle', 'none', 'LineWidth', 3, 'Marker', markerType, ...
                                'MarkerFaceColor', markerFaceColor, 'MarkerEdgeColor', edgeColor, 'MarkerSize', markerSize}); hold on;
                % Average
                plot(barLength, repmat(nanmean(mdl.yBlock(cond,:,block)),1,numel(barLength)), '-', 'Color', lineColor, 'LineWidth', 3, 'HandleVisibility', 'off')

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
                ylim([-50 50])
                xticks(xLimMerged(1):tickCfg.mergedMajor:xLimMerged(2));
                addSkippedTicks(xLimMerged(1), xLimMerged(2), tickCfg.mergedSkip, 'x');
                addSkippedTicks(-50, 50, 10, 'y');
                yline(0,'--','LineWidth',1.5,'Color',.4*[1 1 1],'HandleVisibility','off'); hold on;
                
                % Explicit legend handles for third panel
                % Purple = biasing, gray = masking
                hBiasLegend = plot(nan, nan, 's-', ...
                    'Color', [127, 0, 255]/255, ...
                    'MarkerFaceColor', [127, 0, 255]/255, ...
                    'MarkerEdgeColor', 'k', ...
                    'LineWidth', 3, ...
                    'MarkerSize', 20);
                
                hMaskLegend = plot(nan, nan, 's-', ...
                    'Color', [125, 125, 125]/255, ...
                    'MarkerFaceColor', [125, 125, 125]/255, ...
                    'MarkerEdgeColor', 'k', ...
                    'LineWidth', 3, ...
                    'MarkerSize', 20);
                
                legend([hBiasLegend, hMaskLegend], ...
                    {'Biasing', 'Masking'}, ...
                    'Location', 'southeast', ...
                    'NumColumns', 1, ...
                    'FontSize', 32);
                
                axis square


                deltaBias=mean((rmnan(mdl.yBlock(2,:,block)))-rmnan(mdl.yBlock(3,:,block)));
                deltaMask=mean(rmnan(mdl.yBlock(1,:,block)))-mean(mean(rmnan(mdl.yBlock(2:3,:,block))));
                mdl.deltaBias(block)=deltaBias;
                mdl.deltaMask(block)=deltaMask;
                title(sprintf(['ΔBias_{con-incon} = %.1f%%\n',...
               'ΔMask_{base-opto} = %.1f%%'], ...
               deltaBias, deltaMask), ...
               'Interpreter', 'tex');
            end
        end
        
        % Customize the starting position and spacing
        startPos = [30, 30]; % Starting position in data coordinates
        xSpacing =5; % Horizontal spacing between columns
        ySpacing = 3; % Vertical spacing between rows

        % Call the function to create the table
        subplot(1, 3, 2); ax1=gca;
        %createCustomTable2(ax1, modelTypeStr, mdl.headers, mdl.fittedParams(block,:), startPos, xSpacing, ySpacing);
        
        if plotAverageFlag==1
            [nConditions, ~, nBlocks] = size(behavioralData.gaborContrasts(:, :, clusterBlocks));
            nBlocks=1;
            block=1;
        end
        upFontSize(21, .01);
        %Saving
        if saveFlag
            forceFigureSansSerif(gcf);
            savePDF(savefilename, monkeyName, 1, block, nBlocks)
            % Png/SVG
            %{
            monkey=datastruct(clusterBlocks(block)).monkey;
            date= datastruct(clusterBlocks(block)).date;
            run=datastruct(clusterBlocks(block)).run;
            if ispc
              mainPath='Y:/';
            elseif contains(getenv('HOSTNAME'),'psy.utexas.edu')
              mainPath='/eslab/data/';
            end
            figPath=[mainPath monkey '\Meta\psychometrics\' datastruct(clusterBlocks(block)).chamber '-chamber\' modelTypeStr];
            figName=['\C' num2str(cluster) 'M' datastruct(clusterBlocks(block)).monkeyNo 'D' date 'R' run];
            set(findall(gcf, '-property', 'FontName'), 'FontName', 'SansSerif');                
            set(gcf, 'Renderer', 'painters'); % Use painters for vector graphics
            %print(gcf, [figPath '\png' figName '.png'], '-dpng', '-r600'); % High-res PNG
            %savefig(gcf, [figPath '\fig' figName '.fig']);           % FIG
            %print(gcf, [figPath '\svg' figName '.svg'], '-dsvg');        % SVG
            %}
        end
    end
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

    markerSize = 16;
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