function [figureHandles, distributionData, sourceAudit, statsAudit] = ...
        plotPowerClusterMeanDistributions(mdl, aggregateFits, opts)
% Plot experiment-level means and effects for each aggregate power cluster.

    if nargin < 3 || isempty(opts)
        opts = struct();
    end
    opts = applyDefaults(opts);

    viewNames = {'Horizontal', 'Vertical', 'Merged'};
    viewTitles = {'Horizontal visual stimulus', ...
        'Vertical visual stimulus', 'Merged'};
    conditionLabels = {'Baseline', 'Con-Opto', 'Incon-Opto'};
    conditionColors = [0 0 0; 0.9294 0.1098 0.1373; 0 0.0941 0.6627];
    deltaLabels = {'\DeltaBias', '\DeltaMask'};
    deltaColors = [127 0 255; 125 125 125] ./ 255;
    conditionPairs = [1 2; 2 3; 1 3];
    deltaPairs = [1 1; 2 2; 1 2];
    conditionStatsOpts = performanceStatsOptions(opts);
    deltaStatsOpts = opts.statsOpts;
    deltaStatsOpts.alpha = opts.alpha;

    figureHandles = gobjects(numel(aggregateFits), 1);
    distributionData = repmat(struct(), numel(aggregateFits), 1);
    sourceAudit = emptyDistributionSourceAuditTable();
    statsAudit = emptyDistributionStatsAuditTable();
    if strcmp(opts.statsOpts.testDirection, 'one-sided')
        fprintf(['Distribution statistics: planned one-sided mode is ON. ' ...
            'Holm correction remains within each three-test panel.\n']);
    end
    for aggregateIdx = 1:numel(aggregateFits)
        rows = validateRows(aggregateFits(aggregateIdx).mdlRowIndices, mdl);
        [sourceBlocks, experimentIDs, baselineModes] = auditSessionMetadata( ...
            aggregateFits(aggregateIdx), rows, opts);
        nExperiments = numel(rows);
        clusterTitle = distributionTitle(aggregateFits(aggregateIdx), nExperiments);
        fprintf('\n[%s] %d experiments\n', clusterTitle, nExperiments);

        figureHandles(aggregateIdx) = figure( ...
            'Name', [clusterTitle ' experiment distributions'], ...
            'Color', 'w', ...
            'Visible', opts.figureVisible);
        makeSubplot = @(position) subtightplot(2, 3, position, ...
            [0.13 0.07], [0.23 0.10], [0.10 0.07]);

        distributionData(aggregateIdx).clusterID = ...
            aggregateFits(aggregateIdx).clusterID;
        distributionData(aggregateIdx).mdlRowIndices = rows;
        distributionData(aggregateIdx).nExperiments = nExperiments;
        topAxes = gobjects(1, 3);
        bottomAxes = gobjects(1, 3);
        conditionValuesByView = cell(1, 3);
        deltaValuesByView = cell(1, 3);
        conditionStatsByView = cell(1, 3);
        deltaStatsByView = cell(1, 3);
        conditionTestsByView = cell(1, 3);
        deltaTestsByView = cell(1, 3);

        for viewIdx = 1:3
            meanField = ['meanPsychometric' viewNames{viewIdx}];
            deltaField = ['meanDelta' viewNames{viewIdx}];

            % These are the values saved for each fitted experiment by
            % plotNakaRushtonFit5. They are not recomputed from pooled
            % aggregate curves or aggregate-fit parameters.
            conditionValues = readSessionValues( ...
                mdl, meanField, rows, 3, clusterTitle);
            deltaValues = readSessionValues( ...
                mdl, deltaField, rows, 2, clusterTitle);

            viewField = lower(viewNames{viewIdx});
            distributionData(aggregateIdx).(viewField).conditionSource = meanField;
            distributionData(aggregateIdx).(viewField).deltaSource = deltaField;
            distributionData(aggregateIdx).(viewField).conditions = conditionValues;
            distributionData(aggregateIdx).(viewField).deltas = deltaValues;

            conditionStats = summarizeColumns(conditionValues);
            deltaStats = summarizeColumns(deltaValues);
            conditionMetrics = {'baseline', 'con', 'incon'};
            deltaMetrics = {'deltaBias', 'deltaMask'};
            conditionColumns = { ...
                'column 1 (baseline)', 'column 2 (con)', ...
                'column 3 (incon)'};
            deltaColumns = { ...
                'column 1 (deltaBias)', 'column 2 (deltaMask)'};
            conditionAudit = buildDistributionSourceAuditRows( ...
                aggregateFits(aggregateIdx), 'experimentDistribution', ...
                lower(viewNames{viewIdx}), 'performance', conditionMetrics, ...
                rows, sourceBlocks, experimentIDs, ['mdl.' meanField], ...
                conditionColumns, conditionValues, ...
                isfinite(conditionValues), strings(size(conditionValues)), ...
                baselineModes);
            deltaAudit = buildDistributionSourceAuditRows( ...
                aggregateFits(aggregateIdx), 'deltaDistribution', ...
                lower(viewNames{viewIdx}), 'delta', deltaMetrics, ...
                rows, sourceBlocks, experimentIDs, ['mdl.' deltaField], ...
                deltaColumns, deltaValues, isfinite(deltaValues), ...
                strings(size(deltaValues)), baselineModes);
            assertAuditCounts(conditionAudit, conditionStats.n, ...
                clusterTitle, viewTitles{viewIdx}, conditionMetrics);
            assertAuditCounts(deltaAudit, deltaStats.n, ...
                clusterTitle, [viewTitles{viewIdx} ' deltas'], deltaMetrics);
            sourceAudit = [sourceAudit; conditionAudit; deltaAudit]; %#ok<AGROW>
            printSourceAuditSummary(conditionAudit, clusterTitle);
            printSourceAuditSummary(deltaAudit, clusterTitle);
            conditionValuesByView{viewIdx} = conditionValues;
            deltaValuesByView{viewIdx} = deltaValues;
            conditionStatsByView{viewIdx} = conditionStats;
            deltaStatsByView{viewIdx} = deltaStats;
            distributionData(aggregateIdx).(viewField).conditionStats = conditionStats;
            distributionData(aggregateIdx).(viewField).deltaStats = deltaStats;

            verifyPointCounts(conditionStats.n, nExperiments, ...
                clusterTitle, viewTitles{viewIdx}, conditionLabels);
            verifyPointCounts(deltaStats.n, nExperiments, ...
                clusterTitle, viewTitles{viewIdx}, deltaLabels);
            printMeans(viewTitles{viewIdx}, conditionStats, deltaStats);

            axTop = makeSubplot(viewIdx);
            topAxes(viewIdx) = axTop;
            yline(axTop, 50, '--', 'Color', 0.4 .* [1 1 1], ...
                'LineWidth', 1.5, 'HandleVisibility', 'off');
            conditionPlotOpts = conditionMarkerOptions(opts, conditionColors);
            plotJitteredGroupWithSummary(axTop, conditionValues, ...
                conditionColors, conditionLabels, conditionPlotOpts);
            ylim(axTop, [0 100]);
            styleDistributionAxes(axTop, conditionLabels);
            conditionTests = pairedPanelStats( ...
                conditionValues, conditionPairs, conditionStatsOpts);
            conditionTestsByView{viewIdx} = conditionTests;
            distributionData(aggregateIdx).(viewField).conditionTests = ...
                conditionTests;
            titleHandle = title(axTop, viewTitles{viewIdx}, ...
                'FontWeight', 'normal', 'FontSize', 17);
            set(titleHandle, 'Units', 'normalized', ...
                'Position', [0.5 1.24 0]);
            ylabel(axTop, 'Percent correct', 'FontSize', 17);
            if opts.showSourceAuditText
                addSourceAuditText(axTop, meanField, conditionStats.n);
            end

            axBottom = makeSubplot(viewIdx + 3);
            bottomAxes(viewIdx) = axBottom;
            yline(axBottom, 0, '--', 'Color', [0.35 0.35 0.35], ...
                'LineWidth', 1.2, 'HandleVisibility', 'off');
            deltaPlotOpts = deltaMarkerOptions(opts, deltaColors);
            plotJitteredGroupWithSummary(axBottom, deltaValues, ...
                deltaColors, deltaLabels, deltaPlotOpts);
            setEffectLimits(axBottom, deltaValues);
            styleDistributionAxes(axBottom, deltaLabels);
            deltaTests = oneSampleAndPairedDeltaStats( ...
                deltaValues(:, 1), deltaValues(:, 2), deltaStatsOpts);
            deltaTestsByView{viewIdx} = deltaTests;
            distributionData(aggregateIdx).(viewField).deltaTests = ...
                deltaTests;
            ylabel(axBottom, 'Percentage-point effect size', 'FontSize', 17);
            if opts.showSourceAuditText
                addSourceAuditText(axBottom, deltaField, deltaStats.n);
            end
        end

        conditionBounds = sharedAnnotationBounds( ...
            conditionValuesByView, conditionStatsByView, 50);
        performanceYLim = [0 100];
        deltaBounds = sharedAnnotationBounds( ...
            deltaValuesByView, deltaStatsByView, 0);
        deltaBaseYLim = sharedCurrentYLim(bottomAxes);
        deltaYLim = annotationRowYLim(deltaBaseYLim, deltaBounds, 3, opts);

        for viewIdx = 1:3
            viewField = lower(viewNames{viewIdx});

            ylim(topAxes(viewIdx), performanceYLim);
            tickInfo = applyDistributionYTicks(topAxes(viewIdx), 'performance');
            annotationOpts = sharedAnnotationOptions( ...
                opts, conditionBounds, performanceYLim, 'percentCorrect');
            conditionAnnotation = addSignificanceBrackets( ...
                topAxes(viewIdx), conditionPairs, ...
                conditionTestsByView{viewIdx}, annotationOpts);
            distributionData(aggregateIdx).(viewField).conditionAnnotation = ...
                conditionAnnotation;
            printPanelValidation(clusterTitle, viewTitles{viewIdx}, ...
                conditionLabels, conditionStatsByView{viewIdx}.n, ...
                conditionPairs, conditionTestsByView{viewIdx}, ...
                conditionAnnotation, tickInfo, 50, topAxes(viewIdx));
            conditionFamilyID = sprintf( ...
                '%s|C%d|%s|performance', ...
                datasetLabel(aggregateFits(aggregateIdx)), ...
                aggregateFits(aggregateIdx).clusterID, viewField);
            conditionStatsAudit = buildDistributionStatsAuditRows( ...
                aggregateFits(aggregateIdx), 'experimentDistribution', ...
                'performance', viewField, 'percentCorrect', ...
                conditionFamilyID, conditionTestsByView{viewIdx}, ...
                conditionAnnotation, rows, sourceBlocks, experimentIDs);
            statsAudit = [statsAudit; conditionStatsAudit]; %#ok<AGROW>
            printDistributionStatsAudit(conditionStatsAudit, opts.alpha);

            ylim(bottomAxes(viewIdx), deltaYLim);
            tickInfo = applyDistributionYTicks(bottomAxes(viewIdx), 'delta');
            annotationOpts = sharedAnnotationOptions( ...
                opts, deltaBounds, deltaYLim, 'standard');
            deltaAnnotation = addSignificanceBrackets( ...
                bottomAxes(viewIdx), deltaPairs, ...
                deltaTestsByView{viewIdx}, annotationOpts);
            distributionData(aggregateIdx).(viewField).deltaAnnotation = ...
                deltaAnnotation;
            printPanelValidation(clusterTitle, ...
                [viewTitles{viewIdx} ' deltas'], deltaLabels, ...
                deltaStatsByView{viewIdx}.n, deltaPairs, ...
                deltaTestsByView{viewIdx}, deltaAnnotation, ...
                tickInfo, 0, bottomAxes(viewIdx));
            deltaFamilyID = sprintf('%s|C%d|%s|delta', ...
                datasetLabel(aggregateFits(aggregateIdx)), ...
                aggregateFits(aggregateIdx).clusterID, viewField);
            deltaStatsAudit = buildDistributionStatsAuditRows( ...
                aggregateFits(aggregateIdx), 'deltaDistribution', ...
                'delta', viewField, 'deltaBias/deltaMask', ...
                deltaFamilyID, deltaTestsByView{viewIdx}, ...
                deltaAnnotation, rows, sourceBlocks, experimentIDs);
            statsAudit = [statsAudit; deltaStatsAudit]; %#ok<AGROW>
            printDistributionStatsAudit(deltaStatsAudit, opts.alpha);
        end
        addDistributionTitle(clusterTitle);
    end
end

function opts = applyDefaults(opts)
    defaults = struct( ...
        'figureVisible', 'on', ...
        'jitterWidth', 0.12, ...
        'pointSize', 34, ...
        'meanPointSize', 95, ...
        'pointAlpha', 0.68, ...
        'connectPairedValues', true, ...
        'showIQR', false, ...
        'showNS', true, ...
        'alpha', 0.05, ...
        'annotationLineSpacingFraction', 0.14, ...
        'annotationStarLabelGapFraction', 0.004, ...
        'annotationNSLabelGapFraction', 0.010, ...
        'annotationTopPadFraction', 0.08, ...
        'annotationFinalTopPadFraction', 0.05, ...
        'showSourceAuditText', false, ...
        'experimentIDsByBlock', [], ...
        'baselineModeByBlock', [], ...
        'statsOpts', struct('testDirection', 'two-sided'));
    names = fieldnames(defaults);
    for idx = 1:numel(names)
        if ~isfield(opts, names{idx}) || isempty(opts.(names{idx}))
            opts.(names{idx}) = defaults.(names{idx});
        end
    end
    opts.figureVisible = validatestring(opts.figureVisible, {'on', 'off'});
    if ~isstruct(opts.statsOpts)
        error('plotPowerClusterMeanDistributions:InvalidStatsOpts', ...
            'opts.statsOpts must be a struct.');
    end
    if ~isfield(opts.statsOpts, 'testDirection') || ...
            isempty(opts.statsOpts.testDirection)
        opts.statsOpts.testDirection = 'two-sided';
    end
    opts.statsOpts.testDirection = validatestring( ...
        opts.statsOpts.testDirection, {'two-sided', 'one-sided'});
end

function statsOpts = performanceStatsOptions(opts)
    statsOpts = opts.statsOpts;
    statsOpts.alpha = opts.alpha;
    statsOpts.plannedTails = {'left', 'right', 'right'};
    statsOpts.comparisonLabels = { ...
        'Baseline vs Con-Opto', ...
        'Con-Opto vs Incon-Opto', ...
        'Baseline vs Incon-Opto'};
    statsOpts.plannedAlternativeHypotheses = { ...
        'Con-Opto > Baseline', ...
        'Con-Opto > Incon-Opto', ...
        'Incon-Opto < Baseline'};
    statsOpts.twoSidedAlternativeHypotheses = { ...
        'Baseline ~= Con-Opto', ...
        'Con-Opto ~= Incon-Opto', ...
        'Baseline ~= Incon-Opto'};
end

function label = datasetLabel(item)
    if isfield(item, 'columnTargetLabel') && ...
            ~isempty(item.columnTargetLabel)
        label = char(string(item.columnTargetLabel));
    else
        label = 'unknown';
    end
end

function rows = validateRows(rows, mdl)
    rows = unique(rows(:)', 'stable');
    candidateFields = {'meanPsychometricHorizontal', ...
        'meanPsychometricVertical', 'meanPsychometricMerged'};
    nRows = 0;
    for idx = 1:numel(candidateFields)
        if isfield(mdl, candidateFields{idx})
            nRows = max(nRows, size(mdl.(candidateFields{idx}), 1));
        end
    end
    invalid = ~isfinite(rows) | rows < 1 | rows > nRows | rows ~= round(rows);
    if any(invalid)
        warning('plotPowerClusterMeanDistributions:InvalidRows', ...
            'Ignoring %d invalid mdl row indices.', sum(invalid));
        rows(invalid) = [];
    end
end

function values = readSessionValues(mdl, fieldName, rows, nColumns, clusterTitle)
    values = nan(numel(rows), nColumns);
    if ~isfield(mdl, fieldName)
        warning('plotPowerClusterMeanDistributions:MissingField', ...
            '%s is missing for %s; the corresponding panel will be empty.', ...
            fieldName, clusterTitle);
        return;
    end

    source = mdl.(fieldName);
    if size(source, 2) < nColumns
        warning('plotPowerClusterMeanDistributions:ShortField', ...
            '%s has %d columns, but %d are required for %s.', ...
            fieldName, size(source, 2), nColumns, clusterTitle);
        return;
    end
    values = source(rows, 1:nColumns);
end

function stats = summarizeColumns(values)
    stats.mean = mean(values, 1, 'omitnan');
    stats.n = sum(isfinite(values), 1);
    stats.sem = nan(1, size(values, 2));
    for columnIdx = 1:size(values, 2)
        finiteValues = values(isfinite(values(:, columnIdx)), columnIdx);
        if numel(finiteValues) >= 2
            stats.sem(columnIdx) = std(finiteValues, 0) ./ sqrt(numel(finiteValues));
        elseif numel(finiteValues) == 1
            stats.sem(columnIdx) = 0;
        end
    end
end

function verifyPointCounts(counts, expectedCount, clusterTitle, viewTitle, labels)
    for valueIdx = 1:numel(counts)
        if counts(valueIdx) ~= expectedCount
            warning('plotPowerClusterMeanDistributions:MissingSessionValues', ...
                ['%s, %s, %s: plotting %d of %d experiments because ' ...
                'the remaining session values are missing or nonfinite.'], ...
                clusterTitle, viewTitle, labels{valueIdx}, ...
                counts(valueIdx), expectedCount);
        end
    end
end

function printMeans(viewTitle, conditionStats, deltaStats)
    fprintf(['  %s session means: baseline=%0.3f (n=%d), ' ...
        'con=%0.3f (n=%d), incon=%0.3f (n=%d), ' ...
        'deltaBias=%0.3f (n=%d), deltaMask=%0.3f (n=%d)\n'], ...
        viewTitle, ...
        conditionStats.mean(1), conditionStats.n(1), ...
        conditionStats.mean(2), conditionStats.n(2), ...
        conditionStats.mean(3), conditionStats.n(3), ...
        deltaStats.mean(1), deltaStats.n(1), ...
        deltaStats.mean(2), deltaStats.n(2));
end

function styleDistributionAxes(ax, labels)
    xlim(ax, [0.5, numel(labels) + 0.5]);
    xticks(ax, 1:numel(labels));
    xticklabels(ax, labels);
    xtickangle(ax, 18);
    box(ax, 'off');
    set(ax, 'LineWidth', 2, 'TickDir', 'out', ...
        'TickLength', [0.01 0.01], 'FontName', 'FreeSans', ...
        'FontSize', 17);
end

function plotOpts = conditionMarkerOptions(opts, colors)
    plotOpts = opts;
    plotOpts.groupMarkers = {'o', '^', 'v'};
    plotOpts.groupFaceColors = [1 1 1; colors(2, :); colors(3, :)];
    plotOpts.groupEdgeColors = zeros(3, 3);
    plotOpts.semColor = [0 0 0];
end

function plotOpts = deltaMarkerOptions(opts, colors)
    plotOpts = opts;
    plotOpts.groupMarkers = {'o', 's'};
    plotOpts.groupFaceColors = colors;
    plotOpts.groupEdgeColors = zeros(2, 3);
    plotOpts.semColor = [0 0 0];
end

function bounds = sharedAnnotationBounds(valuesByView, statsByView, referenceValue)
    allValues = referenceValue;
    allUpper = referenceValue;
    allLower = referenceValue;
    for viewIdx = 1:numel(valuesByView)
        values = valuesByView{viewIdx};
        stats = statsByView{viewIdx};
        allValues = [allValues; values(isfinite(values))]; %#ok<AGROW>
        upperSummary = stats.mean + stats.sem;
        lowerSummary = stats.mean - stats.sem;
        allUpper = [allUpper; upperSummary(isfinite(upperSummary))']; %#ok<AGROW>
        allLower = [allLower; lowerSummary(isfinite(lowerSummary))']; %#ok<AGROW>
    end
    bounds = [min([allValues; allLower]), max([allValues; allUpper])];
end

function limits = sharedCurrentYLim(axesHandles)
    limitsByPanel = nan(numel(axesHandles), 2);
    for idx = 1:numel(axesHandles)
        limitsByPanel(idx, :) = ylim(axesHandles(idx));
    end
    limits = [min(limitsByPanel(:, 1)), max(limitsByPanel(:, 2))];
end

function limits = annotationRowYLim(baseLimits, dataBounds, nLanes, opts)
    dataRange = dataBounds(2) - dataBounds(1);
    if ~isfinite(dataRange) || dataRange <= 0
        dataRange = max(abs(dataBounds(2)), 1);
    end
    highestTextY = dataBounds(2) + ...
        opts.annotationTopPadFraction .* dataRange + ...
        (nLanes - 1) .* opts.annotationLineSpacingFraction .* dataRange + ...
        opts.annotationNSLabelGapFraction .* dataRange;
    finalTop = highestTextY + ...
        opts.annotationFinalTopPadFraction .* dataRange;
    limits = [baseLimits(1), max(baseLimits(2), finalTop)];
end

function annotationOpts = sharedAnnotationOptions( ...
        opts, bounds, finalLimits, annotationMode)
    annotationOpts = opts;
    annotationOpts.annotationMode = annotationMode;
    annotationOpts.annotationDataBottom = bounds(1);
    annotationOpts.annotationDataTop = bounds(2);
    annotationOpts.annotationFinalYLim = finalLimits;
    annotationOpts.lineSpacingFraction = ...
        opts.annotationLineSpacingFraction;
    annotationOpts.starLabelGapFraction = ...
        opts.annotationStarLabelGapFraction;
    annotationOpts.nsLabelGapFraction = ...
        opts.annotationNSLabelGapFraction;
    annotationOpts.topPadFraction = opts.annotationTopPadFraction;
    annotationOpts.finalTopPadFraction = ...
        opts.annotationFinalTopPadFraction;
end

function [blockIndices, experimentIDs, baselineModes] = auditSessionMetadata(item, rows, opts)
    blockIndices = nan(numel(rows), 1);
    if isfield(item, 'mdlRowIndices') && ...
            isfield(item, 'sourceBlockIndices')
        itemRows = item.mdlRowIndices(:);
        itemBlocks = item.sourceBlockIndices(:);
        if numel(itemRows) == numel(itemBlocks)
            for rowIdx = 1:numel(rows)
                matchIdx = find(itemRows == rows(rowIdx), 1);
                if ~isempty(matchIdx)
                    blockIndices(rowIdx) = itemBlocks(matchIdx);
                end
            end
        end
    end

    experimentIDs = strings(numel(rows), 1);
    baselineModes = lookupBaselineModes(blockIndices, opts);
    if isempty(opts.experimentIDsByBlock)
        return;
    end
    idsByBlock = opts.experimentIDsByBlock;
    for rowIdx = 1:numel(rows)
        blockIdx = blockIndices(rowIdx);
        if isfinite(blockIdx) && blockIdx == round(blockIdx) && ...
                blockIdx >= 1 && blockIdx <= numel(idsByBlock)
            experimentIDs(rowIdx) = string(idsByBlock(blockIdx));
        end
    end
end

function baselineModes = lookupBaselineModes(blockIndices, opts)
    baselineModes = strings(numel(blockIndices), 1);
    if ~isfield(opts, 'baselineModeByBlock') || isempty(opts.baselineModeByBlock)
        return;
    end
    modesByBlock = string(opts.baselineModeByBlock(:));
    for rowIdx = 1:numel(blockIndices)
        blockIdx = blockIndices(rowIdx);
        if isfinite(blockIdx) && blockIdx == round(blockIdx) && ...
                blockIdx >= 1 && blockIdx <= numel(modesByBlock)
            baselineModes(rowIdx) = modesByBlock(blockIdx);
        end
    end
end

function assertAuditCounts(audit, expectedCounts, clusterTitle, ...
        panelTitle, metricLabels)
    for metricIdx = 1:numel(metricLabels)
        metric = string(metricLabels{metricIdx});
        metricRows = audit.conditionOrMetric == metric;
        validRows = metricRows & audit.isValidForPlot;
        nAudit = sum(validRows);
        if nAudit ~= expectedCounts(metricIdx)
            error('plotPowerClusterMeanDistributions:AuditCountMismatch', ...
                ['%s, %s, %s: audit has %d valid points, but plotting ' ...
                'summary expects %d.'], ...
                clusterTitle, panelTitle, metric, ...
                nAudit, expectedCounts(metricIdx));
        end
        excluded = audit(metricRows & ~audit.isValidForPlot, :);
        for excludedIdx = 1:height(excluded)
            fprintf(['  Excluded source value: %s | %s | sessionRow=%g | ' ...
                'block=%g | reason=%s\n'], ...
                panelTitle, metric, excluded.sessionRowIndex(excludedIdx), ...
                excluded.blockIndex(excludedIdx), ...
                excluded.reasonExcluded(excludedIdx));
        end
    end
end

function printSourceAuditSummary(audit, clusterTitle)
    metrics = unique(audit.conditionOrMetric, 'stable');
    for metricIdx = 1:numel(metrics)
        metricRows = audit.conditionOrMetric == metrics(metricIdx);
        validRows = metricRows & audit.isValidForPlot;
        values = audit.plottedValue(validRows);
        [valueMean, valueSEM, valueMedian, valueMin, valueMax] = ...
            summarizeAuditValues(values);
        exampleRow = find(metricRows, 1);
        figureLabel = char(audit.figureType(exampleRow));
        datasetLabel = char(audit.datasetLabel(exampleRow));
        panelColumnLabel = char(audit.panelColumn(exampleRow));
        panelRowLabel = char(audit.panelRow(exampleRow));
        metricLabel = char(metrics(metricIdx));
        sourceFieldLabel = char(audit.sourceField(exampleRow));
        baselineModeLabel = char(strjoin(unique(audit.baselineMode(metricRows)), ','));
        fprintf(['  Source audit: figure=%s | dataset=%s | cluster=%g | ' ...
            'panel=%s/%s | metric=%s | source=%s | n=%d | BL=%s | ' ...
            'mean=%0.4g | SEM=%0.4g | median=%0.4g | ' ...
            'min=%0.4g | max=%0.4g\n'], ...
            figureLabel, datasetLabel, ...
            audit.powerClusterID(exampleRow), panelColumnLabel, ...
            panelRowLabel, metricLabel, ...
            sourceFieldLabel, sum(validRows), ...
            baselineModeLabel, ...
            valueMean, valueSEM, valueMedian, valueMin, valueMax);
        if ~contains(clusterTitle, sprintf('Power cluster %g', ...
                audit.powerClusterID(exampleRow)))
            warning('plotPowerClusterMeanDistributions:AuditClusterLabel', ...
                'Audit cluster label does not match figure title.');
        end
    end
end

function [valueMean, valueSEM, valueMedian, valueMin, valueMax] = ...
        summarizeAuditValues(values)
    if isempty(values)
        [valueMean, valueSEM, valueMedian, valueMin, valueMax] = ...
            deal(nan);
        return;
    end
    valueMean = mean(values);
    valueMedian = median(values);
    valueMin = min(values);
    valueMax = max(values);
    if numel(values) >= 2
        valueSEM = std(values, 0) ./ sqrt(numel(values));
    else
        valueSEM = 0;
    end
end

function addSourceAuditText(ax, sourceField, counts)
    text(ax, 0.02, 0.98, sprintf('%s | n=%s', sourceField, ...
        mat2str(counts)), ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', ...
        'FontSize', 7, ...
        'Color', [0.25 0.25 0.25], ...
        'Interpreter', 'none', ...
        'Clipping', 'on', ...
        'Tag', 'sourceAuditText');
end

function printPanelValidation(clusterTitle, panelTitle, labels, counts, ...
        pairs, results, annotationInfo, tickInfo, referenceValue, ax)
    fprintf(['  Validation: %s | %s | n/group=%s | ylim=%s | ' ...
        'reference=%g | annotations=%d\n'], ...
        clusterTitle, panelTitle, mat2str(counts), ...
        mat2str(tickInfo.limits, 4), referenceValue, annotationInfo.nDrawn);
    for resultIdx = 1:numel(results)
        if pairs(resultIdx, 1) == pairs(resultIdx, 2)
            comparisonLabel = sprintf('%s vs 0', ...
                labels{pairs(resultIdx, 1)});
        else
            comparisonLabel = sprintf('%s vs %s', ...
                labels{pairs(resultIdx, 1)}, labels{pairs(resultIdx, 2)});
        end
        if annotationInfo.drawnMask(resultIdx)
            annotationStatus = 'drawn';
        else
            annotationStatus = ['skipped: ' ...
                annotationInfo.skippedReasons{resultIdx}];
        end
        fprintf(['    %s: %s, n=%d, raw p=%0.4g, Holm p=%0.4g, ' ...
            '%s, annotation %s, lineY=%0.4g, textY=%0.4g, ' ...
            'lineSpacing=%0.4g, labelGap=%0.4g, fontSize=%g\n'], ...
            comparisonLabel, results(resultIdx).test, ...
            results(resultIdx).n, results(resultIdx).rawP, ...
            results(resultIdx).adjustedP, results(resultIdx).star, ...
            annotationStatus, annotationInfo.pairs(resultIdx).yLine, ...
            annotationInfo.pairs(resultIdx).yText, ...
            annotationInfo.lineSpacing, ...
            annotationInfo.pairs(resultIdx).labelGap, ...
            annotationInfo.pairs(resultIdx).fontSize);
    end
    fprintf(['      annotation geometry: mode=%s, yDataMax=%0.4g, ' ...
        'yRangeBase=%0.4g, starLabelGap=%0.4g, nsLabelGap=%0.4g, ' ...
        'final ylim=%s\n'], ...
        annotationInfo.mode, annotationInfo.yDataMax, ...
        annotationInfo.yRangeBase, annotationInfo.starLabelGap, ...
        annotationInfo.nsLabelGap, mat2str(annotationInfo.finalYLim, 4));
    labelsFound = findall(ax, 'Tag', 'sigLabel');
    for labelIdx = 1:numel(labelsFound)
        fprintf('      sigLabel verify: String=%s, FontSize=%g\n', ...
            labelsFound(labelIdx).String, labelsFound(labelIdx).FontSize);
    end
end

function setEffectLimits(ax, values)
    finiteValues = values(isfinite(values));
    if isempty(finiteValues)
        ylim(ax, [-10 10]);
        return;
    end
    valueRange = max(finiteValues) - min(finiteValues);
    padding = max(5, 0.12 .* max(valueRange, 1));
    limits = [min([finiteValues; 0]) - padding, ...
        max([finiteValues; 0]) + padding];
    if limits(1) == limits(2)
        limits = limits + [-5 5];
    end
    ylim(ax, limits);
end

function titleText = distributionTitle(item, nExperiments)
    titleText = sprintf('Power cluster %d', item.clusterID);
    if isfield(item, 'columnTargetLabel')
        titleText = sprintf('%s, Columns %s', ...
            titleText, item.columnTargetLabel);
    end
    titleText = sprintf('%s (n_{expt}=%d)', titleText, nExperiments);
end

function addDistributionTitle(titleText)
    annotation(gcf, 'textbox', [0.12 0.920 0.76 0.055], ...
        'String', titleText, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', 16, ...
        'Interpreter', 'tex', ...
        'EdgeColor', 'none');
end
