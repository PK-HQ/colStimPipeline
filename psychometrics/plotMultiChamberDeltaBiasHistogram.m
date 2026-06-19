function results = plotMultiChamberDeltaBiasHistogram(opts)
% Plot verified Stage 1 merged experiment-wise deltaBias histograms.
%
% This lightweight analysis uses saved distribution source audits only. It
% does not load raw trials, rerun fits, rerun clustering, run Stage 2
% permutation analysis, or modify existing psychometric outputs.

    if nargin < 1 || isempty(opts)
        opts = struct();
    end

    opts = fillDefaultOptions(opts);
    if ~exist(opts.outputDir, 'dir')
        mkdir(opts.outputDir);
    end

    chamberConfigs = defaultChamberConfigs(opts);
    nChambers = numel(chamberConfigs);
    chamberResults = repmat(emptyChamberResult(), nChambers, 1);
    clusterStats = table();

    fprintf('Verified Stage 1 multi-chamber merged deltaBias histogram\n');
    fprintf('Source: saved distributionSourceAudit tables only. Stage 2 is not run.\n');

    for chamberIdx = 1:nChambers
        cfg = chamberConfigs(chamberIdx);
        fprintf('\nLoading %s-%s audit: %s\n', cfg.monkeyID, cfg.chamber, cfg.auditPath);
        audit = loadAuditTable(cfg.auditPath);
        printAuditColumns(audit);

        [rows090, rowsControl] = extractMergedDeltaBiasRows(audit, cfg);
        chamberResult = summarizeChamberRows(rows090, rowsControl, cfg);
        chamberResults(chamberIdx) = chamberResult;
        clusterStats = [clusterStats; chamberResult.clusterStats]; %#ok<AGROW>
        printChamberSummary(chamberResult);
    end

    displayedValues = collectDisplayedValues(chamberResults);
    binEdges = commonBinEdges(displayedValues, opts.binWidth);
    histogramBins = buildHistogramBinTable(chamberResults, binEdges);
    panelStats = buildPanelStatsTable(chamberResults);

    plotFigure = plotHistogramFigure(chamberResults, binEdges, panelStats, opts);
    saveOutputs(plotFigure, panelStats, histogramBins, clusterStats, opts);
    validateResults(chamberResults, binEdges, histogramBins);

    results = struct();
    results.chamberResults = chamberResults;
    results.clusterStats = clusterStats;
    results.panelStats = panelStats;
    results.histogramBins = histogramBins;
    results.binEdges = binEdges;
    results.outputDir = opts.outputDir;
end

function opts = fillDefaultOptions(opts)
    if ~isfield(opts, 'modelName') || isempty(opts.modelName)
        opts.modelName = 'weibullfreeAll';
    end
    if ~isfield(opts, 'summaryRoot') || isempty(opts.summaryRoot)
        opts.summaryRoot = 'Y:\';
    end
    if ~isfield(opts, 'outputDir') || isempty(opts.outputDir)
        opts.outputDir = fullfile('Y:\users\PK\colStimPipeline', ...
            'outputs', 'multiChamber');
    end
    if ~isfield(opts, 'binWidth') || isempty(opts.binWidth)
        opts.binWidth = 5;
    end
    if ~isfield(opts, 'showRawMarkers') || isempty(opts.showRawMarkers)
        opts.showRawMarkers = false;
    end
    if ~isfield(opts, 'purple090') || isempty(opts.purple090)
        opts.purple090 = [0.55 0.25 0.80];
    end
    if ~isfield(opts, 'orangeControl') || isempty(opts.orangeControl)
        opts.orangeControl = [0.95 0.62 0.12];
    end
end

function configs = defaultChamberConfigs(opts)
    configs = struct( ...
        'animal', {'Chip', 'Chip', 'Pepper'}, ...
        'monkeyID', {'M1', 'M1', 'M2'}, ...
        'chamber', {'L', 'R', 'R'}, ...
        'auditPath', {'', '', ''});

    for idx = 1:numel(configs)
        configs(idx).auditPath = fullfile(opts.summaryRoot, configs(idx).animal, ...
            'Meta', 'summary', sprintf('distributionSourceAudit_%s_%s_%s.mat', ...
            configs(idx).animal, configs(idx).chamber, opts.modelName));
    end
end

function result = emptyChamberResult()
    result = struct( ...
        'animal', '', ...
        'monkeyID', '', ...
        'chamber', '', ...
        'auditPath', '', ...
        'sourceColumns', {{}}, ...
        'clusterStats', table(), ...
        'experimentTable090', table(), ...
        'experimentTableControl', table(), ...
        'retainedValues090', [], ...
        'retainedValuesControl', [], ...
        'includedClusterIDs', [], ...
        'totalClusters', 0, ...
        'controlUnmatchedIDs', strings(0, 1));
end

function audit = loadAuditTable(auditPath)
    if ~isfile(auditPath)
        error('plotMultiChamberDeltaBiasHistogram:MissingAudit', ...
            'Required source audit MAT file is missing: %s', auditPath);
    end

    loaded = load(auditPath);
    varNames = fieldnames(loaded);
    tableVars = {};
    for idx = 1:numel(varNames)
        if istable(loaded.(varNames{idx}))
            tableVars{end + 1} = varNames{idx}; %#ok<AGROW>
        end
    end

    if isempty(tableVars)
        error('plotMultiChamberDeltaBiasHistogram:MissingAuditTable', ...
            'No table variable found inside %s. Variables: %s', ...
            auditPath, strjoin(varNames, ', '));
    end

    preferredNames = {'distributionSourceAudit', 'audit'};
    auditVar = tableVars{1};
    for idx = 1:numel(preferredNames)
        if any(strcmp(tableVars, preferredNames{idx}))
            auditVar = preferredNames{idx};
            break
        end
    end

    audit = loaded.(auditVar);
    fprintf('  using table variable "%s" with %d rows.\n', auditVar, height(audit));
end

function printAuditColumns(audit)
    fprintf('  columns: %s\n', strjoin(audit.Properties.VariableNames, ', '));
end

function [rows090, rowsControl] = extractMergedDeltaBiasRows(audit, cfg)
    requiredVars = {'datasetLabel', 'figureType', 'panelColumn', ...
        'conditionOrMetric', 'powerClusterID', 'experimentID', ...
        'sessionRowIndex', 'blockIndex', 'sourceField', ...
        'sourceSubscriptOrColumn', 'plottedValue', 'isValidForPlot'};
    assertRequiredColumns(audit, requiredVars, cfg);

    figureType = lower(string(audit.figureType));
    panelColumn = lower(string(audit.panelColumn));
    conditionOrMetric = lower(string(audit.conditionOrMetric));
    validMask = logical(audit.isValidForPlot);

    keep = figureType == "deltadistribution" & ...
        panelColumn == "merged" & ...
        conditionOrMetric == "deltabias" & ...
        validMask & isfinite(audit.plottedValue) & ...
        isfinite(audit.powerClusterID);

    rows = audit(keep, :);
    if isempty(rows)
        error('plotMultiChamberDeltaBiasHistogram:NoMergedDeltaBiasRows', ...
            ['No valid merged deltaBias rows found for %s %s. Expected ' ...
            'figureType=deltaDistribution, panelColumn=merged, ' ...
            'conditionOrMetric=deltaBias, isValidForPlot=true.'], ...
            cfg.animal, cfg.chamber);
    end

    labels = normalizeLabel(string(rows.datasetLabel));
    isControl = contains(labels, "45") & contains(labels, "135");
    is090 = contains(labels, "0") & contains(labels, "90") & ~isControl;
    if ~any(is090)
        is090 = ~isControl;
    end

    rows090 = rows(is090, :);
    rowsControl = rows(isControl, :);

    assertUniqueExperiments(rows090, cfg, '0/90');
    if ~isempty(rowsControl)
        assertUniqueExperiments(rowsControl, cfg, '45/135');
    end

    fprintf(['  Stage 1 empirical delta source: %s, %s. These rows are ' ...
        'saved experiment-wise row-2-column-3 merged deltaBias dots, not fits.\n'], ...
        char(rows090.sourceField(1)), char(rows090.sourceSubscriptOrColumn(1)));
    fprintf('  0/90 candidate rows: %d; 45/135 control candidate rows: %d\n', ...
        height(rows090), height(rowsControl));
end

function labels = normalizeLabel(labels)
    labels = lower(labels);
    labels = erase(labels, char(176));
    labels = replace(labels, " ", "");
end

function assertUniqueExperiments(rows, cfg, conditionPair)
    experimentStrings = string(rows.experimentID);
    uniqueExperiments = unique(experimentStrings);
    if numel(uniqueExperiments) ~= height(rows)
        duplicateCounts = zeros(numel(uniqueExperiments), 1);
        for experimentIdx = 1:numel(uniqueExperiments)
            duplicateCounts(experimentIdx) = sum(experimentStrings == uniqueExperiments(experimentIdx));
        end
        duplicateIDs = uniqueExperiments(duplicateCounts > 1);
        error('plotMultiChamberDeltaBiasHistogram:DuplicateExperiments', ...
            'Duplicate %s merged deltaBias rows for %s %s: %s', ...
            conditionPair, cfg.animal, cfg.chamber, strjoin(duplicateIDs, ', '));
    end
end

function assertRequiredColumns(audit, requiredVars, cfg)
    missing = setdiff(requiredVars, audit.Properties.VariableNames);
    if ~isempty(missing)
        error('plotMultiChamberDeltaBiasHistogram:MissingColumns', ...
            'Audit for %s %s is missing required columns: %s', ...
            cfg.animal, cfg.chamber, strjoin(missing, ', '));
    end
end

function result = summarizeChamberRows(rows090, rowsControl, cfg)
    clusterIDs = unique(rows090.powerClusterID(:));
    clusterIDs = clusterIDs(isfinite(clusterIDs));
    clusterStats = table();
    experimentTable090 = table();

    for clusterIdx = 1:numel(clusterIDs)
        clusterID = clusterIDs(clusterIdx);
        clusterRows = rows090(rows090.powerClusterID == clusterID, :);
        values = clusterRows.plottedValue(:);
        finiteValues = values(isfinite(values));
        stats = valueSummaryStats(finiteValues);
        [pValue, reason, testValid] = runSignrank(finiteValues, 'right');
        included = stats.n > 0 && stats.meanDeltaBias > 0 && ...
            isfinite(pValue) && pValue < 0.05;
        if included
            inclusionReason = "mean>0 and one-sided signrank p<0.05";
        else
            inclusionReason = string(reason);
            if inclusionReason == ""
                inclusionReason = "failed mean>0 and one-sided p<0.05";
            end
        end

        clusterStats = [clusterStats; table( ...
            string(cfg.monkeyID), string(cfg.animal), string(cfg.chamber), ...
            clusterID, stats.n, stats.meanDeltaBias, stats.medianDeltaBias, ...
            stats.semDeltaBias, stats.nPositive, stats.nNegative, ...
            stats.nZero, pValue, testValid, included, inclusionReason, ...
            'VariableNames', {'monkeyID', 'animal', 'chamber', ...
            'clusterID', 'n', 'meanDeltaBias', 'medianDeltaBias', ...
            'semDeltaBias', 'nPositive', 'nNegative', 'nZero', ...
            'oneSidedP', 'testValid', 'included', 'inclusionReason'})]; %#ok<AGROW>

        nRows = height(clusterRows);
        experimentTable090 = [experimentTable090; makeExperimentRows( ...
            cfg, clusterRows, '0/90', included, inclusionReason)]; %#ok<AGROW>
    end

    includedClusterIDs = clusterStats.clusterID(clusterStats.included);
    retainedMask090 = experimentTable090.includedInHistogram;
    retainedValues090 = experimentTable090.deltaBiasMerged(retainedMask090);

    experimentTableControl = table();
    retainedValuesControl = [];
    controlUnmatchedIDs = strings(0, 1);
    if strcmp(cfg.monkeyID, 'M2') && strcmp(cfg.chamber, 'R') && ~isempty(rowsControl)
        controlInClusters = rowsControl(ismember(rowsControl.powerClusterID, includedClusterIDs), :);
        retainedIDs090 = string(experimentTable090.experimentID(retainedMask090));
        matchedControl = ismember(string(controlInClusters.experimentID), retainedIDs090);
        unmatchedControl = controlInClusters(~matchedControl, :);
        controlUnmatchedIDs = string(unmatchedControl.experimentID);
        if ~isempty(controlUnmatchedIDs)
            fprintf('  M2-R 45/135 unmatched controls excluded: %s\n', ...
                strjoin(controlUnmatchedIDs, ', '));
        end
        matchedRows = controlInClusters(matchedControl, :);
        experimentTableControl = makeExperimentRows(cfg, matchedRows, '45/135', ...
            true, "matched experiment ID and included 0/90 cluster");
        retainedValuesControl = experimentTableControl.deltaBiasMerged;
    end

    if strcmp(cfg.monkeyID, 'M2') && strcmp(cfg.chamber, 'R')
        fprintf('\nM2-R cluster inclusion verification from saved 0/90 audit rows:\n');
        for rowIdx = 1:height(clusterStats)
            row = clusterStats(rowIdx, :);
            if row.included
                status = 'included';
            else
                status = 'excluded';
            end
            fprintf('  C%d: n=%d, mean=%0.3f, median=%0.3f, p_right=%0.4g, %s (%s)\n', ...
                row.clusterID, row.n, row.meanDeltaBias, row.medianDeltaBias, ...
                row.oneSidedP, status, row.inclusionReason);
        end
        c3Rows = clusterStats(clusterStats.clusterID == 3, :);
        if ~isempty(c3Rows) && ~c3Rows.included
            error('plotMultiChamberDeltaBiasHistogram:M2RC3Excluded', ...
                ['M2-R C3 fails the saved-audit inclusion rule: n=%d, ' ...
                'mean=%0.3f, median=%0.3f, p_right=%0.4g. Reason: %s'], ...
                c3Rows.n, c3Rows.meanDeltaBias, c3Rows.medianDeltaBias, ...
                c3Rows.oneSidedP, c3Rows.inclusionReason);
        end
    end

    result = emptyChamberResult();
    result.animal = cfg.animal;
    result.monkeyID = cfg.monkeyID;
    result.chamber = cfg.chamber;
    result.auditPath = cfg.auditPath;
    result.sourceColumns = rows090.Properties.VariableNames;
    result.clusterStats = clusterStats;
    result.experimentTable090 = experimentTable090;
    result.experimentTableControl = experimentTableControl;
    result.retainedValues090 = retainedValues090;
    result.retainedValuesControl = retainedValuesControl;
    result.includedClusterIDs = includedClusterIDs;
    result.totalClusters = numel(clusterIDs);
    result.controlUnmatchedIDs = controlUnmatchedIDs;
end

function rowsOut = makeExperimentRows(cfg, sourceRows, conditionPair, included, reason)
    nRows = height(sourceRows);
    rowsOut = table( ...
        repmat(string(cfg.monkeyID), nRows, 1), ...
        repmat(string(cfg.animal), nRows, 1), ...
        repmat(string(cfg.chamber), nRows, 1), ...
        repmat(string(conditionPair), nRows, 1), ...
        string(sourceRows.experimentID), ...
        sourceRows.sessionRowIndex, ...
        sourceRows.blockIndex, ...
        sourceRows.powerClusterID, ...
        sourceRows.plottedValue, ...
        repmat(included, nRows, 1), ...
        repmat(included, nRows, 1), ...
        repmat(string(reason), nRows, 1), ...
        'VariableNames', {'monkeyID', 'animal', 'chamber', ...
        'conditionPair', 'experimentID', 'sessionRowIndex', 'blockIndex', ...
        'clusterID', 'deltaBiasMerged', 'clusterIncluded', ...
        'includedInHistogram', 'exclusionReason'});
end

function stats = valueSummaryStats(values)
    values = values(isfinite(values));
    stats = struct();
    stats.n = numel(values);
    if isempty(values)
        stats.meanDeltaBias = NaN;
        stats.medianDeltaBias = NaN;
        stats.semDeltaBias = NaN;
        stats.nPositive = 0;
        stats.nNegative = 0;
        stats.nZero = 0;
        return;
    end
    stats.meanDeltaBias = mean(values);
    stats.medianDeltaBias = median(values);
    if numel(values) > 1
        stats.semDeltaBias = std(values) ./ sqrt(numel(values));
    else
        stats.semDeltaBias = NaN;
    end
    stats.nPositive = sum(values > 0);
    stats.nNegative = sum(values < 0);
    stats.nZero = sum(values == 0);
end

function [pValue, reason, testValid] = runSignrank(values, tail)
    values = values(isfinite(values));
    pValue = NaN;
    reason = "";
    testValid = false;
    if numel(values) < 2
        reason = "insufficient finite data for signrank";
        return;
    end
    if all(values == 0)
        reason = "all values equal zero";
        return;
    end
    try
        pValue = signrank(values, 0, 'tail', tail);
        testValid = isfinite(pValue);
    catch err
        reason = "signrank failed: " + string(err.message);
    end
end

function values = collectDisplayedValues(chamberResults)
    values = [];
    for chamberIdx = 1:numel(chamberResults)
        values = [values; chamberResults(chamberIdx).retainedValues090(:)]; %#ok<AGROW>
        values = [values; chamberResults(chamberIdx).retainedValuesControl(:)]; %#ok<AGROW>
    end
end

function binEdges = commonBinEdges(values, binWidth)
    values = values(isfinite(values));
    if isempty(values)
        error('plotMultiChamberDeltaBiasHistogram:NoRetainedValues', ...
            'No experiments passed the positive-effect cluster filter.');
    end
    lowerEdge = floor(min([values(:); 0]) ./ binWidth) .* binWidth;
    upperEdge = ceil(max([values(:); 0]) ./ binWidth) .* binWidth;
    if lowerEdge == upperEdge
        lowerEdge = lowerEdge - binWidth;
        upperEdge = upperEdge + binWidth;
    end
    binEdges = lowerEdge:binWidth:upperEdge;
    assert(numel(binEdges) >= 2, 'At least two histogram bin edges are required.');
    assert(any(abs(binEdges) < eps), 'Common histogram edges must include zero.');
end

function histogramBins = buildHistogramBinTable(chamberResults, binEdges)
    histogramBins = table();
    for chamberIdx = 1:numel(chamberResults)
        result = chamberResults(chamberIdx);
        histogramBins = [histogramBins; makeBinRows(result, '0/90', ...
            result.retainedValues090, binEdges)]; %#ok<AGROW>
        if ~isempty(result.retainedValuesControl)
            histogramBins = [histogramBins; makeBinRows(result, '45/135', ...
                result.retainedValuesControl, binEdges)]; %#ok<AGROW>
        end
    end
end

function rows = makeBinRows(result, conditionPair, values, binEdges)
    [counts, edges] = histcounts(values, binEdges);
    assert(sum(counts) == numel(values), 'Histogram counts lost values.');
    if numel(values) > 0
        proportions = counts(:) ./ numel(values);
        assert(abs(sum(proportions) - 1) < 1e-12, ...
            'Histogram proportions do not sum to 1.');
    else
        proportions = NaN(numel(counts), 1);
    end
    centers = edges(1:end-1)' + diff(edges(:)) ./ 2;
    nBins = numel(counts);
    rows = table( ...
        repmat(string(result.monkeyID), nBins, 1), ...
        repmat(string(result.animal), nBins, 1), ...
        repmat(string(result.chamber), nBins, 1), ...
        repmat(string(conditionPair), nBins, 1), ...
        (1:nBins)', edges(1:end-1)', edges(2:end)', centers, ...
        counts(:), proportions, ...
        'VariableNames', {'monkeyID', 'animal', 'chamber', ...
        'conditionPair', 'binIndex', 'binLeft', 'binRight', ...
        'binCenter', 'count', 'proportion'});
end

function panelStats = buildPanelStatsTable(chamberResults)
    panelStats = table();
    for chamberIdx = 1:numel(chamberResults)
        result = chamberResults(chamberIdx);
        panelStats = [panelStats; makePanelStatsRow(result, '0/90', ...
            result.retainedValues090, 'right')]; %#ok<AGROW>
        if strcmp(result.monkeyID, 'M2') && strcmp(result.chamber, 'R')
            panelStats = [panelStats; makePanelStatsRow(result, '45/135', ...
                result.retainedValuesControl, 'both')]; %#ok<AGROW>
        end
    end
end

function row = makePanelStatsRow(result, conditionPair, values, tail)
    stats = valueSummaryStats(values);
    [pValue, reason, testValid] = runSignrank(values, tail);
    if reason == ""
        reason = "ok";
    end
    includedText = clusterListText(result.includedClusterIDs);
    row = table( ...
        string(result.monkeyID), string(result.animal), string(result.chamber), ...
        string(conditionPair), string(includedText), result.totalClusters, ...
        stats.n, stats.meanDeltaBias, stats.medianDeltaBias, stats.semDeltaBias, ...
        stats.nPositive, stats.nNegative, stats.nZero, string(tail), ...
        pValue, testValid, string(reason), ...
        'VariableNames', {'monkeyID', 'animal', 'chamber', 'conditionPair', ...
        'includedClusterIDs', 'totalClusters', 'nRetained', ...
        'meanDeltaBias', 'medianDeltaBias', 'semDeltaBias', ...
        'nPositive', 'nNegative', 'nZero', 'testTail', ...
        'postSelectionP', 'testValid', 'testNote'});
end

function figHandle = plotHistogramFigure(chamberResults, binEdges, panelStats, opts)
    binWidth = binEdges(2) - binEdges(1);
    xDisplay = [binEdges(1) - binWidth ./ 2, binEdges(end) + binWidth ./ 2];
    maxProp = 0;
    for chamberIdx = 1:numel(chamberResults)
        maxProp = max(maxProp, maxProportion(chamberResults(chamberIdx).retainedValues090, binEdges));
        maxProp = max(maxProp, maxProportion(chamberResults(chamberIdx).retainedValuesControl, binEdges));
    end
    yMax = min(1, max(0.12, maxProp + 0.10));

    figHandle = figure('Color', 'w', 'Position', [100 100 1450 430]);
    tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
    legendHandles = gobjects(0);

    for chamberIdx = 1:numel(chamberResults)
        result = chamberResults(chamberIdx);
        ax = nexttile;
        hold(ax, 'on');
        if strcmp(result.monkeyID, 'M2') && strcmp(result.chamber, 'R') && ...
                ~isempty(result.retainedValuesControl)
            h = drawGroupedBars(ax, result, binEdges, opts);
            legendHandles = h;
        else
            h = drawSingleBars(ax, result.retainedValues090, binEdges, ...
                opts.purple090, darken(opts.purple090));
            if isempty(legendHandles)
                legendHandles = h;
            end
        end
        xline(ax, 0, '--', 'Color', [0.45 0.45 0.45], ...
            'LineWidth', 1.2, 'HandleVisibility', 'off');
        xline(ax, mean(result.retainedValues090, 'omitnan'), '-', ...
            'Color', darken(opts.purple090), 'LineWidth', 1.8, ...
            'HandleVisibility', 'off');

        if opts.showRawMarkers
            drawRawMarkers(ax, result.retainedValues090, yMax, opts.purple090);
            if ~isempty(result.retainedValuesControl)
                drawRawMarkers(ax, result.retainedValuesControl, yMax, opts.orangeControl);
            end
        end

        title(ax, sprintf('%s %s %s\nn = %d | clusters %s / %d', ...
            result.monkeyID, char(8212), result.chamber, ...
            numel(result.retainedValues090), clusterListText(result.includedClusterIDs), ...
            result.totalClusters), 'Interpreter', 'none');
        annotatePanel(ax, result, panelStats, opts);
        xlabel(ax, 'Merged \DeltaBias (% correct)');
        ylabel(ax, 'Proportion of included experiments');
        xlim(ax, xDisplay);
        ylim(ax, [0 yMax]);
        box(ax, 'off');
        set(ax, 'TickDir', 'out', 'LineWidth', 1);
    end

    if numel(legendHandles) >= 2
        legend(legendHandles, {'0/90', '45/135'}, 'Location', 'northeastoutside');
    end
    sgtitle('Merged experiment-wise \DeltaBias in positive-effect power clusters');
end

function h = drawSingleBars(ax, values, binEdges, faceColor, edgeColor)
    [counts, edges] = histcounts(values, binEdges);
    assert(sum(counts) == numel(values));
    proportions = counts ./ numel(values);
    assert(abs(sum(proportions) - 1) < 1e-12);
    centers = edges(1:end-1) + diff(edges) ./ 2;
    h = bar(ax, centers, proportions, 1.0, ...
        'FaceColor', faceColor, 'FaceAlpha', 1, ...
        'EdgeColor', edgeColor, 'LineWidth', 1.0);
end

function h = drawGroupedBars(ax, result, binEdges, opts)
    [counts090, edges] = histcounts(result.retainedValues090, binEdges);
    [countsControl, ~] = histcounts(result.retainedValuesControl, binEdges);
    assert(sum(counts090) == numel(result.retainedValues090));
    assert(sum(countsControl) == numel(result.retainedValuesControl));
    prop090 = counts090(:) ./ numel(result.retainedValues090);
    propControl = countsControl(:) ./ numel(result.retainedValuesControl);
    assert(abs(sum(prop090) - 1) < 1e-12);
    assert(abs(sum(propControl) - 1) < 1e-12);
    centers = edges(1:end-1)' + diff(edges(:)) ./ 2;
    h = bar(ax, centers, [prop090 propControl], 0.85, 'grouped');
    h(1).FaceColor = opts.purple090;
    h(1).FaceAlpha = 1;
    h(1).EdgeColor = darken(opts.purple090);
    h(2).FaceColor = opts.orangeControl;
    h(2).FaceAlpha = 1;
    h(2).EdgeColor = darken(opts.orangeControl);
end

function maxVal = maxProportion(values, binEdges)
    if isempty(values)
        maxVal = 0;
        return;
    end
    counts = histcounts(values, binEdges);
    maxVal = max(counts ./ numel(values));
end

function annotatePanel(ax, result, panelStats, opts)
    xLimits = xlim(ax);
    yLimits = ylim(ax);
    xText = xLimits(1) + 0.05 .* diff(xLimits);
    yText = yLimits(2) - 0.12 .* diff(yLimits);
    rows = panelStats(panelStats.monkeyID == string(result.monkeyID) & ...
        panelStats.chamber == string(result.chamber), :);
    row090 = rows(rows.conditionPair == "0/90", :);
    if strcmp(result.monkeyID, 'M2') && strcmp(result.chamber, 'R')
        rowControl = rows(rows.conditionPair == "45/135", :);
        text(ax, xText, yText, formatM2AnnotationLine('0/90', row090, 'p+'), ...
            'Color', darken(opts.purple090), 'FontSize', 10, ...
            'FontWeight', 'bold', 'VerticalAlignment', 'top');
        text(ax, xText, yText - 0.10 .* diff(yLimits), ...
            formatM2AnnotationLine('45/135', rowControl, 'p2'), ...
            'Color', darken(opts.orangeControl), 'FontSize', 10, ...
            'FontWeight', 'bold', 'VerticalAlignment', 'top');
    else
        text(ax, xText, yText, sprintf('median = %0.1f pp\np+ %s', ...
            row090.medianDeltaBias, formatPValue(row090.postSelectionP)), ...
            'Color', darken(opts.purple090), 'FontSize', 10, ...
            'FontWeight', 'bold', 'VerticalAlignment', 'top');
    end
end

function txt = formatM2AnnotationLine(label, row, pLabel)
    if isempty(row) || row.nRetained == 0
        txt = sprintf('%s: n=0, median=n/a, %s n/a', label, pLabel);
    else
        txt = sprintf('%s: n=%d, median=%0.1f pp, %s %s', ...
            label, row.nRetained, row.medianDeltaBias, pLabel, ...
            formatPValue(row.postSelectionP));
    end
end

function pText = formatPValue(pValue)
    if ~isfinite(pValue)
        pText = 'n/a';
    elseif pValue < 0.001
        pText = '< 0.001';
    else
        pText = sprintf('= %.3g', pValue);
    end
end

function drawRawMarkers(ax, values, yMax, color)
    if isempty(values)
        return;
    end
    [uniqueVals, ~, groupIdx] = unique(values);
    yBase = 0.035 .* yMax;
    yStep = 0.025 .* yMax;
    for idx = 1:numel(uniqueVals)
        members = find(groupIdx == idx);
        for memberIdx = 1:numel(members)
            plot(ax, uniqueVals(idx), yBase + (memberIdx - 1) .* yStep, ...
                's', 'MarkerSize', 4, 'MarkerFaceColor', color, ...
                'MarkerEdgeColor', darken(color), 'HandleVisibility', 'off');
        end
    end
end

function colorOut = darken(colorIn)
    colorOut = max(0, colorIn .* 0.65);
end

function saveOutputs(figHandle, panelStats, histogramBins, clusterStats, opts)
    pdfPath = fullfile(opts.outputDir, 'multiChamberDeltaBiasHistogram_stage1_verified.pdf');
    pngPath = fullfile(opts.outputDir, 'multiChamberDeltaBiasHistogram_stage1_verified.png');
    figPath = fullfile(opts.outputDir, 'multiChamberDeltaBiasHistogram_stage1_verified.fig');
    panelStatsPath = fullfile(opts.outputDir, 'multiChamberDeltaBiasPanelStats.csv');
    binsPath = fullfile(opts.outputDir, 'multiChamberDeltaBiasHistogramBins_verified.csv');
    clusterPath = fullfile(opts.outputDir, 'multiChamberDeltaBiasClusterStats_verified.csv');

    exportgraphics(figHandle, pdfPath, 'ContentType', 'vector', 'BackgroundColor', 'white');
    exportgraphics(figHandle, pngPath, 'Resolution', 300, 'BackgroundColor', 'white');
    savefig(figHandle, figPath);
    writetable(panelStats, panelStatsPath);
    writetable(histogramBins, binsPath);
    writetable(clusterStats, clusterPath);

    fprintf('\nSaved verified Stage 1 outputs:\n');
    fprintf('  %s\n', pdfPath);
    fprintf('  %s\n', pngPath);
    fprintf('  %s\n', figPath);
    fprintf('  %s\n', panelStatsPath);
    fprintf('  %s\n', binsPath);
    fprintf('  %s\n', clusterPath);
end

function validateResults(chamberResults, binEdges, histogramBins)
    for chamberIdx = 1:numel(chamberResults)
        result = chamberResults(chamberIdx);
        retainedRows = result.experimentTable090(result.experimentTable090.includedInHistogram, :);
        excludedRows = result.experimentTable090(~result.experimentTable090.clusterIncluded, :);
        assert(~any(excludedRows.includedInHistogram), ...
            'Excluded-cluster 0/90 experiment was marked for plotting.');
        retainedIDs = string(retainedRows.experimentID);
        assert(numel(unique(retainedIDs)) == numel(retainedIDs), ...
            'A retained 0/90 experiment appears more than once.');
        validateConditionBins(result, '0/90', result.retainedValues090, binEdges, histogramBins);
        if ~isempty(result.retainedValuesControl)
            validateConditionBins(result, '45/135', result.retainedValuesControl, binEdges, histogramBins);
        end
    end
end

function validateConditionBins(result, conditionPair, values, binEdges, histogramBins)
    counts = histcounts(values, binEdges);
    assert(sum(counts) == numel(values), 'Histogram counts do not sum to retained n.');
    if numel(values) > 0
        proportions = counts(:) ./ numel(values);
        assert(abs(sum(proportions) - 1) < 1e-12, ...
            'Histogram proportions do not sum to 1.');
    else
        proportions = NaN(numel(counts), 1);
    end
    binRows = histogramBins(histogramBins.monkeyID == string(result.monkeyID) & ...
        histogramBins.chamber == string(result.chamber) & ...
        histogramBins.conditionPair == string(conditionPair), :);
    assert(isequal(binRows.count(:), counts(:)), ...
        'Saved histogram counts differ from plotted counts.');
    if numel(values) > 0
        assert(max(abs(binRows.proportion(:) - proportions)) < 1e-12, ...
            'Saved histogram proportions differ from plotted proportions.');
    end
end

function printChamberSummary(result)
    fprintf('\n%s-%s summary\n', result.monkeyID, result.chamber);
    fprintf('  clusters found: %s / %d\n', ...
        clusterListText(result.clusterStats.clusterID), result.totalClusters);
    for rowIdx = 1:height(result.clusterStats)
        row = result.clusterStats(rowIdx, :);
        if row.included
            status = 'included';
        else
            status = 'excluded';
        end
        fprintf(['  C%d: n=%d, mean=%0.3f, median=%0.3f, SEM=%0.3f, ' ...
            'p_right=%0.4g, %s (%s)\n'], ...
            row.clusterID, row.n, row.meanDeltaBias, row.medianDeltaBias, ...
            row.semDeltaBias, row.oneSidedP, status, row.inclusionReason);
    end
    fprintf('  retained 0/90 experiments: %d\n', numel(result.retainedValues090));
    if strcmp(result.monkeyID, 'M2') && strcmp(result.chamber, 'R')
        fprintf('  retained 45/135 controls matched by experiment ID: %d\n', ...
            numel(result.retainedValuesControl));
    end
end

function textValue = clusterListText(clusterIDs)
    clusterIDs = clusterIDs(:)';
    if isempty(clusterIDs)
        textValue = 'none';
        return;
    end
    labels = arrayfun(@(x) sprintf('C%d', x), clusterIDs, 'UniformOutput', false);
    textValue = strjoin(labels, ', ');
end
