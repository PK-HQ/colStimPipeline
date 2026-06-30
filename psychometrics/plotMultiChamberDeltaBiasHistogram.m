function results = plotMultiChamberDeltaBiasHistogram(opts)
% Plot merged experiment-wise deltaBias count histograms.
%
% Total bar height is the number of included experiments in each bin.
% Colored fill is the subset marked by the saved experiment-level
% overallSignificant result. Distribution-level signed-rank p-values against
% zero are reported alongside the distribution median and significant n/N.

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
    experimentStats = table();

    fprintf('Verified multi-chamber merged deltaBias count histograms\n');
    fprintf(['Delta values: saved distributionSourceAudit tables. ' ...
        'Experiment significance: saved permutation summaries when available.\n']);

    for chamberIdx = 1:nChambers
        cfg = chamberConfigs(chamberIdx);
        fprintf('\nLoading %s-%s audit: %s\n', cfg.monkeyID, cfg.chamber, cfg.auditPath);
        audit = loadAuditTable(cfg.auditPath);
        printAuditColumns(audit);

        [rows090, rowsControl] = extractMergedDeltaBiasRows(audit, cfg);
        permutationSource = loadExperimentPermutationSummary(cfg, opts);
        chamberResult = summarizeChamberRows( ...
            rows090, rowsControl, cfg, permutationSource, opts);
        chamberResults(chamberIdx) = chamberResult;
        clusterStats = [clusterStats; chamberResult.clusterStats]; %#ok<AGROW>
        experimentStats = [experimentStats; chamberResult.experimentTable090]; %#ok<AGROW>
        if ~isempty(chamberResult.experimentTableControl)
            experimentStats = [experimentStats; ...
                chamberResult.experimentTableControl]; %#ok<AGROW>
        end
        printChamberSummary(chamberResult);
    end

    combinedResult = buildCombinedResult(chamberResults);

    displayedValues = collectDisplayedValues(chamberResults);
    binEdges = commonBinEdges(displayedValues, opts.binWidth);
    histogramBins = buildHistogramBinTable( ...
        chamberResults, combinedResult, binEdges);
    panelStats = buildPanelStatsTable(chamberResults, combinedResult);

    plotFigure = plotHistogramFigure( ...
        chamberResults, combinedResult, binEdges, panelStats, opts);
    saveOutputs(plotFigure, panelStats, histogramBins, ...
        clusterStats, experimentStats, opts);
    validateResults( ...
        chamberResults, combinedResult, binEdges, histogramBins);

    results = struct();
    results.chamberResults = chamberResults;
    results.combinedResult = combinedResult;
    results.clusterStats = clusterStats;
    results.experimentStats = experimentStats;
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
    if ~isfield(opts, 'combinedColor') || isempty(opts.combinedColor)
        opts.combinedColor = [0.20 0.55 0.75];
    end
    if ~isfield(opts, 'significanceAlpha') || isempty(opts.significanceAlpha)
        % Used only as a fallback when a saved overallSignificant field is absent.
        opts.significanceAlpha = 0.05;
    end
    if ~isfield(opts, 'requireCompleteSignificance') || ...
            isempty(opts.requireCompleteSignificance)
        opts.requireCompleteSignificance = false;
    end
end

function configs = defaultChamberConfigs(opts)
    configs = struct( ...
        'animal', {'Chip', 'Chip', 'Pepper'}, ...
        'monkeyID', {'M1', 'M1', 'M2'}, ...
        'chamber', {'L', 'R', 'R'}, ...
        'auditPath', {'', '', ''}, ...
        'summaryDir', {'', '', ''}, ...
        'significancePath', {'', '', ''}, ...
        'preferredSignificanceModel', {'', '', ''});

    psychometricFileNames = { ...
        'statisticsL-psychometrics43.mat', ...
        'statisticsR-psychometrics16.mat', ...
        'statisticsR-psychometrics38.mat'};

    for idx = 1:numel(configs)
        configs(idx).summaryDir = fullfile(opts.summaryRoot, ...
            configs(idx).animal, 'Meta', 'summary');
        configs(idx).auditPath = fullfile(configs(idx).summaryDir, ...
            sprintf('distributionSourceAudit_%s_%s_%s.mat', ...
            configs(idx).animal, configs(idx).chamber, opts.modelName));

        % The saved model field is chamber-prefixed, for example
        % LweibullfreeAllC1 or RweibullfreeAllC1.
        configs(idx).preferredSignificanceModel = sprintf('%s%sC1', ...
            configs(idx).chamber, opts.modelName);

        expectedPath = fullfile(configs(idx).summaryDir, ...
            psychometricFileNames{idx});
        configs(idx).significancePath = ...
            resolvePsychometricStatisticsPath( ...
            expectedPath, configs(idx).summaryDir, ...
            configs(idx).chamber);
    end
end

function sourcePath = resolvePsychometricStatisticsPath( ...
        expectedPath, summaryDir, chamber)

    if isfile(expectedPath)
        sourcePath = expectedPath;
        return;
    end

    filePattern = sprintf('statistics%s-psychometrics*.mat', chamber);
    candidates = dir(fullfile(summaryDir, filePattern));
    candidates = candidates(~[candidates.isdir]);

    if isempty(candidates)
        sourcePath = expectedPath;
        return;
    end

    % Prefer the newest matching psychometric file when the expected
    % block-count-specific filename has changed.
    [~, newestIdx] = max([candidates.datenum]);
    sourcePath = fullfile(summaryDir, candidates(newestIdx).name);
    warning('plotMultiChamberDeltaBiasHistogram:UsingAlternatePsychometricFile', ...
        ['Expected psychometric file was not found:\n  %s\n' ...
        'Using newest %s-chamber match instead:\n  %s'], ...
        expectedPath, chamber, sourcePath);
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
        'retainedSignificant090', false(0, 1), ...
        'retainedSignificanceAvailable090', false(0, 1), ...
        'retainedPermutationP090', [], ...
        'retainedValuesControl', [], ...
        'retainedSignificantControl', false(0, 1), ...
        'retainedSignificanceAvailableControl', false(0, 1), ...
        'retainedPermutationPControl', [], ...
        'significanceSourcePath', '', ...
        'significanceSourceFields', {{}}, ...
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


function source = loadExperimentPermutationSummary(cfg, opts)
    source = struct('table', table(), 'path', '', ...
        'modelFields', {{}}, 'hasConditionPair', false);

    sourcePath = cfg.significancePath;
    if ~isfile(sourcePath)
        warning('plotMultiChamberDeltaBiasHistogram:MissingPsychometricStatisticsFile', ...
            ['Stored experiment-level significance file is missing for %s-%s:\n' ...
            '  %s\nThe total-count outline will still plot, but the ' ...
            'significant-experiment fill is unavailable.'], ...
            cfg.monkeyID, cfg.chamber, sourcePath);
        return;
    end

    variables = whos('-file', sourcePath);
    if ~any(strcmp({variables.name}, 'mdlStruct'))
        warning('plotMultiChamberDeltaBiasHistogram:MissingMdlStruct', ...
            ['%s contains no mdlStruct. The total-count outline will still ' ...
            'plot, but significance fill is unavailable.'], sourcePath);
        return;
    end

    loaded = load(sourcePath, 'mdlStruct');
    if ~isfield(loaded, 'mdlStruct') || ~isstruct(loaded.mdlStruct)
        warning('plotMultiChamberDeltaBiasHistogram:InvalidMdlStruct', ...
            '%s does not contain a usable mdlStruct.', sourcePath);
        return;
    end

    [summaryTable, modelFields] = ...
        collectPermutationSummariesFromMdlStruct( ...
        loaded.mdlStruct, cfg, opts, sourcePath);

    if isempty(summaryTable)
        modelFieldsPresent = fieldnames(loaded.mdlStruct);
        warning('plotMultiChamberDeltaBiasHistogram:NoPermutationSummary', ...
            ['%s was loaded for %s-%s, but none of its mdlStruct fields ' ...
            'contains deltaBiasPermutationExperimentSummary.\n' ...
            'Preferred model: %s\nFields present: %s\n' ...
            'Rerun psycluster for this chamber with the individual ' ...
            'permutation test enabled, then resave this psychometrics file.'], ...
            sourcePath, cfg.monkeyID, cfg.chamber, ...
            cfg.preferredSignificanceModel, ...
            strjoin(modelFieldsPresent, ', '));
        return;
    end

    source.table = summaryTable;
    source.path = sourcePath;
    source.modelFields = modelFields;
    source.hasConditionPair = any(strlength(summaryTable.conditionPair) > 0);

    nAvailable = sum(summaryTable.significanceAvailable);
    nSignificant = sum(summaryTable.significanceAvailable & ...
        summaryTable.significant);
    fprintf(['  experiment significance source: %s\n' ...
        '  permutation model fields: %s\n' ...
        '  experiment p available: %d/%d; p < %0.3g: %d\n'], ...
        sourcePath, strjoin(modelFields, ', '), ...
        nAvailable, height(summaryTable), opts.significanceAlpha, nSignificant);
end

function [summaryTable, modelFieldsUsed] = ...
        collectPermutationSummariesFromMdlStruct( ...
        mdlStruct, cfg, opts, sourcePath)

    summaryTable = table();
    modelFieldsUsed = {};
    modelFields = fieldnames(mdlStruct);

    % Prefer the exact model field known for this chamber. If the saved
    % file differs, inspect every model field rather than assuming the
    % histogram model name and the permutation-source model are identical.
    preferredIdx = find(strcmpi(modelFields, ...
        cfg.preferredSignificanceModel), 1, 'first');
    if isempty(preferredIdx)
        orderedFields = modelFields;
    else
        orderedFields = [{modelFields{preferredIdx}}; ...
            modelFields(setdiff(1:numel(modelFields), preferredIdx))];
    end

    for fieldIdx = 1:numel(orderedFields)
        modelField = orderedFields{fieldIdx};
        modelValue = mdlStruct.(modelField);
        if ~isstruct(modelValue) || ...
                ~isfield(modelValue, ...
                'deltaBiasPermutationExperimentSummary') || ...
                isempty(modelValue.deltaBiasPermutationExperimentSummary) || ...
                ~istable(modelValue.deltaBiasPermutationExperimentSummary)
            continue;
        end

        currentSummary = standardizePermutationSummary( ...
            modelValue.deltaBiasPermutationExperimentSummary, ...
            modelField, sourcePath, opts.significanceAlpha);
        summaryTable = [summaryTable; currentSummary]; %#ok<AGROW>
        modelFieldsUsed{end + 1} = modelField; %#ok<AGROW>
    end

    if isempty(summaryTable)
        return;
    end

    summaryTable.experimentKey = normalizeExperimentKey( ...
        summaryTable.experimentID);
    duplicateKey = summaryTable.experimentKey + "|" + ...
        string(summaryTable.clusterID) + "|" + summaryTable.conditionPair;
    [uniqueKeys, ~, groupIdx] = unique(duplicateKey, 'stable');
    keepRows = false(height(summaryTable), 1);

    for keyIdx = 1:numel(uniqueKeys)
        rows = find(groupIdx == keyIdx);
        pValues = summaryTable.permutationP(rows);
        sigValues = summaryTable.significant(rows);
        finiteP = pValues(isfinite(pValues));
        if numel(unique(sigValues)) > 1 || ...
                (~isempty(finiteP) && max(finiteP) - min(finiteP) > 1e-12)
            error('plotMultiChamberDeltaBiasHistogram:ConflictingPermutationRows', ...
                'Conflicting permutation summaries for experiment key %s.', ...
                uniqueKeys(keyIdx));
        end
        keepRows(rows(1)) = true;
    end

    summaryTable = summaryTable(keepRows, :);
end

function summary = standardizePermutationSummary( ...
        sourceTable, modelField, sourcePath, alpha)

    variableNames = sourceTable.Properties.VariableNames;
    experimentVar = findVariableName(variableNames, ...
        {'experimentID', 'experimentLabel', 'blockLabel'});
    pVar = findVariableName(variableNames, ...
        {'rawOverallTwoSidedP', 'overallTwoSidedP', ...
        'rawOverallP', 'overallP', 'pValue'});
    significantVar = findVariableName(variableNames, ...
        {'overallSignificant', 'experimentSignificant', ...
        'isSignificant', 'significant'});
    clusterVar = findVariableName(variableNames, ...
        {'powerClusterID', 'clusterID', 'cluster'});
    conditionVar = findVariableName(variableNames, ...
        {'conditionPair', 'datasetLabel', 'conditionLabel'});
    modelRowVar = findVariableName(variableNames, ...
        {'modelRow', 'sessionRowIndex', 'rowIndex', 'blockIndex'});

    if isempty(experimentVar)
        error('plotMultiChamberDeltaBiasHistogram:MissingPermutationExperimentID', ...
            ['%s:%s permutation summary has no experiment-ID column. ' ...
            'Columns: %s'], sourcePath, modelField, ...
            strjoin(variableNames, ', '));
    end
    if isempty(pVar) && isempty(significantVar)
        error('plotMultiChamberDeltaBiasHistogram:MissingPermutationResult', ...
            ['%s:%s permutation summary has neither a p-value nor a ' ...
            'significance column. Columns: %s'], sourcePath, modelField, ...
            strjoin(variableNames, ', '));
    end

    nRows = height(sourceTable);
    experimentID = string(sourceTable.(experimentVar));
    clusterID = nan(nRows, 1);
    if ~isempty(clusterVar)
        clusterID = double(sourceTable.(clusterVar));
    end

    permutationP = nan(nRows, 1);
    if ~isempty(pVar)
        permutationP = double(sourceTable.(pVar));
    end

    significanceAvailable = false(nRows, 1);
    significant = false(nRows, 1);

    % The colored subset must reproduce the already-settled experiment-level
    % result. Prefer the saved overallSignificant flag exactly as stored.
    % Only fall back to a raw-p threshold when that field is unavailable.
    if ~isempty(significantVar)
        significant = logical(sourceTable.(significantVar));
        significanceAvailable(:) = true;
    elseif ~isempty(pVar)
        significanceAvailable = isfinite(permutationP);
        significant(significanceAvailable) = ...
            permutationP(significanceAvailable) < alpha;
        warning('plotMultiChamberDeltaBiasHistogram:NoStoredOverallSignificant', ...
            ['%s:%s has no overallSignificant field. Falling back to ' ...
            'experiment permutation p < %0.3g.'], ...
            sourcePath, modelField, alpha);
    end

    conditionPair = strings(nRows, 1);
    if ~isempty(conditionVar)
        conditionPair = canonicalConditionPair( ...
            string(sourceTable.(conditionVar)));
    end

    modelRow = nan(nRows, 1);
    if ~isempty(modelRowVar)
        modelRow = double(sourceTable.(modelRowVar));
        modelRow = modelRow(:);
    end

    summary = table( ...
        experimentID, clusterID, conditionPair, modelRow, permutationP, ...
        significant, significanceAvailable, ...
        repmat(string(modelField), nRows, 1), ...
        repmat(string(sourcePath), nRows, 1), ...
        'VariableNames', {'experimentID', 'clusterID', ...
        'conditionPair', 'modelRow', 'permutationP', 'significant', ...
        'significanceAvailable', 'sourceModelField', 'sourcePath'});
end

function variableName = findVariableName(variableNames, candidates)
    variableName = '';
    lowerNames = lower(string(variableNames));
    for candidateIdx = 1:numel(candidates)
        matchIdx = find(lowerNames == lower(string(candidates{candidateIdx})), ...
            1, 'first');
        if ~isempty(matchIdx)
            variableName = variableNames{matchIdx};
            return;
        end
    end
end

function conditionPair = canonicalConditionPair(labels)
    normalized = normalizeLabel(labels);
    conditionPair = strings(size(normalized));
    isControl = contains(normalized, "45") & contains(normalized, "135");
    is090 = contains(normalized, "0") & contains(normalized, "90") & ~isControl;
    conditionPair(is090) = "0/90";
    conditionPair(isControl) = "45/135";
end

function keys = normalizeExperimentKey(experimentIDs)
    % The permutation summary uses blockInfo.label (for example
    % "20240202R6"), whereas distributionSourceAudit may retain the full
    % SVG-style identifier (for example "C1M28D20240202R6"). Reduce both
    % forms to the shared date/run key before joining.
    rawIDs = upper(strtrim(string(experimentIDs)));
    keys = strings(size(rawIDs));

    for idx = 1:numel(rawIDs)
        compactID = regexprep(char(rawIDs(idx)), '[^A-Z0-9]', '');
        tokens = regexp(compactID, 'D?(\d{8})R(\d+)', ...
            'tokens', 'once');

        if isempty(tokens)
            keys(idx) = string(compactID);
        else
            keys(idx) = string([tokens{1}, 'R', tokens{2}]);
        end
    end
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

function result = summarizeChamberRows( ...
        rows090, rowsControl, cfg, permutationSource, opts)

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

        experimentTable090 = [experimentTable090; makeExperimentRows( ...
            cfg, clusterRows, '0/90', included, inclusionReason)]; %#ok<AGROW>
    end

    experimentTable090 = attachExperimentSignificance( ...
        experimentTable090, permutationSource, '0/90', cfg, opts);

    includedClusterIDs = clusterStats.clusterID(clusterStats.included);
    retainedMask090 = experimentTable090.includedInHistogram;
    retainedValues090 = experimentTable090.deltaBiasMerged(retainedMask090);
    retainedSignificant090 = ...
        experimentTable090.experimentSignificant(retainedMask090);
    retainedSignificanceAvailable090 = ...
        experimentTable090.experimentSignificanceAvailable(retainedMask090);
    retainedPermutationP090 = ...
        experimentTable090.experimentPermutationP(retainedMask090);

    experimentTableControl = table();
    retainedValuesControl = [];
    retainedSignificantControl = false(0, 1);
    retainedSignificanceAvailableControl = false(0, 1);
    retainedPermutationPControl = [];
    controlUnmatchedIDs = strings(0, 1);

    if strcmp(cfg.monkeyID, 'M2') && strcmp(cfg.chamber, 'R') && ...
            ~isempty(rowsControl)
        controlInClusters = rowsControl( ...
            ismember(rowsControl.powerClusterID, includedClusterIDs), :);
        retainedIDs090 = string( ...
            experimentTable090.experimentID(retainedMask090));
        matchedControl = ismember( ...
            string(controlInClusters.experimentID), retainedIDs090);
        unmatchedControl = controlInClusters(~matchedControl, :);
        controlUnmatchedIDs = string(unmatchedControl.experimentID);
        if ~isempty(controlUnmatchedIDs)
            fprintf('  M2-R 45/135 unmatched controls excluded: %s\n', ...
                strjoin(controlUnmatchedIDs, ', '));
        end

        matchedRows = controlInClusters(matchedControl, :);
        experimentTableControl = makeExperimentRows( ...
            cfg, matchedRows, '45/135', true, ...
            "matched experiment ID and included 0/90 cluster");
        experimentTableControl = attachExperimentSignificance( ...
            experimentTableControl, permutationSource, ...
            '45/135', cfg, opts);

        retainedValuesControl = ...
            experimentTableControl.deltaBiasMerged;
        retainedSignificantControl = ...
            experimentTableControl.experimentSignificant;
        retainedSignificanceAvailableControl = ...
            experimentTableControl.experimentSignificanceAvailable;
        retainedPermutationPControl = ...
            experimentTableControl.experimentPermutationP;
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
            fprintf(['  C%d: n=%d, mean=%0.3f, median=%0.3f, ' ...
                'p_right=%0.4g, %s (%s)\n'], ...
                row.clusterID, row.n, row.meanDeltaBias, ...
                row.medianDeltaBias, row.oneSidedP, ...
                status, row.inclusionReason);
        end
        c3Rows = clusterStats(clusterStats.clusterID == 3, :);
        if ~isempty(c3Rows) && ~c3Rows.included
            error('plotMultiChamberDeltaBiasHistogram:M2RC3Excluded', ...
                ['M2-R C3 fails the saved-audit inclusion rule: n=%d, ' ...
                'mean=%0.3f, median=%0.3f, p_right=%0.4g. Reason: %s'], ...
                c3Rows.n, c3Rows.meanDeltaBias, ...
                c3Rows.medianDeltaBias, c3Rows.oneSidedP, ...
                c3Rows.inclusionReason);
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
    result.retainedSignificant090 = retainedSignificant090;
    result.retainedSignificanceAvailable090 = ...
        retainedSignificanceAvailable090;
    result.retainedPermutationP090 = retainedPermutationP090;
    result.retainedValuesControl = retainedValuesControl;
    result.retainedSignificantControl = retainedSignificantControl;
    result.retainedSignificanceAvailableControl = ...
        retainedSignificanceAvailableControl;
    result.retainedPermutationPControl = retainedPermutationPControl;
    result.significanceSourcePath = permutationSource.path;
    result.significanceSourceFields = permutationSource.modelFields;
    result.includedClusterIDs = includedClusterIDs;
    result.totalClusters = numel(clusterIDs);
    result.controlUnmatchedIDs = controlUnmatchedIDs;
end

function rowsOut = makeExperimentRows( ...
        cfg, sourceRows, conditionPair, included, reason)

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
        nan(nRows, 1), false(nRows, 1), false(nRows, 1), ...
        strings(nRows, 1), strings(nRows, 1), ...
        'VariableNames', {'monkeyID', 'animal', 'chamber', ...
        'conditionPair', 'experimentID', 'sessionRowIndex', ...
        'blockIndex', 'clusterID', 'deltaBiasMerged', ...
        'clusterIncluded', 'includedInHistogram', ...
        'exclusionReason', 'experimentPermutationP', ...
        'experimentSignificant', ...
        'experimentSignificanceAvailable', ...
        'significanceSourceModelField', ...
        'significanceSourcePath'});
end

function rowsOut = attachExperimentSignificance( ...
        rowsOut, permutationSource, conditionPair, cfg, opts)

    if isempty(rowsOut)
        return;
    end
    if isempty(permutationSource.table)
        reportSignificanceCoverage(rowsOut, cfg, conditionPair, opts);
        return;
    end

    summary = permutationSource.table;
    rowKeys = normalizeExperimentKey(rowsOut.experimentID);
    summaryKeys = normalizeExperimentKey(summary.experimentID);
    requestedPair = string(conditionPair);

    for rowIdx = 1:height(rowsOut)
        % Primary join: canonicalized date/run identifier.
        matchMask = summaryKeys == rowKeys(rowIdx);
        matchMask = restrictPermutationMatches(matchMask, summary, ...
            requestedPair);

        availableMatches = find(matchMask & ...
            summary.significanceAvailable);

        % Fallback join: the permutation context stores modelRow, which is
        % the same local row represented by sessionRowIndex in the saved
        % distribution audit. This is used only when the identifier join
        % produced no available result.
        if isempty(availableMatches) && ...
                ismember('modelRow', summary.Properties.VariableNames) && ...
                isfinite(rowsOut.sessionRowIndex(rowIdx))
            rowMatchMask = isfinite(summary.modelRow) & ...
                summary.modelRow == rowsOut.sessionRowIndex(rowIdx);
            rowMatchMask = restrictPermutationMatches(rowMatchMask, ...
                summary, requestedPair);
            availableMatches = find(rowMatchMask & ...
                summary.significanceAvailable);
        end

        if isempty(availableMatches)
            continue;
        end

        pValues = summary.permutationP(availableMatches);
        sigValues = summary.significant(availableMatches);
        finiteP = pValues(isfinite(pValues));
        if numel(unique(sigValues)) > 1 || ...
                (~isempty(finiteP) && max(finiteP) - min(finiteP) > 1e-12)
            error('plotMultiChamberDeltaBiasHistogram:AmbiguousSignificanceMatch', ...
                ['Conflicting permutation summaries matched %s %s ' ...
                'experiment %s.'], cfg.monkeyID, cfg.chamber, ...
                rowsOut.experimentID(rowIdx));
        end

        chosen = availableMatches(1);
        rowsOut.experimentPermutationP(rowIdx) = ...
            summary.permutationP(chosen);
        rowsOut.experimentSignificant(rowIdx) = ...
            summary.significant(chosen);
        rowsOut.experimentSignificanceAvailable(rowIdx) = true;
        rowsOut.significanceSourceModelField(rowIdx) = ...
            summary.sourceModelField(chosen);
        rowsOut.significanceSourcePath(rowIdx) = ...
            summary.sourcePath(chosen);
    end

    reportSignificanceCoverage(rowsOut, cfg, conditionPair, opts);
end

function matchMask = restrictPermutationMatches( ...
        matchMask, summary, requestedPair)

    % Do not constrain this join by powerClusterID. The saved permutation
    % summary was generated while every session belonged to the original
    % fitting cluster C1. The histogram audit stores the later ordered
    % power-effect cluster assignment, so those cluster labels describe
    % different stages of the analysis and are not valid join keys.
    hasSpecificPair = strlength(summary.conditionPair) > 0;
    if requestedPair == "45/135"
        % Never reuse an unlabeled 0/90 experiment result for the separate
        % 45/135 control distribution.
        matchMask = matchMask & ...
            summary.conditionPair == requestedPair;
    else
        matchMask = matchMask & ...
            (~hasSpecificPair | summary.conditionPair == requestedPair);
    end
end

function reportSignificanceCoverage(rows, cfg, conditionPair, opts)
    includedRows = rows.includedInHistogram;
    nIncluded = sum(includedRows);
    nAvailable = sum( ...
        rows.experimentSignificanceAvailable(includedRows));
    nSignificant = sum( ...
        rows.experimentSignificant(includedRows) & ...
        rows.experimentSignificanceAvailable(includedRows));

    fprintf(['  %s-%s %s experiment permutation coverage: ' ...
        '%d/%d available; %d significant at alpha=%0.3g.\n'], ...
        cfg.monkeyID, cfg.chamber, conditionPair, ...
        nAvailable, nIncluded, nSignificant, opts.significanceAlpha);

    if nAvailable < nIncluded
        missingIDs = string(rows.experimentID( ...
            includedRows & ...
            ~rows.experimentSignificanceAvailable));
        warning('plotMultiChamberDeltaBiasHistogram:PartialSignificanceCoverage', ...
            ['%s-%s %s has experiment-level permutation results for %d/%d ' ...
            'included experiments. Missing IDs: %s'], ...
            cfg.monkeyID, cfg.chamber, conditionPair, ...
            nAvailable, nIncluded, strjoin(missingIDs, ', '));

        if opts.requireCompleteSignificance
            error('plotMultiChamberDeltaBiasHistogram:IncompleteSignificance', ...
                ['Complete experiment-level significance was required, but ' ...
                '%s-%s %s has only %d/%d matched results.'], ...
                cfg.monkeyID, cfg.chamber, conditionPair, ...
                nAvailable, nIncluded);
        end
    end
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

function combinedResult = buildCombinedResult(chamberResults)
    combinedResult = emptyChamberResult();
    combinedResult.animal = 'All';
    combinedResult.monkeyID = 'All';
    combinedResult.chamber = 'Combined';

    for chamberIdx = 1:numel(chamberResults)
        result = chamberResults(chamberIdx);
        combinedResult.retainedValues090 = [ ...
            combinedResult.retainedValues090; ...
            result.retainedValues090(:)]; %#ok<AGROW>
        combinedResult.retainedSignificant090 = [ ...
            combinedResult.retainedSignificant090; ...
            logical(result.retainedSignificant090(:))]; %#ok<AGROW>
        combinedResult.retainedSignificanceAvailable090 = [ ...
            combinedResult.retainedSignificanceAvailable090; ...
            logical(result.retainedSignificanceAvailable090(:))]; %#ok<AGROW>
        combinedResult.retainedPermutationP090 = [ ...
            combinedResult.retainedPermutationP090; ...
            result.retainedPermutationP090(:)]; %#ok<AGROW>

        retainedRows = result.experimentTable090( ...
            result.experimentTable090.includedInHistogram, :);
        combinedResult.experimentTable090 = [ ...
            combinedResult.experimentTable090; ...
            retainedRows]; %#ok<AGROW>
    end

    assert(numel(combinedResult.retainedValues090) == ...
        numel(combinedResult.retainedSignificant090));
    assert(numel(combinedResult.retainedValues090) == ...
        numel(combinedResult.retainedSignificanceAvailable090));

    fprintf(['\nCombined 0/90 panel: %d experiments across all three ' ...
        'animal/chamber datasets; %d experiment p-values available; ' ...
        '%d significant.\n'], ...
        numel(combinedResult.retainedValues090), ...
        sum(combinedResult.retainedSignificanceAvailable090), ...
        sum(combinedResult.retainedSignificanceAvailable090 & ...
            combinedResult.retainedSignificant090));
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

function histogramBins = buildHistogramBinTable( ...
        chamberResults, combinedResult, binEdges)

    histogramBins = table();
    for chamberIdx = 1:numel(chamberResults)
        result = chamberResults(chamberIdx);
        histogramBins = [histogramBins; makeBinRows( ...
            result, '0/90', result.retainedValues090, ...
            result.retainedSignificant090, ...
            result.retainedSignificanceAvailable090, ...
            binEdges)]; %#ok<AGROW>
        if ~isempty(result.retainedValuesControl)
            histogramBins = [histogramBins; makeBinRows( ...
                result, '45/135', result.retainedValuesControl, ...
                result.retainedSignificantControl, ...
                result.retainedSignificanceAvailableControl, ...
                binEdges)]; %#ok<AGROW>
        end
    end

    histogramBins = [histogramBins; makeBinRows( ...
        combinedResult, '0/90', combinedResult.retainedValues090, ...
        combinedResult.retainedSignificant090, ...
        combinedResult.retainedSignificanceAvailable090, ...
        binEdges)]; %#ok<AGROW>
end

function rows = makeBinRows( ...
        result, conditionPair, values, significantMask, ...
        significanceAvailable, binEdges)

    values = values(:);
    significantMask = logical(significantMask(:));
    significanceAvailable = logical(significanceAvailable(:));
    assert(numel(values) == numel(significantMask));
    assert(numel(values) == numel(significanceAvailable));

    [counts, edges, binIndex] = histcounts(values, binEdges);
    assert(sum(counts) == numel(values), ...
        'Histogram counts lost values.');

    nBins = numel(counts);
    significantCounts = zeros(nBins, 1);
    testedCounts = zeros(nBins, 1);
    for binIdx = 1:nBins
        inBin = binIndex == binIdx;
        testedCounts(binIdx) = sum( ...
            inBin & significanceAvailable);
        significantCounts(binIdx) = sum( ...
            inBin & significanceAvailable & significantMask);
    end
    assert(all(significantCounts <= testedCounts));
    assert(all(testedCounts <= counts(:)));

    countColumn = counts(:);
    significantFractionOfAll = nan(nBins, 1);
    nonemptyBins = countColumn > 0;
    significantFractionOfAll(nonemptyBins) = ...
        significantCounts(nonemptyBins) ./ countColumn(nonemptyBins);

    significantFractionOfTested = nan(nBins, 1);
    testedBins = testedCounts > 0;
    significantFractionOfTested(testedBins) = ...
        significantCounts(testedBins) ./ testedCounts(testedBins);

    centers = edges(1:end-1)' + diff(edges(:)) ./ 2;
    rows = table( ...
        repmat(string(result.monkeyID), nBins, 1), ...
        repmat(string(result.animal), nBins, 1), ...
        repmat(string(result.chamber), nBins, 1), ...
        repmat(string(conditionPair), nBins, 1), ...
        (1:nBins)', edges(1:end-1)', edges(2:end)', centers, ...
        counts(:), testedCounts, significantCounts, ...
        significantFractionOfAll, significantFractionOfTested, ...
        'VariableNames', {'monkeyID', 'animal', 'chamber', ...
        'conditionPair', 'binIndex', 'binLeft', 'binRight', ...
        'binCenter', 'count', 'testedCount', ...
        'significantCount', 'significantFractionOfAll', ...
        'significantFractionOfTested'});
end

function panelStats = buildPanelStatsTable( ...
        chamberResults, combinedResult)

    panelStats = table();
    for chamberIdx = 1:numel(chamberResults)
        result = chamberResults(chamberIdx);
        panelStats = [panelStats; makePanelStatsRow( ...
            result, '0/90', result.retainedValues090, ...
            result.retainedSignificant090, ...
            result.retainedSignificanceAvailable090, ...
            'right')]; %#ok<AGROW>
    end

    panelStats = [panelStats; makePanelStatsRow( ...
        combinedResult, '0/90', combinedResult.retainedValues090, ...
        combinedResult.retainedSignificant090, ...
        combinedResult.retainedSignificanceAvailable090, ...
        'right')]; %#ok<AGROW>
end

function row = makePanelStatsRow( ...
        result, conditionPair, values, significantMask, ...
        significanceAvailable, tail)

    stats = valueSummaryStats(values);
    [pValue, reason, testValid] = runSignrank(values, tail);
    if reason == ""
        reason = "ok";
    end

    significantMask = logical(significantMask(:));
    significanceAvailable = logical(significanceAvailable(:));
    nExperimentPAvailable = sum(significanceAvailable);
    nExperimentSignificant = sum( ...
        significanceAvailable & significantMask);
    if nExperimentPAvailable > 0
        fractionExperimentSignificant = ...
            nExperimentSignificant ./ nExperimentPAvailable;
    else
        fractionExperimentSignificant = NaN;
    end

    includedText = clusterListText(result.includedClusterIDs);
    row = table( ...
        string(result.monkeyID), string(result.animal), ...
        string(result.chamber), string(conditionPair), ...
        string(includedText), result.totalClusters, ...
        stats.n, stats.meanDeltaBias, stats.medianDeltaBias, ...
        stats.semDeltaBias, stats.nPositive, stats.nNegative, ...
        stats.nZero, string(tail), pValue, testValid, ...
        string(reason), nExperimentPAvailable, ...
        nExperimentSignificant, fractionExperimentSignificant, ...
        'VariableNames', {'monkeyID', 'animal', 'chamber', ...
        'conditionPair', 'includedClusterIDs', 'totalClusters', ...
        'nRetained', 'meanDeltaBias', 'medianDeltaBias', ...
        'semDeltaBias', 'nPositive', 'nNegative', 'nZero', ...
        'testTail', 'postSelectionP', 'testValid', 'testNote', ...
        'nExperimentPAvailable', 'nExperimentSignificant', ...
        'fractionExperimentSignificant'});
end

function figHandle = plotHistogramFigure( ...
        chamberResults, combinedResult, binEdges, panelStats, opts)

    binWidth = binEdges(2) - binEdges(1);
    xDisplay = [binEdges(1) - binWidth ./ 2, ...
        binEdges(end) + binWidth ./ 2];

    allResults = [chamberResults(:); combinedResult];

    figHandle = figure('Color', 'w', ...
        'Position', [60 100 1900 560]);
    tiledlayout(1, 4, 'Padding', 'compact', ...
        'TileSpacing', 'compact');

    legendHandles = gobjects(0);
    legendLabels = {};
    legendAxis = gobjects(1);
    allAnnotationHandles = gobjects(0);

    for panelIdx = 1:numel(allResults)
        result = allResults(panelIdx);
        isCombinedPanel = panelIdx > numel(chamberResults);
        ax = nexttile;
        hold(ax, 'on');

        isControlPanel = ~isCombinedPanel && ...
            strcmp(result.monkeyID, 'M2') && ...
            strcmp(result.chamber, 'R') && ...
            ~isempty(result.retainedValuesControl);

        if isCombinedPanel
            panelColor = opts.combinedColor;
            barHandles = drawSingleCountBars( ...
                ax, result.retainedValues090, ...
                result.retainedSignificant090, ...
                result.retainedSignificanceAvailable090, ...
                binEdges, panelColor);
        elseif isControlPanel
            panelColor = opts.purple090;
            barHandles = drawGroupedCountBars( ...
                ax, result, binEdges, opts);
        else
            panelColor = opts.purple090;
            barHandles = drawSingleCountBars( ...
                ax, result.retainedValues090, ...
                result.retainedSignificant090, ...
                result.retainedSignificanceAvailable090, ...
                binEdges, panelColor);
        end

        if isempty(legendHandles)
            legendAxis = ax;
            legendHandles = [barHandles.total(1), ...
                barHandles.significant(1)];
            legendLabels = {'All experiments', ...
                'Experiment significant'};
        end

        xline(ax, 0, '--', 'Color', [0.45 0.45 0.45], ...
            'LineWidth', 1.2, 'HandleVisibility', 'off');

        maximumCount = maxHistogramCount( ...
            result.retainedValues090, binEdges);
        if isControlPanel
            maximumCount = max(maximumCount, maxHistogramCount( ...
                result.retainedValuesControl, binEdges));
        end
        yHeadroom = max(2, ceil(0.35 .* maximumCount));
        yMax = max(3, maximumCount + yHeadroom);

        xlim(ax, xDisplay);
        ylim(ax, [0 yMax]);

        drawMedianArrow(ax, result.retainedValues090, ...
            0.89 .* yMax, panelColor);

        if opts.showRawMarkers
            drawRawMarkers(ax, result.retainedValues090, ...
                yMax, panelColor);
            if isControlPanel
                drawRawMarkers(ax, result.retainedValuesControl, ...
                    yMax, opts.orangeControl);
            end
        end

        if isCombinedPanel
            title(ax, sprintf('All animals combined\nn = %d', ...
                numel(result.retainedValues090)), ...
                'Interpreter', 'none');
        else
            title(ax, sprintf('%s %s %s\nn = %d | clusters %s / %d', ...
                result.monkeyID, char(8212), result.chamber, ...
                numel(result.retainedValues090), ...
                clusterListText(result.includedClusterIDs), ...
                result.totalClusters), 'Interpreter', 'none');
        end

        annotationHandles = annotatePanel(ax, result, panelStats, opts, panelColor);
        allAnnotationHandles = [allAnnotationHandles; annotationHandles(:)]; %#ok<AGROW>
        xlabel(ax, 'Merged \DeltaBias (% correct)');
        ylabel(ax, 'Experiment count');
        axis(ax, 'square');
        box(ax, 'off');
        set(ax, 'TickDir', 'out', 'LineWidth', 1);

        % Apply the standard panel formatting, then restore the within-panel
        % statistical annotations to the requested fixed 10-point size.
        axes(ax);
        upFontSize(14, .02);
        if ~isempty(allAnnotationHandles)
            validAnnotationHandles = ...
                allAnnotationHandles(isgraphics(allAnnotationHandles));
            if ~isempty(validAnnotationHandles)
                set(validAnnotationHandles, ...
                    'FontSize', 10, 'Color', 'k');
            end
        end
    end

    if ~isempty(legendHandles) && isgraphics(legendAxis)
        legend(legendAxis, legendHandles, legendLabels, ...
            'Location', 'northwest');
    end
    sgtitle(['Merged experiment-wise \DeltaBias in ' ...
        'positive-effect power clusters']);
end

function h = drawSingleCountBars( ...
        ax, values, significantMask, significanceAvailable, ...
        binEdges, faceColor)

    values = values(:);
    significantMask = logical(significantMask(:));
    significanceAvailable = logical(significanceAvailable(:));
    totalCounts = histcounts(values, binEdges);
    significantCounts = histcounts( ...
        values(significanceAvailable & significantMask), binEdges);
    centers = binEdges(1:end-1) + diff(binEdges) ./ 2;

    h.significant = bar(ax, centers, significantCounts, 1.0, ...
        'FaceColor', faceColor, ...
        'FaceAlpha', 0.85, ...
        'EdgeColor', 'none', ...
        'LineWidth', 1.0);

    h.total = bar(ax, centers, totalCounts, 1.0, ...
        'FaceColor', 'none', ...
        'EdgeColor', 'k', ...
        'LineWidth', 1.5);

    assert(all(significantCounts <= totalCounts));
end

function h = drawGroupedCountBars(ax, result, binEdges, opts)
    total090 = histcounts(result.retainedValues090, binEdges);
    totalControl = histcounts( ...
        result.retainedValuesControl, binEdges);
    significant090 = histcounts( ...
        result.retainedValues090( ...
        result.retainedSignificanceAvailable090 & ...
        result.retainedSignificant090), binEdges);
    significantControl = histcounts( ...
        result.retainedValuesControl( ...
        result.retainedSignificanceAvailableControl & ...
        result.retainedSignificantControl), binEdges);

    centers = binEdges(1:end-1)' + diff(binEdges(:)) ./ 2;

    h.significant = bar(ax, centers, ...
        [significant090(:), significantControl(:)], ...
        0.85, 'grouped');
    h.significant(1).FaceColor = opts.purple090;
    h.significant(1).FaceAlpha = 0.85;
    h.significant(1).EdgeColor = 'none';
    h.significant(2).FaceColor = opts.orangeControl;
    h.significant(2).FaceAlpha = 0.85;
    h.significant(2).EdgeColor = 'none';

    h.total = bar(ax, centers, ...
        [total090(:), totalControl(:)], ...
        0.85, 'grouped');
    for handleIdx = 1:numel(h.total)
        h.total(handleIdx).FaceColor = 'none';
        h.total(handleIdx).EdgeColor = 'k';
        h.total(handleIdx).LineWidth = 1.5;
    end

    assert(all(significant090 <= total090));
    assert(all(significantControl <= totalControl));
end

function maxVal = maxHistogramCount(values, binEdges)
    if isempty(values)
        maxVal = 0;
        return;
    end
    counts = histcounts(values, binEdges);
    maxVal = max(counts);
end

function drawMedianArrow(ax, values, yPosition, color)
    values = values(isfinite(values));
    if isempty(values)
        return;
    end
    medianValue = median(values);
    plot(ax, medianValue, yPosition, 'v', ...
        'MarkerSize', 9, ...
        'MarkerFaceColor', color, ...
        'MarkerEdgeColor', 'k', ...
        'LineWidth', 1.0, ...
        'HandleVisibility', 'off');
end

function annotationHandles = annotatePanel(ax, result, panelStats, opts, primaryColor) %#ok<INUSD>
    annotationHandles = gobjects(0);
    rows = panelStats( ...
        panelStats.monkeyID == string(result.monkeyID) & ...
        panelStats.chamber == string(result.chamber), :);
    row090 = rows(rows.conditionPair == "0/90", :);

    % Only the primary experimental distribution is annotated. The 45/135
    % control may remain plotted in M2-R, but it has no separate stats text.
    annotationHandles(end + 1) = text(ax, 0.97, 0.97, ...
        formatPanelAnnotationLine(row090), ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'top', ...
        'Color', 'k', ...
        'FontSize', 10, ...
        'FontWeight', 'bold', ...
        'Interpreter', 'tex', ...
        'HandleVisibility', 'off');

    if row090.nExperimentPAvailable < row090.nRetained
        annotationHandles(end + 1) = text(ax, 0.97, 0.04, ...
            sprintf('Experiment p available: %d/%d', ...
            row090.nExperimentPAvailable, row090.nRetained), ...
            'Units', 'normalized', ...
            'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'bottom', ...
            'Color', 'k', ...
            'FontSize', 10, ...
            'Interpreter', 'none', ...
            'HandleVisibility', 'off');
    end
end

function txt = formatPanelAnnotationLine(row)
    if isempty(row) || row.nRetained == 0
        txt = sprintf('median = n/a\np n/a (n=0/0)');
        return;
    end

    if row.nExperimentPAvailable > 0
        nText = sprintf('n=%d/%d', ...
            row.nExperimentSignificant, ...
            row.nExperimentPAvailable);
    else
        nText = 'n unavailable';
    end

    txt = sprintf('median = %0.1f\np %s (%s)', ...
        row.medianDeltaBias, ...
        formatPValue(row.postSelectionP), nText);
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

function saveOutputs( ...
        figHandle, panelStats, histogramBins, ...
        clusterStats, experimentStats, opts)

    pdfPath = fullfile(opts.outputDir, ...
        'multiChamberDeltaBiasHistogram_counts.pdf');
    pngPath = fullfile(opts.outputDir, ...
        'multiChamberDeltaBiasHistogram_counts.png');
    figPath = fullfile(opts.outputDir, ...
        'multiChamberDeltaBiasHistogram_counts.fig');
    panelStatsPath = fullfile(opts.outputDir, ...
        'multiChamberDeltaBiasPanelStats.csv');
    binsPath = fullfile(opts.outputDir, ...
        'multiChamberDeltaBiasHistogramBins_counts.csv');
    clusterPath = fullfile(opts.outputDir, ...
        'multiChamberDeltaBiasClusterStats.csv');
    experimentPath = fullfile(opts.outputDir, ...
        'multiChamberDeltaBiasExperimentSignificance.csv');

    exportgraphics(figHandle, pdfPath, ...
        'ContentType', 'vector', 'BackgroundColor', 'white');
    exportgraphics(figHandle, pngPath, ...
        'Resolution', 300, 'BackgroundColor', 'white');
    savefig(figHandle, figPath);
    writetable(panelStats, panelStatsPath);
    writetable(histogramBins, binsPath);
    writetable(clusterStats, clusterPath);
    writetable(experimentStats, experimentPath);

    fprintf('\nSaved count-histogram outputs:\n');
    fprintf('  %s\n', pdfPath);
    fprintf('  %s\n', pngPath);
    fprintf('  %s\n', figPath);
    fprintf('  %s\n', panelStatsPath);
    fprintf('  %s\n', binsPath);
    fprintf('  %s\n', clusterPath);
    fprintf('  %s\n', experimentPath);
end

function validateResults( ...
        chamberResults, combinedResult, binEdges, histogramBins)

    for chamberIdx = 1:numel(chamberResults)
        result = chamberResults(chamberIdx);
        retainedRows = result.experimentTable090( ...
            result.experimentTable090.includedInHistogram, :);
        excludedRows = result.experimentTable090( ...
            ~result.experimentTable090.clusterIncluded, :);
        assert(~any(excludedRows.includedInHistogram), ...
            'Excluded-cluster 0/90 experiment was marked for plotting.');
        retainedIDs = string(retainedRows.experimentID);
        assert(numel(unique(retainedIDs)) == numel(retainedIDs), ...
            'A retained 0/90 experiment appears more than once.');

        validateConditionBins( ...
            result, '0/90', result.retainedValues090, ...
            result.retainedSignificant090, ...
            result.retainedSignificanceAvailable090, ...
            binEdges, histogramBins);

        if ~isempty(result.retainedValuesControl)
            validateConditionBins( ...
                result, '45/135', result.retainedValuesControl, ...
                result.retainedSignificantControl, ...
                result.retainedSignificanceAvailableControl, ...
                binEdges, histogramBins);
        end
    end

    validateConditionBins( ...
        combinedResult, '0/90', combinedResult.retainedValues090, ...
        combinedResult.retainedSignificant090, ...
        combinedResult.retainedSignificanceAvailable090, ...
        binEdges, histogramBins);
end

function validateConditionBins( ...
        result, conditionPair, values, significantMask, ...
        significanceAvailable, binEdges, histogramBins)

    values = values(:);
    significantMask = logical(significantMask(:));
    significanceAvailable = logical(significanceAvailable(:));
    [counts, ~, binIndex] = histcounts(values, binEdges);
    assert(sum(counts) == numel(values), ...
        'Histogram counts do not sum to retained n.');

    nBins = numel(counts);
    testedCounts = zeros(nBins, 1);
    significantCounts = zeros(nBins, 1);
    for binIdx = 1:nBins
        inBin = binIndex == binIdx;
        testedCounts(binIdx) = sum( ...
            inBin & significanceAvailable);
        significantCounts(binIdx) = sum( ...
            inBin & significanceAvailable & significantMask);
    end

    binRows = histogramBins( ...
        histogramBins.monkeyID == string(result.monkeyID) & ...
        histogramBins.chamber == string(result.chamber) & ...
        histogramBins.conditionPair == string(conditionPair), :);

    assert(isequal(binRows.count(:), counts(:)), ...
        'Saved histogram counts differ from plotted counts.');
    assert(isequal(binRows.testedCount(:), testedCounts(:)), ...
        'Saved tested counts differ from plotted counts.');
    assert(isequal(binRows.significantCount(:), ...
        significantCounts(:)), ...
        'Saved significant counts differ from plotted counts.');
    assert(all(binRows.significantCount <= binRows.testedCount));
    assert(all(binRows.testedCount <= binRows.count));
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
    fprintf(['  retained 0/90 experiments: %d | permutation p available: ' ...
        '%d | significant: %d\n'], ...
        numel(result.retainedValues090), ...
        sum(result.retainedSignificanceAvailable090), ...
        sum(result.retainedSignificanceAvailable090 & ...
            result.retainedSignificant090));
    if strcmp(result.monkeyID, 'M2') && strcmp(result.chamber, 'R')
        fprintf(['  retained 45/135 controls matched by experiment ID: %d | ' ...
            'permutation p available: %d | significant: %d\n'], ...
            numel(result.retainedValuesControl), ...
            sum(result.retainedSignificanceAvailableControl), ...
            sum(result.retainedSignificanceAvailableControl & ...
                result.retainedSignificantControl));
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
