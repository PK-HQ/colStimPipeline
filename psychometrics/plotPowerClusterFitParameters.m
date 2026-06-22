function [figureHandles, distributionData, sourceAudit, statsAudit] = ...
        plotPowerClusterFitParameters(mdl, aggregateFits, opts)
% Plot individual-session Weibull parameters for horizontal, vertical, and merged fits.

    if nargin < 3 || isempty(opts)
        opts = struct();
    end
    opts = applyDefaults(opts);

    if isfield(mdl, 'signedBX0') && ...
            isfield(opts, 'modelFieldName') && ...
            contains(char(opts.modelFieldName), 'weibullSignedBX0')
        [figureHandles, distributionData, sourceAudit, statsAudit] = ...
            plotSignedBX0ParameterDistribution(mdl, aggregateFits, opts);
        return;
    end

    if isfield(mdl, 'signedX0') && ...
            isfield(opts, 'modelFieldName') && ...
            contains(char(opts.modelFieldName), 'weibullSignedX0')
        [figureHandles, distributionData, sourceAudit, statsAudit] = ...
            plotSignedX0HorizontalDistribution(mdl, aggregateFits, opts);
        return;
    end

    viewNames = {'horizontal', 'vertical', 'merged'};
    viewTitles = {'Horizontal visual stimulus', ...
        'Vertical visual stimulus', 'Merged'};
    fitFields = {'fittedParamsHorizontal', 'fittedParamsVertical', 'fittedParams'};
    headerFields = {'headersHorizontal', 'headersVertical', 'headers'};
    parameterLabels = {'A (%)', 'B (%)', '\alpha', '\beta'};
    conditionLabels = {'Baseline', 'Con-Opto', 'Incon-Opto'};
    conditionColors = [0 0 0; 0.9294 0.1098 0.1373; 0 0.0941 0.6627];
    conditionPairs = [1 2; 2 3; 1 3];
    parameterStatsOpts = performanceStatsOptions(opts);

    sources = cell(1, 3);
    fprintf('\n=== Individual Weibull parameter source audit ===\n');
    for viewIdx = 1:3
        sources{viewIdx} = auditFitSource(mdl, fitFields{viewIdx}, ...
            headerFields{viewIdx}, opts.modelFieldName);
        printSourceAudit(sources{viewIdx}, viewTitles{viewIdx});
    end

    figureHandles = gobjects(numel(aggregateFits), 1);
    distributionData = repmat(struct(), numel(aggregateFits), 1);
    sourceAudit = emptyDistributionSourceAuditTable();
    statsAudit = emptyDistributionStatsAuditTable();
    if strcmp(opts.statsOpts.testDirection, 'one-sided')
        fprintf(['Parameter distribution statistics: planned one-sided ' ...
            'mode is ON. Holm correction remains within each ' ...
            'three-test parameter panel.\n']);
    end
    for aggregateIdx = 1:numel(aggregateFits)
        rows = unique(aggregateFits(aggregateIdx).mdlRowIndices(:)', 'stable');
        sourceBlocks = getSourceBlocks(aggregateFits(aggregateIdx), rows);
        experimentIDs = lookupExperimentIDs(sourceBlocks, opts);
        baselineModes = lookupBaselineModes(sourceBlocks, opts);
        clusterTitle = parameterTitle(aggregateFits(aggregateIdx), numel(rows));

        parameterValues = cell(1, 3);
        clusterAudits = cell(1, 3);
        fprintf('\n[%s]\n', clusterTitle);
        fprintf('  Sessions assigned to cluster: %d\n', numel(rows));
        for viewIdx = 1:3
            [parameterValues{viewIdx}, clusterAudits{viewIdx}] = ...
                extractClusterParameters(sources{viewIdx}, rows, ...
                sourceBlocks, clusterTitle, viewTitles{viewIdx});
            fprintf('  %s valid fits: %d of %d sessions\n', ...
                viewTitles{viewIdx}, clusterAudits{viewIdx}.nValidFits, ...
                numel(rows));
            printExclusions(clusterAudits{viewIdx}, ...
                clusterTitle, viewTitles{viewIdx});
        end

        sharedYLimits = computeSharedYLimits( ...
            parameterValues, parameterLabels, clusterTitle);

        figureHandles(aggregateIdx) = figure( ...
            'Name', [clusterTitle ' parameter distributions'], ...
            'Color', 'w', ...
            'Visible', opts.figureVisible);
        makeSubplot = @(position) subtightplot(4, 3, position, ...
            [0.09 0.055], [0.11 0.09], [0.09 0.05]);

        distributionData(aggregateIdx).clusterID = ...
            aggregateFits(aggregateIdx).clusterID;
        distributionData(aggregateIdx).mdlRowIndices = rows;
        distributionData(aggregateIdx).sourceBlockIndices = sourceBlocks;
        parameterAxes = gobjects(4, 3);
        availablePanels = false(4, 3);
        parameterAnnotations = cell(4, 3);
        parameterTestsByPanel = cell(4, 3);
        parameterCounts = cell(4, 3);
        parameterStatsByPanel = cell(4, 3);
        plottedValuesByPanel = cell(4, 3);
        referenceValues = nan(4, 3);

        for parameterIdx = 1:4
            for viewIdx = 1:3
                ax = makeSubplot((parameterIdx - 1) .* 3 + viewIdx);
                parameterAxes(parameterIdx, viewIdx) = ax;
                if parameterIdx == 1
                    title(ax, viewTitles{viewIdx}, ...
                        'FontWeight', 'normal', 'FontSize', 15);
                end

                viewName = viewNames{viewIdx};
                distributionData(aggregateIdx).(viewName).sourceAudit = ...
                    clusterAudits{viewIdx};
                distributionData(aggregateIdx).(viewName).values = ...
                    parameterValues{viewIdx};

                if isempty(parameterValues{viewIdx})
                    markUnavailablePanel(ax, parameterLabels{parameterIdx});
                    continue;
                end

                values = reshape( ...
                    parameterValues{viewIdx}(:, :, parameterIdx), ...
                    size(parameterValues{viewIdx}, 1), 3);
                if parameterIdx <= 2
                    values = 100 .* values;
                end

                stats = summarizeColumns(values);
                parameterStatsByPanel{parameterIdx, viewIdx} = stats;
                plottedValuesByPanel{parameterIdx, viewIdx} = values;
                distributionData(aggregateIdx).(viewName).stats(parameterIdx) = stats;
                verifyPointCounts(stats.n, clusterAudits{viewIdx}.nValidFits, ...
                    clusterTitle, viewTitles{viewIdx}, ...
                    parameterLabels{parameterIdx}, conditionLabels);
                fprintf('  %s %s dots: baseline=%d, con=%d, incon=%d\n', ...
                    viewTitles{viewIdx}, parameterLabels{parameterIdx}, ...
                    stats.n(1), stats.n(2), stats.n(3));
                printParameterMeans(viewTitles{viewIdx}, ...
                    parameterLabels{parameterIdx}, stats);
                auditValues = reshape( ...
                    clusterAudits{viewIdx}.reconstructedValues(:, :, parameterIdx), ...
                    clusterAudits{viewIdx}.nRequestedSessions, 3);
                if parameterIdx <= 2
                    auditValues = 100 .* auditValues;
                end
                auditReasons = repmat( ...
                    string(clusterAudits{viewIdx}.requestedReasons(:)), 1, 3);
                sourceColumns = parameterSourceColumns( ...
                    sources{viewIdx}.headerMap, parameterIdx);
                parameterAudit = buildDistributionSourceAuditRows( ...
                    aggregateFits(aggregateIdx), 'parameterDistribution', ...
                    viewName, parameterAuditRowLabel(parameterIdx), ...
                    {'baseline', 'con', 'incon'}, ...
                    clusterAudits{viewIdx}.requestedMdlRows, ...
                    clusterAudits{viewIdx}.requestedSourceBlocks, ...
                    experimentIDs, ['mdl.' sources{viewIdx}.fitField], ...
                    sourceColumns, auditValues, isfinite(auditValues), ...
                    auditReasons, baselineModes);
                assertParameterAuditCounts(parameterAudit, stats.n, ...
                    clusterTitle, viewTitles{viewIdx}, ...
                    parameterLabels{parameterIdx}, ...
                    clusterAudits{viewIdx}.nValidFits);
                sourceAudit = [sourceAudit; parameterAudit]; %#ok<AGROW>
                printSourceAuditSummary(parameterAudit);

                referenceValue = parameterReferenceValue(parameterIdx, values);
                yline(ax, referenceValue, '--', ...
                    'Color', 0.4 .* [1 1 1], ...
                    'LineWidth', 1.2, 'HandleVisibility', 'off');
                conditionPlotOpts = conditionMarkerOptions( ...
                    opts, conditionColors);
                plotJitteredGroupWithSummary(ax, values, ...
                    conditionColors, conditionLabels, conditionPlotOpts);
                ylim(ax, sharedYLimits(parameterIdx, :));
                styleParameterAxes(ax, conditionLabels);
                parameterTests = pairedPanelStats( ...
                    values, conditionPairs, parameterStatsOpts);
                distributionData(aggregateIdx).(viewName).parameterTests{parameterIdx} = ...
                    parameterTests;
                parameterTestsByPanel{parameterIdx, viewIdx} = parameterTests;
                parameterCounts{parameterIdx, viewIdx} = stats.n;
                referenceValues(parameterIdx, viewIdx) = referenceValue;
                availablePanels(parameterIdx, viewIdx) = true;
                ylabel(ax, parameterLabels{parameterIdx}, 'FontSize', 15);
                if opts.showSourceAuditText
                    addSourceAuditText(ax, ...
                        ['mdl.' sources{viewIdx}.fitField], stats.n);
                end
            end
        end

        for parameterIdx = 1:4
            viewIndices = find(availablePanels(parameterIdx, :));
            if isempty(viewIndices)
                continue;
            end
            rowBounds = sharedParameterAnnotationBounds( ...
                plottedValuesByPanel(parameterIdx, viewIndices), ...
                parameterStatsByPanel(parameterIdx, viewIndices), ...
                referenceValues(parameterIdx, viewIndices));
            rowYLim = annotationRowYLim( ...
                sharedYLimits(parameterIdx, :), rowBounds, 3, opts);

            for viewIdx = viewIndices
                ax = parameterAxes(parameterIdx, viewIdx);
                ylim(ax, rowYLim);
                tickInfo = applyDistributionYTicks( ...
                    ax, 'parameter');
                annotationOpts = sharedAnnotationOptions( ...
                    opts, rowBounds, rowYLim, 'parameter');
                parameterAnnotation = addSignificanceBrackets( ...
                    ax, conditionPairs, ...
                    parameterTestsByPanel{parameterIdx, viewIdx}, ...
                    annotationOpts);
                viewName = viewNames{viewIdx};
                distributionData(aggregateIdx).(viewName). ...
                    parameterAnnotations{parameterIdx} = parameterAnnotation;
                parameterAnnotations{parameterIdx, viewIdx} = ...
                    parameterAnnotation;
                printPanelValidation(clusterTitle, ...
                    [viewTitles{viewIdx} ' ' parameterLabels{parameterIdx}], ...
                    conditionLabels, parameterCounts{parameterIdx, viewIdx}, ...
                    conditionPairs, ...
                    parameterTestsByPanel{parameterIdx, viewIdx}, ...
                    parameterAnnotations{parameterIdx, viewIdx}, ...
                    tickInfo, referenceValues(parameterIdx, viewIdx), ax);
                validRows = clusterAudits{viewIdx}.validMdlRows(:);
                validBlocks = ...
                    clusterAudits{viewIdx}.validSourceBlocks(:);
                validExperimentIDs = lookupExperimentIDs( ...
                    validBlocks, opts);
                parameterName = parameterAuditRowLabel(parameterIdx);
                familyID = sprintf('%s|C%d|%s|parameter|%s', ...
                    datasetLabel(aggregateFits(aggregateIdx)), ...
                    aggregateFits(aggregateIdx).clusterID, ...
                    viewName, parameterName);
                panelStatsAudit = buildDistributionStatsAuditRows( ...
                    aggregateFits(aggregateIdx), ...
                    'parameterDistribution', parameterName, ...
                    viewName, parameterName, familyID, ...
                    parameterTestsByPanel{parameterIdx, viewIdx}, ...
                    parameterAnnotations{parameterIdx, viewIdx}, ...
                    validRows, validBlocks, validExperimentIDs);
                statsAudit = [statsAudit; panelStatsAudit]; %#ok<AGROW>
                printDistributionStatsAudit(panelStatsAudit, opts.alpha);
            end
        end
        addParameterTitle(clusterTitle);
    end
end

function opts = applyDefaults(opts)
    defaults = struct( ...
        'figureVisible', 'on', ...
        'jitterWidth', 0.13, ...
        'pointSize', 55, ...
        'meanPointSize', 125, ...
        'pointAlpha', 0.68, ...
        'connectPairedValues', true, ...
        'showIQR', false, ...
        'showNS', true, ...
        'alpha', 0.05, ...
        'annotationLineSpacingFraction', 0.16, ...
        'annotationStarLabelGapFraction', 0.004, ...
        'annotationNSLabelGapFraction', 0.010, ...
        'annotationTopPadFraction', 0.08, ...
        'annotationFinalTopPadFraction', 0.05, ...
        'showSourceAuditText', false, ...
        'experimentIDsByBlock', [], ...
        'baselineModeByBlock', [], ...
        'modelFieldName', 'input mdl', ...
        'statsOpts', struct('testDirection', 'two-sided'));
    names = fieldnames(defaults);
    for idx = 1:numel(names)
        if ~isfield(opts, names{idx}) || isempty(opts.(names{idx}))
            opts.(names{idx}) = defaults.(names{idx});
        end
    end
    opts.figureVisible = validatestring(opts.figureVisible, {'on', 'off'});
    if ~isstruct(opts.statsOpts)
        error('plotPowerClusterFitParameters:InvalidStatsOpts', ...
            'opts.statsOpts must be a struct.');
    end
    if ~isfield(opts.statsOpts, 'testDirection') || ...
            isempty(opts.statsOpts.testDirection)
        opts.statsOpts.testDirection = 'two-sided';
    end
    opts.statsOpts.testDirection = validatestring( ...
        opts.statsOpts.testDirection, {'two-sided', 'one-sided'});
end

function [figureHandles, distributionData, sourceAudit, statsAudit] = plotSignedBX0ParameterDistribution(mdl, aggregateFits, opts)
    figureHandles = gobjects(numel(aggregateFits), 1);
    distributionData = repmat(struct(), numel(aggregateFits), 1);
    sourceAudit = emptyDistributionSourceAuditTable();
    statsAudit = emptyDistributionStatsAuditTable();

    fieldNames = {'deltaB', 'deltaX0', 'BHorizontal', 'BVertical', ...
        'X0Horizontal', 'X0Vertical', 'ACon', 'AIncon', ...
        'alphaCon', 'alphaIncon', 'betaCon', 'betaIncon'};
    fieldLabels = {'\DeltaB', '\DeltaX0', 'B_H', 'B_V', ...
        'X0_H', 'X0_V', 'A_{con}', 'A_{incon}', ...
        '\alpha_{con}', '\alpha_{incon}', '\beta_{con}', '\beta_{incon}'};
    fieldColors = [ ...
        0.35 0.10 0.65; 0.35 0.10 0.65; ...
        0.55 0 0; 0 0.05 0.45; ...
        0.55 0 0; 0 0.05 0.45; ...
        0.9294 0.1098 0.1373; 0 0.0941 0.6627; ...
        0.9294 0.1098 0.1373; 0 0.0941 0.6627; ...
        0.9294 0.1098 0.1373; 0 0.0941 0.6627];
    referenceValues = [0 0 50 50 0 0 nan nan nan nan nan nan];

    fprintf('\n=== Signed-BX0 parameter distribution ===\n');
    fprintf('  source field: mdl.signedBX0.* individual-session fitted parameters\n');
    fprintf('  rows available: %d\n', numel(mdl.signedBX0.deltaB));

    for aggregateIdx = 1:numel(aggregateFits)
        rows = unique(aggregateFits(aggregateIdx).mdlRowIndices(:)', 'stable');
        clusterTitle = parameterTitle(aggregateFits(aggregateIdx), numel(rows));
        figureHandles(aggregateIdx) = figure(...
            'Name', [clusterTitle ' signed-BX0 parameter distributions'], ...
            'Color', 'w', ...
            'Visible', opts.figureVisible);
        makeSubplot = @(position) subtightplot(4, 3, position, ...
            [0.08 0.055], [0.10 0.08], [0.09 0.04]);

        distributionData(aggregateIdx).clusterID = aggregateFits(aggregateIdx).clusterID;
        distributionData(aggregateIdx).mdlRowIndices = rows;
        for fieldIdx = 1:numel(fieldNames)
            ax = makeSubplot(fieldIdx);
            hold(ax, 'on');
            values = getSignedBX0FieldValues(mdl, fieldNames{fieldIdx}, rows);
            validValues = values(isfinite(values));
            if isfinite(referenceValues(fieldIdx))
                yline(ax, referenceValues(fieldIdx), '--', ...
                    'Color', 0.4 .* [1 1 1], 'LineWidth', 1.2, ...
                    'HandleVisibility', 'off');
            end
            if ~isempty(validValues)
                jitter = opts.jitterWidth .* (rand(size(validValues)) - 0.5);
                scatter(ax, 1 + jitter, validValues, opts.pointSize, ...
                    fieldColors(fieldIdx, :), 'filled', ...
                    'MarkerFaceAlpha', opts.pointAlpha, ...
                    'MarkerEdgeColor', 'k');
                mu = mean(validValues, 'omitnan');
                sem = std(validValues, 'omitnan') ./ sqrt(sum(isfinite(validValues)));
                errorbar(ax, 1, mu, sem, 'ko', ...
                    'MarkerFaceColor', 'w', 'MarkerSize', 8, ...
                    'LineWidth', 1.5, 'HandleVisibility', 'off');
            else
                text(ax, 0.5, 0.5, 'No valid values', 'Units', 'normalized', ...
                    'HorizontalAlignment', 'center', 'Color', [0.35 0.35 0.35]);
            end
            xlim(ax, [0.5 1.5]);
            set(ax, 'XTick', 1, 'XTickLabel', {fieldLabels{fieldIdx}});
            ylabel(ax, fieldLabels{fieldIdx}, 'Interpreter', 'tex');
            box(ax, 'off');
            axis(ax, 'square');
            distributionData(aggregateIdx).(fieldNames{fieldIdx}) = values;
            fprintf('  %s %s: n=%d mean=%0.4g SEM=%0.4g\n', ...
                clusterTitle, fieldNames{fieldIdx}, numel(validValues), ...
                mean(validValues, 'omitnan'), ...
                std(validValues, 'omitnan') ./ sqrt(max(1, numel(validValues))));
        end
        addParameterTitle([clusterTitle ' | signed-BX0 individual parameters']);
        upFontSize(16, 0.01);
    end
end

function values = getSignedBX0FieldValues(mdl, fieldName, rows)
    values = nan(size(rows));
    if ~isfield(mdl.signedBX0, fieldName)
        warning('plotPowerClusterFitParameters:MissingSignedBX0Field', ...
            'mdl.signedBX0.%s is missing; leaving panel empty.', fieldName);
        return;
    end
    sourceValues = mdl.signedBX0.(fieldName)(:);
    validRows = rows(isfinite(rows) & rows == round(rows) & ...
        rows >= 1 & rows <= numel(sourceValues));
    [~, loc] = ismember(validRows, rows);
    values(loc) = sourceValues(validRows);
end

function [figureHandles, distributionData, sourceAudit, statsAudit] = plotSignedX0HorizontalDistribution(mdl, aggregateFits, opts)
    figureHandles = gobjects(numel(aggregateFits), 1);
    distributionData = repmat(struct(), numel(aggregateFits), 1);
    sourceAudit = emptyDistributionSourceAuditTable();
    statsAudit = emptyDistributionStatsAuditTable();

    x0 = mdl.signedX0.X0Horizontal(:);
    fprintf('\n=== Signed-X0 parameter distribution ===\n');
    fprintf('  source field: mdl.signedX0.X0Horizontal; X0Vertical = -X0Horizontal\n');
    fprintf('  rows available: %d\n', numel(x0));

    for aggregateIdx = 1:numel(aggregateFits)
        rows = unique(aggregateFits(aggregateIdx).mdlRowIndices(:)', 'stable');
        validRows = rows(isfinite(rows) & rows == round(rows) & ...
            rows >= 1 & rows <= numel(x0) & isfinite(x0(rows)));
        values = x0(validRows);
        clusterTitle = parameterTitle(aggregateFits(aggregateIdx), numel(rows));

        figureHandles(aggregateIdx) = figure(...
            'Name', [clusterTitle ' X0 horizontal-opto distribution'], ...
            'Color', 'w', ...
            'Visible', opts.figureVisible);
        ax = axes('Parent', figureHandles(aggregateIdx));
        hold(ax, 'on');
        yline(ax, 0, '--', 'Color', 0.4 .* [1 1 1], ...
            'LineWidth', 1.5, 'HandleVisibility', 'off');
        if ~isempty(values)
            jitter = opts.jitterWidth .* (rand(size(values)) - 0.5);
            scatter(ax, 1 + jitter, values, opts.pointSize, ...
                [0.35 0.10 0.65], 'filled', ...
                'MarkerFaceAlpha', opts.pointAlpha, ...
                'MarkerEdgeColor', 'k');
            mu = mean(values, 'omitnan');
            sem = std(values, 'omitnan') ./ sqrt(sum(isfinite(values)));
            errorbar(ax, 1, mu, sem, 'ko', ...
                'MarkerFaceColor', 'w', ...
                'MarkerSize', 8, ...
                'LineWidth', 1.5, ...
                'HandleVisibility', 'off');
        end
        xlim(ax, [0.5 1.5]);
        set(ax, 'XTick', 1, 'XTickLabel', {'H-Opto'});
        ylabel(ax, 'X0 horizontal-opto (%)');
        title(ax, sprintf('%s | n=%d', clusterTitle, numel(validRows)), ...
            'Interpreter', 'none', 'FontWeight', 'normal');
        axis(ax, 'square');
        upFontSize(18, 0.01);

        distributionData(aggregateIdx).clusterID = aggregateFits(aggregateIdx).clusterID;
        distributionData(aggregateIdx).mdlRowIndices = rows;
        distributionData(aggregateIdx).validMdlRows = validRows;
        distributionData(aggregateIdx).X0Horizontal = values;
        fprintf('  %s: valid X0Horizontal dots %d of %d | mean %.3g | SEM %.3g\n', ...
            clusterTitle, numel(validRows), numel(rows), ...
            mean(values, 'omitnan'), std(values, 'omitnan') ./ sqrt(max(1, sum(isfinite(values)))));
    end
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

function source = auditFitSource(mdl, fitField, headerField, modelFieldName)
    source = struct( ...
        'available', false, ...
        'modelFieldName', modelFieldName, ...
        'fitField', fitField, ...
        'headerField', headerField, ...
        'rawSize', [], ...
        'rawFits', [], ...
        'headers', {{}}, ...
        'headerMap', [], ...
        'failureReason', '');

    if ~isfield(mdl, fitField)
        source.failureReason = sprintf('mdl.%s is missing.', fitField);
        return;
    end
    if ~isfield(mdl, headerField)
        source.rawSize = size(mdl.(fitField));
        source.failureReason = sprintf( ...
            'mdl.%s exists, but mdl.%s is missing.', fitField, headerField);
        return;
    end

    raw = mdl.(fitField);
    source.rawSize = size(raw);
    if size(raw, 2) < 13
        source.failureReason = sprintf( ...
            'mdl.%s has size %s; expected 12 parameters plus AICc.', ...
            fitField, mat2str(source.rawSize));
        return;
    end
    if ndims(raw) >= 3
        raw = raw(:, :, 1);
    end

    source.rawFits = raw;
    source.headers = mdl.(headerField);
    [source.headerMap, mapFailure] = ...
        buildHeaderMap(source.headers, size(raw, 2));
    if ~isempty(mapFailure)
        source.failureReason = mapFailure;
        return;
    end
    source.available = true;
end

function printSourceAudit(source, viewTitle)
    fprintf('%s:\n', viewTitle);
    fprintf('  field: mdl.%s\n', source.fitField);
    if isempty(source.rawSize)
        fprintf('  size: unavailable\n');
    else
        fprintf('  size: %s\n', mat2str(source.rawSize));
    end
    if source.available
        requiredColumns = source.headerMap(:)';
        validRows = all(isfinite(source.rawFits(:, requiredColumns)), 2);
        fprintf('  complete individual fit rows: %d of %d\n', ...
            sum(validRows), size(source.rawFits, 1));
        parameterNames = {'A', 'B', 'alpha', 'beta'};
        for parameterIdx = 1:4
            sourceColumns = parameterSourceColumns( ...
                source.headerMap, parameterIdx);
            fprintf('  %s mapping: baseline=%s | con=%s | incon=%s\n', ...
                parameterNames{parameterIdx}, sourceColumns{1}, ...
                sourceColumns{2}, sourceColumns{3});
        end
    else
        warning('plotPowerClusterFitParameters:FitSourceUnavailable', ...
            '%s parameter source unavailable: %s', ...
            viewTitle, source.failureReason);
    end
end

function [headerMap, failureReason] = buildHeaderMap(headers, nColumns)
    parameterNames = {'A', 'B', '\alpha', '\beta'};
    headerMap = nan(3, 4);
    failureReason = '';
    if ~iscell(headers)
        failureReason = sprintf( ...
            'Parameter headers must be a cell array, but found %s.', ...
            class(headers));
        return;
    end

    headers = headers(:)';
    for columnIdx = 1:min(numel(headers), nColumns)
        header = headers{columnIdx};
        if ~(ischar(header) || (isstring(header) && isscalar(header)))
            continue;
        end
        header = char(header);
        if startsWith(header, 'AUC') || contains(header, 'AICc')
            continue;
        end
        [conditionIdx, parameterIdx] = parseHeader(header, parameterNames);
        if isnan(conditionIdx) || isnan(parameterIdx)
            continue;
        end
        if isfinite(headerMap(conditionIdx, parameterIdx))
            failureReason = sprintf( ...
                'Duplicate parameter header mapping for %s.', header);
            return;
        end
        headerMap(conditionIdx, parameterIdx) = columnIdx;
    end

    if any(~isfinite(headerMap(:)))
        failureReason = ...
            'Could not map all 12 delta-coded Weibull entries from the headers.';
    end
end

function [conditionIdx, parameterIdx] = parseHeader(header, parameterNames)
    conditionIdx = nan;
    if contains(header, 'incon-bl')
        conditionIdx = 3;
    elseif contains(header, 'con-bl')
        conditionIdx = 2;
    elseif contains(header, '^{bl}')
        conditionIdx = 1;
    end

    parameterName = regexprep(header, '^\\Delta', '');
    parameterName = regexprep(parameterName, '\^\{[^}]+\}', '');
    parameterName = strtrim(parameterName);
    parameterIdx = find(strcmp(parameterName, parameterNames), 1);
    if isempty(parameterIdx)
        parameterIdx = nan;
    end
end

function [values, audit] = extractClusterParameters( ...
        source, rows, sourceBlocks, clusterTitle, viewTitle)
    audit = struct( ...
        'fitField', source.fitField, ...
        'rawSize', source.rawSize, ...
        'requestedMdlRows', rows, ...
        'requestedSourceBlocks', sourceBlocks, ...
        'validMdlRows', [], ...
        'validSourceBlocks', [], ...
        'excludedMdlRows', [], ...
        'excludedSourceBlocks', [], ...
        'exclusionReasons', {{}}, ...
        'requestedReasons', {repmat({''}, size(rows))}, ...
        'reconstructedValues', nan(numel(rows), 3, 4), ...
        'nRequestedSessions', numel(rows), ...
        'nValidFits', 0);
    values = [];

    if ~source.available
        audit.excludedMdlRows = rows;
        audit.excludedSourceBlocks = sourceBlocks;
        audit.exclusionReasons = repmat( ...
            {source.failureReason}, size(rows));
        audit.requestedReasons = audit.exclusionReasons;
        return;
    end

    validMask = isfinite(rows) & rows == round(rows) & ...
        rows >= 1 & rows <= size(source.rawFits, 1);
    reasons = repmat({''}, size(rows));
    reasons(~validMask) = {sprintf( ...
        'mdl row is outside %s rows 1:%d', ...
        source.fitField, size(source.rawFits, 1))};

    requiredColumns = source.headerMap(:)';
    for idx = find(validMask)
        rawRow = source.rawFits(rows(idx), :);
        badColumns = requiredColumns(~isfinite(rawRow(requiredColumns)));
        if ~isempty(badColumns)
            validMask(idx) = false;
            reasons{idx} = sprintf( ...
                'nonfinite fitted columns %s', mat2str(badColumns));
        end
    end

    audit.validMdlRows = rows(validMask);
    audit.validSourceBlocks = sourceBlocks(validMask);
    audit.excludedMdlRows = rows(~validMask);
    audit.excludedSourceBlocks = sourceBlocks(~validMask);
    audit.exclusionReasons = reasons(~validMask);
    audit.requestedReasons = reasons;
    audit.nValidFits = sum(validMask);

    if audit.nValidFits == 0
        warning('plotPowerClusterFitParameters:NoValidFits', ...
            '%s, %s: no complete individual fits are available.', ...
            clusterTitle, viewTitle);
        return;
    end

    values = reconstructActualParameters( ...
        source.rawFits(audit.validMdlRows, :), source.headerMap);
    audit.reconstructedValues(validMask, :, :) = values;
end

function values = reconstructActualParameters(rawFits, headerMap)
    nRows = size(rawFits, 1);
    encoded = nan(nRows, 3, 4);
    for conditionIdx = 1:3
        for parameterIdx = 1:4
            encoded(:, conditionIdx, parameterIdx) = ...
                rawFits(:, headerMap(conditionIdx, parameterIdx));
        end
    end

    baseline = reshape(encoded(:, 1, :), nRows, 4);
    conDelta = reshape(encoded(:, 2, :), nRows, 4);
    inconDelta = reshape(encoded(:, 3, :), nRows, 4);

    actualBaseline = [baseline(:, 1), 0.5 .* ones(nRows, 1), baseline(:, 3:4)];
    actualCon = [baseline(:, 1) + conDelta(:, 1), ...
        0.5 + conDelta(:, 2), ...
        baseline(:, 3) + conDelta(:, 3), ...
        baseline(:, 4) + conDelta(:, 4)];
    actualIncon = [baseline(:, 1) + inconDelta(:, 1), ...
        0.5 - conDelta(:, 2), ...
        baseline(:, 3) + inconDelta(:, 3), ...
        baseline(:, 4) + inconDelta(:, 4)];

    values = nan(nRows, 3, 4);
    values(:, 1, :) = reshape(actualBaseline, nRows, 1, 4);
    values(:, 2, :) = reshape(actualCon, nRows, 1, 4);
    values(:, 3, :) = reshape(actualIncon, nRows, 1, 4);
end

function limits = computeSharedYLimits(parameterValues, parameterLabels, clusterTitle)
    defaults = [25, 100, 50, 8];
    limits = nan(4, 2);
    for parameterIdx = 1:4
        allValues = [];
        for viewIdx = 1:numel(parameterValues)
            if isempty(parameterValues{viewIdx})
                continue;
            end
            values = parameterValues{viewIdx}(:, :, parameterIdx);
            if parameterIdx <= 2
                values = 100 .* values;
            end
            allValues = [allValues; values(:)]; %#ok<AGROW>
        end
        allValues = allValues(isfinite(allValues));
        if isempty(allValues)
            limits(parameterIdx, :) = [0 defaults(parameterIdx)];
            continue;
        end

        maxValue = max(allValues);
        yMax = max(1.10 .* maxValue, 0.10 .* defaults(parameterIdx));
        if min(allValues) < 0
            warning('plotPowerClusterFitParameters:NegativeActualParameter', ...
                ['%s, %s: reconstructed actual values include negatives ' ...
                '(minimum %.4g); the requested display axis starts at 0.'], ...
                clusterTitle, parameterLabels{parameterIdx}, min(allValues));
        end
        limits(parameterIdx, :) = [0 yMax];
    end
end

function printExclusions(audit, clusterTitle, viewTitle)
    if isempty(audit.excludedMdlRows)
        fprintf('  %s excluded fits: none\n', viewTitle);
        return;
    end
    warning('plotPowerClusterFitParameters:ExcludedFits', ...
        '%s, %s: excluded %d of %d sessions.', ...
        clusterTitle, viewTitle, numel(audit.excludedMdlRows), ...
        audit.nRequestedSessions);
    for idx = 1:numel(audit.excludedMdlRows)
        fprintf('    mdl row %g, source block %g: %s\n', ...
            audit.excludedMdlRows(idx), ...
            audit.excludedSourceBlocks(idx), ...
            audit.exclusionReasons{idx});
    end
end

function sourceBlocks = getSourceBlocks(item, rows)
    if isfield(item, 'sourceBlockIndices') && ...
            numel(item.sourceBlockIndices) == numel(rows)
        sourceBlocks = item.sourceBlockIndices(:)';
    else
        sourceBlocks = nan(size(rows));
    end
end

function experimentIDs = lookupExperimentIDs(sourceBlocks, opts)
    experimentIDs = strings(numel(sourceBlocks), 1);
    if isempty(opts.experimentIDsByBlock)
        return;
    end
    idsByBlock = opts.experimentIDsByBlock;
    for rowIdx = 1:numel(sourceBlocks)
        blockIdx = sourceBlocks(rowIdx);
        if isfinite(blockIdx) && blockIdx == round(blockIdx) && ...
                blockIdx >= 1 && blockIdx <= numel(idsByBlock)
            experimentIDs(rowIdx) = string(idsByBlock(blockIdx));
        end
    end
end

function baselineModes = lookupBaselineModes(sourceBlocks, opts)
    baselineModes = strings(numel(sourceBlocks), 1);
    if ~isfield(opts, 'baselineModeByBlock') || isempty(opts.baselineModeByBlock)
        return;
    end
    modesByBlock = string(opts.baselineModeByBlock(:));
    for rowIdx = 1:numel(sourceBlocks)
        blockIdx = sourceBlocks(rowIdx);
        if isfinite(blockIdx) && blockIdx == round(blockIdx) && ...
                blockIdx >= 1 && blockIdx <= numel(modesByBlock)
            baselineModes(rowIdx) = modesByBlock(blockIdx);
        end
    end
end

function columns = parameterSourceColumns(headerMap, parameterIdx)
    baselineColumn = headerMap(1, parameterIdx);
    conColumn = headerMap(2, parameterIdx);
    inconColumn = headerMap(3, parameterIdx);
    if parameterIdx == 2
        columns = { ...
            sprintf('fixed 0.5; displayed as 100*0.5 (header col %d unused)', ...
            baselineColumn), ...
            sprintf('100*(0.5 + col%d)', conColumn), ...
            sprintf('100*(0.5 - col%d)', conColumn)};
        return;
    end

    prefix = '';
    if parameterIdx == 1
        prefix = '100*';
    end
    columns = { ...
        sprintf('%scol%d', prefix, baselineColumn), ...
        sprintf('%s(col%d + col%d)', prefix, ...
        baselineColumn, conColumn), ...
        sprintf('%s(col%d + col%d)', prefix, ...
        baselineColumn, inconColumn)};
end

function label = parameterAuditRowLabel(parameterIdx)
    labels = {'A', 'B', 'alpha', 'beta'};
    label = labels{parameterIdx};
end

function assertParameterAuditCounts(audit, expectedCounts, clusterTitle, ...
        viewTitle, parameterLabel, nValidFits)
    metrics = ["baseline", "con", "incon"];
    for metricIdx = 1:numel(metrics)
        metricRows = audit.conditionOrMetric == metrics(metricIdx);
        validRows = metricRows & audit.isValidForPlot;
        nAudit = sum(validRows);
        if nAudit ~= expectedCounts(metricIdx)
            error('plotPowerClusterFitParameters:AuditCountMismatch', ...
                ['%s, %s, %s, %s: audit has %d valid points, but ' ...
                'plotting summary expects %d.'], ...
                clusterTitle, viewTitle, parameterLabel, ...
                metrics(metricIdx), nAudit, expectedCounts(metricIdx));
        end
        if nAudit == 1 && nValidFits > 1
            warning('plotPowerClusterFitParameters:SuspiciousSingleDot', ...
                ['%s, %s, %s, %s has one plotted dot despite %d valid ' ...
                'individual-session fits.'], ...
                clusterTitle, viewTitle, parameterLabel, ...
                metrics(metricIdx), nValidFits);
        end
        excluded = audit(metricRows & ~audit.isValidForPlot, :);
        for excludedIdx = 1:height(excluded)
            fprintf(['  Excluded parameter source: %s | %s | %s | ' ...
                'sessionRow=%g | block=%g | reason=%s\n'], ...
                viewTitle, parameterLabel, metrics(metricIdx), ...
                excluded.sessionRowIndex(excludedIdx), ...
                excluded.blockIndex(excludedIdx), ...
                excluded.reasonExcluded(excludedIdx));
        end
    end
end

function printSourceAuditSummary(audit)
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

function verifyPointCounts(counts, expectedCount, clusterTitle, ...
        viewTitle, parameterLabel, conditionLabels)
    for conditionIdx = 1:numel(counts)
        if counts(conditionIdx) ~= expectedCount
            warning('plotPowerClusterFitParameters:DotCountMismatch', ...
                ['%s, %s, %s, %s: plotted %d dots, but %d valid ' ...
                'individual fits were expected.'], ...
                clusterTitle, viewTitle, parameterLabel, ...
                conditionLabels{conditionIdx}, counts(conditionIdx), expectedCount);
        end
    end
end

function printParameterMeans(viewTitle, parameterLabel, stats)
    fprintf(['  %s %s mean +/- SEM: baseline=%0.4f +/- %0.4f, ' ...
        'con=%0.4f +/- %0.4f, incon=%0.4f +/- %0.4f\n'], ...
        viewTitle, parameterLabel, ...
        stats.mean(1), stats.sem(1), ...
        stats.mean(2), stats.sem(2), ...
        stats.mean(3), stats.sem(3));
end

function styleParameterAxes(ax, labels)
    xlim(ax, [0.5 3.5]);
    xticks(ax, 1:3);
    xticklabels(ax, labels);
    xtickangle(ax, 15);
    box(ax, 'off');
    set(ax, 'LineWidth', 2, 'TickDir', 'out', ...
        'TickLength', [0.01 0.01], 'FontName', 'FreeSans', ...
        'FontSize', 15);
end

function plotOpts = conditionMarkerOptions(opts, colors)
    plotOpts = opts;
    plotOpts.groupMarkers = {'o', '^', 'v'};
    plotOpts.groupFaceColors = [1 1 1; colors(2, :); colors(3, :)];
    plotOpts.groupEdgeColors = zeros(3, 3);
    plotOpts.semColor = [0 0 0];
end

function bounds = sharedParameterAnnotationBounds( ...
        valuesByPanel, statsByPanel, referenceValues)
    allValues = referenceValues(:);
    allUpper = referenceValues(:);
    allLower = referenceValues(:);
    for panelIdx = 1:numel(valuesByPanel)
        values = valuesByPanel{panelIdx};
        stats = statsByPanel{panelIdx};
        allValues = [allValues; values(isfinite(values))]; %#ok<AGROW>
        upperSummary = stats.mean + stats.sem;
        lowerSummary = stats.mean - stats.sem;
        allUpper = [allUpper; upperSummary(isfinite(upperSummary))']; %#ok<AGROW>
        allLower = [allLower; lowerSummary(isfinite(lowerSummary))']; %#ok<AGROW>
    end
    bounds = [min([allValues; allLower]), max([allValues; allUpper])];
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

function referenceValue = parameterReferenceValue(parameterIdx, values)
    if parameterIdx == 2
        finiteValues = values(isfinite(values));
        if isempty(finiteValues) || median(abs(finiteValues), 'omitnan') < 2
            referenceValue = 0.5;
        else
            referenceValue = 50;
        end
    else
        referenceValue = 0;
    end
end

function printPanelValidation(clusterTitle, panelTitle, labels, counts, ...
        pairs, results, annotationInfo, tickInfo, referenceValue, ax)
    fprintf(['  Validation: %s | %s | n/group=%s | ylim=%s | ' ...
        'reference=%g | annotations=%d\n'], ...
        clusterTitle, panelTitle, mat2str(counts), ...
        mat2str(tickInfo.limits, 4), referenceValue, annotationInfo.nDrawn);
    for resultIdx = 1:numel(results)
        if annotationInfo.drawnMask(resultIdx)
            annotationStatus = 'drawn';
        else
            annotationStatus = ['skipped: ' ...
                annotationInfo.skippedReasons{resultIdx}];
        end
        fprintf(['    %s vs %s: %s, n=%d, raw p=%0.4g, ' ...
            'Holm p=%0.4g, %s, annotation %s, lineY=%0.4g, ' ...
            'textY=%0.4g, lineSpacing=%0.4g, labelGap=%0.4g, ' ...
            'fontSize=%g\n'], ...
            labels{pairs(resultIdx, 1)}, labels{pairs(resultIdx, 2)}, ...
            results(resultIdx).test, results(resultIdx).n, ...
            results(resultIdx).rawP, results(resultIdx).adjustedP, ...
            results(resultIdx).star, annotationStatus, ...
            annotationInfo.pairs(resultIdx).yLine, ...
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

function markUnavailablePanel(ax, parameterLabel)
    axis(ax, 'off');
    text(ax, 0.5, 0.5, 'Individual fits unavailable', ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'Color', [0.35 0.35 0.35], ...
        'FontSize', 11);
    text(ax, 0.02, 0.95, parameterLabel, ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', ...
        'FontSize', 11);
end

function titleText = parameterTitle(item, nExperiments)
    titleText = sprintf('Power cluster %d individual Weibull parameters', ...
        item.clusterID);
    if isfield(item, 'columnTargetLabel')
        titleText = sprintf('%s, Columns %s', ...
            titleText, item.columnTargetLabel);
    end
    titleText = sprintf('%s (n_{expt}=%d)', titleText, nExperiments);
end

function addParameterTitle(titleText)
    annotation(gcf, 'textbox', [0.10 0.925 0.80 0.055], ...
        'String', titleText, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', 16, ...
        'Interpreter', 'tex', ...
        'EdgeColor', 'none');
end
