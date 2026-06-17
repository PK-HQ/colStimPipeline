function audit = buildDistributionSourceAuditRows(item, figureType, ...
        panelColumn, panelRow, groupLabels, sessionRows, blockIndices, ...
        experimentIDs, sourceField, sourceColumns, values, validMask, reasons, ...
        baselineModes)
% Build one audit row per candidate individual-session plotted value.

    sessionRows = sessionRows(:);
    blockIndices = normalizeNumericColumn(blockIndices, numel(sessionRows));
    experimentIDs = normalizeStringColumn(experimentIDs, numel(sessionRows));
    if nargin < 14 || isempty(baselineModes)
        baselineModes = strings(numel(sessionRows), 1);
    else
        baselineModes = normalizeStringColumn(baselineModes, numel(sessionRows));
    end
    groupLabels = string(groupLabels(:)');
    sourceColumns = string(sourceColumns(:)');

    nSessions = numel(sessionRows);
    nGroups = numel(groupLabels);
    if ~isequal(size(values), [nSessions nGroups])
        error('buildDistributionSourceAuditRows:ValueSizeMismatch', ...
            'values must have size [%d %d], but found %s.', ...
            nSessions, nGroups, mat2str(size(values)));
    end
    if nargin < 12 || isempty(validMask)
        validMask = isfinite(values);
    end
    if ~isequal(size(validMask), size(values))
        error('buildDistributionSourceAuditRows:ValiditySizeMismatch', ...
            'validMask must match values.');
    end
    if nargin < 13 || isempty(reasons)
        reasons = strings(size(values));
    else
        reasons = string(reasons);
    end
    if ~isequal(size(reasons), size(values))
        error('buildDistributionSourceAuditRows:ReasonSizeMismatch', ...
            'reasons must match values.');
    end

    datasetLabel = normalizeDatasetLabel(item);
    clusterID = item.clusterID;
    nRows = nSessions * nGroups;
    if nRows == 0
        audit = emptyDistributionSourceAuditTable();
        return;
    end

    conditionOrMetric = strings(nRows, 1);
    sessionRowIndex = zeros(nRows, 1);
    blockIndex = zeros(nRows, 1);
    experimentID = strings(nRows, 1);
    baselineMode = strings(nRows, 1);
    sourceSubscriptOrColumn = strings(nRows, 1);
    plottedValue = zeros(nRows, 1);
    isValidForPlot = false(nRows, 1);
    reasonExcluded = strings(nRows, 1);

    outputRow = 0;
    for sessionIdx = 1:nSessions
        for groupIdx = 1:nGroups
            outputRow = outputRow + 1;
            conditionOrMetric(outputRow) = groupLabels(groupIdx);
            sessionRowIndex(outputRow) = sessionRows(sessionIdx);
            blockIndex(outputRow) = blockIndices(sessionIdx);
            experimentID(outputRow) = experimentIDs(sessionIdx);
            baselineMode(outputRow) = baselineModes(sessionIdx);
            sourceSubscriptOrColumn(outputRow) = ...
                sourceColumns(min(groupIdx, numel(sourceColumns)));
            plottedValue(outputRow) = values(sessionIdx, groupIdx);
            isValidForPlot(outputRow) = validMask(sessionIdx, groupIdx);
            if isValidForPlot(outputRow)
                reasonExcluded(outputRow) = "";
            elseif strlength(reasons(sessionIdx, groupIdx)) > 0
                reasonExcluded(outputRow) = reasons(sessionIdx, groupIdx);
            else
                reasonExcluded(outputRow) = "nonfinite source value";
            end
        end
    end

    schema = emptyDistributionSourceAuditTable();
    audit = table( ...
        repmat(datasetLabel, nRows, 1), ...
        repmat(clusterID, nRows, 1), ...
        repmat(string(figureType), nRows, 1), ...
        repmat(string(panelColumn), nRows, 1), ...
        repmat(string(panelRow), nRows, 1), ...
        conditionOrMetric, ...
        sessionRowIndex, ...
        blockIndex, ...
        experimentID, ...
        baselineMode, ...
        repmat(string(sourceField), nRows, 1), ...
        sourceSubscriptOrColumn, ...
        plottedValue, ...
        isValidForPlot, ...
        reasonExcluded, ...
        repmat(clusterID, nRows, 1), ...
        'VariableNames', schema.Properties.VariableNames);
end

function label = normalizeDatasetLabel(item)
    label = "unknown";
    if isfield(item, 'columnTargetGroup')
        if strcmpi(item.columnTargetGroup, 'experimental')
            label = "0/90";
            return;
        elseif strcmpi(item.columnTargetGroup, 'control')
            label = "45/135";
            return;
        end
    end
    if isfield(item, 'columnTargetLabel')
        raw = char(string(item.columnTargetLabel));
        if contains(raw, '45') && contains(raw, '135')
            label = "45/135";
        elseif contains(raw, '0') && contains(raw, '90')
            label = "0/90";
        else
            label = string(raw);
        end
    end
end

function values = normalizeNumericColumn(values, nRows)
    if isempty(values)
        values = nan(nRows, 1);
        return;
    end
    values = values(:);
    if numel(values) ~= nRows
        error('buildDistributionSourceAuditRows:MetadataSizeMismatch', ...
            'Numeric metadata length does not match session rows.');
    end
end

function values = normalizeStringColumn(values, nRows)
    if isempty(values)
        values = strings(nRows, 1);
        return;
    end
    values = string(values(:));
    if numel(values) ~= nRows
        error('buildDistributionSourceAuditRows:MetadataSizeMismatch', ...
            'String metadata length does not match session rows.');
    end
end
