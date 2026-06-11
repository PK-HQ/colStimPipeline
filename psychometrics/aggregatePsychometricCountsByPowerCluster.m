function agg = aggregatePsychometricCountsByPowerCluster(mdl, powerCluster, opts)
% Aggregate psychometric binomial counts across sessions by power cluster.
%
% Each output element corresponds to one unique power-cluster label and
% contains horizontal, vertical, and merged baseline/con/incon summaries.

    if nargin < 3 || isempty(opts)
        opts = struct();
    end
    opts = applyOptionDefaults(opts);

    powerCluster = powerCluster(:);
    nBlocks = numel(powerCluster);
    if nBlocks == 0
        agg = struct([]);
        return;
    end
    validateClusterLabels(powerCluster);

    validClusterLabel = getValidClusterLabels(powerCluster);
    if isnumeric(powerCluster) || islogical(powerCluster)
        clusterIDs = unique(powerCluster(validClusterLabel), 'sorted');
    else
        clusterIDs = unique(powerCluster(validClusterLabel), 'stable');
    end
    viewNames = {'horizontal', 'vertical', 'merged'};
    conditionNames = {'baseline', 'con', 'incon'};

    pointStore = initializePointStore(numel(clusterIDs), viewNames, conditionNames);

    for block = 1:nBlocks
        if ~validClusterLabel(block)
            continue;
        end
        blockClusterID = getClusterLabel(powerCluster, block);
        clusterIdx = findClusterIndex(clusterIDs, blockClusterID);
        blockData = extractBlockData(mdl, block, nBlocks);

        for viewIdx = 1:numel(viewNames)
            viewName = viewNames{viewIdx};
            for conditionIdx = 1:numel(conditionNames)
                conditionName = conditionNames{conditionIdx};
                conditionData = blockData.(viewName).(conditionName);
                pointStore(clusterIdx).(viewName).(conditionName) = appendPoints( ...
                    pointStore(clusterIdx).(viewName).(conditionName), ...
                    conditionData.x, conditionData.y, block);
            end
        end
    end

    agg = repmat(struct(), numel(clusterIDs), 1);
    for clusterIdx = 1:numel(clusterIDs)
        clusterID = getClusterLabel(clusterIDs, clusterIdx);
        agg(clusterIdx).clusterID = clusterID;
        agg(clusterIdx).sessionIDs = findClusterSessions(powerCluster, clusterID);
        agg(clusterIdx).opts = opts;

        for viewIdx = 1:numel(viewNames)
            viewName = viewNames{viewIdx};
            for conditionIdx = 1:numel(conditionNames)
                conditionName = conditionNames{conditionIdx};
                agg(clusterIdx).(viewName).(conditionName) = aggregateCondition( ...
                    pointStore(clusterIdx).(viewName).(conditionName), opts);
            end
        end
    end
end

function opts = applyOptionDefaults(opts)
    if ~isstruct(opts) || ~isscalar(opts)
        error('opts must be a scalar struct.');
    end

    defaults = struct( ...
        'binWidth', 5, ...
        'useMedianForPlot', false, ...
        'minSessionsPerBin', 1);

    optionNames = fieldnames(defaults);
    for ii = 1:numel(optionNames)
        optionName = optionNames{ii};
        if ~isfield(opts, optionName) || isempty(opts.(optionName))
            opts.(optionName) = defaults.(optionName);
        end
    end

    validateattributes(opts.binWidth, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
    validateattributes(opts.useMedianForPlot, {'logical', 'numeric'}, {'scalar'});
    validateattributes(opts.minSessionsPerBin, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'integer', 'positive'});

    opts.useMedianForPlot = logical(opts.useMedianForPlot);
end

function validateClusterLabels(powerCluster)
    if isnumeric(powerCluster) || islogical(powerCluster)
        if any(isinf(double(powerCluster)))
            error('powerCluster cannot contain Inf labels.');
        end
    elseif isstring(powerCluster)
        if any(ismissing(powerCluster))
            error('powerCluster cannot contain missing labels.');
        end
    elseif iscategorical(powerCluster)
        if any(isundefined(powerCluster))
            error('powerCluster cannot contain undefined labels.');
        end
    elseif iscellstr(powerCluster)
        if any(cellfun(@isempty, powerCluster))
            error('powerCluster cannot contain empty labels.');
        end
    else
        error('powerCluster must be numeric, logical, string, categorical, or a cell array of character vectors.');
    end
end

function valid = getValidClusterLabels(powerCluster)
    if isnumeric(powerCluster) || islogical(powerCluster)
        valid = ~isnan(double(powerCluster));
    elseif isstring(powerCluster)
        valid = ~ismissing(powerCluster);
    elseif iscategorical(powerCluster)
        valid = ~isundefined(powerCluster);
    else
        valid = ~cellfun(@isempty, powerCluster);
    end
end

function pointStore = initializePointStore(nClusters, viewNames, conditionNames)
    emptyPoints = struct('x', [], 'y', [], 'session', []);
    pointStore = repmat(struct(), nClusters, 1);

    for clusterIdx = 1:nClusters
        for viewIdx = 1:numel(viewNames)
            viewName = viewNames{viewIdx};
            for conditionIdx = 1:numel(conditionNames)
                conditionName = conditionNames{conditionIdx};
                pointStore(clusterIdx).(viewName).(conditionName) = emptyPoints;
            end
        end
    end
end

function blockData = extractBlockData(mdl, block, nBlocks)
    [savedHorizontal, haveSavedHorizontal] = extractSavedPanelSide( ...
        mdl, 'xPanel1Horizontal', 'yPanel1Horizontal', block, nBlocks);
    [savedVertical, haveSavedVertical] = extractSavedPanelSide( ...
        mdl, 'xPanel1Vertical', 'yPanel1Vertical', block, nBlocks);

    if haveSavedHorizontal && haveSavedVertical
        blockData.horizontal = savedHorizontal;
        blockData.vertical = savedVertical;
    else
        sideData = extractPreMergedSides(mdl, block);
        if haveSavedHorizontal
            blockData.horizontal = savedHorizontal;
        else
            blockData.horizontal = sideData.horizontal;
        end
        if haveSavedVertical
            blockData.vertical = savedVertical;
        else
            blockData.vertical = sideData.vertical;
        end
    end

    blockData.merged = extractMergedData(mdl, block);
end

function [sideData, success] = extractSavedPanelSide(mdl, xField, yField, block, nBlocks)
    sideData = emptyConditionData();
    success = false;

    if ~isfield(mdl, xField) || ~isfield(mdl, yField)
        return;
    end

    xData = mdl.(xField);
    yData = mdl.(yField);

    for conditionIdx = 1:3
        [x, okX] = extractConditionVector(xData, conditionIdx, block, nBlocks);
        [y, okY] = extractConditionVector(yData, conditionIdx, block, nBlocks);
        if ~okX || ~okY
            return;
        end
        [x, y] = cleanPairedData(x, y);
        conditionName = conditionNameFromIndex(conditionIdx);
        sideData.(conditionName) = struct('x', abs(x), 'y', y);
    end

    success = true;
end

function [values, success] = extractConditionVector(data, conditionIdx, block, nBlocks)
    values = [];
    success = false;

    if iscell(data)
        if isequal(size(data), [nBlocks, 3])
            values = data{block, conditionIdx};
        elseif isequal(size(data), [3, nBlocks])
            values = data{conditionIdx, block};
        elseif numel(data) == nBlocks
            blockValue = data{block};
            [values, success] = extractConditionVector(blockValue, conditionIdx, 1, 1);
            return;
        else
            return;
        end
        success = isnumeric(values);
        return;
    end

    if ~isnumeric(data)
        return;
    end

    dataSize = size(data);
    if ismatrix(data) && nBlocks == 1
        if dataSize(1) == 3
            values = data(conditionIdx, :);
        elseif dataSize(2) == 3
            values = data(:, conditionIdx);
        else
            return;
        end
    elseif ndims(data) == 3 && dataSize(1) == 3 && dataSize(3) >= block
        values = squeeze(data(conditionIdx, :, block));
    elseif ndims(data) == 3 && dataSize(1) >= block && dataSize(2) == 3
        values = squeeze(data(block, conditionIdx, :));
    elseif ndims(data) == 3 && dataSize(1) >= block && dataSize(3) == 3
        values = squeeze(data(block, :, conditionIdx));
    else
        return;
    end

    success = true;
end

function sideData = extractPreMergedSides(mdl, block)
    requiredFields = { ...
        'xBaselinePreMerge', 'yBaselinePreMerge', ...
        'xHorizontalOptoPreMerge', 'yHorizontalOptoPreMerge', ...
        'visualTagHorizontalOptoPreMerge', 'congruencyHorizontalOptoPreMerge', ...
        'xVerticalOptoPreMerge', 'yVerticalOptoPreMerge', ...
        'visualTagVerticalOptoPreMerge', 'congruencyVerticalOptoPreMerge'};
    requireFields(mdl, requiredFields, 'pre-merged side data');

    haveBaselineTags = isfield(mdl, 'visualTagBaselinePreMerge');
    if haveBaselineTags
        [xBaseline, yBaseline, baselineTags] = getBlockTriplet( ...
            mdl.xBaselinePreMerge, mdl.yBaselinePreMerge, ...
            mdl.visualTagBaselinePreMerge, block);
    else
        [xBaseline, yBaseline] = getBlockPair( ...
            mdl.xBaselinePreMerge, mdl.yBaselinePreMerge, block);
        baselineTags = [];
    end

    [xH, yH, tagH, congrH] = getBlockOptoData(mdl, block, 'Horizontal');
    [xV, yV, tagV, congrV] = getBlockOptoData(mdl, block, 'Vertical');

    xOpto = [xH, xV];
    yOpto = [yH, yV];
    tagOpto = [tagH, tagV];
    congrOpto = [congrH, congrV];

    if haveBaselineTags
        idxBaseHorizontal = baselineTags == 0;
        idxBaseVertical = baselineTags == 90;
    else
        idxBaseHorizontal = xBaseline <= 0;
        idxBaseVertical = xBaseline >= 0;
    end

    idxConHorizontal = tagOpto == 0 & congrOpto == 1;
    idxInconHorizontal = tagOpto == 0 & congrOpto == -1;
    idxConVertical = tagOpto == 90 & congrOpto == 1;
    idxInconVertical = tagOpto == 90 & congrOpto == -1;

    sideData.horizontal = makeConditionData( ...
        xBaseline(idxBaseHorizontal), yBaseline(idxBaseHorizontal), ...
        xOpto(idxConHorizontal), yOpto(idxConHorizontal), ...
        xOpto(idxInconHorizontal), yOpto(idxInconHorizontal));

    sideData.vertical = makeConditionData( ...
        xBaseline(idxBaseVertical), yBaseline(idxBaseVertical), ...
        xOpto(idxConVertical), yOpto(idxConVertical), ...
        xOpto(idxInconVertical), yOpto(idxInconVertical));
end

function [x, y, visualTag, congruency] = getBlockOptoData(mdl, block, sideName)
    xField = ['x' sideName 'OptoPreMerge'];
    yField = ['y' sideName 'OptoPreMerge'];
    tagField = ['visualTag' sideName 'OptoPreMerge'];
    congrField = ['congruency' sideName 'OptoPreMerge'];

    x = getBlockVector(mdl.(xField), block);
    y = getBlockVector(mdl.(yField), block);
    visualTag = getBlockVector(mdl.(tagField), block);
    congruency = getBlockVector(mdl.(congrField), block);

    nValues = min([numel(x), numel(y), numel(visualTag), numel(congruency)]);
    x = x(1:nValues);
    y = y(1:nValues);
    visualTag = visualTag(1:nValues);
    congruency = congruency(1:nValues);

    valid = isfinite(x) & isfinite(y) & isfinite(visualTag) & isfinite(congruency);
    x = x(valid);
    y = y(valid);
    visualTag = visualTag(valid);
    congruency = congruency(valid);
end

function mergedData = extractMergedData(mdl, block)
    requireFields(mdl, {'xBlock', 'yBlock'}, 'merged data');
    if size(mdl.xBlock, 1) < 3 || size(mdl.yBlock, 1) < 3 || ...
            size(mdl.xBlock, 3) < block || size(mdl.yBlock, 3) < block
        error('mdl.xBlock and mdl.yBlock must have size 3-by-points-by-blocks.');
    end

    mergedData = emptyConditionData();
    for conditionIdx = 1:3
        x = squeeze(mdl.xBlock(conditionIdx, :, block));
        y = squeeze(mdl.yBlock(conditionIdx, :, block));
        [x, y] = cleanPairedData(x, y);
        conditionName = conditionNameFromIndex(conditionIdx);
        mergedData.(conditionName) = struct('x', abs(x), 'y', y);
    end
end

function conditionData = makeConditionData(xBaseline, yBaseline, xCon, yCon, xIncon, yIncon)
    conditionData = emptyConditionData();
    [xBaseline, yBaseline] = cleanPairedData(xBaseline, yBaseline);
    [xCon, yCon] = cleanPairedData(xCon, yCon);
    [xIncon, yIncon] = cleanPairedData(xIncon, yIncon);

    conditionData.baseline = struct('x', abs(xBaseline), 'y', yBaseline);
    conditionData.con = struct('x', abs(xCon), 'y', yCon);
    conditionData.incon = struct('x', abs(xIncon), 'y', yIncon);
end

function conditionData = emptyConditionData()
    emptyXY = struct('x', [], 'y', []);
    conditionData = struct('baseline', emptyXY, 'con', emptyXY, 'incon', emptyXY);
end

function conditionName = conditionNameFromIndex(conditionIdx)
    conditionNames = {'baseline', 'con', 'incon'};
    conditionName = conditionNames{conditionIdx};
end

function [x, y] = getBlockPair(xData, yData, block)
    x = getBlockVector(xData, block);
    y = getBlockVector(yData, block);
    [x, y] = cleanPairedData(x, y);
end

function [x, y, third] = getBlockTriplet(xData, yData, thirdData, block)
    x = getBlockVector(xData, block);
    y = getBlockVector(yData, block);
    third = getBlockVector(thirdData, block);

    nValues = min([numel(x), numel(y), numel(third)]);
    x = x(1:nValues);
    y = y(1:nValues);
    third = third(1:nValues);

    valid = isfinite(x) & isfinite(y) & isfinite(third);
    x = x(valid);
    y = y(valid);
    third = third(valid);
end

function values = getBlockVector(data, block)
    if iscell(data)
        if numel(data) < block
            error('Cell data do not contain block %d.', block);
        end
        values = data{block};
    elseif isnumeric(data)
        if size(data, 1) < block
            error('Numeric data do not contain block %d.', block);
        end
        values = data(block, :);
    else
        error('Expected numeric or cell block data.');
    end
    values = values(:)';
end

function [x, y] = cleanPairedData(x, y)
    x = x(:)';
    y = y(:)';
    nValues = min(numel(x), numel(y));
    x = x(1:nValues);
    y = y(1:nValues);

    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);
end

function requireFields(mdl, fieldNames, description)
    missing = fieldNames(~isfield(mdl, fieldNames));
    if ~isempty(missing)
        error('Cannot reconstruct %s. Missing mdl field(s): %s.', ...
            description, strjoin(missing, ', '));
    end
end

function points = appendPoints(points, x, y, sessionID)
    [x, y] = cleanPairedData(x, y);
    if isempty(x)
        return;
    end
    if any(y < 0 | y > 100)
        error('Percent-correct values must be between 0 and 100.');
    end

    points.x = [points.x, abs(x)];
    points.y = [points.y, y];
    points.session = [points.session, repmat(sessionID, 1, numel(x))];
end

function summary = aggregateCondition(points, opts)
    summary = emptyAggregateSummary(opts);
    if isempty(points.x)
        return;
    end

    x = abs(points.x);
    y = points.y;
    session = points.session;

    binCenter = opts.binWidth .* round(x ./ opts.binWidth);
    zeroTolerance = sqrt(eps) .* max(1, max(x));
    nTrials = repmat(20, size(x));
    nTrials(x <= zeroTolerance) = 40;
    successes = round((y ./ 100) .* nTrials);
    successes = min(max(successes, 0), nTrials);

    uniqueBins = unique(binCenter);
    binSummaries = repmat(struct(), numel(uniqueBins), 1);
    keepBin = false(size(uniqueBins));

    for binIdx = 1:numel(uniqueBins)
        inBin = binCenter == uniqueBins(binIdx);
        sessionIDs = unique(session(inBin), 'stable');
        sessionPct = nan(size(sessionIDs));

        for sessionIdx = 1:numel(sessionIDs)
            inSession = inBin & session == sessionIDs(sessionIdx);
            sessionPct(sessionIdx) = 100 .* sum(successes(inSession)) ./ sum(nTrials(inSession));
        end

        binSummaries(binIdx).x = uniqueBins(binIdx);
        binSummaries(binIdx).successes = sum(successes(inBin));
        binSummaries(binIdx).nTrials = sum(nTrials(inBin));
        binSummaries(binIdx).pctCorrect = 100 .* ...
            binSummaries(binIdx).successes ./ binSummaries(binIdx).nTrials;
        binSummaries(binIdx).nPoints = sum(inBin);
        binSummaries(binIdx).nSessions = numel(sessionIDs);
        binSummaries(binIdx).sessionIDs = sessionIDs;
        binSummaries(binIdx).sessionValues = sessionPct;
        binSummaries(binIdx).sessionMean = mean(sessionPct, 'omitnan');
        binSummaries(binIdx).sessionMedian = median(sessionPct, 'omitnan');
        if numel(sessionPct) > 1
            binSummaries(binIdx).sessionSEM = std(sessionPct, 0, 'omitnan') ./ ...
                sqrt(sum(~isnan(sessionPct)));
        else
            binSummaries(binIdx).sessionSEM = NaN;
        end

        keepBin(binIdx) = binSummaries(binIdx).nSessions >= opts.minSessionsPerBin;
    end

    binSummaries = binSummaries(keepBin);
    if isempty(binSummaries)
        return;
    end

    summary.x = [binSummaries.x];
    summary.pctCorrect = [binSummaries.pctCorrect];
    summary.successes = [binSummaries.successes];
    summary.nTrials = [binSummaries.nTrials];
    summary.nPoints = [binSummaries.nPoints];
    summary.nSessions = [binSummaries.nSessions];
    summary.sessionMean = [binSummaries.sessionMean];
    summary.sessionMedian = [binSummaries.sessionMedian];
    summary.sessionSEM = [binSummaries.sessionSEM];
    summary.sessionIDs = {binSummaries.sessionIDs};
    summary.sessionValues = {binSummaries.sessionValues};

    if opts.useMedianForPlot
        summary.y = summary.sessionMedian;
    else
        summary.y = summary.pctCorrect;
    end
end

function summary = emptyAggregateSummary(opts)
    summary = struct( ...
        'x', [], ...
        'y', [], ...
        'pctCorrect', [], ...
        'successes', [], ...
        'nTrials', [], ...
        'nPoints', [], ...
        'nSessions', [], ...
        'sessionMean', [], ...
        'sessionMedian', [], ...
        'sessionSEM', [], ...
        'sessionIDs', {{}}, ...
        'sessionValues', {{}}, ...
        'binWidth', opts.binWidth, ...
        'useMedianForPlot', opts.useMedianForPlot, ...
        'minSessionsPerBin', opts.minSessionsPerBin);
end

function clusterIdx = findClusterIndex(clusterIDs, clusterID)
    if isnumeric(clusterIDs) || islogical(clusterIDs) || iscategorical(clusterIDs)
        clusterIdx = find(clusterIDs == clusterID, 1);
    elseif isstring(clusterIDs)
        clusterIdx = find(clusterIDs == clusterID, 1);
    else
        clusterIdx = find(strcmp(clusterIDs, clusterID), 1);
    end
end

function sessionIDs = findClusterSessions(powerCluster, clusterID)
    if isnumeric(powerCluster) || islogical(powerCluster) || iscategorical(powerCluster)
        sessionIDs = find(powerCluster == clusterID)';
    elseif isstring(powerCluster)
        sessionIDs = find(powerCluster == clusterID)';
    else
        sessionIDs = find(strcmp(powerCluster, clusterID))';
    end
end

function clusterID = getClusterLabel(clusterIDs, clusterIdx)
    if iscell(clusterIDs)
        clusterID = clusterIDs{clusterIdx};
    else
        clusterID = clusterIDs(clusterIdx);
    end
end
