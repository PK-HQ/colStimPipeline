function [powerEffectCluster, clusterSummary, diagnostics] = clusterOrderedPowerEffect(power, deltaBias, deltaMask, opts)
% Segment sessions into contiguous power bands using power and behavioral effect.

    if nargin < 3 || isempty(deltaMask)
        deltaMask = nan(size(deltaBias));
    end
    if nargin < 4 || isempty(opts)
        opts = struct();
    end
    opts = applyDefaults(opts);

    power = power(:);
    deltaBias = deltaBias(:);
    deltaMask = deltaMask(:);
    if numel(power) ~= numel(deltaBias) || numel(power) ~= numel(deltaMask)
        error('power, deltaBias, and deltaMask must have the same number of elements.');
    end

    nSessions = numel(power);
    powerEffectCluster = nan(nSessions, 1);
    valid = isfinite(power) & power >= 0 & isfinite(deltaBias);
    validIndices = find(valid);

    diagnostics = initializeDiagnostics(opts, validIndices);
    if isempty(validIndices)
        clusterSummary = emptySummary();
        diagnostics.selectionReason = 'No sessions have finite power and deltaBias.';
        return;
    end

    [sortedPower, sortOrder] = sort(power(validIndices), 'ascend');
    sortedIndices = validIndices(sortOrder);
    sortedBias = deltaBias(sortedIndices);
    sortedMask = deltaMask(sortedIndices);
    zLogPower = standardizeVector(log10(sortedPower + eps));
    zDeltaBias = standardizeVector(sortedBias);

    candidateResults = struct([]);
    for kIdx = 1:numel(opts.kCandidates)
        k = opts.kCandidates(kIdx);
        candidate = findBestCandidate(k, sortedPower, zLogPower, ...
            zDeltaBias, opts);
        if candidate.feasible
            candidateResults = [candidateResults; candidate]; %#ok<AGROW>
        end
    end

    if isempty(candidateResults)
        selected = makeSingleClusterCandidate(zLogPower, zDeltaBias);
        diagnostics.selectionReason = sprintf( ...
            'Fewer than %d sessions per group prevented a legal 2- or 3-band split.', ...
            opts.minSessionsPerCluster);
    else
        selected = selectCandidate(candidateResults, opts);
        diagnostics.selectionReason = selected.selectionReason;
    end

    sortedLabels = labelsFromBreaks(numel(sortedPower), selected.breakIndices);
    powerEffectCluster(sortedIndices) = sortedLabels;
    boundaries = boundaryValues(sortedPower, selected.breakIndices);

    diagnostics.method = 'orderedPowerEffect';
    diagnostics.chosenK = selected.k;
    diagnostics.score = selected.score;
    diagnostics.boundaries = boundaries;
    diagnostics.breakIndices = selected.breakIndices;
    diagnostics.sortedSessionIndices = sortedIndices;
    diagnostics.sortedPower = sortedPower;
    diagnostics.sortedDeltaBias = sortedBias;
    diagnostics.sortedDeltaMask = sortedMask;
    diagnostics.candidates = candidateResults;
    clusterSummary = buildSummary(powerEffectCluster, power, deltaBias, ...
        deltaMask, selected.k, selected.score);
end

function opts = applyDefaults(opts)
    defaults = struct( ...
        'minSessionsPerCluster', 3, ...
        'effectWeight', 0.75, ...
        'monotonicPenalty', 2, ...
        'kCandidates', [2 3], ...
        'minRelativeImprovementFor3Clusters', 0.10);

    names = fieldnames(defaults);
    for ii = 1:numel(names)
        name = names{ii};
        if ~isfield(opts, name) || isempty(opts.(name))
            opts.(name) = defaults.(name);
        end
    end

    validateattributes(opts.minSessionsPerCluster, {'numeric'}, ...
        {'scalar', 'integer', 'positive'});
    validateattributes(opts.effectWeight, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'nonnegative'});
    validateattributes(opts.monotonicPenalty, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'nonnegative'});
    validateattributes(opts.kCandidates, {'numeric'}, ...
        {'vector', 'integer', 'positive'});
    validateattributes(opts.minRelativeImprovementFor3Clusters, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'nonnegative'});

    opts.kCandidates = unique(opts.kCandidates(:)', 'stable');
    if any(~ismember(opts.kCandidates, [2 3]))
        error('opts.kCandidates currently supports only k = 2 and k = 3.');
    end
end

function diagnostics = initializeDiagnostics(opts, validIndices)
    diagnostics = struct( ...
        'method', 'orderedPowerEffect', ...
        'opts', opts, ...
        'validSessionIndices', validIndices, ...
        'chosenK', NaN, ...
        'score', NaN, ...
        'boundaries', [], ...
        'breakIndices', [], ...
        'sortedSessionIndices', [], ...
        'sortedPower', [], ...
        'sortedDeltaBias', [], ...
        'sortedDeltaMask', [], ...
        'candidates', struct([]), ...
        'selectionReason', '');
end

function candidate = findBestCandidate(k, sortedPower, zPower, zBias, opts)
    n = numel(sortedPower);
    candidate = emptyCandidate(k);
    if n < k * opts.minSessionsPerCluster
        return;
    end

    legalBreaks = find(sortedPower(1:end-1) < sortedPower(2:end));
    legalBreaks = legalBreaks( ...
        legalBreaks >= opts.minSessionsPerCluster & ...
        legalBreaks <= n - opts.minSessionsPerCluster);
    if isempty(legalBreaks)
        return;
    end

    bestScore = Inf;
    bestBreaks = [];
    if k == 2
        breakSets = num2cell(legalBreaks(:), 2);
    else
        breakSets = {};
        for firstIdx = 1:numel(legalBreaks)
            firstBreak = legalBreaks(firstIdx);
            for secondIdx = firstIdx + 1:numel(legalBreaks)
                secondBreak = legalBreaks(secondIdx);
                groupSizes = [firstBreak, secondBreak - firstBreak, n - secondBreak];
                if all(groupSizes >= opts.minSessionsPerCluster)
                    breakSets{end + 1, 1} = [firstBreak, secondBreak]; %#ok<AGROW>
                end
            end
        end
    end

    for splitIdx = 1:numel(breakSets)
        breaks = breakSets{splitIdx};
        labels = labelsFromBreaks(n, breaks);
        score = segmentationScore(labels, zPower, zBias, opts);
        if score < bestScore
            bestScore = score;
            bestBreaks = breaks;
        end
    end

    if ~isempty(bestBreaks)
        candidate.feasible = true;
        candidate.score = bestScore;
        candidate.breakIndices = bestBreaks;
    end
end

function candidate = emptyCandidate(k)
    candidate = struct( ...
        'k', k, ...
        'feasible', false, ...
        'score', Inf, ...
        'breakIndices', [], ...
        'selectionReason', '');
end

function selected = makeSingleClusterCandidate(zPower, zBias)
    selected = emptyCandidate(1);
    selected.feasible = true;
    selected.score = mean(zPower .^ 2, 'omitnan') + ...
        mean(zBias .^ 2, 'omitnan');
    selected.breakIndices = [];
end

function selected = selectCandidate(candidates, opts)
    kValues = [candidates.k];
    idx2 = find(kValues == 2, 1);
    idx3 = find(kValues == 3, 1);

    if isempty(idx2)
        selected = candidates(idx3);
        selected.selectionReason = 'Only the 3-band solution was feasible.';
        return;
    end

    selected = candidates(idx2);
    selected.selectionReason = 'Selected 2 bands by default.';
    if isempty(idx3)
        selected.selectionReason = 'The 3-band solution was not feasible.';
        return;
    end

    score2 = candidates(idx2).score;
    score3 = candidates(idx3).score;
    relativeImprovement = (score2 - score3) ./ max(abs(score2), eps);
    if relativeImprovement >= opts.minRelativeImprovementFor3Clusters
        selected = candidates(idx3);
        selected.selectionReason = sprintf( ...
            'Selected 3 bands because score improved by %.1f%%.', ...
            100 .* relativeImprovement);
    else
        selected.selectionReason = sprintf( ...
            'Kept 2 bands because the 3-band improvement was %.1f%%.', ...
            100 .* relativeImprovement);
    end
end

function score = segmentationScore(labels, zPower, zBias, opts)
    n = numel(labels);
    powerSSE = 0;
    effectSSE = 0;
    meanBias = nan(1, max(labels));
    for clusterID = 1:max(labels)
        inCluster = labels == clusterID;
        powerSSE = powerSSE + sum( ...
            (zPower(inCluster) - mean(zPower(inCluster))) .^ 2);
        effectSSE = effectSSE + sum( ...
            (zBias(inCluster) - mean(zBias(inCluster))) .^ 2);
        meanBias(clusterID) = mean(zBias(inCluster));
    end

    monotonicViolations = sum(diff(meanBias) < 0);
    score = powerSSE ./ n + opts.effectWeight .* effectSSE ./ n + ...
        opts.monotonicPenalty .* monotonicViolations;
end

function labels = labelsFromBreaks(n, breaks)
    labels = ones(n, 1);
    for breakIdx = 1:numel(breaks)
        labels(breaks(breakIdx) + 1:end) = breakIdx + 1;
    end
end

function boundaries = boundaryValues(sortedPower, breaks)
    boundaries = nan(1, numel(breaks));
    for ii = 1:numel(breaks)
        boundaries(ii) = mean(sortedPower(breaks(ii):breaks(ii) + 1));
    end
end

function z = standardizeVector(values)
    values = values(:);
    scale = std(values, 0, 'omitnan');
    if ~isfinite(scale) || scale == 0
        z = zeros(size(values));
    else
        z = (values - mean(values, 'omitnan')) ./ scale;
    end
end

function summary = buildSummary(labels, power, deltaBias, deltaMask, chosenK, score)
    validClusters = unique(labels(isfinite(labels)))';
    summary = repmat(emptySummaryRow(), numel(validClusters), 1);
    for row = 1:numel(validClusters)
        clusterID = validClusters(row);
        sessionIndices = find(labels == clusterID);
        summary(row).clusterID = clusterID;
        summary(row).nSessions = numel(sessionIndices);
        summary(row).sessionIndices = sessionIndices(:)';
        summary(row).powerMin = min(power(sessionIndices));
        summary(row).powerMean = mean(power(sessionIndices));
        summary(row).powerMax = max(power(sessionIndices));
        summary(row).deltaBiasMin = min(deltaBias(sessionIndices));
        summary(row).deltaBiasMean = mean(deltaBias(sessionIndices));
        summary(row).deltaBiasMax = max(deltaBias(sessionIndices));
        validMask = deltaMask(sessionIndices);
        validMask = validMask(isfinite(validMask));
        if ~isempty(validMask)
            summary(row).deltaMaskMin = min(validMask);
            summary(row).deltaMaskMean = mean(validMask);
            summary(row).deltaMaskMax = max(validMask);
        end
        summary(row).chosenK = chosenK;
        summary(row).score = score;
        summary(row).method = 'orderedPowerEffect';
    end
end

function summary = emptySummary()
    summary = repmat(emptySummaryRow(), 0, 1);
end

function row = emptySummaryRow()
    row = struct( ...
        'clusterID', NaN, ...
        'nSessions', 0, ...
        'sessionIndices', [], ...
        'powerMin', NaN, ...
        'powerMean', NaN, ...
        'powerMax', NaN, ...
        'deltaBiasMin', NaN, ...
        'deltaBiasMean', NaN, ...
        'deltaBiasMax', NaN, ...
        'deltaMaskMin', NaN, ...
        'deltaMaskMean', NaN, ...
        'deltaMaskMax', NaN, ...
        'chosenK', NaN, ...
        'score', NaN, ...
        'method', 'orderedPowerEffect');
end
