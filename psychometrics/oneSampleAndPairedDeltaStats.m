function results = oneSampleAndPairedDeltaStats(deltaBias, deltaMask, opts)
% Test bias and mask against zero, plus their paired difference.

    if nargin < 3 || isempty(opts)
        opts = struct();
    end
    opts = applyDefaults(opts);

    deltaBias = deltaBias(:);
    deltaMask = deltaMask(:);
    template = computePairedDifferenceDiagnostics([], [], 'both', 'both');
    template.groups = [nan nan];
    template.validRows = [];
    template.comparison = '';
    template.alternativeHypothesis = '';
    template.plannedAlternativeHypothesis = '';
    template.testDirection = opts.testDirection;
    results = repmat(template, 3, 1);

    validBias = isfinite(deltaBias);
    results(1) = makeResult(deltaBias(validBias), ...
        zeros(sum(validBias), 1), [1 1], find(validBias)', ...
        '\DeltaBias vs 0', 'right', ...
        '\DeltaBias > 0', '\DeltaBias ~= 0', opts);

    validMask = isfinite(deltaMask);
    results(2) = makeResult(deltaMask(validMask), ...
        zeros(sum(validMask), 1), [2 2], find(validMask)', ...
        '\DeltaMask vs 0', 'both', ...
        '\DeltaMask ~= 0', '\DeltaMask ~= 0', opts);

    validPair = validBias & validMask;
    results(3) = makeResult( ...
        deltaBias(validPair), deltaMask(validPair), ...
        [1 2], find(validPair)', ...
        '\DeltaBias vs \DeltaMask', 'both', ...
        '\DeltaBias ~= \DeltaMask', ...
        '\DeltaBias ~= \DeltaMask', opts);

    rawP = [results.rawP]';
    tTestRawP = [results.pairedTTestRawP]';
    plannedRawP = [results.plannedRawP]';
    adjustedP = holmBonferroni(rawP);
    tTestAdjustedP = holmBonferroni(tTestRawP);
    plannedAdjustedP = holmBonferroni(plannedRawP);
    for resultIdx = 1:3
        results(resultIdx).adjustedP = adjustedP(resultIdx);
        results(resultIdx).pairedTTestAdjustedP = ...
            tTestAdjustedP(resultIdx);
        results(resultIdx).plannedAdjustedP = ...
            plannedAdjustedP(resultIdx);
        results(resultIdx).familySize = 3;
        results(resultIdx).minPossibleAnnotationAdjustedP = min( ...
            1, 3 .* results(resultIdx).minPossibleAnnotationRawP);
        results(resultIdx).star = pValueStars( ...
            adjustedP(resultIdx), opts.alpha);
        results(resultIdx).plannedStar = pValueStars( ...
            plannedAdjustedP(resultIdx), opts.alpha);
    end
end

function result = makeResult(x, y, groups, validRows, comparison, ...
        plannedTail, plannedAlternative, twoSidedAlternative, opts)
    if strcmp(opts.testDirection, 'one-sided')
        annotationTail = plannedTail;
        annotationAlternative = plannedAlternative;
    else
        annotationTail = 'both';
        annotationAlternative = twoSidedAlternative;
    end
    result = computePairedDifferenceDiagnostics( ...
        x, y, annotationTail, plannedTail);
    result.groups = groups;
    result.validRows = validRows;
    result.comparison = comparison;
    result.alternativeHypothesis = annotationAlternative;
    result.plannedAlternativeHypothesis = plannedAlternative;
    result.testDirection = opts.testDirection;
end

function opts = applyDefaults(opts)
    if ~isfield(opts, 'alpha') || isempty(opts.alpha)
        opts.alpha = 0.05;
    end
    if ~isfield(opts, 'testDirection') || isempty(opts.testDirection)
        opts.testDirection = 'two-sided';
    end
    opts.testDirection = validatestring( ...
        opts.testDirection, {'two-sided', 'one-sided'});
end

function star = pValueStars(pValue, alpha)
    if ~isfinite(pValue) || pValue >= alpha
        star = 'n.s.';
    elseif pValue < 1e-4
        star = '****';
    elseif pValue < 1e-3
        star = '***';
    elseif pValue < 1e-2
        star = '**';
    else
        star = '*';
    end
end
