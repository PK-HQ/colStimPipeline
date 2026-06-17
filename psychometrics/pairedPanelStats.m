function results = pairedPanelStats(valuesByGroup, comparisonPairs, opts)
% Run paired comparisons with Holm correction within one plot panel.

    if nargin < 3 || isempty(opts)
        opts = struct();
    end
    opts = applyDefaults(opts, comparisonPairs);

    nComparisons = size(comparisonPairs, 1);
    template = computePairedDifferenceDiagnostics([], [], 'both', 'both');
    template.groups = [nan nan];
    template.validRows = [];
    template.comparison = '';
    template.alternativeHypothesis = '';
    template.plannedAlternativeHypothesis = '';
    template.testDirection = opts.testDirection;
    results = repmat(template, nComparisons, 1);

    rawP = nan(nComparisons, 1);
    tTestRawP = nan(nComparisons, 1);
    plannedRawP = nan(nComparisons, 1);
    for comparisonIdx = 1:nComparisons
        groupA = comparisonPairs(comparisonIdx, 1);
        groupB = comparisonPairs(comparisonIdx, 2);
        valid = isfinite(valuesByGroup(:, groupA)) & ...
            isfinite(valuesByGroup(:, groupB));
        x = valuesByGroup(valid, groupA);
        y = valuesByGroup(valid, groupB);

        plannedTail = opts.plannedTails{comparisonIdx};
        if strcmp(opts.testDirection, 'one-sided')
            annotationTail = plannedTail;
        else
            annotationTail = 'both';
        end
        diagnostic = computePairedDifferenceDiagnostics( ...
            x, y, annotationTail, plannedTail);
        diagnostic.groups = [groupA groupB];
        diagnostic.validRows = find(valid)';
        diagnostic.comparison = opts.comparisonLabels{comparisonIdx};
        diagnostic.plannedAlternativeHypothesis = ...
            opts.plannedAlternativeHypotheses{comparisonIdx};
        if strcmp(annotationTail, 'both')
            diagnostic.alternativeHypothesis = ...
                opts.twoSidedAlternativeHypotheses{comparisonIdx};
        else
            diagnostic.alternativeHypothesis = ...
                diagnostic.plannedAlternativeHypothesis;
        end
        diagnostic.testDirection = opts.testDirection;
        results(comparisonIdx) = diagnostic;
        rawP(comparisonIdx) = diagnostic.rawP;
        tTestRawP(comparisonIdx) = diagnostic.pairedTTestRawP;
        plannedRawP(comparisonIdx) = diagnostic.plannedRawP;
    end

    adjustedP = holmBonferroni(rawP);
    tTestAdjustedP = holmBonferroni(tTestRawP);
    plannedAdjustedP = holmBonferroni(plannedRawP);
    for comparisonIdx = 1:nComparisons
        results(comparisonIdx).adjustedP = adjustedP(comparisonIdx);
        results(comparisonIdx).pairedTTestAdjustedP = ...
            tTestAdjustedP(comparisonIdx);
        results(comparisonIdx).plannedAdjustedP = ...
            plannedAdjustedP(comparisonIdx);
        results(comparisonIdx).familySize = nComparisons;
        results(comparisonIdx).minPossibleAnnotationAdjustedP = min( ...
            1, nComparisons .* ...
            results(comparisonIdx).minPossibleAnnotationRawP);
        results(comparisonIdx).star = pValueStars( ...
            adjustedP(comparisonIdx), opts.alpha);
        results(comparisonIdx).plannedStar = pValueStars( ...
            plannedAdjustedP(comparisonIdx), opts.alpha);
    end
end

function opts = applyDefaults(opts, comparisonPairs)
    nComparisons = size(comparisonPairs, 1);
    if ~isfield(opts, 'alpha') || isempty(opts.alpha)
        opts.alpha = 0.05;
    end
    if ~isfield(opts, 'testDirection') || isempty(opts.testDirection)
        opts.testDirection = 'two-sided';
    end
    opts.testDirection = validatestring( ...
        opts.testDirection, {'two-sided', 'one-sided'});

    if ~isfield(opts, 'plannedTails') || isempty(opts.plannedTails)
        opts.plannedTails = repmat({'both'}, 1, nComparisons);
    end
    if ~isfield(opts, 'comparisonLabels') || isempty(opts.comparisonLabels)
        opts.comparisonLabels = defaultComparisonLabels(comparisonPairs);
    end
    if ~isfield(opts, 'plannedAlternativeHypotheses') || ...
            isempty(opts.plannedAlternativeHypotheses)
        opts.plannedAlternativeHypotheses = ...
            repmat({'two-sided; no planned direction'}, 1, nComparisons);
    end
    if ~isfield(opts, 'twoSidedAlternativeHypotheses') || ...
            isempty(opts.twoSidedAlternativeHypotheses)
        opts.twoSidedAlternativeHypotheses = ...
            cellfun(@(label) [label ' differs'], ...
            opts.comparisonLabels, 'UniformOutput', false);
    end

    validateCellOption(opts.plannedTails, nComparisons, 'plannedTails');
    validateCellOption( ...
        opts.comparisonLabels, nComparisons, 'comparisonLabels');
    validateCellOption(opts.plannedAlternativeHypotheses, ...
        nComparisons, 'plannedAlternativeHypotheses');
    validateCellOption(opts.twoSidedAlternativeHypotheses, ...
        nComparisons, 'twoSidedAlternativeHypotheses');
    for idx = 1:nComparisons
        opts.plannedTails{idx} = validatestring( ...
            opts.plannedTails{idx}, {'both', 'left', 'right'});
    end
end

function labels = defaultComparisonLabels(comparisonPairs)
    labels = cell(1, size(comparisonPairs, 1));
    for idx = 1:size(comparisonPairs, 1)
        labels{idx} = sprintf('Group %d vs Group %d', ...
            comparisonPairs(idx, 1), comparisonPairs(idx, 2));
    end
end

function validateCellOption(value, expectedLength, optionName)
    if ~iscell(value) || numel(value) ~= expectedLength
        error('pairedPanelStats:InvalidOption', ...
            '%s must be a cell array with %d entries.', ...
            optionName, expectedLength);
    end
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
