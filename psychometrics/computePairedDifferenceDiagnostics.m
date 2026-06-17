function result = computePairedDifferenceDiagnostics( ...
        x, y, annotationTail, plannedTail)
% Compute paired-test diagnostics without applying multiple-test correction.

    x = x(:);
    y = y(:);
    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);
    differences = x - y;

    result = emptyResult();
    result.xValues = x;
    result.yValues = y;
    result.differences = differences;
    result.n = numel(differences);
    result.nPositive = sum(differences > 0);
    result.nNegative = sum(differences < 0);
    result.nZero = sum(differences == 0);
    result.meanDifference = mean(differences, 'omitnan');
    result.medianDifference = median(differences, 'omitnan');
    if result.n >= 2
        result.semDifference = std(differences, 0, 'omitnan') ./ sqrt(result.n);
    elseif result.n == 1
        result.semDifference = 0;
    end

    result.annotationTail = annotationTail;
    result.plannedTail = plannedTail;
    result.annotationAlternative = tailDescription(annotationTail);
    result.plannedAlternative = tailDescription(plannedTail);

    result.nWilcoxonNonzero = sum(differences ~= 0);
    [result.minPossibleTwoSidedRawP, ...
        result.minPossibleOneSidedRawP] = ...
        minimumWilcoxonP(result.nWilcoxonNonzero);
    if strcmp(annotationTail, 'both')
        result.minPossibleAnnotationRawP = ...
            result.minPossibleTwoSidedRawP;
    else
        result.minPossibleAnnotationRawP = ...
            result.minPossibleOneSidedRawP;
    end

    [result.rawP, result.test] = ...
        annotationTest(differences, annotationTail);
    [result.pairedTTestRawP, result.tStatistic] = ...
        runTTest(differences, annotationTail);
    result.pairedTTestN = result.n;
    result.plannedRawP = runWilcoxon(differences, plannedTail);
end

function [pValue, testName] = annotationTest(differences, tail)
    [pValue, available] = runWilcoxon(differences, tail);
    if available
        testName = 'Wilcoxon signed-rank';
        return;
    end

    [pValue, ~] = runTTest(differences, tail);
    if isfinite(pValue)
        testName = 'paired t-test fallback';
    elseif numel(differences) < 2
        testName = 'insufficient paired data';
    else
        testName = 'signrank and ttest unavailable';
    end
end

function [pValue, available] = runWilcoxon(differences, tail)
    pValue = nan;
    available = false;
    if numel(differences) < 2
        return;
    end
    if all(differences == 0)
        pValue = 1;
        available = true;
        return;
    end
    if exist('signrank', 'file') ~= 2
        return;
    end
    try
        pValue = signrank(differences, 0, 'tail', tail);
        available = true;
    catch
    end
end

function [pValue, tStatistic] = runTTest(differences, tail)
    pValue = nan;
    tStatistic = nan;
    if numel(differences) < 2 || exist('ttest', 'file') ~= 2
        return;
    end
    try
        [~, pValue, ~, stats] = ttest(differences, 0, 'Tail', tail);
        if isfield(stats, 'tstat')
            tStatistic = stats.tstat;
        end
    catch
    end
end

function [twoSidedP, oneSidedP] = minimumWilcoxonP(nNonzero)
    if nNonzero < 1
        twoSidedP = nan;
        oneSidedP = nan;
        return;
    end
    oneSidedP = 2 .^ (-nNonzero);
    twoSidedP = min(1, 2 .* oneSidedP);
end

function description = tailDescription(tail)
    switch tail
        case 'right'
            description = 'x > y';
        case 'left'
            description = 'x < y';
        otherwise
            description = 'x ~= y';
    end
end

function result = emptyResult()
    result = struct( ...
        'xValues', [], ...
        'yValues', [], ...
        'differences', [], ...
        'n', 0, ...
        'nPositive', 0, ...
        'nNegative', 0, ...
        'nZero', 0, ...
        'meanDifference', nan, ...
        'medianDifference', nan, ...
        'semDifference', nan, ...
        'annotationTail', 'both', ...
        'annotationAlternative', 'x ~= y', ...
        'plannedTail', 'both', ...
        'plannedAlternative', 'x ~= y', ...
        'test', '', ...
        'rawP', nan, ...
        'adjustedP', nan, ...
        'star', 'n.s.', ...
        'familySize', nan, ...
        'nWilcoxonNonzero', 0, ...
        'minPossibleTwoSidedRawP', nan, ...
        'minPossibleOneSidedRawP', nan, ...
        'minPossibleAnnotationRawP', nan, ...
        'minPossibleAnnotationAdjustedP', nan, ...
        'pairedTTestRawP', nan, ...
        'pairedTTestAdjustedP', nan, ...
        'pairedTTestN', 0, ...
        'tStatistic', nan, ...
        'plannedRawP', nan, ...
        'plannedAdjustedP', nan, ...
        'plannedStar', 'n.s.');
end
