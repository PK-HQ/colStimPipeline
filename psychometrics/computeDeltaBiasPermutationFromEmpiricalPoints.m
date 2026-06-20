function result = computeDeltaBiasPermutationFromEmpiricalPoints(conContrast, conPct, inconContrast, inconPct, nPermutations, randomSeed)
% Compute raw uncorrected two-sided permutation tests for empirical deltaBias.
%
% Inputs are direct empirical marker vectors. This helper does not search for
% sessions, dates, model rows, clusters, or files.

    if nargin < 5 || isempty(nPermutations)
        nPermutations = 500;
    end
    if nargin < 6 || isempty(randomSeed)
        randomSeed = 1;
    end
    rng(double(randomSeed), 'twister');

    [conContrast, conPct] = cleanEmpiricalPoints(conContrast, conPct, 'con');
    [inconContrast, inconPct] = cleanEmpiricalPoints(inconContrast, inconPct, 'incon');
    assertUniqueContrasts(conContrast, 'con');
    assertUniqueContrasts(inconContrast, 'incon');

    matchedContrasts = intersect(conContrast, inconContrast, 'stable');
    matchedContrasts = matchedContrasts(:);
    nContrasts = numel(matchedContrasts);

    observedDelta = nan(nContrasts, 1);
    matchedConPct = nan(nContrasts, 1);
    matchedInconPct = nan(nContrasts, 1);
    nCon = nan(nContrasts, 1);
    nIncon = nan(nContrasts, 1);
    conCorrect = nan(nContrasts, 1);
    inconCorrect = nan(nContrasts, 1);
    conReconstructionErrorPct = nan(nContrasts, 1);
    inconReconstructionErrorPct = nan(nContrasts, 1);
    nullDeltaByContrast = nan(nContrasts, nPermutations);

    for contrastIdx = 1:nContrasts
        contrast = matchedContrasts(contrastIdx);
        conIdx = find(conContrast == contrast, 1);
        inconIdx = find(inconContrast == contrast, 1);
        matchedConPct(contrastIdx) = normalizePercent(conPct(conIdx), 'con');
        matchedInconPct(contrastIdx) = normalizePercent(inconPct(inconIdx), 'incon');
        observedDelta(contrastIdx) = matchedConPct(contrastIdx) - matchedInconPct(contrastIdx);

        [nCon(contrastIdx), nIncon(contrastIdx)] = assumedDeltaPermutationTrialCounts(contrast);
        conCorrect(contrastIdx) = round(matchedConPct(contrastIdx) ./ 100 .* nCon(contrastIdx));
        inconCorrect(contrastIdx) = round(matchedInconPct(contrastIdx) ./ 100 .* nIncon(contrastIdx));
        conCorrect(contrastIdx) = min(max(conCorrect(contrastIdx), 0), nCon(contrastIdx));
        inconCorrect(contrastIdx) = min(max(inconCorrect(contrastIdx), 0), nIncon(contrastIdx));

        reconstructedConPct = 100 .* conCorrect(contrastIdx) ./ nCon(contrastIdx);
        reconstructedInconPct = 100 .* inconCorrect(contrastIdx) ./ nIncon(contrastIdx);
        conReconstructionErrorPct(contrastIdx) = reconstructedConPct - matchedConPct(contrastIdx);
        inconReconstructionErrorPct(contrastIdx) = reconstructedInconPct - matchedInconPct(contrastIdx);

        totalCorrect = conCorrect(contrastIdx) + inconCorrect(contrastIdx);
        totalTrials = nCon(contrastIdx) + nIncon(contrastIdx);
        pooledOutcomes = [ones(totalCorrect, 1); zeros(totalTrials - totalCorrect, 1)];
        for permIdx = 1:nPermutations
            permOrder = randperm(totalTrials);
            permConCorrect = sum(pooledOutcomes(permOrder(1:nCon(contrastIdx))));
            permInconCorrect = totalCorrect - permConCorrect;
            nullConPct = 100 .* permConCorrect ./ nCon(contrastIdx);
            nullInconPct = 100 .* permInconCorrect ./ nIncon(contrastIdx);
            nullDeltaByContrast(contrastIdx, permIdx) = nullConPct - nullInconPct;
        end
    end

    nullMedian = nan(nContrasts, 1);
    nullLower95 = nan(nContrasts, 1);
    nullUpper95 = nan(nContrasts, 1);
    rawTwoSidedP = nan(nContrasts, 1);
    rawPositiveOneSidedP = nan(nContrasts, 1);
    for contrastIdx = 1:nContrasts
        nullVals = nullDeltaByContrast(contrastIdx, :);
        nullMedian(contrastIdx) = median(nullVals, 'omitnan');
        nullLower95(contrastIdx) = prctile(nullVals, 2.5);
        nullUpper95(contrastIdx) = prctile(nullVals, 97.5);
        pUpper = (1 + sum(nullVals >= observedDelta(contrastIdx))) ./ (nPermutations + 1);
        pLower = (1 + sum(nullVals <= observedDelta(contrastIdx))) ./ (nPermutations + 1);
        rawTwoSidedP(contrastIdx) = min(1, 2 .* min(pUpper, pLower));
        rawPositiveOneSidedP(contrastIdx) = pUpper;
    end

    if nContrasts > 0
        observedMeanDelta = mean(observedDelta, 'omitnan');
        nullMeanDelta = mean(nullDeltaByContrast, 1, 'omitnan');
        pUpperOverall = (1 + sum(nullMeanDelta >= observedMeanDelta)) ./ (nPermutations + 1);
        pLowerOverall = (1 + sum(nullMeanDelta <= observedMeanDelta)) ./ (nPermutations + 1);
        rawOverallTwoSidedP = min(1, 2 .* min(pUpperOverall, pLowerOverall));
        rawOverallPositiveOneSidedP = pUpperOverall;
        nullMeanMedian = median(nullMeanDelta, 'omitnan');
        nullMeanLower95 = prctile(nullMeanDelta, 2.5);
        nullMeanUpper95 = prctile(nullMeanDelta, 97.5);
    else
        observedMeanDelta = NaN;
        nullMeanDelta = nan(nPermutations, 1);
        rawOverallTwoSidedP = NaN;
        rawOverallPositiveOneSidedP = NaN;
        nullMeanMedian = NaN;
        nullMeanLower95 = NaN;
        nullMeanUpper95 = NaN;
    end

    result = struct();
    result.contrast = matchedContrasts;
    result.conPct = matchedConPct;
    result.inconPct = matchedInconPct;
    result.observedDeltaBias = observedDelta;
    result.nullDeltaByContrast = nullDeltaByContrast;
    result.nullMedian = nullMedian;
    result.nullLower95 = nullLower95;
    result.nullUpper95 = nullUpper95;
    result.rawTwoSidedP = rawTwoSidedP;
    result.rawPositiveOneSidedP = rawPositiveOneSidedP;
    result.significantUncorrected = rawTwoSidedP < 0.05;
    result.nCon = nCon;
    result.nIncon = nIncon;
    result.conCorrect = conCorrect;
    result.inconCorrect = inconCorrect;
    result.conReconstructionErrorPct = conReconstructionErrorPct;
    result.inconReconstructionErrorPct = inconReconstructionErrorPct;
    result.observedMeanDeltaBias = observedMeanDelta;
    result.nullMeanDelta = nullMeanDelta(:);
    result.nullMeanMedian = nullMeanMedian;
    result.nullMeanLower95 = nullMeanLower95;
    result.nullMeanUpper95 = nullMeanUpper95;
    result.rawOverallTwoSidedP = rawOverallTwoSidedP;
    result.rawOverallPositiveOneSidedP = rawOverallPositiveOneSidedP;
    result.overallSignificant = rawOverallTwoSidedP < 0.05;
    result.nPermutations = nPermutations;
    result.multipleComparisonCorrection = 'none';
end

function [x, y] = cleanEmpiricalPoints(x, y, label)
    x = x(:);
    y = y(:);
    keep = isfinite(x) & isfinite(y);
    x = x(keep);
    y = y(keep);
    [x, order] = sort(x);
    y = y(order);
    if isempty(x)
        warning('computeDeltaBiasPermutationFromEmpiricalPoints:EmptyPoints', ...
            'No finite %s empirical marker points found.', label);
    end
end

function assertUniqueContrasts(x, label)
    if numel(unique(x)) ~= numel(x)
        error('computeDeltaBiasPermutationFromEmpiricalPoints:DuplicateContrasts', ...
            'Duplicate %s empirical marker contrasts make exact matching ambiguous.', label);
    end
end

function pct = normalizePercent(value, label)
    pct = value;
    if abs(pct) <= 1.5
        pct = 100 .* pct;
    end
    if pct < -1e-8 || pct > 100 + 1e-8
        error('computeDeltaBiasPermutationFromEmpiricalPoints:InvalidPercent', ...
            '%s empirical value is outside [0, 100] after normalization: %g', label, pct);
    end
end

function [nCon, nIncon] = assumedDeltaPermutationTrialCounts(contrast)
    if abs(contrast) < eps
        nCon = 40;
        nIncon = 40;
    else
        nCon = 20;
        nIncon = 20;
    end
end
