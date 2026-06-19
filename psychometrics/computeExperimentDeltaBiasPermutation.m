function result = computeExperimentDeltaBiasPermutation(conContrast, conPct, inconContrast, inconPct, opts)
% Compute an exact-contrast con-vs-incon empirical deltaBias permutation test.
%
% Inputs are the saved empirical marker arrays from mdl.xBlock/yBlock. Baseline
% values, fitted curves, and tolerance-based contrast matching are not used.

    if nargin < 5 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'nPermutations') || isempty(opts.nPermutations)
        opts.nPermutations = 500;
    end
    nPermutations = opts.nPermutations;

    [conContrast, conPct] = cleanMarkerData(conContrast, conPct, 'con');
    [inconContrast, inconPct] = cleanMarkerData(inconContrast, inconPct, 'incon');

    assertNoDuplicateContrasts(conContrast, 'con');
    assertNoDuplicateContrasts(inconContrast, 'incon');

    matchedContrasts = intersect(conContrast, inconContrast);
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
        observedDelta(contrastIdx) = matchedConPct(contrastIdx) - ...
            matchedInconPct(contrastIdx);

        [nCon(contrastIdx), nIncon(contrastIdx)] = assumedTrialCounts(contrast);
        conCorrect(contrastIdx) = round(matchedConPct(contrastIdx) ./ 100 .* nCon(contrastIdx));
        inconCorrect(contrastIdx) = round(matchedInconPct(contrastIdx) ./ 100 .* nIncon(contrastIdx));

        reconstructedConPct = 100 .* conCorrect(contrastIdx) ./ nCon(contrastIdx);
        reconstructedInconPct = 100 .* inconCorrect(contrastIdx) ./ nIncon(contrastIdx);
        conReconstructionErrorPct(contrastIdx) = reconstructedConPct - ...
            matchedConPct(contrastIdx);
        inconReconstructionErrorPct(contrastIdx) = reconstructedInconPct - ...
            matchedInconPct(contrastIdx);

        totalCorrect = conCorrect(contrastIdx) + inconCorrect(contrastIdx);
        totalTrials = nCon(contrastIdx) + nIncon(contrastIdx);
        pooledOutcomes = [ones(totalCorrect, 1); ...
            zeros(totalTrials - totalCorrect, 1)];

        for permIdx = 1:nPermutations
            permOrder = randperm(totalTrials);
            permConCorrect = sum(pooledOutcomes(permOrder(1:nCon(contrastIdx))));
            permInconCorrect = totalCorrect - permConCorrect;
            nullConPct = 100 .* permConCorrect ./ nCon(contrastIdx);
            nullInconPct = 100 .* permInconCorrect ./ nIncon(contrastIdx);
            nullDeltaByContrast(contrastIdx, permIdx) = nullConPct - nullInconPct;
        end
    end

    nullMeanDelta = mean(nullDeltaByContrast, 1, 'omitnan');
    observedMeanDelta = mean(observedDelta, 'omitnan');
    if isempty(nullMeanDelta) || all(isnan(nullMeanDelta)) || isnan(observedMeanDelta)
        pExperiment = NaN;
    else
        pExperiment = (1 + sum(nullMeanDelta >= observedMeanDelta)) ./ ...
            (nPermutations + 1);
    end

    contrastP = nan(nContrasts, 1);
    nullDeltaMean = nan(nContrasts, 1);
    nullDeltaMedian = nan(nContrasts, 1);
    nullDeltaLower95 = nan(nContrasts, 1);
    nullDeltaUpper95 = nan(nContrasts, 1);
    for contrastIdx = 1:nContrasts
        nullVals = nullDeltaByContrast(contrastIdx, :);
        contrastP(contrastIdx) = (1 + sum(nullVals >= observedDelta(contrastIdx))) ./ ...
            (nPermutations + 1);
        nullDeltaMean(contrastIdx) = mean(nullVals, 'omitnan');
        nullDeltaMedian(contrastIdx) = median(nullVals, 'omitnan');
        nullDeltaLower95(contrastIdx) = prctile(nullVals, 2.5);
        nullDeltaUpper95(contrastIdx) = prctile(nullVals, 97.5);
    end

    result = struct();
    result.matchedContrasts = matchedContrasts;
    result.conPct = matchedConPct;
    result.inconPct = matchedInconPct;
    result.observedDelta = observedDelta;
    result.nCon = nCon;
    result.nIncon = nIncon;
    result.conCorrect = conCorrect;
    result.inconCorrect = inconCorrect;
    result.conReconstructionErrorPct = conReconstructionErrorPct;
    result.inconReconstructionErrorPct = inconReconstructionErrorPct;
    result.nullDeltaMean = nullDeltaMean;
    result.nullDeltaMedian = nullDeltaMedian;
    result.nullDeltaLower95 = nullDeltaLower95;
    result.nullDeltaUpper95 = nullDeltaUpper95;
    result.oneSidedContrastP = contrastP;
    result.observedMeanDelta = observedMeanDelta;
    result.nullMeanDelta = nullMeanDelta(:);
    result.nullMeanDeltaMean = mean(nullMeanDelta, 'omitnan');
    result.nullMeanDeltaMedian = median(nullMeanDelta, 'omitnan');
    result.nullMeanDeltaLower95 = prctile(nullMeanDelta, 2.5);
    result.nullMeanDeltaUpper95 = prctile(nullMeanDelta, 97.5);
    result.oneSidedExperimentP = pExperiment;
    result.experimentSignificant = observedMeanDelta > 0 && ...
        isfinite(pExperiment) && pExperiment < 0.05;
    result.nPermutations = nPermutations;
    result.sourceFields = struct( ...
        'conContrast', 'mdl.xBlock(2,:,sessionRowIndex)', ...
        'conPercentCorrect', 'mdl.yBlock(2,:,sessionRowIndex)', ...
        'inconContrast', 'mdl.xBlock(3,:,sessionRowIndex)', ...
        'inconPercentCorrect', 'mdl.yBlock(3,:,sessionRowIndex)');
end

function [x, y] = cleanMarkerData(x, y, label)
    x = x(:);
    y = y(:);
    keep = isfinite(x) & isfinite(y);
    x = x(keep);
    y = y(keep);
    [x, order] = sort(x);
    y = y(order);
    if isempty(x)
        warning('computeExperimentDeltaBiasPermutation:EmptyMarkers', ...
            'No finite %s empirical marker points found.', label);
    end
end

function assertNoDuplicateContrasts(x, label)
    if numel(unique(x)) ~= numel(x)
        error('computeExperimentDeltaBiasPermutation:DuplicateContrasts', ...
            'Duplicate %s empirical marker contrasts found; exact matching is ambiguous.', ...
            label);
    end
end

function pct = normalizePercent(value, label)
    pct = value;
    if abs(pct) <= 1.5
        pct = 100 .* pct;
    end
    if pct < -1e-8 || pct > 100 + 1e-8
        error('computeExperimentDeltaBiasPermutation:InvalidPercentRange', ...
            '%s empirical value is outside [0, 100] after normalization: %g', ...
            label, pct);
    end
end

function [nCon, nIncon] = assumedTrialCounts(contrast)
    if abs(contrast) < eps
        nCon = 40;
        nIncon = 40;
    else
        nCon = 20;
        nIncon = 20;
    end
end
