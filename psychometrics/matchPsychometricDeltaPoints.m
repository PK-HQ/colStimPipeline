function [xBias, yBias, xMask, yMask, audit] = ...
        matchPsychometricDeltaPoints(sideData, opts)
% Match side-specific psychometric points for delta plots.
%
% Scientific defaults:
%   * horizontal and vertical data are matched separately by the caller
%   * zero contrast is structural and only matches exact zero
%   * nonzero contrasts use exact matches first, then optimal one-to-one
%     nearest matches within opts.contrastPairTolerance
%   * plotted x positions are the mean of matched contrast coordinates

    if nargin < 2 || isempty(opts)
        opts = struct();
    end
    opts = applyDefaults(opts);

    [xBase, yBase] = collapseConditionDataLocal( ...
        sideData.xBaseline, sideData.yBaseline);
    [xCon, yCon] = collapseConditionDataLocal( ...
        sideData.xConOpto, sideData.yConOpto);
    [xIncon, yIncon] = collapseConditionDataLocal( ...
        sideData.xInconOpto, sideData.yInconOpto);

    [biasPairs, biasMethods, biasUnmatched] = pairTwoConditions( ...
        xCon, xIncon, opts.contrastPairTolerance, ...
        opts.nearMatchTolerance);
    xBias = nan(1, size(biasPairs, 1));
    yBias = nan(1, size(biasPairs, 1));
    for idx = 1:size(biasPairs, 1)
        conIdx = biasPairs(idx, 1);
        inconIdx = biasPairs(idx, 2);
        xBias(idx) = mean([xCon(conIdx), xIncon(inconIdx)], 'omitnan');
        yBias(idx) = yCon(conIdx) - yIncon(inconIdx);
    end

    [maskTriples, maskMethods, maskUnmatched] = pairMaskTriples( ...
        xBase, xCon, xIncon, opts.contrastPairTolerance, ...
        opts.nearMatchTolerance);
    xMask = nan(1, size(maskTriples, 1));
    yMask = nan(1, size(maskTriples, 1));
    for idx = 1:size(maskTriples, 1)
        baseIdx = maskTriples(idx, 1);
        conIdx = maskTriples(idx, 2);
        inconIdx = maskTriples(idx, 3);
        xMask(idx) = mean([xBase(baseIdx), xCon(conIdx), ...
            xIncon(inconIdx)], 'omitnan');
        yMask(idx) = yBase(baseIdx) - ...
            mean([yCon(conIdx), yIncon(inconIdx)], 'omitnan');
    end

    audit = buildDeltaAuditRows(sideData.sideName, xBase, xCon, xIncon, ...
        biasPairs, xBias, yBias, biasMethods, biasUnmatched, ...
        maskTriples, xMask, yMask, maskMethods, maskUnmatched, ...
        opts.contrastPairTolerance, opts.nearMatchTolerance);
end

function opts = applyDefaults(opts)
    if ~isfield(opts, 'contrastPairTolerance') || ...
            isempty(opts.contrastPairTolerance)
        opts.contrastPairTolerance = 5;
    end
    if ~isfield(opts, 'nearMatchTolerance') || ...
            isempty(opts.nearMatchTolerance)
        opts.nearMatchTolerance = 10;
    end
end

function [pairs, pairMethods, unmatched] = pairTwoConditions( ...
        xCon, xIncon, tolerance, nearTolerance)
    xCon = xCon(:);
    xIncon = xIncon(:);
    pairs = zeros(0, 2);
    pairMethods = strings(0, 1);
    unmatched = emptyUnmatchedStruct();
    if isempty(xCon) || isempty(xIncon)
        unmatched = addBiasUnmatchedRows(unmatched, xCon, xIncon, ...
            true(size(xCon)), true(size(xIncon)), tolerance, nearTolerance);
        return;
    end

    usedCon = false(numel(xCon), 1);
    usedIncon = false(numel(xIncon), 1);

    [zeroPairs, zeroMethods, usedCon, usedIncon, unmatched] = ...
        matchBiasZeroPoint(xCon, xIncon, usedCon, usedIncon, unmatched, ...
        tolerance, nearTolerance);
    pairs = [pairs; zeroPairs];
    pairMethods = [pairMethods; zeroMethods];

    [exactPairs, exactMethods, usedCon, usedIncon] = matchExactPairs( ...
        xCon, xIncon, usedCon, usedIncon);
    pairs = [pairs; exactPairs];
    pairMethods = [pairMethods; exactMethods];

    remainingCon = find(~usedCon & ~isZeroContrast(xCon));
    remainingIncon = find(~usedIncon & ~isZeroContrast(xIncon));
    tolerancePairs = optimizePairMatches(xCon, xIncon, remainingCon, ...
        remainingIncon, tolerance);
    if ~isempty(tolerancePairs)
        usedCon(tolerancePairs(:, 1)) = true;
        usedIncon(tolerancePairs(:, 2)) = true;
        pairs = [pairs; tolerancePairs]; %#ok<AGROW>
        pairMethods = [pairMethods; ...
            repmat("nearest/tolerance", size(tolerancePairs, 1), 1)]; %#ok<AGROW>
    end

    [pairs, sortIdx] = sortRowsByMeanContrast(pairs, xCon, xIncon);
    pairMethods = pairMethods(sortIdx);
    unmatched = addBiasUnmatchedRows(unmatched, xCon, xIncon, ...
        ~usedCon, ~usedIncon, tolerance, nearTolerance);
end

function [zeroPairs, zeroMethods, usedCon, usedIncon, unmatched] = ...
        matchBiasZeroPoint(xCon, xIncon, usedCon, usedIncon, unmatched, ...
        tolerance, nearTolerance)
    zeroPairs = zeros(0, 2);
    zeroMethods = strings(0, 1);
    conZero = find(isZeroContrast(xCon));
    inconZero = find(isZeroContrast(xIncon));
    nZero = min(numel(conZero), numel(inconZero));
    if nZero > 0
        zeroPairs = [conZero(1:nZero), inconZero(1:nZero)];
        zeroMethods = repmat("exact/zero", nZero, 1);
        usedCon(conZero(1:nZero)) = true;
        usedIncon(inconZero(1:nZero)) = true;
    end
    if numel(conZero) > nZero
        for idx = (nZero + 1):numel(conZero)
            unmatched = appendUnmatched(unmatched, "deltaBias", "con", ...
                NaN, conZero(idx), NaN, "incon zero contrast unavailable", ...
                0, tolerance, nearTolerance, false);
            usedCon(conZero(idx)) = true;
        end
    end
    if numel(inconZero) > nZero
        for idx = (nZero + 1):numel(inconZero)
            unmatched = appendUnmatched(unmatched, "deltaBias", "incon", ...
                NaN, NaN, inconZero(idx), "con zero contrast unavailable", ...
                0, tolerance, nearTolerance, false);
            usedIncon(inconZero(idx)) = true;
        end
    end
end

function [pairs, methods, usedA, usedB] = matchExactPairs( ...
        xA, xB, usedA, usedB)
    pairs = zeros(0, 2);
    methods = strings(0, 1);
    values = unique(xA(~usedA & ~isZeroContrast(xA)), 'stable');
    for valueIdx = 1:numel(values)
        idxA = find(~usedA & xA == values(valueIdx));
        idxB = find(~usedB & xB == values(valueIdx));
        nMatch = min(numel(idxA), numel(idxB));
        if nMatch == 0
            continue;
        end
        newPairs = [idxA(1:nMatch), idxB(1:nMatch)];
        pairs = [pairs; newPairs]; %#ok<AGROW>
        methods = [methods; repmat("exact", nMatch, 1)]; %#ok<AGROW>
        usedA(idxA(1:nMatch)) = true;
        usedB(idxB(1:nMatch)) = true;
    end
end

function pairs = optimizePairMatches(xA, xB, idxA, idxB, tolerance)
    pairs = zeros(0, 2);
    if isempty(idxA) || isempty(idxB)
        return;
    end

    bestPairs = zeros(0, 2);
    bestCount = -Inf;
    bestCost = Inf;
    idxA = idxA(:);
    idxB = idxB(:);
    usedB = false(numel(xB), 1);

    recurse(1, usedB, zeros(0, 2), 0);
    pairs = bestPairs;

    function recurse(pos, usedBLocal, currentPairs, currentCost)
        remainingA = numel(idxA) - pos + 1;
        if size(currentPairs, 1) + remainingA < bestCount
            return;
        end
        if pos > numel(idxA)
            currentCount = size(currentPairs, 1);
            if currentCount > bestCount || ...
                    (currentCount == bestCount && currentCost < bestCost)
                bestCount = currentCount;
                bestCost = currentCost;
                bestPairs = currentPairs;
            end
            return;
        end

        a = idxA(pos);
        candidates = idxB(~usedBLocal(idxB) & ...
            abs(xB(idxB) - xA(a)) <= tolerance);
        [~, order] = sort(abs(xB(candidates) - xA(a)));
        candidates = candidates(order);
        for candidateIdx = 1:numel(candidates)
            b = candidates(candidateIdx);
            nextUsedB = usedBLocal;
            nextUsedB(b) = true;
            recurse(pos + 1, nextUsedB, [currentPairs; a, b], ...
                currentCost + abs(xA(a) - xB(b))); %#ok<AGROW>
        end
        recurse(pos + 1, usedBLocal, currentPairs, currentCost);
    end
end

function [triples, tripleMethods, unmatched] = pairMaskTriples( ...
        xBase, xCon, xIncon, tolerance, nearTolerance)
    xBase = xBase(:);
    xCon = xCon(:);
    xIncon = xIncon(:);
    triples = zeros(0, 3);
    tripleMethods = strings(0, 1);
    unmatched = emptyUnmatchedStruct();
    if isempty(xBase) || isempty(xCon) || isempty(xIncon)
        unmatched = addMaskUnmatchedRows(unmatched, xBase, xCon, xIncon, ...
            true(size(xBase)), true(size(xCon)), true(size(xIncon)), ...
            tolerance, nearTolerance);
        return;
    end

    usedBase = false(numel(xBase), 1);
    usedCon = false(numel(xCon), 1);
    usedIncon = false(numel(xIncon), 1);

    [zeroTriples, zeroMethods, usedBase, usedCon, usedIncon, unmatched] = ...
        matchMaskZeroPoint(xBase, xCon, xIncon, usedBase, usedCon, ...
        usedIncon, unmatched, tolerance, nearTolerance);
    triples = [triples; zeroTriples];
    tripleMethods = [tripleMethods; zeroMethods];

    [exactTriples, exactMethods, usedBase, usedCon, usedIncon] = ...
        matchExactTriples(xBase, xCon, xIncon, usedBase, usedCon, ...
        usedIncon);
    triples = [triples; exactTriples];
    tripleMethods = [tripleMethods; exactMethods];

    remainingBase = find(~usedBase & ~isZeroContrast(xBase));
    remainingCon = find(~usedCon & ~isZeroContrast(xCon));
    remainingIncon = find(~usedIncon & ~isZeroContrast(xIncon));
    toleranceTriples = optimizeTripleMatches(xBase, xCon, xIncon, ...
        remainingBase, remainingCon, remainingIncon, tolerance);
    if ~isempty(toleranceTriples)
        usedBase(toleranceTriples(:, 1)) = true;
        usedCon(toleranceTriples(:, 2)) = true;
        usedIncon(toleranceTriples(:, 3)) = true;
        triples = [triples; toleranceTriples]; %#ok<AGROW>
        tripleMethods = [tripleMethods; ...
            repmat("nearest/tolerance", size(toleranceTriples, 1), 1)]; %#ok<AGROW>
    end

    [triples, sortIdx] = sortTriplesByMeanContrast( ...
        triples, xBase, xCon, xIncon);
    tripleMethods = tripleMethods(sortIdx);
    unmatched = addMaskUnmatchedRows(unmatched, xBase, xCon, xIncon, ...
        ~usedBase, ~usedCon, ~usedIncon, tolerance, nearTolerance);
end

function [zeroTriples, zeroMethods, usedBase, usedCon, usedIncon, ...
        unmatched] = matchMaskZeroPoint(xBase, xCon, xIncon, usedBase, ...
        usedCon, usedIncon, unmatched, tolerance, nearTolerance)
    zeroTriples = zeros(0, 3);
    zeroMethods = strings(0, 1);
    baseZero = find(isZeroContrast(xBase));
    conZero = find(isZeroContrast(xCon));
    inconZero = find(isZeroContrast(xIncon));
    nZero = min([numel(baseZero), numel(conZero), numel(inconZero)]);
    if nZero > 0
        zeroTriples = [baseZero(1:nZero), conZero(1:nZero), ...
            inconZero(1:nZero)];
        zeroMethods = repmat("exact/zero", nZero, 1);
        usedBase(baseZero(1:nZero)) = true;
        usedCon(conZero(1:nZero)) = true;
        usedIncon(inconZero(1:nZero)) = true;
    end

    maxZero = max([numel(baseZero), numel(conZero), numel(inconZero)]);
    for zeroIdx = (nZero + 1):maxZero
        baseIdx = getExtraIndex(baseZero, zeroIdx);
        conIdx = getExtraIndex(conZero, zeroIdx);
        inconIdx = getExtraIndex(inconZero, zeroIdx);
        if isnan(baseIdx)
            reason = "baseline zero contrast unavailable";
        elseif isnan(conIdx)
            reason = "con zero contrast unavailable";
        else
            reason = "incon zero contrast unavailable";
        end
        unmatched = appendUnmatched(unmatched, "deltaMask", "zero", ...
            baseIdx, conIdx, inconIdx, reason, 0, tolerance, ...
            nearTolerance, false);
        if ~isnan(baseIdx)
            usedBase(baseIdx) = true;
        end
        if ~isnan(conIdx)
            usedCon(conIdx) = true;
        end
        if ~isnan(inconIdx)
            usedIncon(inconIdx) = true;
        end
    end
end

function idx = getExtraIndex(indices, pos)
    if pos <= numel(indices)
        idx = indices(pos);
    else
        idx = NaN;
    end
end

function [triples, methods, usedBase, usedCon, usedIncon] = ...
        matchExactTriples(xBase, xCon, xIncon, usedBase, usedCon, ...
        usedIncon)
    triples = zeros(0, 3);
    methods = strings(0, 1);
    values = unique(xBase(~usedBase & ~isZeroContrast(xBase)), 'stable');
    for valueIdx = 1:numel(values)
        idxBase = find(~usedBase & xBase == values(valueIdx));
        idxCon = find(~usedCon & xCon == values(valueIdx));
        idxIncon = find(~usedIncon & xIncon == values(valueIdx));
        nMatch = min([numel(idxBase), numel(idxCon), numel(idxIncon)]);
        if nMatch == 0
            continue;
        end
        newTriples = [idxBase(1:nMatch), idxCon(1:nMatch), ...
            idxIncon(1:nMatch)];
        triples = [triples; newTriples]; %#ok<AGROW>
        methods = [methods; repmat("exact", nMatch, 1)]; %#ok<AGROW>
        usedBase(idxBase(1:nMatch)) = true;
        usedCon(idxCon(1:nMatch)) = true;
        usedIncon(idxIncon(1:nMatch)) = true;
    end
end

function triples = optimizeTripleMatches(xBase, xCon, xIncon, idxBase, ...
        idxCon, idxIncon, tolerance)
    triples = zeros(0, 3);
    if isempty(idxBase) || isempty(idxCon) || isempty(idxIncon)
        return;
    end

    bestTriples = zeros(0, 3);
    bestCount = -Inf;
    bestCost = Inf;
    idxBase = idxBase(:);
    idxCon = idxCon(:);
    idxIncon = idxIncon(:);
    usedCon = false(numel(xCon), 1);
    usedIncon = false(numel(xIncon), 1);

    recurse(1, usedCon, usedIncon, zeros(0, 3), 0);
    triples = bestTriples;

    function recurse(pos, usedConLocal, usedInconLocal, currentTriples, ...
            currentCost)
        remainingBase = numel(idxBase) - pos + 1;
        if size(currentTriples, 1) + remainingBase < bestCount
            return;
        end
        if pos > numel(idxBase)
            currentCount = size(currentTriples, 1);
            if currentCount > bestCount || ...
                    (currentCount == bestCount && currentCost < bestCost)
                bestCount = currentCount;
                bestCost = currentCost;
                bestTriples = currentTriples;
            end
            return;
        end

        baseIdx = idxBase(pos);
        candidates = zeros(0, 3);
        for conLocal = 1:numel(idxCon)
            conIdx = idxCon(conLocal);
            if usedConLocal(conIdx)
                continue;
            end
            for inconLocal = 1:numel(idxIncon)
                inconIdx = idxIncon(inconLocal);
                if usedInconLocal(inconIdx)
                    continue;
                end
                contrasts = [xBase(baseIdx), xCon(conIdx), ...
                    xIncon(inconIdx)];
                spread = max(contrasts) - min(contrasts);
                if spread <= tolerance
                    candidates(end + 1, :) = [conIdx, inconIdx, ...
                        spread]; %#ok<AGROW>
                end
            end
        end
        [~, order] = sort(candidates(:, 3));
        candidates = candidates(order, :);
        for candidateIdx = 1:size(candidates, 1)
            conIdx = candidates(candidateIdx, 1);
            inconIdx = candidates(candidateIdx, 2);
            nextUsedCon = usedConLocal;
            nextUsedIncon = usedInconLocal;
            nextUsedCon(conIdx) = true;
            nextUsedIncon(inconIdx) = true;
            recurse(pos + 1, nextUsedCon, nextUsedIncon, ...
                [currentTriples; baseIdx, conIdx, inconIdx], ...
                currentCost + candidates(candidateIdx, 3)); %#ok<AGROW>
        end
        recurse(pos + 1, usedConLocal, usedInconLocal, currentTriples, ...
            currentCost);
    end
end

function unmatched = addBiasUnmatchedRows(unmatched, xCon, xIncon, ...
        unmatchedCon, unmatchedIncon, tolerance, nearTolerance)
    for idx = find(unmatchedCon(:))'
        if isZeroContrast(xCon(idx))
            continue;
        end
        [nearestDiff, isNear] = nearestBiasDiscrepancy( ...
            xCon(idx), xIncon, tolerance, nearTolerance);
        unmatched = appendUnmatched(unmatched, "deltaBias", "con", NaN, ...
            idx, NaN, "no one-to-one match within 5 contrast points", ...
            nearestDiff, tolerance, nearTolerance, isNear);
    end
    for idx = find(unmatchedIncon(:))'
        if isZeroContrast(xIncon(idx))
            continue;
        end
        [nearestDiff, isNear] = nearestBiasDiscrepancy( ...
            xIncon(idx), xCon, tolerance, nearTolerance);
        unmatched = appendUnmatched(unmatched, "deltaBias", "incon", ...
            NaN, NaN, idx, ...
            "no one-to-one match within 5 contrast points", nearestDiff, ...
            tolerance, nearTolerance, isNear);
    end
end

function unmatched = addMaskUnmatchedRows(unmatched, xBase, xCon, xIncon, ...
        unmatchedBase, unmatchedCon, unmatchedIncon, tolerance, nearTolerance)
    for idx = find(unmatchedBase(:))'
        if isZeroContrast(xBase(idx))
            continue;
        end
        [nearestDiff, isNear] = nearestMaskSpread( ...
            xBase(idx), xCon, xIncon, tolerance, nearTolerance);
        unmatched = appendUnmatched(unmatched, "deltaMask", "baseline", ...
            idx, NaN, NaN, ...
            "no one-to-one match within 5 contrast points", nearestDiff, ...
            tolerance, nearTolerance, isNear);
    end
    for idx = find(unmatchedCon(:))'
        if isZeroContrast(xCon(idx))
            continue;
        end
        [nearestDiff, isNear] = nearestMaskSpread( ...
            xCon(idx), xBase, xIncon, tolerance, nearTolerance);
        unmatched = appendUnmatched(unmatched, "deltaMask", "con", NaN, ...
            idx, NaN, "no one-to-one match within 5 contrast points", ...
            nearestDiff, tolerance, nearTolerance, isNear);
    end
    for idx = find(unmatchedIncon(:))'
        if isZeroContrast(xIncon(idx))
            continue;
        end
        [nearestDiff, isNear] = nearestMaskSpread( ...
            xIncon(idx), xBase, xCon, tolerance, nearTolerance);
        unmatched = appendUnmatched(unmatched, "deltaMask", "incon", ...
            NaN, NaN, idx, ...
            "no one-to-one match within 5 contrast points", nearestDiff, ...
            tolerance, nearTolerance, isNear);
    end
end

function [nearestDiff, isNear] = nearestBiasDiscrepancy( ...
        xValue, xOther, tolerance, nearTolerance)
    xOther = xOther(~isZeroContrast(xOther));
    if isempty(xOther)
        nearestDiff = NaN;
        isNear = false;
        return;
    end
    nearestDiff = min(abs(xOther - xValue));
    isNear = nearestDiff > tolerance && nearestDiff <= nearTolerance;
end

function [nearestDiff, isNear] = nearestMaskSpread( ...
        xValue, xOtherA, xOtherB, tolerance, nearTolerance)
    xOtherA = xOtherA(~isZeroContrast(xOtherA));
    xOtherB = xOtherB(~isZeroContrast(xOtherB));
    if isempty(xOtherA) || isempty(xOtherB)
        nearestDiff = NaN;
        isNear = false;
        return;
    end
    nearestDiff = Inf;
    for idxA = 1:numel(xOtherA)
        for idxB = 1:numel(xOtherB)
            contrasts = [xValue, xOtherA(idxA), xOtherB(idxB)];
            nearestDiff = min(nearestDiff, max(contrasts) - min(contrasts));
        end
    end
    isNear = nearestDiff > tolerance && nearestDiff <= nearTolerance;
end

function unmatched = emptyUnmatchedStruct()
    unmatched = struct( ...
        'deltaMetric', strings(0, 1), ...
        'sourceRole', strings(0, 1), ...
        'baseIdx', nan(0, 1), ...
        'conIdx', nan(0, 1), ...
        'inconIdx', nan(0, 1), ...
        'reason', strings(0, 1), ...
        'nearestDiff', nan(0, 1), ...
        'nearMatch', false(0, 1));
end

function unmatched = appendUnmatched(unmatched, deltaMetric, sourceRole, ...
        baseIdx, conIdx, inconIdx, reason, nearestDiff, tolerance, ...
        nearTolerance, nearMatch)
    if nargin < 11
        nearMatch = nearestDiff > tolerance && nearestDiff <= nearTolerance;
    end
    unmatched.deltaMetric(end + 1, 1) = string(deltaMetric);
    unmatched.sourceRole(end + 1, 1) = string(sourceRole);
    unmatched.baseIdx(end + 1, 1) = baseIdx;
    unmatched.conIdx(end + 1, 1) = conIdx;
    unmatched.inconIdx(end + 1, 1) = inconIdx;
    unmatched.reason(end + 1, 1) = string(reason);
    unmatched.nearestDiff(end + 1, 1) = nearestDiff;
    unmatched.nearMatch(end + 1, 1) = nearMatch;
end

function audit = buildDeltaAuditRows(sideName, xBase, xCon, xIncon, ...
        biasPairs, xBias, yBias, biasMethods, biasUnmatched, ...
        maskTriples, xMask, yMask, maskMethods, maskUnmatched, ...
        tolerance, nearTolerance)
    nBase = numel(xBase);
    nCon = numel(xCon);
    nIncon = numel(xIncon);
    expectedBias = numel(yBias);
    expectedMask = numel(yMask);
    audit = emptyDeltaPointAuditTable();

    for idx = 1:numel(yBias)
        conIdx = biasPairs(idx, 1);
        inconIdx = biasPairs(idx, 2);
        contrasts = [xCon(conIdx), xIncon(inconIdx)];
        audit = [audit; makeDeltaAuditRow(sideName, "deltaBias", ...
            NaN, xCon(conIdx), xIncon(inconIdx), xBias(idx), yBias(idx), ...
            max(contrasts) - min(contrasts), tolerance, nearTolerance, ...
            true, biasMethods(idx), "", nBase, nCon, nIncon, ...
            expectedBias, numel(yBias), false, NaN)]; %#ok<AGROW>
    end
    for idx = 1:numel(biasUnmatched.deltaMetric)
        audit = [audit; makeUnmatchedAuditRow(sideName, xBase, xCon, ...
            xIncon, biasUnmatched, idx, tolerance, nearTolerance, nBase, ...
            nCon, nIncon, expectedBias, numel(yBias))]; %#ok<AGROW>
    end

    for idx = 1:numel(yMask)
        baseIdx = maskTriples(idx, 1);
        conIdx = maskTriples(idx, 2);
        inconIdx = maskTriples(idx, 3);
        contrasts = [xBase(baseIdx), xCon(conIdx), xIncon(inconIdx)];
        audit = [audit; makeDeltaAuditRow(sideName, "deltaMask", ...
            xBase(baseIdx), xCon(conIdx), xIncon(inconIdx), xMask(idx), ...
            yMask(idx), max(contrasts) - min(contrasts), tolerance, ...
            nearTolerance, true, maskMethods(idx), "", nBase, nCon, ...
            nIncon, expectedMask, numel(yMask), false, NaN)]; %#ok<AGROW>
    end
    for idx = 1:numel(maskUnmatched.deltaMetric)
        audit = [audit; makeUnmatchedAuditRow(sideName, xBase, xCon, ...
            xIncon, maskUnmatched, idx, tolerance, nearTolerance, nBase, ...
            nCon, nIncon, expectedMask, numel(yMask))]; %#ok<AGROW>
    end
end

function row = makeUnmatchedAuditRow(sideName, xBase, xCon, xIncon, ...
        unmatched, idx, tolerance, nearTolerance, nBase, nCon, nIncon, ...
        expectedCount, plottedCount)
    baseContrast = NaN;
    conContrast = NaN;
    inconContrast = NaN;
    if ~isnan(unmatched.baseIdx(idx))
        baseContrast = xBase(unmatched.baseIdx(idx));
    end
    if ~isnan(unmatched.conIdx(idx))
        conContrast = xCon(unmatched.conIdx(idx));
    end
    if ~isnan(unmatched.inconIdx(idx))
        inconContrast = xIncon(unmatched.inconIdx(idx));
    end
    row = makeDeltaAuditRow(sideName, unmatched.deltaMetric(idx), ...
        baseContrast, conContrast, inconContrast, NaN, NaN, ...
        unmatched.nearestDiff(idx), tolerance, nearTolerance, false, ...
        "unmatched", unmatched.reason(idx), nBase, nCon, nIncon, ...
        expectedCount, plottedCount, unmatched.nearMatch(idx), ...
        unmatched.nearestDiff(idx));
end

function row = makeDeltaAuditRow(sideName, deltaMetric, baselineContrast, ...
        conContrast, inconContrast, xDelta, deltaValue, maxAbsDiff, ...
        tolerance, nearTolerance, isPlotted, pairingMethod, ...
        reasonExcluded, nBase, nCon, nIncon, expectedCount, plottedCount, ...
        nearMatchOutsideDefaultTolerance, nearMatchMaxAbsDiff)
    schema = emptyDeltaPointAuditTable();
    row = table( ...
        string(""), nan, nan, string(""), string(sideName), ...
        string(deltaMetric), baselineContrast, conContrast, inconContrast, ...
        xDelta, deltaValue, maxAbsDiff, tolerance, nearTolerance, ...
        isPlotted, string(pairingMethod), string(reasonExcluded), ...
        nBase, nCon, nIncon, expectedCount, plottedCount, ...
        nearMatchOutsideDefaultTolerance, nearMatchMaxAbsDiff, ...
        'VariableNames', schema.Properties.VariableNames);
end

function [pairs, sortIdx] = sortRowsByMeanContrast(pairs, xA, xB)
    if isempty(pairs)
        sortIdx = zeros(0, 1);
        return;
    end
    xMean = mean([xA(pairs(:, 1)), xB(pairs(:, 2))], 2, 'omitnan');
    [~, sortIdx] = sort(xMean);
    pairs = pairs(sortIdx, :);
end

function [triples, sortIdx] = sortTriplesByMeanContrast( ...
        triples, xBase, xCon, xIncon)
    if isempty(triples)
        sortIdx = zeros(0, 1);
        return;
    end
    xMean = mean([xBase(triples(:, 1)), xCon(triples(:, 2)), ...
        xIncon(triples(:, 3))], 2, 'omitnan');
    [~, sortIdx] = sort(xMean);
    triples = triples(sortIdx, :);
end

function tf = isZeroContrast(x)
    tf = abs(x) <= 1e-9;
end

function [xOut, yOut] = collapseConditionDataLocal(xIn, yIn)
    xIn = xIn(:)';
    yIn = yIn(:)';
    validIdx = ~isnan(xIn) & ~isnan(yIn);
    xIn = xIn(validIdx);
    yIn = yIn(validIdx);
    [xIn, sortIdx] = sort(abs(xIn));
    yIn = yIn(sortIdx);
    if isempty(xIn)
        xOut = [];
        yOut = [];
        return;
    end

    [xOut, ~, groupIdx] = unique(xIn, 'stable');
    yOut = nan(size(xOut));
    for idx = 1:numel(xOut)
        yOut(idx) = mean(yIn(groupIdx == idx), 'omitnan');
    end
end
