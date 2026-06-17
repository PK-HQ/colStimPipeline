function adjustedP = holmBonferroni(pValues)
% Return Holm-Bonferroni adjusted p-values, preserving NaN entries.

    adjustedP = nan(size(pValues));
    finiteIdx = find(isfinite(pValues));
    if isempty(finiteIdx)
        return;
    end

    finiteP = pValues(finiteIdx);
    [sortedP, sortOrder] = sort(finiteP);
    nTests = numel(sortedP);
    adjustedSorted = nan(size(sortedP));
    runningMaximum = 0;
    for rankIdx = 1:nTests
        candidate = (nTests - rankIdx + 1) .* sortedP(rankIdx);
        runningMaximum = max(runningMaximum, candidate);
        adjustedSorted(rankIdx) = min(1, runningMaximum);
    end

    unsorted = nan(size(finiteP));
    unsorted(sortOrder) = adjustedSorted;
    adjustedP(finiteIdx) = unsorted;
end
