function significantByBlock = testPowerClusterSignificance(clusterByBlock, deltaBiasByBlock, alpha)
% Mark clusters whose delta bias is significantly above zero.
% A cluster passes when either its mean (one-sample t-test) or median
% (signed-rank test) is significant using a one-sided right-tail test.

if nargin < 3 || isempty(alpha)
    alpha = 0.05;
end
significantByBlock = nan(size(clusterByBlock));
clusterIDs = unique(clusterByBlock(isfinite(clusterByBlock)))';

for clusterID = clusterIDs
    values = deltaBiasByBlock( ...
        clusterByBlock == clusterID & isfinite(deltaBiasByBlock));
    if isempty(values)
        continue
    end

    meanSignificant = false;
    medianSignificant = false;
    try
        meanSignificant = ttest(values, 0, ...
            'Tail', 'right', 'Alpha', alpha);
    catch ME
        warning('MetaTable:PowerClusterTTest', ...
            'Mean test failed for cluster %g: %s', clusterID, ME.message);
    end
    try
        [~, medianSignificant] = signrank(values, 0, ...
            'tail', 'right', 'alpha', alpha);
    catch ME
        warning('MetaTable:PowerClusterSignrank', ...
            'Median test failed for cluster %g: %s', clusterID, ME.message);
    end

    significantByBlock(clusterByBlock == clusterID) = ...
        double(logical(meanSignificant) || logical(medianSignificant));
end
end
