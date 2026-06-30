function [labelsByBlock, sourceField] = recomputeMetaTablePowerClusters(bitmapData, mdlStruct, nBlocks)
% Recompute the same ordered power-effect clustering used by psycluster.

labelsByBlock = nan(nBlocks, 1);
sourceField = 'unavailable';
if exist('clusterOrderedPowerEffect', 'file') ~= 2 || ...
        ~isfield(bitmapData, 'meanPowerDensityWithinROI_mWmm2')
    return
end

[blockIDs, deltaBias, deltaMask, modelField] = getModelRows(mdlStruct, nBlocks);
if isempty(blockIDs)
    return
end

power = squeeze(bitmapData.meanPowerDensityWithinROI_mWmm2);
if isvector(power)
    power = power(:);
else
    power = mean(power, 1, 'omitnan')';
end
if numel(power) ~= nBlocks
    return
end

opts = struct( ...
    'minSessionsPerCluster', 3, ...
    'effectWeight', 1, ...
    'monotonicPenalty', 2, ...
    'kCandidates', [2 3], ...
    'minRelativeImprovementFor3Clusters', 0.10);
labelsByBlock(blockIDs) = clusterOrderedPowerEffect( ...
    power(blockIDs), deltaBias, deltaMask, opts);
sourceField = ['recomputed from ' modelField];
end


function [blockIDs, deltaBias, deltaMask, sourceField] = getModelRows(mdlStruct, nBlocks)
blockIDs = [];
deltaBias = [];
deltaMask = [];
sourceField = '';
if ~isstruct(mdlStruct) || isempty(fieldnames(mdlStruct))
    return
end

fields = fieldnames(mdlStruct);
preferred = fields(contains(fields, 'weibullfreeAllC1'));
fields = [preferred; setdiff(fields, preferred, 'stable')];
for ii = 1:numel(fields)
    mdl = mdlStruct.(fields{ii});
    if ~isstruct(mdl) || ~isfield(mdl, 'clusterBlocksIdx') || ...
            ~isfield(mdl, 'deltaBias')
        continue
    end

    rawBlocks = mdl.clusterBlocksIdx(:);
    rawBias = mdl.deltaBias(:);
    nValues = min(numel(rawBlocks), numel(rawBias));
    rawBlocks = rawBlocks(1:nValues);
    rawBias = rawBias(1:nValues);
    valid = isfinite(rawBlocks) & rawBlocks >= 1 & rawBlocks <= nBlocks & ...
        isfinite(rawBias);
    if ~any(valid)
        continue
    end

    rawMask = nan(nValues, 1);
    if isfield(mdl, 'deltaMask')
        mask = mdl.deltaMask(:);
        nMask = min(nValues, numel(mask));
        rawMask(1:nMask) = mask(1:nMask);
    end

    blockIDs = round(rawBlocks(valid));
    deltaBias = rawBias(valid);
    deltaMask = rawMask(valid);
    sourceField = fields{ii};
    return
end
end
