function [deltaBiasByBlock, sourceField] = getMetaTableDeltaBias(mdlStruct, nBlocks)
% Return per-block delta bias, preferring the weibullfreeAll C1 model.

deltaBiasByBlock = nan(nBlocks, 1);
sourceField = 'no psychometric model';
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

    blockIDs = mdl.clusterBlocksIdx(:);
    values = mdl.deltaBias(:);
    nValues = min(numel(blockIDs), numel(values));
    blockIDs = blockIDs(1:nValues);
    values = values(1:nValues);
    valid = isfinite(blockIDs) & blockIDs >= 1 & blockIDs <= nBlocks;
    deltaBiasByBlock(round(blockIDs(valid))) = values(valid);
    sourceField = fields{ii};

    if contains(fields{ii}, 'weibullfreeAllC1')
        return
    end
end
end
