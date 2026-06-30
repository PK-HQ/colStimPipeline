function [labelsByBlock, sourceField] = getSavedMetaTablePowerClusters(mdlStruct, nBlocks)
% Read ordered power-cluster labels saved by psycluster.

labelsByBlock = nan(nBlocks, 1);
sourceField = 'unavailable';
if ~isstruct(mdlStruct) || isempty(fieldnames(mdlStruct))
    return
end
fields = fieldnames(mdlStruct);

% Current output: already aligned to bitmap block index.
for ii = 1:numel(fields)
    name = fields{ii};
    if endsWith(name, 'PowerClusterAggregateLabelsByBlock')
        values = mdlStruct.(name);
        if isnumeric(values) && numel(values) == nBlocks
            labelsByBlock = values(:);
            sourceField = name;
            return
        end
    end
end

% Intermediate output: cluster summaries with explicit block indices.
for ii = 1:numel(fields)
    name = fields{ii};
    if ~endsWith(name, 'PowerClusterAggregateSummary')
        continue
    end
    summary = mdlStruct.(name);
    if ~isstruct(summary)
        continue
    end
    for jj = 1:numel(summary)
        if isfield(summary(jj), 'clusterID') && isfield(summary(jj), 'blockIndices')
            blocks = cleanBlocks(summary(jj).blockIndices, nBlocks);
            labelsByBlock(blocks) = summary(jj).clusterID;
        end
    end
    if any(isfinite(labelsByBlock))
        sourceField = name;
        return
    end
end

% Older aggregate output: sourceBlockIndices within each cluster record.
for ii = 1:numel(fields)
    name = fields{ii};
    if ~endsWith(name, 'PowerClusterAggregate')
        continue
    end
    aggregate = mdlStruct.(name);
    if ~isstruct(aggregate)
        continue
    end
    for jj = 1:numel(aggregate)
        if isfield(aggregate(jj), 'clusterID') && ...
                isfield(aggregate(jj), 'sourceBlockIndices')
            blocks = cleanBlocks(aggregate(jj).sourceBlockIndices, nBlocks);
            labelsByBlock(blocks) = aggregate(jj).clusterID;
        end
    end
    if any(isfinite(labelsByBlock))
        sourceField = name;
        return
    end
end
end


function blocks = cleanBlocks(blocks, nBlocks)
blocks = blocks(:);
blocks = round(blocks(isfinite(blocks) & blocks >= 1 & blocks <= nBlocks));
end
