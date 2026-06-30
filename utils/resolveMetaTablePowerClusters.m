function [labelsByBlock, deltaBiasByBlock, sourceField] = resolveMetaTablePowerClusters(bitmapData, mdlStruct, nBlocks)
% Resolve ordered power clusters for MetaTable rows.

[deltaBiasByBlock, modelField] = getMetaTableDeltaBias(mdlStruct, nBlocks);
[labelsByBlock, sourceField] = getSavedMetaTablePowerClusters(mdlStruct, nBlocks);

haveSavedLabels = any(isfinite(labelsByBlock));
if ~haveSavedLabels
    [labelsByBlock, sourceField] = recomputeMetaTablePowerClusters( ...
        bitmapData, mdlStruct, nBlocks);
end

if strcmp(sourceField, 'unavailable')
    sourceField = modelField;
end
end
