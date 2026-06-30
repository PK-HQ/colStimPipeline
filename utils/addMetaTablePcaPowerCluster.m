function MetaTable = addMetaTablePcaPowerCluster(MetaTable, bitmapData, mdlStruct, datastruct, analysisBlockID, columnsDesired)
% Add PCA-denoised activity and ordered power-cluster metadata.

if nargin < 3 || isempty(mdlStruct), mdlStruct = struct(); end
if nargin < 4 || isempty(datastruct), datastruct = struct(); end
if nargin < 5 || isempty(analysisBlockID), analysisBlockID = 1:height(MetaTable); end
if nargin < 6, columnsDesired = []; end

nBlocks = size(bitmapData.nColumns, 2);
if ismember('blockID', MetaTable.Properties.VariableNames)
    blockIDs = round(MetaTable.blockID(:));
else
    blockIDs = (1:height(MetaTable))';
end
if any(blockIDs < 1 | blockIDs > nBlocks)
    error('MetaTable blockID is outside the %d-block bitmap dataset.', nBlocks);
end

pcaByRow = cell(height(MetaTable), 1);
if isfield(bitmapData, 'pcadenoisedresp')
    pcaSize = size(bitmapData.pcadenoisedresp);
    blockDim = find(pcaSize == nBlocks, 1, 'last');
    if isempty(blockDim)
        warning('MetaTable:PCA', 'Could not identify the block dimension of pcadenoisedresp.');
    else
        for row = 1:height(MetaTable)
            idx = repmat({':'}, 1, ndims(bitmapData.pcadenoisedresp));
            idx{blockDim} = blockIDs(row);
            pcaByRow{row} = squeeze(bitmapData.pcadenoisedresp(idx{:}));
        end
    end
end
MetaTable.opto_PCAdenoisedResp = pcaByRow;
if ismember('opto_bitmapCamSpace', MetaTable.Properties.VariableNames)
    MetaTable = movevars(MetaTable, 'opto_PCAdenoisedResp', 'Before', 'opto_bitmapCamSpace');
end

[clusterByBlock, deltaBiasByBlock, source] = resolveMetaTablePowerClusters( ...
    bitmapData, mdlStruct, nBlocks);
significantByBlock = testPowerClusterSignificance( ...
    clusterByBlock, deltaBiasByBlock, 0.05);

MetaTable.psy_powerClusterID = clusterByBlock(blockIDs);
MetaTable.psy_powerClusterSignificant = significantByBlock(blockIDs);
MetaTable = movevars(MetaTable, ...
    {'psy_powerClusterID', 'psy_powerClusterSignificant'}, 'After', ...
    MetaTable.Properties.VariableNames{end - 2});

nameMap = {
    'opto_PCAdenoisedResp', 'opto.PCAdenoisedResp'
    'psy_powerClusterID', 'psy.powerClusterID'
    'psy_powerClusterSignificant', 'psy.powerClusterSignificant'};
for ii = 1:size(nameMap, 1)
    tf = strcmp(MetaTable.Properties.VariableNames, nameMap{ii, 1});
    MetaTable.Properties.VariableDescriptions{tf} = nameMap{ii, 2};
end

fprintf(['MetaTable: %d/%d rows have power-cluster labels from %s; ' ...
    '%d rows are in a cluster significantly above zero.\n'], ...
    sum(isfinite(MetaTable.psy_powerClusterID)), height(MetaTable), source, ...
    sum(MetaTable.psy_powerClusterSignificant == 1));

saveAugmentedMetaTable(MetaTable, datastruct, analysisBlockID, columnsDesired);
end
