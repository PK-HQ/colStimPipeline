function MetaTable = buildBlockMetadata(behavioralData, bitmapData, columnsDesired, columnsSpread, blockData, datastruct, analysisBlockID, mdlStruct)

if nargin < 5 || isempty(blockData)
    blockData = struct();
end
if nargin < 6 || isempty(datastruct)
    datastruct = struct();
end
if nargin < 7 || isempty(analysisBlockID)
    analysisBlockID = 1:size(bitmapData.nColumns, 2);
end
if nargin < 8 || isempty(mdlStruct)
    mdlStruct = struct();
end

meanCol = mean(bitmapData.nColumns, 1, 'omitnan');
selectedBlocks = 1:numel(meanCol);
nBlocks = numel(selectedBlocks);
nTotalBlocks = numel(meanCol);
[deltaMask, deltaBias] = get_psychometric_deltas(mdlStruct, nTotalBlocks);

rowTemplate = make_empty_metadata_row();
rows = repmat(rowTemplate, nBlocks, 1);

for ii = 1:nBlocks
    blockID = selectedBlocks(ii);

    row = rowTemplate;

    row.blockID = blockID;
    row.meanColumns = meanCol(blockID);

    row = add_experiment_metadata(row, blockData, datastruct, analysisBlockID, blockID, nTotalBlocks);

    % visual / behavioral metadata
    [gaborContrastRows, gaborContrastAll] = get_gabor_contrasts(behavioralData, blockID);
    row.GaborContrast_pc = {gaborContrastAll};
    row.vis_contrast_baseline = {gaborContrastRows{1}};
    row.vis_contrast_conopto = {gaborContrastRows{2}};
    row.vis_contrast_inconopto = {gaborContrastRows{3}};
    row.GaborSize_deg    = get_visual_stim_field(behavioralData, blockID, 'gaborSize', 'Stimulus', 'GaborSize__deg');
    row.GaborSF_cpd      = get_visual_stim_field(behavioralData, blockID, 'gaborSF', 'Stimulus', 'GaborSF_cpd');
    row.GaborOrt_deg     = get_visual_stim_field(behavioralData, blockID, 'gaborOrt', 'Stimulus', 'GaborOrt__deg');
    row.GaborPhs_deg     = get_visual_stim_field(behavioralData, blockID, 'gaborPhs', 'Stimulus', 'GaborPhs__deg');
    [row.GaborX_deg, row.GaborY_deg] = get_gabor_position(behavioralData, blockID);
    row.vis_contrast = {gaborContrastAll};
    row.vis_ort = {row.GaborOrt_deg};
    row.vis_sz = {row.GaborSize_deg};
    row.vis_sf = {row.GaborSF_cpd};
    row.vis_phs = {row.GaborPhs_deg};
    row.vis_pos = {[row.GaborX_deg row.GaborY_deg]};

    % bitmap metadata summaries
    row.gridSize_mean = mean(bitmapData.gridSize(:,blockID), 'omitnan');
    row.nColumns_Mean = mean(bitmapData.nColumns(:,blockID), 'omitnan');
    row.nColumns_min  = min(bitmapData.nColumns(:,blockID), [], 'omitnan');
    row.nColumns_max  = max(bitmapData.nColumns(:,blockID), [], 'omitnan');

    row.sensitivity_mean = mean(bitmapData.sensitivity(:,blockID), 'omitnan');
    row.adaptthresh_mean = mean(bitmapData.adaptthresh(:,blockID), 'omitnan');
    row.pixelsON_mean    = mean(bitmapData.pixelsON(:,blockID), 'omitnan');

    row.meanPowerDensityWithinROI_mWmm2_mean = mean(bitmapData.meanPowerDensityWithinROI_mWmm2(:,:,blockID), 'all', 'omitnan');
    row.meanPowerDensityWithinROI_mWmm2_max  = max(bitmapData.meanPowerDensityWithinROI_mWmm2(:,:,blockID), [], 'all');

    row.totalPowerToOnPixelsWithinROI_mW_mean = mean(bitmapData.totalPowerToOnPixelsWithinROI_mW(:,:,blockID), 'all', 'omitnan');
    row.totalPowerToOnPixelsWithinROI_mW_max  = max(bitmapData.totalPowerToOnPixelsWithinROI_mW(:,:,blockID), [], 'all');

    row.projectorPowerDensity_mWmm2_mean = mean(bitmapData.projectorPowerDensity_mWmm2(:,:,blockID), 'all', 'omitnan');
    row.projectorPowerDensity_mWmm2_max  = max(bitmapData.projectorPowerDensity_mWmm2(:,:,blockID), [], 'all');

    row.temporalDutyCycle_mean = mean(bitmapData.temporalDutyCycle(:,:,blockID), 'all', 'omitnan');
    row = add_bitmap_metadata(row, bitmapData, blockID, nTotalBlocks);
    row.psy_deltaMask = deltaMask(blockID);
    row.psy_deltaBias = deltaBias(blockID);

    rows(ii) = row;
end

MetaTable = struct2table(rows);
MetaTable.Properties.VariableDescriptions = get_export_header_names(MetaTable.Properties.VariableNames);
MetaTable = filter_metatable_by_columns(MetaTable, columnsDesired, columnsSpread);
MetaTable = remove_metatable_output_columns(MetaTable);
MetaTable = move_psychometric_columns_to_end(MetaTable);
save_metatable_outputs(MetaTable, datastruct, analysisBlockID, columnsDesired);

end


function row = make_empty_metadata_row()

row = struct();
row.animal_ID = {[]};
row.animal_hemisphere = {[]};
row.date_camprojcalib = {[]};
row.run_camprojcalib = {[]};
row.date_opto = {[]};
row.run_opto = {[]};
row.date_baseline = {[]};
row.run_baseline = {[]};
row.date_ortmap = {[]};
row.run_ortmap = {[]};
row.date_visfootprint = {[]};
row.run_visfootprint = {[]};
row.blockID = NaN;
row.meanColumns = NaN;
row.GaborContrast_pc = {NaN};
row.vis_contrast_baseline = {NaN};
row.vis_contrast_conopto = {NaN};
row.vis_contrast_inconopto = {NaN};
row.GaborSize_deg = NaN;
row.GaborSF_cpd = NaN;
row.GaborOrt_deg = NaN;
row.GaborPhs_deg = NaN;
row.GaborX_deg = NaN;
row.GaborY_deg = NaN;
row.vis_contrast = {NaN};
row.vis_ort = {NaN};
row.vis_sz = {NaN};
row.vis_sf = {NaN};
row.vis_phs = {NaN};
row.vis_pos = {NaN};
row.opto_480LED = {[]};
row.opto_480ND = {[]};
row.opto_580LED = {[]};
row.opto_580ND = {[]};
row.opto_gridsize = {[]};
row.opto_gamma = {[]};
row.opto_threshSens = {[]};
row.opto_threshAdapt = {[]};
row.opto_ort = {[]};
row.opto_gausscond = {[]};
row.opto_gausslevel = {[]};
row.opto_gaussmax = {[]};
row.opto_transformParams = {[]};
row.opto_bitmapCamSpace = {[]};
row.opto_bitmapProjSpace = {[]};
row.opto_nColumns = {[]};
row.opto_columnArea = {[]};
row.opto_pixelsON = {[]};
row.opto_temporalDutyCycle = {[]};
row.opto_meanPowerDensityWithinROI_mWmm2 = {[]};
row.opto_totalPowerToOnPixelsWithinROI_mW = {[]};
row.opto_projectorPowerDensity_mWmm2 = {[]};
row.gridSize_mean = NaN;
row.nColumns_Mean = NaN;
row.nColumns_min = NaN;
row.nColumns_max = NaN;
row.sensitivity_mean = NaN;
row.adaptthresh_mean = NaN;
row.pixelsON_mean = NaN;
row.meanPowerDensityWithinROI_mWmm2_mean = NaN;
row.meanPowerDensityWithinROI_mWmm2_max = NaN;
row.totalPowerToOnPixelsWithinROI_mW_mean = NaN;
row.totalPowerToOnPixelsWithinROI_mW_max = NaN;
row.projectorPowerDensity_mWmm2_mean = NaN;
row.projectorPowerDensity_mWmm2_max = NaN;
row.temporalDutyCycle_mean = NaN;
row.psy_deltaMask = NaN;
row.psy_deltaBias = NaN;
row.bitmapFile = '';

end


function row = add_experiment_metadata(row, blockData, datastruct, analysisBlockID, blockID, nTotalBlocks)

dsEntry = get_datastruct_entry(datastruct, analysisBlockID, blockID);

row.animal_ID = {get_struct_field(dsEntry, {'monkeyNo', 'monkeyID'})};
row.animal_hemisphere = {get_struct_field(dsEntry, {'chamber', 'hemisphere'})};

row.date_camprojcalib = {get_block_struct_field(blockData, 'alignment_date', blockID, nTotalBlocks)};
row.run_camprojcalib  = {get_block_struct_field(blockData, 'alignment_run', blockID, nTotalBlocks)};

[optoDate, optoRun] = get_opto_date_run(dsEntry, blockData, blockID, nTotalBlocks);
row.date_opto = {optoDate};
row.run_opto  = {optoRun};

row.date_baseline = {get_block_struct_field(blockData, 'baseline_date', blockID, nTotalBlocks)};
row.run_baseline  = {get_block_struct_field(blockData, 'baseline_run', blockID, nTotalBlocks)};

row.date_ortmap = {get_block_struct_field(blockData, 'ortmap_date', blockID, nTotalBlocks)};
row.run_ortmap  = {get_block_struct_field(blockData, 'ortmap_run', blockID, nTotalBlocks)};

rawFootprintPath = get_block_struct_field(blockData, 'gaussfootprint_run', blockID, nTotalBlocks);
[fpDate, fpRun] = parse_vis_footprint_path(rawFootprintPath);
row.date_visfootprint = {fpDate};
row.run_visfootprint  = {fpRun};

end


function row = add_bitmap_metadata(row, bitmapData, blockID, nTotalBlocks)

row.opto_480LED = {get_bitmap_field(bitmapData, 'blueLED', blockID, nTotalBlocks)};
row.opto_480ND = {get_bitmap_field(bitmapData, 'blueND', blockID, nTotalBlocks)};
row.opto_580LED = {get_bitmap_field(bitmapData, 'orangeLED', blockID, nTotalBlocks)};
row.opto_580ND = {get_bitmap_field(bitmapData, 'orangeND', blockID, nTotalBlocks)};
row.opto_gridsize = {get_bitmap_field(bitmapData, 'gridSize', blockID, nTotalBlocks)};
row.opto_gamma = {get_bitmap_field(bitmapData, 'gammaCorrFactor', blockID, nTotalBlocks)};
row.opto_threshSens = {get_bitmap_field(bitmapData, 'sensitivity', blockID, nTotalBlocks)};
row.opto_threshAdapt = {get_bitmap_field(bitmapData, 'adaptthresh', blockID, nTotalBlocks)};
row.opto_ort = {get_bitmap_field(bitmapData, 'orts', blockID, nTotalBlocks)};
row.opto_gausscond = {get_bitmap_field(bitmapData, 'gaussianCond', blockID, nTotalBlocks)};
row.opto_gausslevel = {get_bitmap_field(bitmapData, 'gaussianContourLevel', blockID, nTotalBlocks)};
row.opto_gaussmax = {get_bitmap_field(bitmapData, 'gaussianContourLevelMax', blockID, nTotalBlocks)};
row.opto_transformParams = {get_bitmap_field(bitmapData, 'transformParams', blockID, nTotalBlocks)};
row.opto_bitmapCamSpace = {get_bitmap_field(bitmapData, 'columnarbitmapTFcamspace', blockID, nTotalBlocks)};
row.opto_bitmapProjSpace = {get_bitmap_field(bitmapData, 'columnarbitmapTFprojspace', blockID, nTotalBlocks)};
row.opto_nColumns = {get_bitmap_field(bitmapData, 'nColumns', blockID, nTotalBlocks)};
row.opto_columnArea = {combine_column_areas(get_bitmap_field(bitmapData, 'columnAreas', blockID, nTotalBlocks))};
row.opto_pixelsON = {get_bitmap_field(bitmapData, 'pixelsON', blockID, nTotalBlocks)};
row.opto_temporalDutyCycle = {get_bitmap_field(bitmapData, 'temporalDutyCycle', blockID, nTotalBlocks)};
row.opto_meanPowerDensityWithinROI_mWmm2 = {get_bitmap_field(bitmapData, 'meanPowerDensityWithinROI_mWmm2', blockID, nTotalBlocks)};
row.opto_totalPowerToOnPixelsWithinROI_mW = {get_bitmap_field(bitmapData, 'totalPowerToOnPixelsWithinROI_mW', blockID, nTotalBlocks)};
row.opto_projectorPowerDensity_mWmm2 = {get_bitmap_field(bitmapData, 'projectorPowerDensity_mWmm2', blockID, nTotalBlocks)};

end


function val = get_bitmap_field(bitmapData, fieldName, blockID, nTotalBlocks)

val = [];

if isfield(bitmapData, fieldName)
    val = get_block_slice(bitmapData.(fieldName), blockID, nTotalBlocks);
end

end


function val = get_block_struct_field(S, fieldName, blockID, nTotalBlocks)

val = [];

if isstruct(S) && isfield(S, fieldName)
    val = get_block_slice(S.(fieldName), blockID, nTotalBlocks);
end

end


function val = get_block_slice(allVal, blockID, nTotalBlocks)

sz = size(allVal);

if numel(allVal) == 1
    if iscell(allVal)
        val = allVal{1};
    else
        val = allVal;
    end
    return
end

if isstruct(allVal)
    if numel(allVal) == nTotalBlocks
        val = allVal(blockID);
        return
    end
    blockDim = find_block_dim(sz, nTotalBlocks);
    idx = repmat({':'}, 1, ndims(allVal));
    idx{blockDim} = blockID;
    val = squeeze(allVal(idx{:}));
    return
end

if iscell(allVal)
    if numel(allVal) == nTotalBlocks
        val = allVal{blockID};
        return
    end
    blockDim = find_block_dim(sz, nTotalBlocks);
    idx = repmat({':'}, 1, ndims(allVal));
    idx{blockDim} = blockID;
    tmp = squeeze(allVal(idx{:}));
    if iscell(tmp) && numel(tmp) == 1
        val = tmp{1};
    else
        val = tmp;
    end
    return
end

blockDim = find_block_dim(sz, nTotalBlocks);
idx = repmat({':'}, 1, ndims(allVal));
idx{blockDim} = blockID;
val = squeeze(allVal(idx{:}));

end


function blockDim = find_block_dim(sz, nTotalBlocks)

candidateDims = find(sz == nTotalBlocks);

if isempty(candidateDims)
    error('No dimension has size nBlocks = %d. Size was [%s].', ...
        nTotalBlocks, num2str(sz));
end

blockDim = candidateDims(end);

end


function dsEntry = get_datastruct_entry(datastruct, analysisBlockID, blockID)

dsEntry = struct();

if ~isstruct(datastruct) || isempty(datastruct) || numel(analysisBlockID) < blockID
    return
end

dsIdx = analysisBlockID(blockID);
if dsIdx >= 1 && numel(datastruct) >= dsIdx
    dsEntry = datastruct(dsIdx);
end

end


function val = get_struct_field(S, fieldNames)

val = [];

if ~isstruct(S)
    return
end

for fieldID = 1:numel(fieldNames)
    fieldName = fieldNames{fieldID};
    if isfield(S, fieldName)
        val = S.(fieldName);
        return
    end
end

end


function [optoDate, optoRun] = get_opto_date_run(dsEntry, blockData, blockID, nTotalBlocks)

optoDate = get_struct_field(dsEntry, {'date'});
optoRun = get_struct_field(dsEntry, {'run', 'runNum', 'run_number'});

if isempty(optoDate)
    optoDate = get_block_struct_field(blockData, 'opto_date', blockID, nTotalBlocks);
end

if isempty(optoRun)
    optoRun = get_block_struct_field(blockData, 'opto_run', blockID, nTotalBlocks);
end

end


function [fpDate, fpRun] = parse_vis_footprint_path(rawPath)

fpDate = [];
fpRun = [];

if isempty(rawPath)
    return
end

if isnumeric(rawPath) && all(isnan(rawPath(:)))
    return
end

if iscell(rawPath)
    if isempty(rawPath)
        return
    end
    rawPath = rawPath{1};
    if isnumeric(rawPath) && all(isnan(rawPath(:)))
        return
    end
end

rawPath = char(string(rawPath));

dateTok = regexp(rawPath, '[A-Za-z]+(\d{8})', 'tokens', 'once');
if ~isempty(dateTok)
    fpDate = str2double(dateTok{1});
end

runTok = regexp(rawPath, '[\\/]+run(\d+)[\\/]+', 'tokens', 'once');
if isempty(runTok)
    runTok = regexp(rawPath, 'run(\d+)', 'tokens', 'once');
end

if ~isempty(runTok)
    fpRun = str2double(runTok{1});
end

end


function out = combine_column_areas(val)

if isempty(val) || ~iscell(val)
    out = val;
    return
end

vals = val(:);
nRows = numel(vals);
rowVecs = cell(nRows, 1);
rowLens = zeros(nRows, 1);

for rowID = 1:nRows
    x = vals{rowID};
    if iscell(x) && numel(x) == 1
        x = x{1};
    end
    if isempty(x)
        rowVecs{rowID} = [];
        continue
    end
    if ~isnumeric(x)
        out = val;
        return
    end
    rowVecs{rowID} = x(:)';
    rowLens(rowID) = numel(rowVecs{rowID});
end

maxLen = max(rowLens);
out = NaN(nRows, maxLen);

for rowID = 1:nRows
    if rowLens(rowID) > 0
        out(rowID, 1:rowLens(rowID)) = rowVecs{rowID};
    end
end

end


function [deltaMask, deltaBias] = get_psychometric_deltas(mdlStruct, nTotalBlocks)

deltaMask = NaN(nTotalBlocks, 1);
deltaBias = NaN(nTotalBlocks, 1);

if ~isstruct(mdlStruct) || isempty(fieldnames(mdlStruct))
    return
end

mdlFields = fieldnames(mdlStruct);
for fieldID = 1:numel(mdlFields)
    modelData = mdlStruct.(mdlFields{fieldID});
    if ~isstruct(modelData) || ~isfield(modelData, 'clusterBlocksIdx') || ...
            ~isfield(modelData, 'deltaMask') || ~isfield(modelData, 'deltaBias')
        continue
    end

    blockIDs = modelData.clusterBlocksIdx(:);
    thisDeltaMask = modelData.deltaMask(:);
    thisDeltaBias = modelData.deltaBias(:);
    nVals = min([numel(blockIDs), numel(thisDeltaMask), numel(thisDeltaBias)]);

    for valID = 1:nVals
        blockID = blockIDs(valID);
        if blockID >= 1 && blockID <= nTotalBlocks
            deltaMask(blockID) = thisDeltaMask(valID);
            deltaBias(blockID) = thisDeltaBias(valID);
        end
    end
end

end


function [contrastRows, contrastAll] = get_gabor_contrasts(behavioralData, blockID)

contrastRows = {NaN, NaN, NaN};
contrastAll = NaN;

if isfield(behavioralData, 'gaborContrasts') && ...
        ndims(behavioralData.gaborContrasts) >= 3 && ...
        size(behavioralData.gaborContrasts, 3) >= blockID
    blockContrasts = behavioralData.gaborContrasts(:,:,blockID);
    nRows = min(3, size(blockContrasts, 1));

    for rowID = 1:nRows
        conditionContrasts = blockContrasts(rowID,:);
        conditionContrasts = conditionContrasts(~isnan(conditionContrasts));
        if ~isempty(conditionContrasts)
            contrastRows{rowID} = conditionContrasts(:)';
        end
    end

    combinedContrasts = [contrastRows{:}];
    combinedContrasts = combinedContrasts(~isnan(combinedContrasts));
    if ~isempty(combinedContrasts)
        contrastAll = unique(combinedContrasts, 'stable');
        return
    end
end

if isfield(behavioralData, 'gaborContrast') && ...
        ndims(behavioralData.gaborContrast) >= 4 && ...
        size(behavioralData.gaborContrast, 4) >= blockID
    blockContrasts = behavioralData.gaborContrast(:,:,:,blockID);
    nRows = min(3, size(blockContrasts, 1));

    for rowID = 1:nRows
        conditionContrasts = squeeze(blockContrasts(rowID,:,:));
        conditionContrasts = conditionContrasts(:)';
        conditionContrasts = conditionContrasts(~isnan(conditionContrasts));
        if ~isempty(conditionContrasts)
            contrastRows{rowID} = unique(conditionContrasts, 'stable');
        end
    end

    combinedContrasts = [contrastRows{:}];
    combinedContrasts = combinedContrasts(~isnan(combinedContrasts));
    if ~isempty(combinedContrasts)
        contrastAll = unique(combinedContrasts, 'stable');
        return
    end
end

visualStimContrasts = get_visual_stim_contrasts(behavioralData, blockID);
if ~all(isnan(visualStimContrasts))
    contrastAll = visualStimContrasts;
    return
end

legacyContrasts = get_opto_ts_field(behavioralData, blockID, 'Stimulus', 'GaborContrast__pc');
if ~all(isnan(legacyContrasts))
    contrastAll = legacyContrasts;
end

end


function val = get_visual_stim_contrasts(behavioralData, blockID)

val = NaN;

if ~isfield(behavioralData, 'visualStim') || numel(behavioralData.visualStim) < blockID
    return
end

visualStim = behavioralData.visualStim(blockID);
fieldNames = {'gaborContrast', 'gaborCon'};
for fieldID = 1:numel(fieldNames)
    fieldName = fieldNames{fieldID};
    if isfield(visualStim, fieldName)
        val = format_metadata_value(visualStim.(fieldName));
        if ~all(isnan(val))
            return
        end
    end
end

end


function val = get_visual_stim_field(behavioralData, blockID, visualFieldName, optoGroupName, optoFieldName)

val = NaN;

if isfield(behavioralData, 'visualStim') && ...
        numel(behavioralData.visualStim) >= blockID && ...
        isfield(behavioralData.visualStim(blockID), visualFieldName)
    val = format_metadata_value(behavioralData.visualStim(blockID).(visualFieldName));
    if ~all(isnan(val))
        return
    end
end

val = get_opto_ts_field(behavioralData, blockID, optoGroupName, optoFieldName);

end


function [xDeg, yDeg] = get_gabor_position(behavioralData, blockID)

xDeg = NaN;
yDeg = NaN;

if isfield(behavioralData, 'visualStim') && ...
        numel(behavioralData.visualStim) >= blockID && ...
        isfield(behavioralData.visualStim(blockID), 'gaborPos')
    gaborPos = format_metadata_value(behavioralData.visualStim(blockID).gaborPos);
    if numel(gaborPos) >= 2
        xDeg = gaborPos(1);
        yDeg = gaborPos(2);
        return
    end
end

xDeg = get_opto_ts_field(behavioralData, blockID, 'Stimulus_Position', 'X__deg');
yDeg = get_opto_ts_field(behavioralData, blockID, 'Stimulus_Position', 'Y__deg');

end


function val = get_opto_ts_field(behavioralData, blockID, groupName, fieldName)

val = NaN;

if ~isfield(behavioralData, 'optoTS') || numel(behavioralData.optoTS) < blockID || ...
        ~isfield(behavioralData.optoTS(blockID), 'Header') || ...
        ~isfield(behavioralData.optoTS(blockID).Header, 'ConditionParams') || ...
        ~isfield(behavioralData.optoTS(blockID).Header.ConditionParams, groupName)
    return
end

S = behavioralData.optoTS(blockID).Header.ConditionParams.(groupName);

if isfield(S, fieldName)
    val = format_metadata_value(S.(fieldName));
end

end


function val = format_metadata_value(val)

if isempty(val)
    val = NaN;
elseif ~isscalar(val)
    val = val(:)';
end

end


function MetaTable = filter_metatable_by_columns(MetaTable, columnsDesired, columnsSpread)

if isempty(columnsDesired) || isempty(columnsSpread)
    return
end

lowerBound = columnsDesired - columnsSpread;
upperBound = columnsDesired + columnsSpread;
keepRows = MetaTable.nColumns_Mean >= lowerBound & MetaTable.nColumns_Mean <= upperBound;
MetaTable = MetaTable(keepRows, :);

end


function MetaTable = remove_metatable_output_columns(MetaTable)

columnsToRemove = {
    'GaborContrast_pc'
    'GaborSize_deg'
    'GaborSF_cpd'
    'GaborOrt_deg'
    'GaborPhs_deg'
    'GaborX_deg'
    'GaborY_deg'
    'gridSize_mean'
    'nColumns_Mean'
    'nColumns_min'
    'nColumns_max'
    'sensitivity_mean'
    'adaptthresh_mean'
    'pixelsON_mean'
    'meanPowerDensityWithinROI_mWmm2_mean'
    'meanPowerDensityWithinROI_mWmm2_max'
    'totalPowerToOnPixelsWithinROI_mW_mean'
    'totalPowerToOnPixelsWithinROI_mW_max'
    'projectorPowerDensity_mWmm2_mean'
    'projectorPowerDensity_mWmm2_max'
    'temporalDutyCycle_mean'
    'bitmapFile'
    };

existingColumns = intersect(columnsToRemove, MetaTable.Properties.VariableNames, 'stable');
if ~isempty(existingColumns)
    MetaTable(:, existingColumns) = [];
end

end


function MetaTable = move_psychometric_columns_to_end(MetaTable)

psyNames = {'psy_deltaMask', 'psy_deltaBias'};
psyNames = psyNames(ismember(psyNames, MetaTable.Properties.VariableNames));
if isempty(psyNames)
    return
end

remainingNames = setdiff(MetaTable.Properties.VariableNames, psyNames, 'stable');
if isempty(remainingNames)
    return
end

MetaTable = movevars(MetaTable, psyNames, 'After', remainingNames{end});

end


function save_metatable_outputs(MetaTable, datastruct, analysisBlockID, columnsDesired)

if ~isstruct(datastruct) || isempty(datastruct) || ...
        ~isfield(datastruct, 'monkey') || ~isfield(datastruct, 'chamber')
    return
end

saveEntry = datastruct(1);
if ~isempty(analysisBlockID) && analysisBlockID(1) >= 1 && numel(datastruct) >= analysisBlockID(1)
    saveEntry = datastruct(analysisBlockID(1));
end

monkeyName = saveEntry.monkey;
chamberStr = saveEntry.chamber;
columnsStr = num2str(columnsDesired);
outDir = ['Y:/' monkeyName '/Meta/summary'];

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

matOutPath = [outDir '/metaTable-' chamberStr columnsStr '.mat'];
xlsxOutPath = [outDir '/metaTable-' chamberStr columnsStr '.xlsx'];

save(matOutPath, 'MetaTable');

excelCell = table_to_excel_cell(MetaTable);
if exist(xlsxOutPath, 'file')
    delete(xlsxOutPath);
end
writecell(excelCell, xlsxOutPath, 'Sheet', 'metadata');

end


function excelCell = table_to_excel_cell(MetaTable)

nRows = height(MetaTable);
nCols = width(MetaTable);
headerNames = MetaTable.Properties.VariableDescriptions;

if numel(headerNames) ~= nCols || any(cellfun(@isempty, headerNames))
    headerNames = MetaTable.Properties.VariableNames;
end

excelCell = cell(nRows + 1, nCols);
excelCell(1, :) = headerNames;

for rowID = 1:nRows
    for colID = 1:nCols
        varName = MetaTable.Properties.VariableNames{colID};
        if iscell(MetaTable.(varName))
            rawVal = MetaTable.(varName){rowID};
        else
            rawVal = MetaTable.(varName)(rowID, :);
        end
        excelCell{rowID + 1, colID} = excel_scalarize(rawVal);
    end
end

end


function out = excel_scalarize(x)

if isempty(x)
    out = '';
    return
end

if iscell(x)
    if numel(x) == 1
        out = excel_scalarize(x{1});
        return
    end

    parts = cell(numel(x), 1);
    for cellID = 1:numel(x)
        parts{cellID} = char(string(excel_scalarize(x{cellID})));
    end
    out = strjoin(parts, '; ');
    return
end

if isnumeric(x) || islogical(x)
    if isscalar(x)
        out = x;
        return
    end

    if ndims(x) <= 2
        out = mat2str(x, 6);
        return
    end

    out = sprintf('<%s %s, min=%g, max=%g, mean=%g, nnz=%d>', ...
        class(x), size_vec_to_string(size(x)), min(x(:), [], 'omitnan'), ...
        max(x(:), [], 'omitnan'), mean(double(x(:)), 'omitnan'), nnz(x));
    return
end

if ischar(x)
    out = x;
    return
end

if isstring(x)
    if isscalar(x)
        out = char(x);
    else
        out = strjoin(cellstr(x), '; ');
    end
    return
end

if isstruct(x)
    out = struct_to_string(x);
    return
end

if isobject(x)
    if isprop(x, 'T')
        try
            out = mat2str(x.T, 6);
            return
        catch
        end
    end

    if isprop(x, 'A')
        try
            out = mat2str(x.A, 6);
            return
        catch
        end
    end

    try
        out = char(string(x));
    catch
        out = ['<' class(x) '>'];
    end
    return
end

try
    out = char(string(x));
catch
    out = ['<' class(x) '>'];
end

end


function out = struct_to_string(S)

if numel(S) > 1
    out = sprintf('<%s %s>', class(S), size_vec_to_string(size(S)));
    return
end

fieldNames = fieldnames(S);
parts = cell(numel(fieldNames), 1);

for fieldID = 1:numel(fieldNames)
    fieldName = fieldNames{fieldID};
    parts{fieldID} = sprintf('%s=%s', fieldName, char(string(excel_scalarize(S.(fieldName)))));
end

out = strjoin(parts, '; ');

end


function s = size_vec_to_string(sz)

if isempty(sz)
    s = '[]';
else
    s = ['[' strtrim(sprintf('%d ', sz)) ']'];
end

end


function headerNames = get_export_header_names(variableNames)

headerNames = variableNames;

headerMap = {
    'animal_ID', 'animal.ID'
    'animal_hemisphere', 'animal.hemisphere'
    'date_camprojcalib', 'date.camprojcalib'
    'run_camprojcalib', 'run.camprojcalib'
    'date_opto', 'date.opto'
    'run_opto', 'run.opto'
    'date_baseline', 'date.baseline'
    'run_baseline', 'run.baseline'
    'date_ortmap', 'date.ortmap'
    'run_ortmap', 'run.ortmap'
    'date_visfootprint', 'date.visfootprint'
    'run_visfootprint', 'run.visfootprint'
    'vis_contrast', 'vis.contrast'
    'vis_contrast_baseline', 'vis.contrast_baseline'
    'vis_contrast_conopto', 'vis.contrast_conopto'
    'vis_contrast_inconopto', 'vis.contrast_inconopto'
    'vis_ort', 'vis.ort'
    'vis_sz', 'vis.sz'
    'vis_sf', 'vis.sf'
    'vis_phs', 'vis.phs'
    'vis_pos', 'vis.pos'
    'opto_480LED', 'opto.480LED'
    'opto_480ND', 'opto.480ND'
    'opto_580LED', 'opto.580LED'
    'opto_580ND', 'opto.580ND'
    'opto_gridsize', 'opto.gridsize'
    'opto_gamma', 'opto.gamma'
    'opto_threshSens', 'opto.threshSens'
    'opto_threshAdapt', 'opto.threshAdapt'
    'opto_ort', 'opto.ort'
    'opto_gausscond', 'opto.gausscond'
    'opto_gausslevel', 'opto.gausslevel'
    'opto_gaussmax', 'opto.gaussmax'
    'opto_transformParams', 'opto.transformParams'
    'opto_bitmapCamSpace', 'opto.bitmapCamSpace'
    'opto_bitmapProjSpace', 'opto.bitmapProjSpace'
    'opto_nColumns', 'opto.nColumns'
    'opto_columnArea', 'opto.columnArea'
    'opto_pixelsON', 'opto.pixelsON'
    'opto_temporalDutyCycle', 'opto.temporalDutyCycle'
    'opto_meanPowerDensityWithinROI_mWmm2', 'opto.meanPowerDensityWithinROI_mWmm2'
    'opto_totalPowerToOnPixelsWithinROI_mW', 'opto.totalPowerToOnPixelsWithinROI_mW'
    'opto_projectorPowerDensity_mWmm2', 'opto.projectorPowerDensity_mWmm2'
    'psy_deltaMask', 'psy.deltaMask'
    'psy_deltaBias', 'psy.deltaBias'
    };

for mapID = 1:size(headerMap, 1)
    varIdx = strcmp(variableNames, headerMap{mapID, 1});
    headerNames(varIdx) = headerMap(mapID, 2);
end

end
