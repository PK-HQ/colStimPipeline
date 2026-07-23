function MetaTable = buildBlockMetadata(behavioralData, bitmapData, columnsDesired, columnsSpread, blockData, datastruct, analysisBlockID, mdlStruct, saveOpts)

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
if nargin < 9 || isempty(saveOpts)
    saveOpts = struct();
end

meanCol = mean(bitmapData.nColumns, 1, 'omitnan');
nTotalBlocks = numel(meanCol);
% analysisBlockID is used only as a metadata/datastruct label map (blockID → dsIdx).
% bitmapData and behavioralData are always indexed by local rows 1:nTotalBlocks.
assert(numel(analysisBlockID) >= nTotalBlocks, ...
    'analysisBlockID has %d elements but nTotalBlocks=%d (from bitmapData.nColumns)', ...
    numel(analysisBlockID), nTotalBlocks);
selectedBlocks = 1:nTotalBlocks;
nBlocks = nTotalBlocks;
assert(all(selectedBlocks >= 1 & selectedBlocks <= nTotalBlocks), ...
    'Selected block indices exceed 1:%d', nTotalBlocks);
[deltaMask, deltaBias] = get_psychometric_deltas(mdlStruct, nTotalBlocks);

[powerClusterID, powerClusterSignificant] = ...
    get_power_cluster_metadata(mdlStruct, nTotalBlocks);
psyFull = get_full_psychometrics_struct(mdlStruct, nTotalBlocks, behavioralData);

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
    row.gridSize_mean    = safe_col_mean(bitmapData, 'gridSize',     blockID);
    row.nColumns_Mean    = safe_col_mean(bitmapData, 'nColumns',    blockID);
    row.nColumns_min     = safe_col_min (bitmapData, 'nColumns',    blockID);
    row.nColumns_max     = safe_col_max (bitmapData, 'nColumns',    blockID);
    row.sensitivity_mean = safe_col_mean(bitmapData, 'sensitivity', blockID);
    row.adaptthresh_mean = safe_col_mean(bitmapData, 'adaptthresh', blockID);
    row.pixelsON_mean    = safe_col_mean(bitmapData, 'pixelsON',    blockID);

    if isfield(bitmapData, 'meanPowerDensityWithinROI_mWmm2')
        row.meanPowerDensityWithinROI_mWmm2_mean = mean(bitmapData.meanPowerDensityWithinROI_mWmm2(:,:,blockID), 'all', 'omitnan');
        row.meanPowerDensityWithinROI_mWmm2_max  = max(bitmapData.meanPowerDensityWithinROI_mWmm2(:,:,blockID), [], 'all');
    end
    if isfield(bitmapData, 'totalPowerToOnPixelsWithinROI_mW')
        row.totalPowerToOnPixelsWithinROI_mW_mean = mean(bitmapData.totalPowerToOnPixelsWithinROI_mW(:,:,blockID), 'all', 'omitnan');
        row.totalPowerToOnPixelsWithinROI_mW_max  = max(bitmapData.totalPowerToOnPixelsWithinROI_mW(:,:,blockID), [], 'all');
    end
    if isfield(bitmapData, 'projectorPowerDensity_mWmm2')
        row.projectorPowerDensity_mWmm2_mean = mean(bitmapData.projectorPowerDensity_mWmm2(:,:,blockID), 'all', 'omitnan');
        row.projectorPowerDensity_mWmm2_max  = max(bitmapData.projectorPowerDensity_mWmm2(:,:,blockID), [], 'all');
    end
    if isfield(bitmapData, 'temporalDutyCycle')
        row.temporalDutyCycle_mean = mean(bitmapData.temporalDutyCycle(:,:,blockID), 'all', 'omitnan');
    end
    row = add_bitmap_metadata(row, bitmapData, blockID, nTotalBlocks);
    row.psy_deltaMask = deltaMask(blockID);
    row.psy_deltaBias = deltaBias(blockID);
    
    row.psy_powerClusterID = powerClusterID(blockID);
    row.psy_powerClusterSignificant = powerClusterSignificant(blockID);

    row = add_psychometric_full_data(row, psyFull, blockID);

    if isempty(row.session_baselineSource)
        row.session_baselineSource = infer_baseline_source_from_row(row);
    end

    [bmpHCam, bmpVCam, bmpHProj, bmpVProj] = ...
        get_oriented_bitmaps_for_block(bitmapData, blockID, nTotalBlocks);
    row.bmp_horizontalCamSpace{1}  = bmpHCam;
    row.bmp_verticalCamSpace{1}    = bmpVCam;
    row.bmp_horizontalProjSpace{1} = bmpHProj;
    row.bmp_verticalProjSpace{1}   = bmpVProj;

    fpDate = row.date_visfootprint{1};
    fpRun  = row.run_visfootprint{1};
    hasFPS = ~isempty(fpDate) && ~isempty(fpRun);
    if hasFPS && isnumeric(fpDate) && all(isnan(fpDate(:))), hasFPS = false; end
    if hasFPS && isnumeric(fpRun)  && all(isnan(fpRun(:))),  hasFPS = false; end
    row.session_hasVisFPS = hasFPS;

    rows(ii) = row;
end

MetaTable = struct2table(rows);
MetaTable.Properties.VariableDescriptions = get_export_header_names(MetaTable.Properties.VariableNames);
MetaTable = filter_metatable_by_columns(MetaTable, columnsDesired, columnsSpread);
MetaTable = remove_metatable_output_columns(MetaTable);
MetaTable = move_psychometric_columns_to_end(MetaTable);
save_metatable_outputs(MetaTable, datastruct, analysisBlockID, columnsDesired, saveOpts);

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
row.opto_PCAdenoisedResp = {[]};
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
row.psy_powerClusterID = NaN;
row.psy_powerClusterSignificant = NaN;
row.psy_deltaBiasHorizontal = NaN;
row.psy_deltaMaskHorizontal = NaN;
row.psy_deltaBiasVertical   = NaN;
row.psy_deltaMaskVertical   = NaN;
row.psy_deltaBiasMerged     = NaN;
row.psy_deltaMaskMerged     = NaN;
row.psy_combinedBL          = false;
row.psy_modelField          = '';
row.session_baselineSource  = '';
row.session_hasVisFPS       = false;
row.psy_xBaselinePreMerge               = {[]};
row.psy_yBaselinePreMerge               = {[]};
row.psy_xBaselinePreMergeOrt            = {[]};
row.psy_nTrialsBaselinePreMerge         = {[]};
row.psy_xHorizontalOptoPreMerge         = {[]};
row.psy_yHorizontalOptoPreMerge         = {[]};
row.psy_xHorizontalOptoPreMergeOrt      = {[]};
row.psy_nTrialsHorizontalOptoPreMerge   = {[]};
row.psy_congruencyHorizontalOptoPreMerge = {[]};
row.psy_xVerticalOptoPreMerge           = {[]};
row.psy_yVerticalOptoPreMerge           = {[]};
row.psy_xVerticalOptoPreMergeOrt        = {[]};
row.psy_nTrialsVerticalOptoPreMerge     = {[]};
row.psy_congruencyVerticalOptoPreMerge  = {[]};
row.psy_xBaselineMerged                 = {[]};
row.psy_yBaselineMerged                 = {[]};
row.psy_nTrialsBaselineMerged           = {[]};
row.psy_xConOptoMerged                  = {[]};
row.psy_yConOptoMerged                  = {[]};
row.psy_nTrialsConOptoMerged            = {[]};
row.psy_xInconOptoMerged                = {[]};
row.psy_yInconOptoMerged                = {[]};
row.psy_nTrialsInconOptoMerged          = {[]};
row.bmp_horizontalCamSpace  = {[]};
row.bmp_verticalCamSpace    = {[]};
row.bmp_horizontalProjSpace = {[]};
row.bmp_verticalProjSpace   = {[]};
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
row.opto_transformParams = ...
    {get_bitmap_field( ...
        bitmapData, ...
        'transformParams', ...
        blockID, ...
        nTotalBlocks)};

% Normally 512 x 512 x 2 for one experiment.
row.opto_PCAdenoisedResp = ...
    {get_bitmap_field( ...
        bitmapData, ...
        'pcadenoisedresp', ...
        blockID, ...
        nTotalBlocks)};

row.opto_bitmapCamSpace = ...
    {get_bitmap_field(bitmapData, 'columnarbitmapTFcamspace', blockID, nTotalBlocks)};
row.opto_bitmapProjSpace = {get_bitmap_field(bitmapData, 'columnarbitmapTFprojspace', blockID, nTotalBlocks)};
row.opto_nColumns = {get_bitmap_field(bitmapData, 'nColumns', blockID, nTotalBlocks)};
row.opto_columnArea = {combine_column_areas(get_bitmap_field(bitmapData, 'columnAreas', blockID, nTotalBlocks))};
row.opto_pixelsON = {get_bitmap_field(bitmapData, 'pixelsON', blockID, nTotalBlocks)};
row.opto_temporalDutyCycle = {get_bitmap_field(bitmapData, 'temporalDutyCycle', blockID, nTotalBlocks)};
row.opto_meanPowerDensityWithinROI_mWmm2 = {get_bitmap_field(bitmapData, 'meanPowerDensityWithinROI_mWmm2', blockID, nTotalBlocks)};
row.opto_totalPowerToOnPixelsWithinROI_mW = {get_bitmap_field(bitmapData, 'totalPowerToOnPixelsWithinROI_mW', blockID, nTotalBlocks)};
row.opto_projectorPowerDensity_mWmm2 = {get_bitmap_field(bitmapData, 'projectorPowerDensity_mWmm2', blockID, nTotalBlocks)};

end


function v = safe_col_mean(S, fn, bID)
v = NaN;
if isfield(S, fn)
    col = S.(fn);
    if size(col,2) >= bID
        v = mean(col(:,bID), 'omitnan');
    end
end
end

function v = safe_col_min(S, fn, bID)
v = NaN;
if isfield(S, fn)
    col = S.(fn);
    if size(col,2) >= bID
        v = min(col(:,bID), [], 'omitnan');
    end
end
end

function v = safe_col_max(S, fn, bID)
v = NaN;
if isfield(S, fn)
    col = S.(fn);
    if size(col,2) >= bID
        v = max(col(:,bID), [], 'omitnan');
    end
end
end

function val = get_bitmap_field(bitmapData, fieldName, blockID, nTotalBlocks)

val = [];

if ~isfield(bitmapData, fieldName)
    return
end

try
    val = get_block_slice(bitmapData.(fieldName), blockID, nTotalBlocks);
catch
    % Field exists but block dimension doesn't match nTotalBlocks.
    % Try to index directly if size allows.
    rawField = bitmapData.(fieldName);
    if numel(rawField) >= blockID
        try
            val = get_block_slice(rawField, blockID, numel(rawField));
        catch
        end
    end
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

psyNames = {
    'psy_deltaMask'
    'psy_deltaBias'
    'psy_powerClusterID'
    'psy_powerClusterSignificant'
    };

psyNames = psyNames(ismember(psyNames, ...
    MetaTable.Properties.VariableNames));

if isempty(psyNames)
    return
end

remainingNames = setdiff( ...
    MetaTable.Properties.VariableNames, psyNames, 'stable');

if isempty(remainingNames)
    return
end

MetaTable = movevars( ...
    MetaTable, psyNames, 'After', remainingNames{end});

end


function save_metatable_outputs(MetaTable, datastruct, analysisBlockID, columnsDesired, saveOpts)

if nargin < 5
    saveOpts = struct();
end

% Use custom paths from saveOpts when provided
if isstruct(saveOpts) && isfield(saveOpts, 'matPath') && isfield(saveOpts, 'xlsxPath')
    matOutPath  = saveOpts.matPath;
    xlsxOutPath = saveOpts.xlsxPath;
    outDir = fileparts(matOutPath);
    if ~isempty(outDir) && ~exist(outDir, 'dir')
        mkdir(outDir);
    end
else
    % Auto-derive path from datastruct (original behaviour)
    if ~isstruct(datastruct) || isempty(datastruct) || ...
            ~isfield(datastruct, 'monkey') || ~isfield(datastruct, 'chamber')
        return
    end
    saveEntry = datastruct(1);
    if ~isempty(analysisBlockID) && analysisBlockID(1) >= 1 && numel(datastruct) >= analysisBlockID(1)
        saveEntry = datastruct(analysisBlockID(1));
    end
    monkeyName  = saveEntry.monkey;
    chamberStr  = saveEntry.chamber;
    columnsStr  = num2str(columnsDesired);
    outDir      = ['Y:/' monkeyName '/Meta/summary'];
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end
    matOutPath  = [outDir '/metaTable-' chamberStr columnsStr '.mat'];
    xlsxOutPath = [outDir '/metaTable-' chamberStr columnsStr '.xlsx'];
end

save(matOutPath, 'MetaTable', '-v7.3');

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
        if numel(x) <= 100
            out = mat2str(x, 6);
        else
            out = sprintf('%dx%d %s', size(x,1), size(x,2), class(x));
        end
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
    'opto_PCAdenoisedResp', 'opto.PCAdenoisedResp'
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
    'psy_powerClusterID', 'psy.powerClusterID'
    'psy_powerClusterSignificant', 'psy.powerClusterSignificant'
    'psy_deltaBiasHorizontal', 'psy.deltaBiasHorizontal'
    'psy_deltaMaskHorizontal', 'psy.deltaMaskHorizontal'
    'psy_deltaBiasVertical', 'psy.deltaBiasVertical'
    'psy_deltaMaskVertical', 'psy.deltaMaskVertical'
    'psy_deltaBiasMerged', 'psy.deltaBiasMerged'
    'psy_deltaMaskMerged', 'psy.deltaMaskMerged'
    'psy_combinedBL', 'psy.combinedBL'
    'psy_modelField', 'psy.modelField'
    'psy_xBaselinePreMerge', 'psy.xBaselinePreMerge'
    'psy_yBaselinePreMerge', 'psy.yBaselinePreMerge'
    'psy_xBaselinePreMergeOrt', 'psy.xBaselinePreMergeOrt'
    'psy_nTrialsBaselinePreMerge', 'psy.nTrialsBaselinePreMerge'
    'psy_xHorizontalOptoPreMerge', 'psy.xHorizontalOptoPreMerge'
    'psy_yHorizontalOptoPreMerge', 'psy.yHorizontalOptoPreMerge'
    'psy_xHorizontalOptoPreMergeOrt', 'psy.xHorizontalOptoPreMergeOrt'
    'psy_nTrialsHorizontalOptoPreMerge', 'psy.nTrialsHorizontalOptoPreMerge'
    'psy_congruencyHorizontalOptoPreMerge', 'psy.congruencyHorizontalOptoPreMerge'
    'psy_xVerticalOptoPreMerge', 'psy.xVerticalOptoPreMerge'
    'psy_yVerticalOptoPreMerge', 'psy.yVerticalOptoPreMerge'
    'psy_xVerticalOptoPreMergeOrt', 'psy.xVerticalOptoPreMergeOrt'
    'psy_nTrialsVerticalOptoPreMerge', 'psy.nTrialsVerticalOptoPreMerge'
    'psy_congruencyVerticalOptoPreMerge', 'psy.congruencyVerticalOptoPreMerge'
    'psy_xBaselineMerged', 'psy.xBaselineMerged'
    'psy_yBaselineMerged', 'psy.yBaselineMerged'
    'psy_nTrialsBaselineMerged', 'psy.nTrialsBaselineMerged'
    'psy_xConOptoMerged', 'psy.xConOptoMerged'
    'psy_yConOptoMerged', 'psy.yConOptoMerged'
    'psy_nTrialsConOptoMerged', 'psy.nTrialsConOptoMerged'
    'psy_xInconOptoMerged', 'psy.xInconOptoMerged'
    'psy_yInconOptoMerged', 'psy.yInconOptoMerged'
    'psy_nTrialsInconOptoMerged', 'psy.nTrialsInconOptoMerged'
    'session_baselineSource', 'session.baselineSource'
    'session_hasVisFPS', 'session.hasVisFPS'
    'bmp_horizontalCamSpace', 'bmp.horizontalCamSpace'
    'bmp_verticalCamSpace', 'bmp.verticalCamSpace'
    'bmp_horizontalProjSpace', 'bmp.horizontalProjSpace'
    'bmp_verticalProjSpace', 'bmp.verticalProjSpace'
    };

for mapID = 1:size(headerMap, 1)
    varIdx = strcmp(variableNames, headerMap{mapID, 1});
    headerNames(varIdx) = headerMap(mapID, 2);
end

end


% ---- NEW LOCAL FUNCTIONS -----------------------------------------------

function psyFull = get_full_psychometrics_struct(mdlStruct, nTotalBlocks, behavioralData)
% For every block (1..nTotalBlocks) find the model field that covers it
% and extract all psychometric data.  Returns a struct array [nTotalBlocks x 1].
% behavioralData must contain gaborContrasts [3x12xN] and percentageCorrect [3x12xN].
if nargin < 3
    behavioralData = struct();
end

emptyEntry.modelField               = '';
emptyEntry.combinedBL               = false;
emptyEntry.baselineSource           = '';
emptyEntry.xBaselinePreMerge        = [];
emptyEntry.yBaselinePreMerge        = [];
emptyEntry.xBaselinePreMergeOrt     = [];
emptyEntry.nTrialsBaselinePreMerge  = [];
emptyEntry.xHorizontalOptoPreMerge  = [];
emptyEntry.yHorizontalOptoPreMerge  = [];
emptyEntry.xHorizontalOptoPreMergeOrt       = [];
emptyEntry.nTrialsHorizontalOptoPreMerge    = [];
emptyEntry.congruencyHorizontalOptoPreMerge = [];
emptyEntry.xVerticalOptoPreMerge    = [];
emptyEntry.yVerticalOptoPreMerge    = [];
emptyEntry.xVerticalOptoPreMergeOrt         = [];
emptyEntry.nTrialsVerticalOptoPreMerge      = [];
emptyEntry.congruencyVerticalOptoPreMerge   = [];
emptyEntry.xBaselineMerged          = [];
emptyEntry.yBaselineMerged          = [];
emptyEntry.nTrialsBaselineMerged    = [];
emptyEntry.xConOptoMerged           = [];
emptyEntry.yConOptoMerged           = [];
emptyEntry.nTrialsConOptoMerged     = [];
emptyEntry.xInconOptoMerged         = [];
emptyEntry.yInconOptoMerged         = [];
emptyEntry.nTrialsInconOptoMerged   = [];
emptyEntry.deltaBias                = NaN;
emptyEntry.deltaMask                = NaN;
emptyEntry.deltaBiasHorizontal      = NaN;
emptyEntry.deltaMaskHorizontal      = NaN;
emptyEntry.deltaBiasVertical        = NaN;
emptyEntry.deltaMaskVertical        = NaN;
emptyEntry.deltaBiasMerged          = NaN;
emptyEntry.deltaMaskMerged          = NaN;

psyFull = repmat(emptyEntry, nTotalBlocks, 1);

if ~isstruct(mdlStruct) || isempty(fieldnames(mdlStruct))
    return
end

mdlFields = fieldnames(mdlStruct);

for fIdx = 1:numel(mdlFields)
    fName = mdlFields{fIdx};
    md    = mdlStruct.(fName);
    if ~isstruct(md) || ~isfield(md, 'clusterBlocksIdx')
        continue
    end
    blockIDs = md.clusterBlocksIdx(:);
    nK = numel(blockIDs);

    for kIdx = 1:nK
        bID = blockIDs(kIdx);
        if bID < 1 || bID > nTotalBlocks
            continue
        end

        e = psyFull(bID);
        e.modelField = fName;

        e = set_scalar_field(e, md, 'deltaBias',           kIdx, 1);
        e = set_scalar_field(e, md, 'deltaMask',           kIdx, 1);
        e = set_scalar_field(e, md, 'deltaBiasHorizontal', kIdx, 1);
        e = set_scalar_field(e, md, 'deltaMaskHorizontal', kIdx, 1);
        e = set_scalar_field(e, md, 'deltaBiasVertical',   kIdx, 1);
        e = set_scalar_field(e, md, 'deltaMaskVertical',   kIdx, 1);
        e = set_scalar_field(e, md, 'deltaBiasMerged',     kIdx, 1);
        e = set_scalar_field(e, md, 'deltaMaskMerged',     kIdx, 1);

        e = set_row_field(e, md, 'xBaselinePreMerge',                        kIdx);
        e = set_row_field(e, md, 'yBaselinePreMerge',                        kIdx);
        e = set_row_field_as(e, md, 'visualTagBaselinePreMerge',       'xBaselinePreMergeOrt',       kIdx);
        e = set_row_field(e, md, 'xHorizontalOptoPreMerge',                  kIdx);
        e = set_row_field(e, md, 'yHorizontalOptoPreMerge',                  kIdx);
        e = set_row_field_as(e, md, 'visualTagHorizontalOptoPreMerge', 'xHorizontalOptoPreMergeOrt', kIdx);
        e = set_row_field(e, md, 'congruencyHorizontalOptoPreMerge',         kIdx);
        e = set_row_field(e, md, 'xVerticalOptoPreMerge',                    kIdx);
        e = set_row_field(e, md, 'yVerticalOptoPreMerge',                    kIdx);
        e = set_row_field_as(e, md, 'visualTagVerticalOptoPreMerge',   'xVerticalOptoPreMergeOrt',   kIdx);
        e = set_row_field(e, md, 'congruencyVerticalOptoPreMerge',           kIdx);
        e = set_row_field_as(e, md, 'xBaseline',   'xBaselineMerged',   kIdx);
        e = set_row_field_as(e, md, 'yBaseline',   'yBaselineMerged',   kIdx);
        e = set_row_field_as(e, md, 'xConOpto',    'xConOptoMerged',    kIdx);
        e = set_row_field_as(e, md, 'yConOpto',    'yConOptoMerged',    kIdx);
        e = set_row_field_as(e, md, 'xInconOpto',  'xInconOptoMerged',  kIdx);
        e = set_row_field_as(e, md, 'yInconOpto',  'yInconOptoMerged',  kIdx);
        e = compute_nominal_trial_counts(e);

        if isfield(md, 'combinedBL')
            cbl = md.combinedBL;
            if isvector(cbl) && numel(cbl) >= kIdx
                e.combinedBL = logical(cbl(kIdx));
            elseif isscalar(cbl)
                e.combinedBL = logical(cbl);
            end
        end

        if e.combinedBL
            e.baselineSource = 'same_block_as_opto';
        elseif ~isempty(e.xBaselinePreMerge)
            e.baselineSource = 'separate_block';
        end

        psyFull(bID) = e;
    end
end

% Second pass: fill behavioral arrays from raw data for blocks not covered by any cluster.
% behavioralData.gaborContrasts [3x12xN] = xBlocks, .percentageCorrect [3x12xN] = yBlocks.
hasRaw = isstruct(behavioralData) && ...
         isfield(behavioralData, 'gaborContrasts') && ...
         isfield(behavioralData, 'percentageCorrect');
if hasRaw
    xBlocks = behavioralData.gaborContrasts;
    yBlocks = behavioralData.percentageCorrect;
    nRaw = size(xBlocks, 3);
    for bID = 1:min(nTotalBlocks, nRaw)
        if ~isempty(psyFull(bID).modelField)
            continue
        end
        try
            e = psyFull(bID);
            e = fill_raw_behavioral_entry(e, xBlocks, yBlocks, bID);
            e = compute_nominal_trial_counts(e);
            e = compute_deltas_for_raw_entry(e);
            e.deltaBias = e.deltaBiasMerged;
            e.deltaMask = e.deltaMaskMerged;
            e.modelField = 'behavioralData_raw';
            if isempty(e.baselineSource) && ~isempty(e.xBaselineMerged)
                e.baselineSource = 'separate_block';
            end
            psyFull(bID) = e;
        catch ME
            fprintf('WARNING: fill_raw_behavioral_entry failed for block %d: %s\n', bID, ME.message);
        end
    end
end

end


function e = set_scalar_field(e, md, fname, kIdx, dim)
% Extract scalar at position kIdx from field fname (1-D vector along dim).
if ~isfield(md, fname)
    return
end
val = md.(fname);
if isempty(val)
    return
end
if dim == 1
    vec = val(:);
else
    vec = val(kIdx, :);
end
if numel(vec) >= kIdx
    scalar = vec(kIdx);
    if isnumeric(scalar) || islogical(scalar)
        e.(fname) = scalar;
    end
end
end


function e = set_row_field(e, md, fname, kIdx)
% Extract row kIdx from a 2-D matrix field, stripping trailing NaNs.
if ~isfield(md, fname)
    return
end
mat = md.(fname);
if isempty(mat)
    return
end
if isvector(mat)
    row = mat(:)';
elseif ismatrix(mat) && size(mat,1) >= kIdx
    row = mat(kIdx, :);
else
    return
end
if isnumeric(row)
    row = row(~isnan(row));
end
e.(fname) = row;
end


function row = add_psychometric_full_data(row, psyFull, blockID)
% Copy all psychometric fields from psyFull(blockID) into the metatable row.
if blockID < 1 || blockID > numel(psyFull)
    return
end
e = psyFull(blockID);

row.psy_modelField              = e.modelField;
row.psy_combinedBL              = e.combinedBL;
row.session_baselineSource      = e.baselineSource;

row.psy_xBaselinePreMerge               = {e.xBaselinePreMerge};
row.psy_yBaselinePreMerge               = {e.yBaselinePreMerge};
row.psy_xBaselinePreMergeOrt            = {e.xBaselinePreMergeOrt};
row.psy_nTrialsBaselinePreMerge         = {e.nTrialsBaselinePreMerge};
row.psy_xHorizontalOptoPreMerge         = {e.xHorizontalOptoPreMerge};
row.psy_yHorizontalOptoPreMerge         = {e.yHorizontalOptoPreMerge};
row.psy_xHorizontalOptoPreMergeOrt      = {e.xHorizontalOptoPreMergeOrt};
row.psy_nTrialsHorizontalOptoPreMerge   = {e.nTrialsHorizontalOptoPreMerge};
row.psy_congruencyHorizontalOptoPreMerge = {e.congruencyHorizontalOptoPreMerge};
row.psy_xVerticalOptoPreMerge           = {e.xVerticalOptoPreMerge};
row.psy_yVerticalOptoPreMerge           = {e.yVerticalOptoPreMerge};
row.psy_xVerticalOptoPreMergeOrt        = {e.xVerticalOptoPreMergeOrt};
row.psy_nTrialsVerticalOptoPreMerge     = {e.nTrialsVerticalOptoPreMerge};
row.psy_congruencyVerticalOptoPreMerge  = {e.congruencyVerticalOptoPreMerge};
row.psy_xBaselineMerged                 = {e.xBaselineMerged};
row.psy_yBaselineMerged                 = {e.yBaselineMerged};
row.psy_nTrialsBaselineMerged           = {e.nTrialsBaselineMerged};
row.psy_xConOptoMerged                  = {e.xConOptoMerged};
row.psy_yConOptoMerged                  = {e.yConOptoMerged};
row.psy_nTrialsConOptoMerged            = {e.nTrialsConOptoMerged};
row.psy_xInconOptoMerged                = {e.xInconOptoMerged};
row.psy_yInconOptoMerged                = {e.yInconOptoMerged};
row.psy_nTrialsInconOptoMerged          = {e.nTrialsInconOptoMerged};

row.psy_deltaBiasHorizontal = e.deltaBiasHorizontal;
row.psy_deltaMaskHorizontal = e.deltaMaskHorizontal;
row.psy_deltaBiasVertical   = e.deltaBiasVertical;
row.psy_deltaMaskVertical   = e.deltaMaskVertical;
row.psy_deltaBiasMerged     = e.deltaBiasMerged;
row.psy_deltaMaskMerged     = e.deltaMaskMerged;
if ~isnan(e.deltaBias)
    row.psy_deltaBias = e.deltaBias;
end
if ~isnan(e.deltaMask)
    row.psy_deltaMask = e.deltaMask;
end

end


function [hCam, vCam, hProj, vProj] = get_oriented_bitmaps_for_block(bitmapData, blockID, nTotalBlocks)
% Extract horizontal- and vertical-oriented bitmaps (cam-space and proj-space).
hCam  = [];
vCam  = [];
hProj = [];
vProj = [];

if ~isfield(bitmapData, 'orts')
    return
end

ortsAll = bitmapData.orts;
% shape [1 x 2 x nBlocks]
if ndims(ortsAll) ~= 3 || size(ortsAll,3) < blockID
    return
end
ortsBlock = squeeze(ortsAll(:, :, blockID));  % [1x2] or [2x1] → vector
ortsBlock  = ortsBlock(:);                    % ensure column [2x1]

if numel(ortsBlock) < 2
    return
end

% 0 deg = horizontal, 90 deg = vertical
idxH = find(ortsBlock == 0,  1);
idxV = find(ortsBlock == 90, 1);

if isempty(idxH) && isempty(idxV)
    % fallback: assign by position (index 1 = hor, index 2 = vert)
    idxH = 1;
    idxV = 2;
end

% Camera-space bitmaps: [H x W x 2 x nBlocks]
if isfield(bitmapData, 'columnarbitmapTFcamspace')
    bmpCam = bitmapData.columnarbitmapTFcamspace;
    % find block dimension
    szCam = size(bmpCam);
    blockDimCam = find(szCam == nTotalBlocks, 1, 'last');
    if ~isempty(blockDimCam)
        idxC = repmat({':'}, 1, ndims(bmpCam));
        idxC{blockDimCam} = blockID;
        bmpBlock = bmpCam(idxC{:});   % [H x W x 2] or [H x W] depending on ndims
        bmpBlock = squeeze(bmpBlock);
        % orientation dimension: ndims of squeezed result
        if ndims(bmpBlock) == 3
            % 3rd dim is orientation
            if ~isempty(idxH) && size(bmpBlock,3) >= idxH
                hCam = bmpBlock(:,:,idxH);
            end
            if ~isempty(idxV) && size(bmpBlock,3) >= idxV
                vCam = bmpBlock(:,:,idxV);
            end
        elseif ismatrix(bmpBlock)
            hCam = bmpBlock;
            vCam = bmpBlock;
        end
    end
end

% Projector-space bitmaps
if isfield(bitmapData, 'columnarbitmapTFprojspace')
    bmpProj = bitmapData.columnarbitmapTFprojspace;
    szProj = size(bmpProj);
    blockDimProj = find(szProj == nTotalBlocks, 1, 'last');
    if ~isempty(blockDimProj)
        idxP = repmat({':'}, 1, ndims(bmpProj));
        idxP{blockDimProj} = blockID;
        bmpPBlock = bmpProj(idxP{:});
        bmpPBlock = squeeze(bmpPBlock);
        if ndims(bmpPBlock) == 3
            if ~isempty(idxH) && size(bmpPBlock,3) >= idxH
                hProj = bmpPBlock(:,:,idxH);
            end
            if ~isempty(idxV) && size(bmpPBlock,3) >= idxV
                vProj = bmpPBlock(:,:,idxV);
            end
        elseif ismatrix(bmpPBlock)
            hProj = bmpPBlock;
            vProj = bmpPBlock;
        end
    end
end

end


function e = fill_raw_behavioral_entry(e, xBlocks, yBlocks, blockID)
% Replicates processConditionsBlocks logic for one block using raw behavioral data.
% xBlocks [3 x nContrasts x nBlocks]: row1=baseline, row2=horizontal, row3=vertical.
% yBlocks same shape with percent-correct values.

xBaselineRaw = bmd_rmnan(squeeze(xBlocks(1,:,blockID)));
yBaselineRaw = bmd_rmnan(squeeze(yBlocks(1,:,blockID)));
xHorizontalRaw = bmd_rmnan(squeeze(xBlocks(2,:,blockID)));
yHorizontalRaw = bmd_rmnan(squeeze(yBlocks(2,:,blockID)));
xVerticalRaw = bmd_rmnan(squeeze(xBlocks(3,:,blockID)));
yVerticalRaw = bmd_rmnan(squeeze(yBlocks(3,:,blockID)));

%% Baseline
nBase = numel(xBaselineRaw);
if nBase > 0 && mod(nBase, 2) == 0
    tagBaselineRaw = bmd_make_visual_tag(nBase);
    [xBaseSorted, sortIdx] = sort(xBaselineRaw);
    yBaseSorted = yBaselineRaw(sortIdx);
    tagBaseSorted = tagBaselineRaw(sortIdx);

    % Pre-merge baseline
    [xBPre, yBPre, tagBPre] = bmd_make_pre_merge_pct_correct(xBaseSorted, yBaseSorted, tagBaseSorted, true);
    e.xBaselinePreMerge = xBPre(:)';
    e.yBaselinePreMerge = yBPre(:)';
    e.xBaselinePreMergeOrt = tagBPre(:)';

    % Merged baseline (averaged symmetric halves → percent correct)
    numVal = numel(xBaseSorted);
    xBMerged = mean([fliplr(-xBaseSorted(1:numVal/2)); xBaseSorted(numVal/2+1:end)]);
    yBMerged  = mean([fliplr(100 - yBaseSorted(1:numVal/2)); yBaseSorted(numVal/2+1:end)]);
    e.xBaselineMerged = xBMerged(:)';
    e.yBaselineMerged  = yBMerged(:)';
end

%% Horizontal opto
nH = numel(xHorizontalRaw);
if nH > 0 && mod(nH, 2) == 0
    tagHRaw = bmd_make_visual_tag(nH);
    [xHSorted, sortIdxH] = sort(xHorizontalRaw);
    yHSorted = yHorizontalRaw(sortIdxH);
    tagHSorted = tagHRaw(sortIdxH);

    [xHPre, yHPre, tagHPre] = bmd_make_pre_merge_pct_correct(xHSorted, yHSorted, tagHSorted, false);
    congrH = NaN(size(tagHPre));
    congrH(tagHPre == 0) = 1;
    congrH(tagHPre == 90) = -1;
    e.xHorizontalOptoPreMerge = xHPre(:)';
    e.yHorizontalOptoPreMerge = yHPre(:)';
    e.xHorizontalOptoPreMergeOrt = tagHPre(:)';
    e.congruencyHorizontalOptoPreMerge = congrH(:)';
end

%% Vertical opto
nV = numel(xVerticalRaw);
if nV > 0 && mod(nV, 2) == 0
    tagVRaw = bmd_make_visual_tag(nV);
    [xVSorted, sortIdxV] = sort(xVerticalRaw);
    yVSorted = yVerticalRaw(sortIdxV);
    tagVSorted = tagVRaw(sortIdxV);

    [xVPre, yVPre, tagVPre] = bmd_make_pre_merge_pct_correct(xVSorted, yVSorted, tagVSorted, false);
    congrV = NaN(size(tagVPre));
    congrV(tagVPre == 0) = -1;
    congrV(tagVPre == 90) = 1;
    e.xVerticalOptoPreMerge = xVPre(:)';
    e.yVerticalOptoPreMerge = yVPre(:)';
    e.xVerticalOptoPreMergeOrt = tagVPre(:)';
    e.congruencyVerticalOptoPreMerge = congrV(:)';
end

%% Merged con/incon (requires both H and V)
if nH > 0 && mod(nH, 2) == 0 && nV > 0 && mod(nV, 2) == 0
    contrastNeg = 1:numel(xHorizontalRaw)/2;
    contrastPos = numel(xHorizontalRaw)/2+1:numel(xHorizontalRaw);

    xHSorted2 = sort(xHorizontalRaw);
    yHSorted2 = yHorizontalRaw(argsort_vec(xHorizontalRaw));
    xVSorted2 = sort(xVerticalRaw);
    yVSorted2 = yVerticalRaw(argsort_vec(xVerticalRaw));

    yHCorrect = yHSorted2;
    yHCorrect(contrastNeg) = 100 - yHCorrect(contrastNeg);
    yVCorrect = yVSorted2;
    yVCorrect(contrastNeg) = 100 - yVCorrect(contrastNeg);

    xConOpto  = mean([-fliplr(xHSorted2(contrastNeg)); xVSorted2(contrastPos)]);
    yConOpto  = mean([fliplr(yHCorrect(contrastNeg));  yVCorrect(contrastPos)]);
    xInconOpto = mean([xHSorted2(contrastPos); -fliplr(xVSorted2(contrastNeg))]);
    yInconOpto  = mean([yHCorrect(contrastPos); fliplr(yVCorrect(contrastNeg))]);

    e.xConOptoMerged   = xConOpto(:)';
    e.yConOptoMerged   = yConOpto(:)';
    e.xInconOptoMerged = xInconOpto(:)';
    e.yInconOptoMerged = yInconOpto(:)';
end
end


function idx = argsort_vec(v)
[~, idx] = sort(v);
end


function v = bmd_rmnan(v)
v = v(~isnan(v));
end


function tag = bmd_make_visual_tag(nVal)
tag = NaN(1, nVal);
tag(1:nVal/2) = 0;
tag(nVal/2+1:end) = 90;
end


function [xOut, yOut, tagOut] = bmd_make_pre_merge_pct_correct(xIn, yIn, tagIn, mergeDupZeros)
xIn = xIn(:)'; yIn = yIn(:)'; tagIn = tagIn(:)';
valid = ~isnan(xIn) & ~isnan(yIn) & ~isnan(tagIn);
xIn = xIn(valid); yIn = yIn(valid); tagIn = tagIn(valid);
yCorrect = yIn;
yCorrect(tagIn == 0) = 100 - yCorrect(tagIn == 0);
if mergeDupZeros
    [xOut, yOut, tagOut] = bmd_merge_dup_x_for_baseline(xIn, yCorrect, tagIn);
else
    xOut = xIn; yOut = yCorrect; tagOut = tagIn;
end
end


function [xOut, yOut, tagOut] = bmd_merge_dup_x_for_baseline(xIn, yIn, tagIn)
xIn = xIn(:)'; yIn = yIn(:)'; tagIn = tagIn(:)';
valid = ~isnan(xIn) & ~isnan(yIn) & ~isnan(tagIn);
xIn = xIn(valid); yIn = yIn(valid); tagIn = tagIn(valid);
[xOut, ~, grp] = unique(xIn, 'stable');
yOut = nan(size(xOut)); tagOut = nan(size(xOut));
for ii = 1:numel(xOut)
    sel = grp == ii;
    yOut(ii) = mean(yIn(sel), 'omitnan');
    utags = unique(tagIn(sel));
    utags = utags(~isnan(utags));
    if numel(utags) == 1
        tagOut(ii) = utags;
    else
        tagOut(ii) = 45;
    end
end
end


function e = set_row_field_as(e, md, srcFname, dstFname, kIdx)
% Like set_row_field but reads from srcFname and writes to dstFname.
if ~isfield(md, srcFname)
    return
end
mat = md.(srcFname);
if isempty(mat)
    return
end
if isvector(mat)
    row = mat(:)';
elseif ismatrix(mat) && size(mat,1) >= kIdx
    row = mat(kIdx, :);
else
    return
end
if isnumeric(row)
    row = row(~isnan(row));
end
e.(dstFname) = row;
end


function e = compute_nominal_trial_counts(e)
% Assign nominal trial counts aligned one-for-one with x/y vectors.
% Pre-merge sides: 10 per point, 20 at abs(x)==0 (baseline after zero-merging).
% Merged: 20 per point, 40 at abs(x)==0.
e.nTrialsBaselinePreMerge       = bmd_make_side_weights(e.xBaselinePreMerge);
e.nTrialsHorizontalOptoPreMerge = 10 * ones(1, numel(e.xHorizontalOptoPreMerge));
e.nTrialsVerticalOptoPreMerge   = 10 * ones(1, numel(e.xVerticalOptoPreMerge));
e.nTrialsBaselineMerged         = bmd_make_merged_weights(e.xBaselineMerged);
e.nTrialsConOptoMerged          = bmd_make_merged_weights(e.xConOptoMerged);
e.nTrialsInconOptoMerged        = bmd_make_merged_weights(e.xInconOptoMerged);
end


function w = bmd_make_side_weights(x)
w = 10 * ones(1, numel(x));
w(abs(x) == 0) = 20;
end


function w = bmd_make_merged_weights(x)
w = 20 * ones(1, numel(x));
w(abs(x) == 0) = 40;
end


function e = compute_deltas_for_raw_entry(e)
% Compute all six delta scalars from x/y/ort arrays using nominal weights.
% Formula: deltaBias = conMean - inconMean; deltaMask = baseMean - mean(con,incon).

% Merged deltas
basM = bmd_weighted_mean(e.xBaselineMerged, e.yBaselineMerged, 'merged');
conM = bmd_weighted_mean(e.xConOptoMerged,  e.yConOptoMerged,  'merged');
incM = bmd_weighted_mean(e.xInconOptoMerged,e.yInconOptoMerged,'merged');
e.deltaBiasMerged = conM - incM;
e.deltaMaskMerged = basM - (conM + incM) / 2;

% Horizontal side deltas
ort = e.xBaselinePreMergeOrt;
xBH = e.xBaselinePreMerge(ort == 0);
yBH = e.yBaselinePreMerge(ort == 0);
ortH = e.xHorizontalOptoPreMergeOrt;
xCH = e.xHorizontalOptoPreMerge(ortH == 0);
yCH = e.yHorizontalOptoPreMerge(ortH == 0);
xIH = e.xHorizontalOptoPreMerge(ortH == 90);
yIH = e.yHorizontalOptoPreMerge(ortH == 90);
basH = bmd_weighted_mean(xBH, yBH, 'sideBaseline');
conH = bmd_weighted_mean(xCH, yCH, 'sideOpto');
incH = bmd_weighted_mean(xIH, yIH, 'sideOpto');
e.deltaBiasHorizontal = conH - incH;
e.deltaMaskHorizontal = basH - (conH + incH) / 2;

% Vertical side deltas
xBV = e.xBaselinePreMerge(ort == 90);
yBV = e.yBaselinePreMerge(ort == 90);
ortV = e.xVerticalOptoPreMergeOrt;
xCV = e.xVerticalOptoPreMerge(ortV == 90);
yCV = e.yVerticalOptoPreMerge(ortV == 90);
xIV = e.xVerticalOptoPreMerge(ortV == 0);
yIV = e.yVerticalOptoPreMerge(ortV == 0);
basV = bmd_weighted_mean(xBV, yBV, 'sideBaseline');
conV = bmd_weighted_mean(xCV, yCV, 'sideOpto');
incV = bmd_weighted_mean(xIV, yIV, 'sideOpto');
e.deltaBiasVertical = conV - incV;
e.deltaMaskVertical = basV - (conV + incV) / 2;
end


function m = bmd_weighted_mean(x, y, mode)
x = x(:); y = y(:);
if isempty(x) || isempty(y)
    m = NaN;
    return
end
switch mode
    case 'sideBaseline'
        w = 10 * ones(size(x));
        w(abs(x) == 0) = 20;
    case 'sideOpto'
        w = 10 * ones(size(x));
    case 'merged'
        w = 20 * ones(size(x));
        w(abs(x) == 0) = 40;
    otherwise
        w = ones(size(x));
end
wsum = sum(w);
if wsum == 0
    m = NaN;
else
    m = sum(w .* y) / wsum;
end
end


function src = infer_baseline_source_from_row(row)
src = '';
try
    dOpto = row.date_opto{1};
    rOpto = row.run_opto{1};
    dBase = row.date_baseline{1};
    rBase = row.run_baseline{1};
    if isempty(dOpto) || isempty(dBase)
        return
    end
    if isnumeric(dOpto) && any(isnan(dOpto(:))); return; end
    if isnumeric(dBase) && any(isnan(dBase(:))); return; end
    if isequal(dOpto, dBase) && isequal(rOpto, rBase)
        src = 'same_block_as_opto';
    else
        src = 'separate_block';
    end
catch
end
end
function [clusterByBlock, significantByBlock] = ...
    get_power_cluster_metadata(mdlStruct, nTotalBlocks)
% Read the power-cluster assignments already saved by psycluster.
%
% This function deliberately does NOT rerun clusterOrderedPowerEffect.
% Therefore, the cluster IDs in MetaTable remain identical to those used
% for the existing power-cluster plots.
%
% Significance is tested against zero in the positive direction:
%
%   H0: delta bias <= 0
%   H1: delta bias > 0
%
% A cluster is marked significant when either:
%   1. its mean is significantly greater than zero by a right-tailed
%      one-sample t-test, or
%   2. its median is significantly greater than zero by a right-tailed
%      signed-rank test.
%
% clusterByBlock:
%   1, 2, 3, ... = saved ordered power-cluster assignment
%   NaN          = no saved assignment
%
% significantByBlock:
%   1   = cluster mean or median is significantly greater than zero
%   0   = cluster was tested but was not significantly greater than zero
%   NaN = no saved cluster assignment

clusterByBlock = NaN(nTotalBlocks, 1);
significantByBlock = NaN(nTotalBlocks, 1);

if ~isstruct(mdlStruct) || isempty(fieldnames(mdlStruct))

    warning('MetaTable:NoMdlStruct', ...
        ['mdlStruct is empty. Saved power-cluster assignments ' ...
         'cannot be loaded.']);

    return

end

[modelData, modelField] = ...
    get_preferred_power_cluster_model(mdlStruct);

if isempty(modelField)

    warning('MetaTable:NoPsychometricModel', ...
        ['Could not find a psychometric model containing ' ...
         'clusterBlocksIdx and deltaBias.']);

    return

end

% Model rows are not necessarily identical to all bitmap block rows.
% Map delta bias back to the original bitmap block IDs.
blockIDs = modelData.clusterBlocksIdx(:);
deltaBias = modelData.deltaBias(:);

nValues = min(numel(blockIDs), numel(deltaBias));

blockIDs = blockIDs(1:nValues);
deltaBias = deltaBias(1:nValues);

validRows = ...
    isfinite(blockIDs) & ...
    blockIDs >= 1 & ...
    blockIDs <= nTotalBlocks & ...
    blockIDs == round(blockIDs);

blockIDs = round(blockIDs(validRows));
deltaBias = deltaBias(validRows);

deltaBiasByBlock = NaN(nTotalBlocks, 1);
deltaBiasByBlock(blockIDs) = deltaBias;

% Load only the cluster assignments already saved by psycluster.
clusterByBlock = get_saved_power_cluster_labels( ...
    mdlStruct, modelField, nTotalBlocks);

if all(~isfinite(clusterByBlock))

    warning('MetaTable:NoSavedPowerClusters', ...
        ['No saved power-cluster assignments were found for %s. ' ...
         'Clusters will not be recalculated. Run psycluster with ' ...
         'cluster saving enabled, then rerun MetaTable.'], ...
        modelField);

    return

end

fprintf(['MetaTable: loaded saved power-cluster assignments ' ...
    'associated with %s.\n'], modelField);

clusterIDs = unique( ...
    clusterByBlock(isfinite(clusterByBlock)))';

for clusterID = clusterIDs

    clusterRows = ...
        clusterByBlock == clusterID & ...
        isfinite(deltaBiasByBlock);

    values = deltaBiasByBlock(clusterRows);
    values = values(isfinite(values));

    if isempty(values)
        continue
    end

    hMean = false;
    hMedian = false;

    pMean = NaN;
    pMedian = NaN;

    % Right-tailed t-test:
    % Is the cluster mean delta bias significantly greater than zero?
    try

        [hMean, pMean] = ttest( ...
            values, ...
            0, ...
            'Tail', 'right', ...
            'Alpha', 0.05);

    catch ME

        warning('MetaTable:ClusterTTest', ...
            'T-test failed for cluster %g: %s', ...
            clusterID, ME.message);

    end

    % Right-tailed signed-rank test:
    % Is the cluster median delta bias significantly greater than zero?
    try

        [pMedian, hMedian] = signrank( ...
            values, ...
            0, ...
            'tail', 'right', ...
            'alpha', 0.05);

    catch ME

        warning('MetaTable:ClusterSignrank', ...
            'Signed-rank test failed for cluster %g: %s', ...
            clusterID, ME.message);

    end

    isSignificant = ...
        logical(hMean) || logical(hMedian);

    significantByBlock( ...
        clusterByBlock == clusterID) = ...
        double(isSignificant);

    fprintf([ ...
        'Power cluster %g: n=%d, mean=%.4g, median=%.4g, ' ...
        'right-tailed mean p=%.4g, ' ...
        'right-tailed median p=%.4g, significant=%d\n'], ...
        clusterID, ...
        numel(values), ...
        mean(values, 'omitnan'), ...
        median(values, 'omitnan'), ...
        pMean, ...
        pMedian, ...
        isSignificant);

end

end


function [modelData, modelField] = ...
    get_preferred_power_cluster_model(mdlStruct)
% Prefer the Weibull-free-all C1 model because that is what psycluster
% currently uses to construct the ordered power clusters.

modelData = struct();
modelField = '';

fieldNames = fieldnames(mdlStruct);

% First pass: explicitly prefer weibullfreeAll C1.
for fieldID = 1:numel(fieldNames)

    fieldName = fieldNames{fieldID};
    candidate = mdlStruct.(fieldName);

    if ~isstruct(candidate) || ...
            ~isfield(candidate, 'clusterBlocksIdx') || ...
            ~isfield(candidate, 'deltaBias')
        continue
    end

    if endsWith(fieldName, 'weibullfreeAllC1')
        modelData = candidate;
        modelField = fieldName;
        return
    end

end

% Fallback: use the first model that contains the required fields.
for fieldID = 1:numel(fieldNames)

    fieldName = fieldNames{fieldID};
    candidate = mdlStruct.(fieldName);

    if isstruct(candidate) && ...
            isfield(candidate, 'clusterBlocksIdx') && ...
            isfield(candidate, 'deltaBias')

        modelData = candidate;
        modelField = fieldName;
        return
    end

end

end


function labelsByBlock = get_saved_power_cluster_labels( ...
    mdlStruct, modelField, nTotalBlocks)

labelsByBlock = NaN(nTotalBlocks, 1);

% Example:
% RweibullfreeAllC1
% becomes
% RweibullfreeAllPowerClusterAggregateLabelsByBlock

aggregateBase = regexprep( ...
    modelField, 'C1$', 'PowerClusterAggregate');

labelsField = [aggregateBase 'LabelsByBlock'];

if isfield(mdlStruct, labelsField)

    values = mdlStruct.(labelsField);

    if isnumeric(values) && numel(values) == nTotalBlocks
        labelsByBlock = values(:);
        return
    end

end

% Fallback to the saved cluster summary.
summaryField = [aggregateBase 'Summary'];

if isfield(mdlStruct, summaryField)

    clusterSummary = mdlStruct.(summaryField);

    if isstruct(clusterSummary)

        for summaryID = 1:numel(clusterSummary)

            if ~isfield(clusterSummary(summaryID), 'clusterID') || ...
                    ~isfield(clusterSummary(summaryID), 'blockIndices')
                continue
            end

            blocks = clean_power_cluster_blocks( ...
                clusterSummary(summaryID).blockIndices, ...
                nTotalBlocks);

            labelsByBlock(blocks) = ...
                clusterSummary(summaryID).clusterID;

        end

        if any(isfinite(labelsByBlock))
            return
        end

    end

end

% Older aggregate structure fallback.
aggregateField = aggregateBase;

if isfield(mdlStruct, aggregateField)

    aggregateData = mdlStruct.(aggregateField);

    if isstruct(aggregateData)

        for aggregateID = 1:numel(aggregateData)

            if ~isfield(aggregateData(aggregateID), 'clusterID') || ...
                    ~isfield(aggregateData(aggregateID), ...
                    'sourceBlockIndices')
                continue
            end

            blocks = clean_power_cluster_blocks( ...
                aggregateData(aggregateID).sourceBlockIndices, ...
                nTotalBlocks);

            labelsByBlock(blocks) = ...
                aggregateData(aggregateID).clusterID;

        end

    end

end

end


function blocks = clean_power_cluster_blocks( ...
    blocks, nTotalBlocks)

blocks = blocks(:);

valid = ...
    isfinite(blocks) & ...
    blocks >= 1 & ...
    blocks <= nTotalBlocks & ...
    blocks == round(blocks);

blocks = round(blocks(valid));

end
