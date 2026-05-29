function metaT = buildBlockMetadata(behavioralData, bitmapData, columnsDesired, columnsSpread)

meanCol = mean(bitmapData.nColumns, 1, 'omitnan');
blockIdx = meanCol > columnsDesired-columnsSpread & ...
           meanCol < columnsDesired+columnsSpread;

selectedBlocks = find(blockIdx);
nBlocks = numel(selectedBlocks);

rows = struct([]);

for ii = 1:nBlocks
    blockID = selectedBlocks(ii);

    S = behavioralData.optoTS(blockID).Header.ConditionParams.Stimulus;
    P = behavioralData.optoTS(blockID).Header.ConditionParams.Stimulus_Position;

    row = struct();

    row.blockID = blockID;
    row.meanColumns = meanCol(blockID);

    % visual / behavioral metadata
    row.GaborContrast_pc = getfield_safe(S, 'GaborContrast__pc');
    row.GaborSize_deg    = getfield_safe(S, 'GaborSize__deg');
    row.GaborSF_cpd      = getfield_safe(S, 'GaborSF_cpd');
    row.GaborOrt_deg     = getfield_safe(S, 'GaborOrt__deg');
    row.GaborX_deg       = getfield_safe(P, 'X__deg');
    row.GaborY_deg       = getfield_safe(P, 'Y__deg');

    % bitmap metadata summaries
    row.gridSize_mean = mean(bitmapData.gridSize(:,blockID), 'omitnan');
    row.nColumns_mean = mean(bitmapData.nColumns(:,blockID), 'omitnan');
    row.nColumns_min  = min(bitmapData.nColumns(:,blockID), [], 'omitnan');
    row.nColumns_max  = max(bitmapData.nColumns(:,blockID), [], 'omitnan');

    row.sensitivity_mean = mean(bitmapData.sensitivity(:,blockID), 'omitnan');
    row.adaptthresh_mean = mean(bitmapData.adaptthresh(:,blockID), 'omitnan');
    row.pixelsON_mean    = mean(bitmapData.pixelsON(:,blockID), 'omitnan');

    row.energy_mean = mean(bitmapData.energy(:,:,blockID), 'all', 'omitnan');
    row.energy_max  = max(bitmapData.energy(:,:,blockID), [], 'all');

    row.powerdensity_mean = mean(bitmapData.powerdensity(:,:,blockID), 'all', 'omitnan');
    row.powerdensity_max  = max(bitmapData.powerdensity(:,:,blockID), [], 'all');

    row.timeONPercent_mean = mean(bitmapData.timeONPercent(:,:,blockID), 'all', 'omitnan');

    % optional: filename for heavy bitmap arrays
    row.bitmapFile = sprintf('bitmap_block_%03d.mat', blockID);

    bitmapProjSpace = bitmapData.columnarbitmapTFprojspace(:,:,:,blockID);
    bitmapCamSpace  = bitmapData.columnarbitmapTFcamspace(:,:,:,blockID);
    save(row.bitmapFile, 'bitmapProjSpace', 'bitmapCamSpace');

    rows(ii) = row;
end

metaT = struct2table(rows);

end


function val = getfield_safe(S, fieldName)

if isfield(S, fieldName)
    val = S.(fieldName);
    if ~isscalar(val)
        val = val(:)';
    end
else
    val = NaN;
end

end