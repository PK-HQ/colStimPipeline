function generateAnimalMetaTables()
% Discover Chip and Pepper final-statistics MAT files, build a full
% per-block metatable for each (no column filtering), and save MAT + XLSX.
%
% Output naming: statistics{chamber}-metatable{N}.mat / .xlsx
% Saved to the same Meta/summary folder as the source file.

mainPath = 'Y:/';
sources  = discover_final_statistics(mainPath);

if isempty(sources)
    error('generateAnimalMetaTables:NoSources', ...
        'No statistics*-final*.mat files found under Chip or Pepper Meta/summary.');
end

for srcIdx = 1:numel(sources)
    src = sources(srcIdx);
    fprintf('\n=== Processing: %s ===\n', src.matFile);

    % ---- discover which variables exist in this file ----
    w = whos('-file', src.matFile);
    varNamesInFile = {w.name};

    reqVars = {'behavioralData', 'bitmapData', 'datastruct', ...
               'analysisBlockID', 'mdlStruct'};
    missing = setdiff(reqVars, varNamesInFile);
    if ~isempty(missing)
        warning('generateAnimalMetaTables:MissingVar', ...
            'Source file missing required variables: %s  -- skipping.', ...
            strjoin(missing, ', '));
        continue;
    end

    loadVars = reqVars;
    if ismember('blockData', varNamesInFile)
        loadVars{end+1} = 'blockData';
    else
        fprintf('  NOTE: blockData absent from source — provenance dates will be empty.\n');
    end

    loaded = load(src.matFile, loadVars{:});

    blockData = struct();
    if isfield(loaded, 'blockData')
        blockData = loaded.blockData;
    end

    % ---- build metatable via the extended buildBlockMetadata ----
    % Pass saveOpts so buildBlockMetadata writes directly to the desired paths.
    matOut  = fullfile(src.outDir, [src.stem '.mat']);
    xlsxOut = fullfile(src.outDir, [src.stem '.xlsx']);
    saveOpts = struct('matPath', matOut, 'xlsxPath', xlsxOut);

    MetaTable = buildBlockMetadata( ...
        loaded.behavioralData, loaded.bitmapData, ...
        [], [], ...                      % no column filtering
        blockData, loaded.datastruct, loaded.analysisBlockID, ...
        loaded.mdlStruct, saveOpts);

    % ---- validation ----
    fprintf('\n--- Validation: %s ---\n', src.stem);
    nRows = height(MetaTable);
    nCols = width(MetaTable);
    fprintf('Rows: %d  |  Columns: %d\n', nRows, nCols);

    % Row count vs analysisBlockID
    nAnalyzed = numel(loaded.analysisBlockID);
    if nRows == nAnalyzed
        fprintf('Row count matches analysisBlockID length (%d). OK\n', nAnalyzed);
    else
        fprintf('WARNING: Row count %d != analysisBlockID length %d\n', ...
            nRows, nAnalyzed);
    end

    % Required psychometric columns
    reqPsy = {'psy_xBaselinePreMerge', 'psy_yBaselinePreMerge', ...
              'psy_xHorizontalOptoPreMerge', 'psy_yHorizontalOptoPreMerge', ...
              'psy_xVerticalOptoPreMerge',   'psy_yVerticalOptoPreMerge', ...
              'psy_xBaseline', 'psy_yBaseline', ...
              'psy_xConOpto',  'psy_yConOpto', ...
              'psy_xInconOpto','psy_yInconOpto', ...
              'psy_deltaBias', 'psy_deltaMask', ...
              'psy_deltaBiasHorizontal','psy_deltaMaskHorizontal', ...
              'psy_deltaBiasVertical',  'psy_deltaMaskVertical', ...
              'psy_deltaBiasMerged',    'psy_deltaMaskMerged', ...
              'psy_combinedBL', 'session_baselineSource', ...
              'session_hasVisFPS'};
    existingCols = MetaTable.Properties.VariableNames;
    missingPsy   = setdiff(reqPsy, existingCols);
    if isempty(missingPsy)
        fprintf('All required psychometric columns present.\n');
    else
        fprintf('MISSING psychometric columns: %s\n', strjoin(missingPsy, ', '));
    end

    % MAT bitmap check: first non-empty horizontal bitmap must be numeric
    if ismember('bmp_horizontalCamSpace', existingCols)
        bmpOK = false;
        for rr = 1:nRows
            bmpVal = MetaTable.bmp_horizontalCamSpace{rr};
            if isnumeric(bmpVal) && ~isempty(bmpVal)
                fprintf('bmp_horizontalCamSpace{%d}: %dx%d %s (numeric OK)\n', ...
                    rr, size(bmpVal,1), size(bmpVal,2), class(bmpVal));
                bmpOK = true;
                break;
            end
        end
        if ~bmpOK
            fprintf('WARNING: bmp_horizontalCamSpace is empty or non-numeric for all rows.\n');
        end
    end

    % Excel file check: bitmap cell must be compact summary string
    if exist(xlsxOut, 'file')
        try
            xlData = readcell(xlsxOut, 'Sheet', 'metadata');
            bmpColInXL = find(strcmp(xlData(1,:), 'bmp.horizontalCamSpace'), 1);
            if isempty(bmpColInXL)
                bmpColInXL = find(strcmp(existingCols, 'bmp_horizontalCamSpace'));
                if ~isempty(bmpColInXL)
                    bmpColInXL = bmpColInXL(1);
                end
            end
            if ~isempty(bmpColInXL) && size(xlData,1) >= 2
                cellVal = xlData{2, bmpColInXL};
                if ischar(cellVal) && ~isempty(regexp(cellVal, '^\d+x\d+', 'once'))
                    fprintf('Excel bitmap cell: "%s" (compact summary OK)\n', cellVal);
                else
                    fprintf('Excel bitmap cell class: %s, value: %s\n', ...
                        class(cellVal), char(string(cellVal)));
                end
            end
        catch ME
            fprintf('Excel read check skipped: %s\n', ME.message);
        end
    end

    % Sample rows: x/y length match for baseline and contrast-correct
    sampleRows = unique([1, round(nRows/2), nRows]);
    fprintf('\nSample row xBaseline/yBaseline lengths (must match):\n');
    for rr = sampleRows
        xbl = MetaTable.psy_xBaseline{rr};
        ybl = MetaTable.psy_yBaseline{rr};
        match = numel(xbl) == numel(ybl);
        fprintf('  Row %2d: len=%d/%d  match=%d\n', rr, numel(xbl), numel(ybl), match);
    end
    fprintf('Sample row xHorizOptoPreMerge/yHorizOptoPreMerge lengths:\n');
    for rr = sampleRows
        xh = MetaTable.psy_xHorizontalOptoPreMerge{rr};
        yh = MetaTable.psy_yHorizontalOptoPreMerge{rr};
        match = numel(xh) == numel(yh);
        fprintf('  Row %2d: len=%d/%d  match=%d\n', rr, numel(xh), numel(yh), match);
    end

    % Early combined-BL rows
    if ismember('session_baselineSource', existingCols)
        nCombined  = sum(strcmp(MetaTable.session_baselineSource, 'same_block_as_opto'));
        nSeparate  = sum(strcmp(MetaTable.session_baselineSource, 'separate_block'));
        nUnknown   = nRows - nCombined - nSeparate;
        fprintf('\nBaseline source: same_block=%d  separate=%d  unknown=%d\n', ...
            nCombined, nSeparate, nUnknown);
    end

    % Missing footprint sessions must stay missing (no run=0 invented)
    if ismember('run_visfootprint', existingCols) && ismember('session_hasVisFPS', existingCols)
        nHasFPS  = sum(MetaTable.session_hasVisFPS);
        nNoFPS   = nRows - nHasFPS;
        fprintf('Visual footprint sessions: present=%d  absent=%d\n', nHasFPS, nNoFPS);
        % Check no spurious run=0
        fpRuns = MetaTable.run_visfootprint;
        nRun0  = 0;
        for rr = 1:nRows
            rv = fpRuns{rr};
            if isnumeric(rv) && ~isempty(rv) && any(rv(:) == 0)
                nRun0 = nRun0 + 1;
            end
        end
        if nRun0 > 0
            fprintf('WARNING: %d rows have run_visfootprint==0 (possibly invented).\n', nRun0);
        else
            fprintf('No spurious run_visfootprint=0 entries. OK\n');
        end
    end

    % clusterBlocksIdx coverage: count blocks covered by any model field
    mdlFields = fieldnames(loaded.mdlStruct);
    coveredBlocks = [];
    for fIdx = 1:numel(mdlFields)
        md = loaded.mdlStruct.(mdlFields{fIdx});
        if isstruct(md) && isfield(md, 'clusterBlocksIdx')
            coveredBlocks = union(coveredBlocks, md.clusterBlocksIdx(:)');
        end
    end
    coveredInRange = sum(coveredBlocks >= 1 & coveredBlocks <= nAnalyzed);
    fprintf('clusterBlocksIdx coverage: %d/%d analysis blocks covered.\n', ...
        coveredInRange, nAnalyzed);

    % Source file integrity: whos must match original
    w2 = whos('-file', src.matFile);
    origNames = sort({w.name});
    nowNames  = sort({w2.name});
    if isequal(origNames, nowNames)
        fprintf('Source file integrity: unchanged. OK\n');
    else
        fprintf('WARNING: Source file variable list changed!\n');
    end

    % Output paths
    fprintf('\nOutput MAT:  %s\n', matOut);
    fprintf('Output XLSX: %s\n', xlsxOut);
end

fprintf('\n=== generateAnimalMetaTables done. ===\n');
end


function sources = discover_final_statistics(mainPath)
% Return array of structs: .matFile, .outDir, .stem
animals = {'Chip', 'Pepper'};
sources = struct('matFile', {}, 'outDir', {}, 'stem', {});

for anIdx = 1:numel(animals)
    animal   = animals{anIdx};
    summDir  = fullfile(mainPath, animal, 'Meta', 'summary');
    if ~exist(summDir, 'dir')
        continue;
    end
    listing = dir(fullfile(summDir, 'statistics*-final*.mat'));
    for liIdx = 1:numel(listing)
        fname = listing(liIdx).name;
        tok = regexp(fname, '^(statistics[A-Z])-final(\d+)\.mat$', 'tokens', 'once');
        if isempty(tok)
            continue;
        end
        stem = [tok{1} '-metatable' tok{2}];
        s.matFile = fullfile(summDir, fname);
        s.outDir  = summDir;
        s.stem    = stem;
        sources(end+1) = s; %#ok<AGROW>
    end
end
end
