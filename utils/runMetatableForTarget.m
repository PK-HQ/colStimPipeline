function result = runMetatableForTarget(mainPath, animalName, chamberLetter, srcFilename)
% Generic metatable builder for one animal/chamber combination.
%
% Finds the unique statistics<Chamber>-final*.mat source file,
% calls buildBlockMetadata, validates the output, and returns a summary
% struct.
%
% Signature:
%   result = runMetatableForTarget(mainPath, animalName, chamberLetter)
%   result = runMetatableForTarget(mainPath, animalName, chamberLetter, srcFilename)
%
% When srcFilename is provided:
%   - uses that exact file (no wildcard scan)
%   - asserts the file exists
%   - asserts the name contains chamberLetter
%   - asserts the name does not contain 'metatable'
%
% When srcFilename is omitted or empty:
%   - discovers the unique statistics<Chamber>-final*.mat by wildcard,
%     excluding names containing 'metatable', 'temp', or 'backup'
%   - errors if 0 or >1 matches remain
%
% Output struct fields:
%   animal, chamber, sourcePath, nBlocks, matOutputPath, xlsxOutputPath,
%   success

if nargin < 4
    srcFilename = '';
end

result = struct( ...
    'animal',         animalName, ...
    'chamber',        chamberLetter, ...
    'sourcePath',     '', ...
    'nBlocks',        0, ...
    'matOutputPath',  '', ...
    'xlsxOutputPath', '', ...
    'success',        false);

summaryDir = fullfile(mainPath, animalName, 'Meta', 'summary');

%% 1. Source-file selection
if ~isempty(srcFilename)
    % --- Explicit filename provided ---
    if ~isempty(strfind(lower(srcFilename), 'metatable'))
        error('runMetatableForTarget:ExplicitFilenameIsMetatable', ...
            'Explicit source filename must not contain ''metatable'': %s', srcFilename);
    end
    if isempty(strfind(srcFilename, chamberLetter))
        error('runMetatableForTarget:ChamberMismatch', ...
            'Explicit source filename ''%s'' does not contain chamber letter ''%s''.', ...
            srcFilename, chamberLetter);
    end
    srcFile = fullfile(summaryDir, srcFilename);
    if ~exist(srcFile, 'file')
        error('runMetatableForTarget:SourceNotFound', ...
            'Explicit source file not found: %s', srcFile);
    end
    nStr = regexp(srcFilename, ...
        ['^statistics' chamberLetter '-final(\d+)\.mat$'], 'tokens');
    if isempty(nStr) || isempty(nStr{1})
        error('runMetatableForTarget:UnparsedFilename', ...
            'Cannot extract block count from filename: %s', srcFilename);
    end
    nStr = nStr{1}{1};
else
    % --- Wildcard discovery ---
    listing = dir(fullfile(summaryDir, ...
        ['statistics' chamberLetter '-final*.mat']));
    excludeTerms = {'metatable', 'temp', 'backup'};
    keep = true(numel(listing), 1);
    for ki = 1:numel(listing)
        fname = lower(listing(ki).name);
        for ti = 1:numel(excludeTerms)
            if ~isempty(strfind(fname, excludeTerms{ti}))
                keep(ki) = false;
                break;
            end
        end
    end
    listing = listing(keep);
    if numel(listing) == 0
        error('runMetatableForTarget:NoSourceFile', ...
            'No statistics%s-final*.mat found in: %s', chamberLetter, summaryDir);
    end
    if numel(listing) > 1
        names = strjoin({listing.name}, ', ');
        error('runMetatableForTarget:AmbiguousSource', ...
            'Multiple source files in %s: %s', summaryDir, names);
    end
    srcFile = fullfile(summaryDir, listing(1).name);
    tok = regexp(listing(1).name, ...
        ['^statistics' chamberLetter '-final(\d+)\.mat$'], 'tokens');
    if isempty(tok) || isempty(tok{1})
        error('runMetatableForTarget:UnparsedFilename', ...
            'Cannot extract block count from filename: %s', listing(1).name);
    end
    nStr = tok{1}{1};
end

result.sourcePath = srcFile;

%% 2. Output paths
outMat  = fullfile(summaryDir, ['statistics' chamberLetter '-metatable' nStr '.mat']);
outXlsx = fullfile(summaryDir, ['statistics' chamberLetter '-metatable' nStr '.xlsx']);
result.matOutputPath  = outMat;
result.xlsxOutputPath = outXlsx;

fprintf('--- runMetatableForTarget: %s %s ---\n', animalName, chamberLetter);
fprintf('Source:   %s\n', srcFile);
fprintf('MAT out:  %s\n', outMat);
fprintf('XLSX out: %s\n', outXlsx);

%% 3. Load source variables
fprintf('\nLoading source file...\n');
load(srcFile, 'blockData', 'bitmapData', 'behavioralData', ...
     'analysisBlockID', 'datastruct', 'dataTag', 'mdlStruct');

%% 3b. Validate bitmapData.pcadenoisedresp before building
nTotalBlocks = numel(mean(bitmapData.nColumns, 1, 'omitnan'));
if ~isfield(bitmapData, 'pcadenoisedresp')
    error('runMetatableForTarget:MissingPCAField', ...
        ['bitmapData.pcadenoisedresp absent in %s %s source:\n  %s\n' ...
         '  Regenerate using updated bitmap code (getColumnarBitmapV4).'], ...
        animalName, chamberLetter, srcFile);
end
if ndims(bitmapData.pcadenoisedresp) < 4
    error('runMetatableForTarget:PCAFieldWrongDims', ...
        'bitmapData.pcadenoisedresp has %d dims (need >=4) for %s %s:\n  %s', ...
        ndims(bitmapData.pcadenoisedresp), animalName, chamberLetter, srcFile);
end
pcaActualN = size(bitmapData.pcadenoisedresp, 4);
if pcaActualN ~= nTotalBlocks
    error('runMetatableForTarget:PCAFieldSizeMismatch', ...
        ['bitmapData.pcadenoisedresp 4th dim=%d, expected nTotalBlocks=%d\n' ...
         '  animal=%s  chamber=%s\n  source=%s\n' ...
         '  Regenerate using updated bitmap code (getColumnarBitmapV4).'], ...
        pcaActualN, nTotalBlocks, animalName, chamberLetter, srcFile);
end
fprintf('pcadenoisedresp validated: size [%s]\n', num2str(size(bitmapData.pcadenoisedresp)));

%% 4. Build metatable (save handled inside buildBlockMetadata via saveOpts)
fprintf('Building metatable...\n');
saveOpts = struct('matPath', outMat, 'xlsxPath', outXlsx);
MetaTable = buildBlockMetadata(behavioralData, bitmapData, [], [], ...
    blockData, datastruct, analysisBlockID, mdlStruct, saveOpts);

result.nBlocks = height(MetaTable);
tblCols = MetaTable.Properties.VariableNames;
nRows = height(MetaTable);

%% 5-12. Validation suite
fprintf('\n--- Validation: %s %s ---\n', animalName, chamberLetter);

% V1: row count matches N from filename
nExpected = str2double(nStr);
if nRows == nExpected
    fprintf('  V1 PASS: height = %d\n', nRows);
else
    fprintf('  V1 WARN: height=%d, filename implies %s\n', nRows, nStr);
end

% V2: no duplicate blockIDs
nUniqueIDs = numel(unique(MetaTable.blockID));
if nUniqueIDs == nRows
    fprintf('  V2 PASS: unique blockIDs = %d\n', nUniqueIDs);
else
    fprintf('  V2 WARN: unique blockIDs = %d, height = %d\n', nUniqueIDs, nRows);
end

% V3-V8: nonmissing psychometric y fields
yChecks = { ...
    'psy_yBaselineMerged',         'V3'; ...
    'psy_yConOptoMerged',          'V4'; ...
    'psy_yInconOptoMerged',        'V5'; ...
    'psy_yBaselinePreMerge',       'V6'; ...
    'psy_yHorizontalOptoPreMerge', 'V7'; ...
    'psy_yVerticalOptoPreMerge',   'V8'; ...
};
for fi = 1:size(yChecks, 1)
    col = yChecks{fi, 1};
    tag = yChecks{fi, 2};
    if ~any(strcmp(tblCols, col))
        fprintf('  %s MISS: %s not in table\n', tag, col);
        continue;
    end
    n = sum(~cellfun(@isempty, MetaTable.(col)));
    if n == nRows
        fprintf('  %s PASS: %s = %d/%d\n', tag, col, n, nRows);
    else
        fprintf('  %s WARN: %s = %d/%d\n', tag, col, n, nRows);
    end
end

% V9: delta completeness
deltaFields = { ...
    'psy_deltaBias', 'psy_deltaMask', ...
    'psy_deltaBiasMerged', 'psy_deltaMaskMerged', ...
    'psy_deltaBiasHorizontal', 'psy_deltaMaskHorizontal', ...
    'psy_deltaBiasVertical', 'psy_deltaMaskVertical' ...
};
nDeltaOk = 0; nDeltaTotal = 0;
for di = 1:numel(deltaFields)
    col = deltaFields{di};
    if any(strcmp(tblCols, col))
        nDeltaOk    = nDeltaOk + sum(isfinite(MetaTable.(col)));
        nDeltaTotal = nDeltaTotal + nRows;
    end
end
if nDeltaOk == nDeltaTotal && nDeltaTotal > 0
    fprintf('  V9 PASS: all %d delta values finite\n', nDeltaTotal);
else
    fprintf('  V9 WARN: %d/%d delta values finite\n', nDeltaOk, nDeltaTotal);
end

% V10: expected columns present
expectedCols = { ...
    'blockID', 'psy_modelField', ...
    'psy_xBaselineMerged', 'psy_yBaselineMerged', ...
    'psy_xConOptoMerged', 'psy_yConOptoMerged', ...
    'psy_xInconOptoMerged', 'psy_yInconOptoMerged', ...
    'psy_xBaselinePreMerge', 'psy_yBaselinePreMerge', 'psy_xBaselinePreMergeOrt', ...
    'psy_nTrialsBaselinePreMerge', 'psy_nTrialsBaselineMerged', ...
    'psy_xHorizontalOptoPreMerge', 'psy_yHorizontalOptoPreMerge', ...
    'psy_xHorizontalOptoPreMergeOrt', 'psy_congruencyHorizontalOptoPreMerge', ...
    'psy_xVerticalOptoPreMerge', 'psy_yVerticalOptoPreMerge', ...
    'psy_xVerticalOptoPreMergeOrt', 'psy_congruencyVerticalOptoPreMerge', ...
    'psy_deltaBias', 'psy_deltaBiasMerged', ...
    'psy_deltaBiasHorizontal', 'psy_deltaBiasVertical', ...
    'opto_PCAdenoisedResp', ...
    'bmp_horizontalCamSpace', 'bmp_verticalCamSpace', ...
    'session_hasVisFPS', 'session_baselineSource', ...
};
nMissing = 0;
for ci = 1:numel(expectedCols)
    if ~any(strcmp(tblCols, expectedCols{ci}))
        fprintf('  V10 MISS: %s\n', expectedCols{ci});
        nMissing = nMissing + 1;
    end
end
if nMissing == 0
    fprintf('  V10 PASS: all %d expected columns present\n', numel(expectedCols));
end

% V10b: opto_PCAdenoisedResp must immediately precede opto_bitmapCamSpace
names_mt  = MetaTable.Properties.VariableNames;
pcaColIdx = find(strcmp(names_mt, 'opto_PCAdenoisedResp'));
camColIdx = find(strcmp(names_mt, 'opto_bitmapCamSpace'));
if isempty(pcaColIdx)
    fprintf('  V10b FAIL: opto_PCAdenoisedResp missing from MetaTable\n');
elseif isempty(camColIdx)
    fprintf('  V10b FAIL: opto_bitmapCamSpace missing from MetaTable\n');
elseif pcaColIdx + 1 == camColIdx
    fprintf('  V10b PASS: opto_PCAdenoisedResp col %d immediately before opto_bitmapCamSpace col %d\n', ...
        pcaColIdx, camColIdx);
else
    fprintf('  V10b FAIL: opto_PCAdenoisedResp col %d, opto_bitmapCamSpace col %d (not adjacent)\n', ...
        pcaColIdx, camColIdx);
end

% PCAdenoisedResp blockwise audit
fprintf('\n--- PCAdenoisedResp audit: %s %s ---\n', animalName, chamberLetter);
pcaSrcSz = size(bitmapData.pcadenoisedresp);
fprintf('  source size: [%s]\n', num2str(pcaSrcSz));
if any(strcmp(tblCols, 'opto_PCAdenoisedResp')) && nRows > 0
    firstBID = MetaTable.blockID(1);
    lastBID  = MetaTable.blockID(end);
    srcFirst = bitmapData.pcadenoisedresp(:,:,:,firstBID);
    srcLast  = bitmapData.pcadenoisedresp(:,:,:,lastBID);
    tblFirst = MetaTable.opto_PCAdenoisedResp{1};
    tblLast  = MetaTable.opto_PCAdenoisedResp{end};
    fprintf('  first-block checksum: src=%.6g  tbl=%.6g\n', ...
        sum(double(srcFirst(:)), 'omitnan'), sum(double(tblFirst(:)), 'omitnan'));
    fprintf('  last-block  checksum: src=%.6g  tbl=%.6g\n', ...
        sum(double(srcLast(:)), 'omitnan'), sum(double(tblLast(:)), 'omitnan'));
    pcaOk = 0;
    for bi = 1:nRows
        bID      = MetaTable.blockID(bi);
        srcSlice = bitmapData.pcadenoisedresp(:,:,:,bID);
        tblSlice = MetaTable.opto_PCAdenoisedResp{bi};
        if isequaln(size(tblSlice), size(srcSlice)) && isequaln(tblSlice, srcSlice)
            pcaOk = pcaOk + 1;
        else
            fprintf('  MISMATCH block %d: src size [%s]  tbl size [%s]\n', ...
                bID, num2str(size(srcSlice)), num2str(size(tblSlice)));
        end
    end
    fprintf('  PCAdenoisedResp populated: %d/%d blocks\n', pcaOk, nRows);
    if ~isempty(pcaColIdx) && ~isempty(camColIdx) && pcaColIdx + 1 == camColIdx
        fprintf('  column position: immediately before opto_bitmapCamSpace\n');
    end
else
    fprintf('  opto_PCAdenoisedResp not found in MetaTable\n');
end

% V11: modelField distribution
fprintf('  V11 psy_modelField distribution:\n');
mf = MetaTable.psy_modelField;
uMF = unique(mf);
for i = 1:numel(uMF)
    fprintf('    %-52s  count=%d\n', uMF{i}, sum(strcmp(mf, uMF{i})));
end

% V12: C1 overlap validation using sorted x/y pair comparison
fprintf('  V12 C1 overlap validation:\n');
c1Field = rmt_find_c1_field(mdlStruct, chamberLetter);
if isempty(c1Field)
    fprintf('    No *C1 field with clusterBlocksIdx for chamber %s — skipping.\n', chamberLetter);
else
    fprintf('    Model field: %s\n', c1Field);
    md = mdlStruct.(c1Field);
    c1BlockIDs = md.clusterBlocksIdx(:);
    xRaw = behavioralData.gaborContrasts;
    yRaw = behavioralData.percentageCorrect;

    % x/y pairs: sort by x, apply same permutation to y, then compare.
    % This tolerates any storage-order difference while still detecting
    % genuine value differences.
    xyPairDefs = { ...
        'psy_xBaselineMerged',          'psy_yBaselineMerged',          'xBaselineMerged',          'yBaselineMerged'; ...
        'psy_xConOptoMerged',           'psy_yConOptoMerged',           'xConOptoMerged',           'yConOptoMerged'; ...
        'psy_xInconOptoMerged',         'psy_yInconOptoMerged',         'xInconOptoMerged',         'yInconOptoMerged'; ...
        'psy_xBaselinePreMerge',        'psy_yBaselinePreMerge',        'xBaselinePreMerge',        'yBaselinePreMerge'; ...
        'psy_xHorizontalOptoPreMerge',  'psy_yHorizontalOptoPreMerge',  'xHorizontalOptoPreMerge',  'yHorizontalOptoPreMerge'; ...
        'psy_xVerticalOptoPreMerge',    'psy_yVerticalOptoPreMerge',    'xVerticalOptoPreMerge',    'yVerticalOptoPreMerge'; ...
    };
    % Single fields: sort independently
    singleFieldDefs = { ...
        'psy_xBaselinePreMergeOrt',             'xBaselinePreMergeOrt'; ...
        'psy_xHorizontalOptoPreMergeOrt',       'xHorizontalOptoPreMergeOrt'; ...
        'psy_congruencyHorizontalOptoPreMerge', 'congruencyHorizontalOptoPreMerge'; ...
        'psy_xVerticalOptoPreMergeOrt',         'xVerticalOptoPreMergeOrt'; ...
        'psy_congruencyVerticalOptoPreMerge',   'congruencyVerticalOptoPreMerge'; ...
    };

    TOL = 1e-9;
    nFields = size(xyPairDefs, 1) * 2 + size(singleFieldDefs, 1);
    mismatchCounts = zeros(nFields, 1);
    mismatchLabels = cell(nFields, 1);
    fi = 0;
    for pi = 1:size(xyPairDefs, 1)
        fi = fi + 1; mismatchLabels{fi} = xyPairDefs{pi, 1};
        fi = fi + 1; mismatchLabels{fi} = xyPairDefs{pi, 2};
    end
    for si = 1:size(singleFieldDefs, 1)
        fi = fi + 1; mismatchLabels{fi} = singleFieldDefs{si, 1};
    end

    nTested = 0;
    firstMismatchBlockID = [];
    firstMismatchField   = '';
    firstMismatchDetail  = struct();

    for bi = 1:numel(c1BlockIDs)
        bID = c1BlockIDs(bi);
        rowIdx = find(MetaTable.blockID == bID, 1);
        if isempty(rowIdx)
            fprintf('    WARNING: C1 block %d not in MetaTable\n', bID);
            continue;
        end
        eRaw = rmt_compute_raw_psy(xRaw, yRaw, bID);
        nTested = nTested + 1;

        fi = 0;
        blockHadMismatch = false;

        % --- x/y pairs ---
        for pi = 1:size(xyPairDefs, 1)
            xTblCol = xyPairDefs{pi, 1}; yTblCol = xyPairDefs{pi, 2};
            xRawFld = xyPairDefs{pi, 3}; yRawFld = xyPairDefs{pi, 4};

            fi_x = fi + 1; fi_y = fi + 2; fi = fi + 2;

            if ~any(strcmp(tblCols, xTblCol)) || ~any(strcmp(tblCols, yTblCol))
                continue;
            end
            xT = MetaTable.(xTblCol){rowIdx}; yT = MetaTable.(yTblCol){rowIdx};
            xRv = eRaw.(xRawFld);          yRv = eRaw.(yRawFld);

            % Strip NaN
            validT = ~isnan(xT) & ~isnan(yT);
            xT = xT(validT); yT = yT(validT);
            validR = ~isnan(xRv) & ~isnan(yRv);
            xRv = xRv(validR); yRv = yRv(validR);

            % Sort by x, apply same perm to y
            [xT, ixT] = sort(xT(:)); yT = yT(ixT);
            [xRv, ixR] = sort(xRv(:)); yRv = yRv(ixR);

            xMismatch = numel(xT) ~= numel(xRv) || any(abs(xT - xRv) > TOL);
            yMismatch = numel(yT) ~= numel(yRv) || any(abs(yT(:) - yRv(:)) > TOL);
            if xMismatch; mismatchCounts(fi_x) = mismatchCounts(fi_x) + 1; end
            if yMismatch; mismatchCounts(fi_y) = mismatchCounts(fi_y) + 1; end

            if (xMismatch || yMismatch) && isempty(firstMismatchBlockID)
                blockHadMismatch = true;
                if xMismatch
                    firstMismatchBlockID = bID;
                    firstMismatchField   = xTblCol;
                    d.blockID     = bID;
                    d.modelField  = MetaTable.psy_modelField{rowIdx};
                    d.tblCol      = xTblCol;
                    d.rawFld      = xRawFld;
                    d.xTbl_sorted = xT(:)';
                    d.xRaw_sorted = xRv(:)';
                    d.sizeTbl     = numel(xT);
                    d.sizeRaw     = numel(xRv);
                    d.isRowTbl    = isrow(MetaTable.(xTblCol){rowIdx});
                    d.isRowRaw    = isrow(eRaw.(xRawFld));
                    if numel(xT) == numel(xRv)
                        d.maxAbsDiff = max(abs(xT(:) - xRv(:)));
                        d.setsEqual  = all(ismember(round(xT*1e10)/1e10, round(xRv*1e10)/1e10));
                    else
                        d.maxAbsDiff = Inf;
                        d.setsEqual  = false;
                    end
                    firstMismatchDetail = d;
                end
            end
        end

        % --- single fields ---
        for si = 1:size(singleFieldDefs, 1)
            fi = fi + 1;
            tblCol = singleFieldDefs{si, 1};
            rawFld = singleFieldDefs{si, 2};
            if ~any(strcmp(tblCols, tblCol)); continue; end
            vT = MetaTable.(tblCol){rowIdx};
            vR = eRaw.(rawFld);
            vT = sort(vT(~isnan(vT(:))'));
            vR = sort(vR(~isnan(vR(:))'));
            if numel(vT) ~= numel(vR) || any(abs(vT(:) - vR(:)) > TOL)
                mismatchCounts(fi) = mismatchCounts(fi) + 1;
                if ~blockHadMismatch && isempty(firstMismatchBlockID)
                    firstMismatchBlockID = bID;
                    firstMismatchField   = tblCol;
                end
            end
        end
    end

    fprintf('    Blocks tested: %d / %d\n', nTested, numel(c1BlockIDs));

    % Print diagnostic for first mismatching block
    if ~isempty(firstMismatchBlockID) && isstruct(firstMismatchDetail) && isfield(firstMismatchDetail, 'blockID')
        d = firstMismatchDetail;
        fprintf('\n    --- First mismatch diagnostic ---\n');
        fprintf('    blockID          : %d\n',   d.blockID);
        fprintf('    psy_modelField   : %s\n',   d.modelField);
        fprintf('    field            : %s / %s\n', d.tblCol, d.rawFld);
        fprintf('    size (table/raw) : %d / %d\n', d.sizeTbl, d.sizeRaw);
        fprintf('    isrow (tbl/raw)  : %d / %d\n', d.isRowTbl, d.isRowRaw);
        if numel(d.xTbl_sorted) > 0
            fprintf('    table sorted     : [%s]\n', num2str(d.xTbl_sorted, '%g '));
        else
            fprintf('    table sorted     : []\n');
        end
        if numel(d.xRaw_sorted) > 0
            fprintf('    raw sorted       : [%s]\n', num2str(d.xRaw_sorted, '%g '));
        else
            fprintf('    raw sorted       : []\n');
        end
        fprintf('    sets equal       : %d\n',   d.setsEqual);
        fprintf('    max |diff|       : %g\n',   d.maxAbsDiff);
        fprintf('    ---\n\n');
    end

    anyMismatch = any(mismatchCounts > 0);
    if ~anyMismatch
        fprintf('    V12 PASS: all fields match across %d C1 blocks.\n', nTested);
    else
        for fi = 1:numel(mismatchLabels)
            if mismatchCounts(fi) > 0
                fprintf('    MISMATCH: %-52s  %d blocks\n', mismatchLabels{fi}, mismatchCounts(fi));
            end
        end
        fprintf('    V12 WARN: mismatches found — inspect diagnostic above.\n');
    end
end

% Final delta column order assertion (tail 8)
% Expected order as produced by buildBlockMetadata:
requiredTail = { ...
    'psy_deltaMask',          'psy_deltaBias', ...
    'psy_deltaBiasHorizontal','psy_deltaMaskHorizontal', ...
    'psy_deltaBiasVertical',  'psy_deltaMaskVertical', ...
    'psy_deltaBiasMerged',    'psy_deltaMaskMerged' ...
};
fprintf('\n--- Final delta column order ---\n');
nCols = numel(tblCols);
if nCols >= numel(requiredTail)
    tailCols = tblCols(nCols - numel(requiredTail) + 1 : nCols);
    if isequal(tailCols, requiredTail)
        fprintf('  PASS: final eight psychometric delta columns are in the expected order.\n');
    else
        fprintf('  WARN: tail columns differ.\n');
        fprintf('  Expected: %s\n', strjoin(requiredTail, ', '));
        fprintf('  Actual:   %s\n', strjoin(tailCols, ', '));
    end
else
    fprintf('  WARN: table has only %d columns (need >= %d).\n', ...
        nCols, numel(requiredTail));
end

result.success = true;
fprintf('\nDone: %s %s -> %s\n', animalName, chamberLetter, outMat);
end


% =========================================================================
%  Private helpers
% =========================================================================

function c1Field = rmt_find_c1_field(mdlStruct, chamberLetter)
% Find the first mdlStruct field matching <Chamber>.*C1 that has
% clusterBlocksIdx.
c1Field = '';
if ~isstruct(mdlStruct); return; end
fnames = fieldnames(mdlStruct);
pat = ['^' chamberLetter '.*C1$'];
for fi = 1:numel(fnames)
    fn = fnames{fi};
    if ~isempty(regexp(fn, pat, 'once')) && ...
            isfield(mdlStruct.(fn), 'clusterBlocksIdx')
        c1Field = fn;
        return;
    end
end
end


function e = rmt_compute_raw_psy(xBlocks, yBlocks, blockID)
% Replicate processConditionsBlocks (fitPsyMLE2) for one block.
% Returns struct with same field names as the psyFull entries.
e = struct();
e.xBaselineMerged = []; e.yBaselineMerged = [];
e.xConOptoMerged = []; e.yConOptoMerged = [];
e.xInconOptoMerged = []; e.yInconOptoMerged = [];
e.xBaselinePreMerge = []; e.yBaselinePreMerge = [];
e.xBaselinePreMergeOrt = [];
e.xHorizontalOptoPreMerge = []; e.yHorizontalOptoPreMerge = [];
e.xHorizontalOptoPreMergeOrt = [];
e.congruencyHorizontalOptoPreMerge = [];
e.xVerticalOptoPreMerge = []; e.yVerticalOptoPreMerge = [];
e.xVerticalOptoPreMergeOrt = [];
e.congruencyVerticalOptoPreMerge = [];

xBase = rmt_rmnan(squeeze(xBlocks(1, :, blockID)));
yBase = rmt_rmnan(squeeze(yBlocks(1, :, blockID)));
xH    = rmt_rmnan(squeeze(xBlocks(2, :, blockID)));
yH    = rmt_rmnan(squeeze(yBlocks(2, :, blockID)));
xV    = rmt_rmnan(squeeze(xBlocks(3, :, blockID)));
yV    = rmt_rmnan(squeeze(yBlocks(3, :, blockID)));

%% Baseline
nBase = numel(xBase);
if nBase > 0 && mod(nBase, 2) == 0
    tag = rmt_make_vtag(nBase);
    [xS, si] = sort(xBase); yS = yBase(si); tS = tag(si);
    [xPre, yPre, tPre] = rmt_pre_merge(xS, yS, tS, true);
    e.xBaselinePreMerge    = xPre(:)';
    e.yBaselinePreMerge    = yPre(:)';
    e.xBaselinePreMergeOrt = tPre(:)';
    n = numel(xS);
    % fliplr(-xS(1:n/2)): negate the negative half to get positive contrasts,
    % then flip to ascending order; mean with the positive half = absolute contrast.
    e.xBaselineMerged = mean([fliplr(-xS(1:n/2)); xS(n/2+1:end)]);
    e.yBaselineMerged = mean([fliplr(100 - yS(1:n/2)); yS(n/2+1:end)]);
end

%% Horizontal opto
nH = numel(xH);
if nH > 0 && mod(nH, 2) == 0
    tag = rmt_make_vtag(nH);
    [xS, si] = sort(xH); yS = yH(si); tS = tag(si);
    [xPre, yPre, tPre] = rmt_pre_merge(xS, yS, tS, false);
    congrH = NaN(size(tPre));
    congrH(tPre == 0)  =  1;
    congrH(tPre == 90) = -1;
    e.xHorizontalOptoPreMerge          = xPre(:)';
    e.yHorizontalOptoPreMerge          = yPre(:)';
    e.xHorizontalOptoPreMergeOrt       = tPre(:)';
    e.congruencyHorizontalOptoPreMerge = congrH(:)';
end

%% Vertical opto
nV = numel(xV);
if nV > 0 && mod(nV, 2) == 0
    tag = rmt_make_vtag(nV);
    [xS, si] = sort(xV); yS = yV(si); tS = tag(si);
    [xPre, yPre, tPre] = rmt_pre_merge(xS, yS, tS, false);
    congrV = NaN(size(tPre));
    congrV(tPre == 0)  = -1;
    congrV(tPre == 90) =  1;
    e.xVerticalOptoPreMerge          = xPre(:)';
    e.yVerticalOptoPreMerge          = yPre(:)';
    e.xVerticalOptoPreMergeOrt       = tPre(:)';
    e.congruencyVerticalOptoPreMerge = congrV(:)';
end

%% Merged con/incon
if nH > 0 && mod(nH,2)==0 && nV > 0 && mod(nV,2)==0
    [xHS, si2] = sort(xH); yHS = yH(si2);
    [xVS, si3] = sort(xV); yVS = yV(si3);
    cNeg = 1:nH/2; cPos = nH/2+1:nH;
    yHC = yHS; yHC(cNeg) = 100 - yHC(cNeg);
    yVC = yVS; yVC(cNeg) = 100 - yVC(cNeg);
    e.xConOptoMerged   = mean([-fliplr(xHS(cNeg)); xVS(cPos)])';
    e.yConOptoMerged   = mean([fliplr(yHC(cNeg));  yVC(cPos)])';
    e.xInconOptoMerged = mean([xHS(cPos); -fliplr(xVS(cNeg))])';
    e.yInconOptoMerged = mean([yHC(cPos);  fliplr(yVC(cNeg))])';
end
end


function v = rmt_rmnan(v)
v = v(~isnan(v(:))');
end


function tag = rmt_make_vtag(nVal)
tag = NaN(1, nVal);
tag(1:nVal/2)     = 0;
tag(nVal/2+1:end) = 90;
end


function [xOut, yOut, tagOut] = rmt_pre_merge(xIn, yIn, tagIn, mergeDupZeros)
xIn = xIn(:)'; yIn = yIn(:)'; tagIn = tagIn(:)';
valid = ~isnan(xIn) & ~isnan(yIn) & ~isnan(tagIn);
xIn = xIn(valid); yIn = yIn(valid); tagIn = tagIn(valid);
yCorrect = yIn;
yCorrect(tagIn == 0) = 100 - yCorrect(tagIn == 0);
if mergeDupZeros
    [xOut, yOut, tagOut] = rmt_merge_dup_x(xIn, yCorrect, tagIn);
else
    xOut = xIn; yOut = yCorrect; tagOut = tagIn;
end
end


function [xOut, yOut, tagOut] = rmt_merge_dup_x(xIn, yIn, tagIn)
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
