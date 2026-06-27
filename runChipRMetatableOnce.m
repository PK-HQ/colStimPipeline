function runChipRMetatableOnce()
% Regenerate Chip R metatable (statisticsR-metatable16).
%
% Behavioral curves:
%   C1 cluster (15 blocks): RweibullfreeAllC1 in mdlStruct.
%   Block 16:               behavioralData raw fallback.
%
% Validates C1 overlap: reconstructs x/y arrays from raw behavioral for all
% 15 C1 blocks and verifies they match stored values (tol=1e-9).
%
% Usage:
%   cd('Y:\users\PK\colStimPipeline')
%   runChipRMetatableOnce

%% 1. Path setup
restoredefaultpath;
addpath(genpath('Y:/users/PK/colStimPipeline'));

%% 2. Source and output paths
SRC_MAT  = 'Y:/Chip/Meta/summary/statisticsR-final16.mat';
OUT_MAT  = 'Y:/Chip/Meta/summary/statisticsR-metatable16.mat';
OUT_XLSX = 'Y:/Chip/Meta/summary/statisticsR-metatable16.xlsx';

fprintf('Source : %s\n', SRC_MAT);
fprintf('MAT out: %s\n', OUT_MAT);
fprintf('XLSX out: %s\n', OUT_XLSX);

%% 3. Load source variables
fprintf('\nLoading source file...\n');
load(SRC_MAT, 'blockData', 'bitmapData', 'behavioralData', ...
     'analysisBlockID', 'datastruct', 'dataTag', 'mdlStruct');

%% 4. Generate metatable
fprintf('Building metatable...\n');
saveOpts = struct('matPath', OUT_MAT, 'xlsxPath', OUT_XLSX);
MetaTable = buildBlockMetadata(behavioralData, bitmapData, [], [], ...
    blockData, datastruct, analysisBlockID, mdlStruct, saveOpts);

%% 5. Basic row counts
fprintf('\n--- Row counts ---\n');
fprintf('height(MetaTable)        = %d\n', height(MetaTable));
fprintf('unique blockIDs          = %d\n', numel(unique(MetaTable.blockID)));

%% 6. Nonmissing psychometric y fields
fieldPairs = { ...
    'psy_yBaselineMerged',            'yBaselineMerged'; ...
    'psy_yConOptoMerged',             'yConOptoMerged'; ...
    'psy_yInconOptoMerged',           'yInconOptoMerged'; ...
    'psy_yBaselinePreMerge',          'yBaselinePreMerge'; ...
    'psy_yHorizontalOptoPreMerge',    'yHorizontalOptoPreMerge'; ...
    'psy_yVerticalOptoPreMerge',      'yVerticalOptoPreMerge'; ...
};

fprintf('\n--- Nonmissing psychometric y fields ---\n');
for fi = 1:size(fieldPairs, 1)
    col = fieldPairs{fi, 1};
    if ~any(strcmp(MetaTable.Properties.VariableNames, col))
        fprintf('  %-42s MISSING FROM TABLE\n', col);
        continue
    end
    vals = MetaTable.(col);
    n = sum(~cellfun(@isempty, vals));
    fprintf('  %-42s %d / %d\n', col, n, height(MetaTable));
end

%% 7. modelField distribution
fprintf('\n--- psy_modelField distribution ---\n');
mf = MetaTable.psy_modelField;
uMF = unique(mf);
for i = 1:numel(uMF)
    fprintf('  %-50s  count=%d\n', uMF{i}, sum(strcmp(mf, uMF{i})));
end

%% 8. Delta field coverage
fprintf('\n--- Delta field nonmissing counts ---\n');
deltaFields = {'psy_deltaBias','psy_deltaMask','psy_deltaBiasMerged','psy_deltaMaskMerged', ...
    'psy_deltaBiasHorizontal','psy_deltaMaskHorizontal', ...
    'psy_deltaBiasVertical','psy_deltaMaskVertical'};
for di = 1:numel(deltaFields)
    fn = deltaFields{di};
    if any(strcmp(MetaTable.Properties.VariableNames, fn))
        v = MetaTable.(fn);
        fprintf('  %-40s  %d / %d\n', fn, sum(~isnan(v)), height(MetaTable));
    end
end

%% 9. Validate raw-behavioral reconstruction vs C1 (15 blocks)
fprintf('\n--- Validating raw-behavioral reconstruction vs C1 (15 blocks) ---\n');
C1_FIELD = 'RweibullfreeAllC1';
if ~isstruct(mdlStruct) || ~isfield(mdlStruct, C1_FIELD)
    fprintf('  %s not found — skipping validation.\n', C1_FIELD);
else
    md = mdlStruct.(C1_FIELD);
    if ~isfield(md, 'clusterBlocksIdx')
        fprintf('  clusterBlocksIdx missing — skipping validation.\n');
    else
        c1BlockIDs = md.clusterBlocksIdx(:);
        xRaw = behavioralData.gaborContrasts;
        yRaw = behavioralData.percentageCorrect;

        cmpFields = { ...
            'psy_xBaselineMerged',             'xBaselineMerged'; ...
            'psy_yBaselineMerged',             'yBaselineMerged'; ...
            'psy_xConOptoMerged',              'xConOptoMerged'; ...
            'psy_yConOptoMerged',              'yConOptoMerged'; ...
            'psy_xInconOptoMerged',            'xInconOptoMerged'; ...
            'psy_yInconOptoMerged',            'yInconOptoMerged'; ...
            'psy_xBaselinePreMerge',           'xBaselinePreMerge'; ...
            'psy_yBaselinePreMerge',           'yBaselinePreMerge'; ...
            'psy_xBaselinePreMergeOrt',        'xBaselinePreMergeOrt'; ...
            'psy_xHorizontalOptoPreMerge',     'xHorizontalOptoPreMerge'; ...
            'psy_yHorizontalOptoPreMerge',     'yHorizontalOptoPreMerge'; ...
            'psy_xHorizontalOptoPreMergeOrt',  'xHorizontalOptoPreMergeOrt'; ...
            'psy_congruencyHorizontalOptoPreMerge', 'congruencyHorizontalOptoPreMerge'; ...
            'psy_xVerticalOptoPreMerge',       'xVerticalOptoPreMerge'; ...
            'psy_yVerticalOptoPreMerge',       'yVerticalOptoPreMerge'; ...
            'psy_xVerticalOptoPreMergeOrt',    'xVerticalOptoPreMergeOrt'; ...
            'psy_congruencyVerticalOptoPreMerge', 'congruencyVerticalOptoPreMerge'; ...
        };

        TOL = 1e-9;
        mismatchCounts = zeros(size(cmpFields, 1), 1);
        nTested = 0;

        for bi = 1:numel(c1BlockIDs)
            bID = c1BlockIDs(bi);
            rowIdx = find(MetaTable.blockID == bID, 1);
            if isempty(rowIdx)
                fprintf('  WARNING: C1 block %d not found in MetaTable\n', bID);
                continue
            end
            eRaw = compute_raw_psy_for_block(xRaw, yRaw, bID);
            nTested = nTested + 1;
            for fi = 1:size(cmpFields, 1)
                tblCol = cmpFields{fi, 1};
                rawFld = cmpFields{fi, 2};
                if ~any(strcmp(MetaTable.Properties.VariableNames, tblCol)); continue; end
                vTbl = MetaTable.(tblCol){rowIdx};
                vRaw = eRaw.(rawFld);
                vTbl = vTbl(~isnan(vTbl));
                vRaw = vRaw(~isnan(vRaw));
                if numel(vTbl) ~= numel(vRaw) || any(abs(vTbl(:) - vRaw(:)) > TOL)
                    mismatchCounts(fi) = mismatchCounts(fi) + 1;
                end
            end
        end

        fprintf('  Blocks tested: %d / %d\n', nTested, numel(c1BlockIDs));
        anyMismatch = false;
        for fi = 1:size(cmpFields, 1)
            if mismatchCounts(fi) > 0
                fprintf('  MISMATCH  %-50s  %d blocks\n', cmpFields{fi,1}, mismatchCounts(fi));
                anyMismatch = true;
            end
        end
        if ~anyMismatch
            fprintf('  All %d fields match for all %d C1 blocks. Validation PASSED.\n', ...
                size(cmpFields,1), nTested);
        else
            fprintf('  WARNING: mismatches found — inspect before using.\n');
        end
    end
end

%% 10. Missing expected columns
fprintf('\n--- Missing expected columns ---\n');
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
    'bmp_horizontalCamSpace', 'bmp_verticalCamSpace', ...
    'session_hasVisFPS', 'session_baselineSource', ...
};
tblCols = MetaTable.Properties.VariableNames;
nMissing = 0;
for ci = 1:numel(expectedCols)
    if ~any(strcmp(tblCols, expectedCols{ci}))
        fprintf('  MISSING: %s\n', expectedCols{ci});
        nMissing = nMissing + 1;
    end
end
if nMissing == 0
    fprintf('  All expected columns present.\n');
end

fprintf('\nDone.\n');
end


% =========================================================================
%  Raw-behavioral processing helpers
% =========================================================================

function e = compute_raw_psy_for_block(xBlocks, yBlocks, blockID)
e = struct();
e.xBaselineMerged = []; e.yBaselineMerged = [];
e.xConOptoMerged = []; e.yConOptoMerged = [];
e.xInconOptoMerged = []; e.yInconOptoMerged = [];
e.xBaselinePreMerge = []; e.yBaselinePreMerge = []; e.xBaselinePreMergeOrt = [];
e.xHorizontalOptoPreMerge = []; e.yHorizontalOptoPreMerge = [];
e.xHorizontalOptoPreMergeOrt = []; e.congruencyHorizontalOptoPreMerge = [];
e.xVerticalOptoPreMerge = []; e.yVerticalOptoPreMerge = [];
e.xVerticalOptoPreMergeOrt = []; e.congruencyVerticalOptoPreMerge = [];

xBaselineRaw  = local_rmnan(squeeze(xBlocks(1,:,blockID)));
yBaselineRaw  = local_rmnan(squeeze(yBlocks(1,:,blockID)));
xHorizontalRaw = local_rmnan(squeeze(xBlocks(2,:,blockID)));
yHorizontalRaw = local_rmnan(squeeze(yBlocks(2,:,blockID)));
xVerticalRaw  = local_rmnan(squeeze(xBlocks(3,:,blockID)));
yVerticalRaw  = local_rmnan(squeeze(yBlocks(3,:,blockID)));

nBase = numel(xBaselineRaw);
if nBase > 0 && mod(nBase, 2) == 0
    tagBRaw = local_make_visual_tag(nBase);
    [xBSort, si] = sort(xBaselineRaw);
    yBSort = yBaselineRaw(si);
    tagBSort = tagBRaw(si);
    [xBPre, yBPre, tagBPre] = local_make_pre_merge_pct_correct(xBSort, yBSort, tagBSort, true);
    e.xBaselinePreMerge    = xBPre(:)';
    e.yBaselinePreMerge    = yBPre(:)';
    e.xBaselinePreMergeOrt = tagBPre(:)';
    numVal = numel(xBSort);
    e.xBaselineMerged = mean([fliplr(-xBSort(1:numVal/2)); xBSort(numVal/2+1:end)])';
    e.yBaselineMerged = mean([fliplr(100-yBSort(1:numVal/2)); yBSort(numVal/2+1:end)])';
end

nH = numel(xHorizontalRaw);
if nH > 0 && mod(nH, 2) == 0
    tagHRaw = local_make_visual_tag(nH);
    [xHSort, si] = sort(xHorizontalRaw);
    yHSort = yHorizontalRaw(si); tagHSort = tagHRaw(si);
    [xHPre, yHPre, tagHPre] = local_make_pre_merge_pct_correct(xHSort, yHSort, tagHSort, false);
    congrH = NaN(size(tagHPre));
    congrH(tagHPre == 0) = 1; congrH(tagHPre == 90) = -1;
    e.xHorizontalOptoPreMerge       = xHPre(:)';
    e.yHorizontalOptoPreMerge       = yHPre(:)';
    e.xHorizontalOptoPreMergeOrt    = tagHPre(:)';
    e.congruencyHorizontalOptoPreMerge = congrH(:)';
end

nV = numel(xVerticalRaw);
if nV > 0 && mod(nV, 2) == 0
    tagVRaw = local_make_visual_tag(nV);
    [xVSort, si] = sort(xVerticalRaw);
    yVSort = yVerticalRaw(si); tagVSort = tagVRaw(si);
    [xVPre, yVPre, tagVPre] = local_make_pre_merge_pct_correct(xVSort, yVSort, tagVSort, false);
    congrV = NaN(size(tagVPre));
    congrV(tagVPre == 0) = -1; congrV(tagVPre == 90) = 1;
    e.xVerticalOptoPreMerge         = xVPre(:)';
    e.yVerticalOptoPreMerge         = yVPre(:)';
    e.xVerticalOptoPreMergeOrt      = tagVPre(:)';
    e.congruencyVerticalOptoPreMerge = congrV(:)';
end

if nH > 0 && mod(nH,2)==0 && nV > 0 && mod(nV,2)==0
    cNeg = 1:nH/2; cPos = nH/2+1:nH;
    xHS = sort(xHorizontalRaw); yHS = yHorizontalRaw(local_argsort(xHorizontalRaw));
    xVS = sort(xVerticalRaw);   yVS = yVerticalRaw(local_argsort(xVerticalRaw));
    yHC = yHS; yHC(cNeg) = 100-yHC(cNeg);
    yVC = yVS; yVC(cNeg) = 100-yVC(cNeg);
    e.xConOptoMerged  = mean([-fliplr(xHS(cNeg)); xVS(cPos)])';
    e.yConOptoMerged  = mean([fliplr(yHC(cNeg));  yVC(cPos)])';
    e.xInconOptoMerged = mean([xHS(cPos); -fliplr(xVS(cNeg))])';
    e.yInconOptoMerged = mean([yHC(cPos); fliplr(yVC(cNeg))])';
end
end


function v = local_rmnan(v)
v = v(~isnan(v(:))');
end

function tag = local_make_visual_tag(nVal)
tag = NaN(1, nVal);
tag(1:nVal/2) = 0;
tag(nVal/2+1:end) = 90;
end

function [xOut, yOut, tagOut] = local_make_pre_merge_pct_correct(xIn, yIn, tagIn, mergeDupZeros)
xIn = xIn(:)'; yIn = yIn(:)'; tagIn = tagIn(:)';
valid = ~isnan(xIn) & ~isnan(yIn) & ~isnan(tagIn);
xIn = xIn(valid); yIn = yIn(valid); tagIn = tagIn(valid);
yCorrect = yIn;
yCorrect(tagIn == 0) = 100 - yCorrect(tagIn == 0);
if mergeDupZeros
    [xOut, yOut, tagOut] = local_merge_dup_x_for_baseline(xIn, yCorrect, tagIn);
else
    xOut = xIn; yOut = yCorrect; tagOut = tagIn;
end
end

function [xOut, yOut, tagOut] = local_merge_dup_x_for_baseline(xIn, yIn, tagIn)
xIn = xIn(:)'; yIn = yIn(:)'; tagIn = tagIn(:)';
valid = ~isnan(xIn) & ~isnan(yIn) & ~isnan(tagIn);
xIn = xIn(valid); yIn = yIn(valid); tagIn = tagIn(valid);
[xOut, ~, grp] = unique(xIn, 'stable');
yOut = nan(size(xOut)); tagOut = nan(size(xOut));
for ii = 1:numel(xOut)
    sel = grp == ii;
    yOut(ii) = mean(yIn(sel), 'omitnan');
    utags = unique(tagIn(sel)); utags = utags(~isnan(utags));
    if numel(utags) == 1; tagOut(ii) = utags; else; tagOut(ii) = 45; end
end
end

function idx = local_argsort(v)
[~, idx] = sort(v);
end
