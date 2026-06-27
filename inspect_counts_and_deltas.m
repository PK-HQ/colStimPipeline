% inspect_counts_and_deltas.m
% One-shot inspection: trial count sources, delta definitions, Chip R block 16.
restoredefaultpath;
addpath(genpath('Y:/users/PK/colStimPipeline'));

%% ---- Chip L ----
fprintf('=== CHIP L: all MAT variables ===\n');
wL = whos('-file','Y:/Chip/Meta/summary/statisticsL-final43.mat');
for i=1:numel(wL)
    fprintf('  %-45s %s %s\n', wL(i).name, mat2str(wL(i).size), wL(i).class);
end

fprintf('\n=== CHIP L: full behavioralData fields ===\n');
bL = load('Y:/Chip/Meta/summary/statisticsL-final43.mat','behavioralData');
bd = bL.behavioralData;
fn = fieldnames(bd);
for i=1:numel(fn)
    v = bd.(fn{i});
    if isstruct(v)
        fprintf('  %-30s struct  fields=%d  numel=%d\n', fn{i}, numel(fieldnames(v)), numel(v));
    elseif iscell(v)
        fprintf('  %-30s cell    %s\n', fn{i}, mat2str(size(v)));
    else
        fprintf('  %-30s %-8s %s\n', fn{i}, class(v), mat2str(size(v)));
    end
end

fprintf('\n=== CHIP L: LweibullSignedBX0C1 sub-struct fields ===\n');
mL = load('Y:/Chip/Meta/summary/statisticsL-final43.mat','mdlStruct');
md = mL.mdlStruct.LweibullSignedBX0C1;
fn2 = fieldnames(md);
for i=1:numel(fn2)
    v = md.(fn2{i});
    if isstruct(v)
        fn3 = fieldnames(v);
        fprintf('  %-40s struct  fields=%d\n', fn2{i}, numel(fn3));
        for j=1:numel(fn3)
            vv = v.(fn3{j});
            if isnumeric(vv) || islogical(vv)
                fprintf('    %-38s %-8s %s\n', fn3{j}, class(vv), mat2str(size(vv)));
            end
        end
    elseif isnumeric(v) || islogical(v)
        fprintf('  %-40s %-8s %s\n', fn2{i}, class(v), mat2str(size(v)));
    end
end

fprintf('\n=== CHIP L: sample delta verification (block 4, kIdx=1) ===\n');
md = mL.mdlStruct.LweibullSignedBX0C1;
kIdx = 1; bID = md.clusterBlocksIdx(kIdx);
fprintf('  block=%d (kIdx=%d)\n', bID, kIdx);
fprintf('  stored deltaBiasHorizontal   = %.6f\n', md.deltaBiasHorizontal(kIdx));
fprintf('  stored deltaMaskHorizontal   = %.6f\n', md.deltaMaskHorizontal(kIdx));
fprintf('  stored deltaBiasVertical     = %.6f\n', md.deltaBiasVertical(kIdx));
fprintf('  stored deltaMaskVertical     = %.6f\n', md.deltaMaskVertical(kIdx));
fprintf('  stored deltaBiasMerged       = %.6f\n', md.deltaBiasMerged(kIdx));
fprintf('  stored deltaMaskMerged       = %.6f\n', md.deltaMaskMerged(kIdx));
fprintf('  stored meanBaselineHorizontal= %.6f\n', md.meanBaselineHorizontal(kIdx));
fprintf('  stored meanConOptoHorizontal = %.6f\n', md.meanConOptoHorizontal(kIdx));
fprintf('  stored meanInconOptoHorizontal=%.6f\n', md.meanInconOptoHorizontal(kIdx));
fprintf('  stored meanBaselineVertical  = %.6f\n', md.meanBaselineVertical(kIdx));
fprintf('  stored meanConOptoVertical   = %.6f\n', md.meanConOptoVertical(kIdx));
fprintf('  stored meanInconOptoVertical = %.6f\n', md.meanInconOptoVertical(kIdx));
fprintf('  stored meanBaselineMerged    = %.6f\n', md.meanBaselineMerged(kIdx));
fprintf('  stored meanConOptoMerged     = %.6f\n', md.meanConOptoMerged(kIdx));
fprintf('  stored meanInconOptoMerged   = %.6f\n', md.meanInconOptoMerged(kIdx));

% Try to reproduce from pre-merge y/x/n fields using signed choice data
if isfield(md,'signedBX0') && isfield(md.signedBX0,'nBaselineChoice')
    nb  = md.signedBX0.nBaselineChoice(kIdx,:);
    yb  = md.signedBX0.yBaselineChoice(kIdx,:);
    nh  = md.signedBX0.nHorizontalOptoChoice(kIdx,:);
    yh  = md.signedBX0.yHorizontalOptoChoice(kIdx,:);
    nv  = md.signedBX0.nVerticalOptoChoice(kIdx,:);
    yv  = md.signedBX0.yVerticalOptoChoice(kIdx,:);
    nb = nb(~isnan(nb)); yb = yb(~isnan(yb));
    nh = nh(~isnan(nh)); yh = yh(~isnan(yh));
    nv = nv(~isnan(nv)); yv = yv(~isnan(yv));
    fprintf('  nBaseline  (nonNaN): %s\n', mat2str(nb));
    fprintf('  yBaseline  (nonNaN): %s\n', mat2str(yb,4));
    fprintf('  nHorizontal(nonNaN): %s\n', mat2str(nh));
    fprintf('  nVertical  (nonNaN): %s\n', mat2str(nv));

    % Weighted mean
    meanBH_w = sum(nb(yb>=0) .* yb(yb>=0)) / sum(nb(yb>=0));  % rough approx
    fprintf('  [debug] raw weighted mean baseline (all): %.4f\n', ...
        sum(nb .* yb)/sum(nb));

    % Check if meanPanel1Horizontal contains [H-baseline, H-con, H-incon]
    if isfield(md,'meanPanel1Horizontal')
        fprintf('  meanPanel1Horizontal(kIdx,:) = %s\n', mat2str(md.meanPanel1Horizontal(kIdx,:),6));
        fprintf('  meanPanel1Headers = ');
        disp(md.meanPanel1Headers);
    end
end

%% ---- Chip R ----
fprintf('\n=== CHIP R: all MAT variables ===\n');
wR = whos('-file','Y:/Chip/Meta/summary/statisticsR-final16.mat');
for i=1:numel(wR)
    fprintf('  %-45s %s %s\n', wR(i).name, mat2str(wR(i).size), wR(i).class);
end

fprintf('\n=== CHIP R: behavioralData fields ===\n');
bR = load('Y:/Chip/Meta/summary/statisticsR-final16.mat','behavioralData');
bd2 = bR.behavioralData;
fn = fieldnames(bd2);
for i=1:numel(fn)
    v = bd2.(fn{i});
    if isstruct(v)
        fprintf('  %-30s struct  fields=%d  numel=%d\n', fn{i}, numel(fieldnames(v)), numel(v));
    elseif iscell(v)
        fprintf('  %-30s cell    %s\n', fn{i}, mat2str(size(v)));
    else
        fprintf('  %-30s %-8s %s\n', fn{i}, class(v), mat2str(size(v)));
    end
end

fprintf('\n=== CHIP R: mdlStruct fields ===\n');
mR = load('Y:/Chip/Meta/summary/statisticsR-final16.mat','mdlStruct');
fn = fieldnames(mR.mdlStruct);
for i=1:numel(fn)
    v = mR.mdlStruct.(fn{i});
    if isstruct(v) && isfield(v,'clusterBlocksIdx')
        fprintf('  %-45s has clusterBlocksIdx, nK=%d, blockIDs=%s\n', fn{i}, numel(v.clusterBlocksIdx), mat2str(v.clusterBlocksIdx(:)'));
    elseif isstruct(v)
        fprintf('  %-45s struct, no clusterBlocksIdx\n', fn{i});
    else
        fprintf('  %-45s %s\n', fn{i}, class(v));
    end
end

fprintf('\n=== CHIP R: gaborContrasts block 16 (raw) ===\n');
xR16 = squeeze(bd2.gaborContrasts(:,:,16));
yR16 = squeeze(bd2.percentageCorrect(:,:,16));
for row=1:3
    x = xR16(row,:); x = x(~isnan(x));
    y = yR16(row,:); y = y(~isnan(y));
    fprintf('  row%d: x=%s\n', row, mat2str(x,4));
    fprintf('       y=%s\n', mat2str(y,4));
end

fprintf('\n=== CHIP R: check any C1 cluster for block 16 in R ===\n');
fn = fieldnames(mR.mdlStruct);
for i=1:numel(fn)
    v = mR.mdlStruct.(fn{i});
    if isstruct(v) && isfield(v,'clusterBlocksIdx')
        if any(v.clusterBlocksIdx == 16)
            fprintf('  %s contains block 16\n', fn{i});
        else
            fprintf('  %s does NOT contain block 16 (blocks: %s)\n', fn{i}, mat2str(v.clusterBlocksIdx(:)'));
        end
    end
end

fprintf('\nDone.\n');
