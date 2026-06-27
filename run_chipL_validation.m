restoredefaultpath;
addpath(genpath('Y:/users/PK/colStimPipeline'));
load('Y:/Chip/Meta/summary/statisticsL-final43.mat','blockData','bitmapData','behavioralData','analysisBlockID','datastruct','dataTag','mdlStruct');
saveOpts = struct('matPath','Y:/Chip/Meta/summary/statisticsL-metatable43.mat','xlsxPath','Y:/Chip/Meta/summary/statisticsL-metatable43.xlsx');
MetaTable = buildBlockMetadata(behavioralData, bitmapData, [], [], blockData, datastruct, analysisBlockID, mdlStruct, saveOpts);
fprintf('height=%d\n', height(MetaTable));
fprintf('unique_blockIDs=%d\n', numel(unique(MetaTable.blockID)));
pmy = MetaTable.psy_modelField;
fprintf('nonmissing_modelField=%d\n', sum(~cellfun(@isempty, pmy)));
ybl = MetaTable.psy_yBaseline;
fprintf('nonmissing_yBaseline=%d\n', sum(~cellfun(@isempty, ybl)));
yco = MetaTable.psy_yConOpto;
fprintf('nonmissing_yConOpto=%d\n', sum(~cellfun(@isempty, yco)));
yic = MetaTable.psy_yInconOpto;
fprintf('nonmissing_yInconOpto=%d\n', sum(~cellfun(@isempty, yic)));
yBPre = MetaTable.psy_yBaselinePreMerge;
fprintf('nonmissing_yBaselinePreMerge=%d\n', sum(~cellfun(@isempty, yBPre)));
yHPre = MetaTable.psy_yHorizontalOptoPreMerge;
fprintf('nonmissing_yHorizontalOptoPreMerge=%d\n', sum(~cellfun(@isempty, yHPre)));
yVPre = MetaTable.psy_yVerticalOptoPreMerge;
fprintf('nonmissing_yVerticalOptoPreMerge=%d\n', sum(~cellfun(@isempty, yVPre)));
dbm = MetaTable.psy_deltaBiasMerged;
fprintf('nonmissing_deltaBiasMerged=%d\n', sum(~isnan(dbm)));
mf = MetaTable.psy_modelField;
uniqueMF = unique(mf);
for i=1:numel(uniqueMF)
    n = sum(strcmp(mf, uniqueMF{i}));
    fprintf('modelField=%s count=%d\n', uniqueMF{i}, n);
end
