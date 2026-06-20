function inspection = generateDeltaBiasPermutationInspectionReport(monkeyName, chamberWanted, modelType)
% Generate a saved-data Pepper-R deltaBias permutation inspection report.
%
% This uses existing saved psychometric model/results structures only. It does
% not refit models, rerun clustering, or alter source MAT files.

    if nargin < 1 || isempty(monkeyName)
        monkeyName = 'Pepper';
    end
    if nargin < 2 || isempty(chamberWanted)
        chamberWanted = 'R';
    end
    if nargin < 3 || isempty(modelType)
        modelType = 'weibullfreeAll';
    end

    mainPath = 'Y:';
    summaryDir = fullfile(mainPath, monkeyName, 'Meta', 'summary');
    psychDir = fullfile(mainPath, monkeyName, 'Meta', 'psychometrics');
    summaryFile = findLatestFinalSummary(summaryDir, chamberWanted);
    psychFile = fullfile(psychDir, sprintf('psychfit%s-%s.mat', chamberWanted, modelType));
    if ~isfile(psychFile)
        error('generateDeltaBiasPermutationInspectionReport:MissingPsychFile', ...
            'Missing saved psychfit file: %s', psychFile);
    end

    summaryData = load(summaryFile, 'bitmapData', 'behavioralData', ...
        'analysisBlockID', 'datastruct', 'mdlStruct');
    psychData = load(psychFile, 'mdlStruct');
    if isfield(psychData, 'mdlStruct')
        mdlStruct = psychData.mdlStruct;
    else
        mdlStruct = summaryData.mdlStruct;
    end

    modelField = sprintf('%s%sC1', chamberWanted, modelType);
    if ~isfield(mdlStruct, modelField)
        error('generateDeltaBiasPermutationInspectionReport:MissingModelField', ...
            'Missing mdlStruct.%s in %s', modelField, psychFile);
    end
    clusterMdl = mdlStruct.(modelField);
    nModelRows = size(clusterMdl.fittedParams, 1);
    mappingMethod = '';
    if isfield(clusterMdl, 'clusterBlocksIdx') && ...
            ~isempty(clusterMdl.clusterBlocksIdx)
        clusterBlocks = clusterMdl.clusterBlocksIdx(:)';
        mappingMethod = 'clusterMdl.clusterBlocksIdx';
    elseif isfield(clusterMdl, 'blockIndices') && numel(clusterMdl.blockIndices) == nModelRows
        clusterBlocks = clusterMdl.blockIndices(:)';
        mappingMethod = 'clusterMdl.blockIndices';
    elseif isfield(clusterMdl, 'sourceBlockIndices') && numel(clusterMdl.sourceBlockIndices) == nModelRows
        clusterBlocks = clusterMdl.sourceBlockIndices(:)';
        mappingMethod = 'clusterMdl.sourceBlockIndices';
    else
        clusterBlocks = 1:nModelRows;
        mappingMethod = 'validated fitted-row indices';
        warning('generateDeltaBiasPermutationInspectionReport:MissingClusterBlocksIdx', ...
            ['Saved clusterMdl lacks explicit block indices; using fitted-row ' ...
            'indices after validating row/block dimensions.']);
    end
    validateInspectionMapping(clusterBlocks, nModelRows, summaryData);
    printInspectionMapping(mappingMethod, clusterBlocks, nModelRows, summaryData);

    aggregateFieldCandidates = { ...
        sprintf('%s%sPowerClusterAggregate', chamberWanted, modelType), ...
        sprintf('%s%sC1Aggregate', chamberWanted, modelType)};
    aggregateField = '';
    for candidateIdx = 1:numel(aggregateFieldCandidates)
        if isfield(mdlStruct, aggregateFieldCandidates{candidateIdx})
            aggregateField = aggregateFieldCandidates{candidateIdx};
            break;
        end
    end
    if isempty(aggregateField)
        fprintf('Pepper-R aggregate field: none saved; rebuilding from individual model rows.\n');
        [aggregatePsychometrics, powerClusterLabels] = ...
            rebuildInspectionAggregate(clusterMdl, clusterBlocks, summaryData);
    else
        fprintf('Pepper-R aggregate field: mdlStruct.%s\n', aggregateField);
        aggregatePsychometrics = mdlStruct.(aggregateField);
        powerClusterLabels = resolveSavedClusterLabels(mdlStruct, aggregateField, clusterBlocks);
    end

    outputDir = fullfile('Y:\users\PK\colStimPipeline', 'outputs', ...
        'deltaBiasPermutation');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    plotOpts = struct('showDeltaPermutationStats', true, ...
        'nDeltaPermutations', 500, ...
        'deltaPermutationBaseSeed', 99173, ...
        'saveDeltaPermutationExamples', true, ...
        'deltaPermutationExampleDir', outputDir);
    fprintf('Focused generator plot options: showDeltaPermutationStats=%d | nDeltaPermutations=%d\n', ...
        plotOpts.showDeltaPermutationStats, plotOpts.nDeltaPermutations);

    reportFilename = fullfile(outputDir, 'psychfit-Pepper-R_permStatsInspection.pdf');
    reportState = initializeReportPDFAssembly(reportFilename, monkeyName);
    clusterMdl = ensureInspectionPreMergeFields(clusterMdl);
    xFit = sort(nlinspace(0, 100, 100, 'linear'));
    clusterLabels = makeInspectionClusterLabels(mdlStruct, aggregateField, clusterBlocks, powerClusterLabels);
    fitParams = clusterMdl.fittedParams;
    if ndims(fitParams) == 3
        fitParams = fitParams(:, :, 1);
    end

    [clusterMdl, reportState] = plotNakaRushtonFit5( ...
        summaryData.behavioralData, summaryData.bitmapData, ...
        summaryData.datastruct, summaryData.analysisBlockID, ...
        clusterMdl, fitParams, xFit, monkeyName, clusterBlocks, ...
        0, 1, true, 1, modelType, reportFilename, ...
        clusterLabels, reportState, plotOpts);

    [aggregateFits, aggregateFigures] = ...
        plotAggregatedPowerClusterPsychometrics(aggregatePsychometrics, plotOpts);
    reportState = stageAndCloseInspectionFigures(aggregateFigures, reportState);
    finalPdf = finalizeReportPDFAssembly(reportState);

    auditPaths = saveDeltaBiasPermutationAuditOutputs( ...
        monkeyName, chamberWanted, modelType, clusterMdl, aggregateFits);

    individualPng = fullfile(outputDir, 'PepperR_individualDeltaPermutationExample.png');
    aggregatePng = fullfile(outputDir, 'PepperR_aggregateDeltaPermutationExample.png');
    assertOutputFile(finalPdf, 'inspection PDF');
    assertOutputFile(individualPng, 'individual example PNG');
    assertOutputFile(aggregatePng, 'aggregate example PNG');

    inspection = struct();
    inspection.summaryFile = summaryFile;
    inspection.psychFile = psychFile;
    inspection.pdfPath = finalPdf;
    inspection.auditPaths = auditPaths;
    inspection.individualPng = individualPng;
    inspection.aggregatePng = aggregatePng;
    inspection.mappingMethod = mappingMethod;
    inspection.clusterBlocks = clusterBlocks;
    inspection.nIndividualAnnotated = countIndividualPermutationPanels(clusterMdl);
    inspection.nAggregateAnnotated = countAggregatePermutationPanels(aggregateFits);
    inspection.clusterMdl = clusterMdl;
    inspection.aggregateFits = aggregateFits;
    fprintf('Permutation inspection complete: individual panels %d | aggregate panels %d\n', ...
        inspection.nIndividualAnnotated, inspection.nAggregateAnnotated);
end

function labels = resolveSavedClusterLabels(mdlStruct, aggregateField, clusterBlocks)
    labels = [];
    labelsField = [aggregateField 'Labels'];
    labelsByBlockField = [aggregateField 'LabelsByBlock'];
    if isfield(mdlStruct, labelsField) && numel(mdlStruct.(labelsField)) == numel(clusterBlocks)
        savedLabels = mdlStruct.(labelsField);
        labels = savedLabels(:);
    elseif isfield(mdlStruct, labelsByBlockField)
        labelsByBlock = mdlStruct.(labelsByBlockField);
        if max(clusterBlocks) <= numel(labelsByBlock)
            labels = labelsByBlock(clusterBlocks);
        end
    end
end

function [aggregatePsychometrics, powerClusterLabels] = rebuildInspectionAggregate(clusterMdl, clusterBlocks, summaryData)
    clusterOpts = struct( ...
        'minSessionsPerCluster', 3, ...
        'effectWeight', 0.75, ...
        'monotonicPenalty', 2, ...
        'kCandidates', [2 3], ...
        'minRelativeImprovementFor3Clusters', 0.10);
    aggregateOpts = struct('binWidth', 5, ...
        'useMedianForPlot', false, 'minSessionsPerBin', 1);
    powerAudit = auditPowerMetrics(summaryData.bitmapData, clusterBlocks, ...
        summaryData.analysisBlockID(clusterBlocks), [], []);
    sessionPower = powerAudit.Ptotal_recomputed(:);
    [sessionBias, sessionMask, deltaSource] = resolveInspectionSessionDeltas(clusterMdl);
    fprintf('Pepper-R clustering delta source: %s\n', deltaSource);
    try
        [experimentalByBlock, controlByBlock] = classifyOptostimColumnTargets( ...
            summaryData.bitmapData.orts, max(clusterBlocks));
        experimentalSessionMask = experimentalByBlock(clusterBlocks);
        controlSessionMask = controlByBlock(clusterBlocks);
        unknownSessionMask = ~(experimentalSessionMask | controlSessionMask);
    catch ortError
        warning('generateDeltaBiasPermutationInspectionReport:OrtMaskUnavailable', ...
            ['Could not align bitmapData.orts for saved-data inspection (%s). ' ...
            'Using all sessions for rebuilt aggregate clusters.'], ortError.message);
        experimentalSessionMask = true(size(clusterBlocks(:)));
        controlSessionMask = false(size(clusterBlocks(:)));
        unknownSessionMask = false(size(clusterBlocks(:)));
    end
    if ~any(experimentalSessionMask)
        experimentalSessionMask = true(size(clusterBlocks(:)));
        controlSessionMask = false(size(clusterBlocks(:)));
        unknownSessionMask = false(size(clusterBlocks(:)));
        warning('generateDeltaBiasPermutationInspectionReport:NoOrtExperimentalMask', ...
            'No 0/90 mask found; rebuilding clusters using all sessions.');
    end
    [experimentalCluster, ~, diagnostics] = clusterOrderedPowerEffect( ...
        sessionPower(experimentalSessionMask), ...
        sessionBias(experimentalSessionMask), ...
        sessionMask(experimentalSessionMask), clusterOpts);
    powerClusterLabels = assignOrderedPowerBands(sessionPower, diagnostics.boundaries);
    if any(powerClusterLabels(experimentalSessionMask) ~= experimentalCluster)
        error('generateDeltaBiasPermutationInspectionReport:ClusterRebuildMismatch', ...
            'Rebuilt fixed-boundary cluster labels do not match orderedPowerEffect output.');
    end
    powerClusterLabels(controlSessionMask) = diagnostics.chosenK;
    powerClusterLabels(unknownSessionMask) = NaN;
    fprintf('Rebuilt aggregate power-cluster labels: %s\n', mat2str(powerClusterLabels(:)'));
    clusterMdl = ensureInspectionSidePanelFields(clusterMdl);
    aggregatePsychometrics = aggregatePsychometricCountsByPowerCluster( ...
        clusterMdl, powerClusterLabels, aggregateOpts);
end
function [sessionBias, sessionMask, sourceName] = resolveInspectionSessionDeltas(clusterMdl)
    if isfield(clusterMdl, 'deltaBias') && isfield(clusterMdl, 'deltaMask')
        sessionBias = clusterMdl.deltaBias(:);
        sessionMask = clusterMdl.deltaMask(:);
        sourceName = 'clusterMdl.deltaBias/deltaMask';
        return;
    end
    if isfield(clusterMdl, 'meanDeltaMerged') && size(clusterMdl.meanDeltaMerged, 2) >= 2
        sessionBias = clusterMdl.meanDeltaMerged(:, 1);
        sessionMask = clusterMdl.meanDeltaMerged(:, 2);
        sourceName = 'clusterMdl.meanDeltaMerged(:,1:2)';
        return;
    end
    if isfield(clusterMdl, 'meanPsychometricMerged') && size(clusterMdl.meanPsychometricMerged, 2) >= 3
        baseline = clusterMdl.meanPsychometricMerged(:, 1);
        con = clusterMdl.meanPsychometricMerged(:, 2);
        incon = clusterMdl.meanPsychometricMerged(:, 3);
        sessionBias = con - incon;
        sessionMask = baseline - mean([con incon], 2, 'omitnan');
        sourceName = 'clusterMdl.meanPsychometricMerged condition means';
        return;
    end
    if isfield(clusterMdl, 'xBlock') && isfield(clusterMdl, 'yBlock') && ...
            size(clusterMdl.yBlock, 1) >= 3
        nRows = size(clusterMdl.yBlock, 3);
        baseline = nan(nRows, 1);
        con = nan(nRows, 1);
        incon = nan(nRows, 1);
        for rowIdx = 1:nRows
            baseline(rowIdx) = mean(squeeze(clusterMdl.yBlock(1, :, rowIdx)), 'omitnan');
            con(rowIdx) = mean(squeeze(clusterMdl.yBlock(2, :, rowIdx)), 'omitnan');
            incon(rowIdx) = mean(squeeze(clusterMdl.yBlock(3, :, rowIdx)), 'omitnan');
        end
        sessionBias = con - incon;
        sessionMask = baseline - mean([con incon], 2, 'omitnan');
        sourceName = 'clusterMdl.yBlock condition means';
        return;
    end
    error('generateDeltaBiasPermutationInspectionReport:MissingDeltaSource', ...
        'Could not find individual-session deltaBias/deltaMask source in clusterMdl.');
end
function clusterMdl = ensureInspectionSidePanelFields(clusterMdl)
    required = {'xPanel1Horizontal', 'yPanel1Horizontal', ...
        'xPanel1Vertical', 'yPanel1Vertical'};
    if all(isfield(clusterMdl, required))
        return;
    end
    if ~isfield(clusterMdl, 'xBlock') || ~isfield(clusterMdl, 'yBlock')
        error('generateDeltaBiasPermutationInspectionReport:MissingAggregateSideSource', ...
            'Cannot rebuild aggregate side panels because xBlock/yBlock are missing.');
    end
    clusterMdl.xPanel1Horizontal = clusterMdl.xBlock(1:3, :, :);
    clusterMdl.yPanel1Horizontal = clusterMdl.yBlock(1:3, :, :);
    clusterMdl.xPanel1Vertical = clusterMdl.xBlock(1:3, :, :);
    clusterMdl.yPanel1Vertical = clusterMdl.yBlock(1:3, :, :);
    warning('generateDeltaBiasPermutationInspectionReport:SynthesizedSidePanels', ...
        ['Saved model lacks horizontal/vertical side-panel fields; using ' ...
        'merged xBlock/yBlock as side-panel placeholders for focused ' ...
        'aggregate inspection plotting.']);
end
function clusterMdl = ensureInspectionPreMergeFields(clusterMdl)
    required = {'xBaselinePreMerge', 'yBaselinePreMerge', ...
        'xHorizontalOptoPreMerge', 'yHorizontalOptoPreMerge', ...
        'visualTagHorizontalOptoPreMerge', 'congruencyHorizontalOptoPreMerge', ...
        'xVerticalOptoPreMerge', 'yVerticalOptoPreMerge', ...
        'visualTagVerticalOptoPreMerge', 'congruencyVerticalOptoPreMerge'};
    if all(isfield(clusterMdl, required))
        return;
    end
    if ~isfield(clusterMdl, 'xBlock') || ~isfield(clusterMdl, 'yBlock')
        error('generateDeltaBiasPermutationInspectionReport:MissingPreMergeSource', ...
            'Cannot synthesize pre-merged fields because xBlock/yBlock are missing.');
    end
    xBase = squeezeToRows(clusterMdl.xBlock(1, :, :));
    yBase = squeezeToRows(clusterMdl.yBlock(1, :, :));
    xCon = squeezeToRows(clusterMdl.xBlock(2, :, :));
    yCon = squeezeToRows(clusterMdl.yBlock(2, :, :));
    xIncon = squeezeToRows(clusterMdl.xBlock(3, :, :));
    yIncon = squeezeToRows(clusterMdl.yBlock(3, :, :));
    clusterMdl.xBaselinePreMerge = xBase;
    clusterMdl.yBaselinePreMerge = yBase;
    xH = [-abs(xCon), -abs(xIncon)];
    yH = [yCon, yIncon];
    xV = [abs(xCon), abs(xIncon)];
    yV = [yCon, yIncon];
    tagH = [zeros(size(xCon)), zeros(size(xIncon))];
    congrH = [ones(size(xCon)), zeros(size(xIncon))];
    tagV = [90 .* ones(size(xCon)), 90 .* ones(size(xIncon))];
    congrV = [ones(size(xCon)), zeros(size(xIncon))];
    invalidH = ~isfinite(xH) | ~isfinite(yH);
    invalidV = ~isfinite(xV) | ~isfinite(yV);
    tagH(invalidH) = NaN;
    congrH(invalidH) = NaN;
    tagV(invalidV) = NaN;
    congrV(invalidV) = NaN;
    clusterMdl.xHorizontalOptoPreMerge = xH;
    clusterMdl.yHorizontalOptoPreMerge = yH;
    clusterMdl.visualTagHorizontalOptoPreMerge = tagH;
    clusterMdl.congruencyHorizontalOptoPreMerge = congrH;
    clusterMdl.xVerticalOptoPreMerge = xV;
    clusterMdl.yVerticalOptoPreMerge = yV;
    clusterMdl.visualTagVerticalOptoPreMerge = tagV;
    clusterMdl.congruencyVerticalOptoPreMerge = congrV;
    warning('generateDeltaBiasPermutationInspectionReport:SynthesizedPreMergeFields', ...
        ['Saved model lacks pre-merged panel fields; synthesized focused ' ...
        'inspection placeholders from merged xBlock/yBlock.']);
end

function rows = squeezeToRows(values)
    rows = squeeze(values)';
    if isvector(rows)
        rows = rows(:)';
    end
end
function validateInspectionMapping(clusterBlocks, nModelRows, summaryData)
    assert(numel(clusterBlocks) == nModelRows, ...
        'Cluster block mapping has %d entries for %d model rows.', ...
        numel(clusterBlocks), nModelRows);
    assert(all(isfinite(clusterBlocks)) && all(clusterBlocks >= 1) && ...
        all(abs(clusterBlocks - round(clusterBlocks)) < eps), ...
        'Cluster block mapping contains invalid indices.');
    assert(isfield(summaryData, 'analysisBlockID'), ...
        'Loaded summary data lacks analysisBlockID.');
    assert(max(clusterBlocks) <= numel(summaryData.analysisBlockID), ...
        'Cluster block index exceeds analysisBlockID length.');
end

function printInspectionMapping(mappingMethod, clusterBlocks, nModelRows, summaryData)
    expIDs = summaryData.analysisBlockID(clusterBlocks);
    fprintf('Pepper-R inspection mapping: %s\n', mappingMethod);
    fprintf('model rows: %d\n', nModelRows);
    fprintf('source blocks: %s\n', mat2str(clusterBlocks(:)'));
    fprintf('experiment IDs: %s\n', mat2str(expIDs(:)'));
end

function assertOutputFile(outputPath, label)
    assert(isfile(outputPath), 'Missing %s: %s', label, outputPath);
    fileInfo = dir(outputPath);
    assert(fileInfo.bytes > 0, 'Empty %s: %s', label, outputPath);
    fprintf('%s: %s (%d bytes)\n', label, outputPath, fileInfo.bytes);
end

function n = countIndividualPermutationPanels(clusterMdl)
    n = 0;
    if isfield(clusterMdl, 'deltaBiasPermutationExperimentSummary') && ...
            istable(clusterMdl.deltaBiasPermutationExperimentSummary)
        n = height(clusterMdl.deltaBiasPermutationExperimentSummary);
    end
end

function n = countAggregatePermutationPanels(aggregateFits)
    n = 0;
    for idx = 1:numel(aggregateFits)
        if isfield(aggregateFits(idx), 'merged') && ...
                isfield(aggregateFits(idx).merged, 'deltaBiasPermutationClusterSummary') && ...
                istable(aggregateFits(idx).merged.deltaBiasPermutationClusterSummary) && ...
                height(aggregateFits(idx).merged.deltaBiasPermutationClusterSummary) > 0
            n = n + 1;
        end
    end
end
function summaryFile = findLatestFinalSummary(summaryDir, chamberWanted)
    files = dir(fullfile(summaryDir, sprintf('statistics%s-final*.mat', chamberWanted)));
    if isempty(files)
        error('generateDeltaBiasPermutationInspectionReport:MissingSummary', ...
            'No statistics%s-final*.mat file found in %s', chamberWanted, summaryDir);
    end
    [~, idx] = max([files.datenum]);
    summaryFile = fullfile(files(idx).folder, files(idx).name);
end

function labels = makeInspectionClusterLabels(mdlStruct, aggregateField, clusterBlocks, powerClusterLabels)
    labels = strings(numel(clusterBlocks), 1);
    labels(:) = "";
    if nargin >= 4 && ~isempty(powerClusterLabels)
        labelsByBlock = nan(max(clusterBlocks), 1);
        labelsByBlock(clusterBlocks) = powerClusterLabels(:);
    else
        labelsByBlockField = [aggregateField 'LabelsByBlock'];
        if isempty(aggregateField) || ~isfield(mdlStruct, labelsByBlockField)
            return;
        end
        labelsByBlock = mdlStruct.(labelsByBlockField);
    end
    valid = labelsByBlock(isfinite(labelsByBlock) & labelsByBlock > 0);
    if isempty(valid)
        return;
    end
    nClusters = max(valid);
    for idx = 1:numel(clusterBlocks)
        blockIdx = clusterBlocks(idx);
        if blockIdx <= numel(labelsByBlock) && isfinite(labelsByBlock(blockIdx)) && ...
                labelsByBlock(blockIdx) > 0
            labels(idx) = sprintf('Cluster: %d/%d', labelsByBlock(blockIdx), nClusters);
        end
    end
end

function reportState = stageAndCloseInspectionFigures(figureHandles, reportState)
    for idx = 1:numel(figureHandles)
        if isempty(figureHandles(idx)) || ~isgraphics(figureHandles(idx))
            continue;
        end
        reportState = stageReportPDFPage(reportState, figureHandles(idx));
        close(figureHandles(idx));
    end
end
