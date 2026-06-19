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
    clusterBlocks = clusterMdl.clusterBlocksIdx(:)';
    if isempty(clusterBlocks)
        clusterBlocks = 1:size(clusterMdl.fittedParams, 1);
    end

    aggregateField = sprintf('%s%sC1Aggregate', chamberWanted, modelType);
    if ~isfield(mdlStruct, aggregateField)
        error('generateDeltaBiasPermutationInspectionReport:MissingAggregateField', ...
            'Missing mdlStruct.%s in saved model structures.', aggregateField);
    end
    aggregatePsychometrics = mdlStruct.(aggregateField);

    plotOpts = struct('showDeltaPermutationStats', true, ...
        'nDeltaPermutations', 500, ...
        'deltaPermutationBaseSeed', 99173);

    reportFilename = sprintf('psychometrics/%s-chamber/%s/psychfit-%s-%s-C1_permStatsInspection', ...
        chamberWanted, modelType, chamberWanted, modelType);
    reportState = initializeReportPDFAssembly(reportFilename, monkeyName);
    xFit = sort(nlinspace(0, 100, 100, 'linear'));
    clusterLabels = makeInspectionClusterLabels(mdlStruct, aggregateField, clusterBlocks);
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

    inspection = struct();
    inspection.summaryFile = summaryFile;
    inspection.psychFile = psychFile;
    inspection.pdfPath = finalPdf;
    inspection.auditPaths = auditPaths;
    inspection.clusterMdl = clusterMdl;
    inspection.aggregateFits = aggregateFits;
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

function labels = makeInspectionClusterLabels(mdlStruct, aggregateField, clusterBlocks)
    labels = strings(numel(clusterBlocks), 1);
    labels(:) = "";
    labelsByBlockField = [aggregateField 'LabelsByBlock'];
    if ~isfield(mdlStruct, labelsByBlockField)
        return;
    end
    labelsByBlock = mdlStruct.(labelsByBlockField);
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
