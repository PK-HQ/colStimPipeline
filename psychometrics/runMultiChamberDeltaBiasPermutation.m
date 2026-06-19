function results = runMultiChamberDeltaBiasPermutation(chamberResults, binEdges, opts)
% Run saved-data per-experiment con-vs-incon permutation tests.
%
% Empirical sources are the row-1-column-3 merged psychometric markers saved
% by plotNakaRushtonFit5:
%   con contrasts:       mdl.xBlock(2,:,sessionRowIndex)
%   con percent correct: mdl.yBlock(2,:,sessionRowIndex)
%   incon contrasts:     mdl.xBlock(3,:,sessionRowIndex)
%   incon percent correct: mdl.yBlock(3,:,sessionRowIndex)
% Baseline values, fitted curves, fitted delta curves, and tolerance matching
% are intentionally not used.

    if nargin < 3 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'outputDir') || isempty(opts.outputDir)
        opts.outputDir = fullfile('Y:\users\PK\colStimPipeline', ...
            'outputs', 'multiChamber');
    end
    if ~isfield(opts, 'modelName') || isempty(opts.modelName)
        opts.modelName = 'weibullfreeAll';
    end
    if ~isfield(opts, 'nPermutations') || isempty(opts.nPermutations)
        opts.nPermutations = 500;
    end
    if ~isfield(opts, 'barColor') || isempty(opts.barColor)
        opts.barColor = [0.70 0.45 0.95];
    end
    if ~isfield(opts, 'significantBarColor') || isempty(opts.significantBarColor)
        opts.significantBarColor = [0.42 0.10 0.65];
    end
    if ~isfield(opts, 'meanLineColor') || isempty(opts.meanLineColor)
        opts.meanLineColor = [0.45 0.10 0.70];
    end
    if ~exist(opts.outputDir, 'dir')
        mkdir(opts.outputDir);
    end

    rng(1, 'twister');

    experimentTable = table();
    contrastTable = table();
    nullExperimentMeans = struct([]);
    reconstructionWarnings = table();

    fprintf('\nStage 2 permutation source fields:\n');
    fprintf('  empirical con contrasts: mdl.xBlock(2,:,sessionRowIndex)\n');
    fprintf('  empirical con percent correct: mdl.yBlock(2,:,sessionRowIndex)\n');
    fprintf('  empirical incon contrasts: mdl.xBlock(3,:,sessionRowIndex)\n');
    fprintf('  empirical incon percent correct: mdl.yBlock(3,:,sessionRowIndex)\n');
    fprintf('  experiment identity: distribution source audit experimentID/sessionRowIndex/blockIndex\n');
    fprintf('  chamber identity: configured animal/chamber plus saved psychfit file\n');

    nullIdx = 0;
    for chamberIdx = 1:numel(chamberResults)
        result = chamberResults(chamberIdx);
        fprintf('\nLoading Stage 2 model source for %s-%s: %s\n', ...
            result.monkeyID, result.chamber, result.modelPath);
        mdl = loadChamberModel(result, opts.modelName);

        retainedRows = result.experimentTable(result.experimentTable.includedInHistogram, :);
        for rowIdx = 1:height(retainedRows)
            row = retainedRows(rowIdx, :);
            sessionRow = row.sessionRowIndex;
            if ~isfinite(sessionRow) || sessionRow ~= round(sessionRow) || ...
                    sessionRow < 1 || sessionRow > size(mdl.xBlock, 3)
                error('runMultiChamberDeltaBiasPermutation:InvalidSessionRow', ...
                    ['Invalid sessionRowIndex for %s %s experiment %s: %g. ' ...
                    'mdl.xBlock third dimension has %d rows.'], ...
                    result.animal, result.chamber, string(row.experimentID), ...
                    sessionRow, size(mdl.xBlock, 3));
            end

            conX = squeeze(mdl.xBlock(2, :, sessionRow));
            conY = squeeze(mdl.yBlock(2, :, sessionRow));
            inconX = squeeze(mdl.xBlock(3, :, sessionRow));
            inconY = squeeze(mdl.yBlock(3, :, sessionRow));

            permResult = computeExperimentDeltaBiasPermutation( ...
                conX, conY, inconX, inconY, opts);

            nullIdx = nullIdx + 1;
            nullExperimentMeans(nullIdx).monkeyID = char(row.monkeyID);
            nullExperimentMeans(nullIdx).animal = char(row.animal);
            nullExperimentMeans(nullIdx).chamber = char(row.chamber);
            nullExperimentMeans(nullIdx).experimentID = char(string(row.experimentID));
            nullExperimentMeans(nullIdx).sessionRowIndex = sessionRow;
            nullExperimentMeans(nullIdx).blockIndex = row.blockIndex;
            nullExperimentMeans(nullIdx).powerClusterID = row.clusterID;
            nullExperimentMeans(nullIdx).nullMeanDelta = permResult.nullMeanDelta;

            experimentTable = [experimentTable; buildExperimentRow(row, permResult)]; %#ok<AGROW>
            contrastTable = [contrastTable; buildContrastRows(row, permResult)]; %#ok<AGROW>

            warnRows = reconstructionWarningRows(row, permResult);
            reconstructionWarnings = [reconstructionWarnings; warnRows]; %#ok<AGROW>
        end
    end

    validatePermutationResults(chamberResults, experimentTable, contrastTable, binEdges);
    stage2Figure = plotStage2Histogram(chamberResults, experimentTable, binEdges, opts);
    savePermutationOutputs(stage2Figure, experimentTable, contrastTable, ...
        nullExperimentMeans, reconstructionWarnings, opts);
    printPermutationSummary(chamberResults, experimentTable, reconstructionWarnings);

    results = struct();
    results.experimentTable = experimentTable;
    results.contrastTable = contrastTable;
    results.nullExperimentMeans = nullExperimentMeans;
    results.reconstructionWarnings = reconstructionWarnings;
end

function mdl = loadChamberModel(result, modelName)
    if ~isfile(result.modelPath)
        error('runMultiChamberDeltaBiasPermutation:MissingModelFile', ...
            'Saved psychometric model file is missing: %s', result.modelPath);
    end
    loaded = load(result.modelPath, 'mdlStruct', 'analysisBlockID', 'datastruct', 'dataTag');
    if ~isfield(loaded, 'mdlStruct')
        error('runMultiChamberDeltaBiasPermutation:MissingMdlStruct', ...
            'File does not contain mdlStruct: %s', result.modelPath);
    end

    modelField = sprintf('%s%sC1', result.chamber, modelName);
    if ~isfield(loaded.mdlStruct, modelField)
        error('runMultiChamberDeltaBiasPermutation:MissingModelField', ...
            'mdlStruct is missing %s. Available fields: %s', ...
            modelField, strjoin(fieldnames(loaded.mdlStruct), ', '));
    end
    mdl = loaded.mdlStruct.(modelField);
    if ~isfield(mdl, 'xBlock') || ~isfield(mdl, 'yBlock')
        error('runMultiChamberDeltaBiasPermutation:MissingMarkerFields', ...
            'mdlStruct.%s must contain xBlock and yBlock empirical marker arrays.', ...
            modelField);
    end
    if size(mdl.xBlock, 1) < 3 || size(mdl.yBlock, 1) < 3
        error('runMultiChamberDeltaBiasPermutation:InvalidMarkerShape', ...
            'mdlStruct.%s xBlock/yBlock must have at least 3 condition rows.', ...
            modelField);
    end
    fprintf('  using mdlStruct.%s.xBlock/yBlock, size xBlock=%s, yBlock=%s\n', ...
        modelField, mat2str(size(mdl.xBlock)), mat2str(size(mdl.yBlock)));
end

function rowOut = buildExperimentRow(row, permResult)
    rowOut = table( ...
        string(row.monkeyID), string(row.animal), string(row.chamber), ...
        string(row.experimentID), row.sessionRowIndex, row.blockIndex, ...
        row.clusterID, row.clusterIncluded, ...
        numel(permResult.matchedContrasts), ...
        permResult.observedMeanDelta, ...
        permResult.nullMeanDeltaMean, permResult.nullMeanDeltaMedian, ...
        permResult.nullMeanDeltaLower95, permResult.nullMeanDeltaUpper95, ...
        permResult.oneSidedExperimentP, permResult.experimentSignificant, ...
        permResult.nPermutations, ...
        'VariableNames', {'monkeyID', 'animal', 'chamber', 'experimentID', ...
        'sessionRowIndex', 'blockIndex', 'powerClusterID', 'clusterIncluded', ...
        'nMatchedContrasts', 'observedMeanDeltaBias', ...
        'nullMeanDeltaBiasMean', 'nullMeanDeltaBiasMedian', ...
        'nullMeanDeltaBiasLower95', 'nullMeanDeltaBiasUpper95', ...
        'oneSidedExperimentP', 'experimentSignificant', 'nPermutations'});
end

function contrastRows = buildContrastRows(row, permResult)
    nRows = numel(permResult.matchedContrasts);
    contrastRows = table( ...
        repmat(string(row.monkeyID), nRows, 1), ...
        repmat(string(row.animal), nRows, 1), ...
        repmat(string(row.chamber), nRows, 1), ...
        repmat(string(row.experimentID), nRows, 1), ...
        repmat(row.clusterID, nRows, 1), ...
        permResult.matchedContrasts(:), ...
        permResult.conPct(:), permResult.inconPct(:), ...
        permResult.observedDelta(:), ...
        permResult.nCon(:), permResult.nIncon(:), ...
        permResult.conCorrect(:), permResult.inconCorrect(:), ...
        permResult.conReconstructionErrorPct(:), ...
        permResult.inconReconstructionErrorPct(:), ...
        permResult.nullDeltaMean(:), permResult.nullDeltaMedian(:), ...
        permResult.nullDeltaLower95(:), permResult.nullDeltaUpper95(:), ...
        permResult.oneSidedContrastP(:), ...
        repmat(permResult.nPermutations, nRows, 1), ...
        'VariableNames', {'monkeyID', 'animal', 'chamber', 'experimentID', ...
        'powerClusterID', 'contrast', 'conPct', 'inconPct', ...
        'observedDeltaBias', 'nCon', 'nIncon', 'conCorrect', ...
        'inconCorrect', 'conReconstructionErrorPct', ...
        'inconReconstructionErrorPct', 'nullDeltaMean', ...
        'nullDeltaMedian', 'nullDeltaLower95', 'nullDeltaUpper95', ...
        'oneSidedContrastP', 'nPermutations'});
end

function warnRows = reconstructionWarningRows(row, permResult)
    warnMask = abs(permResult.conReconstructionErrorPct) > 0.25 | ...
        abs(permResult.inconReconstructionErrorPct) > 0.25;
    nRows = sum(warnMask);
    if nRows == 0
        warnRows = table();
        return;
    end
    warnRows = table( ...
        repmat(string(row.monkeyID), nRows, 1), ...
        repmat(string(row.animal), nRows, 1), ...
        repmat(string(row.chamber), nRows, 1), ...
        repmat(string(row.experimentID), nRows, 1), ...
        permResult.matchedContrasts(warnMask), ...
        permResult.conReconstructionErrorPct(warnMask), ...
        permResult.inconReconstructionErrorPct(warnMask), ...
        'VariableNames', {'monkeyID', 'animal', 'chamber', ...
        'experimentID', 'contrast', 'conReconstructionErrorPct', ...
        'inconReconstructionErrorPct'});
    for warnIdx = 1:nRows
        warning('runMultiChamberDeltaBiasPermutation:ReconstructionError', ...
            ['%s-%s experiment %s contrast %g reconstruction error exceeds ' ...
            '0.25 pp: con=%0.3f, incon=%0.3f'], ...
            string(row.monkeyID), string(row.chamber), string(row.experimentID), ...
            warnRows.contrast(warnIdx), ...
            warnRows.conReconstructionErrorPct(warnIdx), ...
            warnRows.inconReconstructionErrorPct(warnIdx));
    end
end

function validatePermutationResults(chamberResults, experimentTable, contrastTable, binEdges)
    for chamberIdx = 1:numel(chamberResults)
        result = chamberResults(chamberIdx);
        retainedRows = result.experimentTable(result.experimentTable.includedInHistogram, :);
        expRows = experimentTable(experimentTable.animal == string(result.animal) & ...
            experimentTable.chamber == string(result.chamber), :);

        assert(height(expRows) == height(retainedRows), ...
            'Every retained experiment must have exactly one permutation result.');
        assert(numel(unique(string(expRows.experimentID))) == height(expRows), ...
            'Duplicate experiment-level permutation rows found.');
        assert(all(expRows.clusterIncluded), ...
            'Permutation result includes an experiment from an excluded cluster.');

        for rowIdx = 1:height(expRows)
            contrastRows = contrastTable(contrastTable.animal == expRows.animal(rowIdx) & ...
                contrastTable.chamber == expRows.chamber(rowIdx) & ...
                contrastTable.experimentID == expRows.experimentID(rowIdx), :);
            assert(height(contrastRows) == expRows.nMatchedContrasts(rowIdx), ...
                'Contrast-level rows must match nMatchedContrasts.');
            assert(max(abs(contrastRows.observedDeltaBias - ...
                (contrastRows.conPct - contrastRows.inconPct))) < 1e-10, ...
                'Observed contrast delta does not equal empirical con minus incon.');
        end

        sigMask = logical(expRows.experimentSignificant);
        sigCounts = histcounts(retainedRows.deltaBiasMerged(sigMask), binEdges);
        nonsigCounts = histcounts(retainedRows.deltaBiasMerged(~sigMask), binEdges);
        stage1Counts = histcounts(retainedRows.deltaBiasMerged, binEdges);
        assert(isequal(sigCounts + nonsigCounts, stage1Counts), ...
            'Stacked Stage 2 bins do not match Stage 1 bins.');
        assert(sum(sigCounts + nonsigCounts) == height(retainedRows), ...
            'Stacked Stage 2 proportions do not sum to retained n.');
    end
end

function figHandle = plotStage2Histogram(chamberResults, experimentTable, binEdges, opts)
    maxProp = 0;
    for chamberIdx = 1:numel(chamberResults)
        result = chamberResults(chamberIdx);
        retainedRows = result.experimentTable(result.experimentTable.includedInHistogram, :);
        expRows = experimentTable(experimentTable.animal == string(result.animal) & ...
            experimentTable.chamber == string(result.chamber), :);
        sigMask = logical(expRows.experimentSignificant);
        nonsigCounts = histcounts(retainedRows.deltaBiasMerged(~sigMask), binEdges);
        sigCounts = histcounts(retainedRows.deltaBiasMerged(sigMask), binEdges);
        if height(retainedRows) > 0
            maxProp = max(maxProp, max((nonsigCounts + sigCounts) ./ height(retainedRows)));
        end
    end
    yMax = min(1, max(0.1, maxProp + 0.08));

    figHandle = figure('Color', 'w', 'Position', [100 100 1450 430]);
    tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
    binCenters = binEdges(1:end-1) + diff(binEdges) ./ 2;

    legendHandles = gobjects(2, 1);
    for chamberIdx = 1:numel(chamberResults)
        result = chamberResults(chamberIdx);
        retainedRows = result.experimentTable(result.experimentTable.includedInHistogram, :);
        expRows = experimentTable(experimentTable.animal == string(result.animal) & ...
            experimentTable.chamber == string(result.chamber), :);
        sigMask = logical(expRows.experimentSignificant);
        nonsigCounts = histcounts(retainedRows.deltaBiasMerged(~sigMask), binEdges);
        sigCounts = histcounts(retainedRows.deltaBiasMerged(sigMask), binEdges);
        if height(retainedRows) > 0
            stackedProps = [nonsigCounts(:), sigCounts(:)] ./ height(retainedRows);
        else
            stackedProps = zeros(numel(binCenters), 2);
        end

        ax = nexttile;
        hold(ax, 'on');
        barHandles = bar(ax, binCenters, stackedProps, 1.0, 'stacked');
        barHandles(1).FaceColor = opts.barColor;
        barHandles(1).FaceAlpha = 0.35;
        barHandles(1).EdgeColor = opts.meanLineColor;
        barHandles(2).FaceColor = opts.significantBarColor;
        barHandles(2).FaceAlpha = 0.65;
        barHandles(2).EdgeColor = opts.meanLineColor;
        if chamberIdx == 1
            legendHandles = barHandles;
        end
        xline(ax, 0, '--', 'Color', [0.45 0.45 0.45], ...
            'LineWidth', 1.2, 'HandleVisibility', 'off');
        if ~isempty(retainedRows.deltaBiasMerged)
            xline(ax, mean(retainedRows.deltaBiasMerged), '-', ...
                'Color', opts.meanLineColor, 'LineWidth', 1.8, ...
                'HandleVisibility', 'off');
            rugY = 0.03 .* yMax;
            plot(ax, retainedRows.deltaBiasMerged, ...
                repmat(rugY, height(retainedRows), 1), '|', ...
                'Color', opts.meanLineColor, 'MarkerSize', 9, ...
                'LineWidth', 1.1, 'HandleVisibility', 'off');
        end
        title(ax, sprintf('%s %s %s\nn = %d | sig = %d | clusters %s / %d', ...
            result.monkeyID, char(8212), result.chamber, height(retainedRows), ...
            sum(sigMask), clusterListText(result.includedClusterIDs), ...
            result.totalClusters), 'Interpreter', 'none');
        xlabel(ax, 'Merged \DeltaBias (% correct)');
        ylabel(ax, 'Proportion of included experiments');
        xlim(ax, [binEdges(1) binEdges(end)]);
        ylim(ax, [0 yMax]);
        box(ax, 'off');
        set(ax, 'TickDir', 'out', 'LineWidth', 1);
    end

    legend(legendHandles, {'experiment p >= 0.05', 'experiment p < 0.05'}, ...
        'Location', 'northeastoutside');
    sgtitle('Merged experiment-wise \DeltaBias by per-experiment permutation significance');
end

function savePermutationOutputs(figHandle, experimentTable, contrastTable, ...
        nullExperimentMeans, reconstructionWarnings, opts)
    experimentPath = fullfile(opts.outputDir, ...
        'multiChamberDeltaBiasPermutation_Experiments.csv');
    contrastPath = fullfile(opts.outputDir, ...
        'multiChamberDeltaBiasPermutation_Contrasts.csv');
    matPath = fullfile(opts.outputDir, ...
        'multiChamberDeltaBiasPermutation.mat');
    pdfPath = fullfile(opts.outputDir, ...
        'multiChamberDeltaBiasHistogram_byExperimentSignificance.pdf');
    pngPath = fullfile(opts.outputDir, ...
        'multiChamberDeltaBiasHistogram_byExperimentSignificance.png');
    figPath = fullfile(opts.outputDir, ...
        'multiChamberDeltaBiasHistogram_byExperimentSignificance.fig');

    writetable(experimentTable, experimentPath);
    writetable(contrastTable, contrastPath);
    save(matPath, 'experimentTable', 'contrastTable', ...
        'nullExperimentMeans', 'reconstructionWarnings');
    exportgraphics(figHandle, pdfPath, 'ContentType', 'vector', ...
        'BackgroundColor', 'white');
    exportgraphics(figHandle, pngPath, 'Resolution', 300, ...
        'BackgroundColor', 'white');
    savefig(figHandle, figPath);

    fprintf('\nSaved Stage 2 outputs:\n');
    fprintf('  %s\n', experimentPath);
    fprintf('  %s\n', contrastPath);
    fprintf('  %s\n', matPath);
    fprintf('  %s\n', pdfPath);
    fprintf('  %s\n', pngPath);
    fprintf('  %s\n', figPath);
end

function printPermutationSummary(chamberResults, experimentTable, reconstructionWarnings)
    fprintf('\nStage 2 concise summary:\n');
    for chamberIdx = 1:numel(chamberResults)
        result = chamberResults(chamberIdx);
        expRows = experimentTable(experimentTable.animal == string(result.animal) & ...
            experimentTable.chamber == string(result.chamber), :);
        fprintf('%s-%s: included clusters %s / %d, n=%d, significant experiments=%d\n', ...
            result.monkeyID, result.chamber, clusterListText(result.includedClusterIDs), ...
            result.totalClusters, height(expRows), sum(expRows.experimentSignificant));
    end
    if isempty(reconstructionWarnings)
        fprintf('  reconstruction-error warnings: none > 0.25 percentage points\n');
    else
        fprintf('  reconstruction-error warnings: %d contrast rows > 0.25 percentage points\n', ...
            height(reconstructionWarnings));
    end
end

function textValue = clusterListText(clusterIDs)
    clusterIDs = clusterIDs(:)';
    if isempty(clusterIDs)
        textValue = 'none';
        return;
    end
    labels = arrayfun(@(x) sprintf('C%d', x), clusterIDs, ...
        'UniformOutput', false);
    textValue = strjoin(labels, ', ');
end
