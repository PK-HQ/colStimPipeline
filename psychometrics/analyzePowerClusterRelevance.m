function relevanceStruct = analyzePowerClusterRelevance( ...
    clusterMdl, powerEffectCluster, bitmapData, opts)
%ANALYZEPOWERCLUSTERRELEVANCE Relate fit parameters and stimulation to behavior.

    if nargin < 4 || isempty(opts)
        opts = struct();
    end

    opts = fillDefaultOpts(opts);
    rngSeed = 97131;

    validateClusterMdl(clusterMdl);

    fittedParams = clusterMdl.fittedParams;
    paramBL = fittedParams(:, 1:4);
    dCon = fittedParams(:, 5:8);
    dIncon = fittedParams(:, 9:12);

    paramNames = {'A', 'B', 'alpha', 'beta'};
    clusterLabels = powerEffectCluster(:);
    deltaBias = clusterMdl.deltaBias(:);
    deltaMask = clusterMdl.deltaMask(:);
    blockIDs = clusterMdl.clusterBlocksIdx(:);

    if numel(clusterLabels) ~= numel(deltaBias)
        error('analyzePowerClusterRelevance:ClusterLabelMismatch', ...
            'powerEffectCluster must have one label per clustered block.');
    end

    parameter = struct();
    parameter.predictorNames = paramNames;
    parameter.baselineParams = paramBL;
    parameter.conParams = paramBL + dCon;
    parameter.inconParams = paramBL + dIncon;
    parameter.biasPredictors = dCon - dIncon;
    parameter.maskPredictors = -(dCon + dIncon) ./ 2;
    parameter.allBlocks = analyzeOutcomePair( ...
        parameter.biasPredictors, deltaBias, ...
        parameter.maskPredictors, deltaMask, paramNames, opts.nBootstrap, ...
        rngSeed);
    parameter.byCluster = analyzeByCluster( ...
        parameter.biasPredictors, deltaBias, ...
        parameter.maskPredictors, deltaMask, paramNames, clusterLabels, ...
        opts.nBootstrap, rngSeed);

    stim = buildStimulationPredictors(bitmapData, blockIDs);
    stimulation = struct();
    stimulation.predictorTable = stim.table;
    stimulation.predictorNames = stim.predictorNames;
    stimulation.familyNames = stim.familyNames;
    stimulation.familyIndex = stim.familyIndex;
    stimulation.allBlocks = analyzeStimulationOutcomes( ...
        stim.values, deltaBias, deltaMask, stim.predictorNames, ...
        stim.familyIndex, stim.familyNames, opts.nBootstrap, rngSeed + 101);
    stimulation.byCluster = analyzeStimulationByCluster( ...
        stim.values, deltaBias, deltaMask, stim.predictorNames, ...
        stim.familyIndex, stim.familyNames, clusterLabels, ...
        opts.nBootstrap, rngSeed + 202);

    summaryTables = struct();
    summaryTables.parameterShapley = buildParameterShapleySummaryTable( ...
        parameter.allBlocks, paramNames);
    summaryTables.stimulationPredictorDiagnostics = stim.diagnostics;
    summaryTables.stimulationFamilyScores = stim.familyScoreTable;
    summaryTables.stimulationFamilyShapley = ...
        buildStimulationFamilyShapleySummaryTable(stimulation.allBlocks);
    fprintf('Stimulation predictor diagnostics:\n');
    disp(summaryTables.stimulationPredictorDiagnostics);
    fprintf('Stimulation family-score diagnostics:\n');
    disp(stim.familyScoreDiagnostics);

    relevanceStruct = struct();
    relevanceStruct.sourceMdlField = opts.sourceMdlField;
    relevanceStruct.aggregateField = opts.aggregateField;
    relevanceStruct.chamberWanted = opts.chamberWanted;
    relevanceStruct.monkeyName = opts.monkeyName;
    relevanceStruct.blockIDs = blockIDs;
    relevanceStruct.powerEffectCluster = clusterLabels;
    relevanceStruct.outcomes = struct('deltaBias', deltaBias, ...
        'deltaMask', deltaMask);
    relevanceStruct.parameter = parameter;
    relevanceStruct.stimulation = stimulation;
    relevanceStruct.summaryTables = summaryTables;
    relevanceStruct.interpretation = { ...
        'Parameter Shapley values partition in-sample multivariable R2.'; ...
        ['Stimulation-family Shapley values use family-level z-score ' ...
         'summaries to reduce over-interpretation of collinear raw predictors.']; ...
        ['In-sample R2 can be optimistic; use leave-one-out predictive R2 ' ...
         'when it is finite.']};

    if opts.makeFigures
        relevanceStruct.figureHandles = makeRelevanceFigures( ...
            relevanceStruct, opts);
    else
        relevanceStruct.figureHandles = [];
    end

    if opts.saveFlag == 1
        relevanceStruct.outputPaths = saveRelevanceOutputs( ...
            relevanceStruct, opts);
    else
        relevanceStruct.outputPaths = struct();
    end
end

function opts = fillDefaultOpts(opts)
    defaults = struct( ...
        'sourceMdlField', '', ...
        'aggregateField', '', ...
        'chamberWanted', '', ...
        'monkeyName', '', ...
        'mainPath', '', ...
        'saveFlag', 0, ...
        'nBootstrap', 1000, ...
        'makeFigures', true, ...
        'includeAllBlocksSummary', true);

    names = fieldnames(defaults);
    for idx = 1:numel(names)
        if ~isfield(opts, names{idx}) || isempty(opts.(names{idx}))
            opts.(names{idx}) = defaults.(names{idx});
        end
    end
end

function validateClusterMdl(clusterMdl)
    requiredFields = {'fittedParams', 'deltaBias', 'deltaMask', ...
        'clusterBlocksIdx'};
    for idx = 1:numel(requiredFields)
        if ~isfield(clusterMdl, requiredFields{idx})
            error('analyzePowerClusterRelevance:MissingField', ...
                'clusterMdl is missing required field "%s".', ...
                requiredFields{idx});
        end
    end

    if size(clusterMdl.fittedParams, 2) < 12
        error('analyzePowerClusterRelevance:InvalidFittedParams', ...
            'clusterMdl.fittedParams must contain at least 12 columns.');
    end

    nRows = size(clusterMdl.fittedParams, 1);
    if numel(clusterMdl.deltaBias) ~= nRows || ...
            numel(clusterMdl.deltaMask) ~= nRows || ...
            numel(clusterMdl.clusterBlocksIdx) ~= nRows
        error('analyzePowerClusterRelevance:RowMismatch', ...
            ['fittedParams, deltaBias, deltaMask, and clusterBlocksIdx ' ...
             'must have matching row counts.']);
    end
end

function result = analyzeOutcomePair(XBias, yBias, XMask, yMask, ...
    predictorNames, nBootstrap, rngSeed)
    result = struct();
    result.bias = analyzePredictorSet( ...
        XBias, yBias, predictorNames, nBootstrap, rngSeed);
    result.mask = analyzePredictorSet( ...
        XMask, yMask, predictorNames, nBootstrap, rngSeed + 17);
end

function byCluster = analyzeByCluster(XBias, yBias, XMask, yMask, ...
    predictorNames, clusterLabels, nBootstrap, rngSeed)
    clusterIDs = unique(clusterLabels(:), 'stable');
    byCluster = repmat(struct('clusterID', [], 'n', [], 'bias', [], ...
        'mask', []), numel(clusterIDs), 1);

    for clusterIdx = 1:numel(clusterIDs)
        clusterID = clusterIDs(clusterIdx);
        keep = clusterLabels == clusterID;
        byCluster(clusterIdx).clusterID = clusterID;
        byCluster(clusterIdx).n = nnz(keep);
        byCluster(clusterIdx).bias = analyzePredictorSet( ...
            XBias(keep, :), yBias(keep), predictorNames, nBootstrap, ...
            rngSeed + 100 * clusterIdx);
        byCluster(clusterIdx).mask = analyzePredictorSet( ...
            XMask(keep, :), yMask(keep), predictorNames, nBootstrap, ...
            rngSeed + 100 * clusterIdx + 31);
    end
end

function result = analyzePredictorSet(X, y, predictorNames, ...
    nBootstrap, rngSeed)
    validRows = all(isfinite(X), 2) & isfinite(y);
    Xv = X(validRows, :);
    yv = y(validRows);
    p = size(Xv, 2);

    result = struct();
    result.n = numel(yv);
    result.validRows = validRows;
    result.predictorNames = predictorNames;
    result.pairwise = pairwiseCorrelations(X, y, predictorNames);

    fitStats = ordinaryR2FromOLS(Xv, yv);
    result.OLS = fitStats;
    result.status = modelStatus(Xv, yv);

    shapleyR2 = shapleyR2Individual(Xv, yv);
    result.shapleyR2 = shapleyR2;
    result.percentOfModelR2 = percentContribution(shapleyR2);
    result.bootstrapPercentCI = bootstrapShapleyPercentCI( ...
        Xv, yv, nBootstrap, rngSeed, []);

    if canComputeLOO(Xv, yv)
        result.leaveOneOutPredictiveR2 = leaveOneOutPredictiveR2(Xv, yv);
    else
        result.leaveOneOutPredictiveR2 = NaN;
    end

    if p == 0
        result.individualPredictorOrder = {};
    else
        [~, orderIdx] = sort(result.percentOfModelR2, 'descend');
        result.individualPredictorOrder = predictorNames(orderIdx);
    end
end

function pairwise = pairwiseCorrelations(X, y, predictorNames)
    p = size(X, 2);
    pairwise = repmat(struct('predictor', '', 'n', 0, ...
        'pearsonR', NaN, 'pearsonP', NaN, ...
        'spearmanRho', NaN, 'spearmanP', NaN), p, 1);

    for predictorIdx = 1:p
        keep = isfinite(X(:, predictorIdx)) & isfinite(y);
        pairwise(predictorIdx).predictor = predictorNames{predictorIdx};
        pairwise(predictorIdx).n = nnz(keep);
        if nnz(keep) >= 3
            [rPearson, pPearson] = corr( ...
                X(keep, predictorIdx), y(keep), 'Type', 'Pearson');
            [rSpearman, pSpearman] = corr( ...
                X(keep, predictorIdx), y(keep), 'Type', 'Spearman');
            pairwise(predictorIdx).pearsonR = rPearson;
            pairwise(predictorIdx).pearsonP = pPearson;
            pairwise(predictorIdx).spearmanRho = rSpearman;
            pairwise(predictorIdx).spearmanP = pSpearman;
        end
    end
end

function status = modelStatus(X, y)
    n = size(X, 1);
    p = size(X, 2);
    status = 'ok';
    if p == 0
        status = 'no valid predictors';
    elseif n < 3
        status = 'too few valid rows';
    elseif n <= p + 1
        status = 'underpowered';
    elseif rank(centerAndScale(X)) < p
        status = 'rank deficient';
    elseif sum((y - mean(y)).^2) <= eps
        status = 'zero outcome variance';
    end
end

function canCompute = canComputeLOO(X, y)
    n = size(X, 1);
    p = size(X, 2);
    canCompute = n > p + 1 && rank(centerAndScale(X)) == p && ...
        sum((y - mean(y)).^2) > eps;
end

function fitStats = ordinaryR2FromOLS(X, y)
    fitStats = struct('R2', NaN, 'coefficients', [], ...
        'intercept', NaN, 'rank', 0);

    if size(X, 1) < 2 || isempty(X) || ...
            sum((y - mean(y)).^2) <= eps
        return
    end

    Xz = centerAndScale(X);
    design = [ones(size(Xz, 1), 1), Xz];
    beta = pinv(design) * y;
    yHat = design * beta;
    totalSS = sum((y - mean(y)).^2);
    residualSS = sum((y - yHat).^2);

    fitStats.R2 = max(0, 1 - residualSS / totalSS);
    fitStats.coefficients = beta(2:end);
    fitStats.intercept = beta(1);
    fitStats.rank = rank(Xz);
end

function Xz = centerAndScale(X)
    if isempty(X)
        Xz = X;
        return
    end

    mu = mean(X, 1, 'omitnan');
    sigma = std(X, 0, 1, 'omitnan');
    sigma(~isfinite(sigma) | sigma == 0) = 1;
    Xz = (X - mu) ./ sigma;
end

function shapleyValues = shapleyR2Individual(X, y)
    p = size(X, 2);
    shapleyValues = nan(1, p);
    if p == 0 || size(X, 1) < 2
        return
    end

    shapleyValues(:) = 0;
    factorialP = factorial(p);
    r2ByMask = nan(2^p, 1);
    r2ByMask(1) = 0;

    for mask = 1:(2^p - 1)
        usePredictor = logical(bitget(mask, 1:p));
        fitStats = ordinaryR2FromOLS(X(:, usePredictor), y);
        r2ByMask(mask + 1) = fitStats.R2;
    end

    for predictorIdx = 1:p
        others = setdiff(1:p, predictorIdx);
        for mask = 0:(2^numel(others) - 1)
            subset = false(1, p);
            subset(others) = logical(bitget(mask, 1:numel(others)));
            subsetSize = nnz(subset);
            withPredictor = subset;
            withPredictor(predictorIdx) = true;
            weight = factorial(subsetSize) * ...
                factorial(p - subsetSize - 1) / factorialP;
            shapleyValues(predictorIdx) = shapleyValues(predictorIdx) + ...
                weight * (r2ForSubset(r2ByMask, withPredictor) - ...
                r2ForSubset(r2ByMask, subset));
        end
    end
end

function value = r2ForSubset(r2ByMask, subset)
    mask = 0;
    for idx = 1:numel(subset)
        if subset(idx)
            mask = bitset(mask, idx);
        end
    end
    value = r2ByMask(mask + 1);
end

function percent = percentContribution(values)
    total = sum(values, 'omitnan');
    if ~isfinite(total) || total <= eps
        percent = nan(size(values));
    else
        percent = 100 .* values ./ total;
    end
end

function rowSums = rowSumOmitNan(values)
    rowSums = sum(values, 2, 'omitnan');
    allMissing = all(~isfinite(values), 2);
    rowSums(allMissing) = NaN;
end

function ci = bootstrapShapleyPercentCI(X, y, nBootstrap, rngSeed, groups)
    p = size(X, 2);
    if isempty(groups)
        nTerms = p;
    else
        nTerms = numel(unique(groups, 'stable'));
    end
    ci = nan(nTerms, 2);
    if nBootstrap <= 0 || size(X, 1) < 3 || p == 0
        return
    end

    stream = RandStream('mt19937ar', 'Seed', rngSeed);
    n = size(X, 1);
    bootValues = nan(nBootstrap, nTerms);
    for bootIdx = 1:nBootstrap
        rowIdx = randi(stream, n, n, 1);
        if isempty(groups)
            values = shapleyR2Individual(X(rowIdx, :), y(rowIdx));
        else
            values = shapleyR2Grouped(X(rowIdx, :), y(rowIdx), groups);
        end
        bootValues(bootIdx, :) = percentContribution(values);
    end

    ci(:, 1) = prctile(bootValues, 2.5, 1)';
    ci(:, 2) = prctile(bootValues, 97.5, 1)';
end

function looR2 = leaveOneOutPredictiveR2(X, y)
    n = size(X, 1);
    yHat = nan(n, 1);

    for rowIdx = 1:n
        train = true(n, 1);
        train(rowIdx) = false;
        XTrain = X(train, :);
        yTrain = y(train);
        mu = mean(XTrain, 1, 'omitnan');
        sigma = std(XTrain, 0, 1, 'omitnan');
        sigma(~isfinite(sigma) | sigma == 0) = 1;
        XTrainZ = (XTrain - mu) ./ sigma;
        xTestZ = (X(rowIdx, :) - mu) ./ sigma;
        beta = pinv([ones(nnz(train), 1), XTrainZ]) * yTrain;
        yHat(rowIdx) = [1, xTestZ] * beta;
    end

    totalSS = sum((y - mean(y)).^2);
    if totalSS <= eps
        looR2 = NaN;
    else
        looR2 = 1 - sum((y - yHat).^2) / totalSS;
    end
end

function stim = buildStimulationPredictors(bitmapData, blockIDs)
    fields = { ...
        'meanPowerDensityWithinROI_mWmm2', 'Energy'; ...
        'totalPowerToOnPixelsWithinROI_mW', 'Energy'; ...
        'projectorPowerDensity_mWmm2', 'Energy'; ...
        'nColumns', 'Spatial'; ...
        'pixelsON', 'Spatial'; ...
        'areaFinalROI', 'Spatial'; ...
        'areaPixelsONWithinROI', 'Spatial'; ...
        'spatialDutyCycleWithinROI', 'Spatial'; ...
        'temporalDutyCycle', 'Temporal'; ...
        'sensitivity', 'QC'; ...
        'adaptthresh', 'QC'};

    nBlocks = numel(blockIDs);
    values = nan(nBlocks, size(fields, 1));
    for fieldIdx = 1:size(fields, 1)
        for rowIdx = 1:nBlocks
            values(rowIdx, fieldIdx) = extractBlockScalar( ...
                bitmapData, fields{fieldIdx, 1}, blockIDs(rowIdx));
        end
    end

    predictorNames = fields(:, 1)';
    predictorFamilies = fields(:, 2)';
    familyNames = unique(predictorFamilies, 'stable');
    familyIndex = nan(1, numel(predictorFamilies));
    for familyIdx = 1:numel(familyNames)
        familyIndex(strcmp(predictorFamilies, familyNames{familyIdx})) = ...
            familyIdx;
    end

    stim = struct();
    stim.values = values;
    stim.predictorNames = predictorNames;
    stim.predictorFamilies = predictorFamilies;
    stim.familyNames = familyNames;
    stim.familyIndex = familyIndex;
    stim.diagnostics = buildStimulationPredictorDiagnostics( ...
        values, predictorNames, predictorFamilies);
    validPredictor = stim.diagnostics.NFiniteRows >= 5 & ...
        stim.diagnostics.NUniqueFiniteValues >= 2 & ...
        isfinite(stim.diagnostics.Std) & stim.diagnostics.Std > 0;
    [familyScores, familyScoreNames, retainedPredictors] = ...
        buildFamilyScorePredictors(values, predictorNames, familyIndex, ...
        familyNames, validPredictor);
    stim.familyScores = familyScores;
    stim.familyScoreNames = familyScoreNames;
    stim.retainedFamilyPredictors = retainedPredictors;
    stim.familyScoreDiagnostics = buildStimulationPredictorDiagnostics( ...
        familyScores, familyScoreNames, familyScoreNames);
    if isempty(familyScoreNames)
        stim.familyScoreTable = table(blockIDs(:), ...
            'VariableNames', {'blockID'});
    else
        stim.familyScoreTable = array2table(familyScores, ...
            'VariableNames', familyScoreNames);
        stim.familyScoreTable.blockID = blockIDs(:);
    end
    stim.table = array2table(values, 'VariableNames', predictorNames);
    stim.table.blockID = blockIDs(:);
end

function value = extractBlockScalar(bitmapData, fieldName, blockID)
    value = NaN;
    if ~isfield(bitmapData, fieldName)
        return
    end

    data = bitmapData.(fieldName);
    if isempty(data) || blockID < 1
        return
    end

    if isscalar(data)
        value = scalarMean(data);
        return
    end

    dims = size(data);
    if isvector(data)
        if numel(data) >= blockID
            value = scalarMean(data(blockID));
        end
    elseif ndims(data) == 2
        if dims(2) >= blockID
            value = scalarMean(data(:, blockID));
        elseif dims(1) >= blockID
            value = scalarMean(data(blockID, :));
        end
    elseif ndims(data) == 3
        if dims(3) >= blockID
            value = scalarMean(data(:, :, blockID));
        elseif dims(2) >= blockID
            value = scalarMean(data(:, blockID, :));
        end
    end
end

function value = scalarMean(x)
    x = double(x(:));
    x = x(isfinite(x));
    if isempty(x)
        value = NaN;
    else
        value = mean(x);
    end
end

function result = analyzeStimulationOutcomes(X, yBias, yMask, ...
    predictorNames, familyIndex, familyNames, nBootstrap, rngSeed)
    predictorFamilies = mapPredictorFamilies(predictorNames, familyIndex, ...
        familyNames);
    rawDiagnostics = buildStimulationPredictorDiagnostics( ...
        X, predictorNames, predictorFamilies);
    validPredictor = rawDiagnostics.NFiniteRows >= 5 & ...
        rawDiagnostics.NUniqueFiniteValues >= 2 & ...
        isfinite(rawDiagnostics.Std) & rawDiagnostics.Std > 0;
    [familyScores, familyScoreNames, retainedPredictors] = ...
        buildFamilyScorePredictors(X, predictorNames, familyIndex, ...
        familyNames, validPredictor);

    result = struct();
    result.bias = analyzeStimulationSet( ...
        familyScores, yBias, familyScoreNames, rawDiagnostics, ...
        retainedPredictors, nBootstrap, rngSeed);
    result.mask = analyzeStimulationSet( ...
        familyScores, yMask, familyScoreNames, rawDiagnostics, ...
        retainedPredictors, nBootstrap, rngSeed + 43);
end

function byCluster = analyzeStimulationByCluster(X, yBias, yMask, ...
    predictorNames, familyIndex, familyNames, clusterLabels, ...
    nBootstrap, rngSeed)
    clusterIDs = unique(clusterLabels(:), 'stable');
    byCluster = repmat(struct('clusterID', [], 'n', [], 'bias', [], ...
        'mask', []), numel(clusterIDs), 1);

    for clusterIdx = 1:numel(clusterIDs)
        keep = clusterLabels == clusterIDs(clusterIdx);
        predictorFamilies = mapPredictorFamilies(predictorNames, ...
            familyIndex, familyNames);
        rawDiagnostics = buildStimulationPredictorDiagnostics( ...
            X(keep, :), predictorNames, predictorFamilies);
        validPredictor = rawDiagnostics.NFiniteRows >= 5 & ...
            rawDiagnostics.NUniqueFiniteValues >= 2 & ...
            isfinite(rawDiagnostics.Std) & rawDiagnostics.Std > 0;
        [familyScores, familyScoreNames, retainedPredictors] = ...
            buildFamilyScorePredictors(X(keep, :), predictorNames, ...
            familyIndex, familyNames, validPredictor);

        byCluster(clusterIdx).clusterID = clusterIDs(clusterIdx);
        byCluster(clusterIdx).n = nnz(keep);
        byCluster(clusterIdx).bias = analyzeStimulationSet( ...
            familyScores, yBias(keep), familyScoreNames, rawDiagnostics, ...
            retainedPredictors, nBootstrap, rngSeed + 100 * clusterIdx);
        byCluster(clusterIdx).mask = analyzeStimulationSet( ...
            familyScores, yMask(keep), familyScoreNames, rawDiagnostics, ...
            retainedPredictors, nBootstrap, rngSeed + 100 * clusterIdx + 29);
    end
end

function result = analyzeStimulationSet(familyScores, y, familyScoreNames, ...
    rawDiagnostics, retainedPredictors, nBootstrap, rngSeed)
    result = analyzePredictorSet(familyScores, y, familyScoreNames, ...
        nBootstrap, rngSeed);
    result.familyNames = familyScoreNames;
    result.familyIndex = 1:numel(familyScoreNames);
    result.retainedRawPredictors = retainedPredictors;
    result.rawPredictorDiagnostics = rawDiagnostics;
    result.groupedShapleyR2 = result.shapleyR2;
    result.groupedPercentOfModelR2 = result.percentOfModelR2;
    result.groupedBootstrapPercentCI = result.bootstrapPercentCI;
end

function predictorFamilies = mapPredictorFamilies(predictorNames, ...
    familyIndex, familyNames)
    predictorFamilies = cell(size(predictorNames));
    for predictorIdx = 1:numel(predictorNames)
        predictorFamilies{predictorIdx} = familyNames{familyIndex(predictorIdx)};
    end
end

function diagnostics = buildStimulationPredictorDiagnostics( ...
    values, predictorNames, predictorFamilies)
    nPredictors = numel(predictorNames);
    nFiniteRows = zeros(nPredictors, 1);
    nUniqueFiniteValues = zeros(nPredictors, 1);
    minValue = nan(nPredictors, 1);
    maxValue = nan(nPredictors, 1);
    stdValue = nan(nPredictors, 1);

    for predictorIdx = 1:nPredictors
        finiteValues = values(:, predictorIdx);
        finiteValues = finiteValues(isfinite(finiteValues));
        nFiniteRows(predictorIdx) = numel(finiteValues);
        nUniqueFiniteValues(predictorIdx) = numel(unique(finiteValues));
        if ~isempty(finiteValues)
            minValue(predictorIdx) = min(finiteValues);
            maxValue(predictorIdx) = max(finiteValues);
            stdValue(predictorIdx) = std(finiteValues, 0);
        end
    end

    diagnostics = table(predictorNames(:), predictorFamilies(:), ...
        nFiniteRows, nUniqueFiniteValues, minValue, maxValue, stdValue, ...
        'VariableNames', {'PredictorName', 'Family', 'NFiniteRows', ...
        'NUniqueFiniteValues', 'Min', 'Max', 'Std'});
end

function [familyScores, familyScoreNames, retainedPredictors] = ...
    buildFamilyScorePredictors(X, predictorNames, familyIndex, ...
    familyNames, validPredictor)
    familyScores = nan(size(X, 1), 0);
    familyScoreNames = {};
    retainedPredictors = struct('family', {}, 'predictorNames', {});

    for familyIdx = 1:numel(familyNames)
        familyPredictors = find(validPredictor(:)' & familyIndex == familyIdx);
        if isempty(familyPredictors)
            continue
        end

        zValues = nan(size(X, 1), numel(familyPredictors));
        for localIdx = 1:numel(familyPredictors)
            predictorValues = X(:, familyPredictors(localIdx));
            finiteRows = isfinite(predictorValues);
            mu = mean(predictorValues(finiteRows));
            sigma = std(predictorValues(finiteRows), 0);
            zValues(finiteRows, localIdx) = ...
                (predictorValues(finiteRows) - mu) ./ sigma;
        end

        score = mean(zValues, 2, 'omitnan');
        scoreDiagnostics = buildStimulationPredictorDiagnostics( ...
            score, familyNames(familyIdx), familyNames(familyIdx));
        if scoreDiagnostics.NFiniteRows < 5 || ...
                scoreDiagnostics.NUniqueFiniteValues < 2 || ...
                ~isfinite(scoreDiagnostics.Std) || ...
                scoreDiagnostics.Std == 0
            continue
        end

        familyScores(:, end + 1) = score; %#ok<AGROW>
        familyScoreNames{end + 1} = familyNames{familyIdx}; %#ok<AGROW>
        retainedPredictors(end + 1).family = familyNames{familyIdx}; %#ok<AGROW>
        retainedPredictors(end).predictorNames = predictorNames(familyPredictors);
    end
end

function shapleyValues = shapleyR2Grouped(X, y, familyIndex)
    familyIDs = unique(familyIndex, 'stable');
    nGroups = numel(familyIDs);
    shapleyValues = nan(1, nGroups);
    if isempty(X) || nGroups == 0 || size(X, 1) < 2
        return
    end

    shapleyValues(:) = 0;
    factorialGroups = factorial(nGroups);
    r2ByMask = nan(2^nGroups, 1);
    r2ByMask(1) = 0;

    for mask = 1:(2^nGroups - 1)
        useGroup = logical(bitget(mask, 1:nGroups));
        usePredictor = ismember(familyIndex, familyIDs(useGroup));
        fitStats = ordinaryR2FromOLS(X(:, usePredictor), y);
        r2ByMask(mask + 1) = fitStats.R2;
    end

    for groupIdx = 1:nGroups
        others = setdiff(1:nGroups, groupIdx);
        for mask = 0:(2^numel(others) - 1)
            subset = false(1, nGroups);
            subset(others) = logical(bitget(mask, 1:numel(others)));
            subsetSize = nnz(subset);
            withGroup = subset;
            withGroup(groupIdx) = true;
            weight = factorial(subsetSize) * ...
                factorial(nGroups - subsetSize - 1) / factorialGroups;
            shapleyValues(groupIdx) = shapleyValues(groupIdx) + ...
                weight * (r2ForSubset(r2ByMask, withGroup) - ...
                r2ForSubset(r2ByMask, subset));
        end
    end
end

function figureHandles = makeRelevanceFigures(relevanceStruct, opts)
    figureHandles = gobjects(0);
    figureHandles(end + 1) = plotParameterShapleyFigure( ...
        relevanceStruct, opts);
    figureHandles(end + 1) = plotStimulationGroupedShapleyFigure( ...
        relevanceStruct, opts);
    figureHandles(end + 1) = plotParameterScatterFigure( ...
        relevanceStruct, opts, 1);
    figureHandles(end + 1) = plotParameterScatterFigure( ...
        relevanceStruct, opts, 2);
end

function shapleyTable = buildParameterShapleySummaryTable(allBlocks, ...
    predictorNames)
    outcome = {'DeltaBias'; 'DeltaMask'};
    n = [allBlocks.bias.n; allBlocks.mask.n];
    status = {allBlocks.bias.status; allBlocks.mask.status};
    totalInSampleR2 = [allBlocks.bias.OLS.R2; allBlocks.mask.OLS.R2];
    leaveOneOutPredictiveR2 = [ ...
        allBlocks.bias.leaveOneOutPredictiveR2; ...
        allBlocks.mask.leaveOneOutPredictiveR2];
    percentValues = [allBlocks.bias.percentOfModelR2; ...
        allBlocks.mask.percentOfModelR2];
    absoluteValues = [allBlocks.bias.shapleyR2; ...
        allBlocks.mask.shapleyR2];
    rowSumPercent = rowSumOmitNan(percentValues);
    rowSumShapleyR2 = rowSumOmitNan(absoluteValues);

    shapleyTable = table(outcome, n, status, totalInSampleR2, ...
        leaveOneOutPredictiveR2, rowSumShapleyR2, rowSumPercent, ...
        'VariableNames', {'Outcome', 'N', 'ModelStatus', ...
        'TotalInSampleR2', 'LeaveOneOutPredictiveR2', ...
        'RowSumShapleyR2', 'RowSumPercent'});
    for predictorIdx = 1:numel(predictorNames)
        absoluteName = matlab.lang.makeValidName( ...
            [predictorNames{predictorIdx}, 'ShapleyR2']);
        percentName = matlab.lang.makeValidName( ...
            [predictorNames{predictorIdx}, 'PercentModelR2']);
        ciLowName = matlab.lang.makeValidName( ...
            [predictorNames{predictorIdx}, 'PercentModelR2CI025']);
        ciHighName = matlab.lang.makeValidName( ...
            [predictorNames{predictorIdx}, 'PercentModelR2CI975']);
        shapleyTable.(absoluteName) = absoluteValues(:, predictorIdx);
        shapleyTable.(percentName) = percentValues(:, predictorIdx);
        shapleyTable.(ciLowName) = [ ...
            allBlocks.bias.bootstrapPercentCI(predictorIdx, 1); ...
            allBlocks.mask.bootstrapPercentCI(predictorIdx, 1)];
        shapleyTable.(ciHighName) = [ ...
            allBlocks.bias.bootstrapPercentCI(predictorIdx, 2); ...
            allBlocks.mask.bootstrapPercentCI(predictorIdx, 2)];
    end
end

function shapleyTable = buildStimulationFamilyShapleySummaryTable(allBlocks)
    familyNames = allBlocks.bias.familyNames;
    outcome = {'DeltaBias'; 'DeltaMask'};
    n = [allBlocks.bias.n; allBlocks.mask.n];
    status = {allBlocks.bias.status; allBlocks.mask.status};
    totalInSampleR2 = [allBlocks.bias.OLS.R2; allBlocks.mask.OLS.R2];
    leaveOneOutPredictiveR2 = [ ...
        allBlocks.bias.leaveOneOutPredictiveR2; ...
        allBlocks.mask.leaveOneOutPredictiveR2];
    percentValues = [allBlocks.bias.groupedPercentOfModelR2; ...
        allBlocks.mask.groupedPercentOfModelR2];
    absoluteValues = [allBlocks.bias.groupedShapleyR2; ...
        allBlocks.mask.groupedShapleyR2];
    rowSumPercent = rowSumOmitNan(percentValues);
    rowSumShapleyR2 = rowSumOmitNan(absoluteValues);

    shapleyTable = table(outcome, n, status, totalInSampleR2, ...
        leaveOneOutPredictiveR2, rowSumShapleyR2, rowSumPercent, ...
        'VariableNames', {'Outcome', 'N', 'ModelStatus', ...
        'TotalInSampleR2', 'LeaveOneOutPredictiveR2', ...
        'RowSumShapleyR2', 'RowSumPercent'});
    for familyIdx = 1:numel(familyNames)
        absoluteName = matlab.lang.makeValidName( ...
            [familyNames{familyIdx}, 'ShapleyR2']);
        percentName = matlab.lang.makeValidName( ...
            [familyNames{familyIdx}, 'PercentModelR2']);
        ciLowName = matlab.lang.makeValidName( ...
            [familyNames{familyIdx}, 'PercentModelR2CI025']);
        ciHighName = matlab.lang.makeValidName( ...
            [familyNames{familyIdx}, 'PercentModelR2CI975']);
        shapleyTable.(absoluteName) = absoluteValues(:, familyIdx);
        shapleyTable.(percentName) = percentValues(:, familyIdx);
        shapleyTable.(ciLowName) = [ ...
            allBlocks.bias.groupedBootstrapPercentCI(familyIdx, 1); ...
            allBlocks.mask.groupedBootstrapPercentCI(familyIdx, 1)];
        shapleyTable.(ciHighName) = [ ...
            allBlocks.bias.groupedBootstrapPercentCI(familyIdx, 2); ...
            allBlocks.mask.groupedBootstrapPercentCI(familyIdx, 2)];
    end
end

function fig = plotParameterScatterFigure(relevanceStruct, opts, outcomeIdx)
    parameter = relevanceStruct.parameter;
    clusterLabels = relevanceStruct.powerEffectCluster;
    clusterIDs = unique(clusterLabels(:), 'stable');
    colors = lines(max(1, numel(clusterIDs)));
    predictorLabels = {'\DeltaA', '\DeltaB', '\Delta\alpha', ...
        '\Delta\beta'};
    xLimits = [-0.15, 0.05; -0.15, 0.30; -10, 10; -4, 4];
    yLimits = [-5, 50; -2, 16];
    outcomes = {'Biasing', 'Masking'};
    outcomeFields = {'bias', 'mask'};
    yLabels = {'\DeltaBias (% correct)', '\DeltaMask (% correct)'};
    XCell = {parameter.biasPredictors, parameter.maskPredictors};
    yCell = {relevanceStruct.outcomes.deltaBias, ...
        relevanceStruct.outcomes.deltaMask};

    fig = figure('Color', 'w', ...
        'Name', ['Parameter-' lower(outcomes{outcomeIdx}) ...
        ' scatter diagnostics']);
    set(fig, 'Position', [100, 100, 720, 680]);

    for predictorIdx = 1:numel(predictorLabels)
        ax = subplot(2, 2, predictorIdx, 'Parent', fig);
        hold(ax, 'on');
        xline(ax, 0, '-', 'Color', 0.70 .* [1 1 1], ...
            'LineWidth', 1.0, 'HandleVisibility', 'off');
        yline(ax, 0, '-', 'Color', 0.70 .* [1 1 1], ...
            'LineWidth', 1.0, 'HandleVisibility', 'off');
        for clusterIdx = 1:numel(clusterIDs)
            keep = clusterLabels == clusterIDs(clusterIdx);
            scatter(ax, XCell{outcomeIdx}(keep, predictorIdx), ...
                yCell{outcomeIdx}(keep), 36, ...
                colors(clusterIdx, :), 's', 'filled', ...
                'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.80);
        end
        addVisualOLSLine(ax, XCell{outcomeIdx}(:, predictorIdx), ...
            yCell{outcomeIdx}, xLimits(predictorIdx, :));
        pairwise = relevanceStruct.parameter.allBlocks.( ...
            outcomeFields{outcomeIdx}).pairwise(predictorIdx);
        title(ax, predictorLabels{predictorIdx}, 'FontName', 'Arial', ...
            'FontSize', 13, 'FontWeight', 'normal');
        xlabel(ax, predictorLabels{predictorIdx}, 'FontName', 'Arial', ...
            'FontSize', 13);
        ylabel(ax, yLabels{outcomeIdx}, 'FontName', 'Arial', ...
            'FontSize', 13);
        text(ax, 0.04, 0.94, sprintf('n=%d\nrho=%.2f\np=%.3g', ...
            pairwise.n, pairwise.spearmanRho, pairwise.spearmanP), ...
            'Units', 'normalized', 'VerticalAlignment', 'top', ...
            'FontName', 'Arial', 'FontSize', 11);
        xlim(ax, xLimits(predictorIdx, :));
        ylim(ax, yLimits(outcomeIdx, :));
        pbaspect(ax, [1 1 1]);
        box(ax, 'off');
        set(ax, 'FontName', 'Arial', 'FontSize', 12, ...
            'LineWidth', 1.0);
    end

    annotation(fig, 'textbox', [0, 0.955, 1, 0.04], ...
        'String', sprintf('%s %s parameter-%s scatter diagnostics', ...
        opts.monkeyName, opts.chamberWanted, lower(outcomes{outcomeIdx})), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'EdgeColor', 'none', 'FontName', 'Arial', 'FontSize', 16, ...
        'FontWeight', 'bold');
end

function fig = plotParameterShapleyFigure(relevanceStruct, opts)
    values = [
        relevanceStruct.parameter.allBlocks.bias.percentOfModelR2;
        relevanceStruct.parameter.allBlocks.mask.percentOfModelR2];
    rowSums = rowSumOmitNan(values);
    looBias = relevanceStruct.parameter.allBlocks.bias.leaveOneOutPredictiveR2;
    looMask = relevanceStruct.parameter.allBlocks.mask.leaveOneOutPredictiveR2;
    fprintf(['Parameter Shapley row-sum sanity check | ' ...
        'DeltaBias=%.3f%% | DeltaMask=%.3f%% | ' ...
        'LOO predictive R2: bias=%.3f, mask=%.3f\n'], ...
        rowSums(1), rowSums(2), looBias, looMask);

    fig = figure('Color', 'w', 'Name', 'Parameter Shapley relevance');
    set(fig, 'Position', [100, 100, 520, 380]);
    ax = axes('Parent', fig);
    imagesc(ax, values);
    axis(ax, 'image');
    caxis(ax, [0, 100]);
    colormap(parula);
    cb = colorbar(ax);
    ylabel(cb, '% of model R^2', 'FontName', 'Arial', 'FontSize', 13);
    xticks(ax, 1:numel(relevanceStruct.parameter.predictorNames));
    xticklabels(ax, {'A', 'B', '\alpha', '\beta'});
    yticks(ax, 1:2);
    yticklabels(ax, {'\DeltaBias', '\DeltaMask'});
    set(ax, 'FontName', 'Arial', 'FontSize', 13, 'LineWidth', 1.0);
    title(ax, sprintf(['Parameter Shapley relevance | n=%d | ' ...
        'in-sample R^2: bias=%.2f, mask=%.2f | ' ...
        'LOO R^2: bias=%.2f, mask=%.2f'], ...
        relevanceStruct.parameter.allBlocks.bias.n, ...
        relevanceStruct.parameter.allBlocks.bias.OLS.R2, ...
        relevanceStruct.parameter.allBlocks.mask.OLS.R2, ...
        looBias, looMask), ...
        'FontName', 'Arial', 'FontSize', 13, 'FontWeight', 'bold');
    addHeatmapText(values);
end

function fig = plotStimulationGroupedShapleyFigure(relevanceStruct, opts)
    familyNames = relevanceStruct.stimulation.allBlocks.bias.familyNames;
    values = [
        relevanceStruct.stimulation.allBlocks.bias.groupedPercentOfModelR2;
        relevanceStruct.stimulation.allBlocks.mask.groupedPercentOfModelR2];
    totalR2 = [relevanceStruct.stimulation.allBlocks.bias.OLS.R2, ...
        relevanceStruct.stimulation.allBlocks.mask.OLS.R2];
    unavailable = isempty(familyNames) || all(~isfinite(values(:))) || ...
        any(~isfinite(totalR2));

    fig = figure('Color', 'w', 'Name', 'Stimulation-family Shapley relevance');
    set(fig, 'Position', [100, 100, 520, 380]);
    ax = axes('Parent', fig);
    if unavailable
        axis(ax, 'off');
        text(ax, 0.5, 0.5, ['Stimulation-family model unavailable: ' ...
            'insufficient complete finite predictors.'], ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontName', 'Arial', 'FontSize', 13);
        return
    end

    imagesc(ax, values);
    axis(ax, 'image');
    caxis(ax, [0, 100]);
    colormap(parula);
    cb = colorbar(ax);
    ylabel(cb, '% of model R^2', 'FontName', 'Arial', 'FontSize', 13);
    xticks(ax, 1:numel(familyNames));
    xticklabels(ax, familyNames);
    yticks(ax, 1:2);
    yticklabels(ax, {'\DeltaBias', '\DeltaMask'});
    set(ax, 'FontName', 'Arial', 'FontSize', 13, 'LineWidth', 1.0);
    title(ax, sprintf(['Stimulation-family Shapley relevance | n=%d | ' ...
        'family-score in-sample R^2: bias=%.2f, mask=%.2f'], ...
        relevanceStruct.stimulation.allBlocks.bias.n, ...
        relevanceStruct.stimulation.allBlocks.bias.OLS.R2, ...
        relevanceStruct.stimulation.allBlocks.mask.OLS.R2), ...
        'FontName', 'Arial', 'FontSize', 13, 'FontWeight', 'bold');
    addHeatmapText(values);
end

function addHeatmapText(values)
    for rowIdx = 1:size(values, 1)
        for colIdx = 1:size(values, 2)
            if isfinite(values(rowIdx, colIdx))
                text(colIdx, rowIdx, sprintf('%.0f', values(rowIdx, colIdx)), ...
                    'HorizontalAlignment', 'center', ...
                    'Color', 'w', 'FontWeight', 'bold', ...
                    'FontName', 'Arial');
            end
        end
    end
end

function outputPaths = saveRelevanceOutputs(relevanceStruct, opts)
    outputDir = fullfile(opts.mainPath, opts.monkeyName, 'Meta', ...
        'psychometrics', 'psycluster-relevance');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    safeBase = regexprep(sprintf('%s_%s_%s_%s', opts.monkeyName, ...
        opts.chamberWanted, opts.aggregateField, opts.sourceMdlField), ...
        '[^A-Za-z0-9_-]', '_');
    outputPaths = struct();
    outputPaths.mat = fullfile(outputDir, [safeBase, '_relevance.mat']);
    figureHandles = relevanceStruct.figureHandles;
    relevanceStructForSave = relevanceStruct;
    relevanceStructForSave.figureHandles = [];
    relevanceStruct = relevanceStructForSave;
    save(outputPaths.mat, 'relevanceStruct');

    figureNames = {'parameterShapley', 'stimulationGroupedShapley', ...
        'parameterBiasingScatterDiagnostic', ...
        'parameterMaskingScatterDiagnostic'};
    outputPaths.figures = struct([]);
    for figIdx = 1:numel(figureHandles)
        baseName = fullfile(outputDir, ...
            sprintf('%s_%s', safeBase, figureNames{figIdx}));
        print(figureHandles(figIdx), [baseName, '.svg'], ...
            '-dsvg', '-painters');
        print(figureHandles(figIdx), [baseName, '.png'], ...
            '-dpng', '-r300');
        outputPaths.figures(figIdx).svg = [baseName, '.svg'];
        outputPaths.figures(figIdx).png = [baseName, '.png'];
    end
end

function addVisualOLSLine(ax, x, y, xLimits)
    validRows = isfinite(x) & isfinite(y);
    if nnz(validRows) < 2 || numel(unique(x(validRows))) < 2
        return
    end

    coefficients = polyfit(x(validRows), y(validRows), 1);
    xFit = linspace(xLimits(1), xLimits(2), 100);
    yFit = polyval(coefficients, xFit);
    plot(ax, xFit, yFit, '-', 'Color', 0.15 .* [1 1 1], ...
        'LineWidth', 1.5, 'HandleVisibility', 'off');
end
