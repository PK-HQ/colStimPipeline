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
    summaryTables.parameterShapleyByCluster = ...
        buildParameterClusterShapleySummaryTable(parameter.byCluster, ...
        paramNames);
    summaryTables.parameterClusterDiagnostic = ...
        buildParameterClusterDiagnosticSummary(parameter, paramNames);
    summaryTables.stimulationPredictorDiagnostics = stim.diagnostics;
    summaryTables.stimulationFamilyScores = stim.familyScoreTable;
    summaryTables.stimulationFamilyScoreDiagnostics = ...
        stim.familyScoreDiagnostics;
    summaryTables.stimulationFamilyShapley = ...
        buildStimulationFamilyShapleySummaryTable(stimulation.allBlocks);
    summaryTables.stimulationFamilyShapleyByCluster = ...
        buildStimulationClusterShapleySummaryTable(stimulation.byCluster);
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
    fitStats = struct('R2', NaN, 'adjustedR2', NaN, ...
        'coefficients', [], 'intercept', NaN, 'rank', 0);

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
    n = size(Xz, 1);
    p = size(Xz, 2);

    fitStats.R2 = max(0, 1 - residualSS / totalSS);
    if n > p + 1
        fitStats.adjustedR2 = 1 - (1 - fitStats.R2) * ...
            (n - 1) / (n - p - 1);
    end
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
        'nColumns', 'Spatial factors'; ...
        'pixelsON', 'Spatial factors'; ...
        'areaFinalROI', 'Spatial factors'; ...
        'areaPixelsONWithinROI', 'Spatial factors'; ...
        'spatialDutyCycleWithinROI', 'Spatial factors'; ...
        'temporalDutyCycle', 'Temporal factors'; ...
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
    validPredictor = validStimulationPredictorMask(stim.diagnostics);
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
    validPredictor = validStimulationPredictorMask(rawDiagnostics);
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
    validPredictor = validStimulationPredictorMask(rawDiagnostics);
    result.excludedRawPredictors = rawDiagnostics(~validPredictor, :);
    result.groupedShapleyR2 = result.shapleyR2;
    result.groupedPercentOfModelR2 = result.percentOfModelR2;
    result.groupedBootstrapPercentCI = result.bootstrapPercentCI;
end

function validPredictor = validStimulationPredictorMask(diagnostics)
    validPredictor = diagnostics.NFiniteRows >= 5 & ...
        diagnostics.NUniqueFiniteValues >= 2 & ...
        isfinite(diagnostics.Std) & diagnostics.Std > 0;
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
    figureHandles(end + 1) = plotParameterClusterSummaryFigure( ...
        relevanceStruct);
    figureHandles(end + 1) = plotStimulationGroupedShapleyFigure( ...
        relevanceStruct, opts);
    figureHandles(end + 1) = plotParameterScatterFigure( ...
        relevanceStruct, opts, 1);
    figureHandles(end + 1) = plotParameterScatterFigure( ...
        relevanceStruct, opts, 2);
end

function shapleyTable = buildParameterShapleySummaryTable(allBlocks, ...
    predictorNames)
    shapleyTable = buildShapleySummaryRows([], allBlocks.bias, ...
        allBlocks.mask, predictorNames, false);
end

function shapleyTable = buildParameterClusterShapleySummaryTable( ...
    byCluster, predictorNames)
    shapleyTable = table();
    for clusterIdx = 1:numel(byCluster)
        clusterTable = buildShapleySummaryRows( ...
            byCluster(clusterIdx).clusterID, byCluster(clusterIdx).bias, ...
            byCluster(clusterIdx).mask, predictorNames, false);
        shapleyTable = [shapleyTable; clusterTable]; %#ok<AGROW>
    end
end

function summaryTable = buildParameterClusterDiagnosticSummary( ...
    parameter, predictorNames)
    p = numel(predictorNames);
    rowLabels = {'All pooled'; 'Power cluster 1'; 'Power cluster 2'; ...
        'Power cluster 3'};
    clusterIDs = [NaN; 1; 2; 3];
    nRows = numel(rowLabels);

    rowStatus = cell(nRows, 1);
    n = nan(nRows, 1);
    biasInSampleR2 = nan(nRows, 1);
    biasPredictiveR2 = nan(nRows, 1);
    biasTopParameter = cell(nRows, 1);
    biasTopParameterPercent = nan(nRows, 1);
    maskInSampleR2 = nan(nRows, 1);
    maskPredictiveR2 = nan(nRows, 1);
    maskTopParameter = cell(nRows, 1);
    maskTopParameterPercent = nan(nRows, 1);

    clusterIDValues = [];
    if ~isempty(parameter.byCluster)
        clusterIDValues = [parameter.byCluster.clusterID];
    end

    for rowIdx = 1:nRows
        if rowIdx == 1
            biasStats = parameter.allBlocks.bias;
            maskStats = parameter.allBlocks.mask;
        else
            matchIdx = find(clusterIDValues == clusterIDs(rowIdx), 1);
            if isempty(matchIdx)
                biasStats = missingParameterOutcomeStats(p);
                maskStats = missingParameterOutcomeStats(p);
            else
                biasStats = parameter.byCluster(matchIdx).bias;
                maskStats = parameter.byCluster(matchIdx).mask;
            end
        end

        n(rowIdx) = min([biasStats.n, maskStats.n]);
        rowStatus{rowIdx} = 'ok';
        if parameterClusterRowUnderpowered(biasStats, maskStats, p)
            rowStatus{rowIdx} = 'underpowered';
        end

        [biasInSampleR2(rowIdx), biasPredictiveR2(rowIdx), ...
            biasTopParameter{rowIdx}, biasTopParameterPercent(rowIdx)] = ...
            parameterOutcomeDiagnosticValues(biasStats, predictorNames);
        [maskInSampleR2(rowIdx), maskPredictiveR2(rowIdx), ...
            maskTopParameter{rowIdx}, maskTopParameterPercent(rowIdx)] = ...
            parameterOutcomeDiagnosticValues(maskStats, predictorNames);
    end

    summaryTable = table(rowLabels, clusterIDs, rowStatus, n, ...
        biasInSampleR2, biasPredictiveR2, biasTopParameter, ...
        biasTopParameterPercent, maskInSampleR2, maskPredictiveR2, ...
        maskTopParameter, maskTopParameterPercent, 'VariableNames', ...
        {'RowLabel', 'ClusterID', 'RowStatus', 'N', ...
        'BiasInSampleR2', 'BiasPredictiveR2', 'BiasTopParameter', ...
        'BiasTopParameterPercent', 'MaskInSampleR2', ...
        'MaskPredictiveR2', 'MaskTopParameter', ...
        'MaskTopParameterPercent'});
end

function stats = missingParameterOutcomeStats(p)
    stats = struct();
    stats.n = 0;
    stats.status = 'underpowered';
    stats.OLS = struct('R2', NaN, 'adjustedR2', NaN);
    stats.leaveOneOutPredictiveR2 = NaN;
    stats.percentOfModelR2 = nan(1, p);
end

function isUnderpowered = parameterClusterRowUnderpowered( ...
    biasStats, maskStats, p)
    isUnderpowered = parameterOutcomeUnderpowered(biasStats, p) || ...
        parameterOutcomeUnderpowered(maskStats, p);
end

function isUnderpowered = parameterOutcomeUnderpowered(stats, p)
    isUnderpowered = ~isfinite(stats.n) || stats.n <= p + 1 || ...
        ~strcmp(stats.status, 'ok');
end

function [inSampleR2, predictiveR2, topParameter, topPercent] = ...
    parameterOutcomeDiagnosticValues(stats, predictorNames)
    inSampleR2 = stats.OLS.R2;
    predictiveR2 = stats.leaveOneOutPredictiveR2;
    topParameter = 'n/a';
    topPercent = NaN;

    percentValues = stats.percentOfModelR2;
    finiteIdx = find(isfinite(percentValues));
    if isempty(finiteIdx)
        return
    end

    [topPercent, localIdx] = max(percentValues(finiteIdx));
    topIdx = finiteIdx(localIdx);
    topParameter = compactParameterLabel(predictorNames{topIdx});
end

function label = compactParameterLabel(parameterName)
    switch parameterName
        case 'A'
            label = '\DeltaA';
        case 'B'
            label = '\DeltaB';
        case 'alpha'
            label = '\Delta\alpha';
        case 'beta'
            label = '\Delta\beta';
        otherwise
            label = ['\Delta', parameterName];
    end
end

function shapleyTable = buildStimulationFamilyShapleySummaryTable(allBlocks)
    familyNames = canonicalStimulationFamilyNames();
    shapleyTable = buildShapleySummaryRows([], allBlocks.bias, ...
        allBlocks.mask, familyNames, true);
end

function shapleyTable = buildStimulationClusterShapleySummaryTable(byCluster)
    shapleyTable = table();
    familyNames = canonicalStimulationFamilyNames();
    for clusterIdx = 1:numel(byCluster)
        clusterTable = buildShapleySummaryRows( ...
            byCluster(clusterIdx).clusterID, byCluster(clusterIdx).bias, ...
            byCluster(clusterIdx).mask, familyNames, true);
        shapleyTable = [shapleyTable; clusterTable]; %#ok<AGROW>
    end
end

function familyNames = canonicalStimulationFamilyNames()
    familyNames = {'Energy', 'Spatial factors', 'Temporal factors', 'QC'};
end

function shapleyTable = buildShapleySummaryRows(clusterID, biasStats, ...
    maskStats, predictorNames, useGroupedFields)
    outcome = {'DeltaBias'; 'DeltaMask'};
    n = [biasStats.n; maskStats.n];
    status = {biasStats.status; maskStats.status};
    totalInSampleR2 = [biasStats.OLS.R2; maskStats.OLS.R2];
    adjustedR2 = [biasStats.OLS.adjustedR2; maskStats.OLS.adjustedR2];
    leaveOneOutPredictiveR2 = [ ...
        biasStats.leaveOneOutPredictiveR2; ...
        maskStats.leaveOneOutPredictiveR2];
    if useGroupedFields
        sourceNames = biasStats.familyNames;
        sourcePercentValues = [biasStats.groupedPercentOfModelR2; ...
            maskStats.groupedPercentOfModelR2];
        sourceAbsoluteValues = [biasStats.groupedShapleyR2; ...
            maskStats.groupedShapleyR2];
        sourceCiBias = biasStats.groupedBootstrapPercentCI;
        sourceCiMask = maskStats.groupedBootstrapPercentCI;
    else
        sourceNames = biasStats.predictorNames;
        sourcePercentValues = [biasStats.percentOfModelR2; ...
            maskStats.percentOfModelR2];
        sourceAbsoluteValues = [biasStats.shapleyR2; maskStats.shapleyR2];
        sourceCiBias = biasStats.bootstrapPercentCI;
        sourceCiMask = maskStats.bootstrapPercentCI;
    end

    nPredictors = numel(predictorNames);
    percentValues = nan(2, nPredictors);
    absoluteValues = nan(2, nPredictors);
    ciLow = nan(2, nPredictors);
    ciHigh = nan(2, nPredictors);
    for predictorIdx = 1:nPredictors
        sourceIdx = find(strcmp(sourceNames, predictorNames{predictorIdx}), 1);
        if isempty(sourceIdx)
            continue;
        end
        if sourceIdx <= size(sourcePercentValues, 2)
            percentValues(:, predictorIdx) = sourcePercentValues(:, sourceIdx);
        end
        if sourceIdx <= size(sourceAbsoluteValues, 2)
            absoluteValues(:, predictorIdx) = sourceAbsoluteValues(:, sourceIdx);
        end
        if sourceIdx <= size(sourceCiBias, 1)
            ciLow(1, predictorIdx) = sourceCiBias(sourceIdx, 1);
            ciHigh(1, predictorIdx) = sourceCiBias(sourceIdx, 2);
        end
        if sourceIdx <= size(sourceCiMask, 1)
            ciLow(2, predictorIdx) = sourceCiMask(sourceIdx, 1);
            ciHigh(2, predictorIdx) = sourceCiMask(sourceIdx, 2);
        end
    end

    approximateAbsoluteFromPercent = ...
        repmat(totalInSampleR2, 1, nPredictors) .* percentValues ./ 100;
    rowSumPercent = rowSumOmitNan(percentValues);
    rowSumShapleyR2 = rowSumOmitNan(absoluteValues);

    if isempty(clusterID)
        shapleyTable = table(outcome, n, status, totalInSampleR2, ...
            adjustedR2, leaveOneOutPredictiveR2, rowSumShapleyR2, ...
            rowSumPercent, 'VariableNames', {'Outcome', 'N', ...
            'ModelStatus', 'TotalInSampleR2', 'AdjustedR2', ...
            'LeaveOneOutPredictiveR2', 'RowSumShapleyR2', ...
            'RowSumPercent'});
    else
        clusterIDColumn = repmat(clusterID, 2, 1);
        shapleyTable = table(clusterIDColumn, outcome, n, status, ...
            totalInSampleR2, adjustedR2, leaveOneOutPredictiveR2, ...
            rowSumShapleyR2, rowSumPercent, 'VariableNames', ...
            {'ClusterID', 'Outcome', 'N', 'ModelStatus', ...
            'TotalInSampleR2', 'AdjustedR2', ...
            'LeaveOneOutPredictiveR2', 'RowSumShapleyR2', ...
            'RowSumPercent'});
    end

    for predictorIdx = 1:numel(predictorNames)
        absoluteName = matlab.lang.makeValidName( ...
            [predictorNames{predictorIdx}, 'ShapleyR2']);
        approximateName = matlab.lang.makeValidName( ...
            [predictorNames{predictorIdx}, 'ApproxContributionR2']);
        percentName = matlab.lang.makeValidName( ...
            [predictorNames{predictorIdx}, 'PercentModelR2']);
        ciLowName = matlab.lang.makeValidName( ...
            [predictorNames{predictorIdx}, 'PercentModelR2CI025']);
        ciHighName = matlab.lang.makeValidName( ...
            [predictorNames{predictorIdx}, 'PercentModelR2CI975']);
        shapleyTable.(absoluteName) = absoluteValues(:, predictorIdx);
        shapleyTable.(approximateName) = ...
            approximateAbsoluteFromPercent(:, predictorIdx);
        shapleyTable.(percentName) = percentValues(:, predictorIdx);
        shapleyTable.(ciLowName) = ciLow(:, predictorIdx);
        shapleyTable.(ciHighName) = ciHigh(:, predictorIdx);
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
    totalR2 = [relevanceStruct.parameter.allBlocks.bias.OLS.R2; ...
        relevanceStruct.parameter.allBlocks.mask.OLS.R2];
    predictiveR2 = [ ...
        relevanceStruct.parameter.allBlocks.bias.leaveOneOutPredictiveR2; ...
        relevanceStruct.parameter.allBlocks.mask.leaveOneOutPredictiveR2];
    rowSums = rowSumOmitNan(values);
    fprintf(['Parameter Shapley row-sum sanity check | ' ...
        'DeltaBias=%.3f%% | DeltaMask=%.3f%% | ' ...
        'predictive R2: bias=%.3f, mask=%.3f\n'], ...
        rowSums(1), rowSums(2), predictiveR2(1), predictiveR2(2));

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
        'in-sample R^2: bias=%.2f, mask=%.2f'], ...
        relevanceStruct.parameter.allBlocks.bias.n, totalR2(1), ...
        totalR2(2)), 'FontName', 'Arial', 'FontSize', 13, ...
        'FontWeight', 'bold');
    xlabel(ax, sprintf('predictive R^2: bias=%.2f, mask=%.2f', ...
        predictiveR2(1), predictiveR2(2)), 'FontName', 'Arial', ...
        'FontSize', 10);
    addHeatmapText(values, totalR2, predictiveR2, true);
end

function fig = plotParameterClusterSummaryFigure(relevanceStruct)
    summaryTable = relevanceStruct.summaryTables.parameterClusterDiagnostic;
    headers = {'', 'n', sprintf('\\DeltaBias\nin-sample R^2'), ...
        sprintf('\\DeltaBias\npredictive R^2'), ...
        sprintf('\\DeltaBias\ntop param'), ...
        sprintf('\\DeltaBias\ntop param %%'), ...
        sprintf('\\DeltaMask\nin-sample R^2'), ...
        sprintf('\\DeltaMask\npredictive R^2'), ...
        sprintf('\\DeltaMask\ntop param'), ...
        sprintf('\\DeltaMask\ntop param %%')};
    colWidths = [1.35, 0.45, 0.95, 0.95, 0.78, 0.88, ...
        0.95, 0.95, 0.78, 0.88];
    xEdges = [0, cumsum(colWidths)];
    tableWidth = xEdges(end);
    nRows = height(summaryTable);

    fig = figure('Color', 'w', ...
        'Name', 'Parameter by-power-cluster summary');
    set(fig, 'Position', [100, 100, 980, 340]);
    ax = axes('Parent', fig, 'Position', [0.04, 0.08, 0.92, 0.76]);
    hold(ax, 'on');
    axis(ax, 'off');
    set(ax, 'XLim', [0, tableWidth], 'YLim', [0, nRows + 1], ...
        'YDir', 'reverse');

    title(ax, ['Parameter relevance by power cluster | ' ...
        'diagnostic summary'], 'FontName', 'Arial', 'FontSize', 13, ...
        'FontWeight', 'bold');

    for colIdx = 1:numel(headers)
        xCenter = mean(xEdges(colIdx:(colIdx + 1)));
        text(ax, xCenter, 0.50, headers{colIdx}, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontName', 'Arial', 'FontSize', 9, 'FontWeight', 'bold', ...
            'Interpreter', 'tex');
    end

    for rowIdx = 1:nRows
        yCenter = rowIdx + 0.50;
        rowColor = [0 0 0];
        rowLabel = summaryTable.RowLabel{rowIdx};
        if strcmp(summaryTable.RowStatus{rowIdx}, 'underpowered')
            rowColor = 0.50 .* [1 1 1];
            rowLabel = sprintf('%s\nunderpowered', rowLabel);
        end

        rowValues = {rowLabel, formatIntegerCell(summaryTable.N(rowIdx)), ...
            formatR2Cell(summaryTable.BiasInSampleR2(rowIdx)), ...
            formatR2Cell(summaryTable.BiasPredictiveR2(rowIdx)), ...
            summaryTable.BiasTopParameter{rowIdx}, ...
            formatPercentCell(summaryTable.BiasTopParameterPercent(rowIdx)), ...
            formatR2Cell(summaryTable.MaskInSampleR2(rowIdx)), ...
            formatR2Cell(summaryTable.MaskPredictiveR2(rowIdx)), ...
            summaryTable.MaskTopParameter{rowIdx}, ...
            formatPercentCell(summaryTable.MaskTopParameterPercent(rowIdx))};

        for colIdx = 1:numel(rowValues)
            xCenter = mean(xEdges(colIdx:(colIdx + 1)));
            horizontalAlignment = 'center';
            if colIdx == 1
                xCenter = xEdges(colIdx) + 0.06;
                horizontalAlignment = 'left';
            end
            text(ax, xCenter, yCenter, rowValues{colIdx}, ...
                'HorizontalAlignment', horizontalAlignment, ...
                'VerticalAlignment', 'middle', 'FontName', 'Arial', ...
                'FontSize', 9, 'Color', rowColor, 'Interpreter', 'tex');
        end
    end

    for rowEdge = 0:(nRows + 1)
        line(ax, [0, tableWidth], [rowEdge, rowEdge], ...
            'Color', 0.82 .* [1 1 1], 'LineWidth', 0.75);
    end
    for colEdge = 1:numel(xEdges)
        line(ax, [xEdges(colEdge), xEdges(colEdge)], [0, nRows + 1], ...
            'Color', 0.88 .* [1 1 1], 'LineWidth', 0.75);
    end
end

function textValue = formatIntegerCell(value)
    if isfinite(value)
        textValue = sprintf('%d', value);
    else
        textValue = 'n/a';
    end
end

function textValue = formatR2Cell(value)
    if isfinite(value)
        textValue = sprintf('%.2f', value);
    else
        textValue = 'n/a';
    end
end

function textValue = formatPercentCell(value)
    if isfinite(value)
        textValue = sprintf('%.0f%%', value);
    else
        textValue = 'n/a';
    end
end

function fig = plotStimulationGroupedShapleyFigure(relevanceStruct, opts)
    familyNames = relevanceStruct.stimulation.allBlocks.bias.familyNames;
    values = [
        relevanceStruct.stimulation.allBlocks.bias.groupedPercentOfModelR2;
        relevanceStruct.stimulation.allBlocks.mask.groupedPercentOfModelR2];
    totalR2 = [relevanceStruct.stimulation.allBlocks.bias.OLS.R2; ...
        relevanceStruct.stimulation.allBlocks.mask.OLS.R2];
    unavailable = isempty(familyNames) || all(~isfinite(values(:))) || ...
        any(~isfinite(totalR2));

    fig = figure('Color', 'w', 'Name', 'Stimulation-family Shapley relevance');
    set(fig, 'Position', [100, 100, 580, 500]);
    ax = axes('Parent', fig, 'Position', [0.16, 0.30, 0.62, 0.56]);
    if unavailable
        axis(ax, 'off');
        text(ax, 0.5, 0.5, ['Stimulation-family model unavailable: ' ...
            'insufficient complete finite predictors.'], ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontName', 'Arial', 'FontSize', 13);
        addStimulationPredictorFooter(fig, relevanceStruct);
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
        relevanceStruct.stimulation.allBlocks.bias.n, totalR2(1), ...
        totalR2(2)), 'FontName', 'Arial', 'FontSize', 13, ...
        'FontWeight', 'bold');
    addHeatmapText(values, totalR2, totalR2, false);
    addStimulationPredictorFooter(fig, relevanceStruct);
end

function addHeatmapText(values, totalR2, predictiveR2, requirePredictiveR2)
    approximateR2 = totalR2 .* values ./ 100;
    for rowIdx = 1:size(values, 1)
        textColor = heatmapTextColor(totalR2(rowIdx), ...
            predictiveR2(rowIdx), requirePredictiveR2);
        for colIdx = 1:size(values, 2)
            if isfinite(values(rowIdx, colIdx))
                text(colIdx, rowIdx, sprintf('%.0f%%\nR^2=%.2f', ...
                    values(rowIdx, colIdx), approximateR2(rowIdx, colIdx)), ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', 'Color', textColor, ...
                    'FontWeight', 'bold', 'FontName', 'Arial', ...
                    'FontSize', 10);
            end
        end
    end
end

function textColor = heatmapTextColor(totalR2, predictiveR2, ...
    requirePredictiveR2)
    if requirePredictiveR2
        strongModel = isfinite(totalR2) && totalR2 >= 0.25 && ...
            isfinite(predictiveR2) && predictiveR2 >= 0.10;
    else
        strongModel = isfinite(totalR2) && totalR2 >= 0.25;
    end

    if strongModel
        textColor = [0 0 0];
    else
        textColor = 0.50 .* [1 1 1];
    end
end

function addStimulationPredictorFooter(fig, relevanceStruct)
    [includedLine, excludedLine] = stimulationPredictorFooterText( ...
        relevanceStruct);
    footerText = includedLine;
    if ~isempty(excludedLine)
        footerText = sprintf('%s\n%s', includedLine, excludedLine);
    end
    annotation(fig, 'textbox', [0.05, 0.02, 0.90, 0.16], ...
        'String', footerText, 'EdgeColor', 'none', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
        'FontName', 'Arial', 'FontSize', 10, 'Interpreter', 'none');
end

function [includedLine, excludedLine] = stimulationPredictorFooterText( ...
    relevanceStruct)
    families = {'Energy', 'Spatial factors', 'Temporal factors', 'QC'};
    retained = relevanceStruct.stimulation.allBlocks.bias.retainedRawPredictors;
    includedParts = cell(1, numel(families));
    for familyIdx = 1:numel(families)
        predictorList = {};
        for retainedIdx = 1:numel(retained)
            if strcmp(retained(retainedIdx).family, families{familyIdx})
                predictorList = retained(retainedIdx).predictorNames;
                break
            end
        end
        if isempty(predictorList)
            predictorText = 'none';
        else
            predictorText = strjoin(predictorList, ', ');
        end
        includedParts{familyIdx} = sprintf('%s = %s', ...
            families{familyIdx}, predictorText);
    end
    includedLine = strjoin(includedParts, ' | ');

    diagnostics = relevanceStruct.stimulation.allBlocks.bias.rawPredictorDiagnostics;
    validPredictor = validStimulationPredictorMask(diagnostics);
    excludedNames = diagnostics.PredictorName(~validPredictor);
    if isempty(excludedNames)
        excludedLine = '';
    else
        excludedLine = ['Excluded from family scores = ', ...
            strjoin(excludedNames(:)', ', ')];
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

    figureNames = {'parameterShapley', 'parameterClusterSummary', ...
        'stimulationGroupedShapley', 'parameterBiasingScatterDiagnostic', ...
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
