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
        ['Stimulation-family Shapley values are preferred over individual ' ...
         'stimulation predictors when predictors are collinear.']; ...
        ['Leave-one-out predictive R2 is reported only when the design has ' ...
         'enough rows and rank.']};

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
    if n < 3
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
        'nColumns', 'Coverage'; ...
        'pixelsON', 'Coverage'; ...
        'areaFinalROI', 'Coverage'; ...
        'areaPixelsONWithinROI', 'Coverage'; ...
        'spatialDutyCycleWithinROI', 'Coverage'; ...
        'temporalDutyCycle', 'Timing'; ...
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

    present = any(isfinite(values), 1);
    values = values(:, present);
    predictorNames = fields(present, 1)';
    predictorFamilies = fields(present, 2)';
    familyNames = unique(predictorFamilies, 'stable');
    familyIndex = nan(1, numel(predictorFamilies));
    for familyIdx = 1:numel(familyNames)
        familyIndex(strcmp(predictorFamilies, familyNames{familyIdx})) = ...
            familyIdx;
    end

    stim = struct();
    stim.values = values;
    stim.predictorNames = predictorNames;
    stim.familyNames = familyNames;
    stim.familyIndex = familyIndex;
    if isempty(predictorNames)
        stim.table = table(blockIDs(:), 'VariableNames', {'blockID'});
    else
        stim.table = array2table(values, 'VariableNames', predictorNames);
        stim.table.blockID = blockIDs(:);
    end
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
    result = struct();
    result.bias = analyzeStimulationSet( ...
        X, yBias, predictorNames, familyIndex, familyNames, ...
        nBootstrap, rngSeed);
    result.mask = analyzeStimulationSet( ...
        X, yMask, predictorNames, familyIndex, familyNames, ...
        nBootstrap, rngSeed + 43);
end

function byCluster = analyzeStimulationByCluster(X, yBias, yMask, ...
    predictorNames, familyIndex, familyNames, clusterLabels, ...
    nBootstrap, rngSeed)
    clusterIDs = unique(clusterLabels(:), 'stable');
    byCluster = repmat(struct('clusterID', [], 'n', [], 'bias', [], ...
        'mask', []), numel(clusterIDs), 1);

    for clusterIdx = 1:numel(clusterIDs)
        keep = clusterLabels == clusterIDs(clusterIdx);
        byCluster(clusterIdx).clusterID = clusterIDs(clusterIdx);
        byCluster(clusterIdx).n = nnz(keep);
        byCluster(clusterIdx).bias = analyzeStimulationSet( ...
            X(keep, :), yBias(keep), predictorNames, familyIndex, ...
            familyNames, nBootstrap, rngSeed + 100 * clusterIdx);
        byCluster(clusterIdx).mask = analyzeStimulationSet( ...
            X(keep, :), yMask(keep), predictorNames, familyIndex, ...
            familyNames, nBootstrap, rngSeed + 100 * clusterIdx + 29);
    end
end

function result = analyzeStimulationSet(X, y, predictorNames, ...
    familyIndex, familyNames, nBootstrap, rngSeed)
    usable = any(isfinite(X), 1);
    X = X(:, usable);
    predictorNames = predictorNames(usable);
    familyIndex = familyIndex(usable);
    [familyIndex, familyNames] = compactFamilyIndex(familyIndex, familyNames);

    result = analyzePredictorSet(X, y, predictorNames, nBootstrap, rngSeed);
    result.familyNames = familyNames;
    result.familyIndex = familyIndex;
    validY = y(result.validRows);
    validX = X(result.validRows, :);
    result.groupedShapleyR2 = shapleyR2Grouped(validX, validY, familyIndex);
    result.groupedPercentOfModelR2 = percentContribution( ...
        result.groupedShapleyR2);
    result.groupedBootstrapPercentCI = bootstrapShapleyPercentCI( ...
        validX, validY, nBootstrap, rngSeed + 67, familyIndex);
end

function [newIndex, newNames] = compactFamilyIndex(familyIndex, familyNames)
    familyIDs = unique(familyIndex, 'stable');
    newIndex = nan(size(familyIndex));
    newNames = cell(1, numel(familyIDs));
    for idx = 1:numel(familyIDs)
        newIndex(familyIndex == familyIDs(idx)) = idx;
        newNames{idx} = familyNames{familyIDs(idx)};
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
    figureHandles(end + 1) = plotParameterPairwiseFigure( ...
        relevanceStruct, opts);
end

function shapleyTable = buildParameterShapleySummaryTable(allBlocks, ...
    predictorNames)
    outcome = {'DeltaBias'; 'DeltaMask'};
    n = [allBlocks.bias.n; allBlocks.mask.n];
    totalModelR2 = [allBlocks.bias.OLS.R2; allBlocks.mask.OLS.R2];
    values = [allBlocks.bias.percentOfModelR2; ...
        allBlocks.mask.percentOfModelR2];
    rowSumPercent = sum(values, 2, 'omitnan');

    shapleyTable = table(outcome, n, totalModelR2, rowSumPercent, ...
        'VariableNames', {'Outcome', 'N', 'TotalModelR2', ...
        'RowSumPercent'});
    for predictorIdx = 1:numel(predictorNames)
        variableName = matlab.lang.makeValidName( ...
            [predictorNames{predictorIdx}, 'PercentModelR2']);
        shapleyTable.(variableName) = values(:, predictorIdx);
    end
end

function fig = plotParameterPairwiseFigure(relevanceStruct, opts)
    parameter = relevanceStruct.parameter;
    clusterLabels = relevanceStruct.powerEffectCluster;
    clusterIDs = unique(clusterLabels(:), 'stable');
    colors = lines(max(1, numel(clusterIDs)));
    predictorNames = parameter.predictorNames;
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

    fig = figure('Color', 'w', 'Name', 'Parameter-behavior relevance');
    set(fig, 'Position', [100, 100, 1050, 1200]);
    tlo = tiledlayout(numel(predictorNames), numel(outcomes), ...
        'TileSpacing', 'compact', 'Padding', 'compact');
    for predictorIdx = 1:numel(predictorNames)
        for outcomeIdx = 1:numel(outcomes)
            ax = nexttile;
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
            title(ax, sprintf('%s - %s', predictorLabels{predictorIdx}, ...
                outcomes{outcomeIdx}), 'FontName', 'Arial', ...
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
            box(ax, 'off');
            set(ax, 'FontName', 'Arial', 'FontSize', 12, ...
                'LineWidth', 1.0);
        end
    end

    title(tlo, sprintf('%s %s parameter-behavior scatter diagnostics', ...
        opts.monkeyName, opts.chamberWanted), 'FontName', 'Arial', ...
        'FontSize', 16, 'FontWeight', 'bold');
end

function fig = plotParameterShapleyFigure(relevanceStruct, opts)
    values = [
        relevanceStruct.parameter.allBlocks.bias.percentOfModelR2;
        relevanceStruct.parameter.allBlocks.mask.percentOfModelR2];
    rowSums = sum(values, 2, 'omitnan');
    fprintf(['Parameter Shapley row-sum sanity check | ' ...
        'DeltaBias=%.3f%% | DeltaMask=%.3f%%\n'], rowSums(1), rowSums(2));

    fig = figure('Color', 'w', 'Name', 'Parameter Shapley relevance');
    set(fig, 'Position', [100, 100, 900, 420]);
    ax = axes('Parent', fig);
    imagesc(ax, values);
    colormap(parula);
    cb = colorbar(ax);
    ylabel(cb, '% of model R^2', 'FontName', 'Arial', 'FontSize', 13);
    xticks(1:numel(relevanceStruct.parameter.predictorNames));
    xticklabels({'A', 'B', '\alpha', '\beta'});
    yticks(1:2);
    yticklabels({'\DeltaBias', '\DeltaMask'});
    set(ax, 'FontName', 'Arial', 'FontSize', 13, 'LineWidth', 1.0);
    title(ax, sprintf(['Parameter Shapley relevance | n=%d | ' ...
        'total R^2: bias=%.2f, mask=%.2f'], ...
        relevanceStruct.parameter.allBlocks.bias.n, ...
        relevanceStruct.parameter.allBlocks.bias.OLS.R2, ...
        relevanceStruct.parameter.allBlocks.mask.OLS.R2), ...
        'FontName', 'Arial', 'FontSize', 15, 'FontWeight', 'bold');
    addHeatmapText(values);
end

function fig = plotStimulationGroupedShapleyFigure(relevanceStruct, opts)
    familyNames = relevanceStruct.stimulation.allBlocks.bias.familyNames;
    if isempty(familyNames)
        fig = figure('Color', 'w', ...
            'Name', 'Stimulation-family Shapley relevance');
        axis off
        text(0.5, 0.5, 'No stimulation predictors available', ...
            'HorizontalAlignment', 'center', 'FontName', 'Arial');
        return
    end
    values = [
        relevanceStruct.stimulation.allBlocks.bias.groupedPercentOfModelR2;
        relevanceStruct.stimulation.allBlocks.mask.groupedPercentOfModelR2];
    fig = figure('Color', 'w', 'Name', 'Stimulation-family Shapley relevance');
    imagesc(values);
    colormap(parula);
    colorbar;
    xticks(1:numel(familyNames));
    xticklabels(familyNames);
    yticks(1:2);
    yticklabels({'Delta bias', 'Delta mask'});
    set(gca, 'FontName', 'Arial');
    title(sprintf(['Stimulation-family Shapley percent model R2 | n=%d | ' ...
        'R2 bias=%.2f mask=%.2f'], ...
        relevanceStruct.stimulation.allBlocks.bias.n, ...
        relevanceStruct.stimulation.allBlocks.bias.OLS.R2, ...
        relevanceStruct.stimulation.allBlocks.mask.OLS.R2), ...
        'FontName', 'Arial');
    addHeatmapText(values);
    sgtitle(sprintf('%s %s %s', opts.monkeyName, opts.chamberWanted, ...
        opts.aggregateField), 'FontName', 'Arial');
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
        'parameterPairwiseDiagnostic'};
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
