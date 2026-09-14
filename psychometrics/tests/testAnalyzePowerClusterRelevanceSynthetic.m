function testAnalyzePowerClusterRelevanceSynthetic
%TESTANALYZEPOWERCLUSTERRELEVANCESYNTHETIC Synthetic dominance checks.

    nBlocks = 36;
    rng(42);

    paramBL = repmat([3, 50, 15, 4], nBlocks, 1);
    dCon = randn(nBlocks, 4);
    dIncon = randn(nBlocks, 4);
    parameterContrast = dCon - dIncon;
    parameterContrast(:, 2) = linspace(-4, 4, nBlocks)' + ...
        0.1 .* randn(nBlocks, 1);
    dCon(:, 2) = parameterContrast(:, 2) ./ 2;
    dIncon(:, 2) = -parameterContrast(:, 2) ./ 2;

    clusterMdl = struct();
    clusterMdl.fittedParams = [paramBL, dCon, dIncon];
    clusterMdl.deltaBias = 5 + 8 .* parameterContrast(:, 2) + ...
        randn(nBlocks, 1);
    clusterMdl.deltaMask = randn(nBlocks, 1);
    clusterMdl.clusterBlocksIdx = (1:nBlocks)';

    energy = linspace(0.1, 7, nBlocks)';
    clusterMdl.deltaBias = clusterMdl.deltaBias + 0.5 .* energy;

    bitmapData = struct();
    bitmapData.meanPowerDensityWithinROI_mWmm2 = reshape( ...
        energy, 1, 1, []);
    bitmapData.totalPowerToOnPixelsWithinROI_mW = reshape( ...
        energy .* 1.2, 1, 1, []);
    bitmapData.projectorPowerDensity_mWmm2 = reshape( ...
        energy .* 0.8, 1, 1, []);
    bitmapData.nColumns = repmat(20, 1, nBlocks);
    bitmapData.pixelsON = repmat(100, 1, nBlocks);
    bitmapData.areaFinalROI = repmat(1, 1, nBlocks);
    bitmapData.areaPixelsONWithinROI = repmat(0.2, 1, nBlocks);
    bitmapData.spatialDutyCycleWithinROI = repmat(0.2, 1, nBlocks);
    bitmapData.temporalDutyCycle = repmat(0.5, 1, nBlocks);
    bitmapData.sensitivity = repmat(1, 1, nBlocks);
    bitmapData.adaptthresh = repmat(0.1, 1, nBlocks);

    clusterLabels = ones(nBlocks, 1);
    opts = struct('saveFlag', 0, 'makeFigures', false, ...
        'nBootstrap', 25, 'sourceMdlField', 'synthetic', ...
        'aggregateField', 'syntheticAggregate');

    relevanceStruct = analyzePowerClusterRelevance( ...
        clusterMdl, clusterLabels, bitmapData, opts);

    [~, dominantParam] = max( ...
        relevanceStruct.parameter.allBlocks.bias.percentOfModelR2);
    assert(dominantParam == 2, ...
        'Synthetic bias should be dominated by predictor B.');

    [~, dominantFamily] = max( ...
        relevanceStruct.stimulation.allBlocks.bias.groupedPercentOfModelR2);
    energyIdx = find(strcmp( ...
        relevanceStruct.stimulation.allBlocks.bias.familyNames, 'Energy'));
    assert(dominantFamily == energyIdx, ...
        'Synthetic stimulation relevance should be dominated by Energy.');
end
