function [out, figureHandles] = plotAggregatedPowerClusterPsychometrics(agg, opts)
% Plot and fit power-cluster aggregate psychometric data.
%
% [out, figureHandles] = plotAggregatedPowerClusterPsychometrics(agg, opts)
%
% agg is produced by aggregatePsychometricCountsByPowerCluster. Each power
% cluster receives one 2-by-3 figure. Fits use pooled binomial counts, not
% least-squares errors on percent-correct values.

    if nargin < 2 || isempty(opts)
        opts = struct();
    end
    opts = applyDefaults(opts);

    validateAggregateInput(agg);

    if isempty(agg)
        out = agg;
        figureHandles = gobjects(0);
        return;
    end

    agg = sortAggregateClusters(agg);

    style = getPlotStyle();
    viewNames = {'horizontal', 'vertical', 'merged'};
    viewTitles = {'Horizontal visual stimulus', ...
        'Vertical visual stimulus', 'Merged'};

    out = agg;
    figureHandles = gobjects(numel(agg), 1);
    for clusterIdx = 1:numel(agg)
        clusterTitle = buildClusterTitle(agg(clusterIdx));
        figureHandles(clusterIdx) = figure( ...
            'Name', clusterTitle, ...
            'Color', 'w', ...
            'Visible', opts.figureVisible);
        panelGap = [0.09 0.055];
        verticalMargins = [0.18 0.18];
        horizontalMargins = [0.15 0.17];
        makeSubplot = @(position) subtightplot(2, 3, position, ...
            panelGap, verticalMargins, horizontalMargins);
        topAxes = gobjects(1, 3);
        bottomAxes = gobjects(1, 3);
        deltaSummaries = repmat(struct('biasing', NaN, 'masking', NaN), 1, 3);
        mergedViewData = struct();
        mergedDeltaData = struct();
        signedBX0Fit = [];
        if strcmp(opts.modelTypeStr, 'weibullSignedBX0')
            signedBX0Fit = fitAggregateSignedBX0(agg(clusterIdx), opts);
            out(clusterIdx).signedBX0 = signedBX0Fit.signedBX0;
        end

        for viewIdx = 1:numel(viewNames)
            viewName = viewNames{viewIdx};
            viewData = agg(clusterIdx).(viewName);

            if strcmp(opts.modelTypeStr, 'weibullSignedBX0')
                fitResult = signedBX0Fit.signedBX0.panelFits.(viewName);
            else
                fitResult = fitAggregateView(viewData, opts);
            end
            conditionMeans = computeConditionMeans(viewData);
            deltaSummary = computeScalarDeltas(conditionMeans);
            deltaSummaries(viewIdx) = deltaSummary;
            deltaData = computeDeltaData(viewData);

            out(clusterIdx).(viewName).fit = fitResult;
            out(clusterIdx).(viewName).aggregateMean = conditionMeans;
            out(clusterIdx).(viewName).deltaBias = deltaSummary.biasing;
            out(clusterIdx).(viewName).deltaMask = deltaSummary.masking;
            out(clusterIdx).(viewName).delta = deltaData;

            axTop = makeSubplot(viewIdx);
            topAxes(viewIdx) = axTop;
            plotPsychometricPanel(axTop, viewData, fitResult, ...
                conditionMeans, viewTitles{viewIdx}, style, opts);

            axBottom = makeSubplot(viewIdx + 3);
            bottomAxes(viewIdx) = axBottom;
            plotDeltaPanel(axBottom, deltaData, fitResult, ...
                deltaSummary, style, opts);

            if strcmp(viewName, 'merged')
                mergedViewData = viewData;
                mergedDeltaData = deltaData;
            end
        end

        axes(bottomAxes(3));
        upFontSize(21, 0.01);
        styleAggregateAxes(topAxes, bottomAxes);

        if isfield(agg(clusterIdx), 'stimulationStats')
            addStimulationStatsText(topAxes(1), ...
                agg(clusterIdx).stimulationStats, agg(clusterIdx));
        end
        for viewIdx = 1:3
            addDeltaSummaryText(bottomAxes(viewIdx), deltaSummaries(viewIdx));
        end
        addFitParameterTable(topAxes(3), out(clusterIdx).merged.fit);
        addFigureTitle(clusterTitle);

        if opts.showDeltaPermutationStats
            axMergedDelta = bottomAxes(3);
            permSeed = stableAggregateDeltaPermutationSeed(...
                opts.deltaPermutationBaseSeed, agg(clusterIdx).clusterID);
            permResult = computeAggregateDeltaBiasPermutation(...
                mergedViewData, mergedDeltaData, opts.nDeltaPermutations, permSeed);
            if isempty(permResult.contrast)
                error('plotAggregatedPowerClusterPsychometrics:NoDeltaPermutationContrasts', ...
                    'No aggregate exact con/incon contrasts for C%d.', agg(clusterIdx).clusterID);
            end
            assertAggregatePermutationMatchesDisplayed(permResult, mergedDeltaData, ...
                agg(clusterIdx).clusterID);
            context = struct('type', 'cluster', ...
                'clusterID', agg(clusterIdx).clusterID);
            fprintf('Aggregate permutation ON | cluster C%d | contrasts %d\n', ...
                agg(clusterIdx).clusterID, numel(permResult.contrast));
            permAudit = addDeltaBiasPermutationVisualization(...
                axMergedDelta, permResult, context);
            assertAggregateDeltaPermutationVisualizationAudit(permAudit, ...
                permResult, agg(clusterIdx).clusterID);
            validateAggregateDeltaPermutationLegend(axMergedDelta, agg(clusterIdx).clusterID);
            saveAggregateDeltaPermutationExampleFigure(figureHandles(clusterIdx), opts);
            out(clusterIdx).merged.deltaBiasPermutation = permResult;
            out(clusterIdx).merged.deltaBiasPermutationClusterContrasts = ...
                permAudit.clusterContrasts;
            out(clusterIdx).merged.deltaBiasPermutationClusterSummary = ...
                permAudit.clusterSummary;
            printAggregateDeltaPermutationSummary(agg(clusterIdx).clusterID, ...
                numel(permResult.contrast), permResult.observedMeanDeltaBias, ...
                permResult.rawOverallTwoSidedP, permResult.overallSignificant, ...
                permResult.rawOverallPositiveOneSidedP);
        end
    end
end

function opts = applyDefaults(opts)
    if ~isstruct(opts) || ~isscalar(opts)
        error('opts must be a scalar struct.');
    end

    defaults = struct( ...
        'showSEM', true, ...
        'modelTypeStr', '', ...
        'figureVisible', 'on', ...
        'xLim', [], ...
        'yLim', [0 100], ...
        'deltaYLim', [-15 45], ...
        'markerFaceAlpha', 0.72, ...
        'deltaMarkerFaceAlpha', 0.62, ...
        'fitGridPoints', 401, ...
        'showDeltaPermutationStats', false, ...
        'nDeltaPermutations', 500, ...
        'deltaPermutationBaseSeed', 99173, ...
        'maxIterations', 5000, ...
        'tileSpacing', 'compact', ...
        'padding', 'compact');

    defaultNames = fieldnames(defaults);
    for ii = 1:numel(defaultNames)
        name = defaultNames{ii};
        if ~isfield(opts, name) || isempty(opts.(name))
            opts.(name) = defaults.(name);
        end
    end

    validateattributes(opts.showSEM, {'logical', 'numeric'}, {'scalar'});
    validateattributes(opts.fitGridPoints, {'numeric'}, ...
        {'scalar', 'integer', 'finite', '>=', 25});
    validateattributes(opts.maxIterations, {'numeric'}, ...
        {'scalar', 'integer', 'finite', 'positive'});
    if ~isempty(opts.xLim)
        validateattributes(opts.xLim, {'numeric'}, ...
            {'vector', 'numel', 2, 'real', 'finite', 'increasing'});
    end
    validateattributes(opts.yLim, {'numeric'}, ...
        {'vector', 'numel', 2, 'real', 'finite', 'increasing'});
    validateattributes(opts.deltaYLim, {'numeric'}, ...
        {'vector', 'numel', 2, 'real', 'finite', 'increasing'});
    validateattributes(opts.markerFaceAlpha, {'numeric'}, ...
        {'scalar', 'real', 'finite', '>=', 0, '<=', 1});
    validateattributes(opts.deltaMarkerFaceAlpha, {'numeric'}, ...
        {'scalar', 'real', 'finite', '>=', 0, '<=', 1});

    opts.showSEM = logical(opts.showSEM);
    opts.showDeltaPermutationStats = logical(opts.showDeltaPermutationStats);
    opts.figureVisible = validatestring(opts.figureVisible, {'on', 'off'});
end

function validateAggregateInput(agg)
    if ~isstruct(agg)
        error('agg must be a struct returned by aggregatePsychometricCountsByPowerCluster.');
    end

    requiredViews = {'horizontal', 'vertical', 'merged'};
    requiredConditions = {'baseline', 'con', 'incon'};
    requiredFields = {'x', 'successes', 'nTrials'};

    for clusterIdx = 1:numel(agg)
        if ~isfield(agg(clusterIdx), 'clusterID')
            error('agg(%d) is missing clusterID.', clusterIdx);
        end
        for viewIdx = 1:numel(requiredViews)
            viewName = requiredViews{viewIdx};
            if ~isfield(agg(clusterIdx), viewName)
                error('agg(%d) is missing view %s.', clusterIdx, viewName);
            end
            for conditionIdx = 1:numel(requiredConditions)
                conditionName = requiredConditions{conditionIdx};
                if ~isfield(agg(clusterIdx).(viewName), conditionName)
                    error('agg(%d).%s is missing condition %s.', ...
                        clusterIdx, viewName, conditionName);
                end
                conditionData = agg(clusterIdx).(viewName).(conditionName);
                missing = requiredFields(~isfield(conditionData, requiredFields));
                if ~isempty(missing)
                    error('agg(%d).%s.%s is missing field(s): %s.', ...
                        clusterIdx, viewName, conditionName, ...
                        strjoin(missing, ', '));
                end
            end
        end
    end
end

function style = getPlotStyle()
    baselineColor = [0 0 0];
    conColor = min([0.9294, 0.1098, 0.1373] .* 1.05, 1);
    inconColor = min([0, 0.0941, 0.6627] .* 1.25, 1);

    style.baseline.color = baselineColor;
    style.baseline.faceColor = [1 1 1];
    style.baseline.errorColor = baselineColor;
    style.baseline.marker = 'o';
    style.baseline.label = 'Baseline';
    style.con.color = conColor;
    style.con.faceColor = conColor;
    style.con.errorColor = conColor;
    style.con.marker = '^';
    style.con.label = 'Con-Opto';
    style.incon.color = inconColor;
    style.incon.faceColor = inconColor;
    style.incon.errorColor = inconColor;
    style.incon.marker = 'v';
    style.incon.label = 'Incon-Opto';
    style.bias.color = [127, 0, 255] ./ 255;
    style.bias.marker = 's';
    style.bias.label = 'Biasing';
    style.mask.color = [125, 125, 125] ./ 255;
    style.mask.marker = 's';
    style.mask.label = 'Masking';
    style.lineWidth = 3;
    style.topMarkerSize = 12;
    style.deltaMarkerSize = 15;
    style.referenceLineWidth = 1.5;
end

function signedFit = fitAggregateSignedBX0(clusterAgg, opts)
    signedChoice = clusterAgg.signedChoice;
    validateSignedAggregateSource(signedChoice);
    data = buildSignedBX0ObjectiveData(signedChoice);
    [initialParams, lb, ub] = getWeibullSignedBX0InitParams();
    [~, objectiveFunction] = getWeibullSignedBX0ModelFuncs(struct());

    starts = makeSignedBX0AggregateStarts(initialParams, lb, ub);
    [fitParams, nLL, exitFlag] = optimizeSignedBX0Starts(...
        objectiveFunction, data, starts, lb, ub, opts, false);
    [noX0FitParams, noX0NLL, noX0ExitFlag] = optimizeSignedBX0Starts(...
        objectiveFunction, data, starts, lb, ub, opts, true);

    nTrials = sum(data.sumBaselineChoice) + ...
        sum(data.sumHorizontalOptoChoice) + ...
        sum(data.sumVerticalOptoChoice);
    aicc = calculateAggregateAICc(nLL, 11, nTrials);
    noX0AICc = calculateAggregateAICc(noX0NLL, 10, nTrials);
    deltaAICcX0 = noX0AICc - aicc;
    weights = aggregateAkaikeWeights([noX0AICc, aicc]);

    if ~isSignedBX0AggregateFitValid(fitParams, objectiveFunction, data)
        error('plotAggregatedPowerClusterPsychometrics:InvalidSignedBX0AggregateFit', ...
            'Aggregate weibullSignedBX0 fit did not produce valid finite predictions.');
    end

    xGridMax = chooseSignedBX0GridMax(signedChoice, opts);
    xGrid = linspace(0, xGridMax, opts.fitGridPoints);
    displayCurves = projectSignedBX0DisplayCurves(xGrid, fitParams);

    signedFit = struct();
    signedFit.success = true;
    signedFit.message = '';
    signedFit.exitFlag = exitFlag;
    signedFit.noX0ExitFlag = noX0ExitFlag;
    signedFit.x = xGrid;
    signedFit.displayCurves = displayCurves;
    signedFit.data = data;
    signedFit.signedBX0 = struct( ...
        'modelVersion', 'fullBeta_slopeCap_v1', ...
        'parameterNames', {{'A_baseline', 'alpha_baseline', 'beta_baseline', ...
            'A_con', 'alpha_con', 'beta_con', ...
            'A_incon', 'alpha_incon', 'beta_incon', ...
            'deltaB', 'deltaX0'}}, ...
        'maxAllowedSlopePctPerContrast', 5.0, ...
        'slopeConstraintActive', true, ...
        'slopeDiagnostics', getWeibullSignedBX0SlopeDiagnostics(fitParams, 5.0), ...
        'noX0SlopeDiagnostics', getWeibullSignedBX0SlopeDiagnostics([noX0FitParams(1:10), 0], 5.0), ...
        'fitParams', fitParams, ...
        'nLL', nLL, ...
        'AICc', aicc, ...
        'noX0FitParams', noX0FitParams(1:10), ...
        'noX0NLL', noX0NLL, ...
        'noX0AICc', noX0AICc, ...
        'deltaAICcX0', deltaAICcX0, ...
        'akaikeWeightBOnly', weights(1), ...
        'akaikeWeightBX0', weights(2), ...
        'BHorizontal', 50 - fitParams(10), ...
        'BVertical', 50 + fitParams(10), ...
        'X0Horizontal', fitParams(11), ...
        'X0Vertical', -fitParams(11), ...
        'deltaB', fitParams(10), ...
        'deltaX0', fitParams(11), ...
        'fitStatus', 'ok', ...
        'nTrials', nTrials, ...
        'sourceAudit', signedChoice.audit, ...
        'kM0', 10, ...
        'kM1', 11);
    signedFit.signedBX0.panelFits = fitAggregateSignedBX0PanelFits(...
        clusterAgg, xGrid, fitParams(11), signedFit.signedBX0, opts);
end

function panelFits = fitAggregateSignedBX0PanelFits(clusterAgg, xGrid, globalDeltaX0, primarySignedBX0, opts)
    viewNames = {'horizontal', 'vertical', 'merged'};
    panelFits = struct();
    for viewIdx = 1:numel(viewNames)
        viewName = viewNames{viewIdx};
        fitResult = fitAggregateSignedBX0Panel(clusterAgg.(viewName), xGrid, ...
            globalDeltaX0, primarySignedBX0.fitParams, opts, viewName);
        fitResult.signedBX0.deltaAICcX0 = primarySignedBX0.deltaAICcX0;
        fitResult.signedBX0.akaikeWeightBX0 = primarySignedBX0.akaikeWeightBX0;
        fitResult.signedBX0.akaikeWeightBOnly = primarySignedBX0.akaikeWeightBOnly;
        fitResult.signedBX0.primaryFitParams = primarySignedBX0.fitParams;
        fitResult.signedBX0.primaryNLL = primarySignedBX0.nLL;
        fitResult.signedBX0.primaryAICc = primarySignedBX0.AICc;
        fitResult.signedBX0.primaryNoX0AICc = primarySignedBX0.noX0AICc;
        panelFits.(viewName) = fitResult;
    end
end

function fitResult = fitAggregateSignedBX0Panel(viewData, xGrid, globalDeltaX0, primaryParams, opts, viewName)
    fitResult = emptyFitResult();
    fitResult.modelType = 'weibullSignedBX0';
    fitResult.message = '';
    if ~hasAggregatePanelData(viewData)
        fitResult.message = 'All three aggregate conditions require populated binomial bins.';
        return;
    end

    [initialParams] = getWeibullSignedBX0InitParams();
    [lb, ub] = getWeibullSignedBX0PanelBounds();
    params0 = initialParams(1:10);
    if numel(primaryParams) >= 10 && all(isfinite(primaryParams(1:10)))
        params0 = primaryParams(1:10);
    end
    params0 = min(max(params0, lb), ub);
    starts = makeSignedBX0AggregatePanelStarts(params0, lb, ub);
    bestNLL = Inf;
    bestParams = nan(1, 10);
    bestExitFlag = NaN;
    for startIdx = 1:size(starts, 1)
        [candidateParams, exitFlag] = optimizeSignedBX0AggregatePanelStart(...
            starts(startIdx, :), lb, ub, viewData, globalDeltaX0, opts);
        candidateNLL = signedBX0AggregatePanelNLL(candidateParams, viewData, globalDeltaX0);
        diagnostics = getSignedBX0AggregatePanelSlopeDiagnostics(candidateParams, globalDeltaX0);
        if isfinite(candidateNLL) && candidateNLL < bestNLL && diagnostics.isValid
            bestNLL = candidateNLL;
            bestParams = candidateParams;
            bestExitFlag = exitFlag;
        end
    end
    if ~isfinite(bestNLL)
        error('plotAggregatedPowerClusterPsychometrics:SignedBX0PanelFitFailed', ...
            'All aggregate signed-BX0 panel starts failed for %s.', viewName);
    end

    validateSignedBX0AggregateBranchIndependence(bestParams, globalDeltaX0, viewName);
    [baselineY, conY, inconY] = predictSignedBX0AggregatePanelCurves(...
        xGrid, bestParams, globalDeltaX0);
    fitResult.success = true;
    fitResult.nLL = bestNLL;
    fitResult.exitFlag = bestExitFlag;
    fitResult.deltaParams = bestParams;
    fitResult.x = xGrid;
    fitResult.jointFitID = ['aggregateSignedBX0Panel_' viewName];
    fitResult.baseline.params = [bestParams(1), 50, bestParams(2), bestParams(3), 0];
    fitResult.con.params = [bestParams(4), 50 + bestParams(10), bestParams(5), bestParams(6), -globalDeltaX0];
    fitResult.incon.params = [bestParams(7), 50 - bestParams(10), bestParams(8), bestParams(9), +globalDeltaX0];
    fitResult.baseline.y = baselineY;
    fitResult.con.y = conY;
    fitResult.incon.y = inconY;
    fitResult.signedBX0 = struct( ...
        'fitParams', bestParams, ...
        'parameterNames', {{'A_baseline', 'alpha_baseline', 'beta_baseline', ...
        'A_con', 'alpha_con', 'beta_con', ...
        'A_incon', 'alpha_incon', 'beta_incon', 'deltaB_panel'}}, ...
        'globalDeltaX0', globalDeltaX0, ...
        'deltaB', bestParams(10), ...
        'sourcePanel', viewName, ...
        'nLL', bestNLL, ...
        'slopeDiagnostics', getSignedBX0AggregatePanelSlopeDiagnostics(bestParams, globalDeltaX0));
end

function ok = hasAggregatePanelData(viewData)
    ok = isfield(viewData, 'baseline') && isfield(viewData, 'con') && ...
        isfield(viewData, 'incon') && ~isempty(viewData.baseline.x) && ...
        ~isempty(viewData.con.x) && ~isempty(viewData.incon.x) && ...
        sum(viewData.baseline.nTrials) > 0 && sum(viewData.con.nTrials) > 0 && ...
        sum(viewData.incon.nTrials) > 0;
end

function starts = makeSignedBX0AggregatePanelStarts(params0, lb, ub)
    starts = repmat(params0(:)', 7, 1);
    starts(2, [3 6 9]) = 1.5;
    starts(3, [3 6 9]) = 5;
    starts(4, 10) = 10;
    starts(5, 10) = -10;
    starts(6, [3 6 9 10]) = [2.5 2.5 2.5 5];
    starts(7, [3 6 9 10]) = [6 6 6 -5];
    starts = min(max(starts, lb), ub);
    starts = unique(starts, 'rows', 'stable');
end

function [params, exitFlag] = optimizeSignedBX0AggregatePanelStart(params0, lb, ub, viewData, globalDeltaX0, opts)
    activeIdx = 1:10;
    u0 = paramsToUnit(params0, lb, ub);
    unitObjective = @(u) signedBX0AggregatePanelNLL(...
        lb + min(max(u(:)', 0), 1) .* (ub - lb), viewData, globalDeltaX0);
    if exist('fmincon', 'file') == 2
        fopts = optimoptions('fmincon', ...
            'Display', 'off', ...
            'MaxIterations', opts.maxIterations, ...
            'MaxFunctionEvaluations', 10 .* opts.maxIterations, ...
            'OptimalityTolerance', 1e-8, ...
            'StepTolerance', 1e-8);
        [uFit, ~, exitFlag] = fmincon(unitObjective, u0(activeIdx), [], [], [], [], ...
            zeros(size(u0(activeIdx))), ones(size(u0(activeIdx))), [], fopts);
    else
        fopts = optimset('Display', 'off', ...
            'MaxIter', opts.maxIterations, ...
            'MaxFunEvals', 10 .* opts.maxIterations, ...
            'TolX', 1e-8, ...
            'TolFun', 1e-8);
        [uFit, ~, exitFlag] = fminsearchbnd(unitObjective, u0(activeIdx), ...
            zeros(size(u0(activeIdx))), ones(size(u0(activeIdx))), fopts);
    end
    params = lb + min(max(uFit(:)', 0), 1) .* (ub - lb);
end

function nLL = signedBX0AggregatePanelNLL(params, viewData, globalDeltaX0)
    if numel(params) ~= 10 || any(~isfinite(params)) || ~isfinite(globalDeltaX0)
        nLL = 1e12;
        return;
    end
    diagnostics = getSignedBX0AggregatePanelSlopeDiagnostics(params, globalDeltaX0);
    if ~diagnostics.isValid
        excess = diagnostics.slopeExcess;
        excess(~isfinite(excess)) = 5.0;
        nLL = 1e12 + 1e6 .* sum(excess .^ 2);
        return;
    end
    [baselineY, conY, inconY] = predictSignedBX0AggregatePanelCurves(...
        [], params, globalDeltaX0, viewData);
    nLL = aggregatePanelConditionNLL(viewData.baseline, baselineY) + ...
        aggregatePanelConditionNLL(viewData.con, conY) + ...
        aggregatePanelConditionNLL(viewData.incon, inconY);
end

function nLL = aggregatePanelConditionNLL(conditionData, prediction)
    probability = min(max(prediction ./ 100, 1e-10), 1 - 1e-10);
    successes = conditionData.successes;
    failures = conditionData.nTrials - successes;
    nLL = -sum(successes .* log(probability) + failures .* log(1 - probability));
end

function [baselineY, conY, inconY] = predictSignedBX0AggregatePanelCurves(xGrid, params, globalDeltaX0, viewData)
    if nargin >= 4 && ~isempty(viewData)
        xBaseline = viewData.baseline.x;
        xCon = viewData.con.x;
        xIncon = viewData.incon.x;
    else
        xBaseline = xGrid;
        xCon = xGrid;
        xIncon = xGrid;
    end
    baselineY = predictShiftedWeibullBranch(abs(xBaseline), params(1), 50, params(2), params(3), 0);
    conY = predictShiftedWeibullBranch(abs(xCon), params(4), 50 + params(10), params(5), params(6), -globalDeltaX0);
    inconY = predictShiftedWeibullBranch(abs(xIncon), params(7), 50 - params(10), params(8), params(9), +globalDeltaX0);
end

function diagnostics = getSignedBX0AggregatePanelSlopeDiagnostics(params, globalDeltaX0)
    maxSlopePctPerContrast = 5.0;
    B_con = 50 + params(10);
    B_incon = 50 - params(10);
    amplitudes = [(100 - params(1)) - 50, ...
        (100 - params(4)) - B_con, ...
        (100 - params(7)) - B_incon];
    alphas = [params(2), params(5), params(8)];
    betas = [params(3), params(6), params(9)];
    slopeValues = nan(1, 3);
    if all(isfinite([amplitudes, alphas, betas, globalDeltaX0])) && ...
            all(amplitudes > 0) && all(alphas > 0) && all(betas > 1)
        slopeValues = getWeibullHalfMaxSlope(amplitudes, alphas, betas);
    end
    diagnostics.maxSlopeValues = slopeValues;
    diagnostics.maxSlopeOverall = max(slopeValues, [], 'omitnan');
    diagnostics.maxAllowedSlopePctPerContrast = maxSlopePctPerContrast;
    diagnostics.slopeConstraintActive = true;
    diagnostics.slopeCapActive = isfinite(diagnostics.maxSlopeOverall) && ...
        diagnostics.maxSlopeOverall > maxSlopePctPerContrast;
    diagnostics.isValid = all(isfinite(slopeValues)) && all(slopeValues <= maxSlopePctPerContrast);
    diagnostics.slopeExcess = max(0, slopeValues - maxSlopePctPerContrast);
end

function validateSignedBX0AggregateBranchIndependence(params, globalDeltaX0, viewName)
% Assert cross-branch independence of the 10-parameter panel prediction.
%
% Parameter groups:
%   baseline 1:3  — must not affect con (curve 2) or incon (curve 3)
%   con      4:6  — must not affect incon (curve 3)
%   incon    7:9  — must not affect con  (curve 2)
%   deltaB   10   — shared by design; excluded from independence checks

    if nargin < 3; viewName = 'unknown'; end

    [lb, ub] = getWeibullSignedBX0PanelBounds();
    c        = linspace(0, 1, 51);
    [base0, con0, incon0] = predictSignedBX0AggregatePanelCurves(c, params, globalDeltaX0);

    step      = max(0.05 .* (ub - lb), 1e-6);   % 5 % of legal range, floored
    indTol    = 1e-9;   % cross-branch contamination threshold
    usefulTol = 1e-6;   % own-branch change required for perturbation to be informative

    % {paramIndices, ownCurveIdx, crossCurveIndices}  (curve 1=base,2=con,3=incon)
    groups = { 1:3, 1, [2 3]; ...
               4:6, 2, 3;     ...
               7:9, 3, 2      };

    violations = {};
    for gi = 1:size(groups, 1)
        paramIdxs = groups{gi, 1};
        ownIdx    = groups{gi, 2};
        crossIdxs = groups{gi, 3};
        for pi = paramIdxs
            pPert = params;
            if params(pi) + step(pi) <= ub(pi)
                pPert(pi) = params(pi) + step(pi);
            else
                pPert(pi) = params(pi) - step(pi);
            end
            pPert(pi) = min(max(pPert(pi), lb(pi)), ub(pi));
            [bP, cP, iP] = predictSignedBX0AggregatePanelCurves(c, pPert, globalDeltaX0);
            dAll = [max(abs(bP(:) - base0(:)), [], 'omitnan'), ...
                    max(abs(cP(:) - con0(:)),   [], 'omitnan'), ...
                    max(abs(iP(:) - incon0(:)), [], 'omitnan')];
            if dAll(ownIdx) < usefulTol; continue; end  % uninformative; skip
            for ci = crossIdxs(:)'
                if dAll(ci) > indTol
                    violations{end+1} = sprintf( ...  %#ok<AGROW>
                        'p%d: contaminates curve %d (d=%.2e, own=%.2e)', ...
                        pi, ci, dAll(ci), dAll(ownIdx));
                end
            end
        end
    end

    if ~isempty(violations)
        fprintf('validateSignedBX0AggregateBranchIndependence FAILED (panel=%s)\n', viewName);
        fprintf('  params=[%s]  deltaX0=%.4g  indTol=%.0e  usefulTol=%.0e\n', ...
            num2str(params, '%.4g '), globalDeltaX0, indTol, usefulTol);
        for vi = 1:numel(violations); fprintf('  %s\n', violations{vi}); end
        error('plotAggregatedPowerClusterPsychometrics:SignedBX0BranchSwitching', ...
            'Aggregate signed-BX0 displayed branches switch con/incon shape parameters.');
    end
end
function validateSignedAggregateSource(signedChoice)
    if isempty(signedChoice) || ~isstruct(signedChoice)
        error('plotAggregatedPowerClusterPsychometrics:MissingSignedBX0AggregateSource', ...
            'weibullSignedBX0 aggregate fitting requires agg(cluster).signedChoice.');
    end
    requiredConditions = {'baseline', 'horizontalOpto', 'verticalOpto'};
    requiredFields = {'x', 'successes', 'nTrials'};
    for conditionIdx = 1:numel(requiredConditions)
        conditionName = requiredConditions{conditionIdx};
        if ~isfield(signedChoice, conditionName)
            error('plotAggregatedPowerClusterPsychometrics:MissingSignedBX0Condition', ...
                'signedChoice is missing physical condition %s.', conditionName);
        end
        conditionData = signedChoice.(conditionName);
        missing = requiredFields(~isfield(conditionData, requiredFields));
        if ~isempty(missing)
            error('plotAggregatedPowerClusterPsychometrics:MissingSignedBX0Fields', ...
                'signedChoice.%s is missing field(s): %s.', ...
                conditionName, strjoin(missing, ', '));
        end
        x = conditionData.x;
        successes = conditionData.successes;
        nTrials = conditionData.nTrials;
        valid = isfinite(x) & isfinite(successes) & isfinite(nTrials) & ...
            nTrials > 0 & successes >= 0 & successes <= nTrials;
        if isempty(x) || ~all(valid)
            error('plotAggregatedPowerClusterPsychometrics:InvalidSignedBX0Counts', ...
                'signedChoice.%s contains invalid or empty signed binomial counts.', ...
                conditionName);
        end
        if ~any(x < 0) || ~any(x > 0)
            warning('plotAggregatedPowerClusterPsychometrics:SignedBX0ContrastCoverage', ...
                'signedChoice.%s does not contain both negative and positive signed contrasts.', ...
                conditionName);
        end
    end
end

function data = buildSignedBX0ObjectiveData(signedChoice)
    data = struct( ...
        'xBaselineChoice', signedChoice.baseline.x, ...
        'sumBaselineChoice', signedChoice.baseline.nTrials, ...
        'successBaselineChoice', signedChoice.baseline.successes, ...
        'xHorizontalOptoChoice', signedChoice.horizontalOpto.x, ...
        'sumHorizontalOptoChoice', signedChoice.horizontalOpto.nTrials, ...
        'successHorizontalOptoChoice', signedChoice.horizontalOpto.successes, ...
        'xVerticalOptoChoice', signedChoice.verticalOpto.x, ...
        'sumVerticalOptoChoice', signedChoice.verticalOpto.nTrials, ...
        'successVerticalOptoChoice', signedChoice.verticalOpto.successes);
end

function starts = makeSignedBX0AggregateStarts(initialParams, lb, ub)
    starts = [ ...
        initialParams; ...
        setSignedBX0Deltas(initialParams, 0, 0); ...
        setSignedBX0Deltas(initialParams, 5, 0); ...
        setSignedBX0Deltas(initialParams, -5, 0); ...
        setSignedBX0Deltas(initialParams, 0, 3); ...
        setSignedBX0Deltas(initialParams, 0, -3); ...
        setSignedBX0Deltas(initialParams, 5, 3); ...
        setSignedBX0Deltas(initialParams, 5, -3); ...
        setSignedBX0Deltas(initialParams, -5, 3); ...
        setSignedBX0Deltas(initialParams, -5, -3)];
    starts = min(max(starts, lb), ub);
    starts = unique(starts, 'rows', 'stable');
end

function params = setSignedBX0Deltas(params, deltaB, deltaX0)
    params(10) = deltaB;
    params(11) = deltaX0;
end

function [bestParams, bestNLL, bestExitFlag] = optimizeSignedBX0Starts(...
        objectiveFunction, data, starts, lb, ub, opts, forceNoX0)
    bestParams = nan(1, 11);
    bestNLL = Inf;
    bestExitFlag = NaN;
    for startIdx = 1:size(starts, 1)
        start = starts(startIdx, :);
        if forceNoX0
            start(11) = 0;
        end
        [candidateParams, exitFlag] = optimizeSignedBX0SingleStart(...
            objectiveFunction, data, start, lb, ub, opts, forceNoX0);
        candidateNLL = objectiveFunction(candidateParams, data);
        if isSignedBX0AggregateFitValid(candidateParams, objectiveFunction, data) && ...
                candidateNLL < bestNLL
            bestParams = candidateParams;
            bestNLL = candidateNLL;
            bestExitFlag = exitFlag;
        end
    end
    if ~isfinite(bestNLL)
        error('plotAggregatedPowerClusterPsychometrics:SignedBX0OptimizationFailed', ...
            'All aggregate weibullSignedBX0 deterministic starts failed.');
    end
end

function [params, exitFlag] = optimizeSignedBX0SingleStart(...
        objectiveFunction, data, start, lb, ub, opts, forceNoX0)
    if forceNoX0
        activeIdx = 1:10;
    else
        activeIdx = 1:11;
    end
    activeLb = lb(activeIdx);
    activeUb = ub(activeIdx);
    u0 = paramsToUnit(start(activeIdx), activeLb, activeUb);
    unitObjective = @(u) objectiveFunction(unitToSignedBX0Params(...
        u, activeIdx, activeLb, activeUb, forceNoX0), data);

    if exist('fmincon', 'file') == 2
        fopts = optimoptions('fmincon', ...
            'Display', 'off', ...
            'MaxIterations', opts.maxIterations, ...
            'MaxFunctionEvaluations', 10 .* opts.maxIterations, ...
            'OptimalityTolerance', 1e-8, ...
            'StepTolerance', 1e-8);
        [uFit, ~, exitFlag] = fmincon(unitObjective, u0, [], [], [], [], ...
            zeros(size(u0)), ones(size(u0)), [], fopts);
    else
        fopts = optimset('Display', 'off', ...
            'MaxIter', opts.maxIterations, ...
            'MaxFunEvals', 10 .* opts.maxIterations, ...
            'TolX', 1e-8, ...
            'TolFun', 1e-8);
        [uFit, ~, exitFlag] = fminsearchbnd(unitObjective, u0, ...
            zeros(size(u0)), ones(size(u0)), fopts);
    end

    params = unitToSignedBX0Params(uFit, activeIdx, activeLb, activeUb, forceNoX0);
end

function u = paramsToUnit(params, lb, ub)
    u = (params - lb) ./ (ub - lb);
    u = min(max(u, 0), 1);
end

function params = unitToSignedBX0Params(u, activeIdx, activeLb, activeUb, forceNoX0)
    params = nan(1, 11);
    params(activeIdx) = activeLb + u .* (activeUb - activeLb);
    if forceNoX0
        params(11) = 0;
    end
end

function ok = isSignedBX0AggregateFitValid(params, objectiveFunction, data)
    ok = numel(params) == 11 && all(isfinite(params));
    if ~ok
        return;
    end
    nLL = objectiveFunction(params, data);
    if ~isfinite(nLL) || nLL >= 1e11
        ok = false;
        return;
    end
    xTest = unique([data.xBaselineChoice, data.xHorizontalOptoChoice, data.xVerticalOptoChoice]);
    if isempty(xTest)
        ok = false;
        return;
    end
    curves = evaluateSignedBX0PhysicalPredictions(xTest, params);
    values = [curves.baseline, curves.horizontalOpto, curves.verticalOpto];
    ok = isreal(values) && all(isfinite(values)) && ...
        all(values >= 0) && all(values <= 100);
end

function curves = evaluateSignedBX0PhysicalPredictions(x, params)
    curves.baseline = weibullSignedBX0Mdl(x, ...
        params(1), params(2), params(3), ...
        params(1), params(2), params(3), 50, 0);
    curves.horizontalOpto = weibullSignedBX0Mdl(x, ...
        params(4), params(5), params(6), ...
        params(7), params(8), params(9), 50 - params(10), params(11));
    curves.verticalOpto = weibullSignedBX0Mdl(x, ...
        params(7), params(8), params(9), ...
        params(4), params(5), params(6), 50 + params(10), -params(11));
end

function aicc = calculateAggregateAICc(nLL, k, n)
    aic = 2 .* k + 2 .* nLL;
    if n > (k + 1)
        aicc = aic + (2 .* k .* (k + 1)) ./ (n - k - 1);
    else
        aicc = Inf;
    end
end

function weights = aggregateAkaikeWeights(aiccValues)
    finiteAICc = aiccValues(isfinite(aiccValues));
    if isempty(finiteAICc)
        weights = nan(size(aiccValues));
        return;
    end
    delta = aiccValues - min(finiteAICc);
    relLike = exp(-0.5 .* delta);
    weights = relLike ./ sum(relLike(isfinite(relLike)));
end

function xGridMax = chooseSignedBX0GridMax(signedChoice, opts)
    if ~isempty(opts.xLim)
        xGridMax = opts.xLim(2);
        return;
    end
    allX = abs([signedChoice.baseline.x, ...
        signedChoice.horizontalOpto.x, signedChoice.verticalOpto.x]);
    if isempty(allX)
        xGridMax = 100;
    else
        xGridMax = max(100, ceil(max(allX) ./ 5) .* 5);
    end
end

function displayCurves = projectSignedBX0DisplayCurves(c, params)
    pBLneg = weibullSignedBX0Mdl(-c, params(1), params(2), params(3), ...
        params(1), params(2), params(3), 50, 0);
    pBLpos = weibullSignedBX0Mdl(+c, params(1), params(2), params(3), ...
        params(1), params(2), params(3), 50, 0);
    pHneg = weibullSignedBX0Mdl(-c, params(4), params(5), params(6), ...
        params(7), params(8), params(9), 50 - params(10), params(11));
    pHpos = weibullSignedBX0Mdl(+c, params(4), params(5), params(6), ...
        params(7), params(8), params(9), 50 - params(10), params(11));
    pVneg = weibullSignedBX0Mdl(-c, params(7), params(8), params(9), ...
        params(4), params(5), params(6), 50 + params(10), -params(11));
    pVpos = weibullSignedBX0Mdl(+c, params(7), params(8), params(9), ...
        params(4), params(5), params(6), 50 + params(10), -params(11));

    displayCurves.horizontal.baseline = 100 - pBLneg;
    displayCurves.horizontal.con = 100 - pHneg;
    displayCurves.horizontal.incon = 100 - pVneg;
    displayCurves.vertical.baseline = pBLpos;
    displayCurves.vertical.con = pVpos;
    displayCurves.vertical.incon = pHpos;
    displayCurves.merged.baseline = 0.5 .* (...
        displayCurves.horizontal.baseline + displayCurves.vertical.baseline);
    displayCurves.merged.con = 0.5 .* (...
        displayCurves.horizontal.con + displayCurves.vertical.con);
    displayCurves.merged.incon = 0.5 .* (...
        displayCurves.horizontal.incon + displayCurves.vertical.incon);
end

function ok = areSignedBX0DisplayCurvesValid(displayCurves)
    viewNames = {'horizontal', 'vertical', 'merged'};
    conditionNames = {'baseline', 'con', 'incon'};
    ok = true;
    for viewIdx = 1:numel(viewNames)
        viewName = viewNames{viewIdx};
        for conditionIdx = 1:numel(conditionNames)
            conditionName = conditionNames{conditionIdx};
            values = displayCurves.(viewName).(conditionName);
            ok = ok && isreal(values) && all(isfinite(values)) && ...
                all(values >= 0) && all(values <= 100);
        end
    end
end
function fitResult = projectSignedBX0FitToView(signedFit, viewName)
    fitResult = emptyFitResult();
    fitResult.modelType = 'weibullSignedBX0';
    fitResult.success = signedFit.success;
    fitResult.message = signedFit.message;
    fitResult.nLL = signedFit.signedBX0.nLL;
    fitResult.exitFlag = signedFit.exitFlag;
    fitResult.deltaParams = signedFit.signedBX0.fitParams;
    fitResult.x = signedFit.x;
    fitResult.signedBX0 = signedFit.signedBX0;
    fitResult.jointFitID = 'aggregateSignedBX0';
    params = signedFit.signedBX0.fitParams;
    fitResult.baseline.params = [params(1), 50, params(2), params(3), 0];
    switch viewName
        case 'horizontal'
            fitResult.con.params = [params(4), 50 - params(10), params(5), params(6), params(11)];
            fitResult.incon.params = [params(7), 50 - params(10), params(8), params(9), params(11)];
        case 'vertical'
            fitResult.con.params = [params(4), 50 + params(10), params(5), params(6), -params(11)];
            fitResult.incon.params = [params(7), 50 + params(10), params(8), params(9), -params(11)];
        otherwise
            fitResult.con.params = [NaN NaN NaN NaN NaN];
            fitResult.incon.params = [NaN NaN NaN NaN NaN];
    end
    fitResult.baseline.y = signedFit.displayCurves.(viewName).baseline;
    fitResult.con.y = signedFit.displayCurves.(viewName).con;
    fitResult.incon.y = signedFit.displayCurves.(viewName).incon;
end
function fitResult = fitAggregateView(viewData, opts)
    conditionNames = {'baseline', 'con', 'incon'};
    fitResult = emptyFitResult();

    haveData = false(1, numel(conditionNames));
    for conditionIdx = 1:numel(conditionNames)
        conditionData = viewData.(conditionNames{conditionIdx});
        haveData(conditionIdx) = numel(conditionData.x) >= 2 && ...
            sum(conditionData.nTrials) > 0;
    end

    if ~all(haveData)
        fitResult.message = 'All three conditions require at least two populated bins.';
        return;
    end

    q0 = [0.08, 15, 3, 0, 0, 0, 0, 0, 0, 0];
    lower = [0, 1e-3, 0.05, -0.20, -0.30, -14.9, -2.95, -0.20, -14.9, -2.95];
    upper = [0.49, 200, 12, 0.20, 0.30, 185, 9, 0.20, 185, 9];

    allX = [viewData.baseline.x, viewData.con.x, viewData.incon.x];
    positiveX = allX(allX > 0);
    if ~isempty(positiveX)
        q0(2) = median(positiveX);
        upper(2) = max(upper(2), 4 .* max(positiveX));
        upper([6 9]) = upper(2);
    end

    objective = @(q) jointBinomialNLL(q, viewData);
    penaltyObjective = @(q) objective(q) + parameterPenalty(q, lower, upper);

    optimizerOptions = optimset( ...
        'Display', 'off', ...
        'MaxIter', opts.maxIterations, ...
        'MaxFunEvals', 10 .* opts.maxIterations, ...
        'TolX', 1e-8, ...
        'TolFun', 1e-8);

    starts = [ ...
        q0; ...
        q0 + [0.04 0 1 0 0 0 0 0 0 0]; ...
        q0 + [0 5 -1 0 0.05 0 0 0 0 0]; ...
        q0 + [0 -5 2 0 -0.05 0 0 0 0 0]];

    bestQ = q0;
    bestNLL = Inf;
    bestExitFlag = NaN;

    for startIdx = 1:size(starts, 1)
        start = min(max(starts(startIdx, :), lower), upper);
        [candidateQ, candidateObjective, exitFlag] = fminsearch( ...
            penaltyObjective, start, optimizerOptions);
        candidateQ = min(max(candidateQ, lower), upper);
        candidateNLL = objective(candidateQ);

        if isfinite(candidateObjective) && candidateNLL < bestNLL
            bestQ = candidateQ;
            bestNLL = candidateNLL;
            bestExitFlag = exitFlag;
        end
    end

    params = unpackParameters(bestQ);
    if ~areParametersValid(params)
        fitResult.message = 'Optimization did not produce valid Weibull parameters.';
        return;
    end

    xMax = max(allX);
    if isempty(opts.xLim)
        xGridMax = max(100, ceil(xMax ./ 5) .* 5);
    else
        xGridMax = opts.xLim(2);
    end
    xGrid = linspace(0, xGridMax, opts.fitGridPoints);

    fitResult.success = true;
    fitResult.message = '';
    fitResult.nLL = bestNLL;
    fitResult.exitFlag = bestExitFlag;
    fitResult.deltaParams = bestQ;
    fitResult.x = xGrid;
    fitResult.baseline.params = params.baseline;
    fitResult.baseline.y = weibullCurve(xGrid, params.baseline);
    fitResult.con.params = params.con;
    fitResult.con.y = weibullCurve(xGrid, params.con);
    fitResult.incon.params = params.incon;
    fitResult.incon.y = weibullCurve(xGrid, params.incon);
end

function fitResult = emptyFitResult()
    emptyCondition = struct('params', nan(1, 4), 'y', []);
    fitResult = struct( ...
        'modelType', '', ...
        'success', false, ...
        'message', '', ...
        'nLL', NaN, ...
        'exitFlag', NaN, ...
        'deltaParams', nan(1, 10), ...
        'x', [], ...
        'baseline', emptyCondition, ...
        'con', emptyCondition, ...
        'incon', emptyCondition);
end

function nLL = jointBinomialNLL(q, viewData)
    params = unpackParameters(q);
    if ~areParametersValid(params)
        nLL = 1e12;
        return;
    end

    conditionNames = {'baseline', 'con', 'incon'};
    nLL = 0;
    for conditionIdx = 1:numel(conditionNames)
        conditionName = conditionNames{conditionIdx};
        conditionData = viewData.(conditionName);
        probability = weibullCurve(conditionData.x, params.(conditionName)) ./ 100;
        probability = min(max(probability, 1e-10), 1 - 1e-10);
        successes = conditionData.successes;
        failures = conditionData.nTrials - successes;
        nLL = nLL - sum(successes .* log(probability) + ...
            failures .* log(1 - probability));
    end
end

function penalty = parameterPenalty(q, lower, upper)
    below = max(lower - q, 0);
    above = max(q - upper, 0);
    penalty = 1e8 .* sum(below .^ 2 + above .^ 2);

    params = unpackParameters(q);
    conditionNames = {'baseline', 'con', 'incon'};
    for conditionIdx = 1:numel(conditionNames)
        p = params.(conditionNames{conditionIdx});
        invalid = [ ...
            max(-p(1), 0), ...
            max(p(1) - 0.49, 0), ...
            max(0.001 - p(2), 0), ...
            max(p(2) - 0.999, 0), ...
            max(1e-3 - p(3), 0), ...
            max(0.05 - p(4), 0), ...
            max(p(2) - (1 - p(1) - 1e-4), 0)];
        penalty = penalty + 1e10 .* sum(invalid .^ 2);
    end
end

function params = unpackParameters(q)
    baseline = [q(1), 0.5, q(2), q(3)];
    con = [q(1) + q(4), 0.5 + q(5), q(2) + q(6), q(3) + q(7)];
    incon = [q(1) + q(8), 0.5 - q(5), q(2) + q(9), q(3) + q(10)];
    params = struct('baseline', baseline, 'con', con, 'incon', incon);
end

function valid = areParametersValid(params)
    conditionNames = {'baseline', 'con', 'incon'};
    valid = true;
    for conditionIdx = 1:numel(conditionNames)
        p = params.(conditionNames{conditionIdx});
        valid = valid && all(isfinite(p)) && ...
            p(1) >= 0 && p(1) <= 0.49 && ...
            p(2) > 0 && p(2) < 1 && ...
            p(3) > 0 && p(4) > 0 && ...
            (1 - p(1)) > p(2);
    end
end

function y = weibullCurve(x, params)
    A = params(1);
    B = params(2);
    alpha = params(3);
    beta = params(4);
    y = 100 .* (B + (1 - exp(-((x ./ alpha) .^ beta))) .* ...
        ((1 - A) - B));
end

function means = computeConditionMeans(viewData)
    means = struct();
    conditionNames = {'baseline', 'con', 'incon'};
    for conditionIdx = 1:numel(conditionNames)
        conditionName = conditionNames{conditionIdx};
        conditionData = viewData.(conditionName);
        totalTrials = sum(conditionData.nTrials);
        if totalTrials > 0
            means.(conditionName) = 100 .* ...
                sum(conditionData.successes) ./ totalTrials;
        else
            means.(conditionName) = NaN;
        end
    end
end

function summary = computeScalarDeltas(means)
    summary.biasing = means.con - means.incon;
    summary.masking = means.baseline - mean([means.con, means.incon], 'omitnan');
end

function deltaData = computeDeltaData(viewData)
    [biasX, idxCon, idxIncon] = intersect( ...
        viewData.con.x, viewData.incon.x, 'stable');
    bias = viewData.con.pctCorrect(idxCon) - ...
        viewData.incon.pctCorrect(idxIncon);
    biasSEM = combineSEM( ...
        getSEM(viewData.con, idxCon), ...
        getSEM(viewData.incon, idxIncon), 1, 1);

    optoMean = 0.5 .* (viewData.con.pctCorrect(idxCon) + ...
        viewData.incon.pctCorrect(idxIncon));
    optoSEM = combineSEM( ...
        getSEM(viewData.con, idxCon), ...
        getSEM(viewData.incon, idxIncon), 0.5, 0.5);

    [maskX, idxBaseline, idxOpto] = intersect( ...
        viewData.baseline.x, biasX, 'stable');
    masking = viewData.baseline.pctCorrect(idxBaseline) - optoMean(idxOpto);
    maskingSEM = combineSEM( ...
        getSEM(viewData.baseline, idxBaseline), optoSEM(idxOpto), 1, 1);

    deltaData = struct( ...
        'biasX', biasX, ...
        'biasing', bias, ...
        'biasSEM', biasSEM, ...
        'maskX', maskX, ...
        'masking', masking, ...
        'maskSEM', maskingSEM);
end

function sem = getSEM(conditionData, indices)
    if isempty(indices)
        sem = [];
    elseif isfield(conditionData, 'sessionSEM') && ...
            numel(conditionData.sessionSEM) >= max(indices)
        sem = conditionData.sessionSEM(indices);
    else
        sem = nan(size(indices));
    end
end

function combined = combineSEM(firstSEM, secondSEM, firstScale, secondScale)
    combined = sqrt((firstScale .* firstSEM) .^ 2 + ...
        (secondScale .* secondSEM) .^ 2);
end

function plotPsychometricPanel(ax, viewData, fitResult, conditionMeans, ...
        titleText, style, opts)
    hold(ax, 'on');
    yline(ax, 50, '--', 'Color', 0.4 .* [1 1 1], ...
        'LineWidth', style.referenceLineWidth, 'HandleVisibility', 'off');

    conditionNames = {'baseline', 'con', 'incon'};
    handles = gobjects(1, numel(conditionNames));
    for conditionIdx = 1:numel(conditionNames)
        conditionName = conditionNames{conditionIdx};
        conditionStyle = style.(conditionName);

        if fitResult.success
            plot(ax, fitResult.x, fitResult.(conditionName).y, ...
                '-', 'Color', conditionStyle.color, ...
                'LineWidth', style.lineWidth, ...
                'HandleVisibility', 'off');
        end
    end

    for conditionIdx = 1:numel(conditionNames)
        conditionName = conditionNames{conditionIdx};
        conditionData = viewData.(conditionName);
        conditionStyle = style.(conditionName);
        handles(conditionIdx) = plotAggregatePoints(ax, conditionData, ...
            conditionStyle, opts, style);
    end

    xLimits = chooseXLimits(viewData, opts);
    xlim(ax, xLimits);
    ylim(ax, opts.yLim);
    plotRightEdgeMeanTicks(ax, ...
        [conditionMeans.baseline, conditionMeans.con, conditionMeans.incon], ...
        {style.baseline.color, style.con.color, style.incon.color}, ...
        style.lineWidth);
    xlabel(ax, 'Gabor contrast (%)');
    ylabel(ax, 'Correct (%)');
    title(ax, titleText, 'FontWeight', 'normal');
    legend(ax, handles, {'Baseline', 'Con-Opto', 'Incon-Opto'}, ...
        'Location', 'southeast');
    box(ax, 'off');
    axis(ax, 'square');
end

function handle = plotAggregatePoints(ax, conditionData, conditionStyle, opts, style)
    if isempty(conditionData.x)
        handle = plot(ax, NaN, NaN, conditionStyle.marker, ...
            'Color', conditionStyle.color, ...
            'MarkerFaceColor', conditionStyle.faceColor, ...
            'MarkerSize', style.topMarkerSize, ...
            'LineWidth', style.lineWidth, ...
            'DisplayName', conditionStyle.label);
        return;
    end

    if opts.showSEM && isfield(conditionData, 'sessionSEM')
        sem = conditionData.sessionSEM;
        errorbar(ax, conditionData.x, conditionData.y, sem, ...
            'LineStyle', 'none', ...
            'Color', conditionStyle.errorColor, ...
            'LineWidth', style.lineWidth, ...
            'CapSize', 0, ...
            'HandleVisibility', 'off');
    end

    handle = scatter(ax, conditionData.x, conditionData.y, ...
        style.topMarkerSize .^ 2, ...
        conditionStyle.faceColor, ...
        conditionStyle.marker, ...
        'filled', ...
        'MarkerFaceColor', conditionStyle.faceColor, ...
        'MarkerEdgeColor', 'k', ...
        'MarkerFaceAlpha', opts.markerFaceAlpha, ...
        'MarkerEdgeAlpha', 1, ...
        'LineWidth', style.lineWidth, ...
        'DisplayName', conditionStyle.label);
end

function plotDeltaPanel(ax, deltaData, fitResult, deltaSummary, style, opts)
    hold(ax, 'on');
    yline(ax, 0, '--', 'Color', 0.4 .* [1 1 1], ...
        'LineWidth', style.referenceLineWidth, 'HandleVisibility', 'off');

    if fitResult.success
        fitBias = fitResult.con.y - fitResult.incon.y;
        fitMask = fitResult.baseline.y - ...
            0.5 .* (fitResult.con.y + fitResult.incon.y);
        plot(ax, fitResult.x, fitBias, '-', ...
            'Color', style.bias.color, 'LineWidth', style.lineWidth, ...
            'HandleVisibility', 'off');
        plot(ax, fitResult.x, fitMask, '-', ...
            'Color', style.mask.color, 'LineWidth', style.lineWidth, ...
            'HandleVisibility', 'off');
    end

    if opts.showSEM && ~isempty(deltaData.biasX)
        errorbar(ax, deltaData.biasX, deltaData.biasing, deltaData.biasSEM, ...
            'LineStyle', 'none', 'Color', style.bias.color, ...
            'LineWidth', style.lineWidth, 'CapSize', 0, 'HandleVisibility', 'off');
    end
    hBias = scatter(ax, deltaData.biasX, deltaData.biasing, ...
        style.deltaMarkerSize .^ 2, style.bias.color, style.bias.marker, ...
        'filled', ...
        'MarkerFaceColor', style.bias.color, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceAlpha', opts.deltaMarkerFaceAlpha, ...
        'MarkerEdgeAlpha', 1, ...
        'LineWidth', style.lineWidth, 'DisplayName', style.bias.label);

    if opts.showSEM && ~isempty(deltaData.maskX)
        errorbar(ax, deltaData.maskX, deltaData.masking, deltaData.maskSEM, ...
            'LineStyle', 'none', 'Color', style.mask.color, ...
            'LineWidth', style.lineWidth, 'CapSize', 0, 'HandleVisibility', 'off');
    end
    hMask = scatter(ax, deltaData.maskX, deltaData.masking, ...
        style.deltaMarkerSize .^ 2, style.mask.color, style.mask.marker, ...
        'filled', ...
        'MarkerFaceColor', style.mask.color, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceAlpha', opts.deltaMarkerFaceAlpha, ...
        'MarkerEdgeAlpha', 1, ...
        'LineWidth', style.lineWidth, 'DisplayName', style.mask.label);

    xLimits = chooseDeltaXLimits(deltaData, fitResult, opts);
    xlim(ax, xLimits);
    ylim(ax, opts.deltaYLim);
    plotRightEdgeMeanTicks(ax, ...
        [deltaSummary.biasing, deltaSummary.masking], ...
        {style.bias.color, style.mask.color}, style.lineWidth);
    xlabel(ax, 'Gabor contrast (%)');
    ylabel(ax, '\DeltaCorrect (%)');
    title(ax, '');
    legend(ax, [hBias, hMask], {'Biasing', 'Masking'}, ...
        'Location', 'southeast');
    box(ax, 'off');
    axis(ax, 'square');
end


function result = computeAggregateDeltaBiasPermutation(viewData, deltaData, nPermutations, randomSeed)
    rng(double(randomSeed), 'twister');
    contrasts = deltaData.biasX(:);
    nContrasts = numel(contrasts);
    nullDeltaByContrast = nan(nContrasts, nPermutations);
    observedDelta = deltaData.biasing(:);
    nContributing = zeros(nContrasts, 1);

    for contrastIdx = 1:nContrasts
        contrast = contrasts(contrastIdx);
        conIdx = find(viewData.con.x == contrast, 1);
        inconIdx = find(viewData.incon.x == contrast, 1);
        if isempty(conIdx) || isempty(inconIdx)
            continue;
        end
        conSessions = viewData.con.sessionIDs{conIdx};
        inconSessions = viewData.incon.sessionIDs{inconIdx};
        conValues = viewData.con.sessionValues{conIdx};
        inconValues = viewData.incon.sessionValues{inconIdx};
        commonSessions = intersect(conSessions, inconSessions, 'stable');
        nContributing(contrastIdx) = numel(commonSessions);
        if isempty(commonSessions)
            continue;
        end
        sessionNull = nan(numel(commonSessions), nPermutations);
        for sessionIdx = 1:numel(commonSessions)
            sessionID = commonSessions(sessionIdx);
            conSessionIdx = find(conSessions == sessionID, 1);
            inconSessionIdx = find(inconSessions == sessionID, 1);
            conPct = conValues(conSessionIdx);
            inconPct = inconValues(inconSessionIdx);
            [nCon, nIncon] = assumedAggregateDeltaTrialCounts(contrast);
            conCorrect = round(conPct ./ 100 .* nCon);
            inconCorrect = round(inconPct ./ 100 .* nIncon);
            conCorrect = min(max(conCorrect, 0), nCon);
            inconCorrect = min(max(inconCorrect, 0), nIncon);
            totalCorrect = conCorrect + inconCorrect;
            totalTrials = nCon + nIncon;
            pooledOutcomes = [ones(totalCorrect, 1); zeros(totalTrials - totalCorrect, 1)];
            for permIdx = 1:nPermutations
                permOrder = randperm(totalTrials);
                permConCorrect = sum(pooledOutcomes(permOrder(1:nCon)));
                permInconCorrect = totalCorrect - permConCorrect;
                sessionNull(sessionIdx, permIdx) = ...
                    100 .* permConCorrect ./ nCon - ...
                    100 .* permInconCorrect ./ nIncon;
            end
        end
        nullDeltaByContrast(contrastIdx, :) = mean(sessionNull, 1, 'omitnan');
    end

    nullMedian = nan(nContrasts, 1);
    nullLower95 = nan(nContrasts, 1);
    nullUpper95 = nan(nContrasts, 1);
    rawTwoSidedP = nan(nContrasts, 1);
    rawPositiveOneSidedP = nan(nContrasts, 1);
    for contrastIdx = 1:nContrasts
        nullVals = nullDeltaByContrast(contrastIdx, :);
        nullMedian(contrastIdx) = median(nullVals, 'omitnan');
        nullLower95(contrastIdx) = prctile(nullVals, 2.5);
        nullUpper95(contrastIdx) = prctile(nullVals, 97.5);
        pUpper = (1 + sum(nullVals >= observedDelta(contrastIdx))) ./ (nPermutations + 1);
        pLower = (1 + sum(nullVals <= observedDelta(contrastIdx))) ./ (nPermutations + 1);
        rawTwoSidedP(contrastIdx) = min(1, 2 .* min(pUpper, pLower));
        rawPositiveOneSidedP(contrastIdx) = pUpper;
    end

    observedMean = mean(observedDelta, 'omitnan');
    nullMean = mean(nullDeltaByContrast, 1, 'omitnan');
    pUpperOverall = (1 + sum(nullMean >= observedMean)) ./ (nPermutations + 1);
    pLowerOverall = (1 + sum(nullMean <= observedMean)) ./ (nPermutations + 1);

    result = struct();
    result.contrast = contrasts;
    result.observedDeltaBias = observedDelta;
    result.nullDeltaByContrast = nullDeltaByContrast;
    result.nullMedian = nullMedian;
    result.nullLower95 = nullLower95;
    result.nullUpper95 = nullUpper95;
    result.rawTwoSidedP = rawTwoSidedP;
    result.rawPositiveOneSidedP = rawPositiveOneSidedP;
    result.significantUncorrected = rawTwoSidedP < 0.05;
    result.nContributingExperiments = nContributing;
    result.observedMeanDeltaBias = observedMean;
    result.nullMeanDelta = nullMean(:);
    result.nullMeanMedian = median(nullMean, 'omitnan');
    result.nullMeanLower95 = prctile(nullMean, 2.5);
    result.nullMeanUpper95 = prctile(nullMean, 97.5);
    result.rawOverallTwoSidedP = min(1, 2 .* min(pUpperOverall, pLowerOverall));
    result.rawOverallPositiveOneSidedP = pUpperOverall;
    result.overallSignificant = result.rawOverallTwoSidedP < 0.05;
    result.nPermutations = nPermutations;
    result.multipleComparisonCorrection = 'none';
end

function [nCon, nIncon] = assumedAggregateDeltaTrialCounts(contrast)
    if abs(contrast) < eps
        nCon = 40;
        nIncon = 40;
    else
        nCon = 20;
        nIncon = 20;
    end
end

function seed = stableAggregateDeltaPermutationSeed(baseSeed, clusterID)
    values = double(char(string(clusterID)));
    if isnumeric(clusterID)
        values = [values, double(clusterID(:)')];
    end
    seed = mod(double(baseSeed) + sum((1:numel(values)) .* values), 2^31 - 1);
    if seed <= 0
        seed = double(baseSeed);
    end
end

function assertAggregatePermutationMatchesDisplayed(permResult, deltaData, clusterID)
    displayedX = deltaData.biasX(:);
    displayedY = deltaData.biasing(:);
    keep = isfinite(displayedX) & isfinite(displayedY);
    displayedX = displayedX(keep);
    displayedY = displayedY(keep);
    [displayedX, order] = sort(displayedX);
    displayedY = displayedY(order);
    if numel(displayedX) ~= numel(permResult.contrast) || ...
            any(abs(displayedX(:) - permResult.contrast(:)) > 1e-9) || ...
            any(abs(displayedY(:) - permResult.observedDeltaBias(:)) > 1e-9)
        error('plotAggregatedPowerClusterPsychometrics:DeltaPermutationMismatch', ...
            ['Aggregate permutation inputs do not reproduce displayed ' ...
            'merged purple deltaBias points for C%d.'], clusterID);
    end
end

function saveAggregateDeltaPermutationExampleFigure(fig, opts)
    if ~isfield(opts, 'saveDeltaPermutationExamples') || ...
            ~opts.saveDeltaPermutationExamples || ...
            ~isfield(opts, 'deltaPermutationExampleDir') || ...
            isempty(opts.deltaPermutationExampleDir)
        return;
    end
    outputPath = fullfile(opts.deltaPermutationExampleDir, ...
        'PepperR_regularAggregate_permStats.png');
    if ~exist(fileparts(outputPath), 'dir')
        mkdir(fileparts(outputPath));
    end
    if isfile(outputPath)
        return;
    end
    exportgraphics(fig, outputPath, 'Resolution', 200);
    info = dir(outputPath);
    assert(isfile(outputPath) && info.bytes > 0, ...
        'Failed to write aggregate delta permutation example PNG: %s', outputPath);
end
function assertAggregateDeltaPermutationVisualizationAudit(audit, permResult, clusterID)
    if ~isfield(audit, 'visualization') || isempty(audit.visualization)
        error('plotAggregatedPowerClusterPsychometrics:MissingDeltaPermutationVisualizationAudit', ...
            'Missing permutation visualization audit for C%d.', clusterID);
    end
    vis = audit.visualization;
    nExpected = numel(permResult.contrast);
    nNullIntervals = getAuditScalarField(vis, 'nNullIntervals', 0);
    nNullMedians = getAuditScalarField(vis, 'nNullMedians', 0);
    nLabels = getAuditScalarField(vis, 'nLabels', 0);
    nOverall = getAuditScalarField(vis, 'nOverall', 0);
    nOverallNullBands = getAuditScalarField(vis, 'nOverallNullBands', 0);
    nOverallNullMedians = getAuditScalarField(vis, 'nOverallNullMedians', 0);

    if nOverallNullBands > 0 || nOverallNullMedians > 0
        if nNullIntervals ~= 0 || nNullMedians ~= 0 || nLabels ~= 0 || ...
                nOverall < 1 || nOverallNullBands ~= 1 || nOverallNullMedians ~= 1
            error('plotAggregatedPowerClusterPsychometrics:DeltaPermutationVisualizationCountMismatch', ...
                ['Overall-only permutation visualization mismatch for C%d: ' ...
                'contrast intervals %d, medians %d, labels %d, overall labels %d, ' ...
                'overall bands %d, overall medians %d.'], clusterID, ...
                nNullIntervals, nNullMedians, nLabels, nOverall, ...
                nOverallNullBands, nOverallNullMedians);
        end
    elseif nNullIntervals ~= nExpected || nNullMedians ~= nExpected || ...
            nLabels ~= nExpected || nOverall < 1
        error('plotAggregatedPowerClusterPsychometrics:DeltaPermutationVisualizationCountMismatch', ...
            ['Permutation visualization count mismatch for C%d: expected %d, ' ...
            'intervals %d, medians %d, labels %d, overall %d.'], ...
            clusterID, nExpected, nNullIntervals, nNullMedians, nLabels, nOverall);
    end

    if isfield(vis, 'labelsShareY') && ~vis.labelsShareY
        error('plotAggregatedPowerClusterPsychometrics:DeltaPermutationLabelYMismatch', ...
            'Permutation labels do not share one y-coordinate for C%d.', clusterID);
    end
    if isfield(vis, 'allHandlesHidden') && ~vis.allHandlesHidden
        error('plotAggregatedPowerClusterPsychometrics:DeltaPermutationHandleVisibility', ...
            'Permutation visualization handles are not hidden from legend for C%d.', clusterID);
    end
    fprintf(['Added permutation visuals | contrast boxes %d | contrast labels %d | ' ...
        'overall bands %d | overall annotations %d\n'], ...
        nNullIntervals, nLabels, nOverallNullBands, nOverall);
end

function value = getAuditScalarField(auditStruct, fieldName, defaultValue)
    if isfield(auditStruct, fieldName) && ~isempty(auditStruct.(fieldName))
        value = auditStruct.(fieldName);
    else
        value = defaultValue;
    end
end
function validateAggregateDeltaPermutationLegend(ax, clusterID)
    legends = findobj(ancestor(ax, 'figure'), 'Type', 'Legend');
    expected = {'Biasing', 'Masking'};
    legendStrings = {};
    for ii = 1:numel(legends)
        candidateStrings = cellstr(string(legends(ii).String));
        if numel(candidateStrings) == 2 && isequal(candidateStrings(:)', expected)
            legendStrings = candidateStrings;
            break;
        end
    end
    if isempty(legendStrings)
        allStrings = cell(size(legends));
        for ii = 1:numel(legends)
            allStrings{ii} = strjoin(cellstr(string(legends(ii).String)), ', ');
        end
        error('plotAggregatedPowerClusterPsychometrics:DeltaPermutationLegendMismatch', ...
            'Legend mismatch for C%d. Expected Biasing/Masking. Found legends: %s', ...
            clusterID, strjoin(allStrings, ' | '));
    end
    fprintf('Aggregate legend validation passed for C%d: {%s, %s}\n', ...
        clusterID, legendStrings{1}, legendStrings{2});
end
function printAggregateDeltaPermutationSummary(clusterID, nContrasts, meanDelta, pValue, significant, positiveOneSidedP)
    if significant
        sigLabel = '*';
    else
        sigLabel = 'n.s.';
    end
    fprintf('cluster C%s | n contrasts %d | mean DeltaBias %.3f | raw two-sided p %.4g | raw positive one-sided p %.4g | %s\n', ...
        char(string(clusterID)), nContrasts, meanDelta, pValue, positiveOneSidedP, sigLabel);
end

function xLimits = chooseXLimits(viewData, opts)
    if ~isempty(opts.xLim)
        xLimits = opts.xLim;
        return;
    end

    allX = [viewData.baseline.x, viewData.con.x, viewData.incon.x];
    if isempty(allX)
        xLimits = [0 100];
    else
        xLimits = [0, max(100, ceil(max(allX) ./ 5) .* 5)];
    end
end

function xLimits = chooseDeltaXLimits(deltaData, fitResult, opts)
    if ~isempty(opts.xLim)
        xLimits = opts.xLim;
    elseif fitResult.success && ~isempty(fitResult.x)
        xLimits = [0, max(fitResult.x)];
    else
        allX = [deltaData.biasX, deltaData.maskX];
        if isempty(allX)
            xLimits = [0 100];
        else
            xLimits = [0, max(100, ceil(max(allX) ./ 5) .* 5)];
        end
    end
end

function plotRightEdgeMeanTicks(ax, meanValues, colors, lineWidth)
    xLimits = xlim(ax);
    xRange = diff(xLimits);
    if ~isfinite(xRange) || xRange <= 0
        return;
    end

    tickX = [xLimits(2) - 0.050 .* xRange, ...
        xLimits(2) - 0.005 .* xRange];
    for valueIdx = 1:numel(meanValues)
        yValue = meanValues(valueIdx);
        if ~isfinite(yValue)
            continue;
        end

        plot(ax, tickX, [yValue, yValue], '-', ...
            'Color', colors{valueIdx}, ...
            'LineWidth', lineWidth, ...
            'Clipping', 'off', ...
            'HandleVisibility', 'off');
    end
end

function label = clusterLabelToString(clusterID)
    if isnumeric(clusterID) || islogical(clusterID)
        label = num2str(clusterID);
    elseif isstring(clusterID)
        label = char(clusterID);
    elseif iscategorical(clusterID)
        label = char(string(clusterID));
    elseif ischar(clusterID)
        label = clusterID;
    else
        label = char(string(clusterID));
    end
end

function titleText = buildClusterTitle(clusterData)
    clusterLabel = clusterLabelToString(clusterData.clusterID);
    titleText = ['Power cluster ' clusterLabel];
    targetText = '';
    if isfield(clusterData, 'columnTargetLabel')
        targetText = [', Columns ' clusterData.columnTargetLabel];
        titleText = [titleText targetText];
    end

    if ~isfield(clusterData, 'powerRange') || ...
            numel(clusterData.powerRange) ~= 2 || ...
            any(~isfinite(clusterData.powerRange))
        return
    end

    powerRange = sort(clusterData.powerRange(:));
    lowPower = formatPowerValue(powerRange(1));
    highPower = formatPowerValue(powerRange(2));
    nExperiments = NaN;
    if isfield(clusterData, 'sourceBlockIndices')
        nExperiments = numel(clusterData.sourceBlockIndices);
    elseif isfield(clusterData, 'sessionIDs')
        nExperiments = numel(clusterData.sessionIDs);
    end

    if isfinite(nExperiments)
        titleText = sprintf( ...
            'Power cluster %s%s (%s-%s mW, n_{expt}=%d)', ...
            clusterLabel, targetText, lowPower, highPower, nExperiments);
    else
        titleText = sprintf('Power cluster %s%s (%s-%s mW)', ...
            clusterLabel, targetText, lowPower, highPower);
    end
end

function valueText = formatPowerValue(value)
    display = formatPowerMetricsForDisplay(NaN, value);
    valueText = display.Ptotal;
end

function styleAggregateAxes(topAxes, bottomAxes)
    for axesIdx = 1:numel(topAxes)
        axes(topAxes(axesIdx));
        xLimits = xlim(topAxes(axesIdx));
        addSkippedTicks(xLimits(1), xLimits(2), diff(xLimits) ./ 8, 'x');
        addSkippedTicks(0, 100, 10, 'y');
        set(topAxes(axesIdx), ...
            'LineWidth', 2, ...
            'TickDir', 'out', ...
            'TickLength', [0.01 0.01], ...
            'FontName', 'FreeSans');

        axes(bottomAxes(axesIdx));
        xLimits = xlim(bottomAxes(axesIdx));
        yLimits = ylim(bottomAxes(axesIdx));
        addSkippedTicks(xLimits(1), xLimits(2), diff(xLimits) ./ 8, 'x');
        setAlternatingTickLabels(bottomAxes(axesIdx), yLimits, 5, 'y');
        set(bottomAxes(axesIdx), ...
            'LineWidth', 2, ...
            'TickDir', 'out', ...
            'TickLength', [0.01 0.01], ...
            'FontName', 'FreeSans');
    end
end

function setAlternatingTickLabels(ax, limits, interval, axisName)
    ticks = limits(1):interval:limits(2);
    labels = strings(size(ticks));
    zeroIdx = find(abs(ticks) < max(eps(max(abs(limits))), 1e-12), 1);
    if isempty(zeroIdx)
        labelIdx = 1:2:numel(ticks);
    else
        labelIdx = mod(1:numel(ticks), 2) == mod(zeroIdx, 2);
    end
    labels(labelIdx) = compose('%g', ticks(labelIdx));

    if strcmp(axisName, 'x')
        set(ax, 'XLim', limits, 'XTick', ticks, 'XTickLabel', labels);
    else
        set(ax, 'YLim', limits, 'YTick', ticks, 'YTickLabel', labels);
    end
end

function sortedAgg = sortAggregateClusters(agg)
    clusterIDs = arrayfun(@(item) double(item.clusterID), agg);
    targetOrder = ones(size(clusterIDs));
    for idx = 1:numel(agg)
        if isfield(agg(idx), 'targetSortOrder')
            targetOrder(idx) = agg(idx).targetSortOrder;
        end
    end
    [~, sortOrder] = sortrows([clusterIDs(:), targetOrder(:)], [1 2]);
    sortedAgg = agg(sortOrder);
end

function addFigureTitle(titleText)
    if exist('suplabel', 'file') == 2
        [~, titleHandle] = suplabel(titleText, 't', [.1 .1 .82 .84]);
        set(titleHandle, ...
            'FontSize', 16, ...
            'FontWeight', 'normal', ...
            'Interpreter', 'tex');
    else
        annotation(gcf, 'textbox', [0.2 0.925 0.6 0.05], ...
            'String', titleText, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 16, ...
            'Interpreter', 'tex', ...
            'EdgeColor', 'none');
    end
end

function addDeltaSummaryText(ax, deltaSummary)
    axPosition = get(ax, 'Position');
    textWidth = 0.58 .* axPosition(3);
    textX = axPosition(1) + 0.5 .* axPosition(3) - 0.5 .* textWidth;
    textY = max(0.001, axPosition(2) - 0.23);
    annotation(gcf, 'textbox', [textX, textY, textWidth, 0.085], ...
        'String', sprintf('biasing: %.1f%%\nmasking: %.1f%%', ...
        deltaSummary.biasing, deltaSummary.masking), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', ...
        'FontSize', 13.5, ...
        'Interpreter', 'tex', ...
        'EdgeColor', 'none', ...
        'BackgroundColor', 'none', ...
        'FitBoxToText', 'off');
end

function addStimulationStatsText(ax, stats, clusterData)
    axPosition = get(ax, 'Position');
    statsRight = axPosition(1) - 0.012;
    statsX = 0.002;
    statsWidth = max(0.02, statsRight - statsX);
    statsHeight = 0.38;
    statsY = axPosition(2) + 0.5 .* axPosition(4) - 0.5 .* statsHeight;

    clusterText = clusterLabelToString(clusterData.clusterID);
    clusterCount = clusterData.clusterID;
    if isfield(clusterData, 'powerClusterCount')
        clusterCount = clusterData.powerClusterCount;
    end
    columnsText = 'n/a';
    if isfield(clusterData, 'columnTargetLabel')
        columnsText = clusterData.columnTargetLabel;
    end
    baselineText = 'n/a';
    if isfield(clusterData, 'baselineModeSummary')
        baselineText = clusterData.baselineModeSummary;
    elseif isfield(clusterData, 'baselineModes')
        baselineText = summarizeBaselineModesForText(clusterData.baselineModes);
    end
    powerDisplay = formatPowerMetricsForDisplay(stats.roiPowerDensity, ...
        stats.totalPower);
    statsText = sprintf([ ...
        'Cluster: %s/%d\n' ...
        'BL: %s\n' ...
        'Columns_{targeted}: %s\n' ...
        'Columns_{n}: %s\n' ...
        'PD_{DMD} %s mW mm^{-2}\n' ...
        'Area_{ROI} %s mm^2\n' ...
        'Area_{ON} %s mm^2\n' ...
        'sDC %s%%\n' ...
        'tDC %s%%\n' ...
        'PD_{ROI} %s mW mm^{-2}\n' ...
        'P_{total} %s mW'], ...
        clusterText, clusterCount, baselineText, columnsText, ...
        formatRange(stats.columns, 1), ...
        formatRange(stats.projectorPowerDensity, 2), ...
        formatRange(stats.areaROI, 2), ...
        formatRange(stats.areaON, 2), ...
        formatRange(stats.spatialDutyCycle, 1), ...
        formatRange(stats.temporalDutyCycle, 1), ...
        powerDisplay.PDROI, ...
        powerDisplay.Ptotal);

    annotation(gcf, 'textbox', [statsX, statsY, statsWidth, statsHeight], ...
        'String', statsText, ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', 10.5, ...
        'Interpreter', 'tex', ...
        'EdgeColor', 'none', ...
        'BackgroundColor', 'none', ...
        'FitBoxToText', 'off');
end

function textValue = summarizeBaselineModesForText(modes)
    modes = string(modes(:));
    modes = modes(strlength(modes) > 0);
    if isempty(modes)
        textValue = 'n/a';
        return;
    end
    uniqueModes = unique(modes, 'stable');
    if numel(uniqueModes) == 1
        textValue = char(uniqueModes);
        return;
    end
    parts = strings(numel(uniqueModes), 1);
    for idx = 1:numel(uniqueModes)
        parts(idx) = sprintf('%s n=%d', ...
            uniqueModes(idx), sum(modes == uniqueModes(idx)));
    end
    textValue = char(strjoin(parts, ', '));
end

function textValue = formatRange(summary, decimalPlaces)
    if any(~isfinite([summary.min, summary.max]))
        textValue = 'n/a';
        return;
    end
    formatString = sprintf('%%.%df-%%.%df', decimalPlaces, decimalPlaces);
    textValue = sprintf(formatString, summary.min, summary.max);
end

function addFitParameterTable(ax, fitResult)
    if ~fitResult.success
        return;
    end
    if isfield(fitResult, 'modelType') && ...
            strcmp(fitResult.modelType, 'weibullSignedBX0')
        addSignedBX0AggregateFitParameterTable(ax, fitResult.signedBX0);
        return;
    end

    parameterHeaders = {'A', 'B', '\alpha', '\beta'};
    rowLabels = {'Baseline', 'Con-Opto', 'Incon-Opto'};
    rowColors = [0 0 0; 0.55 0 0; 0 0.05 0.45];
    tableValues = [ ...
        fitResult.baseline.params; ...
        fitResult.con.params; ...
        fitResult.incon.params];

    axPosition = get(ax, 'Position');
    tableGap = 0.010;
    maxTableRight = 0.992;
    tableX = axPosition(1) + axPosition(3) + tableGap;
    tableWidth = min(0.22, maxTableRight - tableX);
    if tableWidth < 0.18
        tableWidth = 0.18;
        tableX = max(0.01, maxTableRight - tableWidth);
    end

    tableHeight = 0.44 .* axPosition(4);
    tableY = axPosition(2) + 0.5 .* axPosition(4) - 0.5 .* tableHeight;
    rowHeight = tableHeight ./ 4;
    labelWidth = 0.070;
    labelGap = 0.0015;
    parameterGap = 0.0100;
    valueWidth = (tableWidth - labelWidth - labelGap - ...
        3 .* parameterGap) ./ 4;
    fontSize = 10.5;

    addTableCell(tableX, tableY + 3 .* rowHeight, labelWidth, ...
        rowHeight, '', [0 0 0], fontSize, 'bold', 'left');
    for column = 1:4
        xPosition = tableX + labelWidth + labelGap + ...
            (column - 1) .* (valueWidth + parameterGap);
        addTableCell(xPosition, tableY + 3 .* rowHeight, valueWidth, ...
            rowHeight, parameterHeaders{column}, [0 0 0], ...
            fontSize, 'bold', 'center');
    end

    for row = 1:3
        yPosition = tableY + (3 - row) .* rowHeight;
        addTableCell(tableX, yPosition, labelWidth, rowHeight, ...
            rowLabels{row}, rowColors(row,:), fontSize, 'bold', 'left');
        for column = 1:4
            value = tableValues(row, column);
            if column <= 2
                value = 100 .* value;
            end
            valueText = sprintf('%.1f', value);
            xPosition = tableX + labelWidth + labelGap + ...
                (column - 1) .* (valueWidth + parameterGap);
            addTableCell(xPosition, yPosition, valueWidth, rowHeight, ...
                valueText, rowColors(row,:), fontSize, 'normal', 'center');
        end
    end
end

function addSignedBX0AggregateFitParameterTable(ax, signedBX0)
    if ~isfield(signedBX0, 'fitParams') || numel(signedBX0.fitParams) < 10 || ...
            ~isfield(signedBX0, 'globalDeltaX0') || ~isfinite(signedBX0.globalDeltaX0)
        return;
    end

    params = signedBX0.fitParams(1:10);
    deltaB = params(10);
    deltaX0 = signedBX0.globalDeltaX0;
    tableValues = [ ...
        params(1), 50, params(2), params(3), 0; ...
        params(4), 50 + deltaB, params(5), params(6), -deltaX0; ...
        params(7), 50 - deltaB, params(8), params(9), +deltaX0];
    parameterHeaders = {'A', 'B', '\alpha', '\beta', 'X0'};
    rowLabels = {'Baseline', 'Con-Opto', 'Incon-Opto'};
    rowColors = [0 0 0; 0.55 0 0; 0 0.05 0.45];

    axPosition = get(ax, 'Position');
    tableGap = 0.004;
    maxTableRight = 0.988;
    tableX = axPosition(1) + axPosition(3) + tableGap;
    tableWidth = min(0.205, maxTableRight - tableX);
    if tableWidth < 0.185
        tableWidth = 0.185;
        tableX = max(0.01, maxTableRight - tableWidth);
    end

    tableHeight = 0.44 .* axPosition(4);
    tableY = axPosition(2) + 0.5 .* axPosition(4) - 0.5 .* tableHeight;
    rowHeight = tableHeight ./ 5;
    fontSize = 10.5;
    columnX = [0.00, 0.43, 0.56, 0.69, 0.82, 0.94];
    columnWidth = [0.41, 0.105, 0.105, 0.105, 0.105, 0.055];
    columnX = tableX + tableWidth .* columnX;
    columnWidth = tableWidth .* columnWidth;
    if columnX(end) + columnWidth(end) > 0.99
        error('plotAggregatedPowerClusterPsychometrics:SignedBX0TableClipped', ...
            'Aggregate signed-BX0 X0 table column exceeds the normalized figure boundary.');
    end

    addTableCell(columnX(1), tableY + 4 .* rowHeight, columnWidth(1), ...
        rowHeight, '', [0 0 0], fontSize, 'bold', 'left');
    for column = 1:5
        addTableCell(columnX(column + 1), tableY + 4 .* rowHeight, ...
            columnWidth(column + 1), rowHeight, parameterHeaders{column}, ...
            [0 0 0], fontSize, 'bold', 'center');
    end

    for row = 1:3
        yPosition = tableY + (4 - row) .* rowHeight;
        addTableCell(columnX(1), yPosition, columnWidth(1), rowHeight, ...
            rowLabels{row}, rowColors(row,:), fontSize, 'bold', 'left');
        for column = 1:5
            valueText = formatSignedBX0AggregateParameterValue(...
                tableValues(row, column), parameterHeaders{column});
            addTableCell(columnX(column + 1), yPosition, columnWidth(column + 1), ...
                rowHeight, valueText, rowColors(row,:), fontSize, 'normal', 'center');
        end
    end

    if isfield(signedBX0, 'deltaAICcX0') && isfinite(signedBX0.deltaAICcX0) && ...
            signedBX0.deltaAICcX0 >= 8
        aicColor = [0 0 0];
    else
        aicColor = [0.45 0.45 0.45];
    end
    if isfield(signedBX0, 'deltaAICcX0')
        deltaAICcX0 = signedBX0.deltaAICcX0;
    else
        deltaAICcX0 = NaN;
    end
    addTableCell(tableX, tableY, tableWidth, rowHeight, ...
        sprintf('DeltaAICc_X0    %.1f', deltaAICcX0), ...
        aicColor, fontSize, 'normal', 'left');
end
function valueText = formatSignedBX0AggregateParameterValue(value, parameterName)
    if ~isfinite(value)
        valueText = '';
        return;
    end
    if strcmp(parameterName, 'X0')
        valueText = sprintf('%+.1f', value);
    else
        valueText = sprintf('%.1f', value);
    end
end
function addTableCell(xPosition, yPosition, width, height, textValue, ...
        color, fontSize, fontWeight, horizontalAlignment)
    xPosition = max(0, min(0.999, xPosition));
    yPosition = max(0, min(0.999, yPosition));
    width = max(0.001, min(width, 1 - xPosition));
    height = max(0.001, min(height, 1 - yPosition));
    annotation(gcf, 'textbox', [xPosition, yPosition, width, height], ...
        'String', textValue, ...
        'HorizontalAlignment', horizontalAlignment, ...
        'VerticalAlignment', 'middle', ...
        'FontSize', fontSize, ...
        'FontWeight', fontWeight, ...
        'Interpreter', 'tex', ...
        'Color', color, ...
        'EdgeColor', 'none', ...
        'BackgroundColor', 'none', ...
        'FitBoxToText', 'off');
end


