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
        horizontalMargins = [0.12 0.20];
        makeSubplot = @(position) subtightplot(2, 3, position, ...
            panelGap, verticalMargins, horizontalMargins);
        topAxes = gobjects(1, 3);
        bottomAxes = gobjects(1, 3);
        deltaSummaries = repmat(struct('biasing', NaN, 'masking', NaN), 1, 3);

        for viewIdx = 1:numel(viewNames)
            viewName = viewNames{viewIdx};
            viewData = agg(clusterIdx).(viewName);

            fitResult = fitAggregateView(viewData, opts);
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
                viewTitles{viewIdx}, style, opts);

            axBottom = makeSubplot(viewIdx + 3);
            bottomAxes(viewIdx) = axBottom;
            plotDeltaPanel(axBottom, deltaData, fitResult, style, opts);
        end

        axes(bottomAxes(3));
        upFontSize(21, 0.01);
        styleAggregateAxes(topAxes, bottomAxes);

        if isfield(agg(clusterIdx), 'stimulationStats')
            addStimulationStatsText(topAxes(1), ...
                agg(clusterIdx).stimulationStats);
        end
        for viewIdx = 1:3
            addDeltaSummaryText(bottomAxes(viewIdx), deltaSummaries(viewIdx));
        end
        addFitParameterTable(topAxes(3), out(clusterIdx).merged.fit);
        addFigureTitle(clusterTitle);
    end
end

function opts = applyDefaults(opts)
    if ~isstruct(opts) || ~isscalar(opts)
        error('opts must be a scalar struct.');
    end

    defaults = struct( ...
        'showSEM', true, ...
        'figureVisible', 'on', ...
        'xLim', [], ...
        'yLim', [0 100], ...
        'deltaYLim', [-75 75], ...
        'fitGridPoints', 401, ...
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

    opts.showSEM = logical(opts.showSEM);
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
    style.baseline.marker = 'o';
    style.baseline.label = 'Baseline';
    style.con.color = conColor;
    style.con.marker = '^';
    style.con.label = 'Con-Opto';
    style.incon.color = inconColor;
    style.incon.marker = 'v';
    style.incon.label = 'Incon-Opto';
    style.bias.color = [127, 0, 255] ./ 255;
    style.bias.marker = 's';
    style.bias.label = 'Biasing';
    style.mask.color = [125, 125, 125] ./ 255;
    style.mask.marker = 's';
    style.mask.label = 'Masking';
    style.lineWidth = 3;
    style.topMarkerSize = 16;
    style.deltaMarkerSize = 20;
    style.referenceLineWidth = 1.5;
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

function plotPsychometricPanel(ax, viewData, fitResult, titleText, style, opts)
    hold(ax, 'on');
    yline(ax, 50, '--', 'Color', 0.4 .* [1 1 1], ...
        'LineWidth', style.referenceLineWidth, 'HandleVisibility', 'off');

    conditionNames = {'baseline', 'con', 'incon'};
    handles = gobjects(1, numel(conditionNames));
    for conditionIdx = 1:numel(conditionNames)
        conditionName = conditionNames{conditionIdx};
        conditionData = viewData.(conditionName);
        conditionStyle = style.(conditionName);

        handles(conditionIdx) = plotAggregatePoints(ax, conditionData, ...
            conditionStyle, opts, style);
        if fitResult.success
            plot(ax, fitResult.x, fitResult.(conditionName).y, ...
                '-', 'Color', conditionStyle.color, ...
                'LineWidth', style.lineWidth, ...
                'HandleVisibility', 'off');
        end
    end

    xLimits = chooseXLimits(viewData, opts);
    xlim(ax, xLimits);
    ylim(ax, opts.yLim);
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
            'MarkerFaceColor', conditionStyle.color, ...
            'MarkerSize', style.topMarkerSize, ...
            'LineWidth', style.lineWidth, ...
            'DisplayName', conditionStyle.label);
        return;
    end

    if opts.showSEM && isfield(conditionData, 'sessionSEM')
        sem = conditionData.sessionSEM;
        errorbar(ax, conditionData.x, conditionData.y, sem, ...
            'LineStyle', 'none', ...
            'Color', conditionStyle.color, ...
            'LineWidth', style.lineWidth, ...
            'CapSize', 0, ...
            'HandleVisibility', 'off');
    end

    handle = plot(ax, conditionData.x, conditionData.y, ...
        conditionStyle.marker, ...
        'LineStyle', 'none', ...
        'Color', conditionStyle.color, ...
        'MarkerFaceColor', conditionStyle.color, ...
        'MarkerEdgeColor', 'k', ...
        'MarkerSize', style.topMarkerSize, ...
        'LineWidth', style.lineWidth, ...
        'DisplayName', conditionStyle.label);
end

function plotDeltaPanel(ax, deltaData, fitResult, style, opts)
    hold(ax, 'on');
    yline(ax, 0, '--', 'Color', 0.4 .* [1 1 1], ...
        'LineWidth', style.referenceLineWidth, 'HandleVisibility', 'off');

    if opts.showSEM && ~isempty(deltaData.biasX)
        errorbar(ax, deltaData.biasX, deltaData.biasing, deltaData.biasSEM, ...
            'LineStyle', 'none', 'Color', style.bias.color, ...
            'LineWidth', style.lineWidth, 'CapSize', 0, 'HandleVisibility', 'off');
    end
    hBias = plot(ax, deltaData.biasX, deltaData.biasing, style.bias.marker, ...
        'LineStyle', 'none', 'Color', style.bias.color, ...
        'MarkerFaceColor', style.bias.color, 'MarkerEdgeColor', 'k', ...
        'MarkerSize', style.deltaMarkerSize, ...
        'LineWidth', style.lineWidth, 'DisplayName', style.bias.label);

    if opts.showSEM && ~isempty(deltaData.maskX)
        errorbar(ax, deltaData.maskX, deltaData.masking, deltaData.maskSEM, ...
            'LineStyle', 'none', 'Color', style.mask.color, ...
            'LineWidth', style.lineWidth, 'CapSize', 0, 'HandleVisibility', 'off');
    end
    hMask = plot(ax, deltaData.maskX, deltaData.masking, style.mask.marker, ...
        'LineStyle', 'none', 'Color', style.mask.color, ...
        'MarkerFaceColor', style.mask.color, 'MarkerEdgeColor', 'k', ...
        'MarkerSize', style.deltaMarkerSize, ...
        'LineWidth', style.lineWidth, 'DisplayName', style.mask.label);

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

    xLimits = chooseDeltaXLimits(deltaData, fitResult, opts);
    xlim(ax, xLimits);
    ylim(ax, opts.deltaYLim);
    xlabel(ax, 'Gabor contrast (%)');
    ylabel(ax, '\DeltaCorrect (%)');
    title(ax, '');
    legend(ax, [hBias, hMask], {'Biasing', 'Masking'}, ...
        'Location', 'southeast');
    box(ax, 'off');
    axis(ax, 'square');
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
            'Power cluster %s (%s-%s mW mm^{-2}, n_{expt}=%d)', ...
            clusterLabel, lowPower, highPower, nExperiments);
    else
        titleText = sprintf('Power cluster %s (%s-%s mW mm^{-2})', ...
            clusterLabel, lowPower, highPower);
    end
end

function valueText = formatPowerValue(value)
    if abs(value) >= 10
        valueText = sprintf('%.1f', value);
    else
        valueText = sprintf('%.2f', value);
    end
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
        addSkippedTicks(yLimits(1), yLimits(2), diff(yLimits) ./ 10, 'y');
        set(bottomAxes(axesIdx), ...
            'LineWidth', 2, ...
            'TickDir', 'out', ...
            'TickLength', [0.01 0.01], ...
            'FontName', 'FreeSans');
    end
end

function sortedAgg = sortAggregateClusters(agg)
    clusterIDs = arrayfun(@(item) double(item.clusterID), agg);
    [~, sortOrder] = sort(clusterIDs, 'ascend');
    sortedAgg = agg(sortOrder);
end

function addFigureTitle(titleText)
    if exist('suplabel', 'file') == 2
        [~, titleHandle] = suplabel(titleText, 't', [.1 .1 .82 .88]);
        set(titleHandle, ...
            'FontSize', 16, ...
            'FontWeight', 'normal', ...
            'Interpreter', 'tex');
    else
        annotation(gcf, 'textbox', [0.2 0.955 0.6 0.035], ...
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

function addStimulationStatsText(ax, stats)
    axPosition = get(ax, 'Position');
    statsRight = axPosition(1) - 0.012;
    statsX = 0.002;
    statsWidth = max(0.02, statsRight - statsX);
    statsHeight = 0.29;
    statsY = axPosition(2) + 0.5 .* axPosition(4) - 0.5 .* statsHeight;

    statsText = sprintf([ ...
        'cols %s\n' ...
        'PD_{DMD} %s mW/mm^2\n' ...
        'Area_{ROI} %s mm^2\n' ...
        'Area_{ON} %s mm^2\n' ...
        'sDC %s%%\n' ...
        'tDC %s%%\n' ...
        'PD_{ROI} %s mW/mm^2\n' ...
        'P_{total} %s mW'], ...
        formatRange(stats.columns, 1), ...
        formatRange(stats.projectorPowerDensity, 2), ...
        formatRange(stats.areaROI, 2), ...
        formatRange(stats.areaON, 2), ...
        formatRange(stats.spatialDutyCycle, 1), ...
        formatRange(stats.temporalDutyCycle, 1), ...
        formatRange(stats.roiPowerDensity, 2), ...
        formatRange(stats.totalPower, 2));

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
