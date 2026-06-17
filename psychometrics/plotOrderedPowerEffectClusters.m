function figureHandle = plotOrderedPowerEffectClusters(power, deltaBias, clusterLabels, diagnostics, opts)
% Plot ordered power-band assignments and selected boundaries.

    if nargin < 5 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'figureVisible') || isempty(opts.figureVisible)
        opts.figureVisible = 'on';
    end
    if ~isfield(opts, 'controlMask') || isempty(opts.controlMask)
        opts.controlMask = false(size(power));
    end

    power = power(:);
    deltaBias = deltaBias(:);
    clusterLabels = clusterLabels(:);
    controlMask = logical(opts.controlMask(:));
    if numel(controlMask) ~= numel(power)
        error('opts.controlMask must have one value per session.');
    end
    valid = isfinite(power) & isfinite(deltaBias) & isfinite(clusterLabels);
    experimental = valid & ~controlMask;

    figureHandle = figure( ...
        'Name', 'Ordered power-effect clustering', ...
        'Color', 'w', ...
        'Visible', validatestring(opts.figureVisible, {'on', 'off'}));
    ax = axes(figureHandle);
    hold(ax, 'on');

    clusterIDs = unique(clusterLabels(valid))';
    colors = gray(max(numel(clusterIDs), 2));
    colors = colors(round(linspace(1, size(colors, 1), numel(clusterIDs))), :);
    handles = gobjects(1, numel(clusterIDs));
    labels = cell(1, numel(clusterIDs));
    for ii = 1:numel(clusterIDs)
        clusterID = clusterIDs(ii);
        inCluster = experimental & clusterLabels == clusterID;
        handles(ii) = scatter(ax, power(inCluster), deltaBias(inCluster), ...
            90, colors(ii,:), 'filled', ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 2);
        meanPower = mean(power(inCluster));
        meanBias = mean(deltaBias(inCluster));
        powerRange = [min(power(inCluster)), max(power(inCluster))];
        plot(ax, meanPower, meanBias, 'd', ...
            'Color', 'k', ...
            'MarkerFaceColor', colors(ii,:), ...
            'MarkerEdgeColor', 'k', ...
            'MarkerSize', 12, ...
            'LineWidth', 2, ...
            'HandleVisibility', 'off');
        plot(ax, powerRange, [meanBias meanBias], ':', ...
            'Color', 'k', ...
            'LineWidth', 2.5, ...
            'HandleVisibility', 'off');
        labels{ii} = sprintf('C%d, n=%d', clusterID, sum(inCluster));
    end

    control = valid & controlMask;
    if any(control)
        controlHandle = scatter(ax, power(control), deltaBias(control), ...
            115, [1.00 0.62 0.05], 'o', 'filled', ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 2);
        handles(end + 1) = controlHandle;
        labels{end + 1} = sprintf('45^{\\circ}/135^{\\circ} control, n=%d', ...
            sum(control));
        for ii = 1:numel(clusterIDs)
            inControlCluster = control & clusterLabels == clusterIDs(ii);
            if ~any(inControlCluster)
                continue;
            end
            meanPower = mean(power(inControlCluster));
            meanBias = mean(deltaBias(inControlCluster));
            powerRange = [min(power(inControlCluster)), ...
                max(power(inControlCluster))];
            plot(ax, meanPower, meanBias, 'd', ...
                'Color', 'k', ...
                'MarkerFaceColor', [1.00 0.62 0.05], ...
                'MarkerEdgeColor', 'k', ...
                'MarkerSize', 12, ...
                'LineWidth', 2, ...
                'HandleVisibility', 'off');
            plot(ax, powerRange, [meanBias meanBias], ':', ...
                'Color', [0.85 0.42 0], ...
                'LineWidth', 2.5, ...
                'HandleVisibility', 'off');
        end
    end

    for boundary = diagnostics.boundaries
        xline(ax, boundary, '--k', 'LineWidth', 2.5, ...
            'HandleVisibility', 'off');
    end

    xlabel(ax, 'Total delivered stimulation power (mW)');
    ylabel(ax, '\Delta biasing (%)');
    title(ax, sprintf('Ordered power-effect bands, k=%d, score=%.3f', ...
        diagnostics.chosenK, diagnostics.score), ...
        'FontWeight', 'normal');
    if ~isempty(handles)
        legend(ax, handles, labels, 'Location', 'best');
    end
    box(ax, 'off');
    axis(ax, 'square');
    styleDiagnosticAxes(ax, power(valid), deltaBias(valid), false);
end

function styleDiagnosticAxes(ax, xValues, yValues, integerX)
    [xLimits, xInterval] = diagnosticLimits(xValues, integerX);
    [yLimits, yInterval] = positiveBiasLimits(yValues);
    axes(ax);
    addSkippedTicks(xLimits(1), xLimits(2), xInterval, 'x');
    addSkippedTicks(yLimits(1), yLimits(2), yInterval, 'y');
    upFontSize(24, 0.01);
    set(ax, 'LineWidth', 2, 'TickDir', 'out', ...
        'TickLength', [0.01 0.01], 'FontName', 'FreeSans');
end

function [limits, interval] = positiveBiasLimits(values)
    values = values(isfinite(values));
    if isempty(values)
        limits = [0 10];
        interval = 1;
        return;
    end

    upperValue = max([values; 0]);
    interval = niceStep(max(upperValue, 1) ./ 8);
    upperLimit = max(interval, ceil(upperValue ./ interval) .* interval);
    limits = [0 upperLimit];
end

function [limits, interval] = diagnosticLimits(values, integerAxis)
    values = values(isfinite(values));
    if isempty(values)
        limits = [0 1];
        interval = 0.125;
        return;
    end

    valueRange = max(values) - min(values);
    if valueRange == 0
        valueRange = max(abs(values(1)), 1);
    end
    interval = niceStep(valueRange ./ 8);
    limits = [floor(min(values) ./ interval), ...
        ceil(max(values) ./ interval)] .* interval;
    if integerAxis
        interval = max(1, ceil(interval));
        limits = [floor(min(values) ./ interval), ...
            ceil(max(values) ./ interval)] .* interval;
    end
end

function step = niceStep(rawStep)
    exponent = floor(log10(rawStep));
    fraction = rawStep ./ (10 .^ exponent);
    choices = [1 2 2.5 5 10];
    choiceIdx = find(choices >= fraction, 1);
    step = choices(choiceIdx) .* (10 .^ exponent);
end
