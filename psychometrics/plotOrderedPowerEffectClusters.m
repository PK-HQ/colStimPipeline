function figureHandle = plotOrderedPowerEffectClusters(power, deltaBias, clusterLabels, diagnostics, opts)
% Plot ordered power-band assignments and selected boundaries.

    if nargin < 5 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'figureVisible') || isempty(opts.figureVisible)
        opts.figureVisible = 'on';
    end

    power = power(:);
    deltaBias = deltaBias(:);
    clusterLabels = clusterLabels(:);
    valid = isfinite(power) & isfinite(deltaBias) & isfinite(clusterLabels);

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
        inCluster = valid & clusterLabels == clusterID;
        handles(ii) = scatter(ax, power(inCluster), deltaBias(inCluster), ...
            90, colors(ii,:), 'filled', ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 2);
        medianPower = median(power(inCluster));
        medianBias = median(deltaBias(inCluster));
        powerRange = [min(power(inCluster)), max(power(inCluster))];
        plot(ax, medianPower, medianBias, 'd', ...
            'Color', 'k', ...
            'MarkerFaceColor', colors(ii,:), ...
            'MarkerEdgeColor', 'k', ...
            'MarkerSize', 12, ...
            'LineWidth', 2, ...
            'HandleVisibility', 'off');
        plot(ax, powerRange, [medianBias medianBias], ':', ...
            'Color', 'k', ...
            'LineWidth', 2.5, ...
            'HandleVisibility', 'off');
        labels{ii} = sprintf('C%d, n=%d', clusterID, sum(inCluster));
    end

    for boundary = diagnostics.boundaries
        xline(ax, boundary, '--k', 'LineWidth', 2.5, ...
            'HandleVisibility', 'off');
    end

    xlabel(ax, 'Mean power density within ROI (mW/mm^2)');
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
    [yLimits, yInterval] = diagnosticLimits(yValues, false);
    axes(ax);
    addSkippedTicks(xLimits(1), xLimits(2), xInterval, 'x');
    addSkippedTicks(yLimits(1), yLimits(2), yInterval, 'y');
    upFontSize(24, 0.01);
    set(ax, 'LineWidth', 2, 'TickDir', 'out', ...
        'TickLength', [0.01 0.01], 'FontName', 'FreeSans');
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
