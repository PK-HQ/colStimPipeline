function figureHandle = plotDeltaBiasChronology(experimentNumber, deltaBias, clusterLabels, opts)
% Plot behavioral bias in chronological experiment order without clustering.

    if nargin < 3 || isempty(clusterLabels)
        clusterLabels = ones(size(deltaBias));
    end
    if nargin < 4 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'figureVisible') || isempty(opts.figureVisible)
        opts.figureVisible = 'on';
    end
    if ~isfield(opts, 'controlMask') || isempty(opts.controlMask)
        opts.controlMask = false(size(deltaBias));
    end

    experimentNumber = experimentNumber(:);
    deltaBias = deltaBias(:);
    clusterLabels = clusterLabels(:);
    controlMask = logical(opts.controlMask(:));
    if numel(experimentNumber) ~= numel(deltaBias) || ...
            numel(experimentNumber) ~= numel(clusterLabels) || ...
            numel(experimentNumber) ~= numel(controlMask)
        error('All chronology inputs must have the same length.');
    end

    valid = isfinite(experimentNumber) & isfinite(deltaBias) & ...
        isfinite(clusterLabels);
    [experimentNumber, sortOrder] = sort(experimentNumber(valid));
    deltaBias = deltaBias(valid);
    deltaBias = deltaBias(sortOrder);
    clusterLabels = clusterLabels(valid);
    clusterLabels = clusterLabels(sortOrder);
    controlMask = controlMask(valid);
    controlMask = controlMask(sortOrder);
    experimentCount = (1:numel(experimentNumber))';

    figureHandle = figure( ...
        'Name', 'Delta biasing by experiment number', ...
        'Color', 'w', ...
        'Visible', validatestring(opts.figureVisible, {'on', 'off'}));
    ax = axes(figureHandle);
    hold(ax, 'on');
    yline(ax, 0, '--', 'Color', 0.4 .* [1 1 1], ...
        'LineWidth', 2, 'HandleVisibility', 'off');
    clusterIDs = unique(clusterLabels, 'sorted')';
    colors = gray(max(numel(clusterIDs), 2));
    colors = colors(round(linspace(1, size(colors, 1), numel(clusterIDs))), :);
    handles = gobjects(1, numel(clusterIDs));
    labels = cell(1, numel(clusterIDs));
    for clusterIdx = 1:numel(clusterIDs)
        clusterID = clusterIDs(clusterIdx);
        inCluster = clusterLabels == clusterID & ~controlMask;
        handles(clusterIdx) = scatter(ax, experimentCount(inCluster), ...
            deltaBias(inCluster), 110, colors(clusterIdx,:), ...
            'filled', 'Marker', 'o', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 2);
        labels{clusterIdx} = sprintf('C%d', clusterID);
    end
    if any(controlMask)
        handles(end + 1) = scatter(ax, experimentCount(controlMask), ...
            deltaBias(controlMask), 110, [1.00 0.62 0.05], ...
            'filled', 'Marker', 'o', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 2);
        labels{end + 1} = '45^{\circ}/135^{\circ} control';
    end

    xlabel(ax, 'Experiment count');
    ylabel(ax, '\Delta biasing (%)');
    title(ax, '\Delta biasing by experiment chronology', ...
        'FontWeight', 'normal');
    if ~isempty(handles)
        legend(ax, handles, labels, 'Location', 'southeast');
    end
    box(ax, 'off');
    axis(ax, 'square');

    [xLimits, xInterval] = chronologyLimits(experimentCount, true);
    [yLimits, yInterval] = signedBiasLimits(deltaBias);
    axes(ax);
    addSkippedTicks(xLimits(1), xLimits(2), xInterval, 'x');
    setAlternatingTickLabels(ax, yLimits, yInterval);
    upFontSize(24, 0.01);
    set(ax, 'LineWidth', 2, 'TickDir', 'out', ...
        'TickLength', [0.01 0.01], 'FontName', 'FreeSans');
end

function [limits, interval] = signedBiasLimits(values)
    values = values(isfinite(values));
    if isempty(values)
        limits = [-10 10];
        interval = 2;
        return;
    end
    valueRange = max(values) - min(values);
    interval = niceStep(max(valueRange, 1) ./ 10);
    limits = [min(0, floor(min(values) ./ interval) .* interval), ...
        max(0, ceil(max(values) ./ interval) .* interval)];
    if limits(1) == limits(2)
        limits = limits + [-interval interval];
    end
end

function setAlternatingTickLabels(ax, limits, interval)
    ticks = limits(1):interval:limits(2);
    labels = strings(size(ticks));
    zeroIdx = find(abs(ticks) < max(eps(max(abs(limits))), 1e-12), 1);
    if isempty(zeroIdx)
        labelIdx = 1:2:numel(ticks);
    else
        labelIdx = mod(1:numel(ticks), 2) == mod(zeroIdx, 2);
    end
    labels(labelIdx) = compose('%g', ticks(labelIdx));
    set(ax, 'YLim', limits, 'YTick', ticks, 'YTickLabel', labels);
end

function [limits, interval] = chronologyLimits(values, integerAxis)
    if isempty(values)
        limits = [0 1];
        interval = 1;
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
        limits = [max(1, min(values)), ...
            ceil(max(values) ./ interval) .* interval];
    end
end

function step = niceStep(rawStep)
    exponent = floor(log10(rawStep));
    fraction = rawStep ./ (10 .^ exponent);
    choices = [1 2 2.5 5 10];
    choiceIdx = find(choices >= fraction, 1);
    step = choices(choiceIdx) .* (10 .^ exponent);
end
