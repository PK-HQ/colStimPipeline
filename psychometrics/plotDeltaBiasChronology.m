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

    experimentNumber = experimentNumber(:);
    deltaBias = deltaBias(:);
    clusterLabels = clusterLabels(:);
    if numel(experimentNumber) ~= numel(deltaBias) || ...
            numel(experimentNumber) ~= numel(clusterLabels)
        error('experimentNumber, deltaBias, and clusterLabels must have the same length.');
    end

    valid = isfinite(experimentNumber) & isfinite(deltaBias) & ...
        isfinite(clusterLabels);
    [experimentNumber, sortOrder] = sort(experimentNumber(valid));
    deltaBias = deltaBias(valid);
    deltaBias = deltaBias(sortOrder);
    clusterLabels = clusterLabels(valid);
    clusterLabels = clusterLabels(sortOrder);
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
        inCluster = clusterLabels == clusterID;
        handles(clusterIdx) = scatter(ax, experimentCount(inCluster), ...
            deltaBias(inCluster), 110, colors(clusterIdx,:), ...
            'filled', 'Marker', 'o', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 2);
        labels{clusterIdx} = sprintf('C%d', clusterID);
    end

    xlabel(ax, 'Experiment count');
    ylabel(ax, '\Delta biasing (%)');
    title(ax, '\Delta biasing by experiment chronology', ...
        'FontWeight', 'normal');
    if ~isempty(handles)
        legend(ax, handles, labels, 'Location', 'best');
    end
    box(ax, 'off');
    axis(ax, 'square');

    [xLimits, xInterval] = chronologyLimits(experimentCount, true);
    [yLimits, yInterval] = positiveBiasLimits(deltaBias);
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
