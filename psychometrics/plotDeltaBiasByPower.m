function figureHandle = plotDeltaBiasByPower(power, deltaBias, opts)
% Plot behavioral bias against total delivered stimulation power.

    if nargin < 3 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'figureVisible') || isempty(opts.figureVisible)
        opts.figureVisible = 'on';
    end

    power = power(:);
    deltaBias = deltaBias(:);
    if numel(power) ~= numel(deltaBias)
        error('power and deltaBias must have the same length.');
    end
    valid = isfinite(power) & isfinite(deltaBias);

    figureHandle = figure( ...
        'Name', 'Delta biasing by stimulation power', ...
        'Color', 'w', ...
        'Visible', validatestring(opts.figureVisible, {'on', 'off'}));
    ax = axes(figureHandle);
    hold(ax, 'on');
    scatter(ax, power(valid), deltaBias(valid), 110, [0.35 0.35 0.35], ...
        'filled', 'Marker', 'o', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 2);

    xlabel(ax, 'Total delivered stimulation power (mW)');
    ylabel(ax, '\Delta biasing (%)');
    title(ax, '\Delta biasing by stimulation power', ...
        'FontWeight', 'normal');
    box(ax, 'off');
    axis(ax, 'square');

    [xLimits, xInterval] = dataLimits(power(valid));
    [yLimits, yInterval] = positiveBiasLimits(deltaBias(valid));
    axes(ax);
    addSkippedTicks(xLimits(1), xLimits(2), xInterval, 'x');
    addSkippedTicks(yLimits(1), yLimits(2), yInterval, 'y');
    upFontSize(24, 0.01);
    set(ax, 'LineWidth', 2, 'TickDir', 'out', ...
        'TickLength', [0.01 0.01], 'FontName', 'FreeSans');
end

function [limits, interval] = dataLimits(values)
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
    if limits(1) == limits(2)
        limits = limits + [-interval interval];
    end
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

function step = niceStep(rawStep)
    exponent = floor(log10(rawStep));
    fraction = rawStep ./ (10 .^ exponent);
    choices = [1 2 2.5 5 10];
    choiceIdx = find(choices >= fraction, 1);
    step = choices(choiceIdx) .* (10 .^ exponent);
end
