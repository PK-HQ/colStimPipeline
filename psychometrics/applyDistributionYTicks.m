function tickInfo = applyDistributionYTicks(ax, mode)
% Apply clean y ticks while preserving annotation headroom.

    limits = ylim(ax);
    switch lower(mode)
        case 'performance'
            ticks = 0:20:100;
        case 'delta'
            ticks = cleanTicks(limits, true);
        case 'parameter'
            ticks = cleanTicks(limits, false);
        otherwise
            error('applyDistributionYTicks:UnknownMode', ...
                'Unknown distribution tick mode: %s', mode);
    end

    yticks(ax, ticks);
    ytickangle(ax, 0);
    tickInfo = struct('limits', limits, 'ticks', ticks);
end

function ticks = cleanTicks(limits, includeZero)
    span = max(diff(limits), eps);
    step = niceStep(span ./ 6);
    firstTick = ceil(limits(1) ./ step) .* step;
    lastTick = floor(limits(2) ./ step) .* step;
    ticks = firstTick:step:lastTick;

    if includeZero && limits(1) <= 0 && limits(2) >= 0 && ...
            ~any(abs(ticks) < step .* 1e-8)
        ticks = sort([ticks 0]);
    end
    if isempty(ticks)
        ticks = limits;
    end
    ticks(abs(ticks) < step .* 1e-10) = 0;
end

function step = niceStep(rawStep)
    magnitude = 10 .^ floor(log10(max(rawStep, eps)));
    scaled = rawStep ./ magnitude;
    if scaled <= 1
        multiplier = 1;
    elseif scaled <= 2
        multiplier = 2;
    elseif scaled <= 5
        multiplier = 5;
    else
        multiplier = 10;
    end
    step = multiplier .* magnitude;
end
