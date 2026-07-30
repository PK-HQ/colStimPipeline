function figureHandle = plotPowerBehaviorCurves(power, deltaBias, deltaMask, controlMask)
    power = power(:);
    deltaBias = deltaBias(:);
    deltaMask = deltaMask(:);
    if nargin < 4 || isempty(controlMask)
        controlMask = false(size(power));
    else
        controlMask = logical(controlMask(:));
    end
    if numel(power) ~= numel(deltaBias) || numel(power) ~= numel(deltaMask)
        error('plotPowerBehaviorCurves:InputLengthMismatch', ...
            'power, deltaBias, and deltaMask must have the same length.');
    end
    if numel(controlMask) ~= numel(power)
        error('plotPowerBehaviorCurves:ControlMaskLengthMismatch', ...
            'controlMask must align with power.');
    end

    validBias = isfinite(power) & isfinite(deltaBias);
    validMask = isfinite(power) & isfinite(deltaMask);
    if ~any(validBias) && ~any(validMask)
        error('plotPowerBehaviorCurves:NoFiniteData', ...
            'No finite power/behavior pairs are available to plot.');
    end

    figureHandle = figure('Name', 'Power x behavior fitted curves', ...
        'Color', 'w');
    ax = axes(figureHandle);
    hold(ax, 'on');
    yline(ax, 0, '--', 'Color', [0.5 0.5 0.5], ...
        'LineWidth', 1.5, 'HandleVisibility', 'off');

    biasColor =[127, 0, 255]/255;
    maskColor =[125 125 125]/255;
    controlBiasColor = [1.00 0.62 0.05];
    controlMaskColor = [1.00 1.00 1.00];
    if any(validBias)
        fitSaturatingCurve(power(validBias & ~controlMask), deltaBias(validBias & ~controlMask), ...
            biasColor, 1);
        plotPowerBehaviorPoints(ax, power, deltaBias, validBias & ~controlMask, ...
            's', biasColor, 'Biasing: 0/90');
        plotPowerBehaviorPoints(ax, power, deltaBias, validBias & controlMask, ...
            's', controlBiasColor, 'Biasing: 45/135 control');
    end
    if any(validMask)
        fitSaturatingCurve(power(validMask & ~controlMask), deltaMask(validMask & ~controlMask), ...
            maskColor, 1);
        plotPowerBehaviorPoints(ax, power, deltaMask, validMask & ~controlMask, ...
            'o', maskColor, 'Masking: 0/90');
        plotPowerBehaviorPoints(ax, power, deltaMask, validMask & controlMask, ...
            'o', controlMaskColor, 'Masking: 45/135 control');
    end

    xlabel(ax, 'Total power (mW)');
    ylabel(ax, '\Delta Correct (%)');
    title(ax, 'Power x behavior', 'FontWeight', 'normal');
    legend(ax, 'Location', 'best');
    box(ax, 'off');
    axis(ax, 'square');
    currentXLimits = xlim(ax);
    xlim(ax, [0 currentXLimits(2)]);
    upFontSize(20, 0.02);
end


function plotPowerBehaviorPoints(ax, power, behavior, pointMask, marker, color, label)
    if ~any(pointMask)
        return;
    end
    scatter(ax, power(pointMask), behavior(pointMask), 130, ...
        marker, 'MarkerFaceColor', color, ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.5, ...
        'DisplayName', label);
end
