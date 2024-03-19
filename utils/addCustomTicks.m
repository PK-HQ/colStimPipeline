function addCustomTicks(axisType, minValue, maxValue, skipLabelInterval)
    % axisType: 'x' for X-axis, 'y' for Y-axis
    % minX, maxX: Range for ticks
    % interval: Spacing between ticks
    % skipLabelInterval: Interval to skip labels (e.g., 10 to skip every 10th label)

    % Generate the ticks and labels
    ticks = minValue:skipLabelInterval:maxValue;
    labels = arrayfun(@num2str, ticks, 'UniformOutput', false);

    % Determine which labels to skip
    skipIndices = mod(ticks, skipLabelInterval) == 0;

    % Replace skipped labels with an empty string
    labels(skipIndices) = {''};

    % Apply the ticks and labels to the specified axis
    if strcmp(axisType, 'x')
        set(gca, 'XTick', ticks, 'XTickLabel', labels);
        xlim([minValue, maxValue])
        xtickangle(0)
    elseif strcmp(axisType, 'y')
        set(gca, 'YTick', ticks, 'YTickLabel', labels);
        ylim([minValue, maxValue])
        ytickangle(0)
    else
        error('Invalid axisType. Use ''x'' for X-axis or ''y'' for Y-axis.');
    end
end
