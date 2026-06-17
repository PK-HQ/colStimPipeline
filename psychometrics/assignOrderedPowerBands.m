function clusterLabels = assignOrderedPowerBands(power, boundaries)
% Assign powers to ordered bands using fixed boundary values.

    power = power(:);
    boundaries = sort(boundaries(:)');
    clusterLabels = nan(size(power));
    valid = isfinite(power) & power >= 0;
    clusterLabels(valid) = 1;
    for boundary = boundaries
        clusterLabels(valid & power > boundary) = ...
            clusterLabels(valid & power > boundary) + 1;
    end
end
