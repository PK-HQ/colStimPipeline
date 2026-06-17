function stats = summarizeClusterStimulationParameters(bitmapData, blockIndices)
% Summarize stimulation metadata across the sessions in one power band.

    blockIndices = blockIndices(:)';
    values = struct( ...
        'columns', nan(numel(blockIndices), 1), ...
        'projectorPowerDensity', nan(numel(blockIndices), 1), ...
        'areaROI', nan(numel(blockIndices), 1), ...
        'areaON', nan(numel(blockIndices), 1), ...
        'spatialDutyCycle', nan(numel(blockIndices), 1), ...
        'temporalDutyCycle', nan(numel(blockIndices), 1), ...
        'roiPowerDensity', nan(numel(blockIndices), 1), ...
        'totalPower', nan(numel(blockIndices), 1));

    for ii = 1:numel(blockIndices)
        metrics = computePowerMetricsFromSource(bitmapData, blockIndices(ii));
        summary = metrics.summary;
        values.columns(ii) = summary.columns;
        values.projectorPowerDensity(ii) = summary.projectorPowerDensity;
        values.areaROI(ii) = summary.areaROI;
        values.areaON(ii) = summary.areaON;
        values.spatialDutyCycle(ii) = summary.spatialDutyCycleFraction * 100;
        values.temporalDutyCycle(ii) = summary.temporalDutyCycleFraction * 100;
        values.roiPowerDensity(ii) = summary.roiPowerDensityRecomputed;
        values.totalPower(ii) = summary.totalPowerRecomputed;
    end

    stats = struct();
    fields = fieldnames(values);
    for ii = 1:numel(fields)
        stats.(fields{ii}) = summarizeValues(values.(fields{ii}));
    end
end

function summary = summarizeValues(values)
    values = values(isfinite(values));
    summary = struct('min', NaN, 'mean', NaN, 'median', NaN, 'max', NaN);
    if isempty(values)
        return;
    end
    summary.min = min(values);
    summary.mean = mean(values);
    summary.median = median(values);
    summary.max = max(values);
end
