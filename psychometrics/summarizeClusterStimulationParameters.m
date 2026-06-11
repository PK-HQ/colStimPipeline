function stats = summarizeClusterStimulationParameters(bitmapData, blockIndices)
% Summarize stimulation metadata across the sessions in one power band.

    definitions = { ...
        'nColumns', 'columns', 1; ...
        'projectorPowerDensity_mWmm2', 'projectorPowerDensity', 1; ...
        'areaFinalROI', 'areaROI', 1; ...
        'areaPixelsONWithinROI', 'areaON', 1; ...
        'spatialDutyCycleWithinROI', 'spatialDutyCycle', 100; ...
        'temporalDutyCycle', 'temporalDutyCycle', 100; ...
        'meanPowerDensityWithinROI_mWmm2', 'roiPowerDensity', 1; ...
        'totalPowerToOnPixelsWithinROI_mW', 'totalPower', 1};

    stats = struct();
    for ii = 1:size(definitions, 1)
        sourceField = definitions{ii, 1};
        outputField = definitions{ii, 2};
        scale = definitions{ii, 3};
        values = sessionMeans(bitmapData, sourceField, blockIndices) .* scale;
        stats.(outputField) = summarizeValues(values);
    end
end

function values = sessionMeans(bitmapData, fieldName, blockIndices)
    values = nan(numel(blockIndices), 1);
    if ~isfield(bitmapData, fieldName) || isempty(bitmapData.(fieldName))
        return;
    end

    fieldData = bitmapData.(fieldName);
    nBlocks = size(fieldData, ndims(fieldData));
    for ii = 1:numel(blockIndices)
        blockIdx = blockIndices(ii);
        if blockIdx < 1 || blockIdx > nBlocks
            continue;
        end
        indices = repmat({':'}, 1, ndims(fieldData));
        indices{end} = blockIdx;
        blockData = fieldData(indices{:});
        values(ii) = mean(blockData(:), 'omitnan');
    end
end

function summary = summarizeValues(values)
    values = values(isfinite(values));
    summary = struct('min', NaN, 'median', NaN, 'max', NaN);
    if isempty(values)
        return;
    end
    summary.min = min(values);
    summary.median = median(values);
    summary.max = max(values);
end
