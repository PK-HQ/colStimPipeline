function metrics = computePowerMetricsFromSource(bitmapData, blockIdx, opts)
% Recompute canonical optical power metrics from source bitmap fields.
%
% Canonical definitions:
%   sDC   = AreaON / AreaROI
%   PDROI = PDDMD * sDC * tDC
%   Ptotal = PDDMD * AreaON * tDC
%
% Duty-cycle fields are normalized to fractions before multiplication. A
% stored value > 1 is interpreted as percent (14.5 -> 0.145); 1 is 100%.

    if nargin < 3 || isempty(opts)
        opts = struct();
    end
    opts = applyDefaults(opts);

    values = struct();
    values.PDDMD = blockValues(bitmapData, 'projectorPowerDensity_mWmm2', blockIdx);
    values.AreaROI = blockValues(bitmapData, 'areaFinalROI', blockIdx);
    values.AreaON = blockValues(bitmapData, 'areaPixelsONWithinROI', blockIdx);
    values.sDCStored = blockValues(bitmapData, 'spatialDutyCycleWithinROI', blockIdx);
    values.tDCStored = blockValues(bitmapData, 'temporalDutyCycle', blockIdx);
    values.PDROIStored = blockValues(bitmapData, ...
        'meanPowerDensityWithinROI_mWmm2', blockIdx);
    values.PtotalStored = blockValues(bitmapData, ...
        'totalPowerToOnPixelsWithinROI_mW', blockIdx);
    values.Columns = blockValues(bitmapData, 'nColumns', blockIdx);

    [PDDMD, AreaROI, AreaON, sDCStored, tDCStored, PDROIStored, ...
        PtotalStored] = alignVectors(values.PDDMD, values.AreaROI, ...
        values.AreaON, values.sDCStored, values.tDCStored, ...
        values.PDROIStored, values.PtotalStored);

    sDCFromArea = AreaON ./ AreaROI;
    sDCFromArea(~isfinite(sDCFromArea)) = NaN;
    sDCStoredFraction = normalizeDutyCycleFraction(sDCStored);
    tDCFraction = normalizeDutyCycleFraction(tDCStored);

    sDCCanonical = sDCFromArea;
    missingAreaRatio = ~isfinite(sDCCanonical);
    sDCCanonical(missingAreaRatio) = sDCStoredFraction(missingAreaRatio);

    PDROIRecomputed = PDDMD .* sDCCanonical .* tDCFraction;
    PtotalRecomputed = PDDMD .* AreaON .* tDCFraction;
    PtotalFromPDROI = PDROIRecomputed .* AreaROI;

    PDROIDiff = PDROIRecomputed - PDROIStored;
    PtotalDiff = PtotalRecomputed - PtotalStored;
    sDCDiff = sDCStoredFraction - sDCFromArea;
    PtotalIdentityDiff = PtotalRecomputed - PtotalFromPDROI;
    PDROIRelErr = relativeError(PDROIDiff, PDROIStored);
    PtotalRelErr = relativeError(PtotalDiff, PtotalStored);
    sDCRelErr = relativeError(sDCDiff, sDCFromArea);
    PtotalIdentityRelErr = relativeError(PtotalIdentityDiff, PtotalRecomputed);

    PDROIPass = isWithinTolerance(PDROIDiff, PDROIRelErr, opts);
    PtotalPass = isWithinTolerance(PtotalDiff, PtotalRelErr, opts);
    sDCPass = isWithinTolerance(sDCDiff, sDCRelErr, opts) | ...
        isnan(sDCStoredFraction) | isnan(sDCFromArea);
    PtotalIdentityPass = isWithinTolerance( ...
        PtotalIdentityDiff, PtotalIdentityRelErr, opts);

    metrics = struct();
    metrics.blockIdx = blockIdx;
    metrics.source = values;
    metrics.vector = struct( ...
        'PDDMD', PDDMD, ...
        'AreaROI', AreaROI, ...
        'AreaON', AreaON, ...
        'sDCStoredFraction', sDCStoredFraction, ...
        'sDCFromArea', sDCFromArea, ...
        'sDCDiff', sDCDiff, ...
        'sDCRelErr', sDCRelErr, ...
        'sDCCanonical', sDCCanonical, ...
        'tDCFraction', tDCFraction, ...
        'PDROIStored', PDROIStored, ...
        'PDROIRecomputed', PDROIRecomputed, ...
        'PDROIDiff', PDROIDiff, ...
        'PDROIRelErr', PDROIRelErr, ...
        'PtotalStored', PtotalStored, ...
        'PtotalRecomputed', PtotalRecomputed, ...
        'PtotalFromPDROI', PtotalFromPDROI, ...
        'PtotalDiff', PtotalDiff, ...
        'PtotalRelErr', PtotalRelErr, ...
        'PtotalIdentityDiff', PtotalIdentityDiff, ...
        'PtotalIdentityRelErr', PtotalIdentityRelErr);
    metrics.summary = struct( ...
        'columns', mean(values.Columns, 'omitnan'), ...
        'projectorPowerDensity', mean(PDDMD, 'omitnan'), ...
        'areaROI', mean(AreaROI, 'omitnan'), ...
        'areaON', mean(AreaON, 'omitnan'), ...
        'spatialDutyCycleFraction', mean(sDCCanonical, 'omitnan'), ...
        'spatialDutyCycleStoredFraction', mean(sDCStoredFraction, 'omitnan'), ...
        'temporalDutyCycleFraction', mean(tDCFraction, 'omitnan'), ...
        'roiPowerDensityStored', mean(PDROIStored, 'omitnan'), ...
        'roiPowerDensityRecomputed', mean(PDROIRecomputed, 'omitnan'), ...
        'totalPowerStored', mean(PtotalStored, 'omitnan'), ...
        'totalPowerRecomputed', mean(PtotalRecomputed, 'omitnan'), ...
        'totalPowerFromPDROI', mean(PtotalFromPDROI, 'omitnan'), ...
        'maxAbsSDCDiff', max(abs(sDCDiff), [], 'omitnan'), ...
        'maxRelSDCError', max(sDCRelErr, [], 'omitnan'), ...
        'maxAbsPDROIDiff', max(abs(PDROIDiff), [], 'omitnan'), ...
        'maxRelPDROIError', max(PDROIRelErr, [], 'omitnan'), ...
        'maxAbsPtotalDiff', max(abs(PtotalDiff), [], 'omitnan'), ...
        'maxRelPtotalError', max(PtotalRelErr, [], 'omitnan'), ...
        'maxAbsPtotalIdentityDiff', max(abs(PtotalIdentityDiff), [], 'omitnan'), ...
        'maxRelPtotalIdentityError', max(PtotalIdentityRelErr, [], 'omitnan'), ...
        'sDCPass', all(sDCPass), ...
        'PDROIPass', all(PDROIPass | isnan(PDROIStored)), ...
        'PtotalPass', all(PtotalPass | isnan(PtotalStored)), ...
        'PtotalIdentityPass', all(PtotalIdentityPass));
    metrics.pass = metrics.summary.sDCPass && ...
        metrics.summary.PDROIPass && metrics.summary.PtotalPass && ...
        metrics.summary.PtotalIdentityPass;
end

function opts = applyDefaults(opts)
    defaults = struct('absoluteTolerance', 1e-9, 'relativeTolerance', 1e-6);
    names = fieldnames(defaults);
    for idx = 1:numel(names)
        if ~isfield(opts, names{idx}) || isempty(opts.(names{idx}))
            opts.(names{idx}) = defaults.(names{idx});
        end
    end
end

function values = blockValues(S, fieldName, blockIdx)
    values = NaN;
    if ~isstruct(S) || ~isfield(S, fieldName) || isempty(S.(fieldName)) || ...
            ~isfinite(blockIdx)
        return;
    end
    data = S.(fieldName);
    if isempty(data)
        return;
    end
    nDims = ndims(data);
    nBlocks = size(data, nDims);
    if blockIdx < 1 || blockIdx > nBlocks
        return;
    end
    indices = repmat({':'}, 1, nDims);
    indices{nDims} = blockIdx;
    blockData = data(indices{:});
    values = double(blockData(:));
end

function varargout = alignVectors(varargin)
    lengths = cellfun(@numel, varargin);
    targetLength = max(lengths);
    if isempty(targetLength) || targetLength < 1
        targetLength = 1;
    end
    varargout = cell(size(varargin));
    for idx = 1:numel(varargin)
        x = double(varargin{idx}(:));
        if isempty(x)
            x = NaN;
        end
        if numel(x) == targetLength
            varargout{idx} = x;
        elseif numel(x) == 1
            varargout{idx} = repmat(x, targetLength, 1);
        else
            padded = nan(targetLength, 1);
            padded(1:min(numel(x), targetLength)) = ...
                x(1:min(numel(x), targetLength));
            varargout{idx} = padded;
        end
    end
end

function fraction = normalizeDutyCycleFraction(values)
    fraction = double(values);
    percentMask = isfinite(fraction) & abs(fraction) > 1;
    fraction(percentMask) = fraction(percentMask) ./ 100;
end

function relErr = relativeError(diffValue, reference)
    denom = max(abs(reference), eps);
    relErr = abs(diffValue) ./ denom;
    relErr(~isfinite(reference)) = NaN;
end

function pass = isWithinTolerance(absDiff, relErr, opts)
    pass = abs(absDiff) <= opts.absoluteTolerance | ...
        relErr <= opts.relativeTolerance;
    pass(~isfinite(absDiff) & ~isfinite(relErr)) = true;
end
