function display = formatPowerMetricsForDisplay(PDROI, Ptotal)
% Display-only formatting for optical power metrics.
%
% Raw analysis values must stay full precision. This helper only converts
% values to text for figure annotations/tables/audit display fields.

    display = struct();
    display.PDROI = formatMetricValueOrRange(PDROI, 4);
    display.Ptotal = formatMetricValueOrRange(Ptotal, 3);
    display.PDROIUnit = 'mW mm^{-2}';
    display.PtotalUnit = 'mW';
    display.PDROIRoundsPositiveToZero = roundsPositiveToZero(PDROI, 4);
    display.PtotalRoundsPositiveToZero = roundsPositiveToZero(Ptotal, 3);
    display.anyRoundsPositiveToZero = ...
        display.PDROIRoundsPositiveToZero || ...
        display.PtotalRoundsPositiveToZero;
end

function textValue = formatMetricValueOrRange(value, decimalPlaces)
    if isstruct(value)
        textValue = formatSummaryRange(value, decimalPlaces);
        return;
    end

    value = double(value);
    finiteValues = value(isfinite(value));
    if isempty(finiteValues)
        textValue = 'n/a';
        return;
    end

    if numel(finiteValues) == 1
        textValue = sprintf(['%0.' num2str(decimalPlaces) 'f'], ...
            finiteValues);
    else
        textValue = sprintf('%s-%s', ...
            sprintf(['%0.' num2str(decimalPlaces) 'f'], min(finiteValues)), ...
            sprintf(['%0.' num2str(decimalPlaces) 'f'], max(finiteValues)));
    end
end

function textValue = formatSummaryRange(summary, decimalPlaces)
    values = [summary.min, summary.max];
    values = values(isfinite(values));
    if isempty(values)
        textValue = 'n/a';
        return;
    end
    fmt = ['%0.' num2str(decimalPlaces) 'f'];
    if numel(values) == 1 || abs(values(1) - values(end)) <= eps(max(abs(values)))
        textValue = sprintf(fmt, values(1));
    else
        textValue = sprintf('%s-%s', sprintf(fmt, min(values)), ...
            sprintf(fmt, max(values)));
    end
end

function tf = roundsPositiveToZero(value, decimalPlaces)
    if isstruct(value)
        fields = {'min', 'mean', 'median', 'max'};
        tf = false;
        for idx = 1:numel(fields)
            if isfield(value, fields{idx})
                tf = tf || roundsPositiveToZero(value.(fields{idx}), ...
                    decimalPlaces);
            end
        end
        return;
    end

    value = double(value);
    finitePositive = value(isfinite(value) & value > 0);
    if isempty(finitePositive)
        tf = false;
        return;
    end
    roundedValue = round(finitePositive, decimalPlaces);
    tf = any(roundedValue == 0);
end
