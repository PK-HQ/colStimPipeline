function summaryTable = summarizeSignedBX0FitDiagnostics(mdl, experimentID)
%SUMMARIZESIGNEDBX0FITDIAGNOSTICS Print compact weibullSignedBX0 fit diagnostics.
% The table uses stored individual-session M1 fit diagnostics and does not
% refit or alter model results.

if ~isfield(mdl, 'signedBX0')
    error('summarizeSignedBX0FitDiagnostics:MissingSignedBX0', ...
        'Input mdl does not contain mdl.signedBX0.');
end
s = mdl.signedBX0;
if nargin < 2 || isempty(experimentID)
    nRows = size(s.fitParams, 1);
    experimentID = arrayfun(@(idx) sprintf('session_%d', idx), ...
        (1:nRows)', 'UniformOutput', false);
elseif ischar(experimentID)
    experimentID = cellstr(experimentID);
elseif isstring(experimentID)
    experimentID = cellstr(experimentID(:));
else
    experimentID = experimentID(:);
end

nRows = size(s.fitParams, 1);
experimentID = localPadIDs(experimentID, nRows);
sessionRow = (1:nRows)';

maxSlopeOverall = localField(s, 'maxSlopeOverall', nRows, NaN);
slopeCapActive = localField(s, 'slopeCapActive', nRows, false);
noX0MaxSlopeOverall = localField(s, 'noX0MaxSlopeOverall', nRows, NaN);
noX0SlopeCapActive = localField(s, 'noX0SlopeCapActive', nRows, false);
boundHit = localField(s, 'boundHit', nRows, false);
exitFlag = localField(s, 'exitFlag', nRows, NaN);
nStarts = localField(s, 'nStarts', nRows, NaN);
bestStartIndex = localField(s, 'bestStartIndex', nRows, NaN);
deltaAICcX0 = localField(s, 'deltaAICcX0', nRows, NaN);
akaikeWeightBX0 = localField(s, 'akaikeWeightBX0', nRows, NaN);
deltaB = localField(s, 'deltaB', nRows, NaN);
deltaX0 = localField(s, 'deltaX0', nRows, NaN);

boundHitFields = repmat({''}, nRows, 1);
if isfield(s, 'boundHitFields')
    for idx = 1:min(nRows, numel(s.boundHitFields))
        fields = s.boundHitFields{idx};
        if iscell(fields)
            boundHitFields{idx} = strjoin(fields, ',');
        elseif ischar(fields)
            boundHitFields{idx} = fields;
        end
    end
end

summaryTable = table(sessionRow, experimentID, deltaB, deltaX0, ...
    deltaAICcX0, akaikeWeightBX0, maxSlopeOverall, slopeCapActive, ...
    noX0MaxSlopeOverall, noX0SlopeCapActive, boundHit, boundHitFields, ...
    nStarts, bestStartIndex, exitFlag);

fprintf('weibullSignedBX0 diagnostics: %d sessions, %d slope-cap hits, %d bound hits.\n', ...
    nRows, sum(logical(slopeCapActive)), sum(logical(boundHit)));
if any(logical(slopeCapActive))
    fprintf('Slope-cap sessions: %s\n', mat2str(sessionRow(logical(slopeCapActive))'));
end
if any(logical(boundHit))
    fprintf('Bound-hit sessions: %s\n', mat2str(sessionRow(logical(boundHit))'));
end
end

function values = localField(s, fieldName, nRows, defaultValue)
if isfield(s, fieldName)
    values = s.(fieldName);
    values = values(:);
else
    values = repmat(defaultValue, nRows, 1);
end
if numel(values) < nRows
    values(end+1:nRows, 1) = defaultValue;
end
values = values(1:nRows);
end

function ids = localPadIDs(ids, nRows)
ids = ids(:);
if numel(ids) < nRows
    for idx = numel(ids)+1:nRows
        ids{idx, 1} = sprintf('session_%d', idx);
    end
end
ids = ids(1:nRows);
end