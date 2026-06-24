function results = validateWeibullSignedBX0Aggregate()
%VALIDATEWEIBULLSIGNEDBX0AGGREGATE Minimal aggregate-path checks for signedBX0.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repoRoot, 'psychometrics', 'weibullModel'));

agg = localSyntheticAggregate();
opts = struct('fitGridPoints', 200, 'maxIterations', 200, 'xLim', [0 100], ...
    'figureVisible', 'off', 'showDeltaPermutationStats', false);
[out, figureHandles] = plotAggregatedPowerClusterPsychometrics(agg, opts);
if ~isempty(figureHandles)
    close(figureHandles(ishandle(figureHandles)));
end

assert(isfield(out(1), 'signedBX0'), 'Aggregate output missing signedBX0 fit.');
assert(strcmp(out(1).signedBX0.modelVersion, 'fullBeta_slopeCap_v1'), ...
    'Aggregate signedBX0 fit used unexpected modelVersion.');
assert(numel(out(1).signedBX0.fitParams) == 11, ...
    'Aggregate M1 must store 11 parameters.');
assert(numel(out(1).signedBX0.noX0FitParams) == 10, ...
    'Aggregate M0 must store 10 parameters.');
assert(out(1).signedBX0.kM0 == 10 && out(1).signedBX0.kM1 == 11, ...
    'Aggregate AICc parameter counts must be 10/11.');
assert(isfield(out(1).signedBX0, 'slopeDiagnostics'), ...
    'Aggregate output missing slope diagnostics.');

results = struct('passed', true, 'signedBX0', out(1).signedBX0);
fprintf('weibullSignedBX0 aggregate validation passed: M0=10, M1=11.\n');
end

function agg = localSyntheticAggregate()
x = [-60 -30 -15 0 15 30 60];
n = 80 .* ones(size(x));
agg = struct();
agg(1).clusterID = 1;
agg(1).signedChoice.baseline = localCounts(x, n, [52 54 49 50 51 53 55]);
agg(1).signedChoice.horizontalOpto = localCounts(x, n, [12 25 38 44 66 78 88]);
agg(1).signedChoice.verticalOpto = localCounts(x, n, [9 20 32 56 63 76 90]);
agg(1).signedChoice.audit = table();
agg(1).horizontal = localView();
agg(1).vertical = localView();
agg(1).merged = localView();
end

function out = localCounts(x, n, pct)
out = struct('x', x, 'nTrials', n, 'successes', round(n .* pct ./ 100));
end

function view = localView()
x = [0 15 30 60];
n = [80 80 80 80];
view.baseline = struct('x', x, 'nTrials', n, 'successes', round(n .* [50 55 65 80] ./ 100));
view.con = struct('x', x, 'nTrials', n, 'successes', round(n .* [55 65 78 90] ./ 100));
view.incon = struct('x', x, 'nTrials', n, 'successes', round(n .* [45 52 62 74] ./ 100));
end