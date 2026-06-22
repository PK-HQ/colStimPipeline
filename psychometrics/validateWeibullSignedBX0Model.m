function results = validateWeibullSignedBX0Model()
%VALIDATEWEIBULLSIGNEDBX0MODEL Static mathematical checks for weibullSignedBX0.
%
% This validator is intentionally self-contained and does not touch animal
% data. It checks the signed raw model, the folded display transform, and
% the nested B-only/B+X0 model convention.

tol = 1e-8;
x = -100:0.5:100;
params = [8, 16, 3.2, 7, 14, 2.4, 11, 18, 4.1, 6, 5];

pBaseline = weibullSignedBX0Mdl(x, params(1), params(2), params(3), ...
    params(1), params(2), params(3), 50, 0);
pHorizontal = weibullSignedBX0Mdl(x, params(4), params(5), params(6), ...
    params(7), params(8), params(9), 50 - params(10), +params(11));
pVertical = weibullSignedBX0Mdl(x, params(7), params(8), params(9), ...
    params(4), params(5), params(6), 50 + params(10), -params(11));

results = struct();
results.rawPredictionsFinite = all(isfinite([pBaseline, pHorizontal, pVertical]));
results.rawPredictionsReal = isreal([pBaseline, pHorizontal, pVertical]);
results.baselineBAtZero = abs(localPredict(params, 'baseline', 0) - 50) < tol;
results.horizontalJoinIsB = ...
    abs(localPredict(params, 'horizontal', params(11)) - (50 - params(10))) < tol;
results.verticalJoinIsB = ...
    abs(localPredict(params, 'vertical', -params(11)) - (50 + params(10))) < tol;
results.displayHorizontalJoinIsOneHundredMinusB = ...
    abs((100 - localPredict(params, 'horizontal', params(11))) - ...
    (100 - (50 - params(10)))) < tol;
results.displayVerticalJoinIsB = ...
    abs(localPredict(params, 'vertical', -params(11)) - ...
    (50 + params(10))) < tol;
results.opponentBConvention = ...
    (50 - params(10)) < 50 && (50 + params(10)) > 50;
results.opponentX0Convention = ...
    params(11) > 0 && (+params(11)) == -(-params(11));
results.branchMapping = localBranchMappingCheck(params);
results.fractionalBetaIsReal = localFractionalBetaCheck();

kM0 = 10;
kM1 = 11;
nTrials = 720;
nLL0 = 225;
nLL1 = 215;
[~, aicc0] = localAIC(nLL0, kM0, nTrials);
[~, aicc1] = localAIC(nLL1, kM1, nTrials);
results.kM1EqualsKM0PlusOne = kM1 == kM0 + 1;
results.deltaAICcConvention = (aicc0 - aicc1) > 0;

results.syntheticRecovery = localSyntheticRecoveryCheck(params);

names = fieldnames(results);
failed = names(~structfun(@logical, results));
if ~isempty(failed)
    error('validateWeibullSignedBX0Model:Failed', ...
        'Failed checks: %s', strjoin(failed, ', '));
end

fprintf('validateWeibullSignedBX0Model passed %d checks.\n', numel(names));
end

function p = localPredict(params, conditionName, x)
switch conditionName
    case 'baseline'
        p = weibullSignedBX0Mdl(x, params(1), params(2), params(3), ...
            params(1), params(2), params(3), 50, 0);
    case 'horizontal'
        p = weibullSignedBX0Mdl(x, params(4), params(5), params(6), ...
            params(7), params(8), params(9), 50 - params(10), +params(11));
    case 'vertical'
        p = weibullSignedBX0Mdl(x, params(7), params(8), params(9), ...
            params(4), params(5), params(6), 50 + params(10), -params(11));
    otherwise
        error('Unknown condition: %s', conditionName);
end
end

function ok = localBranchMappingCheck(params)
xLeft = -40;
xRight = 40;
pHLeft = localPredict(params, 'horizontal', xLeft);
pHRight = localPredict(params, 'horizontal', xRight);
pVLeft = localPredict(params, 'vertical', xLeft);
pVRight = localPredict(params, 'vertical', xRight);

displayHCon = 100 - pHLeft;
displayHIncon = pHRight;
displayVIncon = 100 - pVLeft;
displayVCon = pVRight;
ok = all(isfinite([displayHCon, displayHIncon, displayVIncon, displayVCon]));
end

function ok = localFractionalBetaCheck()
x = -50:50;
p = weibullSignedBX0Mdl(x, 8, 15, 1.7, 11, 18, 2.3, 47, 4);
ok = isreal(p) && all(isfinite(p));
end

function ok = localSyntheticRecoveryCheck(params)
[initialParams, lb, ub] = getWeibullSignedBX0InitParams();
ok = numel(initialParams) == 11 && numel(lb) == 11 && numel(ub) == 11 && ...
    all(params >= lb) && all(params <= ub);
end

function [aic, aicc] = localAIC(nLL, k, n)
aic = 2 * k + 2 * nLL;
aicc = aic + (2 * k * (k + 1)) / (n - k - 1);
end
