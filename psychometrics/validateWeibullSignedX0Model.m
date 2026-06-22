function results = validateWeibullSignedX0Model()
%VALIDATEWEIBULLSIGNEDX0MODEL Lightweight math checks for signed-X0 Weibull.
%
% This validation is intentionally synthetic and does not run the full
% empirical pipeline.

x = -100:0.5:100;
shared = [0.10, 15, 3];
params0 = [shared, shared, shared, 0];
params5 = [shared, shared, shared, 5];

[model, objectiveFunction] = getWeibullSignedX0ModelFuncs(struct());

p0 = weibullSignedX0Mdl(x, shared(1), shared(2), shared(3), 0);
pH5 = weibullSignedX0Mdl(x, shared(1), shared(2), shared(3), 5);
pV5 = weibullSignedX0Mdl(x, shared(1), shared(2), shared(3), -5);
fold0 = model.baseline(abs(x), params0);
fold5 = model.baseline(abs(x), params5);

results = struct();
results.baselineCrossesAtZero = ...
    abs(weibullSignedX0Mdl(0, shared(1), shared(2), shared(3), 0) - 50) < 1e-10;
results.horizontalCrossesAtPlusFive = ...
    abs(weibullSignedX0Mdl(5, shared(1), shared(2), shared(3), 5) - 50) < 1e-10;
results.verticalCrossesAtMinusFive = ...
    abs(weibullSignedX0Mdl(-5, shared(1), shared(2), shared(3), -5) - 50) < 1e-10;
results.realFinite = all(isreal(p0)) && all(isfinite(p0)) && ...
    all(isreal(pH5)) && all(isfinite(pH5)) && ...
    all(isreal(pV5)) && all(isfinite(pV5));
results.mirroredSymmetry = max(abs(p0 - (100 - fliplr(p0)))) < 1e-8;
results.derivedCurvesFinite = all(isfinite(fold0.pcntrl)) && ...
    all(isfinite(fold5.pc)) && all(isfinite(fold5.pic));
results.derivedConInconSeparate = mean(fold5.pc - fold5.pic, 'omitnan') > 0;

data = makeSyntheticData(params5, x);
nLL0 = objectiveFunction(params0, data);
nLL1 = objectiveFunction(params5, data);
kM0 = 9;
kM1 = 10;
nTrials = sum([data.sumBaselineChoice, ...
    data.sumHorizontalOptoChoice, data.sumVerticalOptoChoice]);
[~, aicc0] = localAIC(nLL0, kM0, nTrials);
[~, aicc1] = localAIC(nLL1, kM1, nTrials);
results.identicalObservationCount = nTrials == sum([ ...
    data.sumBaselineChoice, data.sumHorizontalOptoChoice, ...
    data.sumVerticalOptoChoice]);
results.kM1EqualsKM0PlusOne = kM1 == kM0 + 1;
results.deltaAICcConvention = (aicc0 - aicc1) > 0;

fields = fieldnames(results);
for idx = 1:numel(fields)
    assert(results.(fields{idx}), ...
        'validateWeibullSignedX0Model:%s', fields{idx});
end
fprintf('validateWeibullSignedX0Model: all synthetic checks passed.\n');
end

function data = makeSyntheticData(params, x)
xUse = x(1:20:end);
[model] = getWeibullSignedX0ModelFuncs(struct());
pred = model.signed(xUse, params);
data.xBaselineChoice = xUse;
data.xHorizontalOptoChoice = xUse;
data.xVerticalOptoChoice = xUse;
% Use scaled synthetic counts so the AICc convention check has enough
% information to overcome the one-extra-parameter penalty for M1.
data.sumBaselineChoice = 1000 .* ones(size(xUse));
data.sumBaselineChoice(xUse == 0) = 2000;
data.sumHorizontalOptoChoice = 1000 .* ones(size(xUse));
data.sumVerticalOptoChoice = 1000 .* ones(size(xUse));
data.successBaselineChoice = round(pred.pBaseline ./ 100 .* ...
    data.sumBaselineChoice);
data.successHorizontalOptoChoice = round(pred.pHorizontalOpto ./ 100 .* ...
    data.sumHorizontalOptoChoice);
data.successVerticalOptoChoice = round(pred.pVerticalOpto ./ 100 .* ...
    data.sumVerticalOptoChoice);
end

function [aic, aicc] = localAIC(nLL, k, n)
aic = 2 .* k + 2 .* nLL;
if n > k + 1
    aicc = aic + (2 .* k .* (k + 1)) ./ (n - k - 1);
else
    aicc = Inf;
end
end
