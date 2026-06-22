function [mdl, objectiveFunction] = getWeibullSignedX0ModelFuncs(common)
%GETWEIBULLSIGNEDX0MODELFUNCS Model and likelihood for signed-X0 Weibull.
%
% Parameter order:
%   1  A_baseline
%   2  alpha_baseline
%   3  beta_baseline
%   4  A_horizontalOpto
%   5  alpha_horizontalOpto
%   6  beta_horizontalOpto
%   7  A_verticalOpto
%   8  alpha_verticalOpto
%   9  beta_verticalOpto
%   10 deltaX0

if isfield(common, 'epsilon')
    epsilon = common.epsilon;
else
    epsilon = 1e-10;
end

mdl = struct();
mdl.signed = @(x, params) signedPredictions(x, params);
mdl.baseline = @(x, params) foldedPredictions(x, params);
mdl.con = @(x, params) foldedPredictions(x, params);
mdl.incon = @(x, params) foldedPredictions(x, params);

objectiveFunction = @(params, data) signedBinomialNLL( ...
    params, data, epsilon);
end

function y = signedPredictions(x, params)
y = struct();
y.pBaseline = weibullSignedX0Mdl(x, params(1), params(2), params(3), 0);
y.pHorizontalOpto = weibullSignedX0Mdl( ...
    x, params(4), params(5), params(6), params(10));
y.pVerticalOpto = weibullSignedX0Mdl( ...
    x, params(7), params(8), params(9), -params(10));
end

function y = foldedPredictions(x, params)
% Convert signed choice-vertical predictions into the existing folded
% percent-correct display convention. c is contrast magnitude.
c = abs(x);

pBLneg = weibullSignedX0Mdl(-c, params(1), params(2), params(3), 0);
pBLpos = weibullSignedX0Mdl(+c, params(1), params(2), params(3), 0);

pHneg = weibullSignedX0Mdl(-c, params(4), params(5), params(6), params(10));
pHpos = weibullSignedX0Mdl(+c, params(4), params(5), params(6), params(10));

pVneg = weibullSignedX0Mdl(-c, params(7), params(8), params(9), -params(10));
pVpos = weibullSignedX0Mdl(+c, params(7), params(8), params(9), -params(10));

y.pcntrl = 0.5 .* ((100 - pBLneg) + pBLpos);
y.pc = 0.5 .* ((100 - pHneg) + pVpos);
y.pic = 0.5 .* (pHpos + (100 - pVneg));
end

function nLL = signedBinomialNLL(params, data, epsilon)
pBaseline = weibullSignedX0Mdl(data.xBaselineChoice, ...
    params(1), params(2), params(3), 0) ./ 100;
pHorizontal = weibullSignedX0Mdl(data.xHorizontalOptoChoice, ...
    params(4), params(5), params(6), params(10)) ./ 100;
pVertical = weibullSignedX0Mdl(data.xVerticalOptoChoice, ...
    params(7), params(8), params(9), -params(10)) ./ 100;

pBaseline = min(max(pBaseline, epsilon), 1 - epsilon);
pHorizontal = min(max(pHorizontal, epsilon), 1 - epsilon);
pVertical = min(max(pVertical, epsilon), 1 - epsilon);

nLL = ...
    conditionNLL(data.successBaselineChoice, ...
        data.sumBaselineChoice, pBaseline) + ...
    conditionNLL(data.successHorizontalOptoChoice, ...
        data.sumHorizontalOptoChoice, pHorizontal) + ...
    conditionNLL(data.successVerticalOptoChoice, ...
        data.sumVerticalOptoChoice, pVertical);
end

function nLL = conditionNLL(successes, nTrials, probability)
failures = nTrials - successes;
nLL = -sum(successes .* log(probability) + ...
    failures .* log(1 - probability));
end
