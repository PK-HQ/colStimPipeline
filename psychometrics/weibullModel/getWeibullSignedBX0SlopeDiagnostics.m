function diagnostics = getWeibullSignedBX0SlopeDiagnostics(params, maxSlopePctPerContrast)
%GETWEIBULLSIGNEDBX0SLOPEDIAGNOSTICS Six physical half-slopes for BX0.

if nargin < 2 || isempty(maxSlopePctPerContrast)
    maxSlopePctPerContrast = 5.0;
end
params = params(:)';
if numel(params) ~= 11
    error('getWeibullSignedBX0SlopeDiagnostics:InvalidParameterCount', ...
        'Expected 11 weibullSignedBX0 parameters; got %d.', numel(params));
end

A_BL = params(1);
alpha_BL = params(2);
beta_BL = params(3);
A_con = params(4);
alpha_con = params(5);
beta_con = params(6);
A_incon = params(7);
alpha_incon = params(8);
beta_incon = params(9);
deltaB = params(10);

B_H = 50 - deltaB;
B_V = 50 + deltaB;

ampBaselineLeft = 50 - A_BL;
ampBaselineRight = 50 - A_BL;
ampHLeft = B_H - A_con;
ampHRight = (100 - A_incon) - B_H;
ampVLeft = B_V - A_incon;
ampVRight = (100 - A_con) - B_V;

amplitudes = [ampBaselineLeft, ampBaselineRight, ...
    ampHLeft, ampHRight, ampVLeft, ampVRight];
alphas = [alpha_BL, alpha_BL, alpha_con, alpha_incon, alpha_incon, alpha_con];
betas = [beta_BL, beta_BL, beta_con, beta_incon, beta_incon, beta_con];

slopeValues = nan(1, 6);
positiveAmplitude = all(isfinite(amplitudes)) && all(amplitudes > 0);
validShape = all(isfinite(alphas)) && all(isfinite(betas)) && ...
    all(alphas > 0) && all(betas > 1);
if positiveAmplitude && validShape
    slopeValues = getWeibullHalfMaxSlope(amplitudes, alphas, betas);
end

diagnostics = struct();
diagnostics.amplitudes = amplitudes;
diagnostics.maxSlopeBaselineLeft = slopeValues(1);
diagnostics.maxSlopeBaselineRight = slopeValues(2);
diagnostics.maxSlopeHorizontalLeft = slopeValues(3);
diagnostics.maxSlopeHorizontalRight = slopeValues(4);
diagnostics.maxSlopeVerticalLeft = slopeValues(5);
diagnostics.maxSlopeVerticalRight = slopeValues(6);
diagnostics.maxSlopeValues = slopeValues;
diagnostics.maxSlopeOverall = max(slopeValues, [], 'omitnan');
diagnostics.maxAllowedSlopePctPerContrast = maxSlopePctPerContrast;
diagnostics.slopeConstraintActive = true;
diagnostics.slopeCapActive = isfinite(diagnostics.maxSlopeOverall) && ...
    diagnostics.maxSlopeOverall >= 0.99 .* maxSlopePctPerContrast;
diagnostics.positiveAmplitude = positiveAmplitude;
diagnostics.validShape = validShape;
diagnostics.isValid = positiveAmplitude && validShape && ...
    all(isfinite(slopeValues)) && all(slopeValues <= maxSlopePctPerContrast);
diagnostics.slopeExcess = max(0, slopeValues - maxSlopePctPerContrast);
end