function maxSlope = getWeibullHalfMaxSlope(amplitude, alpha, beta)
%GETWEIBULLHALFMAXSLOPE Analytic maximum slope of one Weibull half.
%
% The transition is amplitude .* (1 - exp(-(distance ./ alpha).^beta)).
% maxSlope is in percentage points per one percentage point of contrast.

if any(~isfinite(amplitude(:))) || any(~isfinite(alpha(:))) || ...
        any(~isfinite(beta(:)))
    error('getWeibullHalfMaxSlope:NonFiniteInput', ...
        'amplitude, alpha, and beta must be finite.');
end
if any(alpha(:) <= 0)
    error('getWeibullHalfMaxSlope:InvalidAlpha', ...
        'alpha must be positive.');
end
if any(beta(:) <= 1)
    error('getWeibullHalfMaxSlope:InvalidBeta', ...
        'beta must be greater than 1 for a finite analytic maximum slope.');
end

q = (beta - 1) ./ beta;
maxSlope = abs(amplitude) .* (beta ./ alpha) .* q .^ q .* exp(-q);
end