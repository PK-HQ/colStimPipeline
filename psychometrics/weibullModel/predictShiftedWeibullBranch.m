function y = predictShiftedWeibullBranch(contrastMagnitude, A, B, alpha, beta, X0)
%PREDICTSHIFTEDWEIBULLBRANCH Positive-axis shifted Weibull branch.
% Inputs and output are in percent-correct units. The branch is flat at B
% through X0, then rises toward 100-A using only its own shape parameters.

values = [A, B, alpha, beta, X0];
if any(~isfinite(values)) || alpha <= 0 || beta <= 1 || A >= 100 - B
    y = nan(size(contrastMagnitude));
    return;
end

c = max(0, contrastMagnitude);
distance = max(c - X0, 0);
y = B + ((100 - A) - B) .* ...
    (1 - exp(-(distance ./ alpha) .^ beta));

if ~isreal(y) || any(~isfinite(y(:))) || any(y(:) < 0) || any(y(:) > 100)
    y(:) = NaN;
end
end