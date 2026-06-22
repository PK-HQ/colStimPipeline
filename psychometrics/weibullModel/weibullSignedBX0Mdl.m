function p = weibullSignedBX0Mdl(x, ALeft, alphaLeft, betaLeft, ...
        ARight, alphaRight, betaRight, B, X0)
%WEIBULLSIGNEDBX0MDL Signed choice-vertical Weibull with fitted B and X0.
%
% x is signed contrast. Negative x is horizontal evidence, positive x is
% vertical evidence, and p is percent probability of choosing vertical.
% The two raw branches meet continuously at (X0, B).

x = x(:)';
leftIdx = x < X0;
rightIdx = ~leftIdx;
p = nan(size(x));

if any(leftIdx)
    distanceLeft = abs(x(leftIdx) - X0);
    p(leftIdx) = B - (B - ALeft) .* ...
        (1 - exp(-(distanceLeft ./ alphaLeft) .^ betaLeft));
end

if any(rightIdx)
    distanceRight = abs(x(rightIdx) - X0);
    p(rightIdx) = B + ((100 - ARight) - B) .* ...
        (1 - exp(-(distanceRight ./ alphaRight) .^ betaRight));
end

% Numerical likelihood safety only. Full-precision fitted values remain
% unchanged; the clamp prevents log(0) in binomial likelihood evaluation.
p = min(max(p, eps), 100 - eps);
end
