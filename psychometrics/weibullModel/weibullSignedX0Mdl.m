function p = weibullSignedX0Mdl(x, A, alpha, beta, X0)
%WEIBULLSIGNEDX0MDL Signed mirrored Weibull choice-vertical prediction.
%
% x is signed Gabor contrast. Negative x is horizontal evidence, positive x
% is vertical evidence, and p is percent probability of choosing vertical.
% The decision boundary is fixed by X0, so p(X0) == 50.

z = (x - X0) ./ alpha;

p = 100 .* ( ...
    0.5 + ...
    (0.5 - A) .* sign(z) .* ...
    (1 - exp(-abs(z) .^ beta)) ...
    );

% Floating-point safety only. Parameters are bounded by the fitter.
p = min(max(p, eps), 100 - eps);
end
