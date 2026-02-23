function [mdl, p] = fitSaturatingCurve(x, y, linecolor, robustFlag)
% One-phase rectangular hyperbola: y = A + (B*x)/(K+x)
% Returns mdl (fit on scaled x) and p = [A B K] back-scaled.

if nargin < 3, linecolor = 'k';      end
if nargin < 4, robustFlag = false;   end

% --- SCALE X to tame the Jacobian ---------------------------------------
xScale = max(x);              % put data into roughly [0-1]
xs     = x ./ xScale;

% --- SMART INITIALS ------------------------------------------------------
A0 = min(y);
B0 = max(y) - A0;
[~,iHalf] = min(abs(y - (A0 + 0.5*B0)));        % closest to half-max
K0 = xs(iHalf);
b0 = [A0, B0, max(K0, 0.05)];                   % avoid K0 = 0

% --- MODEL ---------------------------------------------------------------
fh = @(b,x) b(1) + (b(2).*x) ./ (b(3) + x);     % hyperbola

% --- OPTIONS: more iterations + robust weights --------------------------
opts             = statset('nlinfit');
opts.MaxIter     = 5000;                        % raise the ceiling:contentReference[oaicite:4]{index=4}
if robustFlag
    opts.Robust  = 'on';                        % bisquare weighting:contentReference[oaicite:5]{index=5}
end

tbl = table(xs(:), y(:), 'VariableNames', {'x','y'});
mdl = fitnlm(tbl, fh, b0, 'Options', opts);     % fit on scaled x

% --- RESCALE COEFFICIENTS (read-only workaround) ------------------------
c          = mdl.Coefficients.Estimate;
p          = [c(1), c(2), c(3)*xScale];         % A, B, K in original units

% --- PLOT ---------------------------------------------------------------
xx = linspace(min(x), max(x), 200);
yy = fh([c(1:2).', c(3)], xx./xScale);          % note transpose!
hold on
plot(xx, yy, '-', 'Color', linecolor, 'LineWidth', 2);
end
