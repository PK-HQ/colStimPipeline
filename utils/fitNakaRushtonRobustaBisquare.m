function fitParams = fitNakaRushtonRobustaBisquare(x, y, initParam)
    % Remove NaNs
    x = x(~isnan(y));
    y = y(~isnan(y));

    % Define the Naka-Rushton function
    modelFunc = @(b, x) (b(1) .* (x.^b(2)) ./ (max(b(3).^b(2) + x.^b(2), 1e-8))) + b(4);

    % Adjust the tuning constant to make the method more sensitive to outliers
    tuningConst = 4.685; % Lowering this value can increase sensitivity to outliers

    % Initial parameter guess and bounds
    initialParams = (initParam.min + initParam.max) / 2;
    lb = [initParam.min(1), initParam.min(2), initParam.min(3), initParam.min(4)];
    ub = [initParam.max(1), initParam.max(2), initParam.max(3), initParam.max(4)];

    % Optimization options
    opts = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off');

    % Nonlinear constraint to ensure rmax + beta <= 100
    nonlincon = @(b) deal([], max(0, b(1) + b(4) - 100));

    % Perform an initial fit to estimate residuals
    [initialFitParams, ~, ~, exitflag, output] = fmincon(@(b) sum((y - modelFunc(b, x)).^2), initialParams, [], [], [], [], lb, ub, nonlincon, opts);

    % Calculate residuals from the initial fit
    residuals = y - modelFunc(initialFitParams, x);
    madResiduals = mad(residuals, 1); % Median Absolute Deviation of residuals

    % Calculate bisquare weights for each residual
    weights = ((abs(residuals) / (tuningConst * madResiduals)) < 1) .* (1 - (residuals / (tuningConst * madResiduals)).^2).^2;

    % Robust fitting: Minimize the weighted sum of squared errors
    robustLossFunc = @(b) sum(weights .* (y - modelFunc(b, x)).^2);
    [fittedParams, ~, ~, exitflag, output] = fmincon(robustLossFunc, initialFitParams, [], [], [], [], lb, ub, nonlincon, opts);

    % Evaluate the model with the fitted parameters
    fittedY = modelFunc(fittedParams, x);

    % Calculate weighted residuals for R-squared calculation
    finalResiduals = y - fittedY;
    finalWeights = ((abs(finalResiduals) / (tuningConst * mad(finalResiduals, 1))) < 1) .* (1 - (finalResiduals / (tuningConst * mad(finalResiduals, 1))).^2).^2;
    weightedSSRes = sum(finalWeights .* finalResiduals.^2);
    weightedSSTot = sum(finalWeights .* (y - mean(y)).^2);

    % Calculate pseudo-R-squared based on weighted residuals
    pseudoR2 = 1 - weightedSSRes / weightedSSTot;

    % Construct output with best fit parameters and goodness-of-fit
    fitParams = struct('rmax', fittedParams(1), 'exponent', fittedParams(2), 'c50', fittedParams(3), 'beta', fittedParams(4), 'R2', pseudoR2, 'ExitFlag', exitflag, 'Output', output);
end
