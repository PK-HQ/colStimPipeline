function fitParams = fitNakaRushtonRobustIRLS(x, y, initParam)
    % Remove NaNs
    x = x(~isnan(y));
    y = y(~isnan(y));

    % Define the Naka-Rushton function
    modelFunc = @(b, x) b(1) .* (x.^b(2)) ./ (max(b(3).^b(2) + x.^b(2), 1e-8)) + b(4);

    % Initial parameter guess and bounds
    initialParams = (initParam.min + initParam.max) / 2;
    lb = [initParam.min(1), initParam.min(2), initParam.min(3), initParam.min(4)];
    ub = [initParam.max(1), initParam.max(2), initParam.max(3), initParam.max(4)];

    % Optimization options
    opts = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off');

    % Nonlinear constraint to ensure rmax + beta <= 100
    nonlincon = @(b) deal([], max(0, b(1) + b(4) - 100));

    % Iterative fitting process
    numIterations = 5; % Number of iterations for reweighting
    weights = ones(size(y)); % Initial weights
    for iter = 1:numIterations
        % Weighted loss function
        weightedLossFunc = @(b) sum(weights .* (y - modelFunc(b, x)).^2);
        
        % Perform the fitting
        [fittedParams, ~, ~, exitflag, output] = fmincon(weightedLossFunc, initialParams, [], [], [], [], lb, ub, nonlincon, opts);

        % Update initial parameters for the next iteration
        initialParams = fittedParams;

        % Calculate new residuals
        residuals = y - modelFunc(fittedParams, x);

        % Update weights - Huber-like weighting function
        k = 1.5 * mad(residuals, 1); % Scale factor (MAD = median absolute deviation)
        weights = (abs(residuals) < k) + (abs(residuals) >= k) .* k ./ abs(residuals);
    end

    % Evaluate the model with the final fitted parameters
    fittedY = modelFunc(fittedParams, x);

    % Calculate goodness-of-fit (R^2)
    ssRes = sum((y - fittedY).^2);
    ssTot = sum((y - mean(y)).^2);
    rSquared = 1 - (ssRes / ssTot);

    % Construct output with best fit parameters and goodness-of-fit
    fitParams = struct('rmax', fittedParams(1), 'exponent', fittedParams(2), 'c50', fittedParams(3), 'beta', fittedParams(4), 'R2', rSquared, 'ExitFlag', exitflag, 'Output', output);
end
