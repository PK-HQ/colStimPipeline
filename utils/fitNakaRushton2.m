function fitParams = fitNakaRushton2(x, y, initParam)
    % Define Naka-Rushton function with additional safeguards
    epsilon = 1e-8; % Small constant to prevent division by zero
    modelFunc = @(b, x) (b(1) .* (x.^b(2))) ./ (b(3).^b(2) + x.^b(2) + epsilon) + b(4);

    % Adjust initial guess if necessary
    initialGuess = [10, 2, 45, 10]; % Example: Choose based on your data's expected behavior

    % Set options for lsqcurvefit
    options = optimoptions('lsqcurvefit', 'Display', 'iter', 'Algorithm', 'trust-region-reflective');

    % Lower and upper bounds
    lb = initParam.min;
    ub = initParam.max;

    % Check model function at initial guess for debugging
    initialModelOutput = modelFunc(initialGuess, x);
    if any(isnan(initialModelOutput) | isinf(initialModelOutput))
        error('Initial model output contains NaN or Inf values. Adjust the initial guess or model function.');
    end

    % Fit model using lsqcurvefit
    [bestParams,~,residual,~,~,~,J] = lsqcurvefit(modelFunc, initialGuess, x(~isnan(y)), y(~isnan(y)), lb, ub, options);

    % Simplified LogLikelihood calculation
    sigma2 = mean(residual.^2); % Estimate of the variance
    logLikelihood = -0.5 * numel(y) * log(2 * pi * sigma2) - 0.5 * sum((residual.^2) / sigma2);

    % Construct output with best fit parameters
    fitParams = struct('rmax', bestParams(1), 'exponent', bestParams(2), 'c50', bestParams(3), 'beta', bestParams(4), 'LogLikelihood', logLikelihood);
end
