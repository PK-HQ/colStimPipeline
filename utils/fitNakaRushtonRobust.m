function fitParams = fitNakaRushtonRobust(x, y, initParam)
    x=x(~isnan(y));
    y=y(~isnan(y));
    % Define Naka-Rushton function
    modelFunc = @(b, x) (b(1) .* (x.^b(2)) ./ (max(b(3).^b(2) + x.^b(2),1e-8))) + b(4);

    % Custom loss function for robust fitting
    lossFunc = @(b) sum((y - modelFunc(b, x)).^2); % Simple SSE for demonstration
    
    % Set bounds for the parameters
    lb = [initParam.min(1), initParam.min(2), initParam.min(3), initParam.min(4)];
    ub = [initParam.max(1), initParam.max(2), initParam.max(3), initParam.max(4)];
    
    % Initial guess for the parameters (midpoint of the bounds)
    initialParams = (lb + ub) / 2;
    
    % Optimization options
    opts = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off');

    % Perform the fitting
    [fittedParams,~,exitflag,output] = fmincon(lossFunc, initialParams, [], [], [], [], lb, ub, [], opts);
    
    % Compute log-likelihood as an approximation based on the loss function
    logLikelihood = -0.5 * lossFunc(fittedParams); % Simplified log-likelihood
    
    % Construct output with best fit parameters
    fitParams = struct('rmax', fittedParams(1), 'exponent', fittedParams(2), 'c50', fittedParams(3), 'beta', fittedParams(4), 'LogLikelihood', logLikelihood, 'ExitFlag', exitflag, 'Output', output);
end
