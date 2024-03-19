function fitParams = fitNakaRushtonRobusta(x, y, initParam)
    x = x(~isnan(y));
    y = y(~isnan(y));

    % Define Naka-Rushton function
    modelFunc = @(b, x) (b(1) .* (x.^b(2)) ./ (max(b(3).^b(2) + x.^b(2), 1e-8))) + b(4);

    % Custom loss function for robust fitting
    lossFunc = @(b) sum((y - modelFunc(b, x)).^2); % Simple SSE for demonstration
    
    % Nonlinear constraint to ensure rmax + beta <= 100
    nonlincon = @(b) deal([], max(0, b(1) + b(4) - 100)); % Second output deals with inequality constraints

    % Set bounds for the parameters
    lb = [initParam.min(1), initParam.min(2), initParam.min(3), initParam.min(4)];
    ub = [initParam.max(1), initParam.max(2), initParam.max(3), initParam.max(4)];
    
    % Initial guess for the parameters (midpoint of the bounds)
    initialParams = (lb + ub) / 2;
    
    % Optimization options
    opts = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off');

    % Perform the fitting with nonlinear constraint
    [fittedParams,~,exitflag,output] = fmincon(lossFunc, initialParams, [], [], [], [], lb, ub, nonlincon, opts);
    
    % Evaluate the model with the fitted parameters
    fittedY = modelFunc(fittedParams, x);

    % Calculate R-squared as a goodness-of-fit measure
    ssRes = sum((y - fittedY).^2);
    ssTot = sum((y - mean(y)).^2);
    rSquared = 1 - (ssRes / ssTot);

    % Construct output with best fit parameters and goodness-of-fit
    fitParams = struct('rmax', fittedParams(1), 'exponent', fittedParams(2), 'c50', fittedParams(3), 'beta', fittedParams(4), 'R2', rSquared, 'ExitFlag', exitflag, 'Output', output);
end
