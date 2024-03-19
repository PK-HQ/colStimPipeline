function fitParams = fitNakaRushton(x, y, initParam)    
    x=x(~isnan(y));
    y=y(~isnan(y));
    % Define Naka-Rushton function with safeguards against Inf/NaN values: response = rmax .* (x^n) / (c50^n + x^n) + B
    epsilon = 1e-8; % Small value to avoid division by zero
    modelFunc = @(b, x) (b(1) .* (x.^b(2)) ./ (b(3).^b(2) + x.^b(2) + epsilon)) + b(4);
    
    % Initialize best fit tracking variables
    bestModel = [];
    bestLogLikelihood = -Inf;
    bestParams = [];
    
    % Parameter grid setup
    numSteps = 3; % Example grid density
    rmaxGrid = linspace(initParam.min(1), initParam.max(1), numSteps);
    exponentGrid = linspace(initParam.min(2), initParam.max(2), numSteps);
    c50Grid = linspace(initParam.min(3), initParam.max(3), numSteps);
    betaGrid = linspace(initParam.min(4), initParam.max(4), numSteps);
    
    % Attempt fits over the parameter grid
    for rmax = rmaxGrid
        for exponent = exponentGrid
            for c50 = c50Grid
                for beta = betaGrid
                    initialParams = [rmax, exponent, c50, beta];

                    try
                        % Enforce limits
                        if rmax + beta > 100 || rmax >100
                            continue
                        end
                        if exponent < initParam.min(2) || exponent > initParam.max(2) || exponent < 1 || exponent >10
                            continue
                        end
                        % Prepare data table for fitnlm
                        tbl = table(x(:), y(:), 'VariableNames', {'x', 'y'});
                        opts = statset('nlinfit');
                        opts.RobustWgtFun = 'bisquare'; % Enable robust fitting
                        opts.Display = 'off';
                        % Attempt to fit model with current parameters
                        mdl = fitnlm(tbl, @(b,x)modelFunc(b,x), initialParams, 'Options', opts);
                        
                        % Simplified log-likelihood calculation
                        logLikelihood = -0.5 * mdl.SSE; % Note: Real log-likelihood computation might differ
                        
                        % Update best model if this is the best fit so far
                        if logLikelihood > bestLogLikelihood
                            bestModel = mdl;
                            bestLogLikelihood = logLikelihood;
                            bestParams = mdl.Coefficients.Estimate;
                        end
                    catch ME
                        %fprintf('Fit attempt failed with parameters: [%f, %f, %f, %f]\n', initialParams);
                        %fprintf('Reason: %s\n', ME.message);
                    end
                end
            end
        end
    end
    
    % Verify successful fit
    if isempty(bestModel)
        error('No successful fit was achieved.');
    else
        % Construct output with best fit parameters
        fitParams = struct('rmax', bestParams(1), 'exponent', bestParams(2), 'c50', bestParams(3), 'beta', bestParams(4), 'LogLikelihood', bestLogLikelihood);
    end
    %{
    paramsText = sprintf(...
     ['{\\bfParameters}\n',...
     'R_{max}: %.1f\n', ...
     'n: %.1f',... 
     'C_{50}: %.1f\n', ...
      '\\beta: %.1f\n',...
      'L: %.1f'], ...
      fitParams.rmax, fitParams.exponent, fitParams.c50, fitParams.beta, fitParams.LogLikelihood);
      text(.3,-.1,paramsText,'Units','normalized', 'VerticalAlignment','top', 'HorizontalAlignment','left')
    %}
end