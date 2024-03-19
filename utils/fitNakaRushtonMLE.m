function [fitParams, negLogLikelihood] = fitNakaRushtonMLE(x, y, initParam)
    % Naka-Rushton function
    modelFunc = @(b, x) (b(1) .* x.^b(2)) ./ (b(3).^b(2) + x.^b(2)) + b(4);
    
    % Negative log-likelihood function
    negLogLikelihoodFunc = @(params) -sum(log(normpdf(y, modelFunc(params, x), params(5))));
    
    % Initial parameter guesses and bounds
    initialParams = [initParam.rmax, initParam.n, initParam.c50, initParam.beta, initParam.sigma];
    lb = [0, 0, 0, 0, 0]; % Lower bounds, assuming all parameters are non-negative
    ub = [Inf, Inf, Inf, Inf, Inf]; % Upper bounds, assuming no upper limit
    
    % Minimize the negative log-likelihood
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
    [fitParams, negLogLikelihood] = fmincon(negLogLikelihoodFunc, initialParams, [], [], [], [], lb, ub, [], options);
end
