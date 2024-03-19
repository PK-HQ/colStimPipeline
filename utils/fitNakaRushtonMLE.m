function [fitParams, negLogLikelihood] = fitNakaRushtonMLE(x, yCells, initParams)
    % Convert inputs to ensure compatibility
    x = double(x);
    yCells = cellfun(@double, yCells, 'UniformOutput', false);

    % Define the custom negative log-likelihood function
    negLogLikelihoodFunc = @(params) -sum(cellfun(@(y, xi) logLikelihoodPerContrast(params, xi, y), yCells, num2cell(x)));

    % Optimization options
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
    
    % Minimize the negative log-likelihood with bounds
    initGuess=(initParams.min + initParams.max)/2;
    [outputParams, negLogLikelihood] = fmincon(negLogLikelihoodFunc, initGuess, [], [], [], [], initParams.min, initParams.max, [], options);
    fitParams.rmax=outputParams(1);
    fitParams.exponent=outputParams(2);
    fitParams.c50=outputParams(3);
    fitParams.beta=outputParams(4);
end

function ll = logLikelihoodPerContrast(params, xi, yi)
    % Guard against division by zero or undefined behavior for xi = 0
    xi(xi == 0) = eps; % Small positive value close to zero
    
    % Naka-Rushton function calculation
    Rmax = params(1);
    n = params(2);
    C50 = params(3);
    offset = params(4);
    predicted = (Rmax .* xi.^n) ./ (C50.^n + xi.^n) + offset;
    
    % Standard deviation for yi observations
    sigmaYi = std(yi);
    if sigmaYi == 0, sigmaYi = eps; end % Avoid division by zero in likelihood

    % Compute log-likelihood assuming normally distributed residuals
    ll = -0.5 * sum(((yi - predicted).^2) / sigmaYi^2 + log(2 * pi * sigmaYi^2));
end
