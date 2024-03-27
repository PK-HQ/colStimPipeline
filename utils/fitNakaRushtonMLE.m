function [fitParams, negLogLikelihood] = fitNakaRushtonMLE(x, yCells, initParams)
    % Assume more informed initial guesses and bounds are set here
    yMeans = cellfun(@mean, cellfun(@double, yCells, 'UniformOutput', false));
    
    %% Adjust initial guesses and bounds based on the data (rmax, exp, c50, beta)
    % guess initial params
    [initRmax, initExp, initC50, initBeta] = guessInitialParams(x, yMeans);
    initGuess = [initRmax, initExp, initC50, initBeta];
    bounds.min = [initRmax * 0.9, 0, initC50*.8, initBeta*.6]; % Example adjustment
    bounds.max = [initRmax * 2, 4, initC50*2, initBeta*1.1]; % Example adjustment

    % Define the custom negative log-likelihood function
    negLogLikelihoodFunc = @(params) -sum(cellfun(@(y, xi) logLikelihoodPerContrast(params, xi, y), yCells, num2cell(x)));

    % Optimization options
    options = optimoptions('fmincon', 'Display', 'none', 'Algorithm', 'sqp','MaxIterations',10000);

    % Minimize the negative log-likelihood with adjusted bounds
    [outputParams, negLogLikelihood] = fmincon(negLogLikelihoodFunc, initGuess, [], [], [], [], bounds.min, bounds.max, [], options);

    % Pseudo R2
    pseudoR2 = NaN;%calculateNagelkerkeR2(yCells, negLogLikelihood);

    % Assign fitted parameters
    fitParams = struct('rmax', outputParams(1), 'exponent', outputParams(2), 'c50', outputParams(3), 'beta', outputParams(4), 'pseudoR2', pseudoR2);
end

function [initRmax, initExp, initC50, initBeta] = guessInitialParams(x, yMeans)
    % guess NR: rmax (delta between mean responses)
    initRmax = max(yMeans) - min(yMeans);
    
    % guess NR: exp (intermediate between min and max exponent)
    initExp = mean([0 4]);
    
    % guess NR: c50 (midpoint between the delta of mean responses)
    closest_index = find(abs(yMeans - (initRmax/2 + min(yMeans))) == min(abs(yMeans - (initRmax/2 + min(yMeans)))));
    flanking_indices = closest_index + [0, 1];
    initC50 = mean(x(flanking_indices));
    
    % guess vertical offset: beta (lowest mean response)
    initBeta = min(yMeans);
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

function pseudoR2 = calculateNagelkerkeR2(yCells, negLogLikelihood)
    % Convert negative log-likelihood of the fitted model to likelihood
    Lm = exp(-negLogLikelihood);

    % Estimate the null model likelihood (L0)
    yAll = cell2mat(yCells); % Combine all observations into a single array
    meanResponse = mean(yAll);
    L0 = exp(-sum(cellfun(@(yi) sum(logLikelihoodConstantModel(yi, meanResponse)), yCells)));

    % Calculate Cox and Snell's R^2
    n = numel(yAll); % Total number of observations
    R2_CoxSnell = 1 - (L0 / Lm)^(2 / n);

    % Adjust with Nagelkerke's formula to get Nagelkerke's R^2
    pseudoR2 = R2_CoxSnell / (1 - L0^(2 / n));
end
