function [fitParams, negLogLikelihood, pseudoR2] = fitNakaRushtonMLE(x, yCells, fixedParams)
    % Assume more informed initial guesses and bounds are set here
    yMeans = cellfun(@mean, cellfun(@double, yCells, 'UniformOutput', false));
    
    % Define bounds and initial guesses based on the data
    [initRmax, initExp, initC50, initBeta] = guessInitialParams(x, yMeans);
    bounds.min = [initRmax * 0.9, 0, initC50*.8, initBeta*.6]; 
    bounds.max = [initRmax * 2, 5, initC50*1.2, initBeta*1.1]; 

    % Optimization options
    options = optimoptions('fmincon', 'Display', 'none', 'Algorithm', 'sqp', 'MaxIterations', 20000, 'TolFun', 1e-8, 'TolX', 1e-8);
    
    if ~isempty(fixedParams)
        fprintf('Constrained fitting with fixed values... ')
        % Adjust the fitting process based on fixedParams
        fixedIndices = ~isnan(fixedParams); % Indices of fixed parameters
        toFitIndices = isnan(fixedParams); % Indices of parameters to fit

        % Adjust initial guesses and bounds for parameters to fit
        initGuessToFit = [initRmax, initExp, initC50, initBeta];
        initGuessToFit=initGuessToFit(toFitIndices);
        boundsToFit.min = bounds.min(toFitIndices);
        boundsToFit.max = bounds.max(toFitIndices);

        if ~isempty(initGuessToFit)
        
            % Adjust negLogLikelihood function for partial fitting
            negLogLikelihoodFunc = @(fitParams) negLogLikelihoodWithFixedParams(fitParams, fixedParams, fixedIndices, toFitIndices, x, yCells);
    
            % Nonlinear constraint to ensure Rmax + Beta <= 100, considering fixedParams
            nonlcon = @(fitParams) deal([], constraintRmaxPlusBeta(fitParams, fixedParams, toFitIndices));
            
            % Perform fitting only for parameters that are NaN in fixedParams
            [fittedParams, negLogLikelihood] = fmincon(negLogLikelihoodFunc, initGuessToFit, [], [], [], [], boundsToFit.min, boundsToFit.max, nonlcon, options);
            
            % Merge fitted parameters with non-NaN entries in fixedParams
            outputParams = fixedParams;
            outputParams(toFitIndices) = fittedParams;
        else
            outputParams=fixedParams;
        end
    else
        fprintf('Free fitting with no fixed values... ')
        % If fixedParams is empty, follow your existing fitting process
        initGuess = [initRmax, initExp, initC50, initBeta];
        nonlcon = @(params) deal([], max(0, params(1) + params(4) - 100)); % Existing constraint
        [outputParams, negLogLikelihood] = fmincon(negLogLikelihoodFunc, initGuess, [], [], [], [], bounds.min, bounds.max, nonlcon, options);
    end

    % Compute pseudo R2 as before
    pseudoR2 = NaN;%calculateNagelkerkeR2(yCells, negLogLikelihood);

    % Assign fitted parameters
    fitParams = struct('rmax', outputParams(1), 'exponent', outputParams(2), 'c50', outputParams(3), 'beta', outputParams(4), 'pseudoR2', pseudoR2);
    
    % Check that fixedParams are obeyed
    checkTmp=vertcat(outputParams)==fixedParams;
    if sum(checkTmp==fixedIndices)==numel(fixedIndices)
        fprintf('Constrains are obeyed\n')
    else
        fprintf('Constrains are not obeyed\n')
    end
end

function [initRmax, initExp, initC50, initBeta] = guessInitialParams(x, yMeans)
    % guess NR: rmax (delta between mean responses)
    initRmax = max(yMeans) - min(yMeans);
    
    % guess NR: exp (intermediate between min and max exponent)
    initExp = mean([3]);
    
    % guess NR: c50 (midpoint between the delta of mean responses)
    initC50 = computeInitC50(yMeans, initRmax, x);
    
    % guess vertical offset: beta (lowest mean response)
    initBeta = min(yMeans);
end

function initC50 = computeInitC50(yMeans, initRmax, x)
    absoluteDelta = abs(yMeans - (initRmax / 2 + min(yMeans)));
    closest_index = find(absoluteDelta == min(absoluteDelta));

    % Adjust for cases with multiple equally close values
    if numel(closest_index) > 2
        closest_index = closest_index([1, 2]); % Take first two if more than two indices are equally close
    elseif numel(closest_index) == 2
        % Already ideal scenario, take them as they are
    else
        % For a single closest index, attempt to add the next index if it doesn't exceed bounds
        if closest_index < numel(yMeans)
            closest_index = [closest_index, closest_index + 1];
        % If the closest index is the last element, consider the previous one as well, if possible
        elseif closest_index > 1
            closest_index = [closest_index - 1, closest_index];
        end
    end

    % Compute the mean of the x values at the adjusted closest indices
    initC50 = mean(x(closest_index));
end

function [c, ceq] = constraintRmaxPlusBeta(fitParams, fullParams, toFitIndices)
    % Update the full parameter set with fitted parameters
    fullParams(toFitIndices) = fitParams;
    
    % Nonlinear inequality constraint ensures Rmax + Beta <= 100
    % Since 'fmincon' minimizes under the condition c <= 0,
    % we rearrange Rmax + Beta <= 100 to Rmax + Beta - 100 <= 0.
    c = fullParams(1) + fullParams(4) - 100;
    
    % No equality constraints
    ceq = [];
end


% Define a sub-function to adjust the negLogLikelihood function to incorporate fixed parameters
function negLL = negLogLikelihoodWithFixedParams(fitParams, fixedParams, fixedIndices, toFitIndices, x, yCells)
    params = fixedParams;
    params(toFitIndices) = fitParams; % Insert fitted parameters into their correct positions
    % Calculate negative log-likelihood using the combined parameter vector
    negLL = -sum(cellfun(@(y, xi) logLikelihoodPerContrast(params, xi, y), yCells, num2cell(x)));
end

function ll = logLikelihoodPerContrast(params, xi, yi)
    % Guard against division by zero or undefined behavior for xi = 0
    xi(xi == 0) = eps; % Small positive value close to zero
    
    % Naka-Rushton function calculation
    Rmax = 100-params(4);
    n = params(2);
    C50 = params(3);
    beta = params(4);
    predicted = (Rmax .* xi.^n) ./ (C50.^n + xi.^n) + beta;
    
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

function ll = logLikelihoodConstantModel(yi, meanResponse)
    sigmaYi = std(yi);
    if sigmaYi == 0, sigmaYi = eps; end % Avoid division by zero in likelihood calculation
    % Compute log-likelihood for a constant model
    ll = -0.5 * numel(yi) * log(2 * pi * sigmaYi^2) + sum(((yi - meanResponse).^2) / (2 * sigmaYi^2));
end

