function [bestParams, fitInfo] = fitBayesianModel(xDataBaseline, yDataBaseline, xDataOpto, yDataOpto, initialGuess, varargin)
% FITBAYESIANMODEL - Fit normalization model parameters to psychometric data
%
% Usage:
%   [bestParams, fitInfo] = fitBayesianModel(xDataBaseline, yDataBaseline, xDataOpto, yDataOpto, initialGuess)
%   [bestParams, fitInfo] = fitBayesianModel(..., 'OptionName', OptionValue)
%
% Inputs:
%   xDataBaseline - Contrast levels for baseline condition
%   yDataBaseline - Baseline condition performance data
%   xDataOpto - Contrast levels for opto condition
%   yDataOpto - Optogenetic condition performance data
%   initialGuess - Initial parameter values [a, b, l, w, g0, n, rmx, e, o]

    % Persistent variables for tracking optimization progress
    persistent bestNLL iterCount;
    
    % Default parameter bounds if not provided
    defaultLB = [0.01, 0.01, 0.01, 0.01, 40, 4, 10, 0.3, 10];
    defaultUB = [1.0, 1.0, 1.0, 1.0, 100, 10, 100, 0.7, 100];
    
    % Parse inputs
    p = inputParser;
    addParameter(p, 'LowerBounds', defaultLB, @isnumeric);
    addParameter(p, 'UpperBounds', defaultUB, @isnumeric);
    addParameter(p, 'ModelOptions', struct('nTrials', 2000, 'useGPU', true, 'smartAllocation', true), @isstruct);
    addParameter(p, 'OptimMethod', 'GlobalThenLocal', @ischar);
    addParameter(p, 'UseParallel', false, @islogical);
    addParameter(p, 'Verbose', true, @islogical);
    addParameter(p, 'MaxIterations', 1000, @isnumeric);
    addParameter(p, 'StrictBoundsEnforcement', true, @islogical); % New parameter for strict bounds enforcement
    parse(p, varargin{:});
    
    lb = p.Results.LowerBounds;
    ub = p.Results.UpperBounds;
    modelOptions = p.Results.ModelOptions;
    optimMethod = p.Results.OptimMethod;
    useParallel = p.Results.UseParallel;
    verbose = p.Results.Verbose;
    maxIter = p.Results.MaxIterations;
    strictBounds = p.Results.StrictBoundsEnforcement;
    
    % Reset persistent variables
    bestNLL = Inf;
    iterCount = 0;
    
    % Helper function to enforce bounds
    function params = enforceParameterBounds(params)
        for i = 1:length(params)
            if params(i) < lb(i) || params(i) > ub(i)
                params(i) = max(min(params(i), ub(i)), lb(i));
            end
        end
    end
    
    % Ensure initial parameters are within bounds
    initialGuess = enforceParameterBounds(initialGuess);
    
    % Initialize output structure
    fitInfo = struct();
    fitInfo.methodUsed = optimMethod;
    fitInfo.initialParams = initialGuess;
    fitInfo.lowerBounds = lb;
    fitInfo.upperBounds = ub;
    fitInfo.modelOptions = modelOptions;
    
    % Convert percentage correct to trial counts (assuming 20 trials per condition)
    [sumBaselineTrials, successBaselineTrials] = convertToTrialCounts(xDataBaseline, yDataBaseline);
    [sumOptoTrials, successOptoTrials] = convertToTrialCounts(xDataOpto, yDataOpto);
    
    % Define objective function (negative log-likelihood)
    function nll = objectiveNLL(params)
        % Ensure parameters are within bounds - strictly enforce if enabled
        if strictBounds
            params = enforceParameterBounds(params);
        end
        
        % Get model predictions
        try
            % Add an extra bounds check before model evaluation
            boundedParams = enforceParameterBounds(params);
            
            predBaseline = normMdlSimOptimized(xDataBaseline, boundedParams, 'baseline', modelOptions);
            predOpto = normMdlSimOptimized(xDataOpto, boundedParams, 'opto', modelOptions);
            
            % Ensure predictions are within valid range for probability
            predBaseline = max(0.001, min(0.999, predBaseline/100));
            predOpto = max(0.001, min(0.999, predOpto/100));
            
            % Calculate negative log-likelihood
            nllBaseline = -sum(successBaselineTrials .* log(predBaseline) + ...
                         (sumBaselineTrials - successBaselineTrials) .* log(1 - predBaseline));
            nllOpto = -sum(successOptoTrials .* log(predOpto) + ...
                     (sumOptoTrials - successOptoTrials) .* log(1 - predOpto));
            
            nll = nllBaseline + nllOpto;
            
            % Check for invalid values and handle them
            if ~isfinite(nll)
                nll = 1e10; % Return a large value for invalid parameters
            end
            
            % Display progress if verbose
            if nll < bestNLL
                bestNLL = nll;
                if verbose
                    fprintf('Iteration %d: NLL = %.4f, [a=%.2f b=%.2f l=%.2f w=%.2f g0=%.1f n=%.1f rmx=%.1f e=%.2f o=%.1f]\n', ...
                        iterCount, nll, boundedParams(1), boundedParams(2), boundedParams(3), boundedParams(4), boundedParams(5), ...
                        boundedParams(6), boundedParams(7), boundedParams(8), boundedParams(9));
                end
            end
            iterCount = iterCount + 1;
            
        catch ME
            % If error occurs, return large value
            if verbose
                warning('Error in objective function: %s', ME.message);
            end
            nll = 1e10;
        end
    end
    
    % Create wrapper functions for optimizers that enforce bounds
    function [x, fval, exitflag, output] = fminconWrapper(fun, x0, A, b, Aeq, beq, lb, ub, nonlcon, options)
        % Wrapper for fmincon with bounds enforcement
        function y = boundedFun(x)
            if strictBounds
                x = enforceParameterBounds(x);
            end
            y = fun(x);
        end
        
        [x, fval, exitflag, output] = fmincon(@boundedFun, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);
        x = enforceParameterBounds(x); % Ensure final result is within bounds
    end
    
    % Start timing
    if verbose
        fprintf('Starting parameter fitting using %s method...\n', optimMethod);
        tic;
    end
    
    % Set display option based on verbose flag
    if verbose
        display_opt = 'iter';
    else
        display_opt = 'off';
    end
    
    % Perform optimization based on selected method
    switch optimMethod
        case 'GridSearch'
            if verbose
                fprintf('Running grid search for key parameters...\n');
            end
            % Define grid search values for key parameters
            g0_values = [50, 75, 100];
            n_values = [4, 6, 8];
            e_values = [0.4, 0.5, 0.6];
            
            % Initialize tracking
            bestNLL = Inf;
            bestGridParams = initialGuess;
            
            % Loop through parameter combinations
            totalIters = length(g0_values) * length(n_values) * length(e_values);
            counter = 0;
            
            for g0_idx = 1:length(g0_values)
                for n_idx = 1:length(n_values)
                    for e_idx = 1:length(e_values)
                        counter = counter + 1;
                        if verbose
                            fprintf('Grid search: %d/%d combinations\n', counter, totalIters);
                        end
                        
                        % Create test parameters
                        testParams = initialGuess;
                        testParams(5) = g0_values(g0_idx);  % g0
                        testParams(6) = n_values(n_idx);    % n
                        testParams(8) = e_values(e_idx);    % e
                        
                        % Ensure test parameters are within bounds
                        testParams = enforceParameterBounds(testParams);
                        
                        % Calculate NLL
                        nll = objectiveNLL(testParams);
                        
                        % Update best parameters if better
                        if nll < bestNLL
                            bestNLL = nll;
                            bestGridParams = testParams;
                            if verbose
                                fprintf('  New best: NLL = %.4f\n', bestNLL);
                            end
                        end
                    end
                end
            end
            
            bestParams = bestGridParams;
            fitInfo.gridSearchNLL = bestNLL;
            
        case 'Global'
            % Global optimization using Genetic Algorithm
            if verbose
                fprintf('Running global optimization with GA...\n');
            end
            
            % Setup GA options
            gaOptions = optimoptions('ga', 'Display', display_opt, ...
                'PopulationSize', 50, ...
                'MaxGenerations', ceil(maxIter/50), ...
                'UseParallel', useParallel, ...
                'MaxStallGenerations', 20, ...
                'FunctionTolerance', 1e-6);
            
            % Create constraint function to enforce bounds
            nonlcon = [];
            
            % Run GA with strictly enforced bounds
            [bestParams, nll] = ga(@objectiveNLL, length(initialGuess), [], [], [], [], lb, ub, nonlcon, gaOptions);
            
            % Double-check bounds compliance
            bestParams = enforceParameterBounds(bestParams);
            
            fitInfo.finalNLL = nll;
            
        case 'Local'
            % Local optimization using fmincon with interior-point
            if verbose
                fprintf('Running local optimization with fmincon...\n');
            end
            
            % Setup fmincon options
            options = optimoptions('fmincon', 'Algorithm', 'interior-point', ...
                               'Display', display_opt, ...
                               'MaxFunctionEvaluations', maxIter*2, ...
                               'MaxIterations', maxIter, ...
                               'OptimalityTolerance', 1e-6);
            
            % Run fmincon with strict bounds enforcement via our wrapper
            [bestParams, nll] = fminconWrapper(@objectiveNLL, initialGuess, [], [], [], [], lb, ub, [], options);
            
            % Ensure parameters are within bounds
            bestParams = enforceParameterBounds(bestParams);
            
            fitInfo.finalNLL = nll;
            
        case 'GlobalThenLocal'
            % First run grid search for a good starting point
            if verbose
                fprintf('Running initial grid search...\n');
            end
            
            % Define grid search values for key parameters
            g0_values = [50, 75, 100];
            n_values = [4, 6, 8];
            e_values = [0.4, 0.5, 0.6];
            
            % Initialize tracking
            bestNLL = Inf;
            bestGridParams = initialGuess;
            
            % Loop through parameter combinations
            totalIters = length(g0_values) * length(n_values) * length(e_values);
            counter = 0;
            
            for g0_idx = 1:length(g0_values)
                for n_idx = 1:length(n_values)
                    for e_idx = 1:length(e_values)
                        counter = counter + 1;
                        if verbose && mod(counter, 5) == 0
                            fprintf('Grid search: %d/%d combinations\n', counter, totalIters);
                        end
                        
                        % Create test parameters
                        testParams = initialGuess;
                        testParams(5) = g0_values(g0_idx);  % g0
                        testParams(6) = n_values(n_idx);    % n
                        testParams(8) = e_values(e_idx);    % e
                        
                        % Ensure test parameters are within bounds
                        testParams = enforceParameterBounds(testParams);
                        
                        % Calculate NLL
                        nll = objectiveNLL(testParams);
                        
                        % Update best parameters if better
                        if nll < bestNLL
                            bestNLL = nll;
                            bestGridParams = testParams;
                            if verbose
                                fprintf('  New best: NLL = %.4f\n', bestNLL);
                            end
                        end
                    end
                end
            end
            
            % Set grid search results as initial point for local search
            initialGuess = bestGridParams;
            fitInfo.gridSearchNLL = bestNLL;
            fitInfo.gridSearchParams = bestGridParams;
            
            % Run local optimization with fmincon
            if verbose
                fprintf('\nRunning local optimization from grid search result...\n');
            end
            
            % Setup fmincon options
            options = optimoptions('fmincon', 'Algorithm', 'interior-point', ...
                               'Display', display_opt, ...
                               'MaxFunctionEvaluations', maxIter*2, ...
                               'MaxIterations', maxIter, ...
                               'OptimalityTolerance', 1e-6);
            
            % Run fmincon with strict bounds enforcement via our wrapper
            [bestParams, nll] = fminconWrapper(@objectiveNLL, initialGuess, [], [], [], [], lb, ub, [], options);
            
            % Ensure parameters are within bounds
            bestParams = enforceParameterBounds(bestParams);
            
            fitInfo.finalNLL = nll;

        otherwise
            error('Unknown optimization method: %s', optimMethod);
    end
    
    % Final verification that parameters are within bounds
    bestParams = enforceParameterBounds(bestParams);
    
    % Calculate model statistics
    totalTrials = sum(sumBaselineTrials) + sum(sumOptoTrials);
    numParams = length(bestParams);
    
    % Calculate AIC and BIC
    fitInfo.finalNLL = objectiveNLL(bestParams); % Recalculate with enforced bounds
    aic = 2 * numParams + 2 * fitInfo.finalNLL;
    if totalTrials > numParams + 1
        aicc = aic + (2 * numParams * (numParams + 1)) / (totalTrials - numParams - 1);
    else
        aicc = Inf;
    end
    bic = log(totalTrials) * numParams + 2 * fitInfo.finalNLL;
    
    % Calculate final predictions
    predBaseline = normMdlSimOptimized(xDataBaseline, bestParams, 'baseline', modelOptions);
    predOpto = normMdlSimOptimized(xDataOpto, bestParams, 'opto', modelOptions);
    
    % Calculate R-squared
    ssTotal_baseline = sum((yDataBaseline - mean(yDataBaseline)).^2);
    ssResid_baseline = sum((yDataBaseline - predBaseline).^2);
    rSquared_baseline = 1 - (ssResid_baseline / ssTotal_baseline);
    
    ssTotal_opto = sum((yDataOpto - mean(yDataOpto)).^2);
    ssResid_opto = sum((yDataOpto - predOpto).^2);
    rSquared_opto = 1 - (ssResid_opto / ssTotal_opto);
    
    % Store fitting info
    fitInfo.aic = aic;
    fitInfo.aicc = aicc;
    fitInfo.bic = bic;
    fitInfo.rSquared_baseline = rSquared_baseline;
    fitInfo.rSquared_opto = rSquared_opto;
    fitInfo.predictedBaseline = predBaseline;
    fitInfo.predictedOpto = predOpto;
    
    % End timing
    if verbose
        elapsedTime = toc;
        fprintf('\nOptimization completed in %.2f seconds\n', elapsedTime);
        fprintf('Best params: [a=%.3f b=%.3f l=%.3f w=%.3f g0=%.1f n=%.1f rmx=%.1f e=%.3f o=%.1f]\n', ...
            bestParams(1), bestParams(2), bestParams(3), bestParams(4), bestParams(5), ...
            bestParams(6), bestParams(7), bestParams(8), bestParams(9));
        fprintf('AIC: %.2f, AICc: %.2f, BIC: %.2f\n', aic, aicc, bic);
        fprintf('R-squared (baseline): %.4f, R-squared (opto): %.4f\n', ...
            rSquared_baseline, rSquared_opto);
        
        % Additional diagnostics for bounds
        for i = 1:length(bestParams)
            if bestParams(i) <= lb(i)*1.01 || bestParams(i) >= ub(i)*0.99
                fprintf('Warning: Parameter %d (%.4f) is near boundary [%.4f, %.4f]\n', ...
                    i, bestParams(i), lb(i), ub(i));
            end
        end
    end
end

function [sumTrials, successTrials] = convertToTrialCounts(x, y)
    % Convert percentage correct to trial counts
    % Assumes 20 trials per condition by default
    
    % Make sure x and y are row vectors
    x = x(:)';
    y = y(:)';
    
    % Create trial counts (assuming 20 trials per condition)
    sumTrials = 20 * ones(size(x));
    
    % Convert percentage to counts
    successTrials = (y / 100) .* sumTrials;
end