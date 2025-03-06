function [bestParams, fitInfo] = fitBayesianModel(xDataBaseline, yDataBaseline, xDataOpto, yDataOpto, initialGuess, varargin)
% FITBAYESIANMODEL - Fit normalization model parameters to psychometric data
%
% Usage:
%   [bestParams, fitInfo] = fitBayesianModel(xData, yDataBaseline, yDataOpto, initialGuess)
%   [bestParams, fitInfo] = fitBayesianModel(..., 'OptionName', OptionValue)
%
% Inputs:
%   xData - Contrast levels
%   yDataBaseline - Baseline condition performance data
%   yDataOpto - Optogenetic condition performance data
%   initialGuess - Initial parameter values [a, b, l, w, g0, n, rmx, e, o]

    % Persistent variables for tracking optimization progress
    persistent bestNLL iterCount;
    
    % Default parameter bounds if not provided
    defaultLB = [0.01, 0.01, 0.01, 0.01, 25, 2, 10, 0.2, 10];
    defaultUB = [1.0, 1.0, 1.0, 1.0, 150, 10, 100, 0.8, 100];
    
    % Parse inputs
    p = inputParser;
    addParameter(p, 'LowerBounds', defaultLB, @isnumeric);
    addParameter(p, 'UpperBounds', defaultUB, @isnumeric);
    addParameter(p, 'ModelOptions', struct('nTrials', 2000, 'useGPU', true, 'smartAllocation', true), @isstruct);
    addParameter(p, 'OptimMethod', 'GlobalThenLocal', @ischar);
    addParameter(p, 'UseParallel', false, @islogical);
    addParameter(p, 'Verbose', true, @islogical);
    addParameter(p, 'MaxIterations', 1000, @isnumeric);
    parse(p, varargin{:});
    
    lb = p.Results.LowerBounds;
    ub = p.Results.UpperBounds;
    modelOptions = p.Results.ModelOptions;
    optimMethod = p.Results.OptimMethod;
    useParallel = p.Results.UseParallel;
    verbose = p.Results.Verbose;
    maxIter = p.Results.MaxIterations;
    
    % Reset persistent variables
    bestNLL = Inf;
    iterCount = 0;
    
    % Ensure initial parameters are within bounds
    initialGuess = max(min(initialGuess, ub), lb);
    
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
        % Get model predictions
        try
            predBaseline = normMdlSimOptimized(xDataBaseline, params, 'baseline', modelOptions);
            predOpto = normMdlSimOptimized(xDataOpto, params, 'opto', modelOptions);
            
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
                        iterCount, nll, params(1), params(2), params(3), params(4), params(5), ...
                        params(6), params(7), params(8), params(9));
                end
            end
            iterCount = iterCount + 1;
            
        catch ME
            % If error occurs, return large value
            %warning('Error in objective function: %s', ME.message);
            nll = 1e10;
        end
    end
    
    % Start timing
    if verbose
        fprintf('Starting parameter fitting using %s method...\n', optimMethod);
        tic;
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
            gaOptions = optimoptions('ga', 'Display', 'iter', ...
                'PopulationSize', 50, ...
                'MaxGenerations', ceil(maxIter/50), ...
                'UseParallel', useParallel, ...
                'MaxStallGenerations', 20, ...
                'FunctionTolerance', 1e-6);
            
            if ~verbose
                gaOptions.Display = 'off';
            end
            
            % Run GA
            [bestParams, nll] = ga(@objectiveNLL, length(initialGuess), [], [], [], [], lb, ub, [], gaOptions);
            fitInfo.finalNLL = nll;
            
        case 'Local'
            % Local optimization using fminsearch with bounds
            if verbose
                fprintf('Running local optimization with fminsearchbnd...\n');
            end
            
            % Setup fminsearch options
            options = optimset('fminsearch');
            options = optimset(options, 'MaxIter', maxIter);
            options = optimset(options, 'MaxFunEvals', maxIter * 2);
            options = optimset(options, 'TolX', 1e-6);
            options = optimset(options, 'TolFun', 1e-6);
            
            if verbose
                options = optimset(options, 'Display', 'iter');
            else
                options = optimset(options, 'Display', 'off');
            end
            
            % Run bounded fminsearch
            [bestParams, nll] = fminsearchbnd(@objectiveNLL, initialGuess, lb, ub, options);
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
            
            % Run local optimization with fminsearchbnd
            if verbose
                fprintf('\nRunning local optimization from grid search result...\n');
            end
            
            % Setup fminsearch options
            options = optimset('fminsearch');
            options = optimset(options, 'MaxIter', maxIter);
            options = optimset(options, 'MaxFunEvals', maxIter * 2);
            options = optimset(options, 'TolX', 1e-6);
            options = optimset(options, 'TolFun', 1e-6);
            
            if verbose
                options = optimset(options, 'Display', 'iter');
            else
                options = optimset(options, 'Display', 'off');
            end
            
            % Run bounded fminsearch
            [bestParams, nll] = fminsearchbnd(@objectiveNLL, initialGuess, lb, ub, options);
            fitInfo.finalNLL = nll;
            
        otherwise
            error('Unknown optimization method: %s', optimMethod);
    end
    
    % Calculate model statistics
    totalTrials = sum(sumBaselineTrials) + sum(sumOptoTrials);
    numParams = length(bestParams);
    
    % Calculate AIC and BIC
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
    end
end

% --- Subfunctions ---

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