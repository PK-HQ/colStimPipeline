function [bestParams, fitInfo] = fitBayesianModel(xDataBaseline, yDataBaseline, xDataOpto, yDataOpto, initialGuess, varargin)
% FITBAYESIANMODEL - Fit normalization model parameters to psychometric data
%
% IMPROVEMENTS SUMMARY:
% 1. Enhanced global optimization with Particle Swarm Optimization (PSO)
% 2. Multi-start local optimization to avoid local minima
% 3. Adaptive parameter space exploration with logarithmic sampling
% 4. Improved convergence criteria and numerical stability
% 5. Better parameter tracking during optimization
% 6. More robust error handling and recovery strategies
%
% Usage:
%   [bestParams, fitInfo] = fitBayesianModel(xDataBaseline, yDataBaseline, xDataOpto, yDataOpto, initialGuess)
%   [bestParams, fitInfo] = fitBayesianModel(..., 'OptionName', OptionValue)
%
% Inputs:
%   xDataBaseline - Contrast levels for baseline condition
%   yDataBaseline - Baseline condition performance data (% correct)
%   xDataOpto     - Contrast levels for opto condition
%   yDataOpto     - Optogenetic condition performance data (% correct)
%   initialGuess  - Initial parameter values [a, b, l, w, g0, n, rmx, e, o]
%
% Optional Name-Value Parameters:
%   'LowerBounds'    - Lower parameter bounds (default: conservative values)
%   'UpperBounds'    - Upper parameter bounds (default: conservative values)
%   'ModelOptions'   - Options for normMdlSim* (struct)
%   'OptimMethod'    - Optimization method ('PSO_Multi', 'GlobalThenLocal', etc.)
%   'UseParallel'    - Whether to use parallel processing (logical)
%   'Verbose'        - Display detailed progress (logical)
%   'MaxIterations'  - Maximum iterations for optimization
%   'TolFun'         - Function tolerance for convergence
%   'TolX'           - Parameter tolerance for convergence
%   'NumRestarts'    - Number of optimization restarts (multi-start)
%
% Outputs:
%   bestParams - Fitted model parameters
%   fitInfo    - Structure with fitting information (NLL, AIC, BIC, R², etc.)

    % Persistent variables for tracking optimization progress
    persistent bestNLL iterCount paramHistory;
    
    % Default parameter bounds if not provided
    defaultLB = [0.01, 0.01, 0.01, 0.01, 40, 3, 50, 0.3, 1];
    defaultUB = [0.99, 0.99, 0.99, 0.99, 120, 12, 150, 0.7, 100];
    
    % Parse inputs with enhanced options
    p = inputParser;
    addParameter(p, 'LowerBounds', defaultLB, @isnumeric);
    addParameter(p, 'UpperBounds', defaultUB, @isnumeric);
    addParameter(p, 'ModelOptions', struct('nTrials', 5000, 'useGPU', true, 'smartAllocation', true), @isstruct);
    addParameter(p, 'OptimMethod', 'PSO_Multi', @ischar);
    addParameter(p, 'UseParallel', false, @islogical);
    addParameter(p, 'Verbose', true, @islogical);
    addParameter(p, 'MaxIterations', 2000, @isnumeric);
    addParameter(p, 'TolFun', 1e-8, @isnumeric);
    addParameter(p, 'TolX', 1e-8, @isnumeric);
    addParameter(p, 'NumRestarts', 3, @isnumeric);
    parse(p, varargin{:});
    
    lb = p.Results.LowerBounds;
    ub = p.Results.UpperBounds;
    modelOptions = p.Results.ModelOptions;
    optimMethod = p.Results.OptimMethod;
    useParallel = p.Results.UseParallel;
    verbose = p.Results.Verbose;
    maxIter = p.Results.MaxIterations;
    tolFun = p.Results.TolFun;
    tolX = p.Results.TolX;
    numRestarts = p.Results.NumRestarts;
    
    % Model
    minBound=1e-10;
    maxBound=100-minBound;
    mdlBaseline = @(x, params) max(minBound, min(maxBound, normMdlSimClaude(x, params, 'baseline',modelOptions))); % zero opto input
    mdlOpto = @(x, params) max(minBound, min(maxBound, normMdlSimClaude(x, params, 'opto',modelOptions))); % nonzero opto input

    % Reset persistent variables
    bestNLL = Inf;
    iterCount = 0;
    paramHistory = [];
    
    % Ensure initial parameters are within bounds
    initialGuess = max(min(initialGuess, ub), lb);
    
    % Initialize output structure
    fitInfo = struct();
    fitInfo.methodUsed = optimMethod;
    fitInfo.initialParams = initialGuess;
    fitInfo.mdlBaseline=mdlBaseline;
    fitInfo.mdlOpto=mdlOpto;
    fitInfo.lowerBounds = lb;
    fitInfo.upperBounds = ub;
    fitInfo.modelOptions = modelOptions;
    
    % Convert percentage correct to trial counts (assuming 20 trials per condition)
    [sumBaselineTrials, successBaselineTrials] = convertToTrialCounts(xDataBaseline, yDataBaseline);
    [sumOptoTrials, successOptoTrials] = convertToTrialCounts(xDataOpto, yDataOpto);
    
    % Define objective function (negative log-likelihood) with enhanced error handling
    function [nll, gradient] = objectiveNLL(params)
        % Track iterations
        iterCount = iterCount + 1;
        
        % Store parameter history
        paramHistory(end+1,:) = params;
        
        % Get model predictions
        try
            predBaseline = mdlBaseline(xDataBaseline, params);
            predOpto = mdlOpto(xDataOpto, params);
            
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
            
            % Display progress if verbose and significant improvement
            if nll < bestNLL * 0.99  % Only report >1% improvements
                bestNLL = nll;
                if verbose
                    fprintf('Iteration %d: NLL = %.4f, [a=%.3f b=%.3f l=%.3f w=%.3f g0=%.1f n=%.1f rmx=%.1f e=%.3f o=%.1f]\n', ...
                        iterCount, nll, params(1), params(2), params(3), params(4), params(5), ...
                        params(6), params(7), params(8), params(9));
                end
            end
            
            % Periodically report progress
            if verbose && mod(iterCount, 50) == 0
                fprintf('... Iteration %d: Current NLL = %.4f (best: %.4f)\n', iterCount, nll, bestNLL);
            end
            
            % Compute numerical gradient if requested
            if nargout > 1
                gradient = zeros(size(params));
                h = 1e-5;  % Step size for numerical differentiation
                
                for i = 1:length(params)
                    paramsPlus = params;
                    paramsPlus(i) = paramsPlus(i) + h;
                    
                    % Ensure within bounds
                    paramsPlus = max(min(paramsPlus, ub), lb);
                    
                    predBaselinePlus = mdlBaseline(xDataBaseline, paramsPlus, modelOptions);
                    predOptoPlus = mdlOpto(xDataOpto, paramsPlus, modelOptions);
                    
                    predBaselinePlus = max(0.001, min(0.999, predBaselinePlus/100));
                    predOptoPlus = max(0.001, min(0.999, predOptoPlus/100));
                    
                    nllBaselinePlus = -sum(successBaselineTrials .* log(predBaselinePlus) + ...
                                     (sumBaselineTrials - successBaselineTrials) .* log(1 - predBaselinePlus));
                    nllOptoPlus = -sum(successOptoTrials .* log(predOptoPlus) + ...
                                 (sumOptoTrials - successOptoTrials) .* log(1 - predOptoPlus));
                    
                    nllPlus = nllBaselinePlus + nllOptoPlus;
                    
                    if ~isfinite(nllPlus)
                        nllPlus = 1e10;
                    end
                    
                    gradient(i) = (nllPlus - nll) / h;
                end
            end
            
        catch ME
            % If error occurs, return large value
            nll = 1e10;
            if nargout > 1
                gradient = zeros(size(params));
            end
            
            if verbose && mod(iterCount, 100) == 0
                fprintf('Error in objective function: %s\n', ME.message);
            end
        end
    end
    
    % Start timing
    if verbose
        fprintf('Starting parameter fitting using %s method...\n', optimMethod);
        tic;
    end
    
    % Perform optimization based on selected method
    switch optimMethod
        case 'PSO'
            % Particle Swarm Optimization (global)
            if verbose
                fprintf('Running PSO global optimization...\n');
            end
            
            % Setup PSO options
            options = optimoptions('particleswarm', ...
                'Display', ifelse(verbose, 'iter', 'off'), ...
                'SwarmSize', 100, ...
                'MaxIterations', maxIter, ...
                'FunctionTolerance', tolFun, ...
                'UseParallel', useParallel, ...
                'MaxStallIterations', 50);
            
            % Run PSO
            [bestParams, nll] = particleswarm(@objectiveNLL, length(initialGuess), lb, ub, options);
            fitInfo.finalNLL = nll;
            
        case 'PSO_Multi'
            % Enhanced approach with PSO followed by multiple local optimization runs
            if verbose
                fprintf('Running PSO global search followed by multi-start local optimization...\n');
            end
            
            % Run PSO first
            options = optimoptions('particleswarm', ...
                'Display', ifelse(verbose, 'iter', 'off'), ...
                'SwarmSize', 50, ...
                'MaxIterations', ceil(maxIter/4), ...
                'FunctionTolerance', tolFun*10, ...
                'UseParallel', useParallel, ...
                'MaxStallIterations', 20);
            
            % Run PSO to get a good starting point
            [psoParams, psoNLL] = particleswarm(@objectiveNLL, length(initialGuess), lb, ub, options);
            fitInfo.psoNLL = psoNLL;
            fitInfo.psoParams = psoParams;
            
            if verbose
                fprintf('PSO completed: NLL = %.4f\n', psoNLL);
                fprintf('Running multi-start local optimization...\n');
            end
            
            % Now run multiple local optimizations from different starting points
            bestOptimNLL = Inf;
            bestOptimParams = psoParams;
            
            % Generate starting points: PSO result, initial guess, and perturbed points
            startingPoints = zeros(numRestarts, length(initialGuess));
            startingPoints(1,:) = psoParams;  % PSO result
            
            % Random perturbations of the best PSO result
            for i = 2:numRestarts
                % Generate each parameter as random point between bounds
                % For some parameters, use log-scaled sampling
                for j = 1:length(initialGuess)
                    if j >= 5 && j <= 7  % g0, n, rmx - log scale sampling
                        logLB = log(lb(j));
                        logUB = log(ub(j));
                        startingPoints(i,j) = exp(logLB + (logUB-logLB)*rand());
                    else
                        startingPoints(i,j) = lb(j) + (ub(j)-lb(j))*rand();
                    end
                end
            end
            
            % Setup fminsearch options
            options = optimset('fminsearch');
            options = optimset(options, 'MaxIter', maxIter);
            options = optimset(options, 'MaxFunEvals', maxIter * 2);
            options = optimset(options, 'TolX', tolX);
            options = optimset(options, 'TolFun', tolFun);
            options = optimset(options, 'Display', ifelse(verbose, 'iter', 'off'));
            
            % Run local optimization from each starting point
            for i = 1:numRestarts
                if verbose
                    fprintf('Local optimization run %d/%d...\n', i, numRestarts);
                end
                
                try
                    [localParams, localNLL] = fminsearchbnd(@objectiveNLL, startingPoints(i,:), lb, ub, options);
                    
                    if localNLL < bestOptimNLL
                        bestOptimNLL = localNLL;
                        bestOptimParams = localParams;
                        if verbose
                            fprintf('New best local solution found: NLL = %.4f\n', bestOptimNLL);
                        end
                    end
                catch ME
                    if verbose
                        fprintf('Local optimization %d failed: %s\n', i, ME.message);
                    end
                end
            end
            
            % Use the best result
            bestParams = bestOptimParams;
            fitInfo.finalNLL = bestOptimNLL;
            fitInfo.numLocalRuns = numRestarts;
            
        case 'Global'
            % Global optimization using Genetic Algorithm
            if verbose
                fprintf('Running global optimization with GA...\n');
            end
            
            % Setup GA options
            gaOptions = optimoptions('ga', ...
                'Display', ifelse(verbose, 'iter', 'off'), ...
                'PopulationSize', 50, ...
                'MaxGenerations', ceil(maxIter/50), ...
                'UseParallel', useParallel, ...
                'MaxStallGenerations', 20, ...
                'FunctionTolerance', tolFun);
            
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
            options = optimset(options, 'TolX', tolX);
            options = optimset(options, 'TolFun', tolFun);
            options = optimset(options, 'Display', ifelse(verbose, 'iter', 'off'));
            
            % Run bounded fminsearch
            [bestParams, nll] = fminsearchbnd(@objectiveNLL, initialGuess, lb, ub, options);
            fitInfo.finalNLL = nll;
            
        case 'GlobalThenLocal'
            if verbose
                fprintf('Running initial grid search...\n');
            end
            
            % More extensive grid search with logarithmic spacing for some parameters
            % Focus on most influential parameters
            g0_values = logspace(log10(lb(5)), log10(ub(5)), 4);
            n_values = linspace(lb(6), ub(6), 4);
            e_values = linspace(lb(8), ub(8), 3);
            o_values = logspace(log10(max(1, lb(9))), log10(ub(9)), 3);
            
            % Initialize tracking
            bestNLL = Inf;
            bestGridParams = initialGuess;
            
            % Precompute total grid search iterations
            totalIters = length(g0_values) * length(n_values) * length(e_values) * length(o_values);
            
            if verbose
                fprintf('Grid search: testing %d parameter combinations...\n', totalIters);
            end
            
            % Using a more efficient grid search implementation
            counter = 0;
            bestFound = false;
            
            % Vectorize the grid search to avoid nested loops
            [G0, N, E, O] = ndgrid(g0_values, n_values, e_values, o_values);
            
            for i = 1:numel(G0)
                counter = counter + 1;
                
                % Create test parameters
                testParams = initialGuess;
                testParams(5) = G0(i);  % g0
                testParams(6) = N(i);   % n
                testParams(8) = E(i);   % e
                testParams(9) = O(i);   % o
                
                % Calculate NLL
                nll = objectiveNLL(testParams);
                
                % Update best parameters if better
                if nll < bestNLL
                    bestNLL = nll;
                    bestGridParams = testParams;
                    bestFound = true;
                    
                    if verbose && mod(counter, 10) == 0
                        fprintf('  Grid point %d/%d: Found better parameters, NLL = %.4f\n', ...
                            counter, totalIters, bestNLL);
                    end
                end
                
                % Show periodic progress
                if verbose && mod(counter, ceil(totalIters/10)) == 0
                    fprintf('  Grid search progress: %d/%d combinations (%.1f%%)\n', ...
                        counter, totalIters, (counter/totalIters)*100);
                end
            end
            
            % If grid search didn't find better parameters, use initial guess
            if ~bestFound
                bestGridParams = initialGuess;
                if verbose
                    fprintf('Grid search did not find better parameters than initial guess\n');
                end
            end
            
            % Set grid search results for refinement
            initialGuess = bestGridParams;
            fitInfo.gridSearchNLL = bestNLL;
            fitInfo.gridSearchParams = bestGridParams;
            
            if verbose
                fprintf('\nRunning local optimization from grid search result...\n');
            end
            
            % Setup fminsearch options for more thorough search
            options = optimset('fminsearch');
            options = optimset(options, 'MaxIter', maxIter);
            options = optimset(options, 'MaxFunEvals', maxIter * 2);
            options = optimset(options, 'TolX', tolX);
            options = optimset(options, 'TolFun', tolFun);
            options = optimset(options, 'Display', ifelse(verbose, 'iter', 'off'));
            
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
    predBaseline = mdlBaseline(xDataBaseline, bestParams);
    predOpto = mdlOpto(xDataOpto, bestParams);
    
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
    fitInfo.iterations = iterCount;
    fitInfo.optimTime = toc; % Store optimization time
    
    % Calculate parameter changes from initial values
    fitInfo.paramChange = bestParams - initialGuess;
    fitInfo.percentChange = 100 * (bestParams - initialGuess) ./ abs(initialGuess);
    
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
        
        % Report parameter changes
        fprintf('\nParameter changes from initial values:\n');
        paramNames = {'a', 'b', 'l', 'w', 'g0', 'n', 'rmx', 'e', 'o'};
        for i = 1:length(initialGuess)
            change = bestParams(i) - initialGuess(i);
            percentChange = (change / abs(initialGuess(i))) * 100;
            fprintf('  %s: %.4f → %.4f (%.2f%%)\n', paramNames{i}, initialGuess(i), bestParams(i), percentChange);
        end
    end
end

% Helper Functions

% Function to handle conditional expression (replacement for a?b:c in other languages)
function result = ifelse(condition, trueVal, falseVal)
    if condition
        result = trueVal;
    else
        result = falseVal;
    end
end

function [sumTrials, successTrials] = convertToTrialCounts(x, y)
    % Convert percentage correct to trial counts
    % Assumes 20 trials per contrast level by default
    
    % Handle empty inputs
    if isempty(x) || isempty(y)
        sumTrials = [];
        successTrials = [];
        return;
    end
    
    % Remove NaN values
    validIdx = ~isnan(x) & ~isnan(y);
    x = x(validIdx);
    y = y(validIdx);
    
    % Make sure x and y are row vectors
    x = x(:)';
    y = y(:)';
    
    % Create trial counts (assuming 20 trials per condition)
    sumTrials = 20 * ones(size(x));
    
    % Convert percentage to success counts
    successTrials = (y / 100) .* sumTrials;
end