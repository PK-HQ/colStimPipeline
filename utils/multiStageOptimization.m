function [bestParams, fitInfo] = multiStageOptimization(xDataBaseline, yDataBaseline, xDataOpto, yDataOpto, initialGuess, lb, ub, modelOptions, varargin)
% MULTISTAGEOPTIMIZATION - Three-stage optimization for normalization model fitting
%
% Implements a three-stage optimization strategy:
% 1. Global exploration using PSO or GA with parameter scaling
% 2. Local refinement using multiple starting points
% 3. Final polishing with tight convergence criteria
%
% Usage:
%   [bestParams, fitInfo] = multiStageOptimization(xDataBaseline, yDataBaseline, xDataOpto, yDataOpto, initialGuess, lb, ub, modelOptions)
%   [bestParams, fitInfo] = multiStageOptimization(..., 'OptionName', OptionValue)
%
% Optional Name-Value Parameters:
%   'UseParallel'    - Whether to use parallel processing (logical, default: false)
%   'Verbose'        - Display detailed progress (logical, default: true)
%   'NumRestarts'    - Number of optimization restarts in stage 2 (default: 5)
%   'GlobalMethod'   - Global optimization method: 'PSO', 'GA' (default: 'PSO')
%   'PreserveTiming' - Track timing performance (default: true)
%
% Outputs:
%   bestParams - Fitted model parameters
%   fitInfo    - Structure with fitting information and optimization details

    % Persistent variables for tracking optimization progress
    persistent bestNLL iterCount paramHistory stageTimes;
    
    % Parse inputs
    p = inputParser;
    addParameter(p, 'UseParallel', false, @islogical);
    addParameter(p, 'Verbose', true, @islogical);
    addParameter(p, 'NumRestarts', 5, @isnumeric);
    addParameter(p, 'GlobalMethod', 'PSO', @ischar);
    addParameter(p, 'PreserveTiming', true, @islogical);
    parse(p, varargin{:});
    
    useParallel = p.Results.UseParallel;
    verbose = p.Results.Verbose;
    numRestarts = p.Results.NumRestarts;
    globalMethod = p.Results.GlobalMethod;
    preserveTiming = p.Results.PreserveTiming;
    
    % Initialize timing tracking
    if preserveTiming
        stageTimes = struct('sensitivity', 0, 'global', 0, 'local', 0, 'polish', 0);
    end
    
    % Define model functions
    minBound = 1e-10;
    maxBound = 100 - minBound;
    mdlBaseline = @(x, params) max(minBound, min(maxBound, normMdlSimClaude(x, params, 'baseline', modelOptions)));
    mdlOpto = @(x, params) max(minBound, min(maxBound, normMdlSimClaude(x, params, 'opto', modelOptions)));
    
    % Start with parameter sensitivity analysis to guide optimization
    if verbose
        fprintf('=== Stage 0: Parameter Sensitivity Analysis ===\n');
        tic;
    end
    
    % Convert percentage correct to trial counts
    [sumBaselineTrials, successBaselineTrials] = convertToTrialCounts(xDataBaseline, yDataBaseline);
    [sumOptoTrials, successOptoTrials] = convertToTrialCounts(xDataOpto, yDataOpto);
    
    % Reset tracking variables
    bestNLL = Inf;
    iterCount = 0;
    paramHistory = [];
    
    % Ensure initial parameters are within bounds
    initialGuess = max(min(initialGuess, ub), lb);
    
    % Define the full objective function (for non-scaled parameters)
    function [nll, gradient] = objectiveNLLFull(params)
        % Track iterations
        iterCount = iterCount + 1;
        
        % Store parameter history
        paramHistory(end+1,:) = params;
        
        try
            % Get model predictions
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
                    
                    predBaselinePlus = mdlBaseline(xDataBaseline, paramsPlus);
                    predOptoPlus = mdlOpto(xDataOpto, paramsPlus);
                    
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
    
    % Define scaled objective function for global optimization
    function nll = objectiveNLLScaled(scaledParams)
        % Unscale parameters for model evaluation
        unscaledParams = unscaleParameters(scaledParams, lb, ub);
        
        % Call the unscaled objective function
        nll = objectiveNLLFull(unscaledParams);
    end
    
    % Run quick sensitivity analysis to identify most important parameters
    [sensitivity, paramRanking] = analyzeParameterSensitivity(xDataBaseline, yDataBaseline, xDataOpto, yDataOpto, initialGuess, lb, ub, modelOptions);
    
    if preserveTiming
        stageTimes.sensitivity = toc;
    end
    
    % Stage 1: Global exploration
    if verbose
        fprintf('\n=== Stage 1: Global Exploration ===\n');
        tic;
    end
    
    % Scale initial guess for global optimization
    scaledInitialGuess = scaleParameters(initialGuess, lb, ub);
    
    % Define bounds for scaled parameters (always [0,1])
    scaledLB = zeros(size(lb));
    scaledUB = ones(size(ub));
    
    % Global optimization method selection
    switch globalMethod
        case 'PSO'
            if verbose
                fprintf('Running PSO global optimization with scaled parameters...\n');
            end
            
            % Setup PSO options with larger swarm for better exploration
            psoOptions = optimoptions('particleswarm', ...
                'Display', ifelse(verbose, 'iter', 'off'), ...
                'SwarmSize', 100, ...  % Increased swarm size
                'MaxIterations', 200, ...
                'FunctionTolerance', 1e-6, ...
                'UseParallel', useParallel, ...
                'InitialSwarmMatrix', repmat(scaledInitialGuess, 100, 1), ... % Use scaled initial guess
                'MaxStallIterations', 30);
            
            % Run PSO with scaled parameters
            [scaledGlobalParams, globalNLL] = particleswarm(@objectiveNLLScaled, length(scaledInitialGuess), scaledLB, scaledUB, psoOptions);
            
            % Unscale parameters
            globalParams = unscaleParameters(scaledGlobalParams, lb, ub);
            
        case 'GA'
            if verbose
                fprintf('Running GA global optimization with scaled parameters...\n');
            end
            
            % Setup GA options
            gaOptions = optimoptions('ga', ...
                'Display', ifelse(verbose, 'iter', 'off'), ...
                'PopulationSize', 100, ...
                'MaxGenerations', 100, ...
                'FunctionTolerance', 1e-6, ...
                'UseParallel', useParallel, ...
                'InitialPopulationMatrix', scaledInitialGuess, ... % Use scaled initial guess
                'MaxStallGenerations', 20);
            
            % Run GA with scaled parameters
            [scaledGlobalParams, globalNLL] = ga(@objectiveNLLScaled, length(scaledInitialGuess), [], [], [], [], scaledLB, scaledUB, [], gaOptions);
            
            % Unscale parameters
            globalParams = unscaleParameters(scaledGlobalParams, lb, ub);
            
        otherwise
            error('Unknown global optimization method: %s', globalMethod);
    end
    
    if verbose
        fprintf('Global optimization completed: NLL = %.4f\n', globalNLL);
        fprintf('Global best parameters: [a=%.3f b=%.3f l=%.3f w=%.3f g0=%.1f n=%.1f rmx=%.1f e=%.3f o=%.1f]\n', ...
            globalParams(1), globalParams(2), globalParams(3), globalParams(4), globalParams(5), ...
            globalParams(6), globalParams(7), globalParams(8), globalParams(9));
    end
    
    if preserveTiming
        stageTimes.global = toc;
    end
    
    % Stage 2: Local refinement from multiple starting points
    if verbose
        fprintf('\n=== Stage 2: Local Refinement ===\n');
        tic;
    end
    
    % Setup starting points: 
    % 1. Global result
    % 2. Original initial guess
    % 3. Additional points with perturbations
    startingPoints = zeros(numRestarts, length(initialGuess));
    startingPoints(1,:) = globalParams;        % Global optimization result
    startingPoints(2,:) = initialGuess;        % Original initial guess
    
    % Generate additional points with intelligent perturbation based on sensitivity
    for i = 3:numRestarts
        newPoint = globalParams;
        
        % Perturb more sensitive parameters more significantly
        for j = 1:length(paramRanking)
            paramIdx = paramRanking(j);
            
            % Linear decay of perturbation amount based on sensitivity ranking
            perturbFactor = 0.5 * (numRestarts - j) / numRestarts;
            
            % Add random perturbation based on parameter range
            if paramIdx >= 5 && paramIdx <= 7  % Log-scale params (g0, n, rmx)
                % Log-scale perturbation
                logParam = log(newPoint(paramIdx));
                logLB = log(lb(paramIdx));
                logUB = log(ub(paramIdx));
                logRange = logUB - logLB;
                
                % Apply perturbation in log space
                perturbation = (rand() * 2 - 1) * perturbFactor * logRange;
                newLogValue = logParam + perturbation;
                
                % Convert back and ensure within bounds
                newPoint(paramIdx) = exp(max(logLB, min(logUB, newLogValue)));
            else
                % Linear scale perturbation
                paramRange = ub(paramIdx) - lb(paramIdx);
                perturbation = (rand() * 2 - 1) * perturbFactor * paramRange;
                newPoint(paramIdx) = max(lb(paramIdx), min(ub(paramIdx), newPoint(paramIdx) + perturbation));
            end
        end
        
        startingPoints(i,:) = newPoint;
    end
    
    % Setup options for local optimization stage
    localOptions = optimset('fminsearch');
    localOptions = optimset(localOptions, 'MaxIter', 2000);
    localOptions = optimset(localOptions, 'MaxFunEvals', 4000);
    localOptions = optimset(localOptions, 'TolX', 1e-6);
    localOptions = optimset(localOptions, 'TolFun', 1e-6);
    localOptions = optimset(localOptions, 'Display', ifelse(verbose, 'iter', 'off'));
    
    % Run local optimization from each starting point
    localResults = cell(numRestarts, 2);  % Store [params, nll]
    
    for i = 1:numRestarts
        if verbose
            fprintf('Local optimization run %d/%d...\n', i, numRestarts);
        end
        
        try
            % Run bounded fminsearch
            [localParams, localNLL] = fminsearchbnd(@objectiveNLLFull, startingPoints(i,:), lb, ub, localOptions);
            
            % Store results
            localResults{i,1} = localParams;
            localResults{i,2} = localNLL;
            
            if verbose
                fprintf('  Run %d completed: NLL = %.4f\n', i, localNLL);
            end
        catch ME
            if verbose
                fprintf('  Run %d failed: %s\n', i, ME.message);
            end
            localResults{i,1} = startingPoints(i,:);
            localResults{i,2} = 1e10;
        end
    end
    
    % Find best local result
    bestLocalNLL = Inf;
    bestLocalIdx = 0;
    
    for i = 1:numRestarts
        if localResults{i,2} < bestLocalNLL
            bestLocalNLL = localResults{i,2};
            bestLocalIdx = i;
        end
    end
    
    bestLocalParams = localResults{bestLocalIdx,1};
    
    if verbose
        fprintf('Best local optimization: Run %d, NLL = %.4f\n', bestLocalIdx, bestLocalNLL);
        fprintf('Local best parameters: [a=%.3f b=%.3f l=%.3f w=%.3f g0=%.1f n=%.1f rmx=%.1f e=%.3f o=%.1f]\n', ...
            bestLocalParams(1), bestLocalParams(2), bestLocalParams(3), bestLocalParams(4), bestLocalParams(5), ...
            bestLocalParams(6), bestLocalParams(7), bestLocalParams(8), bestLocalParams(9));
    end
    
    if preserveTiming
        stageTimes.local = toc;
    end
    
    % Stage 3: Final polishing with tighter convergence criteria
    if verbose
        fprintf('\n=== Stage 3: Final Polishing ===\n');
        tic;
    end
    
    % Setup options for final polishing
    polishOptions = optimset('fminsearch');
    polishOptions = optimset(polishOptions, 'MaxIter', 5000);
    polishOptions = optimset(polishOptions, 'MaxFunEvals', 10000);
    polishOptions = optimset(polishOptions, 'TolX', 1e-8);
    polishOptions = optimset(polishOptions, 'TolFun', 1e-8);
    polishOptions = optimset(polishOptions, 'Display', ifelse(verbose, 'iter', 'off'));
    
    % Run final polish optimization
    try
        [bestParams, finalNLL] = fminsearchbnd(@objectiveNLLFull, bestLocalParams, lb, ub, polishOptions);
        
        if verbose
            fprintf('Final polishing completed: NLL = %.4f\n', finalNLL);
            fprintf('Final parameters: [a=%.3f b=%.3f l=%.3f w=%.3f g0=%.1f n=%.1f rmx=%.1f e=%.3f o=%.1f]\n', ...
                bestParams(1), bestParams(2), bestParams(3), bestParams(4), bestParams(5), ...
                bestParams(6), bestParams(7), bestParams(8), bestParams(9));
        end
    catch ME
        if verbose
            fprintf('Final polishing failed: %s\n', ME.message);
            fprintf('Using best local parameters instead.\n');
        end
        bestParams = bestLocalParams;
        finalNLL = bestLocalNLL;
    end
    
    if preserveTiming
        stageTimes.polish = toc;
    end
    
    % Calculate model statistics
    totalTrials = sum(sumBaselineTrials) + sum(sumOptoTrials);
    numParams = length(bestParams);
    
    % Calculate AIC and BIC
    aic = 2 * numParams + 2 * finalNLL;
    if totalTrials > numParams + 1
        aicc = aic + (2 * numParams * (numParams + 1)) / (totalTrials - numParams - 1);
    else
        aicc = Inf;
    end
    bic = log(totalTrials) * numParams + 2 * finalNLL;
    
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
    
    % Store results in fitInfo structure
    fitInfo = struct();
    fitInfo.finalNLL = finalNLL;
    fitInfo.globalNLL = globalNLL;
    fitInfo.localNLL = bestLocalNLL;
    fitInfo.globalParams = globalParams;
    fitInfo.localParams = bestLocalParams;
    fitInfo.aic = aic;
    fitInfo.aicc = aicc;
    fitInfo.bic = bic;
    fitInfo.mdlBaseline = mdlBaseline;
    fitInfo.mdlOpto = mdlOpto;
    fitInfo.rSquared_baseline = rSquared_baseline;
    fitInfo.rSquared_opto = rSquared_opto;
    fitInfo.predictedBaseline = predBaseline;
    fitInfo.predictedOpto = predOpto;
    fitInfo.iterations = iterCount;
    fitInfo.paramHistory = paramHistory;
    fitInfo.sensitivity = sensitivity;
    fitInfo.paramRanking = paramRanking;
    
    % Parameter changes from initial values
    fitInfo.paramChangeFromInitial = bestParams - initialGuess;
    fitInfo.percentChangeFromInitial = 100 * (bestParams - initialGuess) ./ abs(initialGuess);
    
    % Also track changes between optimization stages
    fitInfo.paramChangeGlobalToLocal = bestLocalParams - globalParams;
    fitInfo.paramChangeLocalToFinal = bestParams - bestLocalParams;
    
    % Store timing information if available
    if preserveTiming
        fitInfo.stageTimes = stageTimes;
        fitInfo.totalTime = sum(struct2array(stageTimes));
    end
    
    % Display summary if verbose
    if verbose
        fprintf('\n=== Optimization Complete ===\n');
        fprintf('Initial NLL: %.4f\n', objectiveNLLFull(initialGuess));
        fprintf('Global optimization NLL: %.4f\n', globalNLL);
        fprintf('Local optimization NLL: %.4f\n', bestLocalNLL);
        fprintf('Final NLL: %.4f\n', finalNLL);
        fprintf('AIC: %.2f, AICc: %.2f, BIC: %.2f\n', aic, aicc, bic);
        fprintf('R-squared (baseline): %.4f, R-squared (opto): %.4f\n', rSquared_baseline, rSquared_opto);
        
        % Show parameter comparison
        fprintf('\nParameter comparison through optimization stages:\n');
        paramNames = {'a', 'b', 'l', 'w', 'g0', 'n', 'rmx', 'e', 'o'};
        fprintf('%-5s %-10s %-10s %-10s %-10s\n', 'Param', 'Initial', 'Global', 'Local', 'Final');
        for i = 1:length(paramNames)
            fprintf('%-5s %-10.4f %-10.4f %-10.4f %-10.4f\n', ...
                paramNames{i}, initialGuess(i), globalParams(i), bestLocalParams(i), bestParams(i));
        end
        
        % Show timing information
        if preserveTiming
            fprintf('\nTime breakdown (seconds):\n');
            fprintf('Sensitivity analysis: %.2f\n', stageTimes.sensitivity);
            fprintf('Global optimization: %.2f\n', stageTimes.global);
            fprintf('Local refinement: %.2f\n', stageTimes.local);
            fprintf('Final polishing: %.2f\n', stageTimes.polish);
            fprintf('Total time: %.2f\n', fitInfo.totalTime);
        end
    end
end

function result = ifelse(condition, trueVal, falseVal)
    % Simple utility function for conditional expression
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