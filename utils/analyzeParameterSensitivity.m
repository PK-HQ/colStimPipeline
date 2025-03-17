function [sensitivity, paramRanking] = analyzeParameterSensitivity(xBaseline, yBaseline, xOpto, yOpto, initialParams, lb, ub, modelOptions)
% ANALYZEPARAMETERSENSITIVITY - Analyze the sensitivity of model parameters
%
% Inputs:
%   xBaseline, yBaseline - Baseline condition data
%   xOpto, yOpto - Opto condition data
%   initialParams - Initial parameter values [a, b, l, w, g0, n, rmx, e, o]
%   lb, ub - Lower and upper bounds for parameters
%   modelOptions - Options for normMdlSim* (struct)
%
% Outputs:
%   sensitivity - Matrix of sensitivity scores for each parameter
%   paramRanking - Parameters ranked by sensitivity (most to least)

    % Model functions
    minBound = 1e-10;
    maxBound = 100 - minBound;
    mdlBaseline = @(x, params) max(minBound, min(maxBound, normMdlSimClaude(x, params, 'baseline', modelOptions)));
    mdlOpto = @(x, params) max(minBound, min(maxBound, normMdlSimClaude(x, params, 'opto', modelOptions)));
    
    % Convert to trial counts for likelihood calculation
    [sumBaselineTrials, successBaselineTrials] = convertToTrialCounts(xBaseline, yBaseline);
    [sumOptoTrials, successOptoTrials] = convertToTrialCounts(xOpto, yOpto);
    
    % Define negative log-likelihood function
    function nll = objectiveNLL(params)
        try
            % Get model predictions
            predBaseline = mdlBaseline(xBaseline, params);
            predOpto = mdlOpto(xOpto, params);
            
            % Ensure predictions are within valid range for probability
            predBaseline = max(0.001, min(0.999, predBaseline/100));
            predOpto = max(0.001, min(0.999, predOpto/100));
            
            % Calculate negative log-likelihood
            nllBaseline = -sum(successBaselineTrials .* log(predBaseline) + ...
                         (sumBaselineTrials - successBaselineTrials) .* log(1 - predBaseline));
            nllOpto = -sum(successOptoTrials .* log(predOpto) + ...
                     (sumOptoTrials - successOptoTrials) .* log(1 - predOpto));
            
            nll = nllBaseline + nllOpto;
            
            % Check for invalid values
            if ~isfinite(nll)
                nll = 1e10; % Return a large value for invalid parameters
            end
        catch
            nll = 1e10; % Return a large value for errors
        end
    end
    
    % Parameter names for reporting
    paramNames = {'a', 'b', 'l', 'w', 'g0', 'n', 'rmx', 'e', 'o'};
    numParams = length(initialParams);
    
    % Calculate baseline NLL with initial parameters
    baselineNLL = objectiveNLL(initialParams);
    fprintf('Baseline NLL with initial parameters: %.4f\n', baselineNLL);
    
    % Determine step sizes for each parameter (based on range)
    paramRanges = ub - lb;
    stepSizes = paramRanges * 0.1; % 10% steps
    
    % Initialize sensitivity matrix
    % Each row will represent one parameter
    % Columns will be [-2step, -1step, baseline, +1step, +2step]
    sensitivity = zeros(numParams, 5);
    relativeChanges = zeros(numParams, 4);
    
    fprintf('Running parameter sensitivity analysis...\n');
    
    % For each parameter, vary it while keeping others constant
    for i = 1:numParams
        testParams = initialParams;
        
        % Test -2 steps
        testParams(i) = max(lb(i), initialParams(i) - 2*stepSizes(i));
        sensitivity(i, 1) = objectiveNLL(testParams);
        
        % Test -1 step
        testParams(i) = max(lb(i), initialParams(i) - stepSizes(i));
        sensitivity(i, 2) = objectiveNLL(testParams);
        
        % Baseline (already calculated, just store)
        sensitivity(i, 3) = baselineNLL;
        
        % Test +1 step
        testParams(i) = min(ub(i), initialParams(i) + stepSizes(i));
        sensitivity(i, 4) = objectiveNLL(testParams);
        
        % Test +2 steps
        testParams(i) = min(ub(i), initialParams(i) + 2*stepSizes(i));
        sensitivity(i, 5) = objectiveNLL(testParams);
        
        % Calculate relative changes
        for j = 1:4
            relativeChanges(i, j) = abs((sensitivity(i, j+1) - sensitivity(i, j)) / baselineNLL);
        end
        
        fprintf('Parameter %s: NLLs = [%.2f, %.2f, %.2f, %.2f, %.2f]\n', ...
            paramNames{i}, sensitivity(i, 1), sensitivity(i, 2), ...
            sensitivity(i, 3), sensitivity(i, 4), sensitivity(i, 5));
    end
    
    % Calculate overall sensitivity score for each parameter (average relative change)
    sensitivityScores = mean(relativeChanges, 2);
    
    % Rank parameters by sensitivity
    [~, paramRanking] = sort(sensitivityScores, 'descend');
    
    % Display results
    fprintf('\nParameter sensitivity ranking (most to least sensitive):\n');
    for i = 1:numParams
        idx = paramRanking(i);
        fprintf('%d. %s (score: %.4f)\n', i, paramNames{idx}, sensitivityScores(idx));
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