function [bestParams, fitInfo] = multiStageOptimization(xDataBaseline, yDataBaseline, xDataOpto, yDataOpto, initialGuess, lb, ub, modelOptions, varargin)
% MULTISTAGEOPTIMIZATION - Simple optimization for normalization model fitting using fminsearchbnd
%
% Inputs:
%   xDataBaseline, yDataBaseline - Baseline condition data
%   xDataOpto, yDataOpto - Opto condition data
%   initialGuess - Initial parameter values [a, b, l, w, g0, n, rmx, e, o]
%   lb, ub - Lower and upper bounds for parameters
%   modelOptions - Options for normMdlSimClaude
%   varargin - Optional name-value pairs:
%       'Verbose' - Display progress (logical, default: true)
%
% Outputs:
%   bestParams - Fitted model parameters
%   fitInfo - Structure with fitting information

    % Parse inputs
    p = inputParser;
    addParameter(p, 'Verbose', true, @islogical);
    parse(p, varargin{:});
    verbose = p.Results.Verbose;
    
    % Initialize tracking variables
    paramHistory = initialGuess;
    bestNLL = Inf;
    iterCount = 0;
    
    % Define model functions
    minBound = 1e-10;
    maxBound = 100 - minBound;
    mdlBaseline = @(x, params) max(minBound, min(maxBound, normMdlSimClaude(x, params, 'baseline', modelOptions)));
    mdlOpto = @(x, params) max(minBound, min(maxBound, normMdlSimClaude(x, params, 'opto', modelOptions)));
    
    % Convert percentage correct to trial counts
    [sumBaselineTrials, successBaselineTrials] = convertToTrialCounts(xDataBaseline, yDataBaseline);
    [sumOptoTrials, successOptoTrials] = convertToTrialCounts(xDataOpto, yDataOpto);
    
    % Objective function (negative log-likelihood)
    function [nll, gradient] = objectiveNLL(params)
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
            
            % Handle invalid values
            if ~isfinite(nll)
                nll = 1e10;
            end
            
            % Track best NLL and display progress
            if nll < bestNLL
                bestNLL = nll;
                if verbose
                    fprintf('Iteration %d: NLL = %.4f, [a=%.3f b=%.3f l=%.3f w=%.3f g0=%.1f n=%.1f rmx=%.1f e=%.3f o=%.1f]\n', ...
                        iterCount, nll, params(1), params(2), params(3), params(4), params(5), ...
                        params(6), params(7), params(8), params(9));
                end
            end
            
            % Compute gradient if requested
            if nargout > 1
                gradient = zeros(size(params));
                h = 1e-5;
                
                for i = 1:length(params)
                    paramsPlus = params;
                    paramsPlus(i) = min(ub(i), params(i) + h);
                    nllPlus = objectiveNLL(paramsPlus);
                    gradient(i) = (nllPlus - nll) / h;
                end
            end
            
        catch ME
            nll = 1e10;
            if nargout > 1
                gradient = zeros(size(params));
            end
            if verbose
                fprintf('Error in objective function: %s\n', ME.message);
            end
        end
    end
    
    % Optimization options
    options = optimset('fminsearch');
    options = optimset(options, 'MaxIter', 5000);
    options = optimset(options, 'MaxFunEvals', 10000);
    options = optimset(options, 'TolX', 1e-8);
    options = optimset(options, 'TolFun', 1e-8);
    options = optimset(options, 'Display', ifelse(verbose, 'iter', 'off'));
    
    % Run optimization
    [bestParams, finalNLL] = fminsearchbnd(@objectiveNLL, initialGuess, lb, ub, options);
    
    % Calculate model statistics
    totalTrials = sum(sumBaselineTrials) + sum(sumOptoTrials);
    numParams = length(bestParams);
    
    % Calculate information criteria
    aic = 2 * numParams + 2 * finalNLL;
    aicc = aic + (2 * numParams * (numParams + 1)) / (totalTrials - numParams - 1);
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
    
    % Store results
    fitInfo = struct();
    fitInfo.finalNLL = finalNLL;
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
    
    % Parameter changes from initial values
    fitInfo.paramChangeFromInitial = bestParams - initialGuess;
    fitInfo.percentChangeFromInitial = 100 * (bestParams - initialGuess) ./ abs(initialGuess);
    
    % Display summary if verbose
    if verbose
        fprintf('\n=== Optimization Complete ===\n');
        fprintf('Final NLL: %.4f\n', finalNLL);
        fprintf('AIC: %.2f, AICc: %.2f, BIC: %.2f\n', aic, aicc, bic);
        fprintf('R-squared (baseline): %.4f, R-squared (opto): %.4f\n', rSquared_baseline, rSquared_opto);
        
        % Show parameter changes
        fprintf('\nParameter changes from initial values:\n');
        paramNames = {'a', 'b', 'l', 'w', 'g0', 'n', 'rmx', 'e', 'o'};
        for i = 1:length(paramNames)
            fprintf('%-5s: %.4f → %.4f (%.2f%%)\n', ...
                paramNames{i}, initialGuess(i), bestParams(i), fitInfo.percentChangeFromInitial(i));
        end
    end
end

function result = ifelse(condition, trueVal, falseVal)
    if condition
        result = trueVal;
    else
        result = falseVal;
    end
end

function [sumTrials, successTrials] = convertToTrialCounts(x, y)
    if isempty(x) || isempty(y)
        sumTrials = [];
        successTrials = [];
        return;
    end
    
    validIdx = ~isnan(x) & ~isnan(y);
    x = x(validIdx);
    y = y(validIdx);
    
    x = x(:)';
    y = y(:)';
    
    sumTrials = 20 * ones(size(x));
    successTrials = (y / 100) .* sumTrials;
end