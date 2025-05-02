function [BayesModel] = fitBayesianModelMLE_normMultiStart(xData, successData, sumData, modelType)
% fitBayesianModelMLE_normMultiStart
% -------------------------------------------------------------------------
% Example function that parallels the structure of fitBayesianModelMLE but:
%   1) Uses a Normal-based psychometric function,
%   2) Employs MultiStart for a more robust global optimization,
%   3) Generates heuristic initial guesses from the data.
%
% INPUTS:
%   xData{1},    xData{2}      -- stimulus intensities for baseline/opto
%   successData{1}, successData{2} -- # successes for baseline/opto
%   sumData{1}, sumData{2}     -- total # trials for baseline/opto
%   modelType                    -- string, e.g. 'normalFreeAll'
%                                   (used for consistency with your existing code)
%
% OUTPUT:
%   BayesModel -- struct containing:
%       .params       : best-fit parameters
%       .negLogLik    : negative log-likelihood at best-fit
%       .exitflag     : exit flag from optimization
%       .optimizerOutput : additional info from MultiStart/fmincon
%       .modelType    : echo of input
%       .paramLabels  : cell array describing each parameter
%
% Author: ChatGPT Example
% -------------------------------------------------------------------------

    % ---------------------------------------------------------------------
    %  Data unpacking / check
    % ---------------------------------------------------------------------
    xBaseline       = xData{1};
    xOpto           = xData{2};
    successBaseline = successData{1};
    successOpto     = successData{2};
    sumBaseline     = sumData{1};
    sumOpto         = sumData{2};

    if length(xData) ~= 2
        error('This example assumes exactly 2 conditions (baseline & opto).');
    end

    if ~exist('modelType','var') || isempty(modelType)
        modelType = 'normalFreeAll';  % default
    end

    switch lower(modelType)
        % =================================================================
        %   CASE: Normal-based free (baseline, opto), 6-parameter version
        % =================================================================
        case 'normalfreeall'
            % The typical parameter indexing might be:
            %   params = [muB, sigmaB, lapseB, muO, sigmaO, lapseO]
            %
            % Where:
            %   muB = threshold (mean) for Baseline
            %   sigmaB = slope-related (std dev) for Baseline
            %   lapseB = lapse for Baseline
            %   muO = threshold for Opto
            %   sigmaO = slope for Opto
            %   lapseO = lapse for Opto

            % 1) Define the Normal psychometric
            normalFun = @(x, mu, sigma) 0.5 + 0.5 * erf( (x - mu)./(sqrt(2).*sigma ) );

            % 2) Define the model for baseline and opto
            mdlBaseline = @(x, p) (1 - p(3)) .* normalFun(x, p(1), abs(p(2))) + 0.5 .* p(3);
            mdlOpto     = @(x, p) (1 - p(6)) .* normalFun(x, p(4), abs(p(5))) + 0.5 .* p(6);

            % 3) Negative log-likelihood function
            epsilon = 1e-8;  % small value to avoid log(0)
            objectiveFunction = @(params) negLogLik(params, ...
                xBaseline, successBaseline, sumBaseline, ...
                xOpto, successOpto, sumOpto, ...
                mdlBaseline, mdlOpto, epsilon);

            % 4) Get heuristic initial guesses based on data
            initGuess = getInitialGuessNormal( ...
                xBaseline, successBaseline, sumBaseline, ...
                xOpto, successOpto, sumOpto );

            % 5) Parameter bounds
            minX = min([xBaseline(:); xOpto(:)]);
            maxX = max([xBaseline(:); xOpto(:)]);

            lb = [ minX,  0.0001, 0,       minX,  0.0001, 0     ];
            ub = [ maxX,  100,    0.4,     maxX,  100,    0.4   ];

            % 6) Create an optimization problem for fmincon
            opts = optimoptions('fmincon', ...
                'Display','none', ...        % 'iter' if you want output
                'Algorithm','interior-point',...
                'MaxFunctionEvaluations',1e4,...
                'MaxIterations',1e4);

            problem = createOptimProblem('fmincon', ...
                'objective', objectiveFunction, ...
                'x0',        initGuess, ...
                'lb',        lb, ...
                'ub',        ub, ...
                'options',   opts);

            % 7) Use MultiStart for more robust fits
            ms = MultiStart('Display','off','UseParallel',false);
            nStarts = 20; % increase if needed
            [bestParams, fval, exitflag, output] = run(ms, problem, nStarts);

            % 8) Package results into the BayesModel struct
            BayesModel.params           = bestParams;
            BayesModel.negLogLik        = fval;
            BayesModel.exitflag         = exitflag;
            BayesModel.optimizerOutput  = output;
            BayesModel.modelType        = modelType;
            BayesModel.paramLabels      = {'muB','sigmaB','lapseB','muO','sigmaO','lapseO'};

        % =================================================================
        % Add other model cases here if needed (Weibull, logistic, etc.)
        % =================================================================
        otherwise
            error('Unknown modelType: %s', modelType);
    end

end


%% ========================================================================
% Subfunction: Negative Log-Likelihood
%% ========================================================================
function nll = negLogLik(params, ...
    xBaseline, successBaseline, sumBaseline, ...
    xOpto, successOpto, sumOpto, ...
    mdlBaseline, mdlOpto, epsilon)

    % Baseline predicted probability
    pBaseline = mdlBaseline(xBaseline, params);
    % clip to avoid log(0)
    pBaseline = max(epsilon, min(1-epsilon, pBaseline));

    % Opto predicted probability
    pOpto = mdlOpto(xOpto, params);
    pOpto = max(epsilon, min(1-epsilon, pOpto));

    % Log-likelihoods
    llBaseline = successBaseline .* log(pBaseline) + ...
                 (sumBaseline - successBaseline) .* log(1 - pBaseline);

    llOpto = successOpto .* log(pOpto) + ...
             (sumOpto - successOpto) .* log(1 - pOpto);

    % Negative log-likelihood
    nll = -sum(llBaseline) - sum(llOpto);

end


%% ========================================================================
% Subfunction: Heuristic Initial Guesses for Normal Model
%% ========================================================================
function initParams = getInitialGuessNormal( ...
    xBaseline, successBaseline, sumBaseline, ...
    xOpto, successOpto, sumOpto )

    % -------------------------------------------------------
    % Estimate approximate threshold & slope for Baseline
    % -------------------------------------------------------
    pBaseline = successBaseline ./ sumBaseline;
    [muB, sigmaB] = estimateThresholdSlope(xBaseline, pBaseline);

    % -------------------------------------------------------
    % Estimate approximate threshold & slope for Opto
    % -------------------------------------------------------
    pOpto = successOpto ./ sumOpto;
    [muO, sigmaO] = estimateThresholdSlope(xOpto, pOpto);

    % A small default lapse
    lapseB = 0.02;
    lapseO = 0.02;

    % params = [muB, sigmaB, lapseB, muO, sigmaO, lapseO]
    initParams = [ muB, sigmaB, lapseB, muO, sigmaO, lapseO ];

end


%% ========================================================================
% Subfunction: estimateThresholdSlope
% Heuristic approach: 
%  - threshold ~ x-value near 50% correct
%  - slope ~ x-distance between 25% and 75% correct
%% ========================================================================
function [approxThresh, approxSlope] = estimateThresholdSlope(xVals, pVals)

    % Sort data by x
    [xSorted, iSort] = sort(xVals);
    pSorted = pVals(iSort);

    % Approx threshold = first x where p >= 0.5
    halfIdx = find(pSorted >= 0.5, 1, 'first');
    if isempty(halfIdx)
        % fallback: median x
        approxThresh = median(xVals);
    else
        approxThresh = xSorted(halfIdx);
    end

    % Approx slope = difference in x between 0.25 and 0.75 correct
    pct25idx = find(pSorted >= 0.25, 1, 'first');
    pct75idx = find(pSorted >= 0.75, 1, 'first');

    if ~isempty(pct25idx) && ~isempty(pct75idx)
        x25 = xSorted(pct25idx);
        x75 = xSorted(pct75idx);
        approxSlope = max( (x75 - x25), 0.0001 );
    else
        approxSlope = 1; % fallback if data doesn't span 25%-75%
    end

end
