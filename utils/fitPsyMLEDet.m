function [mdl, mdlAvg] = fitPsyMLEDet(xBlocks, yBlocks, modelType, solverOption)
% fitPsyMLE fits the psychometric model to the data.
%
% This version fixes the random seed for reproducibility and uses a modular 
% optimization subfunction (optimizeParametersDet) that lets you choose between
% multiple restarts using fmincon or a global search approach.
%
% INPUTS:
%   xBlocks, yBlocks : Data arrays (as in your original code)
%   modelType        : A string specifying the model type (e.g., 'bill')
%   solverOption     : (Optional) 'multipleRestarts' (default) or 'globalSearch'
%
% OUTPUTS:
%   mdl   : Structure containing fitted model results
%   mdlAvg: Averaged model structure (as in your original code)
%
% Example:
%   [mdl, mdlAvg] = fitPsyMLE(xBlocks, yBlocks, 'bill', 'multipleRestarts');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Fix the Random Seed for Reproducibility
    rng(42, 'twister');

    if nargin < 4 || isempty(solverOption)
        solverOption = 'multipleRestarts';
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Initialization and Model Configuration
    [nConditions, nContrasts, nBlocks] = size(xBlocks);
    mdl = struct(); mdlAvg = struct();
    manualFittingFlag = 0;
    
    config = mdlConfig();
    modelConfig = config.models.(modelType);
    
    % Process all blocks data (using your helper function)
    [xBaselineAll, yBaselineAll, xOptoAll, yOptoAll] = processConditionsBlocks(xBlocks, yBlocks, modelType);
    
    % Get averaged data (if needed)
    xBaselineAverage = rmnan(reshape(xBaselineAll, 1, numel(xBaselineAll)));
    yBaselineAverage = rmnan(reshape(yBaselineAll, 1, numel(yBaselineAll)));
    xOptoAverage = rmnan(reshape(xOptoAll, 1, numel(xOptoAll)));
    yOptoAverage = rmnan(reshape(yOptoAll, 1, numel(yOptoAll)));
    
    fprintf('Fitting model, independent parameters per session...\n');

    % Get model functions (the likelihood and model prediction)
    modelFunc = modelConfig.getModelFuncs;  % Get the function handle
    [mdlStruct, objectiveFunction] = modelFunc(config.common);  % Call it
    
    % Setup optimization options (using fminsearch options in the original code)
    % (These options can be changed or replaced with fmincon options in the optimizer.)
    opts = optimset('fminsearch');
    opts = optimset(opts, 'TolX', 1e-4);
    opts = optimset(opts, 'TolFun', 1e-4);
    opts = optimset(opts, 'Display', 'iter');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Loop over blocks and Fit Model
    for block = 1:nBlocks
        % Get block data
        xBaseline = rmnan(xBaselineAll(block,:));
        yBaseline = rmnan(yBaselineAll(block,:));

        % Process opto data: split into incongruent and congruent halves
        numOpto = numel(rmnan(xOptoAll(block,:)));
        inconIdx = 1:floor(numOpto/2);
        conIdx = floor(numOpto/2)+1 : numOpto;
        xInconOpto = -rmnan(xOptoAll(block,inconIdx));
        yInconOpto = -rmnan(yOptoAll(block,inconIdx));
        xConOpto   = rmnan(xOptoAll(block,conIdx));
        yConOpto   = rmnan(yOptoAll(block,conIdx));

        % Get model parameters (initial parameters and bounds)
        initParamFunc = modelConfig.getInitParams;  
        [initialParams, lb, ub] = initParamFunc();
        switch modelType
            case 'bill'
                mdl.headers = {'a','b','l','w','g0','n','rmx','e','o', 'AICc'};
        end        
        lb = max(lb, eps);
        ub(isinf(ub)) = 1e10;
        
        % Get trial counts from the data
        [sumBaseline, successBaseline] = convert2counts(xBaseline, yBaseline);
        [sumInconOpto, successInconOpto] = convert2counts(xInconOpto, yInconOpto);
        [sumConOpto, successConOpto] = convert2counts(xConOpto, yConOpto);

        % Package data for the objective function
        data = struct(...
            'xBaseline', xBaseline, ...
            'yBaseline', yBaseline, ...
            'xInconOpto', xInconOpto, ...
            'yInconOpto', yInconOpto, ...
            'xConOpto', xConOpto, ...
            'yConOpto', yConOpto, ...
            'sumBaseline', sumBaseline, ...
            'successBaseline', successBaseline, ...
            'sumInconOpto', sumInconOpto, ...
            'successInconOpto', successInconOpto, ...
            'sumConOpto', sumConOpto, ...
            'successConOpto', successConOpto);
        
        % Fit the model using either manual fit or automated optimization
        if manualFittingFlag
            fittedParams = initialParams;
            nLL = objectiveFunction(fittedParams, data);
        else
            % Use the modular optimization subfunction
            nRestarts = 5;  % You can adjust the number of restarts here
            [fittedParams, nLL] = optimizeParametersDet(objectiveFunction, initialParams, lb, ub, opts, data, solverOption, nRestarts);
        end
        
        % Calculate metrics (AICc)
        n = sum([sumBaseline sumInconOpto]);
        k = sum(~isnan(fittedParams));
        [~, aicc, ~] = calculateAIC(nLL, k, n);
        
        % Save results (using your helper function)
        mdl = saveModelResults(mdl, block, xBaseline, yBaseline, xInconOpto, yInconOpto, xConOpto, yConOpto, ...
            mdlStruct, objectiveFunction, ub, lb, initialParams, ...
            opts, fittedParams, aicc, manualFittingFlag);
        
        fprintf('Fitted block #%d\n', block)
    end
    fprintf('Done!\n\n');
end

%% Modular Optimization Subfunction
function [bestParams, bestNegLogLik] = optimizeParametersDet(objectiveFunction, initialParams, lb, ub, opts, data, solverOption, nRestarts)
% This function selects the optimization method based on solverOption.
% Allowed options:
%   'multipleRestarts' - Runs fmincon from several random starting points.
%   'globalSearch'     - Uses GlobalSearch with fmincon.
%
    if strcmpi(solverOption, 'multipleRestarts')
        bestNegLogLik = Inf;
        bestParams = [];
        paramDim = numel(lb);
        % Set fmincon options for fmincon (you can modify these if needed)
        fminconOpts = optimoptions('fmincon',...
            'Algorithm','interior-point',...
            'Display','iter',...
            'MaxIterations',1000,...
            'MaxFunctionEvaluations',5000,...
            'TolFun',1e-4,...
            'TolX',1e-4);
        for i = 1:nRestarts
            % Generate a random starting point within bounds
            p0 = lb + rand(1, paramDim).*(ub - lb);
            [pFit, fVal] = fmincon(@(p) objectiveFunction(p, data), p0, [], [], [], [], lb, ub, [], fminconOpts);
            if fVal < bestNegLogLik
                bestNegLogLik = fVal;
                bestParams = pFit;
            end
        end
    elseif strcmpi(solverOption, 'globalSearch')
        fminconOpts = optimoptions('fmincon',...
            'Algorithm','interior-point',...
            'Display','iter',...
            'MaxIterations',1000,...
            'MaxFunctionEvaluations',5000,...
            'TolFun',1e-4,...
            'TolX',1e-4);
        problem = createOptimProblem('fmincon',...
            'x0', initialParams, ...
            'objective', @(p) objectiveFunction(p, data), ...
            'lb', lb, ...
            'ub', ub, ...
            'options', fminconOpts);
        gs = GlobalSearch('Display','iter');
        [bestParams, bestNegLogLik] = run(gs, problem);
    else
        error('Unknown solver option: %s. Use "multipleRestarts" or "globalSearch".', solverOption);
    end
end

%% --- Helper Functions (as in your original file) ---

function [fittedParams, nLL] = fitParameters(objectiveFunction, params, lb, ub, opts, data)
    % (Not used in this modified version; kept here for reference.)
    objFuncWrapped = @(p) objectiveFunction(p, data);
    fittedParams = fminsearchbnd(objFuncWrapped, params, lb, ub, opts);
    nLL = objFuncWrapped(fittedParams);
end

function [aic, aicc, bic] = calculateAIC(nLL, k, n)
    aic = 2 * k + 2 * nLL;
    if n > k + 1
        aicc = aic + (2 * k * (k + 1)) / (n - k - 1);
    else
        aicc = Inf;
    end
    bic = log(n) * k + 2 * nLL;
end

function mdl = saveModelResults(mdl, block, xBaseline, yBaseline, xInconOpto, yInconOpto, xConOpto, yConOpto, ...
    mdlStruct, objectiveFunction, ub, lb, params, options, fittedParams, aicc, manualFittingFlag)
    mdl.xBaseline(block,:) = padArray(xBaseline, 12, 2, NaN);
    mdl.yBaseline(block,:) = padArray(yBaseline, 12, 2, NaN);
    mdl.xInconOpto(block,:) = padArray(xInconOpto, 12, 2, NaN);
    mdl.yInconOpto(block,:) = padArray(yInconOpto, 12, 2, NaN);
    mdl.xConOpto(block,:) = padArray(xConOpto, 12, 2, NaN);
    mdl.yConOpto(block,:) = padArray(yConOpto, 12, 2, NaN);
    mdl.mdlBaseline = mdlStruct;
    mdl.mdlOpto = mdlStruct;
    mdl.objectiveFunction = objectiveFunction;
    mdl.ub = ub;
    mdl.lb = lb;
    mdl.params(block,:) = params;
    mdl.options = options;
    mdl.fittedParams(block,:,manualFittingFlag+1) = [fittedParams aicc];
end

function padded = padArray(arr, desiredLength, dim, padValue)
    if length(arr) < desiredLength
        padded = [arr, repmat(padValue, 1, desiredLength - length(arr))];
    else
        padded = arr(1:desiredLength);
    end
end

function out = rmnan(arr)
    out = arr(~isnan(arr));
end

function [xBaselineAll, yBaselineAll, xOptoAll, yOptoAll] = processConditionsBlocks(xBlocks, yBlocks, modelType)
    nBlocks = size(xBlocks, 3);
    maxBaseLen = 0;
    maxOptoLen = 0;
    for block = 1:nBlocks
        xBaseline = rmnan(squeeze(xBlocks(1,:,block)));
        xHorizontal = rmnan(squeeze(xBlocks(2,:,block)));
        xVertical = rmnan(squeeze(xBlocks(3,:,block)));
        maxBaseLen = max(maxBaseLen, numel(xBaseline));
        maxOptoLen = max(maxOptoLen, numel(xHorizontal) + numel(xVertical));
    end
    xBaselineAll = NaN(nBlocks, maxBaseLen);
    yBaselineAll = NaN(nBlocks, maxBaseLen);
    xOptoAll = NaN(nBlocks, maxOptoLen);
    yOptoAll = NaN(nBlocks, maxOptoLen);
    
    if strcmp(modelType, 'bill')
        for block = 1:nBlocks
            [xBaseline, sortIdx] = sort(rmnan(squeeze(xBlocks(1,:,block))));
            yBaseline = rmnan(squeeze(yBlocks(1,:,block)));
            yBaseline = yBaseline(sortIdx);
            numVal = numel(xBaseline);
            xBaseline = mean([fliplr(-xBaseline(1:numVal/2)); xBaseline(numVal/2+1:end)]);
            yBaseline = mean([fliplr(100 - yBaseline(1:numVal/2)); yBaseline(numVal/2+1:end)]);
            xBaselineAll(block, 1:numel(xBaseline)) = xBaseline;
            yBaselineAll(block, 1:numel(yBaseline)) = yBaseline;
            
            [xHorizontal, sortIdxH] = sort(rmnan(squeeze(xBlocks(2,:,block))));
            [xVertical, sortIdxV] = sort(rmnan(squeeze(xBlocks(3,:,block))));
            contrastNeg = 1:floor(numel(xHorizontal)/2);
            contrastPos = floor(numel(xHorizontal)/2)+1:numel(xHorizontal);
            
            yHorizontal = rmnan(squeeze(yBlocks(2,:,block)));
            yHorizontal = yHorizontal(sortIdxH);
            yHorizontal(contrastNeg) = 100 - yHorizontal(contrastNeg);
            yVertical = rmnan(squeeze(yBlocks(3,:,block)));
            yVertical = yVertical(sortIdxV);
            yVertical(contrastNeg) = 100 - yVertical(contrastNeg);
            
            xConOpto = mean([-fliplr(xHorizontal(contrastNeg)); xVertical(contrastPos)]);
            yConOpto = mean([fliplr(yHorizontal(contrastNeg)); yVertical(contrastPos)]);
            xInconOpto = mean([xHorizontal(contrastPos); -fliplr(xVertical(contrastNeg))]);
            yInconOpto = mean([yHorizontal(contrastPos); fliplr(yVertical(contrastNeg))]);
            
            optoLen = numel([-xInconOpto, xConOpto]);
            xOptoAll(block, 1:optoLen) = [-xInconOpto, xConOpto];
            yOptoAll(block, 1:optoLen) = [-yInconOpto, yConOpto];
        end
    else
        for block = 1:nBlocks
            [xBaseline, sortIdx] = sort(rmnan(squeeze(xBlocks(1,:,block))));
            yBaseline = rmnan(squeeze(yBlocks(1,:,block)));
            yBaseline = yBaseline(sortIdx);
            numVal = numel(xBaseline);
            yBaselineAvg = mean([fliplr(100 - yBaseline(1:numVal/2)); yBaseline(numVal/2+1:end)]);
            yBaseline = [100 - fliplr(yBaselineAvg), yBaselineAvg];
            [xBaseline, yBaseline] = mergeZeros(xBaseline, yBaseline);
            xBaselineAll(block, 1:numel(xBaseline)) = xBaseline;
            yBaselineAll(block, 1:numel(yBaseline)) = yBaseline;
            
            [xHorizontal, sortIdx] = sort(-1 * rmnan(squeeze(xBlocks(2,:,block))));
            yHorizontal = rmnan(squeeze(yBlocks(2,:,block)));
            yHorizontal = 100 - yHorizontal(sortIdx);
            
            [xVertical, sortIdx] = sort(rmnan(squeeze(xBlocks(3,:,block))));
            yVertical = rmnan(squeeze(yBlocks(3,:,block)));
            yVertical = yVertical(sortIdx);
            
            xOpto = mean([xVertical; xHorizontal]);
            yOpto = mean([yVertical; yHorizontal]);
            [xOpto, yOpto] = mergeZeros(xOpto, yOpto);
            
            optoLen = numel(xOpto);
            xOptoAll(block, 1:optoLen) = xOpto;
            yOptoAll(block, 1:optoLen) = yOpto;
        end
    end
end

function [xProcessed, yProcessed] = mergeZeros(x, y)
    zeroIndices = find(x == 0);
    switch length(zeroIndices)
        case 0
            xProcessed = x;
            yProcessed = y;
        case 2
            xProcessed = x;
            xProcessed(zeroIndices(2)) = [];
            yProcessed = y;
            avgValue = mean(y(zeroIndices));
            yProcessed(zeroIndices) = avgValue;
            yProcessed(zeroIndices(2)) = [];
        otherwise
            error('Unexpected number of zeros. Expected 0 or 2 zeros.');
    end
end

%% Dummy mdlConfig Function (as in your original code)
function config = mdlConfig()
    config.common = struct(...
        'minBound', 1e-10, ...
        'maxBound', 100 - 1e-10, ...
        'betaBL', 50, ...
        'lowerAsymptoteBL', 0.5, ...
        'epsilon', 1e-8 );
    config.models = struct();
    config.models.bill = struct(...
        'getInitParams', @getBillInitParams, ...
        'getModelFuncs', @getBillModelFuncs, ...
        'headers', {'a','b','l','w','g0','n','rmx','e','o','AICc'} );
end

%% Dummy getBillInitParams Function (as in your original code)
function [initialParams, lb, ub] = getBillInitParams()
    initialParams = [0.1, 0.13, 0, 0, 60, 3.5, 100, 0.4, 20];
    lb = [0, 0, 0, 0, 30, 1, 70, 0, 0];
    ub = [0.3, 0.3, 0.3, 0.3, 100, 6, 100, 1, 30];
end

%% Dummy getBillModelFuncs Function (as in your original code)
function [mdlStruct, objectiveFunction] = getBillModelFuncs(common)
    options = struct('nTrials', 20000);
    mdlStruct = [];  % Placeholder for your model structure
    sumAll = @(x) sum(x(:));
    epsilon = 1e-10;
    objectiveFunction = @(params, data) ...
        - sumAll(data.successBaseline .* log(max(columnarBayesMdl(data.xBaseline, params).pcntrl/100, epsilon)) + ...
                 (data.sumBaseline - data.successBaseline) .* log(max(1 - columnarBayesMdl(data.xBaseline, params).pcntrl/100, epsilon))) ...
        - sumAll(data.successInconOpto .* log(max(columnarBayesMdl(data.xInconOpto, params).pic/100, epsilon)) + ...
                 (data.sumInconOpto - data.successInconOpto) .* log(max(1 - columnarBayesMdl(data.xInconOpto, params).pic/100, epsilon))) ...
        - sumAll(data.successConOpto .* log(max(columnarBayesMdl(data.xConOpto, params).pc/100, epsilon)) + ...
                 (data.sumConOpto - data.successConOpto) .* log(max(1 - columnarBayesMdl(data.xConOpto, params).pc/100, epsilon)));
end

%% Dummy Model Function for Simulation (as in your original code)
function yPredicted = columnarBayesMdl(x, params)
    % This is a placeholder for your simulation-based model.
    % Replace with your actual implementation.
    a = params(1); b = params(2); l = params(3); w = params(4);
    g0 = params(5); n = params(6); rmx = params(7);
    e = params(8); o_stim = params(9);
    
    % Example: compute a simple sigmoid as percent correct
    pc = 1./(1+exp(-(x-50)/10));
    pic = pc;
    yPredicted.pcntrl = pc*100;
    yPredicted.pc = pc*100;
    yPredicted.pic = pic*100;
end
