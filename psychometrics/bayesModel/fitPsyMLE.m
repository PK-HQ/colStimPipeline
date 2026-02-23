function [mdl, mdlAvg] = fitPsyMLE(xBlocks, yBlocks, modelType, plotLine)
    % Initialize
    [nConditions, nContrasts, nBlocks] = size(xBlocks);
    mdl = struct(); mdlAvg = struct();
    manualFittingFlag = 0;
    
    % Get model configuration
    config = mdlConfig();
    modelConfig = config.models.(modelType);
    
    % Process all blocks data
    [xBaselineAll, yBaselineAll, xOptoAll, yOptoAll] = processConditionsBlocks(xBlocks, yBlocks, modelType);
    
    % Get averaged data
    xBaselineAverage = rmnan(reshape(xBaselineAll, 1, numel(xBaselineAll)));
    yBaselineAverage = rmnan(reshape(yBaselineAll, 1, numel(yBaselineAll)));
    xOptoAverage = rmnan(reshape(xOptoAll, 1, numel(xOptoAll)));
    yOptoAverage = rmnan(reshape(yOptoAll, 1, numel(yOptoAll)));
    
    fprintf('Fitting model, independent parameters per session...')
    for block = 1:nBlocks
        % Get block data
        xBaseline = rmnan(xBaselineAll(block,:));
        yBaseline = rmnan(yBaselineAll(block,:));

        inconIdx=1:numel(rmnan(xOptoAll(block,:)))/2;
        conIdx=numel(rmnan(xOptoAll(block,:)))/2 +1 : numel(rmnan(xOptoAll(block,:)));
        xInconOpto = -rmnan(xOptoAll(block,inconIdx));
        yInconOpto = -rmnan(yOptoAll(block,inconIdx));
        xConOpto = rmnan(xOptoAll(block,conIdx));
        yConOpto = rmnan(yOptoAll(block,conIdx));

        % Get model parameters
        initParamFunc = modelConfig.getInitParams;  % Get the function handle
        [initialParams, lb, ub] = initParamFunc();  % Call it
        switch modelType
            case 'bill'
                mdl.headers={'a','b','l','w','g0','n','rmx','e','o', 'AICc'};
            case 'weibullfreeAll'
                mdl.headers={'A^{bl}','B^{bl}','\alpha^{bl}','\beta^{bl}', ...
                   '\DeltaA^{con-bl}','\DeltaB^{con-bl}','\Delta\alpha^{con-bl}','\Delta\beta^{con-bl}', ...
                   '\DeltaA^{incon-bl}','\DeltaB^{incon-bl}','\Delta\alpha^{incon-bl}','\Delta\beta^{incon-bl}', ...
                   'AICc^{total}','AUC^{bl}','AUC^{con-bl}','AUC^{incon-bl}','AUC^{con-incon}'};
        end        
        % Adjust bounds
        lb = max(lb, eps);
        ub(isinf(ub)) = 1e10;
        
        % Setup optimization options
        if plotLine==1
            maxIterVal=1000;
            maxFunEvalsVal=1000;
        else
            maxIterVal=1;
            maxFunEvalsVal=1;
        end
        opts = optimset('fminsearch');
        opts = optimset(opts, 'MaxIter', maxIterVal);
        opts = optimset(opts, 'MaxFunEvals', maxFunEvalsVal);
        opts = optimset(opts, 'TolX', 1e-6);
        opts = optimset(opts, 'TolFun', 1e-6);
        opts = optimset(opts, 'Display', 'iter');

        % Get model functions
        modelFunc = modelConfig.getModelFuncs;  % Get the function handle
        [mdlStruct, objectiveFunction] = modelFunc(config.models.(modelType));  % Call it
        
        % Get trial counts
        [sumBaseline, successBaseline] = convert2counts(xBaseline, yBaseline);
        [sumInconOpto, successInconOpto] = convert2counts(xInconOpto, yInconOpto);
        [sumConOpto, successConOpto] = convert2counts(xConOpto, yConOpto);

        % Package data for objective function
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
            'successInconOpto', successInconOpto,...
            'sumConOpto', sumConOpto, ...
            'successConOpto', successConOpto);
        
        % Fit model
        if manualFittingFlag
            fittedParams = initialParams;
            nLL = objectiveFunction(fittedParams, data);
        else
            rng(42 + block, 'twister');   % different PSO swarm per block
            solverOption = 'globalSearch';
            % Use the modular optimization subfunction
            %nRestarts = 5;  % You can adjust the number of restarts here
            %[fittedParams, nLL] = globalSolvers(objectiveFunction, initialParams, lb, ub, opts, data, solverOption, nRestarts);            
            [fittedParams, nLL] = fitParameters(objectiveFunction, initialParams, lb, ub, opts, data);
        end
        fittedParams
        
        
        % Calculate metrics
        n = sum([sumBaseline sumInconOpto]);
        k = sum(~isnan(fittedParams));
        [~, aicc, ~] = calculateAIC(nLL, k, n);
        
        % Save results
        mdl = saveModelResults(mdl, block, xBaseline, yBaseline, xInconOpto, yInconOpto, xConOpto, yConOpto, ...
            mdlStruct, objectiveFunction, ub, lb, initialParams, ...
            opts, fittedParams, aicc, manualFittingFlag);
        
        fprintf('Fitted block #%.0f', block)
    end
    fprintf('Done!\n\n')
    fprintf('\n=========\n\n')
end

%% Helper Functions
function [fittedParams, nLL] = fitParameters(objectiveFunction, params0, lb, ub, opts, data)

    % ----- Unit-cube mapping -----
    toParams = @(u) lb + u(:)'.*(ub - lb);
    toUnit   = @(p) (p(:)' - lb)./(ub - lb);

    % Clamp initial
    u0 = toUnit(params0);
    u0 = min(max(u0, 0), 1);

    havePSO = exist('particleswarm','file') == 2;
    haveFMC = exist('fmincon','file') == 2;
    havePAT = exist('patternsearch','file') == 2;

    % ----- Stage budgets (tune these) -----
    % Stage 1: cheap baseline shaping
    budget1 = struct('trialScale', 10, 'minTrialsPerLevel', 300, 'useGPU', true, 'seed', 1, 'balanceTrialTypes', true);
    % Stage 2: medium for opto parameters
    budget2 = struct('trialScale', 25, 'minTrialsPerLevel', 500, 'useGPU', true, 'seed', 2, 'balanceTrialTypes', true);
    % Stage 3: expensive final refinement
    budget3 = struct('trialScale', 25, 'minTrialsPerLevel', 500,'useGPU', true, 'seed', 3, 'balanceTrialTypes', true);

    % Optional tiny bound barrier (off by default)
    data.useBoundBarrier = false;

    % ----- Param indices in your order: [a b l w g0 n rmx e o] -----
    idx_base = [1 3 5 6 7];      % a, l, g0, n, rmx
    idx_opto = [2 9];        % b, w, e, o

    % ===== Stage 1: baseline-only local fit (CPU/GPU deterministic) =====
    data_base = data;
    data_base.sumConOpto(:) = 0;    data_base.successConOpto(:) = 0;
    data_base.sumInconOpto(:)= 0;   data_base.successInconOpto(:)= 0;
    data_base.simBudget = budget1;

    u1 = u0;
    obj1 = @(uSub) objectiveFunction(toParams(localSetSub(u1, idx_base, uSub)), data_base);

    uSub0 = u1(idx_base);
    [uSub1, ~] = fminsearchbnd(obj1, uSub0, zeros(size(uSub0)), ones(size(uSub0)), opts);
    u1 = localSetSub(u1, idx_base, uSub1);

    % ===== Stage 2: opto-only PSO on [b,w,e,o] (seeded swarm + hybrid refine) =====
    data_opto = data;
    data_opto.sumBaseline(:) = 0;   data_opto.successBaseline(:) = 0;
    data_opto.simBudget = budget2;

    u2 = u1;
    % ---- FIX shared w,e across blocks ----
    % These should be passed in (shared), but if not, fall back to initial values
    if isfield(data,'w_shared') && ~isempty(data.w_shared)
        w_fix = data.w_shared;
    else
        w_fix = params0(4);
    end
    
    if isfield(data,'e_shared') && ~isempty(data.e_shared)
        e_fix = data.e_shared;
    else
        e_fix = params0(8);
    end
    
    p2tmp = toParams(u2);
    p2tmp(4) = w_fix;   % w fixed
    p2tmp(8) = e_fix;   % e fixed
    u2 = min(max(toUnit(p2tmp), 0), 1);

    obj2 = @(uSub) objectiveFunction(toParams(localSetSub(u2, idx_opto, uSub)), data_opto);

    uSub0 = u2(idx_opto);

    % Seeded swarm around current params in unit space
    swarmSize = 150;     % 150–300
    maxIters  = 150;     % 150–300
    sigma     = 0.08;    % jitter radius in unit space (tune 0.05–0.15)

    initSwarm = min(max(uSub0 + sigma*randn(swarmSize, numel(uSub0)), 0), 1);
    initSwarm(1,:) = uSub0; % include exact current point

    if havePSO
        % Hybrid local solver for Stage 2
        hybrid = [];
        if haveFMC
            fopts = optimoptions('fmincon', ...
                'Algorithm','interior-point', ...
                'Display','off', ...
                'MaxFunctionEvaluations', 2e4, ...
                'FiniteDifferenceType','central');
            hybrid = {@fmincon, fopts};
        elseif havePAT
            popts2 = optimoptions('patternsearch', ...
                'Display','off', ...
                'UseCompletePoll', true);
            hybrid = {@patternsearch, popts2};
        end

        psopts = optimoptions('particleswarm', ...
            'SwarmSize', swarmSize, ...
            'MaxIterations', maxIters, ...
            'UseParallel', false, ...
            'Display', 'iter', ...
            'InitialSwarmMatrix', initSwarm);

        if ~isempty(hybrid)
            psopts = optimoptions(psopts, 'HybridFcn', hybrid);
        end

        [uSub2, ~] = particleswarm(obj2, numel(uSub0), zeros(size(uSub0)), ones(size(uSub0)), psopts);
    else
        % Fallback if PSO not available
        [uSub2, ~] = fminsearchbnd(obj2, uSub0, zeros(size(uSub0)), ones(size(uSub0)), opts);
    end

    u2 = localSetSub(u2, idx_opto, uSub2);

    % ===== Stage 3: joint local refinement (expensive budget) =====
    data_all = data;
    data_all.simBudget = budget3;

    obj3 = @(u) objectiveFunction(toParams(u), data_all);

    if haveFMC
        fopts = optimoptions('fmincon', ...
            'Algorithm','interior-point', ...
            'Display','iter', ...
            'MaxFunctionEvaluations', 5e4, ...
            'FiniteDifferenceType','central');
        u3 = fmincon(obj3, u2, [],[],[],[], zeros(size(u2)), ones(size(u2)), [], fopts);
    else
        u3 = fminsearchbnd(obj3, u2, zeros(size(u2)), ones(size(u2)), opts);
    end

    fittedParams = toParams(u3);
    nLL = objectiveFunction(fittedParams, data_all);
end

function u = localSetSub(u, idx, uSub)
    u(idx) = uSub;
end

function mdl = saveModelResults(mdl, block, xBaseline, yBaseline, xInconOpto, yInconOpto, xConOpto, yConOpto, ...
    mdlStruct, objectiveFunction, ub, lb, params, options, ...
    fittedParams, aicc, manualFittingFlag)

    
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

function [sumY, successY] = convert2counts(x,y)
    x = rmnan(x);
    y = rmnan(y);

    % Experiment design:
    %   40 trials at 0% contrast
    %   20 trials at all other contrasts
    sumY = 20 * ones(size(x));
    sumY(x==0) = 40;

    successY = (y / 100) .* sumY;
end


function [aic, aicc, bic] = calculateAIC(nLL, k, n)
    aic = 2 * k + 2 * nLL;
    if n > k + 1
        aicc = aic + (2 * k * (k + 1)) / (n - k - 1);
    else
        aicc = Inf;
    end
    bic = log(n) * k + 2 * nLL;
    [aic, aicc, bic];
end

function [xBaselineAll, yBaselineAll, xOptoAll, yOptoAll] = processConditionsBlocks(xBlocks, yBlocks, modelType)
    % Arrays
    % xBlock/yBlock (3,:,blockNo): Baseline; Horizontal opto; Vertical opto
    nBlocks = size(xBlocks, 3);
    maxBaseLen = 0;
    maxOptoLen = 0;

    % First pass: determine maximum lengths needed for preallocation
    for block = 1:nBlocks
        xBaseline = rmnan(squeeze(xBlocks(1,:,block)));
        xHorizontal = rmnan(squeeze(xBlocks(2,:,block)));
        xVertical = rmnan(squeeze(xBlocks(3,:,block)));

        maxBaseLen = max(maxBaseLen, numel(xBaseline));
        maxOptoLen = max(maxOptoLen, numel(xHorizontal) + numel(xVertical));
    end

    % Initialize arrays with NaNs
    xBaselineAll = NaN(nBlocks, maxBaseLen);
    yBaselineAll = NaN(nBlocks, maxBaseLen);
    xOptoAll = NaN(nBlocks, maxOptoLen);
    yOptoAll = NaN(nBlocks, maxOptoLen);

    switch modelType
        case {'bill','weibullfreeAll'}
            % Second pass: process each block and store results
            for block = 1:nBlocks
                % Process baseline: sort X and Y (% vertical) with same index, mean Y across both 
                % -ve (*-1 to convert to % correct) and +ve X (already is % correct). 
                [xBaseline, sortIdx] = sort(rmnan(squeeze(xBlocks(1,:,block)))); %palindromic -ve to +ve
                yBaseline = rmnan(squeeze(yBlocks(1,:,block)));
                yBaseline = yBaseline(sortIdx);
                numVal = numel(xBaseline);
                xBaseline = mean([fliplr(-xBaseline(1:numVal/2)); xBaseline(numVal/2+1:end)]);
                yBaseline = mean([fliplr(100 - yBaseline(1:numVal/2)); yBaseline(numVal/2+1:end)]);
                % Insert into preallocated arrays
                xBaselineAll(block, 1:numel(xBaseline)) = xBaseline;
                yBaselineAll(block, 1:numel(yBaseline)) = yBaseline;

                % Sort
                [xHorizontal, sortIdxH]=sort(rmnan(squeeze(xBlocks(2,:,block)))); %palindromic -ve to +ve
                [xVertical, sortIdxV] = sort(rmnan(squeeze(xBlocks(3,:,block))));
                contrastNeg = 1:numel(xHorizontal)/2; contrastPos = numel(xHorizontal)/2+1:numel(xHorizontal);

                yHorizontal = rmnan(squeeze(yBlocks(2,:,block))); %get data
                yHorizontal=yHorizontal(sortIdxH); %sort with same idx as x
                yHorizontal(contrastNeg) = 100 - yHorizontal(contrastNeg);
                yVertical = rmnan(squeeze(yBlocks(3,:,block)));
                yVertical = yVertical(sortIdxV);
                yVertical(contrastNeg) = 100 - yVertical(contrastNeg);

                xConOpto=mean([-fliplr(xHorizontal(contrastNeg)); xVertical(contrastPos)]);
                yConOpto=mean([fliplr(yHorizontal(contrastNeg)); yVertical(contrastPos)]);

                xInconOpto=mean([xHorizontal(contrastPos); -fliplr(xVertical(contrastNeg))]);
                yInconOpto=mean([yHorizontal(contrastPos); fliplr(yVertical(contrastNeg))]);

                % Insert into preallocated arrays
                optoLen = numel([-xInconOpto, xConOpto]);
                xOptoAll(block, 1:optoLen) = [-xInconOpto, xConOpto];
                yOptoAll(block, 1:optoLen) = [-yInconOpto, yConOpto];
            end
        otherwise
            % Second pass: process each block and store results
            for block = 1:nBlocks
                % Process baseline: sort X and Y (% vertical) with same index, mean Y across both -ve and
                % +ve X, convert back to % vertical.
                [xBaseline, sortIdx] = sort(rmnan(squeeze(xBlocks(1,:,block))));
                yBaseline = rmnan(squeeze(yBlocks(1,:,block)));
                yBaseline = yBaseline(sortIdx);
                numVal = numel(xBaseline);
                yBaselineAvg = mean([fliplr(100 - yBaseline(1:numVal/2)); yBaseline(numVal/2+1:end)]);
                yBaseline = [100 - fliplr(yBaselineAvg), yBaselineAvg];
                [xBaseline, yBaseline] = mergeZeros(xBaseline, yBaseline);
        
                % Insert into preallocated arrays
                xBaselineAll(block, 1:numel(xBaseline)) = xBaseline;
                yBaselineAll(block, 1:numel(yBaseline)) = yBaseline;
        
                % Process horizontal
                [xHorizontal, sortIdx] = sort(-1 * rmnan(squeeze(xBlocks(2,:,block))));
                yHorizontal = rmnan(squeeze(yBlocks(2,:,block)));
                yHorizontal = 100 - yHorizontal(sortIdx);
        
                % Process vertical
                [xVertical, sortIdx] = sort(rmnan(squeeze(xBlocks(3,:,block))));
                yVertical = rmnan(squeeze(yBlocks(3,:,block)));
                yVertical = yVertical(sortIdx);
        
                % Calculate mean and merge for Opto
                xOpto = mean([xVertical; xHorizontal]);
                yOpto = mean([yVertical; yHorizontal]);
                [xOpto, yOpto] = mergeZeros(xOpto, yOpto);
        
                % Insert into preallocated arrays
                optoLen = numel(xOpto);
                xOptoAll(block, 1:optoLen) = xOpto;
                yOptoAll(block, 1:optoLen) = yOpto;
            end
    end
end

function [xProcessed, yProcessed] = mergeZeros(x, y)
    % Subfunction to merge zeros and process y based on zero positions in x
    zeroIndices = find(x == 0); % Find indices of zeros in x
    
    switch length(zeroIndices)
        case 0
            % If no zeros are found, return the arrays as they are
            xProcessed = x;
            yProcessed = y;
            
        case 2
            % If two zeros are found, process the arrays
            % Removing one zero from x
            xProcessed = x;
            xProcessed(zeroIndices(2)) = []; % Remove the second zero
            
            % Averaging corresponding elements in y and removing one
            yProcessed = y;
            avgValue = mean(y(zeroIndices));
            yProcessed(zeroIndices) = avgValue; % Replace both indices with their average
            yProcessed(zeroIndices(2)) = []; % Remove the second element
            
        otherwise
            error('Unexpected number of zeros. Expected 0 or 2 zeros.');
    end
end

function config = mdlConfig()
    % Define common parameters used across models
    %{
    config.common = struct(...
        'minBound', 1e-10, ...
        'maxBound', 100 - 1e-10, ...
        'betaBL', 50, ...
        'lowerAsymptoteBL', 0.5, ...
        'epsilon', 1e-8 ...
    );
    %}
    
    % Define all model configurations
    config.models = struct();
    
    % Bill Model
    config.models.bill = struct(...
        'getInitParams', @getBillInitParams, ...
        'getModelFuncs', @getBillModelFuncs, ...
        'headers', {'a','b','l','w','g0','n','rmx','e','o', 'AICc'} ...
    );

    % Weibull Free All Model
    config.models.weibullfreeAll = struct(...
        'getInitParams', @getWeibullInitParams, ...
        'getModelFuncs', @getWeibullModelFuncs, ...
        'headers', {'A^{bl}','B^{bl}','\alpha^{bl}','\beta^{bl}', ...
                   '\DeltaA^{con-bl}','\DeltaB^{con-bl}','\Delta\alpha^{con-bl}','\Delta\beta^{con-bl}', ...
                   '\DeltaA^{incon-bl}','\DeltaB^{incon-bl}','\Delta\alpha^{incon-bl}','\Delta\beta^{incon-bl}', ...
                   'AICc^{total}','AUC^{bl}','AUC^{con-bl}','AUC^{incon-bl}','AUC^{con-incon}'} ...
    );
    
end
