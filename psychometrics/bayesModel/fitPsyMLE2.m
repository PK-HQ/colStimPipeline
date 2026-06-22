function [mdl, mdlAvg] = fitPsyMLE2(xBlocks, yBlocks, modelType, plotLine)
    % Initialize
    [nConditions, nContrasts, nBlocks] = size(xBlocks);
    mdl = struct(); mdlAvg = struct();
    manualFittingFlag = 0;
    
    % Get model configuration
    config = mdlConfig();
    modelConfig = config.models.(modelType);
    
    % Process all blocks data, with premerged
    [xBaselineAll, yBaselineAll, xOptoAll, yOptoAll, ...
         xBaselinePreAll, yBaselinePreAll, tagBaselinePreAll, ...
         xHorizontalOptoPreAll, yHorizontalOptoPreAll, tagHorizontalOptoPreAll, congrHorizontalOptoPreAll, ...
         xVerticalOptoPreAll, yVerticalOptoPreAll, tagVerticalOptoPreAll, congrVerticalOptoPreAll] = processConditionsBlocks(xBlocks, yBlocks, modelType);

    if strcmp(modelType, 'weibullfreeAll')
        mdl.fittedParamsHorizontal = nan(nBlocks, 13);
        mdl.fittedParamsVertical = nan(nBlocks, 13);
        mdl.fitStatusHorizontal = repmat({''}, nBlocks, 1);
        mdl.fitStatusVertical = repmat({''}, nBlocks, 1);
    elseif strcmp(modelType, 'weibullSignedX0')
        mdl.signedX0.fitParams = nan(nBlocks, 10);
        mdl.signedX0.nLL = nan(nBlocks, 1);
        mdl.signedX0.AICc = nan(nBlocks, 1);
        mdl.signedX0.noX0FitParams = nan(nBlocks, 9);
        mdl.signedX0.noX0NLL = nan(nBlocks, 1);
        mdl.signedX0.noX0AICc = nan(nBlocks, 1);
        mdl.signedX0.deltaAICcX0 = nan(nBlocks, 1);
        mdl.signedX0.akaikeWeightX0 = nan(nBlocks, 1);
        mdl.signedX0.akaikeWeightNoX0 = nan(nBlocks, 1);
        mdl.signedX0.X0Horizontal = nan(nBlocks, 1);
        mdl.signedX0.X0Vertical = nan(nBlocks, 1);
        mdl.signedX0.fitStatus = repmat({''}, nBlocks, 1);
    elseif strcmp(modelType, 'weibullSignedBX0')
        mdl.signedBX0.fitParams = nan(nBlocks, 11);
        mdl.signedBX0.nLL = nan(nBlocks, 1);
        mdl.signedBX0.AICc = nan(nBlocks, 1);
        mdl.signedBX0.noX0FitParams = nan(nBlocks, 10);
        mdl.signedBX0.noX0NLL = nan(nBlocks, 1);
        mdl.signedBX0.noX0AICc = nan(nBlocks, 1);
        mdl.signedBX0.deltaAICcX0 = nan(nBlocks, 1);
        mdl.signedBX0.akaikeWeightBX0 = nan(nBlocks, 1);
        mdl.signedBX0.akaikeWeightBOnly = nan(nBlocks, 1);
        mdl.signedBX0.BHorizontal = nan(nBlocks, 1);
        mdl.signedBX0.BVertical = nan(nBlocks, 1);
        mdl.signedBX0.X0Horizontal = nan(nBlocks, 1);
        mdl.signedBX0.X0Vertical = nan(nBlocks, 1);
        mdl.signedBX0.deltaB = nan(nBlocks, 1);
        mdl.signedBX0.deltaX0 = nan(nBlocks, 1);
        mdl.signedBX0.fitStatus = repmat({''}, nBlocks, 1);
    end
    
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

        % pre-merged: signed contrast, percent correct, visual tag, congruency tag
        xBaselinePre = rmnan(xBaselinePreAll(block,:));
        yBaselinePre = rmnan(yBaselinePreAll(block,:));
        tagBaselinePre = rmnan(tagBaselinePreAll(block,:));
        
        xHorizontalOptoPre = rmnan(xHorizontalOptoPreAll(block,:));
        yHorizontalOptoPre = rmnan(yHorizontalOptoPreAll(block,:));
        tagHorizontalOptoPre = rmnan(tagHorizontalOptoPreAll(block,:));
        congrHorizontalOptoPre = rmnan(congrHorizontalOptoPreAll(block,:));
        
        xVerticalOptoPre = rmnan(xVerticalOptoPreAll(block,:));
        yVerticalOptoPre = rmnan(yVerticalOptoPreAll(block,:));
        tagVerticalOptoPre = rmnan(tagVerticalOptoPreAll(block,:));
        congrVerticalOptoPre = rmnan(congrVerticalOptoPreAll(block,:));

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
                mdl.headersHorizontal = mdl.headers;
                mdl.headersVertical = mdl.headers;
            case 'weibullSignedX0'
                mdl.headers={'A^{bl}','\alpha^{bl}','\beta^{bl}', ...
                   'A^{HOpto}','\alpha^{HOpto}','\beta^{HOpto}', ...
                   'A^{VOpto}','\alpha^{VOpto}','\beta^{VOpto}', ...
                   '\DeltaX0','AICc^{X0}'};
            case 'weibullSignedBX0'
                mdl.headers={'A^{bl}','\alpha^{bl}','\beta^{bl}', ...
                   'A^{con}','\alpha^{con}','\beta^{con}', ...
                   'A^{incon}','\alpha^{incon}','\beta^{incon}', ...
                   '\DeltaB','\DeltaX0','AICc^{BX0}'};
        end        
        % Adjust bounds
        if ~ismember(modelType, {'weibullSignedX0', 'weibullSignedBX0'})
            lb = max(lb, eps);
        end
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
        opts = optimset(opts, 'Display', 'off');

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
        if ismember(modelType, {'weibullSignedX0', 'weibullSignedBX0'})
            signedData = buildSignedChoiceData(xBlocks, yBlocks, block);
            data = mergeStructs(data, signedData);
        end
        
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
            if ismember(modelType, {'weibullSignedX0', 'weibullSignedBX0'})
                [fittedParams, nLL] = fitParametersSimple(...
                    objectiveFunction, initialParams, lb, ub, opts, data);
            else
                [fittedParams, nLL] = fitParameters(objectiveFunction, initialParams, lb, ub, opts, data);
            end
        end
        % Calculate metrics
        if ismember(modelType, {'weibullSignedX0', 'weibullSignedBX0'})
            n = sum([data.sumBaselineChoice, ...
                data.sumHorizontalOptoChoice, data.sumVerticalOptoChoice]);
        else
            n = sum([sumBaseline sumInconOpto]);
        end
        k = sum(~isnan(fittedParams));
        [~, aicc, ~] = calculateAIC(nLL, k, n);

        signedFitSummary = [];
        if strcmp(modelType, 'weibullSignedX0')
            signedFitSummary = fitSignedNoX0Comparison(...
                objectiveFunction, initialParams, lb, ub, opts, data, ...
                fittedParams, nLL, aicc, n);
        elseif strcmp(modelType, 'weibullSignedBX0')
            signedFitSummary = fitSignedBX0NoX0Comparison(...
                objectiveFunction, initialParams, lb, ub, opts, data, ...
                fittedParams, nLL, aicc, n);
        end
        
        % Save results
        mdl = saveModelResults(mdl, block, ...
            xBaseline, yBaseline, ...
            xInconOpto, yInconOpto, ...
            xConOpto, yConOpto, ...
            xBaselinePre, yBaselinePre, tagBaselinePre, ...
            xHorizontalOptoPre, yHorizontalOptoPre, tagHorizontalOptoPre, congrHorizontalOptoPre, ...
            xVerticalOptoPre, yVerticalOptoPre, tagVerticalOptoPre, congrVerticalOptoPre, ...
            mdlStruct, objectiveFunction, ub, lb, initialParams, ...
            opts, fittedParams, aicc, manualFittingFlag);

        if strcmp(modelType, 'weibullSignedX0')
            mdl = saveSignedX0Results(mdl, block, signedFitSummary);
            mdl = saveSignedX0SourceData(mdl, block, data);
        elseif strcmp(modelType, 'weibullSignedBX0')
            mdl = saveSignedBX0Results(mdl, block, signedFitSummary);
            mdl = saveSignedBX0SourceData(mdl, block, data);
        end

        if strcmp(modelType, 'weibullfreeAll')
            sideData = buildSideFitData( ...
                xBaselinePre, yBaselinePre, ...
                xHorizontalOptoPre, yHorizontalOptoPre, ...
                tagHorizontalOptoPre, congrHorizontalOptoPre, ...
                xVerticalOptoPre, yVerticalOptoPre, ...
                tagVerticalOptoPre, congrVerticalOptoPre);

            rng(10000 + block, 'twister');
            [mdl.fittedParamsHorizontal(block,:), mdl.fitStatusHorizontal{block}] = ...
                fitSideWeibullSession(sideData.horizontal, ...
                objectiveFunction, fittedParams, lb, ub, opts);
            if ~strcmp(mdl.fitStatusHorizontal{block}, 'ok')
                warning('fitPsyMLE2:HorizontalFitFailed', ...
                    'Horizontal Weibull fit failed for local session row %d: %s', ...
                    block, mdl.fitStatusHorizontal{block});
            end

            rng(20000 + block, 'twister');
            [mdl.fittedParamsVertical(block,:), mdl.fitStatusVertical{block}] = ...
                fitSideWeibullSession(sideData.vertical, ...
                objectiveFunction, fittedParams, lb, ub, opts);
            if ~strcmp(mdl.fitStatusVertical{block}, 'ok')
                warning('fitPsyMLE2:VerticalFitFailed', ...
                    'Vertical Weibull fit failed for local session row %d: %s', ...
                    block, mdl.fitStatusVertical{block});
            end
        end
        
        fprintf('Fitted block #%.0f', block)
    end
    fprintf('Done!\n\n')
    fprintf('\n=========\n\n')
end

function sideData = buildSideFitData( ...
        xBaseline, yBaseline, ...
        xHorizontalOpto, yHorizontalOpto, tagHorizontalOpto, congrHorizontalOpto, ...
        xVerticalOpto, yVerticalOpto, tagVerticalOpto, congrVerticalOpto)

    xOpto = [xHorizontalOpto, xVerticalOpto];
    yOpto = [yHorizontalOpto, yVerticalOpto];
    tagOpto = [tagHorizontalOpto, tagVerticalOpto];
    congrOpto = [congrHorizontalOpto, congrVerticalOpto];

    sideData.horizontal = makeSideFitData( ...
        xBaseline(xBaseline <= 0), yBaseline(xBaseline <= 0), ...
        xOpto(tagOpto == 0 & congrOpto == 1), ...
        yOpto(tagOpto == 0 & congrOpto == 1), ...
        xOpto(tagOpto == 0 & congrOpto == -1), ...
        yOpto(tagOpto == 0 & congrOpto == -1));

    sideData.vertical = makeSideFitData( ...
        xBaseline(xBaseline >= 0), yBaseline(xBaseline >= 0), ...
        xOpto(tagOpto == 90 & congrOpto == 1), ...
        yOpto(tagOpto == 90 & congrOpto == 1), ...
        xOpto(tagOpto == 90 & congrOpto == -1), ...
        yOpto(tagOpto == 90 & congrOpto == -1));
end

function data = makeSideFitData( ...
        xBaseline, yBaseline, xConOpto, yConOpto, xInconOpto, yInconOpto)

    [data.xBaseline, data.yBaseline] = cleanSideXY(xBaseline, yBaseline);
    [data.xConOpto, data.yConOpto] = cleanSideXY(xConOpto, yConOpto);
    [data.xInconOpto, data.yInconOpto] = cleanSideXY(xInconOpto, yInconOpto);
end

function [x, y] = cleanSideXY(x, y)
    x = x(:)';
    y = y(:)';
    valid = isfinite(x) & isfinite(y);
    x = abs(x(valid));
    y = y(valid);
    [x, sortIdx] = sort(x);
    y = y(sortIdx);
end

function [fitRow, status] = fitSideWeibullSession( ...
        sideData, objectiveFunction, initialParams, lb, ub, opts)

    fitRow = nan(1, 13);
    status = 'ok';
    conditionLengths = [numel(sideData.xBaseline), ...
        numel(sideData.xConOpto), numel(sideData.xInconOpto)];
    if any(conditionLengths < 2)
        status = sprintf( ...
            'insufficient points [baseline=%d, con=%d, incon=%d]', ...
            conditionLengths(1), conditionLengths(2), conditionLengths(3));
        return;
    end

    [sumBaseline, successBaseline] = ...
        convertSideCounts(sideData.xBaseline, sideData.yBaseline, 'baseline');
    [sumConOpto, successConOpto] = ...
        convertSideCounts(sideData.xConOpto, sideData.yConOpto, 'opto');
    [sumInconOpto, successInconOpto] = ...
        convertSideCounts(sideData.xInconOpto, sideData.yInconOpto, 'opto');

    data = struct( ...
        'xBaseline', sideData.xBaseline, ...
        'yBaseline', sideData.yBaseline, ...
        'xConOpto', sideData.xConOpto, ...
        'yConOpto', sideData.yConOpto, ...
        'xInconOpto', sideData.xInconOpto, ...
        'yInconOpto', sideData.yInconOpto, ...
        'sumBaseline', sumBaseline, ...
        'successBaseline', successBaseline, ...
        'sumConOpto', sumConOpto, ...
        'successConOpto', successConOpto, ...
        'sumInconOpto', sumInconOpto, ...
        'successInconOpto', successInconOpto);

    try
        [fittedParams, nLL] = fitParameters( ...
            objectiveFunction, initialParams, lb, ub, opts, data);
        if numel(fittedParams) ~= 12 || any(~isfinite(fittedParams))
            status = 'optimizer returned incomplete or nonfinite parameters';
            return;
        end

        n = sum([sumBaseline, sumInconOpto]);
        k = sum(isfinite(fittedParams));
        [~, aicc, ~] = calculateAIC(nLL, k, n);
        fitRow = [fittedParams, aicc];
    catch fitError
        status = fitError.message;
    end
end

function [sumY, successY] = convertSideCounts(x, y, conditionType)
    switch conditionType
        case 'baseline'
            % Split baseline panels: 20 trials at zero, 10 elsewhere.
            sumY = 10 .* ones(size(x));
            sumY(x == 0) = 20;
        case 'opto'
            % Split con/incon panels: 10 trials at every contrast.
            sumY = 10 .* ones(size(x));
        otherwise
            error('Unknown side-fit condition type: %s', conditionType);
    end
    successY = (y ./ 100) .* sumY;
end

function data = buildSignedChoiceData(xBlocks, yBlocks, block)
    xBaseline = rmnan(squeeze(xBlocks(1,:,block)));
    yBaseline = rmnan(squeeze(yBlocks(1,:,block)));
    xHorizontal = rmnan(squeeze(xBlocks(2,:,block)));
    yHorizontal = rmnan(squeeze(yBlocks(2,:,block)));
    xVertical = rmnan(squeeze(xBlocks(3,:,block)));
    yVertical = rmnan(squeeze(yBlocks(3,:,block)));

    validateSignedXY(xBaseline, yBaseline, 'baseline', block);
    validateSignedXY(xHorizontal, yHorizontal, 'horizontal opto', block);
    validateSignedXY(xVertical, yVertical, 'vertical opto', block);

    [sumBaseline, successBaseline] = convertSignedChoiceCounts(...
        xBaseline, yBaseline, 'baseline');
    [sumHorizontal, successHorizontal] = convertSignedChoiceCounts(...
        xHorizontal, yHorizontal, 'opto');
    [sumVertical, successVertical] = convertSignedChoiceCounts(...
        xVertical, yVertical, 'opto');

    data = struct(...
        'xBaselineChoice', xBaseline, ...
        'yBaselineChoice', yBaseline, ...
        'xHorizontalOptoChoice', xHorizontal, ...
        'yHorizontalOptoChoice', yHorizontal, ...
        'xVerticalOptoChoice', xVertical, ...
        'yVerticalOptoChoice', yVertical, ...
        'sumBaselineChoice', sumBaseline, ...
        'successBaselineChoice', successBaseline, ...
        'sumHorizontalOptoChoice', sumHorizontal, ...
        'successHorizontalOptoChoice', successHorizontal, ...
        'sumVerticalOptoChoice', sumVertical, ...
        'successVerticalOptoChoice', successVertical);
end

function validateSignedXY(x, y, conditionName, block)
    if numel(x) ~= numel(y)
        error('fitPsyMLE2:SignedLengthMismatch', ...
            'Block %d %s signed x/y lengths differ (%d vs %d).', ...
            block, conditionName, numel(x), numel(y));
    end
    if any(y < 0 | y > 100)
        error('fitPsyMLE2:SignedProbabilityRange', ...
            'Block %d %s choice percentages must be within 0-100.', ...
            block, conditionName);
    end
    if ~any(x < 0) || ~any(x > 0)
        warning('fitPsyMLE2:SignedContrastCoverage', ...
            'Block %d %s does not contain both negative and positive signed contrasts.', ...
            block, conditionName);
    end
end

function [sumY, successY] = convertSignedChoiceCounts(x, y, conditionType)
    switch conditionType
        case 'baseline'
            % Signed baseline keeps the two zero points separate: 20 trials
            % at each signed zero entry, 10 trials at each nonzero entry.
            sumY = 10 .* ones(size(x));
            sumY(x == 0) = 20;
        case 'opto'
            % Horizontal- and vertical-opto signed rows have 10 trials per
            % signed contrast point, including zero.
            sumY = 10 .* ones(size(x));
        otherwise
            error('Unknown signed condition type: %s', conditionType);
    end
    successY = round((y ./ 100) .* sumY);
end

function out = mergeStructs(out, extra)
    names = fieldnames(extra);
    for idx = 1:numel(names)
        out.(names{idx}) = extra.(names{idx});
    end
end

function summary = fitSignedNoX0Comparison(objectiveFunction, initialParams, lb, ub, opts, data, fittedParams, nLL, aicc, nTrials)
    objectiveM0 = @(q9, fitData) objectiveFunction([q9(:)', 0], fitData);
    [noX0Params, noX0NLL] = fitParametersSimple(...
        objectiveM0, initialParams(1:9), lb(1:9), ub(1:9), opts, data);
    [~, noX0AICc] = calculateAIC(noX0NLL, 9, nTrials);

    deltaAICc = noX0AICc - aicc;
    weights = akaikeWeights([noX0AICc, aicc]);

    summary = struct(...
        'fitParams', fittedParams, ...
        'nLL', nLL, ...
        'AICc', aicc, ...
        'noX0FitParams', noX0Params, ...
        'noX0NLL', noX0NLL, ...
        'noX0AICc', noX0AICc, ...
        'deltaAICcX0', deltaAICc, ...
        'akaikeWeightNoX0', weights(1), ...
        'akaikeWeightX0', weights(2), ...
        'X0Horizontal', fittedParams(10), ...
        'X0Vertical', -fittedParams(10), ...
        'fitStatus', 'ok');
end

function summary = fitSignedBX0NoX0Comparison(objectiveFunction, initialParams, lb, ub, opts, data, fittedParams, nLL, aicc, nTrials)
    objectiveM0 = @(q10, fitData) objectiveFunction([q10(:)', 0], fitData);
    [noX0Params, noX0NLL] = fitParametersSimple(...
        objectiveM0, initialParams(1:10), lb(1:10), ub(1:10), opts, data);
    [~, noX0AICc] = calculateAIC(noX0NLL, 10, nTrials);

    deltaAICc = noX0AICc - aicc;
    weights = akaikeWeights([noX0AICc, aicc]);

    summary = struct(...
        'fitParams', fittedParams, ...
        'nLL', nLL, ...
        'AICc', aicc, ...
        'noX0FitParams', noX0Params, ...
        'noX0NLL', noX0NLL, ...
        'noX0AICc', noX0AICc, ...
        'deltaAICcX0', deltaAICc, ...
        'akaikeWeightBOnly', weights(1), ...
        'akaikeWeightBX0', weights(2), ...
        'BHorizontal', 50 - fittedParams(10), ...
        'BVertical', 50 + fittedParams(10), ...
        'X0Horizontal', fittedParams(11), ...
        'X0Vertical', -fittedParams(11), ...
        'deltaB', fittedParams(10), ...
        'deltaX0', fittedParams(11), ...
        'fitStatus', 'ok');
end

function weights = akaikeWeights(aiccValues)
    finiteAICc = aiccValues(isfinite(aiccValues));
    if isempty(finiteAICc)
        weights = nan(size(aiccValues));
        return;
    end
    delta = aiccValues - min(finiteAICc);
    relLike = exp(-0.5 .* delta);
    weights = relLike ./ sum(relLike(isfinite(relLike)));
end

function mdl = saveSignedX0Results(mdl, block, summary)
    mdl.signedX0.fitParams(block,:) = summary.fitParams;
    mdl.signedX0.nLL(block) = summary.nLL;
    mdl.signedX0.AICc(block) = summary.AICc;
    mdl.signedX0.noX0FitParams(block,:) = summary.noX0FitParams;
    mdl.signedX0.noX0NLL(block) = summary.noX0NLL;
    mdl.signedX0.noX0AICc(block) = summary.noX0AICc;
    mdl.signedX0.deltaAICcX0(block) = summary.deltaAICcX0;
    mdl.signedX0.akaikeWeightX0(block) = summary.akaikeWeightX0;
    mdl.signedX0.akaikeWeightNoX0(block) = summary.akaikeWeightNoX0;
    mdl.signedX0.X0Horizontal(block) = summary.X0Horizontal;
    mdl.signedX0.X0Vertical(block) = summary.X0Vertical;
    mdl.signedX0.fitStatus{block} = summary.fitStatus;
end

function mdl = saveSignedBX0Results(mdl, block, summary)
    mdl.signedBX0.fitParams(block,:) = summary.fitParams;
    mdl.signedBX0.nLL(block) = summary.nLL;
    mdl.signedBX0.AICc(block) = summary.AICc;
    mdl.signedBX0.noX0FitParams(block,:) = summary.noX0FitParams;
    mdl.signedBX0.noX0NLL(block) = summary.noX0NLL;
    mdl.signedBX0.noX0AICc(block) = summary.noX0AICc;
    mdl.signedBX0.deltaAICcX0(block) = summary.deltaAICcX0;
    mdl.signedBX0.akaikeWeightBX0(block) = summary.akaikeWeightBX0;
    mdl.signedBX0.akaikeWeightBOnly(block) = summary.akaikeWeightBOnly;
    mdl.signedBX0.BHorizontal(block) = summary.BHorizontal;
    mdl.signedBX0.BVertical(block) = summary.BVertical;
    mdl.signedBX0.X0Horizontal(block) = summary.X0Horizontal;
    mdl.signedBX0.X0Vertical(block) = summary.X0Vertical;
    mdl.signedBX0.deltaB(block) = summary.deltaB;
    mdl.signedBX0.deltaX0(block) = summary.deltaX0;
    mdl.signedBX0.fitStatus{block} = summary.fitStatus;

    params = summary.fitParams;
    mdl.signedBX0.ABaseline(block) = params(1);
    mdl.signedBX0.alphaBaseline(block) = params(2);
    mdl.signedBX0.betaBaseline(block) = params(3);
    mdl.signedBX0.ACon(block) = params(4);
    mdl.signedBX0.alphaCon(block) = params(5);
    mdl.signedBX0.betaCon(block) = params(6);
    mdl.signedBX0.AIncon(block) = params(7);
    mdl.signedBX0.alphaIncon(block) = params(8);
    mdl.signedBX0.betaIncon(block) = params(9);
end

function mdl = saveSignedBX0SourceData(mdl, block, data)
    mdl.signedBX0.xBaselineChoice(block,:) = ...
        padArray(data.xBaselineChoice, 24, 2, NaN);
    mdl.signedBX0.yBaselineChoice(block,:) = ...
        padArray(data.yBaselineChoice, 24, 2, NaN);
    mdl.signedBX0.nBaselineChoice(block,:) = ...
        padArray(data.sumBaselineChoice, 24, 2, NaN);
    mdl.signedBX0.successBaselineChoice(block,:) = ...
        padArray(data.successBaselineChoice, 24, 2, NaN);

    mdl.signedBX0.xHorizontalOptoChoice(block,:) = ...
        padArray(data.xHorizontalOptoChoice, 24, 2, NaN);
    mdl.signedBX0.yHorizontalOptoChoice(block,:) = ...
        padArray(data.yHorizontalOptoChoice, 24, 2, NaN);
    mdl.signedBX0.nHorizontalOptoChoice(block,:) = ...
        padArray(data.sumHorizontalOptoChoice, 24, 2, NaN);
    mdl.signedBX0.successHorizontalOptoChoice(block,:) = ...
        padArray(data.successHorizontalOptoChoice, 24, 2, NaN);

    mdl.signedBX0.xVerticalOptoChoice(block,:) = ...
        padArray(data.xVerticalOptoChoice, 24, 2, NaN);
    mdl.signedBX0.yVerticalOptoChoice(block,:) = ...
        padArray(data.yVerticalOptoChoice, 24, 2, NaN);
    mdl.signedBX0.nVerticalOptoChoice(block,:) = ...
        padArray(data.sumVerticalOptoChoice, 24, 2, NaN);
    mdl.signedBX0.successVerticalOptoChoice(block,:) = ...
        padArray(data.successVerticalOptoChoice, 24, 2, NaN);
end

function mdl = saveSignedX0SourceData(mdl, block, data)
    mdl.signedX0.xBaselineChoice(block,:) = ...
        padArray(data.xBaselineChoice, 24, 2, NaN);
    mdl.signedX0.yBaselineChoice(block,:) = ...
        padArray(data.yBaselineChoice, 24, 2, NaN);
    mdl.signedX0.nBaselineChoice(block,:) = ...
        padArray(data.sumBaselineChoice, 24, 2, NaN);
    mdl.signedX0.successBaselineChoice(block,:) = ...
        padArray(data.successBaselineChoice, 24, 2, NaN);

    mdl.signedX0.xHorizontalOptoChoice(block,:) = ...
        padArray(data.xHorizontalOptoChoice, 24, 2, NaN);
    mdl.signedX0.yHorizontalOptoChoice(block,:) = ...
        padArray(data.yHorizontalOptoChoice, 24, 2, NaN);
    mdl.signedX0.nHorizontalOptoChoice(block,:) = ...
        padArray(data.sumHorizontalOptoChoice, 24, 2, NaN);
    mdl.signedX0.successHorizontalOptoChoice(block,:) = ...
        padArray(data.successHorizontalOptoChoice, 24, 2, NaN);

    mdl.signedX0.xVerticalOptoChoice(block,:) = ...
        padArray(data.xVerticalOptoChoice, 24, 2, NaN);
    mdl.signedX0.yVerticalOptoChoice(block,:) = ...
        padArray(data.yVerticalOptoChoice, 24, 2, NaN);
    mdl.signedX0.nVerticalOptoChoice(block,:) = ...
        padArray(data.sumVerticalOptoChoice, 24, 2, NaN);
    mdl.signedX0.successVerticalOptoChoice(block,:) = ...
        padArray(data.successVerticalOptoChoice, 24, 2, NaN);
end
function [fittedParams, nLL] = fitParametersSimple(objectiveFunction, params0, lb, ub, opts, data)
    toParams = @(u) lb + u(:)' .* (ub - lb);
    toUnit = @(p) (p(:)' - lb) ./ (ub - lb);
    u0 = min(max(toUnit(params0), 0), 1);
    obj = @(u) objectiveFunction(toParams(min(max(u, 0), 1)), data);

    haveFMC = exist('fmincon','file') == 2;
    if haveFMC
        fopts = optimoptions('fmincon', ...
            'Algorithm','interior-point', ...
            'Display','off', ...
            'MaxFunctionEvaluations', max(5000, 20 .* opts.MaxFunEvals), ...
            'FiniteDifferenceType','central');
        uFit = fmincon(obj, u0, [], [], [], [], zeros(size(u0)), ...
            ones(size(u0)), [], fopts);
    else
        [uFit, ~] = fminsearchbnd(obj, u0, zeros(size(u0)), ...
            ones(size(u0)), opts);
    end

    fittedParams = toParams(min(max(uFit, 0), 1));
    nLL = objectiveFunction(fittedParams, data);
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
            'Display', 'off', ...
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
            'Display','off', ...
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

function mdl = saveModelResults(mdl, block, ...
    xBaseline, yBaseline, ...
    xInconOpto, yInconOpto, ...
    xConOpto, yConOpto, ...
    xBaselinePre, yBaselinePre, tagBaselinePre, ...
    xHorizontalOptoPre, yHorizontalOptoPre, tagHorizontalOptoPre, congrHorizontalOptoPre, ...
    xVerticalOptoPre, yVerticalOptoPre, tagVerticalOptoPre, congrVerticalOptoPre, ...
    mdlStruct, objectiveFunction, ub, lb, params, options, ...
    fittedParams, aicc, manualFittingFlag)

    % Merged data used for fitting
    mdl.xBaseline(block,:) = padArray(xBaseline, 12, 2, NaN);
    mdl.yBaseline(block,:) = padArray(yBaseline, 12, 2, NaN);
    mdl.xInconOpto(block,:) = padArray(xInconOpto, 12, 2, NaN);
    mdl.yInconOpto(block,:) = padArray(yInconOpto, 12, 2, NaN);
    mdl.xConOpto(block,:) = padArray(xConOpto, 12, 2, NaN);
    mdl.yConOpto(block,:) = padArray(yConOpto, 12, 2, NaN);

    % Pre-merged/split data for plotting
    % x = signed Gabor contrast
    % y = percent correct
    % visualTag: 0 = horizontal/V0 visual tag, 90 = vertical/V90 visual tag
    % congruency: 1 = congruent, -1 = incongruent, NaN = not applicable/baseline

    mdl.xBaselinePreMerge(block,:) = padArray(xBaselinePre, 12, 2, NaN);
    mdl.yBaselinePreMerge(block,:) = padArray(yBaselinePre, 12, 2, NaN);
    mdl.visualTagBaselinePreMerge(block,:) = padArray(tagBaselinePre, 12, 2, NaN);

    mdl.xHorizontalOptoPreMerge(block,:) = padArray(xHorizontalOptoPre, 12, 2, NaN);
    mdl.yHorizontalOptoPreMerge(block,:) = padArray(yHorizontalOptoPre, 12, 2, NaN);
    mdl.visualTagHorizontalOptoPreMerge(block,:) = padArray(tagHorizontalOptoPre, 12, 2, NaN);
    mdl.congruencyHorizontalOptoPreMerge(block,:) = padArray(congrHorizontalOptoPre, 12, 2, NaN);

    mdl.xVerticalOptoPreMerge(block,:) = padArray(xVerticalOptoPre, 12, 2, NaN);
    mdl.yVerticalOptoPreMerge(block,:) = padArray(yVerticalOptoPre, 12, 2, NaN);
    mdl.visualTagVerticalOptoPreMerge(block,:) = padArray(tagVerticalOptoPre, 12, 2, NaN);
    mdl.congruencyVerticalOptoPreMerge(block,:) = padArray(congrVerticalOptoPre, 12, 2, NaN);

    % Fitting data
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

function [xBaselineAll, yBaselineAll, xOptoAll, yOptoAll, ...
          xBaselinePreAll, yBaselinePreAll, tagBaselinePreAll, ...
          xHorizontalOptoPreAll, yHorizontalOptoPreAll, tagHorizontalOptoPreAll, congrHorizontalOptoPreAll, ...
          xVerticalOptoPreAll, yVerticalOptoPreAll, tagVerticalOptoPreAll, congrVerticalOptoPreAll] = processConditionsBlocks(xBlocks, yBlocks, modelType)

    % Arrays
    % xBlocks/yBlocks rows:
    %   1 = baselineRaw
    %   2 = horizontalOptoRaw
    %   3 = verticalOptoRaw
    %
    % Each row was generated upstream as [V0, V90].
    % Therefore:
    %   first half  = horizontal/V0 visual tag
    %   second half = vertical/V90 visual tag
    %
    % For pre-merged plotting:
    %   x stays signed contrast
    %   y becomes percent correct
    %   visualTag: 0 = horizontal/V0, 90 = vertical/V90
    %   congruency: 1 = congruent, -1 = incongruent

    nBlocks = size(xBlocks, 3);
    maxBaseLen = 0;
    maxOptoLen = 0;
    maxSingleLen = 0;

    % First pass: determine maximum lengths needed for preallocation
    for block = 1:nBlocks
        xBaseline   = rmnan(squeeze(xBlocks(1,:,block)));
        xHorizontal = rmnan(squeeze(xBlocks(2,:,block)));
        xVertical   = rmnan(squeeze(xBlocks(3,:,block)));

        maxBaseLen = max(maxBaseLen, numel(xBaseline));
        maxOptoLen = max(maxOptoLen, numel(xHorizontal) + numel(xVertical));
        maxSingleLen = max([maxSingleLen, numel(xBaseline), numel(xHorizontal), numel(xVertical)]);
    end

    % Merged/fitting arrays
    xBaselineAll = NaN(nBlocks, maxBaseLen);
    yBaselineAll = NaN(nBlocks, maxBaseLen);
    xOptoAll = NaN(nBlocks, maxOptoLen);
    yOptoAll = NaN(nBlocks, maxOptoLen);

    % Pre-merged/split plotting arrays
    xBaselinePreAll = NaN(nBlocks, maxSingleLen);
    yBaselinePreAll = NaN(nBlocks, maxSingleLen);
    tagBaselinePreAll = NaN(nBlocks, maxSingleLen);

    xHorizontalOptoPreAll = NaN(nBlocks, maxSingleLen);
    yHorizontalOptoPreAll = NaN(nBlocks, maxSingleLen);
    tagHorizontalOptoPreAll = NaN(nBlocks, maxSingleLen);
    congrHorizontalOptoPreAll = NaN(nBlocks, maxSingleLen);

    xVerticalOptoPreAll = NaN(nBlocks, maxSingleLen);
    yVerticalOptoPreAll = NaN(nBlocks, maxSingleLen);
    tagVerticalOptoPreAll = NaN(nBlocks, maxSingleLen);
    congrVerticalOptoPreAll = NaN(nBlocks, maxSingleLen);

    switch modelType
        case {'bill','weibullfreeAll','weibullSignedX0','weibullSignedBX0'}

            for block = 1:nBlocks

                %% =========================
                %  Baseline
                %  =========================

                xBaselineRaw = rmnan(squeeze(xBlocks(1,:,block)));
                yBaselineRaw = rmnan(squeeze(yBlocks(1,:,block)));

                nBase = numel(xBaselineRaw);
                if mod(nBase, 2) ~= 0
                    error('Baseline raw vector has odd length in block %.0f. Expected [V0, V90] halves.', block);
                end

                visualTagBaselineRaw = makeVisualTagVector(nBase);

                % Sort x, y, and visual tag together
                [xBaseline, sortIdx] = sort(xBaselineRaw);
                yBaseline = yBaselineRaw(sortIdx);
                visualTagBaseline = visualTagBaselineRaw(sortIdx);

                % Pre-merged baseline for plotting:
                % signed contrast, percent correct.
                % Baseline duplicate zeros should be averaged.
                [xBaselinePre, yBaselinePre, tagBaselinePre] = ...
                    makePreMergePercentCorrect(xBaseline, yBaseline, visualTagBaseline, true);

                xBaselinePreAll(block, 1:numel(xBaselinePre)) = xBaselinePre;
                yBaselinePreAll(block, 1:numel(yBaselinePre)) = yBaselinePre;
                tagBaselinePreAll(block, 1:numel(tagBaselinePre)) = tagBaselinePre;

                % Original merged baseline used for fitting
                numVal = numel(xBaseline);
                xBaselineMerged = mean([fliplr(-xBaseline(1:numVal/2)); xBaseline(numVal/2+1:end)]);
                yBaselineMerged = mean([fliplr(100 - yBaseline(1:numVal/2)); yBaseline(numVal/2+1:end)]);

                xBaselineAll(block, 1:numel(xBaselineMerged)) = xBaselineMerged;
                yBaselineAll(block, 1:numel(yBaselineMerged)) = yBaselineMerged;


                %% =========================
                %  Horizontal opto
                %  =========================

                xHorizontalRaw = rmnan(squeeze(xBlocks(2,:,block)));
                yHorizontalRaw = rmnan(squeeze(yBlocks(2,:,block)));

                nHorizontal = numel(xHorizontalRaw);
                if mod(nHorizontal, 2) ~= 0
                    error('Horizontal opto raw vector has odd length in block %.0f. Expected [V0, V90] halves.', block);
                end

                visualTagHorizontalRaw = makeVisualTagVector(nHorizontal);

                % Sort x, y, and visual tag together
                [xHorizontal, sortIdxH] = sort(xHorizontalRaw);
                yHorizontal = yHorizontalRaw(sortIdxH);
                visualTagHorizontal = visualTagHorizontalRaw(sortIdxH);

                % Pre-merged horizontal opto:
                % DO NOT merge duplicate x=0 points.
                [xHorizontalOptoPre, yHorizontalOptoPre, tagHorizontalOptoPre] = ...
                    makePreMergePercentCorrect(xHorizontal, yHorizontal, visualTagHorizontal, false);

                % Horizontal opto congruency:
                % visual V0/horizontal tag = congruent
                % visual V90/vertical tag = incongruent
                congrHorizontalOptoPre = NaN(size(tagHorizontalOptoPre));
                congrHorizontalOptoPre(tagHorizontalOptoPre == 0) = 1;
                congrHorizontalOptoPre(tagHorizontalOptoPre == 90) = -1;

                xHorizontalOptoPreAll(block, 1:numel(xHorizontalOptoPre)) = xHorizontalOptoPre;
                yHorizontalOptoPreAll(block, 1:numel(yHorizontalOptoPre)) = yHorizontalOptoPre;
                tagHorizontalOptoPreAll(block, 1:numel(tagHorizontalOptoPre)) = tagHorizontalOptoPre;
                congrHorizontalOptoPreAll(block, 1:numel(congrHorizontalOptoPre)) = congrHorizontalOptoPre;


                %% =========================
                %  Vertical opto
                %  =========================

                xVerticalRaw = rmnan(squeeze(xBlocks(3,:,block)));
                yVerticalRaw = rmnan(squeeze(yBlocks(3,:,block)));

                nVertical = numel(xVerticalRaw);
                if mod(nVertical, 2) ~= 0
                    error('Vertical opto raw vector has odd length in block %.0f. Expected [V0, V90] halves.', block);
                end

                visualTagVerticalRaw = makeVisualTagVector(nVertical);

                % Sort x, y, and visual tag together
                [xVertical, sortIdxV] = sort(xVerticalRaw);
                yVertical = yVerticalRaw(sortIdxV);
                visualTagVertical = visualTagVerticalRaw(sortIdxV);

                % Pre-merged vertical opto:
                % DO NOT merge duplicate x=0 points.
                [xVerticalOptoPre, yVerticalOptoPre, tagVerticalOptoPre] = ...
                    makePreMergePercentCorrect(xVertical, yVertical, visualTagVertical, false);

                % Vertical opto congruency:
                % visual V0/horizontal tag = incongruent
                % visual V90/vertical tag = congruent
                congrVerticalOptoPre = NaN(size(tagVerticalOptoPre));
                congrVerticalOptoPre(tagVerticalOptoPre == 0) = -1;
                congrVerticalOptoPre(tagVerticalOptoPre == 90) = 1;

                xVerticalOptoPreAll(block, 1:numel(xVerticalOptoPre)) = xVerticalOptoPre;
                yVerticalOptoPreAll(block, 1:numel(yVerticalOptoPre)) = yVerticalOptoPre;
                tagVerticalOptoPreAll(block, 1:numel(tagVerticalOptoPre)) = tagVerticalOptoPre;
                congrVerticalOptoPreAll(block, 1:numel(congrVerticalOptoPre)) = congrVerticalOptoPre;


                %% =========================
                %  Existing merged con/incon fitting data
                %  =========================

                contrastNeg = 1:numel(xHorizontal)/2;
                contrastPos = numel(xHorizontal)/2+1:numel(xHorizontal);

                % Convert negative/V0 side to percent correct
                yHorizontalCorrect = yHorizontal;
                yHorizontalCorrect(contrastNeg) = 100 - yHorizontalCorrect(contrastNeg);

                yVerticalCorrect = yVertical;
                yVerticalCorrect(contrastNeg) = 100 - yVerticalCorrect(contrastNeg);

                % Congruent:
                % horizontal opto + V0, vertical opto + V90
                xConOpto = mean([-fliplr(xHorizontal(contrastNeg)); xVertical(contrastPos)]);
                yConOpto = mean([fliplr(yHorizontalCorrect(contrastNeg)); yVerticalCorrect(contrastPos)]);

                % Incongruent:
                % horizontal opto + V90, vertical opto + V0
                xInconOpto = mean([xHorizontal(contrastPos); -fliplr(xVertical(contrastNeg))]);
                yInconOpto = mean([yHorizontalCorrect(contrastPos); fliplr(yVerticalCorrect(contrastNeg))]);

                optoLen = numel([-xInconOpto, xConOpto]);
                xOptoAll(block, 1:optoLen) = [-xInconOpto, xConOpto];
                yOptoAll(block, 1:optoLen) = [-yInconOpto, yConOpto];
            end

        otherwise
            error('This updated pre-merge/tag logic is currently implemented for bill, weibullfreeAll, weibullSignedX0, and weibullSignedBX0 only.');
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

    % Signed decision-boundary Weibull with fitted opponent X0
    config.models.weibullSignedX0 = struct(...
        'getInitParams', @getWeibullSignedX0InitParams, ...
        'getModelFuncs', @getWeibullSignedX0ModelFuncs, ...
        'headers', {'A^{bl}','\alpha^{bl}','\beta^{bl}', ...
                   'A^{HOpto}','\alpha^{HOpto}','\beta^{HOpto}', ...
                   'A^{VOpto}','\alpha^{VOpto}','\beta^{VOpto}', ...
                   '\DeltaX0','AICc^{X0}'} ...
    );

    % Production signed Weibull with opponent B and signed decision boundary.
    config.models.weibullSignedBX0 = struct(...
        'getInitParams', @getWeibullSignedBX0InitParams, ...
        'getModelFuncs', @getWeibullSignedBX0ModelFuncs, ...
        'headers', {'A^{bl}','\alpha^{bl}','\beta^{bl}', ...
                   'A^{con}','\alpha^{con}','\beta^{con}', ...
                   'A^{incon}','\alpha^{incon}','\beta^{incon}', ...
                   '\DeltaB','\DeltaX0','AICc^{BX0}'} ...
    );
    
end
function [xOut, yOut] = mergeDuplicateXForPlot(xIn, yIn)
    % Merge duplicate x-values by averaging y-values.
    % Intended for pre-merged plotting data, especially duplicate 0 contrast.
    %
    % Example:
    %   x = [-30 -15 0 0 15 30]
    %   y = [ 10  30 70 80 90 95]
    % becomes:
    %   x = [-30 -15 0 15 30]
    %   y = [ 10  30 75 90 95]

    xIn = xIn(:)';
    yIn = yIn(:)';

    validIdx = ~isnan(xIn) & ~isnan(yIn);
    xIn = xIn(validIdx);
    yIn = yIn(validIdx);

    [xOut, ~, groupIdx] = unique(xIn, 'stable');

    yOut = nan(size(xOut));
    for ii = 1:numel(xOut)
        yOut(ii) = mean(yIn(groupIdx == ii), 'omitnan');
    end
end

function visualTag = makeVisualTagVector(nVal)
    % Upstream analyzeBlockPsychometrics stores raw vectors as [V0, V90].
    % Therefore:
    %   first half  = V0 / horizontal visual tag = 0
    %   second half = V90 / vertical visual tag = 90

    if mod(nVal, 2) ~= 0
        error('Expected even number of raw contrast values because data should be [V0, V90].');
    end

    visualTag = NaN(1, nVal);
    visualTag(1:nVal/2) = 0;
    visualTag(nVal/2+1:end) = 90;
end

function [xOut, yCorrectOut, visualTagOut] = makePreMergePercentCorrect(xIn, yIn, visualTagIn, mergeDuplicateBaselineZeros)
    % Convert raw pre-merged data to:
    %   x = signed contrast
    %   y = percent correct
    %
    % visualTag:
    %   0  = horizontal/V0 visual tag
    %   90 = vertical/V90 visual tag
    %
    % For visualTag == 0, raw y is percent vertical report,
    % so percent correct = 100 - y.
    %
    % For visualTag == 90, percent correct = y.

    xIn = xIn(:)';
    yIn = yIn(:)';
    visualTagIn = visualTagIn(:)';

    validIdx = ~isnan(xIn) & ~isnan(yIn) & ~isnan(visualTagIn);
    xIn = xIn(validIdx);
    yIn = yIn(validIdx);
    visualTagIn = visualTagIn(validIdx);

    yCorrect = yIn;
    yCorrect(visualTagIn == 0) = 100 - yCorrect(visualTagIn == 0);

    if mergeDuplicateBaselineZeros
        [xOut, yCorrectOut, visualTagOut] = mergeDuplicateXAndTagsForBaseline(xIn, yCorrect, visualTagIn);
    else
        xOut = xIn;
        yCorrectOut = yCorrect;
        visualTagOut = visualTagIn;
    end
end

function [xOut, yOut, tagOut] = mergeDuplicateXAndTagsForBaseline(xIn, yIn, tagIn)
    % For baseline only:
    % merge duplicate x-values, especially the two x=0 points.
    % y is averaged.
    %
    % For tag:
    %   if merged tags differ, use 45 to mean "merged V0/V90 baseline".
    %   if same, keep the tag.

    xIn = xIn(:)';
    yIn = yIn(:)';
    tagIn = tagIn(:)';

    validIdx = ~isnan(xIn) & ~isnan(yIn) & ~isnan(tagIn);
    xIn = xIn(validIdx);
    yIn = yIn(validIdx);
    tagIn = tagIn(validIdx);

    [xOut, ~, groupIdx] = unique(xIn, 'stable');

    yOut = nan(size(xOut));
    tagOut = nan(size(xOut));

    for ii = 1:numel(xOut)
        thisIdx = groupIdx == ii;
        yOut(ii) = mean(yIn(thisIdx), 'omitnan');

        uniqueTags = unique(tagIn(thisIdx));
        uniqueTags = uniqueTags(~isnan(uniqueTags));

        if numel(uniqueTags) == 1
            tagOut(ii) = uniqueTags;
        else
            tagOut(ii) = 45; % merged baseline V0/V90, mainly for x=0
        end
    end
end
