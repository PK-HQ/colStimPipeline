function [mdl, mdlAvg] = fitBayesianModelMLE(xBlocks, yBlocks, modelType)
% Modified version of fitNakaRushtonMLE3 to use the normalization model
% Designed to be a drop-in replacement with the same interface

    % Define the number of blocks and conditions
    [nConditions, nContrasts, nBlocks] = size(xBlocks);
    mdl = []; mdlAvg = [];
    
    % Process data across all blocks, get average
    [xBaselineAll, yBaselineAll, xOptoAll, yOptoAll] = processConditionsBlocks(xBlocks, yBlocks);
    
    % Get average data for model fitting
    xBaselineAverage = rmnan(reshape(xBaselineAll, 1, numel(xBaselineAll)));
    yBaselineAverage = rmnan(reshape(yBaselineAll, 1, numel(yBaselineAll)));
    xOptoAverage = rmnan(reshape(xOptoAll, 1, numel(xOptoAll)));
    yOptoAverage = rmnan(reshape(yOptoAll, 1, numel(yOptoAll)));
    
    fprintf('Fitting normalization model, independent parameters per session...\n');
    
    % Set up headers and parameter storage based on model type
    switch modelType
        case 'bill'
            mdl.headers = {'a', 'b', 'l', 'w', 'g0', 'n', 'rmx', 'e', 'o', 'AICc'};
            initialParams = [0.5, 0.5, 0.1, 0.1, 75, 6, 50, 0.5, 50];
            lb = [0.01, 0.01, 0.01, 0.01, 40, 4, 10, 0.3, 10];
            ub = [1.0, 1.0, 1.0, 1.0, 100, 10, 100, 0.7, 100];
        otherwise
            error('Unsupported model type: %s. Only "bill" model is supported.', modelType);
    end
    
    % Setup model options
    modelOptions = struct('nTrials', 2000, 'useGPU', true, 'smartAllocation', true, 'verbose', false);
    
    % Fit for each block
    for block = 1%:nBlocks
        xBaseline = rmnan(xBaselineAll(block, :));
        yBaseline = rmnan(yBaselineAll(block, :));
        
        xOpto = rmnan(xOptoAll(block, :));
        yOpto = rmnan(yOptoAll(block, :));
        
        % Skip blocks with insufficient data
        if length(xBaseline) < 3 || length(xOpto) < 3
            fprintf('Skipping block %d due to insufficient data\n', block);
            continue;
        end
        
        % Fit the model using our new function
        [bestParams, fitInfo] = fitBayesianModel(xBaseline, yBaseline, xOpto, yOpto, initialParams, ...
            'LowerBounds', lb, 'UpperBounds', ub, 'ModelOptions', modelOptions, ...
            'OptimMethod', 'GlobalThenLocal', 'Verbose', true);
        
        % Store results - Fix: Use proper padding approach
        % Create padded arrays with NaN
        padded_xBaseline = NaN(1, 12);
        padded_yBaseline = NaN(1, 12);
        padded_xOpto = NaN(1, 12);
        padded_yOpto = NaN(1, 12);
        
        % Copy actual data (up to 12 elements)
        len_xBaseline = min(length(xBaseline), 12);
        len_yBaseline = min(length(yBaseline), 12);
        len_xOpto = min(length(xOpto), 12);
        len_yOpto = min(length(yOpto), 12);
        
        padded_xBaseline(1:len_xBaseline) = xBaseline(1:len_xBaseline);
        padded_yBaseline(1:len_yBaseline) = yBaseline(1:len_yBaseline);
        padded_xOpto(1:len_xOpto) = xOpto(1:len_xOpto);
        padded_yOpto(1:len_yOpto) = yOpto(1:len_yOpto);

        % Store mdl
        minBound = 1e-10;
        maxBound = 100 - minBound;
        mdlBaseline = @(x, params) max(minBound, min(maxBound, normMdlSimOptimized(x, params, 'baseline'))); % zero opto input
        mdlOpto = @(x, params) max(minBound, min(maxBound, normMdlSimOptimized(x, params, 'opto'))); % nonzero opto input
        
        mdl.mdlBaseline = mdlBaseline;
        mdl.mdlOpto = mdlOpto;

        % Store in mdl structure
        mdl.xBaseline(block, :) = padded_xBaseline;
        mdl.yBaseline(block, :) = padded_yBaseline;
        mdl.xOpto(block, :) = padded_xOpto;
        mdl.yOpto(block, :) = padded_yOpto;
        
        % Save parameter values with AICc
        mdl.fittedParams(block, :, 1) = [bestParams, fitInfo.aicc];
        mdl.r2Baseline(block) = fitInfo.rSquared_baseline;
        mdl.r2Opto(block) = fitInfo.rSquared_opto;
        
        fprintf('Fitted block #%d, AICc = %.2f, R² (baseline) = %.4f, R² (opto) = %.4f\n', ...
            block, fitInfo.aicc, fitInfo.rSquared_baseline, fitInfo.rSquared_opto);
    end
    
    % Fit average data across all blocks
    %{
    if ~isempty(xBaselineAverage) && ~isempty(xOptoAverage)
        [avgParams, avgFitInfo] = fitBayesianModel(xBaselineAverage, yBaselineAverage, yOptoAverage, initialParams, ...
            'LowerBounds', lb, 'UpperBounds', ub, 'ModelOptions', modelOptions, ...
            'OptimMethod', 'GlobalThenLocal', 'Verbose', true);
        
        mdlAvg.headers = mdl.headers;
        mdlAvg.fittedParams(1, :, 1) = [avgParams, avgFitInfo.aicc];
        mdlAvg.r2Baseline = avgFitInfo.rSquared_baseline;
        mdlAvg.r2Opto = avgFitInfo.rSquared_opto;
    end
    %}
    % Display results
    fprintf('\nFitting Results:\n');
    displayTable(mdl.fittedParams(:, :, 1), mdl.headers);
    fprintf('\n=========\n\n');
end


function [xBaselineAll, yBaselineAll, xOptoAll, yOptoAll] = processConditionsBlocks(xBlocks, yBlocks)
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

    % Second pass: process each block and store results
    for block = 1:nBlocks
        % Process baseline
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

function [sumY, successY] = convert2counts(x,y)
x=rmnan(x); y=rmnan(y);
% Total trials
sumY =  [(20*ones(1,sum(x<0))), (20*ones(1,sum(x==0))),...
    (20*ones(1,sum(x>0)))]; % trials per contrast level
% Convert percent correct into success counts
successY = (y / 100) .* sumY;
end

function [aic, aicc, bic] = calculateAIC(nLL, k, n)
    % Calculate AIC
    aic = 2 * k + 2 * nLL;  % since nll is negative log-likelihood
    
    % Calculate Corrected AIC (AICc)
    if n > k + 1
        aicc = aic + (2 * k * (k + 1)) / (n - k - 1);
    else
        aicc = Inf; % AICc is undefined when n <= k + 1
    end

    % Calculate BIC
    bic = log(n) * k + 2 * nLL;

    % Return AIC, Corrected AIC (AICc), and BIC
    [aic, aicc, bic];
end

