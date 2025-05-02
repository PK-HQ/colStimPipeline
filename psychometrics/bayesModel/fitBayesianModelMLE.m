function [mdl, mdlAvg] = fitBayesianModelMLE(xBlocks, yBlocks, modelType)
% FITBAYESIANMODELMLE - Enhanced model fitting with multi-stage optimization
%
% A complete rewrite of the fitting function to implement:
% 1. Parameter sensitivity analysis to identify influential parameters
% 2. Parameter scaling for better optimization
% 3. Multi-stage optimization (global, local, polish)
% 4. Improved performance monitoring and reporting
%
% Inputs:
%   xBlocks - Contrast levels for all conditions and blocks
%   yBlocks - Performance data for all conditions and blocks
%   modelType - Model type to fit (only 'bill' currently supported)
%
% Outputs:
%   mdl - Structure with fitted model information
%   mdlAvg - Structure with fitted model to average data across blocks

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
    
    fprintf('Fitting normalization model with multi-stage optimization...\n');
    
    % Set up headers and parameter storage based on model type
    switch modelType
        case 'bill'
            mdl.headers = {'a', 'b', 'l', 'w', 'g0', 'n', 'rmx', 'e', 'o', 'AICc'};
            
            % Parameters with appropriate ranges based on model understanding
            initialParams = [0.2, 0.2, 0.2, 0.2, 75, 4, 100, 0.5, 10];
            
            % Lower and upper bounds - carefully chosen based on model physics
            lb = [0.01, 0.01, 0.01, 0.01, 50, 2, 50, 0.3, 0];
            ub = [0.99, 0.99, 0.99, 0.99, 120, 10, 150, 0.7, 50];
        otherwise
            error('Unsupported model type: %s. Only "bill" model is supported.', modelType);
    end
    
    % Setup model options - adjust trials for faster sensitivity analysis
    initialModelOptions = struct('nTrials', 1000, 'useGPU', true, 'smartAllocation', false, 'verbose', false);
    fullModelOptions = struct('nTrials', 5000, 'useGPU', true, 'smartAllocation', false, 'verbose', false);
    
    % Fit for each block
    for block = nBlocks
        xBaseline = rmnan(xBaselineAll(block, :));
        yBaseline = rmnan(yBaselineAll(block, :));
        
        xOpto = rmnan(xOptoAll(block, :));
        yOpto = rmnan(yOptoAll(block, :));
        
        % Skip blocks with insufficient data
        if length(xBaseline) < 3 || length(xOpto) < 3
            fprintf('Skipping block %d due to insufficient data\n', block);
            continue;
        end
        
        fprintf('\n==== FITTING BLOCK %d/%d ====\n', block, nBlocks);
        
        %{
        % Optional: Perform quick parameter sensitivity analysis with fewer trials
        [sensitivity, paramRanking] = analyzeParameterSensitivity(xBaseline, yBaseline, xOpto, yOpto, initialParams, lb, ub, initialModelOptions);
        
        % Display sensitivity analysis results
        fprintf('Parameters ranked by sensitivity (most to least):\n');
        paramNames = {'a', 'b', 'l', 'w', 'g0', 'n', 'rmx', 'e', 'o'};
        for i = 1:length(paramRanking)
            fprintf('  %d. %s\n', i, paramNames{paramRanking(i)});
        end
        %}
        
        [bestParams, fitInfo] = multiStageOptimization(xBaseline, yBaseline, xOpto, yOpto, initialParams, lb, ub, fullModelOptions, ...
            'Verbose', true);
        
        % Store results - Use proper padding approach
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
        
        % Store model functions
        mdl.mdlBaseline = fitInfo.mdlBaseline;
        mdl.mdlOpto = fitInfo.mdlOpto;

        % Store data in mdl structure
        mdl.xBaseline(block, :) = padded_xBaseline;
        mdl.yBaseline(block, :) = padded_yBaseline;
        mdl.xOpto(block, :) = padded_xOpto;
        mdl.yOpto(block, :) = padded_yOpto;
        
        % Save parameter values with AICc
        mdl.fittedParams(block, :, 1) = [bestParams, fitInfo.aicc];
        mdl.r2Baseline(block) = fitInfo.rSquared_baseline;
        mdl.r2Opto(block) = fitInfo.rSquared_opto;
        
        % Store additional fitting information
        mdl.fitInfoByBlock{block} = fitInfo;
        
        fprintf('Fitted block #%d, AICc = %.2f, R² (baseline) = %.4f, R² (opto) = %.4f\n', ...
            block, fitInfo.aicc, fitInfo.rSquared_baseline, fitInfo.rSquared_opto);
    end
    
    % Optionally fit average data across all blocks
    if ~isempty(xBaselineAverage) && ~isempty(xOptoAverage) && nargout > 1
        fprintf('\n==== FITTING AVERAGE DATA ACROSS BLOCKS ====\n');
        
        [avgParams, avgFitInfo] = multiStageOptimization(xBaselineAverage, yBaselineAverage, xOptoAverage, yOptoAverage, initialParams, lb, ub, fullModelOptions, ...
            'Verbose', true);
        
        mdlAvg.headers = mdl.headers;
        mdlAvg.fittedParams(1, :, 1) = [avgParams, avgFitInfo.aicc];
        mdlAvg.r2Baseline = avgFitInfo.rSquared_baseline;
        mdlAvg.r2Opto = avgFitInfo.rSquared_opto;
        mdlAvg.fitInfo = avgFitInfo;
    end
    
    % Display summary results
    fprintf('\nFitting Results Summary:\n');
    displayTable(mdl.fittedParams(:, :, 1), mdl.headers);
    fprintf('\n=========\n\n');
end

% --- Helper Functions Below ---

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

function y = rmnan(x)
    y = x;
    y(isnan(y)) = [];
end