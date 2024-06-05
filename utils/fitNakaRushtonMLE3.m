function mdl = fitNakaRushtonMLE3(xBlocks, yBlocks, modelType)
% Define the number of blocks and conditions
[nConditions, nContrasts, nBlocks] = size(xBlocks);
%result = zeros(nBlocks, 5); % For Beta, exponent, and C50 for each condition and block
mdl=[];
%% Model
baselineModelFlag=0;

fprintf('Fitting model, independent parameters per session...')
for block = 1:nBlocks
    % Extract and sort data for baseline
    [xBaseline, sortIdx] = sort(rmnan(squeeze(xBlocks(1,:,block))));
    yBaseline = rmnan(squeeze(yBlocks(1,:,block))); yBaseline=yBaseline(sortIdx);
    numVal=numel(xBaseline);
    yBaselineAvg=mean([fliplr(100-yBaseline(1:numVal/2)); yBaseline(numVal/2+1:end)]); % average both sides
    yBaseline=[100-fliplr(yBaselineAvg), yBaselineAvg];  % average both sides and mirror
    [xBaseline, yBaseline] = mergeZeros(xBaseline, yBaseline);
    
    % Extract and sort data for horizontal (changed into con and incon)
    [xHorizontal,sortIdx]=sort(-1*rmnan(squeeze(xBlocks(2,:,block)))); 
    yHorizontal = rmnan(squeeze(yBlocks(2,:,block))); yHorizontal=100-yHorizontal(sortIdx);
    % Extract and sort data for vertical (changed into con and incon)
    [xVertical, sortIdx] = sort(rmnan(squeeze(xBlocks(3,:,block))));
    yVertical = rmnan(squeeze(yBlocks(3,:,block))); yVertical=yVertical(sortIdx);
    
    % Merge
    xOpto=mean([xVertical;xHorizontal]);
    yOpto=mean([yVertical;yHorizontal]);
    [xOpto, yOpto] = mergeZeros(xOpto, yOpto);

    % Define a common parameter vector [n, C50, Rmax, beta, delta] 
    % Assuming delta (lateral shift) is a new shared parameter
    [initialParams, lb, ub] = defineInitialParams(xHorizontal, yHorizontal, modelType);

    % Unified objective function
    options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'off', 'TolFun', 1e-14, 'TolX', 1e-14);
    mdl = calculateTotalError(mdl, initialParams, xBaseline, yBaseline, xOpto, yOpto, lb, ub,...
        options, modelType, baselineModelFlag, block);
end
fprintf('Done!\n\n')
headers = {'Beta_{opto}', 'n_{Baseline}', 'C_{50}', '\DeltaX', 'n_{con-opto}', 'n_{incon-opto}', 'AICc'};
displayTable(mdl.fittedParams(:,:,1), headers)
fprintf('\n=========\n\n')

%% Reference model (same parameters across sessions)
%{
initialParamsRef=nanmean(mdl.result(:,1:end-1),1);
fprintf('Fitting model, median of parameters across sessions...')
if baselineModelFlag
    for block = 1:nBlocks
        % Extract data for each condition
        xBaseline = rmnan(squeeze(xBlocks(1,:,block)));
        yBaseline = rmnan(squeeze(yBlocks(1,:,block)));
        xHorizontal = rmnan(squeeze(xBlocks(2,:,block)));
        yHorizontal = rmnan(squeeze(yBlocks(2,:,block)));
        xVertical = rmnan(squeeze(xBlocks(3,:,block)));
        yVertical = rmnan(squeeze(yBlocks(3,:,block)));
    
        % Define a common parameter vector [n, C50, Rmax, beta, delta] 
        % Assuming delta (lateral shift) is a new shared parameter
        [~, lb, ub] = defineInitialParams(xHorizontal, yHorizontal, modelType);
    
        % Unified objective function
        options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'off', 'TolFun', 1e-14, 'TolX', 1e-14);
        mdl = calculateTotalError(mdl, initialParamsRef, xBaseline, yBaseline, xHorizontal, yHorizontal, xVertical, yVertical, lb, ub,...
            options, modelType, baselineModelFlag, block);
    end
    fprintf('Done!\n\n')
    headers = {'Beta_{opto}', 'n_{Baseline}', 'C_{50}', '\DeltaX', 'n_{con-opto}', 'n_{incon-opto}', 'AICc'};
    displayTable(mdl.result(:,:,2), headers)
    fprintf('\n=========\n\n')
end
%}
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

function [initialParams, lb, ub] = defineInitialParams(x, y, hypothesisType)
    % Beta
    initBeta = min(y);
    % Exponent
    initExp = 3; % mean([2:4]); % Assuming an intermediate exponent value
    % Rmax
    initRmax=max(y) - min(y);
    % C50
    tmp = abs(y - (initRmax / 2 + min(y)));
    [~, closestIdx] = min(tmp);
    midpointIdx = unique([max(1, closestIdx-1), closestIdx, min(length(x), closestIdx+1)]);
    initC50 = 10; %mean(x(midpointIdx));
    % Delta-x (translation)
    initDeltaX = 0;

    % Duplicate params to allow free fitting, depending on hypothesis
    switch hypothesisType
        case 'base' 
            % Compile into initial fitting params
            initialParams = [NaN, initExp, initC50, NaN, initExp, initExp];
        
            % Define bounds for the params
            lb =  [    0,       0,         5,      -30,     0,       0]; % Lower bounds
            ub = [100,       10,     100,       30      10,     10]; % Upper bounds

        case 'beta' 
            % Compile into initial fitting params
            initialParams = [initBeta, initExp, initC50, NaN, NaN, NaN];
        
            % Define bounds for the params
            lb =  [    0,       2,         5,      -30,     0,     0]; % Lower bounds
            ub = [100,       6,     100,       30,     1,     1]; % Upper bounds

        case 'deltax' % translation and thus free delta-x (other params shared)
            % Compile into initial fitting params
            initialParams = [NaN, initExp, initC50, initDeltaX, NaN, NaN]; % Example: n, C50, Rmax, beta, delta
        
            % Define bounds for the params
            lb =  [    0,       2,         5,      -30,      0,     0]; % Lower bounds
            ub = [100,       6,     100,       30,      1,     1]; % Upper bounds
    end
end

function mdl = calculateTotalError(mdl, params, xBaseline, yBaseline, xOpto, yOpto, lb, ub,...
    options, modelType, baselineModelFlag, block)
    % Total trials and success counts
    [sumBaseline, successBaseline] = convert2counts(xBaseline, yBaseline);
    [sumOpto, successOpto] = convert2counts(xOpto, yOpto);

    % Define Naka-Rushton functions adjusted for potential lateral shift
    betaBL=50;
    % Duplicate params to allow free fitting, depending on hypothesis
    switch modelType
        case 'base'
            % Params = [beta exp C50 NaN]
            % other params shared, beta=50
            nakaRushtonBaseline = @(x, params) ((100-betaBL) .* (x).^params(2)) ./ (params(3).^params(2) + (x).^params(2)) + betaBL;
            nakaRushtonCon = @(x, params) ((100-betaBL) .* (x).^params(5)) ./ (params(3).^params(5) + (x).^params(5)) + betaBL;
            nakaRushtonIncon = @(x, params) (betaBL .* (x).^params(6)) ./ (params(3).^params(6) + (x).^params(6)) + (100-betaBL);
            
        case 'beta'
            % Params = [beta exp C50 NaN]
            % other params shared, beta=50
            %{
            nakaRushtonBaseline = @(x, params) ((100-betaBL) .* (x).^params(2)) ./ (params(3).^params(2) + (x).^params(2)) + betaBL;
            nakaRushtonCon = @(x, params) ((100-params(1)) .* (x).^params(2)) ./ (params(3).^params(2) + (x).^params(2)) + params(1);
            nakaRushtonIncon = @(x, params) (params(1) .* (x).^params(2)) ./ (params(3).^params(2) + (x).^params(2)) + (100-params(1));
            %}
            
            minBound = 1e-10;
            maxBound = 100 - minBound;
            betaBL = 50; % Example value, adjust as needed
            
            nakaRushtonBaseline = @(x, params) max(minBound, min(maxBound, (sign(x) .* ((100 - betaBL) .* abs(x).^params(2)) ./ (params(3).^params(2) + abs(x).^params(2))) + betaBL));
            %nakaRushtonOpto = @(x, params) max(minBound, min(maxBound, (sign(x) .* ((100 - params(1)) .* abs(x).^params(2)) ./ (params(3).^params(2) + abs(x).^params(2))) + params(1)));
            
            nakaRushtonOpto = @(x, params) (...
                ((x~=0) + 0.5 * (x==0)) .* ...
                (x<=0) .* max(minBound, min(maxBound, (sign(x) .* ((params(1)) .* abs(x).^params(2)) ./ (params(3).^params(2) + abs(x).^params(2))) + params(1))) + ...
                (x>=0) .* max(minBound, min(maxBound, (sign(x) .* ((100-params(1)) .* abs(x).^params(2)) ./ (params(3).^params(2) + abs(x).^params(2))) + params(1)))...
                );

        case 'deltax'
            minBound = 1e-10;
            maxBound = 100 - minBound;
            betaBL = 50; % Example value, adjust as needed
            
            nakaRushtonBaseline = @(x, params) max(minBound, min(maxBound, (sign(x) .* ((100 - betaBL) .* abs(x).^params(2)) ./ (params(3).^params(2) + abs(x).^params(2))) + betaBL));
            %nakaRushtonOpto = @(x, params) max(minBound, min(maxBound, (sign(x + params(4)) .* ((100 - betaBL) .* abs(x + params(4)).^params(2)) ./ (params(3).^params(2) + abs(x + params(4)).^params(2))) + betaBL));
            nakaRushtonOpto = @(x, params) (...
                ((x~=0) + 0.5 * (x==0)) .* ...
                (x<=0) .* max(minBound, min(maxBound, (sign(x + params(4)) .* ((betaBL) .* abs(x + params(4)).^params(2)) ./ (params(3).^params(2) + abs(x + params(4)).^params(2))) + betaBL)) + ...
                (x>=0) .* max(minBound, min(maxBound, (sign(x + params(4)) .* ((100 - betaBL) .* abs(x + params(4)).^params(2)) ./ (params(3).^params(2) + abs(x + params(4)).^params(2))) + betaBL))...
                );

    end

    % Calculate errors for all conditions
    objectiveFunction = @(params) ...
            -sum(successBaseline .* log(nakaRushtonBaseline(xBaseline, params)/100) + ...   % baseline neg-LogLikelihood
                       (sumBaseline - successBaseline) .* log(1 - nakaRushtonBaseline(xBaseline, params)/100)) + ...
            -sum(successOpto .* log(nakaRushtonOpto(xOpto, params)/100) + ...                         % con neg-LogLikelihood
                    (sumOpto - successOpto) .* log(1 - nakaRushtonOpto(xOpto, params)/100));

    % Testing
    %scatter(xBaseline, nakaRushtonBaseline(xBaseline,params)); xlim([-50 50]); ylim([0 150])

    if baselineModelFlag
        % Calculate nLL
        fittedParams=params;
        nLL=objectiveFunction(fittedParams); 
    else
        % Solve the optimization problem
        [fittedParams, nLL] = fmincon(objectiveFunction, params, [], [], [], [], lb, ub, [], options);
    end
    % Sum of errors from all conditions
    n = sum([sumBaseline sumOpto]);  % Total number of data points
    k = sum(~isnan(fittedParams));  % Number of parameters in the model
    [~, aicc, ~] = calculateAIC(nLL, k, n);

    % Save the input, equations, output 
    mdl.xBaseline(block,:)=padArray(xBaseline, 12, 2, NaN);
    mdl.yBaseline(block,:)=padArray(yBaseline, 12, 2, NaN);
    mdl.xOpto(block,:)=padArray(xOpto, 12, 2, NaN);
    mdl.yOpto(block,:)=padArray(yOpto, 12, 2, NaN);
    mdl.nakaRushtonBaseline=nakaRushtonBaseline;
    mdl.nakaRushtonOpto=nakaRushtonOpto;
    mdl.objectiveFunction=objectiveFunction;
    mdl.ub=ub;
    mdl.lb=lb;
    mdl.params(block,:)=params;
    mdl.options=options;
    mdl.fittedParams(block,:,baselineModelFlag+1)=[fittedParams aicc];
end

function [sumY, successY] = convert2counts(x,y)
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
