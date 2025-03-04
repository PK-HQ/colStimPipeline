function [mdl,mdlAvg]= fitNakaRushtonMLE3(xBlocks, yBlocks, modelType)
% Define the number of blocks and conditions
[nConditions, nContrasts, nBlocks] = size(xBlocks);
%result = zeros(nBlocks, 5); % For Beta, exponent, and C50 for each condition and block
mdl=[];mdlAvg=[];
%% Model
baselineModelFlag=0;
% Process data across all blocks, get average

[xBaselineAll, yBaselineAll, xOptoAll, yOptoAll] = processConditionsBlocks(xBlocks, yBlocks);
%{
[xBaselineAverage, yBaselineAverage] = binData3(xBaselineAll, yBaselineAll);
[xOptoAverage, yOptoAverage] = binData3(xOptoAll, yOptoAll);
%}
xBaselineAverage=rmnan(reshape(xBaselineAll,1,numel(xBaselineAll)));
yBaselineAverage=rmnan(reshape(yBaselineAll,1,numel(yBaselineAll)));
xOptoAverage=rmnan(reshape(xOptoAll,1,numel(xOptoAll)));
yOptoAverage=rmnan(reshape(yOptoAll,1,numel(yOptoAll)));
% Solver params
options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'off', 'TolFun', 1e-14, 'TolX', 1e-14);
%options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off', 'TolFun', 1e-14, 'TolX', 1e-14);
fprintf('Fitting model, independent parameters per session...')
for block = 1%:3%round(nBlocks/2)
    xBaseline=rmnan(xBaselineAll(block,:));
    yBaseline=rmnan(yBaselineAll(block,:));
    
    xOpto=rmnan(xOptoAll(block,:));
    yOpto=rmnan(yOptoAll(block,:));
    % Define a common parameter vector [n, C50, Rmax, beta, delta] 
    % Assuming delta (lateral shift) is a new shared parameter
    
    [initialParams, lb, ub, mdl.headers] = defineInitialParams(xOpto(xOpto<0), yOpto(yOpto<0), modelType);

    lb = max(lb, eps); % Set lower bounds to a small positive value to avoid zero
    ub(isinf(ub)) = 1e10; % Set a practical upper limit if infinity is used

    % Fitting for block
    % Create optimization options with increased iterations and improved tolerance
    options = optimset('fminsearch');
    options = optimset(options, 'MaxIter', 100000);      % Increase maximum iterations (default: 200*numberOfVariables)
    options = optimset(options, 'MaxFunEvals', 20000);  % Increase function evaluations (default: 200*numberOfVariables)
    options = optimset(options, 'TolX', 1e-6);          % Tighter tolerance on parameter values
    options = optimset(options, 'TolFun', 1e-6);        % Tighter tolerance on function value
    %options = optimset(options, 'Display', 'iter');     % Show iteration progress
    
    mdl = calculateTotalError(mdl, initialParams, xBaseline, yBaseline, xOpto, yOpto, lb, ub,...
        options, modelType, baselineModelFlag, block);
end
%{
% Fitting for binned average
opts=[];
%opts = optimset('fminsearch');
%opts.TolX = 1e-8;  % Tolerance for parameter updates
%opts.TolFun = 1e-8; % Tolerance for objective function value changes
[initialParams, lb, ub, mdlAvg.headers] = defineInitialParams(xOptoAverage(xOpto<0), yOptoAverage(yOpto<0), modelType);
%Options for fitting
% Customize the options for a thorough search
options = optimset('fminsearch');
options.Display = 'off';        % Full verbose output
options.MaxIter = 100000;        % Maximum number of iterations
options.MaxFunEvals = 200000;    % Maximum number of function evaluations
options.TolFun = 1e-8;           % Function tolerance for convergence
options.TolX = 1e-8;             % Parameter tolerance for convergence
mdlAvg = calculateTotalError(mdlAvg, initialParams, xBaselineAverage, yBaselineAverage, xOptoAverage, yOptoAverage, lb, ub,...
    options, modelType, baselineModelFlag, 1);
%}
fprintf('Done!\n\n')
displayTable(mdl.fittedParams(:,:,1), mdl.headers(1:end))%mdl.headers(1:end-4))
%displayTable(mdlAvg.fittedParams(:,:,1), mdl.headers(1:end))%mdlAvg.headers(1:end-4))
fprintf('\n=========\n\n')

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

function [xBaseline, yBaseline, xOpto, yOpto] = processConditions(xBlocks, yBlocks, block)
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
    
end
% [xAverage, yAverage, counts] = binData(xOpto, yOpto)

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

function [initialParams, lb, ub, headers] = defineInitialParams(x, y, hypothesisType)
    % Beta
    initBeta = 50;%min(y);
    % Exponent
    initExp = 3; % mean([2:4]); % Assuming an intermediate exponent value
    % Rmax
    initRmax=50;
    % C50
    tmp = abs(y - (initRmax / 2 + min(y)));
    [~, closestIdx] = min(tmp);
    initC50 = 15; %mean(x(midpointIdx));
    % Delta-x (translation)
    initDeltaX = 0;

    % Duplicate params to allow free fitting, depending on hypothesis
    switch hypothesisType
        case 'base' 
            % Compile into initial fitting params
            initialParams = [NaN, initExp, initC50, NaN, initExp, initExp];
        
            % Define bounds for the params
            lb =  [    0,       2,         5,      -30,     0,       0]; % Lower bounds
            ub = [100,       6,     100,       30      1,       1]; % Upper bounds
            headers = {'Beta_opto', 'n_Baseline', 'C_50', '\Delta_x', 'n_con-opto', 'n_incon-opto', 'AICc'};
        case 'weibull'
            headers = {'C50','Beta','lowerAsy','lapse','AICc'};
            initialParams = [15, 5, .5, 0.1];
            lb = [10 .1 .2 0];
            ub = [70 8 .8 .3];
        case 'weibullfreeL'
            % C50 - threshold contrast
            % Beta - slope
            % lowerAsy = lower asymptote, joins 2 curves at x=0
            % lapse = 1-lapse is the upper asymptote that controls the
            % extreme 
            headers = {'C50','Beta','lowerAsy','lapse1','lapse2','AICc'};
            initialParams = [15, 5, .5, 0.1, 0.1];
            lb = [10 .1 .2 0 0];
            ub = [70 8 .8 .3 .3];
        case 'weibullfreeLS'
            % C50 - threshold contrast
            % Beta - slope
            % lowerAsy = lower asymptote, joins 2 curves at x=0
            % lapse = 1-lapse is the upper asymptote that controls the
            % extreme 
            headers = {'C_{50}','Beta_{incon}','B_{con-incon}','A_{incon}',...
                'Beta_{con}','A_{con}','AICc'};
            initialParams = [15, 3, .5, 0.1,  5, 0.1];
            lb = [10 1 .2 0 1 0]; %lb = [10 .1 .2 0 .1 0];
            ub = [70 20 .8 .2 20 .2];
        case 'weibullfreeAll'
            % C50 - threshold contrast
            % Beta - slope
            % lowerAsy = lower asymptote, joins 2 curves at x=0
            % lapse = 1-lapse is the upper asymptote that controls the
            % extreme 
            headers = {'A^{bl}','B^{bl}','\alpha^{bl}','\beta^{bl}',...
                '\DeltaA^{con-bl}','\DeltaB^{con-bl}','\Delta\alpha^{con-bl}','\Delta\beta^{con-bl}',...
                '\DeltaA^{incon-bl}','\DeltaB^{incon-bl}','\Delta\alpha^{incon-bl}','\Delta\beta^{incon-bl}',...
                'AICc^{total}','AUC^{bl}','AUC^{con-bl}','AUC^{incon-bl}','AUC^{con-incon}'};
            initialParams = [.1 .5 15  3,...
                                      .1 .5 15 3,...
                                      .1 nan 15 3];
            lb = [0 .5 10 1,...
                    0 .2 10 1,...
                    0 nan 10 1]; 
            ub = [.2 .5 50 8,...
                     .2 .8 50 8,...
                     .2 nan 50 8];

        case 'weibull-beta'
            % C50 - threshold contrast
            % Beta - slope
            % lowerAsy = lower asymptote, joins 2 curves at x=0
            % lapse = 1-lapse is the upper asymptote that controls the
            % extreme 
            headers = {'A^{bl}','B^{bl}','C_{50}^{bl}','\beta^{bl}',...
                '\DeltaA^{con-bl}','\DeltaB^{con-bl}','\DeltaC_{50}^{con-bl}','\Delta\beta^{con-bl}',...
                '\DeltaA^{incon-bl}','\DeltaB^{incon-bl}','\DeltaC_{50}^{incon-bl}','\Delta\beta^{incon-bl}','AICc'};
            initialParams = [.1 .5 15  3,...
                                      .1 .3 20 3,...
                                      .1 .3 20 3];
            lb = [0 .5 10 1,...
                    0 .2 10 1,...
                    0 .2 10 1]; 
            ub = [.2 .5 50 8,...
                     .2 .8 50 8,...
                     .2 .8 50 8];

        case 'bill'
            % MODIFIED: Bill model with expanded parameter range and additional zero-contrast parameters
            % Initial parameters for the 'bill' model: [a, b, l, w, g0, n, rmx, e, o, zero_con, zero_incon]
            initialParams = [.5,.5,.1,.1, ... %a,  b,   l,   w
                50, 4, 5,... % g0,   n,   rmx
                .5, 10]; % e o
            
            % Parameters: [a,  b,   l,   w,    g0,   n,   rmx,   e,    o]
            %                |     |     |     |       |       |       |        |      | 
            %                |     |     |     |       |       |       |        |      |     
            %                |     |     |     |       |       |       |        |      |        
            %                |     |     |     |       |       |       |        |      Opto-stim effective contrast
            %                |     |     |     |       |       |       |        Relative weight of optostim effects on excitation
            %                |     |     |     |       |       |       Max response
            %                |     |     |     |       |      Spiking exponent
            %                |     |     |     |     Normalization constant
            %                |     |     |    Weight of ortho-tuned visual input norm
            %                |     |   Weight of iso-tuned visual input norm
            %                |   Weight of ortho-tuned visual input
            %              Weight of iso-tuned visual input

            % MODIFIED: Expanded bounds for the parameters
            lb = [1e-8,1e-8,1e-8,1e-8, ... %a,  b,   l,   w
                40, 4, 0,... % g0,   n,   rmx
                .5, 0]; % e o
            % Parameters: [a,  b,   l,   w,    g0,   n,   rmx,   e,    o]

            % MODIFIED: Expanded upper bounds 
            ub = [1,1,1,1, ... %a,  b,   l,   w
                100, 12, 5,... % g0,   n,   rmx
                .6, 100]; % e o
            % Parameters: [a,  b,   l,   w,    g0,   n,   rmx,   e,    o]
            
            % MODIFIED: Added zero_con and zero_incon to headers
            headers = {'a','b','l','w','g0','n','rmx','e','o', 'AICc'};

        case 'beta' 
            % Compile into initial fitting params
            initialParams = [initBeta, initExp, initC50, NaN, NaN, NaN];
        
            % Define bounds for the params
            lb =  [    0,       2,         5,      -30,     0,     0]; % Lower bounds
            ub = [100,       6,     100,       30,     1,     1]; % Upper bounds
            headers = {'Beta_opto', 'n_Baseline', 'C_50', '\Delta_x', 'n_conopto', 'n_inconopto', 'AICc'};

        case '\Deltax' % translation and thus free delta-x (other params shared)
            % Compile into initial fitting params
            initialParams = [NaN, initExp, initC50, initDeltaX, NaN, NaN]; % Example: n, C50, Rmax, beta, delta
        
            % Define bounds for the params
            lb =  [    0,       2,         5,      -30,      0,     0]; % Lower bounds
            ub = [100,       6,     100,       30,      1,     1]; % Upper bounds
            headers = {'Beta_opto', 'n_Baseline', 'C_50', '\Delta_x', 'n_con-opto', 'n_incon-opto', 'AICc'};

    end
end

function mdl = calculateTotalError(mdl, params, xBaseline, yBaseline, xOpto, yOpto, lb, ub,...
    opts, modelType, baselineModelFlag, block)
    % Total trials and success counts
    [sumBaseline, successBaseline] = convert2counts(xBaseline, yBaseline);
    [sumOpto, successOpto] = convert2counts(xOpto, yOpto);

    % Define Naka-Rushton functions adjusted for potential lateral shift
    switch modelType
        case 'base'
            % Params = [beta exp C50 NaN]
            % other params shared, beta=50
            minBound = 1e-10;
            maxBound = 100 - minBound;
            betaBL = 50; % Example value, adjust as needed
            
            mdlBaseline = @(x, params) max(minBound, min(maxBound, (sign(x) .* ((100 - betaBL) .* abs(x).^params(2)) ./ (params(3).^params(2) + abs(x).^params(2))) + betaBL));
            
            mdlOpto = @(x, params) (...
                ((x~=0) + 0.5 * (x==0)) .* ...
                (x<=0) .* max(minBound, min(maxBound, (sign(x) .* ((params(1)) .* abs(x).^params(2)) ./ (params(3).^params(2) + abs(x).^params(2))) + params(1))) + ...
                (x>=0) .* max(minBound, min(maxBound, (sign(x) .* ((100-params(1)) .* abs(x).^params(2)) ./ (params(3).^params(2) + abs(x).^params(2))) + params(1)))...
                );
            
            % Calculate errors for all conditions
            objectiveFunction = @(params) ...
                    -sum(successBaseline .* log(mdlBaseline(xBaseline, params)/100) + ...   % baseline neg-LogLikelihood
                               (sumBaseline - successBaseline) .* log(1 - mdlBaseline(xBaseline, params)/100));
        case 'weibull'
            % WEIBULL MODEL
            % Params = [Alpha Beta Gamma Lapse]
            lowerAsymptoteBL=.5;
            % Weibull baseline model: Y = (1-exp(-(X/Alpha).^Beta))*(Gamma-Lapse)+1-Gamma;
            mdlBaseline = @(x, params) 100*((x==0)*lowerAsymptoteBL) + ... % midpoint
                [fliplr(100-(100 .* ((x>=0) .* (sign(x) .* (1 - exp(-((x) ./ params(1)).^params(2))) * ((1-lowerAsymptoteBL) - params(4)) + 1 - (1-lowerAsymptoteBL))))) + ...
                ((100 .* (x>=0) .* (sign(x) .* (1 - exp(-((x) ./ params(1)).^params(2))) * (lowerAsymptoteBL - params(4)) + 1 - lowerAsymptoteBL))-100.* (x>=0))];
            
            % Weibull  model
            mdlOpto = @(x, params) 100*((x==0)*params(3)) + ... % midpoint
                [fliplr(100-(100 .* ((x>=0) .* (sign(x) .* (1 - exp(-((x) ./ params(1)).^params(2))) * ((1-params(3)) - params(4)) + 1 - (1-params(3)))))) + ... % x<=0
                ((100 .* (x>=0) .* (sign(x) .* (1 - exp(-((x) ./ params(1)).^params(2))) * (params(3) - params(4)) + 1 - params(3)))-100.* (x>=0))]; % x>=0
            
            % Objective function for Weibull model
            objectiveFunction = @(params) ...
                -sum(successBaseline .* log(mdlBaseline(xBaseline, params)/100) + ...   % Baseline neg-LogLikelihood
                       (sumBaseline - successBaseline) .* log(1 - mdlBaseline(xBaseline, params)/100)) + ...
                -sum(successOpto .* log(mdlOpto(xOpto, params)/100) + ...               % Opto neg-LogLikelihood
                       (sumOpto - successOpto) .* log(1 - mdlOpto(xOpto, params)/100));                        
        case 'weibullfreeL'
            % WEIBULL MODEL
            % Params = [Alpha Beta Gamma Lapse1 Lapse2]
            lowerAsymptoteBL=.5;
            % Weibull baseline model, upper asymptote is mean of opto
            % upper asymptotes
            mdlBaseline = @(x, params) 100*((x==0)*lowerAsymptoteBL) + ... % midpoint
                [fliplr(100-(100 .* ((x>=0) .* (sign(x) .* (1 - exp(-((x) ./ params(1)).^params(2))) * ((1-lowerAsymptoteBL) - mean([params(4) params(5)])) + 1 - (1-lowerAsymptoteBL))))) + ...
                ((100 .* (x>=0) .* (sign(x) .* (1 - exp(-((x) ./ params(1)).^params(2))) * (lowerAsymptoteBL - mean([params(4) params(5)])) + 1 - lowerAsymptoteBL))-100.* (x>=0))];
            
            % Weibull opto model, independent upper asymptotes
            mdlOpto = @(x, params) 100*((x==0)*params(3)) + ... % midpoint
                [fliplr(100-(100 .* ((x>=0) .* (sign(x) .* (1 - exp(-((x) ./ params(1)).^params(2))) * ((1-params(3)) - params(4)) + 1 - (1-params(3)))))) + ... % x<=0
                ((100 .* (x>=0) .* (sign(x) .* (1 - exp(-((x) ./ params(1)).^params(2))) * (params(3) - params(5)) + 1 - params(3)))-100.* (x>=0))]; % x>=0
            
            % Objective function for Weibull model
            objectiveFunction = @(params) ...
                -sum(successBaseline .* log(mdlBaseline(xBaseline, params)/100) + ...   % Baseline neg-LogLikelihood
                       (sumBaseline - successBaseline) .* log(1 - mdlBaseline(xBaseline, params)/100)) + ...
                -sum(successOpto .* log(mdlOpto(xOpto, params)/100) + ...               % Opto neg-LogLikelihood
                       (sumOpto - successOpto) .* log(1 - mdlOpto(xOpto, params)/100));                        
        case 'weibullfreeLS'
            % WEIBULL MODEL
            % Params = [Alpha Beta1 Gamma Lapse1 Beta1 Lapse2]
            lowerAsymptoteBL=.5;
            % Weibull baseline model, upper asymptote is mean of opto
            % upper asymptotes
            mdlBaseline = @(x, params) 100*((x==0)*lowerAsymptoteBL) + ... % midpoint
                [fliplr(100-(100 .* ((x>=0) .* (sign(x) .* (1 - exp(-((x) ./ params(1)).^mean([params(2) params(5)]))) * ((1-lowerAsymptoteBL) - mean([params(4) params(6)])) + 1 - (1-lowerAsymptoteBL))))) + ...
                ((100 .* (x>=0) .* (sign(x) .* (1 - exp(-((x) ./ params(1)).^mean([params(2) params(5)]))) * (lowerAsymptoteBL - mean([params(4) params(6)])) + 1 - lowerAsymptoteBL))-100.* (x>=0))];
            
            % Weibull opto model, independent upper asymptotes
            mdlOpto = @(x, params) 100*((x==0)*params(3)) + ... % midpoint
                [fliplr(100-(100 .* ((x>=0) .* (sign(x) .* (1 - exp(-((x) ./ params(1)).^params(2))) * ((1-params(3)) - params(4)) + 1 - (1-params(3)))))) + ... % x<=0
                ((100 .* (x>=0) .* (sign(x) .* (1 - exp(-((x) ./ params(1)).^params(5))) * (params(3) - params(6)) + 1 - params(3)))-100.* (x>=0))]; % x>=0
            
            % Objective function for Weibull model
            objectiveFunction = @(params) ...
                -sum(successBaseline .* log(mdlBaseline(xBaseline, params)/100) + ...   % Baseline neg-LogLikelihood
                       (sumBaseline - successBaseline) .* log(1 - mdlBaseline(xBaseline, params)/100)) + ...
                -sum(successOpto .* log(mdlOpto(xOpto, params)/100) + ...               % Opto neg-LogLikelihood
                       (sumOpto - successOpto) .* log(1 - mdlOpto(xOpto, params)/100));                    
        case 'weibullfreeAll'
            % WEIBULL MODEL
            % Params = [Alpha Beta1 Gamma Lapse1 Beta1 Lapse2]
            % Weibull baseline model: Y = (1-exp(-(X/C50).^Beta))*(Gamma-Lapse)+1-Gamma;

            lowerAsymptoteBL=.5;
            % Weibull baseline model, upper asymptote is mean of opto
            % upper asymptotes
            mdlBaseline = @(x, params) 100*((x==0)*lowerAsymptoteBL) + ... % midpoint
                [fliplr(100-(100 .* ((x>=0) .* (sign(x) .* (1 - exp(-((x) ./ params(3)).^params(4))) * ((1-lowerAsymptoteBL) - params(1)) + 1 - (1-lowerAsymptoteBL))))) + ...
                ((100 .* (x>=0) .* (sign(x) .* (1 - exp(-((x) ./ params(3)).^params(4))) * (lowerAsymptoteBL - params(1)) + 1 - lowerAsymptoteBL))-100.* (x>=0))];
            
            % Weibull opto model, independent upper asymptotes
            mdlOpto = @(x, params) 100*((x==0)*params(6)) + ... % midpoint
                [fliplr(100-(100 .* ((x>=0) .* (sign(x) .* (1 - exp(-((x) ./ params(7)).^params(8))) * ((1-params(6)) - params(5)) + 1 - (1-params(6)))))) + ... % x<=0
                ((100 .* (x>=0) .* (sign(x) .* (1 - exp(-((x) ./ params(11)).^params(12))) * (params(6) - params(9)) + 1 - params(6)))-100.* (x>=0))]; % x>=0
            
            epsilon = 1e-8; % Small value to avoid log(0)
            
            objectiveFunction = @(params) ...
                -sum(successBaseline .* log((mdlBaseline(xBaseline, params) + epsilon)/100) + ... 
                     (sumBaseline - successBaseline) .* log(1 - (mdlBaseline(xBaseline, params) + epsilon)/100)) + ...
                -sum(successOpto .* log((mdlOpto(xOpto, params) + epsilon)/100) + ...
                     (sumOpto - successOpto) .* log(1 - (mdlOpto(xOpto, params) + epsilon)/100));
            
        case 'weibull-beta'
            % WEIBULL MODEL
            % Params = [Alpha Beta1 Gamma Lapse1 Beta1 Lapse2]
            % Weibull baseline model: Y = (1-exp(-(X/C50).^Beta))*(Gamma-Lapse)+1-Gamma;

            lowerAsymptoteBL=.5;
            % Weibull baseline model, upper asymptote is mean of opto
            %{
            % upper asymptotes
            mdlBaseline = @(x, params) 100*((x==0)*lowerAsymptoteBL) + ... % midpoint
                [fliplr(100-(100 .* ((x>=0) .* (sign(x) .* (1 - exp(-((x) ./ params(3)).^params(4))) * ((1-lowerAsymptoteBL) - params(1)) + 1 - (1-lowerAsymptoteBL))))) + ...
                ((100 .* (x>=0) .* (sign(x) .* (1 - exp(-((x) ./ params(3)).^params(4))) * (lowerAsymptoteBL - params(1)) + 1 - lowerAsymptoteBL))-100.* (x>=0))];
            
            % Weibull opto model, independent upper asymptotes
            mdlOpto = @(x, params) 100*((x==0)*params(6)) + ... % midpoint
                [fliplr(100-(100 .* ((x>=0) .* (sign(x) .* (1 - exp(-((x) ./ params(7)).^params(8))) * ((1-params(6)) - params(5)) + 1 - (1-params(6)))))) + ... % x<=0
                ((100 .* (x>=0) .* (sign(x) .* (1 - exp(-((x) ./ params(11)).^params(8))) * (params(6) - params(9)) + 1 - params(6)))-100.* (x>=0))]; % x>=0
            
            % Objective function for Weibull model
            objectiveFunction = @(params) ...
                -sum(successBaseline .* log(mdlBaseline(xBaseline, params)/100) + ...   % Baseline neg-LogLikelihood
                       (sumBaseline - successBaseline) .* log(1 - mdlBaseline(xBaseline, params)/100)) + ...
                -sum(successOpto .* log(mdlOpto(xOpto, params)/100) + ...               % Opto neg-LogLikelihood
                       (sumOpto - successOpto) .* log(1 - mdlOpto(xOpto, params)/100));         
            %}

            mdlBaseline = @(x, params) 100*((x==0)*lowerAsymptoteBL) + ... % midpoint
                [fliplr(100-(100 .* ((x>=0) .* (sign(x) .* (1 - exp(-((x) ./ params(3)).^params(4))) * ((1-lowerAsymptoteBL) - params(1)) + 1 - (1-lowerAsymptoteBL))))) + ...
                ((100 .* (x>=0) .* (sign(x) .* (1 - exp(-((x) ./ params(3)).^params(4))) * (lowerAsymptoteBL - params(1)) + 1 - lowerAsymptoteBL))-100.* (x>=0))];
            
            % Weibull opto model, independent upper asymptotes
            mdlOpto = @(x, params) 100*((x==0)*params(6)) + ... % midpoint
                [fliplr(100-(100 .* ((x>=0) .* (sign(x) .* (1 - exp(-((x) ./ params(7)).^params(8))) * ((1-params(6)) - params(5)) + 1 - (1-params(6)))))) + ... % x<=0
                ((100 .* (x>=0) .* (sign(x) .* (1 - exp(-((x) ./ params(11)).^params(8))) * (params(6) - params(9)) + 1 - params(6)))-100.* (x>=0))]; % x>=0
            
            epsilon = 1e-8; % Small value to avoid log(0)
            
            objectiveFunction = @(params) ...
                -sum(successBaseline .* log((mdlBaseline(xBaseline, params) + epsilon)/100) + ... 
                     (sumBaseline - successBaseline) .* log(1 - (mdlBaseline(xBaseline, params) + epsilon)/100)) + ...
                -sum(successOpto .* log((mdlOpto(xOpto, params) + epsilon)/100) + ...
                     (sumOpto - successOpto) .* log(1 - (mdlOpto(xOpto, params) + epsilon)/100));

        case 'bill'
            minBound = 1e-10;
            maxBound = 100 - minBound;
            % GPU/ParPool options
            options = struct('useGPU', true, 'useParallel', false, 'nTrials', 1000, 'verbose', false);
            %mdlBaseline = @(x, params) max(minBound, min(maxBound, normMdlSimOptimized(x, params,'baseline',options))); % zero opto input
            %mdlOpto = @(x, params) max(minBound, min(maxBound, normMdlSimOptimized(x, params,'opto',options))); % nonzero opto input
            
            mdlBaseline = @(x, params) max(minBound, min(maxBound, normMdlSimOptimized(x, params,'baseline'))); % zero opto input
            mdlOpto = @(x, params) max(minBound, min(maxBound, normMdlSimOptimized(x, params,'opto'))); % nonzero opto input
            
            % Define NLL objective function
            objectiveFunction = @(params) ...
                -sum(successBaseline .* log(mdlBaseline(xBaseline, params) / 100) + ...
                     (sumBaseline - successBaseline) .* log(1 - mdlBaseline(xBaseline, params) / 100)) + ...
                -sum(successOpto .* log(mdlOpto(xOpto, params) / 100) + ...
                     (sumOpto - successOpto) .* log(1 - mdlOpto(xOpto, params) / 100));

        case 'beta'
            % Params = [beta exp C50 NaN]
            % other params shared, beta=50
           
            minBound = 1e-10;
            maxBound = 100 - minBound;
            betaBL = 50; % Example value, adjust as needed
            
            mdlBaseline = @(x, params) max(minBound, min(maxBound, (sign(x) .* ((100 - betaBL) .* abs(x).^params(2)) ./ (params(3).^params(2) + abs(x).^params(2))) + betaBL));
            %mdlOpto = @(x, params) max(minBound, min(maxBound, (sign(x) .* ((100 - params(1)) .* abs(x).^params(2)) ./ (params(3).^params(2) + abs(x).^params(2))) + params(1)));
            
            mdlOpto = @(x, params) (...
                ((x~=0) + 0.5 * (x==0)) .* ...
                (x<=0) .* max(minBound, min(maxBound, (sign(x) .* ((params(1)) .* abs(x).^params(2)) ./ (params(3).^params(2) + abs(x).^params(2))) + params(1))) + ...
                (x>=0) .* max(minBound, min(maxBound, (sign(x) .* ((100-params(1)) .* abs(x).^params(2)) ./ (params(3).^params(2) + abs(x).^params(2))) + params(1)))...
                );

            % Calculate errors for all conditions
            objectiveFunction = @(params) ...
                    -sum(successBaseline .* log(mdlBaseline(xBaseline, params)/100) + ...   % baseline neg-LogLikelihood
                               (sumBaseline - successBaseline) .* log(1 - mdlBaseline(xBaseline, params)/100)) + ...
                    -sum(successOpto .* log(mdlOpto(xOpto, params)/100) + ...                         % con neg-LogLikelihood`
                            (sumOpto - successOpto) .* log(1 - mdlOpto(xOpto, params)/100));
                        
            %{
            mdlBaseline = @(x, params) ((100-betaBL) .* (x).^params(2)) ./ (params(3).^params(2) + (x).^params(2)) + betaBL;
            nakaRushtonCon = @(x, params) ((100-params(1)) .* (x).^params(2)) ./ (params(3).^params(2) + (x).^params(2)) + params(1);
            nakaRushtonIncon = @(x, params) (params(1) .* (x).^params(2)) ./ (params(3).^params(2) + (x).^params(2)) + (100-params(1));
            %}

        case '\Deltax'
            minBound = 1e-10;
            maxBound = 100 - minBound;
            betaBL = 50; % Example value, adjust as needed
            
            mdlBaseline = @(x, params) max(minBound, min(maxBound, (sign(x) .* ((100 - betaBL) .* abs(x).^params(2)) ./ (params(3).^params(2) + abs(x).^params(2))) + betaBL));
            %mdlOpto = @(x, params) max(minBound, min(maxBound, (sign(x + params(4)) .* ((100 - betaBL) .* abs(x + params(4)).^params(2)) ./ (params(3).^params(2) + abs(x + params(4)).^params(2))) + betaBL));
            mdlOpto = @(x, params) (...
                ((x~=0) + 0.5 * (x==0)) .* ...
                (x<=0) .* max(minBound, min(maxBound, (sign(x + params(4)) .* ((betaBL) .* abs(x + params(4)).^params(2)) ./ (params(3).^params(2) + abs(x + params(4)).^params(2))) + betaBL)) + ...
                (x>=0) .* max(minBound, min(maxBound, (sign(x + params(4)) .* ((100 - betaBL) .* abs(x + params(4)).^params(2)) ./ (params(3).^params(2) + abs(x + params(4)).^params(2))) + betaBL))...
                );

            % Calculate errors for all conditions
            objectiveFunction = @(params) ...
                    -sum(successBaseline .* log(mdlBaseline(xBaseline, params)/100) + ...   % baseline neg-LogLikelihood
                               (sumBaseline - successBaseline) .* log(1 - mdlBaseline(xBaseline, params)/100)) + ...
                    -sum(successOpto .* log(mdlOpto(xOpto, params)/100) + ...                         % con neg-LogLikelihood
                            (sumOpto - successOpto) .* log(1 - mdlOpto(xOpto, params)/100));
    end

    % Solve
    solverType='FMSB';

    if baselineModelFlag
        % Calculate nLL
        fittedParams=params;
        nLL=objectiveFunction(fittedParams); 
    else
        switch solverType
            case 'fmincon'
                % Solve the optimization problem
                % Test the objective function value for the initial guess
                initialObjectiveValue = objectiveFunction(params);
                [fittedParams, nLL] = fmincon(objectiveFunction, params, [], [], [], [], lb, ub, [], opts);
            case 'FMSB'
                [fittedParams]=fminsearchbnd(@(params) objectiveFunction(params),...
                    params,...
                    lb,...
                    ub, opts);
                nLL=objectiveFunction(fittedParams);
        end
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
    mdl.mdlBaseline=mdlBaseline;
    mdl.mdlOpto=mdlOpto;
    mdl.objectiveFunction=objectiveFunction;
    mdl.ub=ub;
    mdl.lb=lb;
    mdl.params(block,:)=params;
    mdl.options=opts;
    mdl.fittedParams(block,:,baselineModelFlag+1)=[fittedParams aicc];
    fprintf('Fitted block #%.0f',block)
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

