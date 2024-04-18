function params = fitNakaRushtonMLE(xBlocks, yBlocks, objFunc, hypothesisType, fixedParams)
    %% Fits the 3D matrices of x and y (3 conditions x 6 contrast levels x 16 blocks) with a 
    % Naka-Rushton function, and with MLE (binomial dist assumption)
    % Input: 
    %       x = contrast values (range 0:100%), 3 conditions x 6 contrast levels x 16 blocks
    %       y = percentage correct (range 0:100%), 3 conditions x 6 contrast levels x 16 blocks
    % Output:
    %       params = [β, n, C50]

    % Define the number of blocks and conditions
    [nConditions, nContrasts, nBlocks] = size(xBlocks);
    
    % Initialize output
    params = zeros(nBlocks, nConditions, 5); % For Beta, exponent, and C50 for each condition and block
    
    % Options for fmincon
    options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off', 'TolFun', 1e-14, 'TolX', 1e-14);
    
    % Loop over blocks
    for b = 1:nBlocks
        % Baseline condition fitting with fixed Beta and Rmax
        x = rmnan(squeeze(xBlocks(1,:,b))); % Contrast values for baseline
        y = rmnan(squeeze(yBlocks(1,:,b))); % Percentage correct for baseline
        params_bl = fitNakaRushtonBaseline(x, y, objFunc, options, hypothesisType, fixedParams);
        
        % Fit con-opto and incon-opto conditions simultaneously
        xCon = rmnan(squeeze(xBlocks(2,:,b))); % Contrast values for con-opto
        yCon = rmnan(squeeze(yBlocks(2,:,b))); % Percentage correct for con-opto
        xIncon = rmnan(squeeze(xBlocks(3,:,b))); % Contrast values for incon-opto
        yIncon = rmnan(squeeze(yBlocks(3,:,b))); % Percentage correct for incon-opto
        [params_con, params_incon] = fitNakaRushtonOptos(xCon, yCon, xIncon, yIncon, objFunc, options, hypothesisType, fixedParams);
        
        params(b,1,:)= params_bl;
        params(b,2,:) = params_con;
        params(b,3,:) = params_incon;
    end
end

  
%% Subfunctions
function params = fitNakaRushtonBaseline(x, y, objFunc, options, hypothesisType, fixedParams)
    % N trials
    N =  [(40*ones(1,sum(x==0))),...
        (20*ones(1,sum(x>0)))]; % trials per contrast level
    
    % Convert percent correct into success counts
    y_counts = (y / 100) .* N;

    % Define the fixed beta and Rmax
    beta = 50;
    Rmax = 50;
    
    % Fixed params
    if sum(~isnan(fixedParams))==0
        % Define the Naka-Rushton function with fitted params = [n, C50] (fixed Rmax, beta)
        nakaRushton = @(x, params) (Rmax .* x.^params(1)) ./ (params(2).^params(1) + x.^params(1)) + beta;
    else 
        initExp=fixedParams(1,2);
        % Define the Naka-Rushton function with fitted params = [n, C50, Rmax, beta] (beta is shared across both opto conditions)
        % Define the Naka-Rushton function with fitted params = [n, C50] (fixed Rmax, beta)
        nakaRushton = @(x, params) (Rmax .* x.^initExp) ./ (params(2).^initExp + x.^initExp) + beta;
    end

    % Initial guesses for the parameters [n, C50]
    initialParams = [2, 10];  % These should be adjusted based on your data
    % Constraints
    lb = [2,      5];  % Lower bounds: no negative values for exponent and C50
    ub = [6, 100];  % Upper bounds: allowing exponent and C50 to be any positive value
    
    % Objective function
    switch objFunc
        case {'SSE'}
            objectiveFunction = @(params) sum((nakaRushton(x, params) - y).^2, 'omitnan');
        case {'MLE'}
            % Negative Log-Likelihood function
            % NLL=-y*log(R(c)) - (N-y)*log(1?R(c))
            objectiveFunction = @(params) -sum(y_counts .* log(nakaRushton(x, params)/100) + ...
                               (N - y_counts) .* log(1 - nakaRushton(x, params)/100));
    end

    % Solve the optimization problem
    [bestParams, nLL] = fmincon(objectiveFunction, initialParams, [], [], [], [], lb, ub, [], options);
       
    % Example usage within your main fitting loop or after it
    n = sum(N);  % Total number of data points
    k = length(initialParams);  % Number of parameters in the model
    [aic, aicc, bic] = calculateAIC(nLL, k, n);

    % Prepare the output with fixed beta included
    % Output params = [beta, n, C50]
    params = [beta, bestParams(1), bestParams(2), 0, aicc];
end

function [params_con, params_incon] = fitNakaRushtonOptos(xCon, yCon, xIncon, yIncon, objFunc, options, hypothesisType, fixedParams)
    % N trials
    N_con = [(40*ones(1,sum(xCon==0))),...
        (20*ones(1,sum(xCon>0)))]; % trials per contrast level
    N_incon = N_con; 

    % Convert percent correct into success counts
    yCon_counts = (yCon / 100) .* N_con;
    yIncon_counts = (yIncon / 100) .* N_incon;

    % Estimate initial parameters for both conditions (assuming these functions exst and are correct)
    [initBetaCon, initExpCon, initC50Con, ~] = guessInitialParams(xCon, yCon);
    [~, initExpIncon, initC50Incon, ~] = guessInitialParams(xIncon, yIncon);
    
    % Select hypothesis
    switch hypothesisType
        case 'base'

            % Fixed params
            if sum(~isnan(fixedParams))==0
                % Define the Naka-Rushton function with fitted params = [n, C50, Rmax, beta] (beta is shared across both opto conditions)
                nakaRushtonCon = @(x, params) ((100-params(1)) .* x.^params(2)) ./ (params(3).^params(2) + x.^params(2)) + params(1); %params([1:3])
                nakaRushtonIncon = @(x, params) ((100-(100-params(1))) .* x.^params(4)) ./ (params(5).^params(4) + x.^params(4)) + (100-params(1)); % params([1 4 5]), Using 1-beta_con for beta_incon
            else 
                initExpCon=fixedParams(2,2);
                initExpIncon=fixedParams(3,2);
                % Define the Naka-Rushton function with fitted params = [n, C50, Rmax, beta] (beta is shared across both opto conditions)
                nakaRushtonCon = @(x, params) ((100-params(1)) .* x.^initExpCon) ./ (params(3).^initExpCon + x.^initExpCon) + params(1); %params([1:3])
                nakaRushtonIncon = @(x, params) ((100-(100-params(1))) .* x.^initExpIncon) ./ (params(5).^initExpIncon + x.^initExpIncon) + (100-params(1)); % params([1 4 5]), Using 1-beta_con for beta_incon
            end
            
            % Initial parameter guesses now include only one beta value, for con condition
            % The incon condition's beta is calculated as 1 - beta_con
            initialParams = [initBetaCon, initExpCon, initC50Con, initExpIncon, initC50Incon];

            % Define bounds for the parameters
            lb = [0,      2,     5,      2,       5]; % Lower bounds
            ub = [100, 6, 100,      6,   100]; % Upper bounds

        case 'goris' 
            % Define the Naka-Rushton function with fitted params = [beta, nCon, c50Con, nIncon, c50Incon, xDelta] (beta is shared across both opto conditions)
            nakaRushtonCon = @(x, params) ((100-params(1)) .* (x - params(6)).^params(2)) ./ (params(3).^params(2) + (x - params(6)).^params(2)) + params(1); %params([1:3])
            nakaRushtonIncon = @(x, params) ((100-(100-params(1))) .* (x + params(6)).^params(4)) ./ (params(5).^params(4) + (x + params(6)).^params(4)) + (100-params(1)); % params([1 4 5]), Using 1-beta_con for beta_incon
            
            % Initial parameter guesses now include only one beta value, for con condition
            % The incon condition's beta is calculated as 1 - beta_con
            initxDelta=0;
            initialParams = [initBetaCon, initExpCon, initC50Con, initExpIncon, initC50Incon, initxDelta];

            % Define bounds for the parameters
            lb = [0,      2,     5,      2,       5,     0]; % Lower bounds
            ub = [100, 6, 100,      6,   100,    30]; % Upper bounds
    end

    % Objective function
    switch objFunc
        case {'SSE'}
            objectiveFunction = @(params) sum(((nakaRushtonCon(xCon, params) - yCon).^2), 'omitnan') + ...
                                          sum(((nakaRushtonIncon(xIncon, params) - yIncon).^2), 'omitnan');
        case {'MLE'}
             % Negative Log-Likelihood function combining both conditions
             % NLL=?y*log(R(c)) ? (N?y)*log(1?R(c))
            objectiveFunction = @(params) ...
                -sum(yCon_counts .* log(nakaRushtonCon(xCon, params)/100) + ...
                (N_con - yCon_counts) .* log(1 - nakaRushtonCon(xCon, params)/100)) ...
                -sum(yIncon_counts .* log(nakaRushtonIncon(xIncon, params)/100) + ...
                (N_incon - yIncon_counts) .* log(1 - nakaRushtonIncon(xIncon, params)/100));
    end

    % Solve the optimization problem
    [bestParams,nLL] = fmincon(objectiveFunction, initialParams, [], [], [], [], lb, ub, [], options);
    
    % Example usage within your main fitting loop or after it
    n = sum(N_con);  % Total number of data points
    k = length(initialParams);  % Number of parameters in the model
    [aic, aicc, bic] = calculateAIC(nLL, k, n);
    
    switch hypothesisType
        case 'base'
            % Extract fitted parameters for each condition, ensuring beta_incon = 1 - beta_con
            params_con = [bestParams(1), bestParams(2), bestParams(3), 0, aicc];
            params_incon = [100 - bestParams(1), bestParams(4), bestParams(5), 0, aicc];
        case 'goris'
            % Extract fitted parameters for each condition, ensuring beta_incon = 1 - beta_con
            params_con = [bestParams(1), bestParams(2), bestParams(3), bestParams(6), aicc];
            params_incon = [100 - bestParams(1), bestParams(4), bestParams(5), bestParams(6), aicc];
    end
end

function [initBeta, initExp, initC50, initRmax] = guessInitialParams(x, y)
    initBeta = min(y);
    initRmax = max(y) - min(y);
    initExp = mean([2:4]); % Assuming an intermediate exponent value
    initC50 = computeInitC50(y, initRmax, x);
end

function initC50 = computeInitC50(y, initRmax, x)
    absoluteDelta = abs(y - (initRmax / 2 + min(y)));
    [~, closest_index] = min(absoluteDelta);
    closest_index = unique([max(1, closest_index-1), closest_index, min(length(x), closest_index+1)]);
    initC50 = mean(x(closest_index));
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
