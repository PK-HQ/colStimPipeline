function params = fitNakaRushtonV2(xBlocks, yBlocks)
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
    params = zeros(nBlocks, nConditions, 3); % For Beta, exponent, and C50 for each condition and block
    
    % Loop over blocks
    for b = 1:nBlocks
        % Baseline condition fitting with fixed Beta and Rmax
        xi = squeeze(xBlocks(1,:,b)); % Contrast values for baseline
        yi = squeeze(yBlocks(1,:,b)); % Percentage correct for baseline
        params(b,1,:) = fitWithFixedBetaRmax(xi, yi);
        
        % Fit con-opto and incon-opto conditions simultaneously
        xi_con = squeeze(xBlocks(2,:,b)); % Contrast values for con-opto
        yi_con = squeeze(yBlocks(2,:,b)); % Percentage correct for con-opto
        xi_incon = squeeze(xBlocks(3,:,b)); % Contrast values for incon-opto
        yi_incon = squeeze(yBlocks(3,:,b)); % Percentage correct for incon-opto
        [params_con, params_incon] = fitConInconOpto(xi_con, yi_con, xi_incon, yi_incon);
        
        params(b,2,:) = params_con;
        params(b,3,:) = params_incon;
    end
end

%% Subfunctions
function params = fitWithFixedBetaRmax(xi, yi)
    % Define the fixed beta and Rmax
    beta = 50;
    Rmax = 50;
    
    % Define the Naka-Rushton function with fixed beta and Rmax
    % Here, params = [n, C50]
    nakaRushton = @(xi, params) (Rmax .* xi.^params(1)) ./ (params(2).^params(1) + xi.^params(1)) + beta;
    
    % Objective function to minimize: sum of squared differences
    objectiveFunction = @(params) sum((nakaRushton(xi, params) - yi).^2, 'omitnan');
    
    % Initial guesses for the parameters [n, C50]
    initialParams = [2, 10];  % These should be adjusted based on your data
    
    % No constraints on parameters, but you could add bounds if necessary
    lb = [1, 5];  % Lower bounds: no negative values for exponent and C50
    ub = [6, 100];  % Upper bounds: allowing exponent and C50 to be any positive value
    
    % Options for fmincon
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point');
    
    % Solve the optimization problem
    [bestParams, ~] = fmincon(objectiveFunction, initialParams, [], [], [], [], lb, ub, [], options);
    
    % Prepare the output with fixed beta included
    % Output params = [beta, n, C50]
    params = [beta, bestParams(1), bestParams(2)];
end

function [params_con, params_incon] = fitConInconOpto(xi_con, yi_con, xi_incon, yi_incon)
    % Use initial parameter estimation functions as defined previously

    % Estimate initial parameters for both conditions (assuming these functions exist and are correct)
    [initBetaCon, initExpCon, initC50Con, ~] = guessInitialParams(xi_con, yi_con);
    [~, initExpIncon, initC50Incon, ~] = guessInitialParams(xi_incon, yi_incon);
    
    % Initial parameter guesses now include only one beta value, for con condition
    % The incon condition's beta is calculated as 1 - beta_con
    initialParams = [initBetaCon, initExpCon, initC50Con, initExpIncon, initC50Incon];
    
    % Define bounds for the parameters
    lb = [0, 1, 5, -100, 5]; % Lower bounds
    ub = [100, 10, 100, 10, 100]; % Upper bounds

    % Objective function incorporating shared beta
    nakaRushtonCon = @(xi, params) ((100-params(1)) .* xi.^params(2)) ./ (params(3).^params(2) + xi.^params(2)) + params(1);
    nakaRushtonIncon = @(xi, params) ((100-(100-params(1))) .* xi.^params(4)) ./ (params(5).^params(4) + xi.^params(4)) + (100-params(1)); % Using 1-beta_con for beta_incon
    
    objectiveFunction = @(params) sum(((nakaRushtonCon(xi_con, params) - yi_con).^2), 'omitnan') + ...
                                  sum(((nakaRushtonIncon(xi_incon, params) - yi_incon).^2), 'omitnan');

    % Options for fmincon
    options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off', ...
        'MaxIterations', 20000, 'TolFun', 1e-10, 'TolX', 1e-10);
    % Solve the optimization problem
    bestParams = fmincon(objectiveFunction, initialParams, [], [], [], [], lb, ub, [], options);
    
    % Extract fitted parameters for each condition, ensuring beta_incon = 1 - beta_con
    params_con = [bestParams(1), bestParams(2), bestParams(3)];
    params_incon = [100 - bestParams(1), bestParams(4), bestParams(5)];
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
