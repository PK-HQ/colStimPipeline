function fitNonlinearModelWithNROptimization(clusterPsy, predictors, predictorNames, method, includeInteractions, includeQuadratic)
    % Fits a nonlinear model using Lasso or Elastic Net with optimized Naka-Rushton terms
    %
    % Args:
    %   clusterPsy: Dependent variable (vector, e.g., behavior performance)
    %   predictors: Matrix of independent variables (each column is a predictor)
    %   predictorNames: Cell array of predictor names (must match the number of predictors)
    %   method: String specifying the regression method ('lasso', 'elasticnet')
    %   includeInteractions: Boolean, whether to include interaction terms
    %   includeQuadratic: Boolean, whether to include quadratic terms
    %
    % Returns:
    %   Prints the model summary, selected terms, coefficients, and goodness-of-fit metrics.
    %   Includes plots for diagnostics and displays the final equation.

    % Validate inputs
    if size(predictors, 1) ~= numel(clusterPsy)
        error('The number of rows in predictors must match the length of clusterPsy.');
    end
    if size(predictors, 2) ~= numel(predictorNames)
        error('The number of predictor names must match the number of columns in predictors.');
    end
    if ~ismember(method, {'lasso', 'elasticnet'})
        error('Method must be one of: ''lasso'', ''elasticnet''.');
    end

    % Remove NaN values
    validIdx = ~isnan(clusterPsy) & all(~isnan(predictors), 2);
    clusterPsy = clusterPsy(validIdx);
    predictors = predictors(validIdx, :);

    % **Apply Z-Score Scaling**
    predictors = (predictors - mean(predictors, 1)) ./ std(predictors, 0, 1);

    % Optimize NR parameters
    fprintf('Optimizing Naka-Rushton parameters...\n');
    [optimizedPredictors, nrParams] = optimizeNakaRushton(predictors);

    fprintf('Optimized NR Parameters: Rmax=%.2f, n=%.2f, k=%.2f\n', nrParams(1), nrParams(2), nrParams(3));

    % Initialize expanded predictors and names
    expandedPredictors = optimizedPredictors;
    expandedPredictorNames = predictorNames;

    % Include interaction terms if flag is true
    if includeInteractions
        numPredictors = size(optimizedPredictors, 2);
        interactionTerms = [];
        interactionNames = {};
        for i = 1:numPredictors
            for j = i+1:numPredictors % Only include interactions, no self-products
                interactionTerms = [interactionTerms, optimizedPredictors(:, i) .* optimizedPredictors(:, j)];
                interactionNames = [interactionNames, {sprintf('%s*%s', predictorNames{i}, predictorNames{j})}];
            end
        end
        expandedPredictors = [expandedPredictors, interactionTerms];
        expandedPredictorNames = [expandedPredictorNames, interactionNames];
    end

    % Include quadratic terms if flag is true
    if includeQuadratic
        quadraticTerms = optimizedPredictors.^2;
        quadraticNames = strcat(predictorNames, '^2');
        expandedPredictors = [expandedPredictors, quadraticTerms];
        expandedPredictorNames = [expandedPredictorNames, quadraticNames];
    end

    % Lasso or Elastic Net parameters
    alpha = 1; % Lasso by default
    if strcmp(method, 'elasticnet')
        alpha = 0.5; % Elastic Net mixing parameter
    end

    % Perform Lasso or Elastic Net regression
    [B, FitInfo] = lasso(expandedPredictors, clusterPsy, 'Standardize', 1, 'MCReps', 100, 'CV', 10, 'Alpha', alpha);

    % Optimal model coefficients
    idxLambdaMinMSE = FitInfo.IndexMinMSE;
    bestCoefficients = B(:, idxLambdaMinMSE);
    intercept = FitInfo.Intercept(idxLambdaMinMSE);

    % Predict values
    predicted = intercept + expandedPredictors * bestCoefficients;

    % Calculate goodness-of-fit metrics
    MAE = mean(abs(predicted - clusterPsy));
    RMSE = sqrt(mean((predicted - clusterPsy).^2));
    R2 = 1 - sum((predicted - clusterPsy).^2) / sum((clusterPsy - mean(clusterPsy)).^2);

    % Display results
    fprintf('\nGoodness-of-Fit Metrics:\n');
    fprintf('MAE: %.4f\n', MAE);
    fprintf('RMSE: %.4f\n', RMSE);
    fprintf('R^2: %.4f\n', R2);

    % Additional diagnostic plots...

end

function [optimizedPredictors, nrParams] = optimizeNakaRushton(predictors)
    % Optimizes the Naka-Rushton parameters for each predictor
    %
    % Args:
    %   predictors: Matrix of independent variables (each column is a predictor)
    %
    % Returns:
    %   optimizedPredictors: Transformed predictors after applying optimized NR function
    %   nrParams: Optimized NR parameters [Rmax, n, k]

    numPredictors = size(predictors, 2);
    optimizedPredictors = zeros(size(predictors));
    nrParams = zeros(3, numPredictors);

    for i = 1:numPredictors
        x = predictors(:, i);
        % Initial parameter guesses
        initialParams = [1, 2, 0.3]; % Rmax, n, k
        lb = [0, 0, 0]; % Lower bounds
        ub = [10, 10, 1]; % Upper bounds

        % Objective function to minimize the sum of squared errors
        objective = @(params) sum((applyNakaRushton(x, params) - x).^2);

        % Optimize parameters using fmincon
        %opts = optimoptions('fmincon', 'Display', 'none');
        params = fmincon(objective, initialParams, [], [], [], [], lb, ub, []);

        % Store optimized parameters and transform predictor
        nrParams(:, i) = params;
        optimizedPredictors(:, i) = applyNakaRushton(x, params);
    end
end

function response = applyNakaRushton(x, params)
    % Applies the Naka-Rushton function to a predictor
    %
    % Args:
    %   x: Predictor variable (vector)
    %   params: Parameters [Rmax, n, k]
    %
    % Returns:
    %   response: Transformed predictor

    Rmax = params(1);
    n = params(2);
    k = params(3);
    response = (Rmax * x.^n) ./ (x.^n + k^n);
end
