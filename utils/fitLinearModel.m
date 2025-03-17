function fitLinearModel(clusterPsy, predictors, predictorNames, method)
    % Fits a linear model using ordinary regression, Lasso, or Elastic Net
    %
    % Args:
    %   clusterPsy: Dependent variable (vector, e.g., behavior performance)
    %   predictors: Matrix of independent variables (each column is a predictor)
    %   predictorNames: Cell array of predictor names (must match the number of predictors)
    %   method: String specifying the regression method ('linear', 'lasso', 'elasticnet')
    %
    % Returns:
    %   Prints the model summary, coefficients, p-values, and goodness-of-fit metrics.

    % Validate inputs
    if size(predictors, 1) ~= numel(clusterPsy)
        error('The number of rows in predictors must match the length of clusterPsy.');
    end
    if size(predictors, 2) ~= numel(predictorNames)
        error('The number of predictor names must match the number of columns in predictors.');
    end
    if ~ismember(method, {'linear', 'lasso', 'elasticnet'})
        error('Method must be one of: ''linear'', ''lasso'', ''elasticnet''.');
    end

    % Remove NaN values
    validIdx = ~isnan(clusterPsy) & all(~isnan(predictors), 2);
    clusterPsy = clusterPsy(validIdx);
    predictors = predictors(validIdx, :);

    switch method
        case 'linear'
            % Prepare data table for linear model
            dataTable = array2table([predictors, clusterPsy], ...
                'VariableNames', [predictorNames, {'Behavior'}]);

            % Dynamically construct the formula
            predictorFormula = strjoin(predictorNames, ' + ');
            formula = sprintf('Behavior ~ %s', predictorFormula);

            % Fit the linear model
            mdl = fitlm(dataTable, formula);

            % Display model summary
            disp(mdl);

            % Extract coefficients and p-values
            coefficients = mdl.Coefficients.Estimate;
            pValues = mdl.Coefficients.pValue;

            % Predict values
            predicted = predict(mdl, dataTable);

            % Calculate goodness-of-fit metrics
            MAE = mean(abs(predicted - clusterPsy));
            RMSE = sqrt(mean((predicted - clusterPsy).^2));
            R2 = mdl.Rsquared.Ordinary;

            % Print coefficients, p-values, and metrics
            fprintf('\nLinear Regression Coefficients and p-values:\n');
            equation = sprintf('Behavior = %.2f', coefficients(1));
            for i = 2:numel(coefficients)
                fprintf('%s: %.2f (p=%.2g)\n', ...
                    mdl.CoefficientNames{i}, coefficients(i), pValues(i));
                equation = sprintf('%s + %.2f*%s', equation, coefficients(i), predictorNames{i-1});
            end

            fprintf('\nFinal Equation:\n%s\n', equation);
            fprintf('\nGoodness-of-Fit Metrics:\n');
            fprintf('MAE: %.4f\n', MAE);
            fprintf('RMSE: %.4f\n', RMSE);
            fprintf('R^2: %.4f\n', R2);

        case {'lasso', 'elasticnet'}
            % Lasso or Elastic Net parameters
            alpha = 1; % Lasso by default
            if strcmp(method, 'elasticnet')
                alpha = 0.5; % Elastic Net mixing parameter
            end

            % Perform Lasso or Elastic Net regression
            [B, FitInfo] = lasso(predictors, clusterPsy, 'CV', 10, 'Alpha', alpha);

            % Optimal model coefficients
            idxLambdaMinMSE = FitInfo.IndexMinMSE;
            bestCoefficients = B(:, idxLambdaMinMSE);
            intercept = FitInfo.Intercept(idxLambdaMinMSE);

            % Predict values
            predicted = intercept + predictors * bestCoefficients;

            % Calculate goodness-of-fit metrics
            MAE = mean(abs(predicted - clusterPsy));
            RMSE = sqrt(mean((predicted - clusterPsy).^2));
            R2 = 1 - sum((predicted - clusterPsy).^2) / sum((clusterPsy - mean(clusterPsy)).^2);

            % Bootstrapping for p-value estimation
            nBoot = 1000; % Number of bootstrap iterations
            bootCoefficients = zeros(nBoot, size(predictors, 2));
            for i = 1:nBoot
                resampleIdx = randi(size(predictors, 1), size(predictors, 1), 1);
                resampledPredictors = predictors(resampleIdx, :);
                resampledClusterPsy = clusterPsy(resampleIdx);
                [BResample, ~] = lasso(resampledPredictors, resampledClusterPsy, ...
                    'Lambda', FitInfo.Lambda(idxLambdaMinMSE), 'Alpha', alpha);
                bootCoefficients(i, :) = BResample';
            end

            % Compute p-values
            pValues = mean(abs(bootCoefficients) >= abs(bestCoefficients'), 1);

            % Print coefficients, p-values, and metrics
            fprintf('\n%s Regression Coefficients and p-values:\n', method);
            equation = sprintf('Behavior = %.2f', intercept);
            for i = 1:numel(bestCoefficients)
                fprintf('%s: %.2f (p=%.2g)\n', ...
                    predictorNames{i}, bestCoefficients(i), pValues(i));
                if bestCoefficients(i) ~= 0
                    equation = sprintf('%s + %.2f*%s', equation, bestCoefficients(i), predictorNames{i});
                end
            end

            fprintf('\nFinal Equation:\n%s\n', equation);
            fprintf('\nGoodness-of-Fit Metrics:\n');
            fprintf('MAE: %.4f\n', MAE);
            fprintf('RMSE: %.4f\n', RMSE);
            fprintf('R^2: %.4f\n', R2);

            % Plot cross-validation results
            lassoPlot(B, FitInfo, 'PlotType', 'CV');
            title(sprintf('%s Regression Cross-Validation', method));
    end
end
