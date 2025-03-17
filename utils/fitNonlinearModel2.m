function fitNonlinearModel2(clusterPsy, predictors, predictorNames, method, includeInteractions, includeQuadratic, flagAnalyticalTools, applyLog, applyExp, applyReciprocal, applyNakaRushton)
    % Fits a nonlinear model using Lasso or Elastic Net with interaction terms and transformations
    %
    % Args:
    %   clusterPsy: Dependent variable (vector, e.g., behavior performance)
    %   predictors: Matrix of independent variables (each column is a predictor)
    %   predictorNames: Cell array of predictor names (must match the number of predictors)
    %   method: String specifying the regression method ('lasso', 'elasticnet')
    %   includeInteractions: Boolean, whether to include interaction terms
    %   includeQuadratic: Boolean, whether to include quadratic terms
    %   flagAnalyticalTools: Boolean, whether to show diagnostic plots
    %   applyLog: Boolean, whether to apply logarithmic transformation
    %   applyExp: Boolean, whether to apply exponential transformation
    %   applyReciprocal: Boolean, whether to apply reciprocal transformation
    %   applyNakaRushton: Boolean, whether to apply Naka-Rushton transformation
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

    % Initialize expanded predictors and names
    expandedPredictors = predictors;
    expandedPredictorNames = predictorNames;

    % Apply logarithmic transformation if flag is true
    if applyLog
        logTerms = log(abs(predictors) + 1); % Avoid log(0)
        logNames = strcat('log_', predictorNames);
        expandedPredictors = [expandedPredictors, logTerms];
        expandedPredictorNames = [expandedPredictorNames, logNames];
    end

    % Apply exponential transformation if flag is true
    if applyExp
        expTerms = exp(predictors);
        expNames = strcat('exp_', predictorNames);
        expandedPredictors = [expandedPredictors, expTerms];
        expandedPredictorNames = [expandedPredictorNames, expNames];
    end

    % Apply reciprocal transformation if flag is true
    if applyReciprocal
        reciprocalTerms = 1 ./ (predictors + 1e-6); % Avoid division by zero
        reciprocalNames = strcat('recip_', predictorNames);
        expandedPredictors = [expandedPredictors, reciprocalTerms];
        expandedPredictorNames = [expandedPredictorNames, reciprocalNames];
    end

    % Apply Naka-Rushton transformation if flag is true
    if applyNakaRushton
        % Define Naka-Rushton parameters (adjust as necessary)
        Rmax = 1; % Maximum response (scaling factor)
        n = 2;    % Exponent controlling steepness
        k = 0.3;  % Semi-saturation constant
        
        nakaRushtonTerms = (Rmax .* (predictors .^ n)) ./ ((predictors .^ n) + k^n);
        nakaRushtonNames = strcat('naka_', predictorNames);
        expandedPredictors = [expandedPredictors, nakaRushtonTerms];
        expandedPredictorNames = [expandedPredictorNames, nakaRushtonNames];
    end

    % Include interaction terms if flag is true
    if includeInteractions
        numPredictors = size(predictors, 2);
        interactionTerms = [];
        interactionNames = {};
        for i = 1:numPredictors
            for j = i+1:numPredictors % Only include interactions, no self-products
                interactionTerms = [interactionTerms, predictors(:, i) .* predictors(:, j)];
                interactionNames = [interactionNames, {sprintf('%s*%s', predictorNames{i}, predictorNames{j})}];
            end
        end
        expandedPredictors = [expandedPredictors, interactionTerms];
        expandedPredictorNames = [expandedPredictorNames, interactionNames];
    end

    % Include quadratic terms if flag is true
    if includeQuadratic
        quadraticTerms = predictors.^2;
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
    [B, FitInfo] = lasso(expandedPredictors, clusterPsy, 'Standardize', 1, 'MCReps', 1, 'CV', 10, 'Alpha', alpha);

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
    
    % Select non-zero coefficients
    selectedTerms = expandedPredictorNames(bestCoefficients ~= 0);
    selectedCoefficients = bestCoefficients(bestCoefficients ~= 0);

    % Construct the final equation as a single line
    equation = sprintf('Behavior = %.2f', intercept);
    for i = 1:numel(selectedTerms)
        equation = sprintf('%s + %.1f*%s', equation, selectedCoefficients(i), selectedTerms{i});
    end

    % Display results
    fprintf('\n%s Nonlinear Regression Summary:\n', method);
    fprintf('Intercept: %.2f\n', intercept);
    fprintf('\nSelected Terms and Coefficients:\n');
    for i = 1:numel(selectedTerms)
        fprintf('%s: %.2f\n', selectedTerms{i}, selectedCoefficients(i));
    end

    fprintf('\nFinal Equation:\n%s\n', equation);

    fprintf('\nGoodness-of-Fit Metrics:\n');
    fprintf('MAE: %.4f\n', MAE);
    fprintf('RMSE: %.4f\n', RMSE);
    fprintf('R^2: %.4f\n', R2);

    % Predicted vs Actual Plot
    figure('Name', 'Predicted vs. Actual Behavior');
    scatter(clusterPsy, predicted, 100, 'ko', 'LineWidth', 2.5, 'markerFaceColor','w');
    hold on;
    plot([min(clusterPsy), max(clusterPsy)], [min(clusterPsy), max(clusterPsy)], 'k--', 'LineWidth', 2.5);
    xlabel('Actual Behavior');
    ylabel('Predicted Behavior');
    title('Predicted vs. Actual Behavior (Evaluate if points align along the diagonal)');
    xlim([0 100]);
    ylim([0 100]);
    text(50, 10, equation, 'FontSize', 8, 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Interpreter', 'none');
    text(50, 90, sprintf('R^2 = %.2f', R2), 'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');upFontSize(18,.005)

    switch flagAnalyticalTools
        case {'true'}
            % Plot cross-validation results
            lassoPlot(B, FitInfo, 'PlotType', 'CV');
            title(sprintf('%s Nonlinear Regression Cross-Validation', method));
            xlabel('\lambda (log scale)');
            ylabel('Cross-Validation MSE (Look for minimum error, check stability of \\lambda)');
            set(gcf, 'Name', 'Cross-Validation MSE');upFontSize(18,.01)

            % Residual Plot
            residuals = clusterPsy - predicted;
            figure('Name', 'Residuals vs Predicted Behavior');
            scatter(predicted, residuals, 'filled');
            xlabel('Predicted Behavior');
            ylabel('Residuals');
            title('Residuals vs. Predicted (Check for random distribution)');upFontSize(18,.01)

            % Regularization Path
            lassoPlot(B, FitInfo, 'PlotType', 'Lambda', 'XScale', 'log');
            title('Regularization Path (Look for stability of coefficients as \\lambda changes)');
            set(gcf, 'Name', 'Regularization Path');upFontSize(18,.01)
    end
end
