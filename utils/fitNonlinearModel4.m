function fitNonlinearModel4(clusterPsy, predictors, predictorNames, method, includeInteractions, transformType, predictorSelection, flagDiagnostics, optimizationGoal)
    % Fits a nonlinear model using Lasso or Elastic Net with transformations
    %
    % Args:
    %   clusterPsy: Dependent variable (vector, e.g., behavior performance)
    %   predictors: Matrix of independent variables (each column is a predictor)
    %   predictorNames: Cell array of predictor names (must match the number of predictors)
    %   method: String specifying the regression method ('lasso', 'elasticnet')
    %   includeInteractions: Boolean, whether to include interaction terms
    %   transformType: String, {'nakarushton', 'weibull'}, transformation type
    %   predictorSelection: String, {'original', 'original+transformed', 'transformed'}
    %   flagAnalyticalTools: Boolean, whether to show diagnostic plots
    %   optimizationGoal: String, {'R2', 'RMSE'} to optimize the transformations
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
    if ~ismember(predictorSelection, {'original', 'original+transformed', 'transformed'})
        error('Predictor selection must be ''original'', ''original+transformed'', or ''transformed''.');
    end
    if ~ismember(optimizationGoal, {'R2', 'RMSE'})
        error('Optimization goal must be ''R2'' or ''RMSE''.');
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

    % Optimization settings
    opts = optimset('fminsearch');
    opts.Display = 'off';
    opts.MaxIter = 100000;
    opts.MaxFunEvals = 200000;
    opts.TolFun = 1e-8;
    opts.TolX = 1e-8;

    % Apply transformations with optimization
    switch transformType
        case 'nakarushton'
            % Define bounds and initial parameters
            initialParams = [1, 2, 0.3]; % [Rmax, n, k]
            lb = [1e-8, 1e-8, 1e-8];              % Lower bounds
            ub = [10, 10, 1];            % Upper bounds
            % Optimize
            nakaParams = fminsearchbnd(@(params) optimizeTransform(params, clusterPsy, predictors, 'nakarushton', optimizationGoal), ...
                                       initialParams, lb, ub, opts);
            [~, nakaRushtonTerms] = optimizeTransform(nakaParams, clusterPsy, predictors, 'nakarushton', optimizationGoal);
            nakaRushtonNames = strcat('\it{f}\rm(', predictorNames,')');
            expandedPredictors = [expandedPredictors, nakaRushtonTerms];
            expandedPredictorNames = [expandedPredictorNames, nakaRushtonNames];

        case 'weibull'
            % Define bounds and initial parameters
            initialParams = [1, 0.3, 2, 2.5]; % [A, B, alpha, beta]
            lb = [1e-8, 1e-8, 1e-8, 1e-8];              % Lower bounds
            ub = [10, 1, 10, 20];               % Upper bounds            
            % Optimize
            weibullParams = fminsearchbnd(@(params) optimizeTransform(params, clusterPsy, predictors, 'weibull', optimizationGoal), ...
                                          initialParams, lb, ub, opts);
            [~, weibullTerms] = optimizeTransform(weibullParams, clusterPsy, predictors, 'weibull', optimizationGoal);
            weibullNames = strcat('\it{f}\rm(', predictorNames, ')');
            expandedPredictors = [expandedPredictors, weibullTerms];
            expandedPredictorNames = [expandedPredictorNames, weibullNames];
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

    % Select predictors based on the predictorSelection flag
    switch predictorSelection
        case 'original'
            finalPredictors = predictors;
            finalPredictorNames = predictorNames;
        case 'original+transformed'
            finalPredictors = expandedPredictors;
            finalPredictorNames = expandedPredictorNames;
        case 'transformed'
            finalPredictors = expandedPredictors(:, size(predictors, 2)+1:end);
            finalPredictorNames = expandedPredictorNames(size(predictors, 2)+1:end);
    end

    % Lasso or Elastic Net parameters
    alpha = 1; % Lasso by default
    if strcmp(method, 'elasticnet')
        alpha = 0.5; % Elastic Net mixing parameter
    end

    %% Perform Lasso or Elastic Net regression
    [B, FitInfo] = lasso(finalPredictors, clusterPsy, 'Standardize', 1, 'MCReps', 1000, 'CV', 10, 'Alpha', alpha);

    % Optimal model coefficients
    idxLambdaMinMSE = FitInfo.IndexMinMSE;
    bestCoefficients = B(:, idxLambdaMinMSE);
    intercept = FitInfo.Intercept(idxLambdaMinMSE);

    % Predict values
    predicted = intercept + finalPredictors * bestCoefficients;

    % Calculate goodness-of-fit metrics
    MAE = mean(abs(predicted - clusterPsy));
    RMSE = sqrt(mean((predicted - clusterPsy).^2));
    R2 = 1 - sum((predicted - clusterPsy).^2) / sum((clusterPsy - mean(clusterPsy)).^2);

    % Select non-zero coefficients
    selectedTerms = finalPredictorNames(bestCoefficients ~= 0);
    selectedCoefficients = bestCoefficients(bestCoefficients ~= 0);

    %% Construct the final equation based on predictorSelection
    switch predictorSelection
        case {'original', 'transformed'}
            % Single-line equation for 'original' or 'transformed'
            equation = sprintf('Behavior = %.2f', intercept);
            for i = 1:numel(selectedTerms)
                equation = sprintf('%s + %.1f*%s', equation, selectedCoefficients(i), selectedTerms{i});
            end
        case 'original+transformed'
            % Two-line equation for 'original+transformed'
            originalTerms = selectedTerms(1:size(predictors, 2));
            originalCoefficients = selectedCoefficients(1:size(predictors, 2));
            transformedTerms = selectedTerms(size(predictors, 2)+1:end);
            transformedCoefficients = selectedCoefficients(size(predictors, 2)+1:end);

            % First line: Original predictors
            equation = sprintf('Predicted = %.1f', intercept);
            for i = 1:numel(originalTerms)
                equation = sprintf('%s + %.1f*%s', equation, originalCoefficients(i), originalTerms{i});
            end

            % Second line: Transformed predictors
            transformedEquation = '';
            for i = 1:numel(transformedTerms)
                transformedEquation = sprintf('%s + %.1f*%s', transformedEquation, transformedCoefficients(i), transformedTerms{i});
            end
    end

    % Display results
    fprintf('\n%s Nonlinear Regression Summary:\n', method);
    fprintf('Intercept: %.2f\n', intercept);
    fprintf('\nSelected Terms and Coefficients:\n');
    for i = 1:numel(selectedTerms)
        fprintf('%s: %.2f\n', selectedTerms{i}, selectedCoefficients(i));
    end
    fprintf('\nFinal Equation:\n%s\n%s\n', equation, transformedEquation);

    fprintf('\nGoodness-of-Fit Metrics:\n');
    fprintf('MAE: %.4f\n', MAE);
    fprintf('RMSE: %.4f\n', RMSE);
    fprintf('R^2: %.4f\n', R2);

    %% Predicted vs Actual Plot
    figure('Name', 'Predicted vs. Actual Behavior');
    plot([0 100], [0 100], '--','Color',[.65 .65 .65], 'LineWidth', 2); hold on;

    %plot([min(clusterPsy), max(clusterPsy)], [min(clusterPsy), max(clusterPsy)], '--','Color',[.65 .65 .65], 'LineWidth', 2.5); hold on;
    scatter(clusterPsy, predicted, 100, 'ko', 'LineWidth', 2.5, 'markerFaceColor','w'); hold on;
    xlabel('True Performance (%)');
    ylabel('Predicted Performance (%)');
    xlim([0 100]); addSkippedTicks(0,100,12.5,'x')
    ylim([0 100]); addSkippedTicks(0,100,12.5,'y')

    % Add R² value to the plot
    text(5, 22, sprintf('R^2 = %.2f', R2), 'FontSize', 12, 'Color', 'k', 'HorizontalAlignment', 'left', 'FontWeight', 'bold');
    upFontSize(24, .01); axis square
    legend({'Reference line','True vs. Predicted'},'Location','northwest','FontSize', 18)

    fontSizeEqn=12;
    % Add Equation
    if strcmp(predictorSelection, 'original+transformed')
        text(5, 14, equation, 'FontSize', fontSizeEqn, 'Color', 'k', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Interpreter', 'none');
        text(15, 9, transformedEquation, 'FontSize', fontSizeEqn, 'Color', 'k', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Interpreter', 'Tex');
    else
        text(5, 7, equation, 'FontSize', fontSizeEqn, 'Color', 'k', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Interpreter', 'Tex');
    end    
    % Add function parameters to the plot if applicable
    if strcmp(transformType, 'weibull')
        text(5, 4.5, sprintf('\\it{f}\\rm_{Weibull}: A=%.1f, B=%.1f, \\alpha=%.1f, \\beta=%.1f', weibullParams(1), weibullParams(2), weibullParams(3), weibullParams(4)), ...
            'FontSize', fontSizeEqn, 'Color', 'k', 'HorizontalAlignment', 'left', 'Interpreter', 'tex');
    elseif strcmp(transformType, 'nakarushton')
        text(5, 4.5, sprintf('\\it{f}\\rm_{Naka-Rushton}: R_{max}=%.1f, n=%.1f, \\kappa=%.1f', nakaParams(1), nakaParams(2), nakaParams(3)), ...
            'FontSize', fontSizeEqn, 'Color', 'k', 'HorizontalAlignment', 'left', 'Interpreter', 'tex');
    end

    %% Diagnostic plots
    if flagDiagnostics==true
        % Plot cross-validation results
        lassoPlot(B, FitInfo, 'PlotType', 'CV');
        title(sprintf('%s Nonlinear Regression Cross-Validation', method));
        xlabel('\lambda (log scale)');
        ylabel('Cross-Validation MSE (Look for minimum error, check stability of \\lambda)');

        % Residual Plot
        residuals = clusterPsy - predicted;
        figure('Name', 'Residuals vs Predicted Behavior');
        scatter(predicted, residuals, 'filled');
        xlabel('Predicted Behavior');
        ylabel('Residuals');
        title('Residuals vs. Predicted (Check for random distribution)');

        % Regularization Path
        lassoPlot(B, FitInfo, 'PlotType', 'Lambda', 'XScale', 'log');
        title('Regularization Path (Look for stability of coefficients as \\lambda changes)');
    end
    
end

function [metric, transformed] = optimizeTransform(params, clusterPsy, predictors, transformType, optimizationGoal)
    % Apply Naka-Rushton or Weibull transformation and compute the optimization metric
    if strcmp(transformType, 'nakarushton')
        Rmax = params(1);
        n = params(2);
        k = params(3);
        transformed = (Rmax .* (predictors .^ n)) ./ ((predictors .^ n) + k^n);
    elseif strcmp(transformType, 'weibull')
        A = params(1);
        B = params(2);
        alpha = params(3);
        beta = params(4);
        
        % Ensure predictors are non-negative and clamped to a minimum value
        predictors = max(predictors, 1e-8); % Clamp values to prevent invalid operations
        
        % Apply Weibull transformation
        transformed = A .* (1 - exp(-((predictors ./ B) .^ alpha))) .^ beta;
    end

    % Ensure transformed values remain real
    if any(imag(transformed) ~= 0)
        disp('Complex values detected in transformation. Clamping and using real part only.')
        transformed = real(transformed); % Use only the real part
    end

    % Predict behavior using simple linear fit and compute RMSE or R²
    mdl = fitlm(transformed, clusterPsy);
    if strcmp(optimizationGoal, 'R2')
        metric = -mdl.Rsquared.Ordinary; % Negate to maximize R²
    else
        metric = mdl.RMSE; % Minimize RMSE
    end
end
