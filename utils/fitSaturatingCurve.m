function fitSaturatingCurve(x, y, linecolor, hAx, subaxes)
    % Define the logistic function
    logisticFunc = @(b,x) b(1) ./ (1 + exp(-b(2)*(x-b(3))));
    
    % Initial guesses for parameters [L, k, x0]
    initialGuess = [max(y), 1, median(x)];
    
    % Fit the model
    fittedModel = fitnlm(x, y, logisticFunc, initialGuess);
    
    % Get the fitted coefficients
    coefficients = fittedModel.Coefficients{:, 'Estimate'};
    
    % Generate x values for plotting the fitted curve
    xFit = linspace(min(x), max(x), 200);
    
    % Evaluate the fitted model
    yFit = logisticFunc(coefficients, xFit);
    
    % Plot the original data
    axes(hAx(subaxes));
    hold on;
    % Overlay the fitted curve
    plot(xFit, yFit, '--', 'Color', linecolor, 'LineWidth', 2); % Fitted curve in red    
    hold off;
end
