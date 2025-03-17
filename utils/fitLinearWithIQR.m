function fitLinearWithIQR(x, y, linecolor)
    % Function to fit a standard linear line to data after removing outliers using IQR
    %
    % Args:
    %   x: Independent variable data (vector)
    %   y: Dependent variable data (vector)
    %   linecolor: Color of the fitted line (e.g., [0, 0, 1] for blue)

    % Remove NaN values
    validIdx = ~isnan(x) & ~isnan(y);
    x = x(validIdx);
    y = y(validIdx);

    % Check if there are enough points to perform outlier removal
    if numel(x) < 3
        warning('Not enough points to perform outlier removal (requires at least 3 points).');
        return;
    end

    % Identify outliers using IQR
    q1 = quantile(y, 0.25); % 1st quartile
    q3 = quantile(y, 0.75); % 3rd quartile
    iqr = q3 - q1; % Interquartile range
    lowerBound = q1 - 1.5 * iqr; % Lower bound
    upperBound = q3 + 1.5 * iqr; % Upper bound

    % Filter inliers
    inlierIdx = (y >= lowerBound) & (y <= upperBound);
    xInliers = x(inlierIdx);
    yInliers = y(inlierIdx);

    % Check if enough inliers remain after outlier removal
    if numel(xInliers) < 3
        warning('Not enough inliers to perform standard fitting (requires at least 3 points).');
        return;
    end

    % Perform standard linear fitting on inliers
    p = polyfit(xInliers, yInliers, 1);

    % Generate x values for the fitted line
    xFit = linspace(min(xInliers), max(xInliers), 200);

    % Compute the fitted y values
    yFit = polyval(p, xFit);

    % Calculate R^2 for the standard linear fit
    yPred = polyval(p, xInliers); % Predicted values for inliers
    ssRes = sum((yInliers - yPred).^2); % Residual sum of squares
    ssTot = sum((yInliers - mean(yInliers)).^2); % Total sum of squares
    R2 = 1 - (ssRes / ssTot); % Coefficient of determination

    % Plot the original data
    hold on;
    scatter(x, y, 50, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.8, 0.8, 0.8], 'HandleVisibility', 'off'); % Original data (gray)

    % Highlight inliers
    scatter(xInliers, yInliers, 50, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', linecolor, 'HandleVisibility', 'off');

    % Plot the fitted line
    plot(xFit, yFit, '-', 'Color', linecolor, 'LineWidth', 2);

    % Annotate with R^2
    text(mean(xInliers), max(yInliers) * 0.9, sprintf('R^2 = %.2f', R2), 'Color', linecolor, 'FontSize', 12);

end
