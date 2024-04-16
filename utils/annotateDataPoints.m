function annotateDataPoints(x, y, counts, markerFaceColor)
    % Calculate luminance of the markerFaceColor
    % Assume markerFaceColor is an RGB triplet [R, G, B]
    luminance = 0.299 * markerFaceColor(1) + 0.587 * markerFaceColor(2) + 0.114 * markerFaceColor(3);

    % Choose text color based on luminance
    if luminance < 0.5
        textColor = 'w'; % Use white text for dark backgrounds
    else
        textColor = 'k'; % Use black text for light backgrounds
    end

    % Annotate each data point with its count
    for i = 1:length(x)
        if ~isnan(y(i))
            text(x(i), y(i), sprintf('%d', counts(i)), ...
                'VerticalAlignment', 'middle', ...
                'HorizontalAlignment', 'center', ...
                'FontSize', 8, ...
                'FontWeight', 'bold', ...
                'Color', textColor,...
                'HandleVisibility','on'); % Apply chosen text color
        end
    end
end
