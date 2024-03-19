function plotFittedLine(fittedX, fittedY, lineColor)
    % Plot fitted line on the given axes with specified styles and colors
    % ax: Axes handle where the plot will be drawn
    % fittedX: X values of the fitted line
    % fittedY: Y values of the fitted line
    % lineStyle: Style of the line (e.g., '-', '--', ':')
    % lineColor: Color of the line (e.g., 'r', 'g', [0.5, 0.5, 0.5])
    % markerStyle: Marker style (e.g., 'o', '*', 'none')
    % markerColor: Marker color (e.g., 'b', 'k', [1, 0, 0])

    % Ensure plotting in the given axes
    lineStyle='-';
    if strcmp(lineColor,'k') || isequal(lineColor,[0 0 0])
        markerFaceColor='w';
    else
        markerFaceColor=lineColor;
    end
    % Plot the line with specified styles and colors
    hLine=plot(fittedX, fittedY, 'LineStyle', lineStyle, 'Color', lineColor, 'LineWidth', 2.5, ...
         'Marker', 'none', 'MarkerEdgeColor', 'k', ...
         'MarkerFaceColor', markerFaceColor);
     
     % flip to bottom
chH = get(gca,'Children');
uistack(hLine,'bottom');
end
