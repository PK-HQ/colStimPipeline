function shadePlot(xBounds, yBounds, Color, Alpha)
% Define the boundaries of the rectangle
xMin = xBounds(1); % Set your xMin
xMax = xBounds(2); % Set your xMax
yMin = yBounds(1); % Set your yMin
yMax = yBounds(2); % Set your yMax

% Coordinates of the rectangle's corners
xCoords = [xMin, xMax, xMax, xMin];
yCoords = [yMin, yMin, yMax, yMax];

% Plot your existing data (for example purposes)
ax=gca;
hold on; % Keep the plot so the rectangle is added
% Fill the rectangle with shading and set transparency
fill(xCoords, yCoords, Color, 'FaceAlpha', Alpha, 'EdgeColor','none','HandleVisibility','off');
end
