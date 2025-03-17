function desiredPos=setFig(varargin)
posInput=~isempty(varargin);
if posInput
    % Get figure size and position
    desiredPos=varargin{1};
    % Create a new figure with the same size
    newFigure = figure; % Create a new figure (popped out)
    set(newFigure, 'Position', desiredPos); % Set the position to match the previous figure
else
    % Get current figure's position and size
    h = figure; % Create a figure (popped out)
    desiredPos = get(h, 'Position'); % Get the position and size of the figure window
end
end
