function upFontSize(varargin)
if length(varargin)==2
    size=varargin{1};
    tickLength=varargin{2};
else
    size=24;
    tickLength=0.01;
end
    %INCREASE OVERALL FONT SIZE AND ADJUST FONT WEIGHT
    fontSizeSmall=12;
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', size) % Increase font size globally
    set(findall(gcf, '-property', 'FontWeight'), 'FontWeight', 'Normal') % Set font weight to normal globally

    set(gca, 'LineWidth', 2) % Increase axis line width
    set(gca, 'TickDir', 'out') % Set tick direction outwards
    set(gca, 'TickLength', [tickLength, tickLength]) % Adjust tick length
    set(gcf, 'Color', 'w') % Set figure background color to white
    box off % Turn off the box surrounding the plot

    % After making global changes, find the legend and reset its properties
    leg = findobj(gcf, 'Type', 'Legend');
    if ~isempty(leg)
        % Set the legend's font size back to default or a preferred size
        % MATLAB doesn't have a 'default' keyword, so you might need to specify a number
        set(leg, 'FontSize', fontSizeSmall); % Example: set back to MATLAB's default size or choose your preferred size
        
        % For legend title, if you have one and want to adjust its size separately
        % MATLAB versions R2017a and newer support legend titles directly
        % For older versions, this might not apply
        if isprop(leg, 'Title')
            titleProp = get(leg(end), 'Title');
            set(titleProp, 'FontSize', fontSizeSmall); % Adjust legend title font size if needed
        end
    end
    
    % same for figure text
    hText=findobj(gcf,'Type','Text');
    if ~isempty(hText)
        % Set the legend's font size back to default or a preferred size
        % MATLAB doesn't have a 'default' keyword, so you might need to specify a number
        set(hText, 'FontSize', fontSizeSmall); % Example: set back to MATLAB's default size or choose your preferred size
    end
    
    % Optionally, if you want to remove the upper and right axis lines (spines) explicitly
    ax = gca; % Get current axis handle
    ax.Box = 'off'; % Turn off the box to remove top and right lines
    ax.XAxis.TickDirection = 'out'; % Ensure X-axis ticks are outward
    ax.YAxis.TickDirection = 'out'; % Ensure Y-axis ticks are outward
end
