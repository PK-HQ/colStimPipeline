% Function to add vertical label to the left of each row
function addVerticalLabel(hAx, row, col, nCols, label)
    if col==1
        % Get position of the first subplot in the row
        subplotPos = get(hAx((row - 1) * nCols + 1), 'Position');
        
        % Position for the label (to the left of the subplot)
        xPos = subplotPos(1) - 0.035; % Adjust this value as needed
        yPos = subplotPos(2) + subplotPos(4)/2 - length(label)*.00085; % Centered vertically in the subplot
    
        % Add the text to the figure, not the axes
        fig = gcf;
        annotation(fig, 'textbox', [xPos, yPos, 0.01, 0.01], 'String', label, ...
                   'Rotation', 90, 'HorizontalAlignment', 'center', ...
                   'VerticalAlignment', 'middle', 'EdgeColor', 'none');
    end
end