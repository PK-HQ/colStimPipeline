function interactiveFittingGUI()
    %% Parameter bounds and initial guess
    lb = [0,   0,   0,   0,   0.1, 0.5,   1,   0.0,   0 ];
    ub = [1,   1,   1,   1,   300,   4,   600, 1.0, 50];
    p0 = [0.45, 0.56, 0.04, 0.59,  172, 1.9, 400, 0.39,  18];
    paramNames = {'a','b','l','w','G0','n','rmx','e','o'};
    
    % Data points (example)
    data.xBlock = [0     8    11    15    40; ...
                   0     7    10    13    40; ...
                   0     7    10    13    40];
    data.yBlock = [55    40    80    95   100; ...
                   65    70    85   100   100; ...
                   40    45    65    95   100];
               
    % Define contrast vector and simulation options
    cmin = 0; cmax = 100;
    x = linspace(cmin, cmax, 100);
    options = struct('nTrials', 100000);
    
    %% Create uifigure and main plot
    % Create uifigure with white background
    fig = uifigure('Name', 'Interactive Fitting (columnarBayesMdl)', ...
        'Position', [100 100 800 600], 'Color', 'w');
    
    % Create uiaxes with font size 16
    ax = uiaxes(fig, 'Position', [50 200 700 350], 'FontSize', 16,'LineWidth',2.5);

    hold(ax, 'on');
    % Colors for data: red, blue, black (matching your compareMdls)
    lineCol = [1 1 1; 0.9759 0.1153 0.1442; 0 0.1176 0.8284];
    
    % Plot data points for each condition
    for cond = 1:3
        scatter(ax, data.xBlock(cond,:), data.yBlock(cond,:), 200, ...
            'LineWidth', 2.5, 'Marker', 'o', 'MarkerFaceColor', lineCol(cond,:), ...
            'MarkerEdgeColor', 'k');
    end
    
    % Compute initial curves using columnarBayesMdl with p0
    curves = columnarBayesMdl(x, p0, options);
    h_baseline = plot(ax, x, curves.pcntrl, 'k', 'LineWidth', 2);
    h_con      = plot(ax, x, curves.pc, 'r', 'LineWidth', 2);
    h_incon    = plot(ax, x, curves.pic, 'b', 'LineWidth', 2);
    
    xlabel(ax, 'Gabor Contrast (%)');
    ylabel(ax, 'Correct (%)');
    title(ax, 'Interactive Fitting (columnarBayesMdl)');
    axis(ax, 'square'); ylim(ax, [0 100]);
    %upFontSize(16,.01)
    %% Create slider panel for bottom row (simulate subplot(4,3,10:12))
    % We'll create a uipanel occupying the bottom ~25% of the figure.
    sliderPanel = uipanel(fig, 'Position', [50 20 700 150], 'BorderColor','none','BorderWidth',0, 'BackgroundColor','w');
    
    % Arrange 9 sliders in a 3x3 grid within sliderPanel.
    nSliders = 9; nCols = 9; nRows = 1;
    panelPos = sliderPanel.Position; % [left bottom width height]
    sliderWidth = 2.75 * panelPos(3) / nCols;
    sliderHeight = 1.5* panelPos(4) / nRows;
    
    sliderHandles = gobjects(nSliders,1);
    labelHandles = gobjects(nSliders,1);
    currentParams = p0;
    
    for i = 1:nSliders
        col = mod(i-1, nCols);
        row = floor((i-1) / nCols);
        % Define a margin within each grid cell for aesthetics
        marginX = 0.3 * sliderWidth;
        marginY = 0.2 * sliderHeight;
        % Position relative to sliderPanel (x, y, width, height)
        pos = [col * sliderWidth + marginX, ...
               (nRows - row - 1)*sliderHeight + marginY, ...
               sliderWidth - 2*marginX, sliderHeight - 2*marginY];
        % Create vertical uislider; set Limits from lb/ub and initial Value from p0
        sliderHandles(i) = uislider(sliderPanel, ...
            'Position', pos, ...
            'Orientation', 'vertical', ...
            'Limits', [lb(i) ub(i)], ...
            'Value', p0(i));
        % Create a uilabel to display the parameter name and value just below the slider
        labelPos = [pos(1)-30, pos(2)+150, pos(3), 15];
        labelHandles(i) = uilabel(sliderPanel, 'Position', labelPos, ...
            'Text', sprintf('%s: %.2f', paramNames{i}, p0(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 14, 'BackgroundColor','w');
        % Set the ValueChangedFcn callback for each slider
        sliderHandles(i).ValueChangedFcn = @(sld,event) sliderCallback();
    end
    
    %% Callback function: update parameter vector and plot
    function sliderCallback()
        for j = 1:nSliders
            currentParams(j) = sliderHandles(j).Value;
            labelHandles(j).Text = sprintf('%s: %.2f', paramNames{j}, currentParams(j));
        end
        newCurves = columnarBayesMdl(x, currentParams, options);
        set(h_baseline, 'YData', newCurves.pcntrl);
        set(h_con, 'YData', newCurves.pc);
        set(h_incon, 'YData', newCurves.pic);
        drawnow;
    end
    %upFontSize(16,.01)
end
