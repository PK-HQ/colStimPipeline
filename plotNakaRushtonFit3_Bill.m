function mdl = plotmdlFit3_Bill(behavioralData, bitmapData, datastruct, analysisBlockID,...
    mdl, fitParams, x, monkeyName, clusterBlocks, plotAverageFlag,...
    saveFlag, cluster, modelTypeStr, savefilename)
   
    % Number of blocks and conditions
    [nConditions, ~, nBlocks] = size(behavioralData.gaborContrasts(:, :, clusterBlocks));
    if plotAverageFlag == 1
        nBlocks = 1;
    end

    for block = 1:nBlocks
        % Initialize figure
        figure('Name', ['Block #', datastruct(analysisBlockID(block)).date]);
        subplot(3, 2, [1:4]);
        hold on;
        
        % Reference lines
        yline(50, '--', 'LineWidth', 1.5, 'Color', 0.4 * [1 1 1], 'HandleVisibility', 'off'); 
        xline(0, '--', 'LineWidth', 1.5, 'Color', 0.4 * [1 1 1], 'HandleVisibility', 'off'); 

        % Loop through conditions
        for cond = 1:3
            % Extract fitted parameters
            alpha = fitParams(block, 1);  % Weight of iso-tuned visual input
            beta = fitParams(block, 2);   % Weight of ortho-tuned visual input
            lambda = fitParams(block, 3); % Weight of iso-tuned inputs in normalization pool
            omega = fitParams(block, 4);  % Weight of cross-orientation inputs in normalization pool
            g0 = fitParams(block, 5);     % Normalization constant
            n = fitParams(block, 6);      % Spiking exponent
            rmax = fitParams(block, 7);   % Max response
            eta = fitParams(block, 8);    % Relative weight of optostim effects
            o = fitParams(block, 9);      % Opto-stim effective contrast

            % Baseline prediction (no optostim)
            if cond == 1
                xPlot=x;
                predictedCurve = mdl.mdlBaseline(xPlot, [alpha, beta, lambda, omega, g0, n, rmax, eta, 0]);
                lineColor = [0 0 0]; % Black for baseline
                markerFaceColor = [1 1 1]; 
                markerType = 'o';
                xBlock = rmnan(mdl.xBaseline(block, :));
                yBlock = rmnan(mdl.yBaseline(block, :));
            end

            % Congruent optostim prediction
            if cond == 2
                xPlot=x(x >= 0);
                predictedCurve = mdl.mdlOpto(xPlot, [alpha, beta, lambda, omega, g0, n, rmax, eta, o]);
                lineColor = [0.9294, 0.1098, 0.1373] * 1.05; % Red for congruent
                markerFaceColor = lineColor;
                markerType = '^';
                xBlock = rmnan(mdl.xOpto(block, :));
                yBlock = rmnan(mdl.yOpto(block, :));
                xBlock = xBlock(xBlock >= 0);
                yBlock = yBlock(xBlock >= 0);
            end

            % Incongruent optostim prediction
            if cond == 3
                xPlot=x(x <= 0);
                predictedCurve = mdl.mdlOpto(xPlot, [alpha, beta, lambda, omega, g0, n, rmax, eta, o]);
                lineColor = [0, 0.0941, 0.6627] * 1.25; % Blue for incongruent
                markerFaceColor = lineColor;
                markerType = 'v';
                xBlock = rmnan(mdl.xOpto(block, :));
                yBlock = rmnan(mdl.yOpto(block, :));
                xBlock = xBlock(xBlock <= 0);
                yBlock = yBlock(xBlock <= 0);
            end

            % Plot predicted curves
            plot(xPlot, predictedCurve, ...
                'Color', lineColor, 'LineWidth', 3, 'HandleVisibility', 'on'); hold on;

            % Plot data points
            scatter(xBlock, yBlock, 50, 'Marker', markerType, 'MarkerFaceColor', markerFaceColor, 'MarkerEdgeColor', 'k'); hold on;
        end

        % Set axis limits and labels
        ylim([0 100]);
        yticks(0:25:100);
        xlim([-100 100]);
        xticks(-100:25:100);
        xlabel('Gabor contrast (%)');
        ylabel('% correct');
        title(['Bill Model Fit: Block ' num2str(block)]);

        % Save the figure
        if saveFlag
            saveas(gcf, [savefilename, '_block', num2str(block), '.png']);
        end
    end
end

