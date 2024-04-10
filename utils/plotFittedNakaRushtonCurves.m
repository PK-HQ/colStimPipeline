function plotFittedNakaRushtonCurves(behavioralData, fitParams, x)
    % Number of blocks and conditions
    [nBlocks, nConditions, ~] = size(fitParams);
    
    % Define contrast range for plotting smooth curves
    xPlot = linspace(min(x(:)), max(x(:)), 100);
    
    for b = 1:nBlocks
        figure('Name', ['Block ', num2str(b)]);
        hold on;
        
        for c = 1:nConditions
            % Extract fitted parameters for current condition and block
            beta = fitParams(b, c, 1);
            n = fitParams(b, c, 2);
            C50 = fitParams(b, c, 3);
            
            % Define the Naka-Rushton functions with the fitted parameters
            switch c
                case 1 % Baseline
                    Rmax = 100 - beta; % Assuming Rmax is derived from beta
                    predictedCurve = (Rmax .* xPlot.^n) ./ (C50^n + xPlot.^n) + beta;
                    lineColor = 'k'; % Black for baseline
                    markerFaceColor = [1 1 1]; 
                    edgeColor = 'k';
                    markerType='o';
                case 2 % Con-Opto
                    predictedCurve = ((100-beta) .* xPlot.^n) ./ (C50^n + xPlot.^n) + beta;
                    lineColor = [0.9294, 0.1098, 0.1373]*1.05; % Red for con-opto
                    markerFaceColor = lineColor;
                    edgeColor = 'k';
                    markerType='^';
                case 3 % Incon-Opto
                    predictedCurve = ((100-(beta)) .* xPlot.^n) ./ (C50^n + xPlot.^n) + (beta);
                    lineColor = [0, 0.0941, 0.6627]*1.25; % Blue for incon-opto
                    markerFaceColor = lineColor;
                    edgeColor = 'k';
                    markerType='v';
                otherwise
                    lineColor = 'g'; % Fallback lineColor
                    markerFaceColor = lineColor;
                    edgeColor = 'k';
                    markerType='.';
            end
            
            % Plot NR fit
            plot(xPlot, predictedCurve, 'Color', lineColor, 'LineWidth', 2); hold on;
            
            % Plot scatterplot
            xBlock=squeeze(behavioralData.gaborContrast(c,3,:, b)); xBlock=xBlock(~isnan(xBlock));
            yBlock=squeeze(behavioralData.percentageCorrect(c,3,:,b)); yBlock=yBlock(~isnan(yBlock));
            semY=std(yBlock', 'omitnan') / sqrt(length(yBlock));
            counts=size(yBlock,1);
            if size(yBlock,2) > 1 % if no sem, don't shade
                patchSaturationVal = 0.2;
            else
                semY=zeros(counts,1);
                patchSaturationVal = 0;
            end

                markerSize = 20;
                shadedErrorBar(xBlock, yBlock, semY, 'patchSaturation', patchSaturationVal, 'lineprops', ...
                               {'Color', lineColor, 'LineStyle', 'none', 'LineWidth', 2, 'Marker', markerType, ...
                                'MarkerFaceColor', markerFaceColor, 'MarkerEdgeColor', 'k', 'MarkerSize', markerSize});

                % Assuming markerFaceColor is defined, e.g., markerFaceColor = [0.2, 0.5, 0.7];
                % Convert markerFaceColor to grayscale using luminance to decide on text lineColor
                %annotateDataPoints(xBlock, yBlock, counts, markerFaceColor);
            
            % Customizations for readability
            xlabel('Contrast (%)');
            ylabel('Predicted Response (%)');
            title(['Block ', num2str(b)]);
            %grid off;
        end
        
        % Adding legend after plotting to ensure it covers all conditions
        legend({'Baseline', 'Con-Opto', 'Incon-Opto'}, 'Location', 'best');
        axis square
        yline(50,'--','LineWidth',2,'Color',.3*[1 1 1])
        hold off;
        xlim([0 100]); ylim([0 100])
        upFontSize(24, 0.01)
    end
end
