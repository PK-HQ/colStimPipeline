function plotNakaRushtonFit3(behavioralData, bitmapData, datastruct, analysisBlockID, mdl, fitParams, x, monkeyName, saveFlag, savefilename)
   
    % Number of blocks and conditions
    [nConditions, ~, nBlocks] = size(behavioralData.gaborContrasts);
        
    for block = 1:nBlocks
        figure('Name', ['Block #', datastruct(analysisBlockID(block)).date]);
        hold on;
        yline(50,'--','LineWidth',1.5,'Color',.4*[1 1 1],'HandleVisibility','off'); hold on;
        xline(0,'--','LineWidth',1.5,'Color',.4*[1 1 1],'HandleVisibility','off'); hold on
        for cond = 1:nConditions
            % Extract fitted parameters for current condition and block
            beta = fitParams(block, 1);
            n = fitParams(block, 2);
            C50 = fitParams(block, 3);
            if cond==1
                deltax = 0;
            elseif cond==2
                deltax = fitParams(block, end-2);
            elseif cond==3
                deltax = fitParams(block, end-1);
            end
            aicc = fitParams(block, end);
            
            % Define the Naka-Rushton functions with the fitted parameters
            betaBL=50;

            switch cond
                case 1 % Baseline
                    xPlot=x;
                    beta=50;
                    Rmax = 100 - beta; % Assuming Rmax is derived from beta
                    predictedCurve = mdl.nakaRushtonBaseline(xPlot,fitParams(block,1:end-1)); %(Rmax .* (xPlot-deltax).^n) ./ (C50^n + (xPlot-deltax).^n) + beta;
                    lineColor = [0 0 0]; % Black for baseline
                    markerFaceColor = [1 1 1]; 
                    edgeColor = 'k';
                    markerType='o';
                    
                    xBlock=rmnan(mdl.xBaseline(block,:));
                    yBlock=rmnan(mdl.yBaseline(block,:));

                case 2 % Con-Opto
                    xPlot=x(x>=0);
                    predictedCurve = mdl.nakaRushtonOpto(xPlot,fitParams(block,1:end-1)); %((100-beta) .* (xPlot-deltax).^n) ./ (C50^n + (xPlot-deltax).^n) + beta;
                    lineColor = [0.9294, 0.1098, 0.1373]*1.05; % Red for con-opto
                    markerFaceColor = lineColor;
                    edgeColor = 'k';
                    markerType='^';
                    xBlock=rmnan(mdl.xOpto(block,:));
                    yBlock=rmnan(mdl.yOpto(block,:));
                    [~,idx]=find(xBlock>=0);
                    xBlock=xBlock(idx);
                    yBlock=yBlock(idx);
                case 3 % Incon-Opto
                    xPlot=x(x<=0);
                    predictedCurve = mdl.nakaRushtonOpto(xPlot,fitParams(block,1:end-1)); %(100-(100-(beta)) .* (xPlot+deltax).^n) ./ (C50^n + (xPlot+deltax).^n) + (100-beta);
                    lineColor = [0, 0.0941, 0.6627]*1.25; % Blue for incon-opto
                    markerFaceColor = lineColor;
                    edgeColor = 'k';
                    markerType='v';
                    xBlock=rmnan(mdl.xOpto(block,:));
                    yBlock=rmnan(mdl.yOpto(block,:));
                    [~,idx]=find(xBlock<0);
                    xBlock=xBlock(idx);
                    yBlock=yBlock(idx);
                otherwise
                    lineColor = 'g'; % Fallback lineColor
                    markerFaceColor = lineColor;
                    edgeColor = 'g';
                    markerType='.';
            end
            %% Plots
            % Plot line fit
            plot(xPlot, predictedCurve, 'Color', lineColor, 'LineWidth', 3, 'HandleVisibility', 'on'); hold on;
            
            % Plot scatterplot
            semY=std(yBlock, 'omitnan') / sqrt(length(yBlock));
            counts=size(yBlock,2);
            if size(yBlock,1) > 1 % if no sem, don't shade
                patchSaturationVal = 0.2;
            else
                semY=zeros(counts,1);
                patchSaturationVal = 0;
            end
            
            % Add data points and shaded error bar
            markerSize = 25;
            patchSaturationVal=1;
            shadedErrorBar(xBlock, yBlock, semY', 'patchSaturation', patchSaturationVal, 'lineprops', ...
                           {'Color', lineColor, 'LineStyle', 'none', 'LineWidth', 3, 'Marker', markerType, ...
                            'MarkerFaceColor', markerFaceColor, 'MarkerEdgeColor', edgeColor, 'MarkerSize', markerSize}); hold on;

            % Add datapoints's count annotation
            zeroConstrastPoint=xBlock==0;
            if sum(zeroConstrastPoint)>0
                nTrials=[20*ones(1,numel(zeroConstrastPoint)) 20*ones(1,numel(~zeroConstrastPoint))];
            else
                 nTrials=[20*ones(1,numel(~zeroConstrastPoint))];
            end
            annotateDataPoints(xBlock, yBlock, nTrials, markerFaceColor); hold on;
            
            xlabel('Gabor contrast (%)');

        end
        if isfield(bitmapData, 'adjustedSPD_uW') & ~isempty(bitmapData.adjustedSPD_uW(:,:,block))
            bitmapSPD=squeeze(bitmapData.adjustedSPD_uW(:,:,block));
        else
            bitmapSPD=[0 0];
        end
        
        title({[datastruct(analysisBlockID(block)).date 'R' datastruct(analysisBlockID(block)).run],...
            ['Adj. SPD: ' num2str(mean(bitmapSPD),2) ' \pm ' num2str(std(bitmapSPD),1) ' \muW']});
        % Labels etc
        axis square
        hold off;
        addSkippedTicks(-100,100,12.5,'x'); yticks([0:25:100]); xlim([-100 100]); ylim([0 100]); 
        % Adding legend after plotting to ensure it covers all conditions
        moveLines()
        h2 = get(gca,'Children');
        legend(h2([end-2:end]), {'Baseline', 'Con-Opto', 'Incon-Opto'}, 'Location', 'northwest',...
            'NumColumns',1,'FontSize',20);
        upFontSize(32, 0.01)

        % Add params table
        header={'Beta','n','C50','DeltaX','AICc'};
        %header = {'\Beta';'n';'C_{50}';'\DeltaX_{H-opto}';'\DeltaX_{V-opto}'};
        if isnan(fitParams(block,1))
            dat=round([50 fitParams(block,[2:4 end])],1);
        else
            dat=round([fitParams(block,1) fitParams(block,[2:4 end])],1);
        end
        T = array2table(dat,'VariableNames',header);
        % Get the table in string form.
        TString = evalc('disp(T)');
        % Use TeX Markup for bold formatting and underscores.
        TString = strrep(TString,'<strong>','\bf');
        TString = strrep(TString,'</strong>','\rm');
        TString = strrep(TString,'_','\_');
        % Get a fixed-width font.
        FixedWidth = get(0,'FixedWidthFontName');
        % Output the table using the annotation command.
        annotation(gcf,'Textbox','String',TString,'Interpreter','Tex','FontName',FixedWidth,'FontSize',12,'EdgeColor','none',...
            'Units','Normalized','Position',[.72 -.85 1 1]);
        
        % Dual axes
        ax = gca;
        xx = abs([-100:12.5:100]);  
        cellArray = cellstr(num2str(xx(:))); cellArray(2:2:end) = {''}; xticklabels(cellArray)
        yyaxis left; ylabel('Correct, incongruent (%)'); 
        ax = gca;
        ax.YTickLabel = flipud(ax.YTickLabel);
        yyaxis right; ylabel('Correct, congruent (%)'); ylim([0 100]); yticks([0:25:100]); set(gca,'ycolor','k') 
        
        % Add text for biasing
        xOffset=.569;
        yOffset=.02;
        text('Units', 'Normalized', 'Position', [1 1]-[xOffset yOffset], 'string', 'More bias', 'color', 'k','FontWeight','bold', 'Fontsize',14)
        
        savePDF(savefilename, monkeyName, saveFlag, block, nBlocks)

    end
end

function moveLines()
% Get the current handles of all children
h = get(gca, 'Children');

% Initialize an array to hold the indices of lines to move
indicesToMove = [];

% Loop through each child handle
for idx = 1:length(h)
    % Check if the handle is a line object and isline
    if isa(h(idx), 'matlab.graphics.chart.primitive.Line') && ~strcmp(get(h(idx), 'LineStyle'),'none')
        % Add the index to the list
        indicesToMove = [indicesToMove, idx];
    end
end

% Extract the elements you want to move
elementsToMove = h(fliplr(indicesToMove));

% Remove these elements from the original array
h(indicesToMove) = [];

% Concatenate the removed elements at the end of the array
h = [h; elementsToMove];

% Set the new order of children
set(gca, 'Children', h);

end