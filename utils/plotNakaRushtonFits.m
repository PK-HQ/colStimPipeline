function plotNakaRushtonFits(behavioralData, bitmapData, datastruct, analysisBlockID, fitParams, x, saveFlag, savefilename, monkeyName)
   
    % Number of blocks and conditions
    [nBlocks, nConditions, ~] = size(fitParams);
    
    % Define contrast range for plotting smooth curves
    xPlot = linspace(min(x(:)), max(x(:)), 100);
    
    for block = 1:nBlocks
        figure('Name', ['Block #', datastruct(analysisBlockID(block)).date]);
        hold on;
        yline(50,'--','LineWidth',2,'Color',.3*[1 1 1],'HandleVisibility','off'); hold on;
        for cond = 1:nConditions
            % Extract fitted parameters for current condition and block
            beta = fitParams(block, cond, 1);
            n = fitParams(block, cond, 2);
            C50 = fitParams(block, cond, 3);
            deltax = fitParams(block, cond, 4);
            aicc = fitParams(block, cond, 5);
            
            % Define the Naka-Rushton functions with the fitted parameters
            switch cond
                case 1 % Baseline
                    Rmax = 100 - beta; % Assuming Rmax is derived from beta
                    predictedCurve = (Rmax .* xPlot.^n) ./ (C50^n + xPlot.^n) + beta;
                    lineColor = [0 0 0]; % Black for baseline
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
            %% Plots
            % Plot Naka-Rushton fit
            plot(xPlot, predictedCurve, 'Color', lineColor, 'LineWidth', 3, 'HandleVisibility', 'on'); hold on;
            
            % Plot scatterplot
            xBlock=squeeze(behavioralData.gaborContrast(cond,3,:, block)); xBlock=xBlock(~isnan(xBlock));
            yBlock=squeeze(behavioralData.percentageCorrect(cond,3,:,block)); yBlock=yBlock(~isnan(yBlock));
            semY=std(yBlock', 'omitnan') / sqrt(length(yBlock));
            counts=size(yBlock,1);
            if size(yBlock,2) > 1 % if no sem, don't shade
                patchSaturationVal = 0.2;
            else
                semY=zeros(counts,1);
                patchSaturationVal = 0;
            end
            markerSize = 25;
            patchSaturationVal=1;
            shadedErrorBar(xBlock, yBlock, semY, 'patchSaturation', patchSaturationVal, 'lineprops', ...
                           {'Color', lineColor, 'LineStyle', 'none', 'LineWidth', 3, 'Marker', markerType, ...
                            'MarkerFaceColor', markerFaceColor, 'MarkerEdgeColor', edgeColor, 'MarkerSize', markerSize}); hold on;

            % Annotate markers
            nTrials=[40*ones(1,sum(xBlock==0)) 20*ones(1,sum(xBlock>0))];
            annotateDataPoints(xBlock, yBlock, nTrials, markerFaceColor); hold on;
            
            % Display fitParams
            lineOrder=[2 1 3 4];
            drawParamText([beta, n, C50, deltax, aicc], lineOrder, lineColor, cond)

            % Customizations for readability
            xlabel('Absolute gabor contrast (%)');
            ylabel('Percent correct (%)');

        end
        bitmapSPD=squeeze(bitmapData.adjustedSPD_uW(:,:,block));
        title({[datastruct(analysisBlockID(block)).date 'R' datastruct(analysisBlockID(block)).run],...
            ['Adj. SPD: ' num2str(mean(bitmapSPD),2) ' \pm ' num2str(std(bitmapSPD),1) ' \muW']});
        % Labels etc
        axis square
        hold off;
        xlim([0 100]); ylim([0 100]); addSkippedTicks(0,100,12.5,'x');addSkippedTicks(0,100,12.5,'y');
        % Adding legend after plotting to ensure it covers all conditions
        moveLines()
        h2 = get(gca,'Children');
        legend(h2([end-2:end]), {'Baseline', 'Con-Opto', 'Incon-Opto'}, 'Location', 'southeast',...
            'NumColumns',1,'FontSize',20);
        %h = get(gca,'Children');
        %resortHandles(20, 'down')
        upFontSize(32, 0.01)
        % Save summary of all blocks
        savePDF(savefilename, monkeyName, saveFlag, block)
    end
end

function drawParamText(fitParams, lineOrder, lineColor, lineNo)
% Plot stuff
xpos=.025;
ypos=.18;
% Split the fit params
beta=fitParams(1);
rmax=100-fitParams(1);
exponent=fitParams(2);
c50=fitParams(3);
deltax=fitParams(4);
aicc=fitParams(5);
pseudoR2=[];
if  lineNo==1 % Add header
    % State fit parameters
    headerText = sprintf(...
     [' {\\bf\\beta}        {\\bfR_{max}}        {\\bfC_{50}}          {\\bfn}          {\\bf\\Deltax}         {\\bfAICc}']);
     text(xpos,ypos,headerText,'Units','normalized', 'VerticalAlignment','top', 'HorizontalAlignment','left','Color','k')
end
if lineNo<=3
    % Add param text
    newlineStr=repmat('\n',1,lineOrder(lineNo));
    paramsText=sprintf(...
    [newlineStr '%.0f         %.0f            %02.f          %0.1f         %02.f           %03.f'], ...
    beta, rmax, c50, exponent, deltax, aicc);
    text(xpos,ypos-.02,paramsText,'Units','normalized', 'VerticalAlignment','top', 'HorizontalAlignment','left','Color',lineColor);
    upFontSize(24,0.01)
end
end

function resortHandles(linesToMove, position)
% Get the current handles of all children
h = get(gca, 'Children');

% Extract the elements you want to move to the bottom
elementsToMove = h(linesToMove);

% Remove these elements from the original array
h(linesToMove) = [];

% Concatenate the removed elements at the end of the array
switch strcmp(position,'down')
    case 1
        h = [h; elementsToMove];
    case 0
        h = [elementsToMove;h];
end

% Set the new order of children
set(gca, 'Children', h);
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