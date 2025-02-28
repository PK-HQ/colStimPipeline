function mdl=plotNakaRushtonFit3(behavioralData, bitmapData, datastruct, analysisBlockID,...
    mdl, fitParams, x, monkeyName, clusterBlocks, plotAverageFlag,...
    saveFlag, cluster, modelTypeStr, savefilename)
   
    % Number of blocks and conditions
    [nConditions, ~, nBlocks] = size(behavioralData.gaborContrasts(:, :, clusterBlocks));
    if plotAverageFlag==1
        nBlocks=1;
    end
    for block = 1:nBlocks
        % Init figure
        make_it_tight = true;
        hmarg=.15;
        subplot = @(m,n,p) subtightplot (m, n, p, [0.1 0.1], [hmarg hmarg], [0.1 0.1]);
        if ~make_it_tight,  clear subplot;  end
       
        figure('Name', ['Block #', datastruct(analysisBlockID(block)).date]);
        subplot(1,2,1)

        hold on;
        yline(50,'--','LineWidth',1.5,'Color',.4*[1 1 1],'HandleVisibility','off'); hold on;
        xline(0,'--','LineWidth',1.5,'Color',.4*[1 1 1],'HandleVisibility','off'); hold on
        for cond = 1:4%nConditions
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
                    idx=x>=0;
                    xPlot = x(idx);
                    predictedCurve = mdl.mdlBaseline(xPlot, fitParams(block, 1:end-1)); 
                    lineColor = [0 0 0]; % Black for baseline
                    markerFaceColor = [1 1 1]; 
                    edgeColor = 'k';
                    markerType = 'o';
                    xBlock = rmnan(mdl.xBaseline(block, :));
                    yBlock = rmnan(mdl.yBaseline(block, :));
                    [~, idx] = find(xBlock >= 0);
                    xBlock = xBlock(idx);
                    yBlock = yBlock(idx);
                    mdl.thresholdContrast(block, 1) = getThreshold(xPlot, predictedCurve, 70);

                case 2 % Con-Opto
                    if strcmp(modelTypeStr,'bill')
                        xPlot = x; % Positive contrasts for congruent condition
                        predictedCurve = mdl.mdlOpto(xPlot, fitParams(block, 1:end-1)); 
                        xPlot = x(x >= 0); % Positive contrasts for congruent condition
                        predictedCurve=predictedCurve(x >= 0);
                    else
                        xPlot = x(x >= 0); % Positive contrasts for congruent condition
                        predictedCurve = mdl.mdlOpto(xPlot, fitParams(block, 1:end-1)); 
                    end
                    lineColor = [0.9294, 0.1098, 0.1373] * 1.05; % Red for con-opto
                    markerFaceColor = lineColor;
                    edgeColor = 'k';
                    markerType = '^';
                    xBlock = rmnan(mdl.xOpto(block, :));
                    yBlock = rmnan(mdl.yOpto(block, :));
                    [~, idx] = find(xBlock >= 0);
                    xBlock = xBlock(idx);
                    yBlock = yBlock(idx);
            
                case 3 % Incon-Opto
                    if strcmp(modelTypeStr,'bill')
                        xPlot = x; % Positive contrasts for congruent condition
                        predictedCurve = mdl.mdlOpto(xPlot, fitParams(block, 1:end-1)); 
                        xPlot = x(x <= 0); % Neg contrasts for congruent condition
                        predictedCurve=predictedCurve(x <= 0);
                    else
                        xPlot = x(x <= 0);
                        predictedCurve = 100-mdl.mdlOpto(xPlot, fitParams(block, 1:end-1)); 
                        xPlot = -1.* x(x <= 0); % Neg contrasts for congruent condition

                    end
                    lineColor = [0, 0.0941, 0.6627] * 1.25; % Blue for incon-opto
                    markerFaceColor = lineColor;
                    edgeColor = 'k';
                    markerType = 'v';
                    idx=find(rmnan(mdl.xOpto(block, :))<1);
                    xBlock = rmnan(mdl.xOpto(block, :));
                    yBlock = rmnan(mdl.yOpto(block, :));
                    [~, idx] = find(xBlock < 0);
                    xBlock = -1.*xBlock(idx);
                    yBlock = 100-yBlock(idx);
                    
                case 4
                    xBlock = rmnan(mdl.xOpto(block, :));
                    yBlock = rmnan(mdl.yOpto(block, :));
                    [~, idxCon] = find(xBlock >= 0);
                    [~, idxIncon] = find(xBlock <= 0);

                    xBlock = xBlock(idxCon);
                    yBlock = yBlock(idxCon) - fliplr(100-yBlock(idxIncon));

                    markerType = 'square';

                    lineColor =[127, 0, 255]/255; % Gray for combined case
                    markerFaceColor =lineColor;
                    edgeColor = 'k';
                    
                    
                otherwise
                    lineColor = 'g'; % Fallback lineColor
                    markerFaceColor = lineColor;
                    edgeColor = 'g';
                    markerType='.';
            end
            %% Plots
            if cond<=3
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
                elseif sum(zeroConstrastPoint)>0 && plotAverageFlag==1
                    nTrials=[nBlocks*ones(1,numel(zeroConstrastPoint)) nBlocks*ones(1,numel(~zeroConstrastPoint))];
                elseif sum(zeroConstrastPoint)==0
                     nTrials=[20*ones(1,numel(~zeroConstrastPoint))];
                elseif sum(zeroConstrastPoint)==0 && plotAverageFlag==1
                     nTrials=[nBlocks*ones(1,numel(~zeroConstrastPoint))];
                end
                annotateDataPoints(xBlock, yBlock, nTrials, markerFaceColor); hold on;

            end
        end
        if isfield(bitmapData, 'energy') & ~isempty(bitmapData.energy(:,:,clusterBlocks(block)))
            if plotAverageFlag % for plotting the average across all blocks
                bitmapSPD=squeeze(bitmapData.energy(:,:,clusterBlocks));
                bitmapColumns=bitmapData.nColumns(:,clusterBlocks)';
            else
                bitmapSPD=squeeze(bitmapData.energy(:,:,clusterBlocks(block)));
                bitmapColumns=bitmapData.nColumns(:,block)';
            end
        else
            bitmapSPD=[0 0];
        end
        
        bitmapColumnhv=nanmean(bitmapColumns,1);
        bitmapSPDhv=nanmean(bitmapSPD,1);
        bitmapSPDmean=nanmean(bitmapSPD,'all');
        bitmapSPDstd=nanstd(bitmapSPD,[],'all');
        if plotAverageFlag==1
            title({[modelTypeStr ', cluster ' num2str(cluster) ' average'],...
                ['Energy: ' num2str(bitmapSPDhv(1),2) ' & ' num2str(bitmapSPDhv(2),2) ' mW (' num2str(bitmapSPDmean,2) ' \pm ' num2str(bitmapSPDstd,1) ' mW)',...
                ', Columns: ' num2str(bitmapColumnhv(1),2) ' & ' num2str(bitmapColumnhv(2),2)]});
        else
            title({[datastruct(analysisBlockID(block)).date 'R' datastruct(analysisBlockID(block)).run],...
                ['Energy: ' num2str(bitmapSPDhv(1),2) ' & ' num2str(bitmapSPDhv(2),2) ' mW, ',...
                'Columns: ' num2str(bitmapColumnhv(1),2) ' & ' num2str(bitmapColumnhv(2),2)]});
        end

        % Labels etc
        %axis square
        hold off;
        addSkippedTicks(-100,100,12.5,'x'); yticks([0:25:100]); xlim([0 100]); ylim([0 100]); axis square
        % Adding legend after plotting to ensure it covers all conditions
        moveLines()
        h2 = get(gca,'Children');
        legend(h2([end-2:end]), {'Baseline', 'Con-Opto', 'Incon-Opto'}, 'Location', 'southeast',...
            'NumColumns',1,'FontSize',20);
        upFontSize(32, 0.01)

        
        % Add text for biasing
        xOffset=.53;
        yOffset=.05;
        %text('Units', 'Normalized', 'Position', [1 1]-[xOffset yOffset], 'string', 'More biasing', 'color', 'k','FontWeight','bold', 'Fontsize',14)

       % Dual axes
        ax = gca;
        xx = abs([-100:12.5:100]);  
        cellArray = cellstr(num2str(xx(:))); cellArray(2:2:end) = {''}; xticklabels(cellArray)

        ax = gca;
        ylabel('Correct %'); ylim([0 100]); yticks([0:10:100]); set(gca,'ycolor','k') 
        xlabel('Gabor contrast (%)');

        % Plot inset for difference plot
        if cond==4
            subplot(1,2,2)
            
            plot(xBlock, yBlock, 'Color', lineColor, 'LineWidth',2, ...
                'Marker', markerType, 'MarkerSize', 15, 'MarkerFaceColor', markerFaceColor, 'MarkerEdgeColor', edgeColor);

            ylabel('\Delta%');

            yline(0,'--','LineWidth',1.5,'Color',.4*[1 1 1],'HandleVisibility','off'); hold on
            xline(0,'--','LineWidth',1.5,'Color',.4*[1 1 1],'HandleVisibility','off'); hold on
            
            % Set limits if necessary (optional)
            box off
            xlim([0 100]);
            ylim([0 100]);
            addSkippedTicks(-20, 100,10,'y')
            addSkippedTicks(-0,100,12.5,'x')
            set(gca,'FontSize', 10); % Smaller font size for inset        
            set(gca,'linewidth',2)
            upFontSize(32, 0.01)            
            % Add text for biasing
            xOffset=.6;
            yOffset=.05;
            text('Units', 'Normalized', 'Position', [1 1]-[xOffset yOffset], 'string', 'Correct biasing', 'color', 'k','FontWeight','normal', 'Fontsize',14)
            title('Con-Incon optostim', 'FontWeight','normal')
            % Add text for biasing
            axis square
        end
        
        
        % Add params table
        if isnan(fitParams(block,1))
            dat=round([50 fitParams(block,[2:4 end])],1);
        else
            dat=round([fitParams(block,1) fitParams(block,[2:4 end])],1);
            headers=[mdl.headers(1) mdl.headers([2:4 end])];
        end
        T = array2table(dat,'VariableNames',headers);
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
           'Units','Normalized','Position',[.18 -.65 1 1]);
        
        if plotAverageFlag==1
            [nConditions, ~, nBlocks] = size(behavioralData.gaborContrasts(:, :, clusterBlocks));
            nBlocks=1;
            block=1;
        end
        savePDF(savefilename, monkeyName, 1, block, nBlocks)

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


%{
if billplot==1
    clf
    plot(xPlot, predictedCurve, 'Color', lineColor, 'LineWidth', 3, 'HandleVisibility', 'on'); hold on;
end
%}
end

function thresholdContrast=getThreshold(xData,yData,target);
% Find threshold value

% Find the index of the closest value
[~, index] = min(abs(yData - target));

% Get the closest value
thresholdContrast = xData(index);
end