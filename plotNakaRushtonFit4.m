function mdl=plotNakaRushtonFit4(behavioralData, bitmapData, datastruct, analysisBlockID,...
    mdl, fitParams, x, monkeyName, clusterBlocks, plotAverageFlag,...
    saveFlag, cluster, modelTypeStr, savefilename)
    endIdx=size(mdl.headers,2);
    % Number of blocks and conditions
    [nConditions, ~, nBlocks] = size(behavioralData.gaborContrasts(:, :, clusterBlocks));
    if plotAverageFlag==1
        nBlocks=1;
    end
    for block = nBlocks
        % Init figure
        dat=[];
        make_it_tight = true;
        hmarg=.15;
        subplot = @(m,n,p) subtightplot (m, n, p, [0.1 0.1], [hmarg hmarg], [0.1 0.1]);
        if ~make_it_tight,  clear subplot;  end
       
        figure('Name', ['Block #', datastruct(analysisBlockID(clusterBlocks(block))).date]);
        subplot(1,2,1)

        hold on;
        yline(50,'--','LineWidth',1.5,'Color',.4*[1 1 1],'HandleVisibility','off'); hold on;
        %xline(0,'--','LineWidth',1.5,'Color',.4*[1 1 1],'HandleVisibility','off'); hold on
        for cond = 1:nConditions+1
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
                    xPlot = x;
                    predictedCurve = mdl.mdlBaseline(xPlot, fitParams(block, 1:end-1)); 
                    lineColor = [0 0 0]; % Black for baseline
                    markerFaceColor = [1 1 1]; 
                    edgeColor = 'k';
                    markerType = 'o';
                    xBlock = rmnan(mdl.xBaseline(block, :));
                    yBlock = rmnan(mdl.yBaseline(block, :));
                    [~, idx] = find(xBlock>=0);
                    xBlock = xBlock(idx);
                    yBlock = yBlock(idx);
                    mdl.thresholdContrast(block, 1) = getThreshold(xPlot, predictedCurve, 70);
                    predictedCurve = replaceNanSections(predictedCurve);
                    mdl.xFitted(cond,:,block)=padArray(xPlot,400,2,nan);
                    mdl.yFitted(cond,:,block)=padArray(predictedCurve,400,2,nan);
                    mdl.xBlock(cond,:,block)=padArray(xBlock, 120, 2, nan);
                    mdl.yBlock(cond,:,block)=padArray(yBlock, 120, 2, nan);
                    %get AUC
                    idxx=mdl.xFitted(cond,:,block)>=0 & ~isnan(mdl.yFitted(cond,:,block));
                    contr=mdl.xFitted(cond,idxx,block);
                    curve=mdl.yFitted(cond,idxx,block);
                    mdl.fittedParams(block,endIdx) = trapz(contr, curve) / (max(contr) - min(contr));
                case 2 % Con-Opto
                    if strcmp(modelTypeStr,'bill')
                        xPlot = x; % Positive contrasts for congruent condition
                        predictedCurve = mdl.mdlOpto(xPlot, fitParams(block, 1:end-1)); 
                        xPlot = x(x >= 0); % Positive contrasts for congruent condition
                        predictedCurve=predictedCurve(x >= 0);
                    else
                        xPlot = x; % Positive contrasts for congruent condition
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
                    predictedCurve = replaceNanSections(predictedCurve);
                    mdl.xFitted(cond,:,block)=padArray(xPlot,400,2,nan);
                    mdl.yFitted(cond,:,block)=padArray(predictedCurve,400,2,nan);
                    mdl.xBlock(cond,:,block)=padArray(xBlock, 120, 2, nan);
                    mdl.yBlock(cond,:,block)=padArray(yBlock, 120, 2, nan);
                    %get AUC
                    idxx=mdl.xFitted(cond,:,block)>=0 & ~isnan(mdl.yFitted(cond,:,block));
                    contr=mdl.xFitted(cond,idxx,block);
                    curve=mdl.yFitted(cond,idxx,block);
                    mdl.fittedParams(block,endIdx+1) = trapz(contr, curve) / (max(contr) - min(contr));
                case 3 % Incon-Opto
                    if strcmp(modelTypeStr,'bill')
                        xPlot = x; % Positive contrasts for congruent condition
                        predictedCurve = 100-mdl.mdlOpto(xPlot, fitParams(block, 1:end-1)); 
                        xPlot = -1*x(x <= 0); % Neg contrasts for congruent condition
                        predictedCurve=predictedCurve(x <= 0);
                    else
                        xPlot = x;
                        predictedCurve = 100-mdl.mdlOpto(xPlot, fitParams(block, 1:end-1)); 
                        xPlot = -1.* xPlot; % Neg contrasts for congruent condition
                    end
                    lineColor = [0, 0.0941, 0.6627] * 1.25; % Blue for incon-opto
                    markerFaceColor = lineColor;
                    edgeColor = 'k';
                    markerType = 'v';
                    idx=find(rmnan(mdl.xOpto(block, :))<1);
                    xBlock = rmnan(mdl.xOpto(block, :));
                    yBlock = rmnan(mdl.yOpto(block, :));
                    [~, idx] = find(xBlock <= 0);
                    xBlock = -1.*xBlock(idx);
                    yBlock = 100-yBlock(idx);
                    predictedCurve = replaceNanSections(predictedCurve);
                    mdl.xFitted(cond,:,block)=padArray(xPlot,400,2,nan);
                    mdl.yFitted(cond,:,block)=padArray(predictedCurve,400,2,nan);
                    mdl.xBlock(cond,:,block)=padArray(xBlock, 120, 2, nan);
                    mdl.yBlock(cond,:,block)=padArray(yBlock, 120, 2, nan);

                case 4
                    idxPos=mdl.xFitted(2,:, block)>=0;
                    idxNeg=mdl.xFitted(3,:, block)>=0;

                    
                    xBlock = mdl.xBlock(2,:, block);
                    yBlock = rmnan(mdl.yBlock(2,:, block))-rmnan(fliplr(mdl.yBlock(3,:, block)));

                    xPlot = mdl.xFitted(2,idxPos, block);
                    predictedCurve = mdl.yFitted(2,idxPos, block)-fliplr(mdl.yFitted(3,idxNeg, block));

                    markerType = 'square';
                    ylimMax=50;%roundup(max(abs(yBlock(:))),10);

                    lineColor =[127, 0, 255]/255; % Gray for combined case
                    markerFaceColor =lineColor;
                    edgeColor = 'k';
                    %TEMPORARY FIX
                    predictedCurve = replaceNanSections(predictedCurve);
                    mdl.xFitted(cond,:,block)=padArray([nan(size(xPlot,1),size(xPlot,2)), xPlot],400,2,nan);
                    mdl.yFitted(cond,:,block)=padArray([nan(size(predictedCurve,1),size(predictedCurve,2)), predictedCurve],400,2,nan);
                    mdl.xBlock(cond,:,block)=padArray(xBlock, 120, 2, nan);
                    mdl.yBlock(cond,:,block)=padArray(yBlock, 120, 2, nan);
                    %get AUC
                    idxx=mdl.xFitted(cond,:,block)>=0 & ~isnan(mdl.yFitted(cond,:,block));
                    contr=mdl.xFitted(cond,idxx,block);
                    curve=mdl.yFitted(cond,idxx,block);
                    mdl.fittedParams(block,endIdx+3) = trapz(contr, curve) / (max(contr) - min(contr));
                    
                    mdl.fittedParams(block,endIdx+2) = mdl.fittedParams(block,end-3) - mdl.fittedParams(block,end);
                otherwise
                    lineColor = 'g'; % Fallback lineColor
                    markerFaceColor = lineColor;
                    edgeColor = 'g';
                    markerType='.';
            end


            %% Plots
            if cond<=3
                % Plot line fit
                plot(mdl.xFitted(cond,:,block), mdl.yFitted(cond,:,block), 'Color', lineColor, 'LineWidth', 3, 'HandleVisibility', 'on'); hold on;
                upFontSize(32, 0.01); legend; axis square

                % Plot scatterplot
                semY=std(mdl.yBlock(cond,:,block), 'omitnan') / sqrt(length(yBlock));
                counts=size(mdl.yBlock(cond,:,block),2);
                if size(mdl.yBlock(cond,:,block),1) > 1 % if no sem, don't shade
                    patchSaturationVal = 0.2;
                else
                    semY=zeros(counts,1);
                    patchSaturationVal = 0;
                end

                %{
                % Add data points and shaded error bar
                markerSize = 25;
                patchSaturationVal=1;
                shadedErrorBar(mdl.xBlock(cond,:,block), mdl.yBlock(cond,:,block), semY', 'patchSaturation', patchSaturationVal, 'lineprops', ...
                               {'Color', lineColor, 'LineStyle', 'none', 'LineWidth', 3, 'Marker', markerType, ...
                                'MarkerFaceColor', markerFaceColor, 'MarkerEdgeColor', edgeColor, 'MarkerSize', markerSize}); hold on;
                %}

                % Add data points and shaded error bar
                markerSize = 16;
                patchSaturationVal=1;
                % Data points
                shadedErrorBar(mdl.xBlock(cond,:,block), mdl.yBlock(cond,:,block), semY', 'patchSaturation', patchSaturationVal, 'lineprops', ...
                               {'Color', lineColor, 'LineStyle', 'none', 'LineWidth', 3, 'Marker', markerType, ...
                                'MarkerFaceColor', markerFaceColor, 'MarkerEdgeColor', edgeColor, 'MarkerSize', markerSize}); hold on;
                % Average
                barLength=45:50;
                plot(barLength, repmat(nanmean(mdl.yBlock(cond,:,block)),1,numel(barLength)), '-', 'Color', lineColor, 'LineWidth', 3, 'HandleVisibility', 'off')
                %plot(100, nanmean(mdl.yBlock(cond,:,block)),'Color', lineColor, 'LineStyle', 'none', 'LineWidth', 2.5, 'Marker', markerType, ...
                %                'MarkerFaceColor', markerFaceColor, 'MarkerEdgeColor', edgeColor, 'MarkerSize', markerSize+6); hold on;
            

                % Add datapoints's count annotation
                zeroConstrastPoint=xBlock==0;
                if sum(zeroConstrastPoint)>0 % data contains 0-point
                    nTrials=[20*ones(1,numel(zeroConstrastPoint)) 20*ones(1,numel(~zeroConstrastPoint))];
                    if cond>1
                            nTrials=[40*zeroConstrastPoint + 20*~zeroConstrastPoint];
                    end
                elseif sum(zeroConstrastPoint)>0 && plotAverageFlag==1 % data contains 0-point & average
                    nTrials=[nBlocks*ones(1,numel(zeroConstrastPoint)) nBlocks*ones(1,numel(~zeroConstrastPoint))];
                    if cond>1
                            nTrials=[40*zeroConstrastPoint + 20*~zeroConstrastPoint];
                    end
                elseif sum(zeroConstrastPoint)==0
                     nTrials=[20*ones(1,numel(~zeroConstrastPoint))];
                elseif sum(zeroConstrastPoint)==0 && plotAverageFlag==1
                     nTrials=[nBlocks*ones(1,numel(~zeroConstrastPoint))];
                end
                %annotateDataPoints(mdl.xBlock(cond,:,block), mdl.yBlock(cond,:,block), nTrials, markerFaceColor); hold on;
            
                
                if cond==3
                    [combinedBLStr,mdl.combinedBL(block)]=checkMergedBaseline(datastruct,analysisBlockID,block);

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
                        %{
                        title({[datastruct(analysisBlockID(clusterBlocks(block))).date 'R' datastruct(analysisBlockID(block)).run ' (' combinedBLStr ')'],...
                            ['Energy: ' num2str(bitmapSPDhv(1),2) ' & ' num2str(bitmapSPDhv(2),2) ' mW, ',...
                            'Columns: ' num2str(bitmapColumnhv(1),2) ' & ' num2str(bitmapColumnhv(2),2)]});
                        %}
                    end
                    % Labels etc
                    %axis square
                    xlim([0 50]); ylim([0 100]); xticks(0:12.5:100);addSkippedTicks(0,50,5,'x'); addSkippedTicks(0,100,10,'y'); axis square
                    % Adding legend after plotting to ensure it covers all conditions
                    moveLines()
                    h2 = get(gca,'Children');
                    legend
                    legend(h2([end-2:end]), {'Baseline', 'Con-Opto', 'Incon-Opto','Con-Incon O'}, 'Location', 'east',...
                       'NumColumns',1,'FontSize',32);
                    upFontSize(32, 0.01)
            
                    
                    % Add text for biasing
                    xOffset=.53;
                    yOffset=.05;
                    %text('Units', 'Normalized', 'Position', [1 1]-[xOffset yOffset], 'string', 'More biasing', 'color', 'k','FontWeight','bold', 'Fontsize',14)
                    ax = gca;
                    ylabel('Correct (%)'); set(gca,'ycolor','k') 
                    xlabel('Gabor contrast (%)');
                end
            elseif cond==4
                %subplot(1,2,2)
                % Plot line fit
                plot(mdl.xFitted(cond,:,block), mdl.yFitted(cond,:,block), 'Color', lineColor, 'LineWidth', 3, 'HandleVisibility', 'on'); hold on;
                upFontSize(32, 0.01); legend; axis square

                % Plot scatterplot
                semY=std(mdl.yBlock(cond,:,block), 'omitnan') / sqrt(length(yBlock));
                counts=size(mdl.yBlock(cond,:,block),2);
                if size(mdl.yBlock(cond,:,block),1) > 1 % if no sem, don't shade
                    patchSaturationVal = 0.2;
                else
                    semY=zeros(counts,1);
                    patchSaturationVal = 0;
                end

                % Add data points and shaded error bar
                markerSize = 16;
                patchSaturationVal=1;
                % Data points
                shadedErrorBar(mdl.xBlock(cond,:,block), mdl.yBlock(cond,:,block), semY', 'patchSaturation', patchSaturationVal, 'lineprops', ...
                               {'Color', lineColor, 'LineStyle', 'none', 'LineWidth', 3, 'Marker', markerType, ...
                                'MarkerFaceColor', markerFaceColor, 'MarkerEdgeColor', edgeColor, 'MarkerSize', markerSize}); hold on;
                % Average
                plot(100, nanmean(mdl.yBlock(cond,:,block)),'Color', lineColor, 'LineStyle', 'none', 'LineWidth', 2.5, 'Marker', markerType, ...
                                'MarkerFaceColor', markerFaceColor, 'MarkerEdgeColor', edgeColor, 'MarkerSize', markerSize+6, 'HandleVisibility', 'off'); hold on;

                % Add datapoints's count annotation
                zeroConstrastPoint=xBlock==0;
                if sum(zeroConstrastPoint)>0
                    nTrials=[20*ones(1,numel(zeroConstrastPoint)) 20*ones(1,numel(~zeroConstrastPoint))];
                    if cond>1
                            nTrials=[40*zeroConstrastPoint + 20*~zeroConstrastPoint];
                    end
                elseif sum(zeroConstrastPoint)>0 && plotAverageFlag==1
                    nTrials=[nBlocks*ones(1,numel(zeroConstrastPoint)) nBlocks*ones(1,numel(~zeroConstrastPoint))];
                    if cond>1
                            nTrials=[40*zeroConstrastPoint + 20*~zeroConstrastPoint];
                    end
                elseif sum(zeroConstrastPoint)==0
                     nTrials=[20*ones(1,numel(~zeroConstrastPoint))];
                elseif sum(zeroConstrastPoint)==0 && plotAverageFlag==1
                     nTrials=[nBlocks*ones(1,numel(~zeroConstrastPoint))];
                end
                %annotateDataPoints(mdl.xBlock(cond,:,block), mdl.yBlock(cond,:,block), nTrials, markerFaceColor); hold on;
            end
        end
        
        % Customize the starting position and spacing
        startPos = [30, 30]; % Starting position in data coordinates
        xSpacing =16; % Horizontal spacing between columns
        ySpacing = 6; % Vertical spacing between rows

        % Call the function to create the table
        subplot(1, 2, 1); ax1=gca;
        %createCustomTable2(ax1, modelTypeStr, mdl.headers, mdl.fittedParams(block,:), startPos, xSpacing, ySpacing);
        
        if plotAverageFlag==1
            [nConditions, ~, nBlocks] = size(behavioralData.gaborContrasts(:, :, clusterBlocks));
            nBlocks=1;
            block=1;
        end
        upFontSize(21, .01);
        %Saving
        %savePDF(savefilename, monkeyName, 1, block, nBlocks)
        % Png/SVG
        monkey=datastruct(analysisBlockID(clusterBlocks(block))).monkey;
        date= datastruct(analysisBlockID(clusterBlocks(block))).date;
        run=datastruct(analysisBlockID(clusterBlocks(block))).run;
        if ispc
          mainPath='Y:/';
        elseif contains(getenv('HOSTNAME'),'psy.utexas.edu')
          mainPath='/eslab/data/';
        end
        figFilename=[mainPath monkey '\Meta\psychometrics\L-chamber\weibullfreeAll\C' num2str(cluster) 'M28D' date 'R' run 'projection'];
        set(findall(gcf, '-property', 'FontName'), 'FontName', 'SansSerif');                
        set(gcf, 'Renderer', 'painters'); % Use painters for vector graphics
        print(gcf, [figFilename '.png'], '-dpng', '-r600'); % High-res PNG
        savefig(gcf, [figFilename '.fig']);           % FIG
        print(gcf, [figFilename '.svg'], '-dsvg');        % SVG
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
    plot(mdl.xFitted(cond,:,block), mdl.yFitted(cond,:,block), 'Color', lineColor, 'LineWidth', 3, 'HandleVisibility', 'on'); hold on;
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

function createCustomTable(ax, headers, fitParams, startPos, xSpacing, ySpacing)
    % Create a custom table on axes using text objects with subscripts.
    %
    % Parameters:
    %   ax - The axes handle where the table will be drawn.
    %   headers - A cell array of header names with subscripts.
    %   fitParams - A 1xN numeric array containing the values.
    %   startPos - Starting position of the table ([x, y]) in data coordinates.
    %   xSpacing - Horizontal spacing between columns.
    %   ySpacing - Vertical spacing between rows.

    % Ensure fitParams and headers are the same length
    if length(headers) ~= length(fitParams)
        error('Headers and fitParams must have the same length.');
    end
    fitParams=round(fitParams,1);
    % Plot headers
    for i = 1:length(headers)
        text(startPos(1) + (i - 1) * xSpacing, startPos(2), ...
             headers{i}, ...
             'Interpreter', 'tex', ...
             'HorizontalAlignment', 'center', ...
             'FontWeight', 'bold', ...
             'Parent', ax); % Use the specified axes
    end

    % Plot values
    for i = 1:length(fitParams)
        text(startPos(1) + (i - 1) * xSpacing, startPos(2) - ySpacing, ...
             sprintf('%.1f', fitParams(i)), ...
             'HorizontalAlignment', 'center', ...
             'Parent', ax); % Use the specified axes
    end
end

function createCustomTable2(ax, modelTypeStr, headers, fitParams, startPos, xSpacing, ySpacing)
    % Create a custom table on axes using text objects with subscripts.
    %
    % Parameters:
    %   ax - The axes handle where the table will be drawn.
    %   headers - A cell array of header names with subscripts.
    %   fitParams - A 1xN numeric array containing the values.
    %   startPos - Starting position of the table ([x, y]) in data coordinates.
    %   xSpacing - Horizontal spacing between columns.
    %   ySpacing - Vertical spacing between rows.
    
    % Exclude AICc for now
    % Strings to remove
    stringsToRemove = {'AICc^{total}', 'AUC^{con-incon}'};

    % Find the indices of the strings to remove
    idxToRemove = find(ismember(headers, stringsToRemove));
    headers(idxToRemove)=[];
    fitParams(idxToRemove)=[];

    % Resort for weibullfreeAll
    if strcmp(modelTypeStr,'weibullfreeAll')
        % Re-sort values to have rows = bl; con; incon.
        sortIdx=[1:4 13, 9:12 14, 5:8 15];
        % Sort
        fitParams=fitParams(sortIdx);
        headers=headers(sortIdx);
        
        % Indices for parameters
        idxA=1:5:numel(headers);
        idxB=idxA+1;
        idxAlpha=idxA+2;
        idxBeta=idxA+3;
        idxAUC=idxA+4;
        % Convert A & B to 100%
        fitParams([idxA,idxB]) = fitParams([idxA,idxB]) * 100;

        % A: Convert A to 100-A, calculate delta from baseline
        fitParams(idxA) = 100-fitParams(idxA);
        fitParams(idxA(2:3)) = fitParams(idxA(2:3)) - fitParams(idxA(1));
        % B: Flip, calculate delta from baseline
        fitParams(idxB(3))=fitParams(idxB(3))-50;
        fitParams(idxB(2))=-fitParams(idxB(3));
        % Alpha: calculate delta from baseline
        fitParams(idxAlpha(2:3)) = fitParams(idxAlpha(2:3)) - fitParams(idxAlpha(1));
        % Beta: calculate delta from baseline
        fitParams(idxBeta(2:3)) = fitParams(idxBeta(2:3)) - fitParams(idxBeta(1));
        % AUC: calculate delta from baseline
        fitParams(idxAUC(2:3)) = fitParams(idxAUC(2:3)) - fitParams(idxAUC(1));
        1;
    elseif strcmp(modelTypeStr,'weibull-beta')
        sortIdx=[1:4, 9:12, 5:8];
        % C50
        fitParams(10)=1-fitParams(6)-.5;
        fitParams(6)=fitParams(6)-.5;
        % Sort
        
        fitParams([2 6 10])=fitParams([2 6 10])*100;
        fitParams([1 5 9])=fitParams([1 5 9])*100;
        
        %A
        fitParams([5 9]) = fitParams([5 9]) - fitParams(1);
        %B
        fitParams([7 11]) = fitParams([7 11]) - fitParams(3);
        %Beta
        fitParams([8 12]) = fitParams([12 12]) - fitParams(4);
    end

    %fitParams = round(fitParams, 3, 'significant');

    % Determine the number of items per row
    numRows = 3;
    itemsPerRow = floor(length(headers) / numRows);

    % Plot headers and values
    for row = 1:numRows
        % Extract headers and fitParams for the current row
        startIdx = (row - 1) * itemsPerRow + 1;
        endIdx = min(row * itemsPerRow, length(headers));
        currentHeaders = headers(startIdx:endIdx);
        currentParams = fitParams(startIdx:endIdx);

        % Plot headers for the current row
        for i = 1:length(currentHeaders)
            text(startPos(1) + (i - 1) * xSpacing, ...
                 startPos(2) - (row - 1) * ySpacing-(5 * (row-1)), ...
                 currentHeaders{i}, ...
                 'Interpreter', 'tex', ...
                 'HorizontalAlignment', 'center', ...
                 'FontWeight', 'bold', ...
                 'Parent', ax); % Use the specified axes
        end
        
        for i = 1:length(currentParams)
            value = currentParams(i);
            if 1==1% Check if the value is an integer
                %formattedValue = sprintf('%.0f', value); % Format integers with .0

                % Force 2 significant figures while ensuring consistent formatting
                formattedValue = sprintf('%.2f', round(value, 2)); 
            else
            end
            text(startPos(1) + (i - 1) * xSpacing, ...
                 startPos(2) - row * ySpacing - (5 * (row - 1)), ...
                 formattedValue, ...
                 'HorizontalAlignment', 'center', ...
                 'Parent', ax); % Use the specified axes
        end


        
        %{
        for i = 1:length(currentParams)
            text(startPos(1) + (i - 1) * xSpacing, ...
                 startPos(2) - row * ySpacing-(5 * (row-1)), ...
                 sprintf('%.2f', currentParams(i)), ...
                 'HorizontalAlignment', 'center', ...
                 'Parent', ax); % Use the specified axes
        end
        %}
    end
end
function y = replaceNanSections(y)
    % Replace NaN sections in an array for a Weibull curve.
    % If NaN sections are flanked by equal values, replace them with those values.
    % For NaN sections at the end of the curve, use the last valid plateau value
    % if the last two valid values are similar within 0.01.
    y=real(y);
    % Ensure input is a row vector for easier indexing
    if iscolumn(y)
        y = y';
    end

    % Identify NaN sections
    isnanIndices = isnan(y);
    nanStarts = find(diff([0 isnanIndices]) == 1); % Start indices of NaN sections
    nanEnds = find(diff([isnanIndices 0]) == -1);  % End indices of NaN sections

    % Iterate through each NaN section
    for i = 1:length(nanStarts)
        startIdx = nanStarts(i);
        endIdx = nanEnds(i);

        % Case 1: Regular NaN sections with left and right flanks
        if startIdx > 1 && endIdx < length(y)
            leftFlank = y(startIdx - 1);
            rightFlank = y(endIdx + 1);

            % Replace NaNs if left and right flanks are equal (or nearly equal)
            if abs(leftFlank - rightFlank) <= 0.01
                y(startIdx:endIdx) = leftFlank;
            end
        % Case 2: NaN sections at the end of the curve
        elseif startIdx > 1 && endIdx == length(y)
            % Check the last two valid values to the left of the NaN section
            leftFlank1 = y(startIdx - 1);
            leftFlank2 = y(startIdx - 2);

            % Replace NaNs if the last two values are similar
            if abs(leftFlank1 - leftFlank2) <= 0.01
                y(startIdx:endIdx) = leftFlank1;
            end
        end
    end
end


