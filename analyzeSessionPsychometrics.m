function analyzeSessionPsychometrics(behavioralData, bitmapData, datastruct, analysisBlockID, saveFlag)
% Naka-Rushton params (Rmax, Exponent, C50, Beta)
initParams.min = [10, .1, 10, 20];
initParams.max = [80, 4, 50, 80];

lineColor={[0 0 0],[0.9294, 0.1098, 0.1373]*1.05,[0, 0.0941, 0.6627]*1.25,[0.4471, 0.0353, 0.7176]*1.1};%{'k','#ED1C23','#0018A9','#7209B7'};
markerType={'o','^','v','diamond'};
lineOrder=[2 1 3 4];xpos=0;ypos=-.20;


%% Plot summary of all blocks
columnRanges=[1 40; 1 5; 6 19; 20 40];

for columnRange=1:size(columnRanges,1)
    % Get blocks within column range
    [selectedBlocks, selectedBlockIDs]=powerfilter(bitmapData, analysisBlockID, 'fixed',columnRanges(columnRange,:));
    % Setup subplot axes
    figure('Name',sprintf('Psych. summary (%.0f-%.0f)', columnRanges(columnRange,1),columnRanges(columnRange,2)))
    nRows=1;nCols=3;
    gap=.05;marginV=.01;marginH=.075;
    [hAx,~]=tight_subplot(nRows,nCols,[gap gap], [marginV+.27 marginV+.13], [marginH marginH]);
    for plotID=1:3 % for H-stim, V-stim or combined
        axes(hAx(plotID));
        for lineNo=1:4 % plot BL/H/V/Diff-opto
            if lineNo==1
                %line([0 0],[0 100],'Color',.6*[1 1 1],'LineStyle','--','LineWidth',.25,'HandleVisibility','off'); hold on
                line([0 100],[50 50],'Color',.6*[1 1 1], 'LineStyle','--','LineWidth',.25,'HandleVisibility','off'); hold on
            end
            x=squeeze(behavioralData.gaborContrast(lineNo,plotID,:,selectedBlocks));
            y=squeeze(behavioralData.percentageCorrect(lineNo,plotID,:,selectedBlocks));

            %SEM plot
            [binnedX, meanY, binnedY]=plotSEM(x,y,lineColor{lineNo},markerType{lineNo}); hold on;

            if lineNo<=3
                [fitParams] = fitNakaRushtonMLE(binnedX, binnedY, initParams);
                %fitParams = fitNakaRushtonRobustaBisquare(binnedX, meanY, initParam);
                % NR fit plot
                behavioralData.fit(lineNo,plotID)=fitParams; % saving params
                %offwarning
                modelFunc = @(b, x) (b(1) .* (x.^b(2)) ./ (b(3).^b(2) + x.^b(2))) + b(4); %remove epsilon
                xplot= 0:1:100;
                fitParams.fittedMeanY=modelFunc([fitParams.rmax, fitParams.exponent, fitParams.c50, fitParams.beta],xplot);
                plotFittedLine(xplot, fitParams.fittedMeanY, lineColor{lineNo}); hold on
            elseif lineNo==4
                [fitParams] = fitNakaRushtonMLE(binnedX, 100-binnedY, initParams);
                %fitParams = fitNakaRushtonRobustaBisquare(binnedX, 100-meanY, initParam);
                % NR fit plot
                behavioralData.fit(lineNo,plotID)=fitParams; % saving params
                %offwarning
                modelFunc = @(b, x) (b(1) .* (x.^b(2)) ./ (b(3).^b(2) + x.^b(2))) + b(4); %remove epsilon
                xplot= 0:1:100;
                fitParams.fittedMeanY=modelFunc([fitParams.rmax, fitParams.exponent, fitParams.c50, fitParams.beta],xplot);
                fitParams.fittedMeanY=100-fitParams.fittedMeanY;
                plotFittedLine(xplot, fitParams.fittedMeanY, lineColor{lineNo}); hold on
            end

            if  lineNo==1 % Add header
                % State fit parameters
                headerText = sprintf(...
                 ['{\\bfParameters}        {\\bf\\beta}        {\\bfR_{max}}        {\\bfC_{50}}         {\\bfn}          {\\bfR_{robust}^{2}}']);
                 text(xpos,ypos,headerText,'Units','normalized', 'VerticalAlignment','top', 'HorizontalAlignment','left','Color','k')
            end

                % Add param text
                newlineStr=repmat('\n',1,lineOrder(lineNo));
                paramsText=sprintf(...
                [newlineStr '                           %.0f        %.0f            %.0f          %02.f          %0.2f'], ...
                fitParams.beta, fitParams.rmax, fitParams.c50, fitParams.exponent, fitParams.R2);
                text(xpos,ypos-.02,paramsText,'Units','normalized', 'VerticalAlignment','top', 'HorizontalAlignment','left','Color',lineColor{lineNo});
                upFontSize(24,0.01)
        end
        upFontSize(24,0.01)
        addSkippedTicks(0,100,5,'y')
        addSkippedTicks(0,100,12.5,'x')
    end
    % Subplot #1
    axes(hAx(1)); title('Vertical stimuli','FontWeight','normal'); ylabel('Correct (%)'); xlabel('Absolute gabor contrast (%)'); 
    % legend order
    legend({'Difference','Incongruent','Congruent','Baseline'},'Location', 'east', 'FontSize',18);
    h = get(gca,'Children'); set(gca,'Children',[h(1:end-4); h([end end-1 end-2 end-3])]);    
    % Subplot #2    
    axes(hAx(2)); title('Horizontal stimuli','FontWeight','normal');
    % Subplot #3
    axes(hAx(3)); title('Combined','FontWeight','normal');
    
    suplabel(sprintf('Session summary (n_{col}=%.0f - %.0f, n_{blocks}=%.0f)',columnRanges(columnRange,1), columnRanges(columnRange,2),...
        numel(selectedBlocks)),'t',[.1 .1 .82 .84])
    upFontSize(24,0.01)
    if columnRange==1 && saveFlag==1
        export_fig('Y:/Chip/Meta/biasingSeries/biasingSeriesSummaryRclean.pdf','-pdf','-nocrop');
    elseif columnRange>1 && saveFlag==1
        export_fig('Y:/Chip/Meta/biasingSeries/biasingSeriesSummaryRclean.pdf','-pdf','-append','-nocrop');
    end
end

%dupe to meetings
copyfile('Y:/Chip/Meta/biasingSeries/biasingSeriesSummaryRclean.pdf', 'Y:/users/PK/Eyal/meetings/biasingSeries/biasingSeriesSummaryRclean.pdf')

%% Plot individual blocks
% Filter by power
[selectedBlocks, selectedBlockIDs]=powerfilter(bitmapData, analysisBlockID, 'fixed');
% define session IDs, session meta info for filename
currentSessID=selectedBlockIDs;
dsCurrentSess=datastruct(currentSessID); %fixed, to get ort map projected onto
dsReferenceSess=datastruct(vertcat(dsCurrentSess.referenceBlockNo));
alignmentSession=datastruct(vertcat(dsCurrentSess.alignmentBlockNo)).date;

for blockNo=1:numel(selectedBlocks)
    selectedBlock=selectedBlocks(blockNo);
    % Get savefilename
     filenameStructCurrent=generateFilenames(dsCurrentSess(blockNo));
    % Setup subplot axes
    figure('Name','Psych. summary')
    nRows=1;nCols=3;
    gap=.05;marginV=.01;marginH=.075;
    [hAx,~]=tight_subplot(nRows,nCols,[gap gap], [marginV+.27 marginV+.13], [marginH marginH]);
    for plotID=1:3 % for H-stim, V-stim or combined
        axes(hAx(plotID));
        for lineNo=1:4 % plot BL/H/V/Diff-opto
            if lineNo==1
                %line([0 0],[0 100],'Color',.6*[1 1 1],'LineStyle','--','LineWidth',.25,'HandleVisibility','off'); hold on
                line([0 100],[50 50],'Color',.6*[1 1 1], 'LineStyle','--','LineWidth',.25,'HandleVisibility','off'); hold on
            end
            x=squeeze(behavioralData.gaborContrast(lineNo,plotID,:,selectedBlocks(blockNo)));
            y=squeeze(behavioralData.percentageCorrect(lineNo,plotID,:,selectedBlocks(blockNo)));

            %SEM plot
            [binnedX, meanY]=plotSEM(x,y,lineColor{lineNo},markerType{lineNo}); hold on;

            if lineNo<=3
                fitParams = fitNakaRushtonRobustaBisquare(binnedX, meanY, initParams);
                % NR fit plot
                behavioralData.fit(lineNo,plotID)=fitParams; % saving params
                %offwarning
                modelFunc = @(b, x) (b(1) .* (x.^b(2)) ./ (b(3).^b(2) + x.^b(2))) + b(4); %remove epsilon
                xplot= 0:1:100;
                fitParams.fittedMeanY=modelFunc([fitParams.rmax, fitParams.exponent, fitParams.c50, fitParams.beta],xplot);
                plotFittedLine(xplot, fitParams.fittedMeanY, lineColor{lineNo}); hold on
            elseif lineNo==4
                %{
                fitParams = fitNakaRushtonRobustIRLS(binnedX, 100-meanY, initParam);
                % NR fit plot
                behavioralData.fit(lineNo,plotID)=fitParams; % saving params
                %offwarning
                modelFunc = @(b, x) (b(1) .* (x.^b(2)) ./ (b(3).^b(2) + x.^b(2))) + b(4); %remove epsilon
                xplot= 0:1:100;
                fitParams.fittedMeanY=modelFunc([fitParams.rmax, fitParams.exponent, fitParams.c50, fitParams.beta],xplot);
                fitParams.fittedMeanY=100-fitParams.fittedMeanY;
                plotFittedLine(xplot, fitParams.fittedMeanY, lineColor{lineNo}); hold on
                %}
            end

            if  lineNo==1 % Add header
                % State fit parameters
                headerText = sprintf(...
                 ['{\\bfParameters}        {\\bf\\beta}        {\\bfR_{max}}        {\\bfC_{50}}         {\\bfn}          {\\bfR_{robust}^{2}}']);
                 text(xpos,ypos,headerText,'Units','normalized', 'VerticalAlignment','top', 'HorizontalAlignment','left','Color','k')
            end
            if lineNo<=3
                % Add param text
                newlineStr=repmat('\n',1,lineOrder(lineNo));
                paramsText=sprintf(...
                [newlineStr '                           %.0f        %.0f            %.0f          %02.f          %0.2f'], ...
                fitParams.beta, fitParams.rmax, fitParams.c50, fitParams.exponent, fitParams.R2);
                text(xpos,ypos-.02,paramsText,'Units','normalized', 'VerticalAlignment','top', 'HorizontalAlignment','left','Color',lineColor{lineNo});
                upFontSize(24,0.01)
            end
        end
        upFontSize(24,0.01)
        addSkippedTicks(0,100,5,'y')
        addSkippedTicks(0,100,12.5,'x')
    end
    % Subplot #1
    axes(hAx(1)); title('Vertical stimuli','FontWeight','normal'); ylabel('Correct (%)'); xlabel('Absolute gabor contrast (%)'); 
    % legend order
    legend({'Incongruent','Congruent','Baseline'},'Location', 'east', 'FontSize',18);
    h = get(gca,'Children'); set(gca,'Children',h([11 9 10 1:8]));    
    % Subplot #2    
    axes(hAx(2)); title('Horizontal stimuli','FontWeight','normal');
    % Subplot #3
    axes(hAx(3)); title('Combined','FontWeight','normal');
    upFontSize(24,0.01)
    
    % Supertitle
    [~,h]=suplabel({[filenameStructCurrent.date 'R' filenameStructCurrent.run ... % Block date and number
        ' (' num2str(bitmapData.nColumns(1,selectedBlock)) ' & ' num2str(bitmapData.nColumns(2,selectedBlock)) ' column @ '... % No. of columns
        num2str(mean(bitmapData.adjustedSPD_uW(1,:,selectedBlock)),2) ' ' ... % Surface Power density mean +- stdev
        char(0177) ' ' num2str(std(bitmapData.adjustedSPD_uW(1,:,selectedBlock)),2) ' ' char(181) 'W mm^{-2})']},'t',[.1 .1 .87 .83]);
    upFontSize(24,0.01)
    
    % Save
    if blockNo==1 && saveFlag==1
        export_fig('Y:/Chip/Meta/biasingSeries/biasingSeriesR.pdf','-pdf','-nocrop');
    elseif blockNo>1 && saveFlag==1
        export_fig('Y:/Chip/Meta/biasingSeries/biasingSeriesR.pdf','-pdf','-append','-nocrop');
    end
end


% dupe to meetings
copyfile('Y:/Chip/Meta/biasingSeries/biasingSeriesR.pdf', 'Y:/users/PK/Eyal/meetings/biasingSeries/biasingSeriesR.pdf')

end