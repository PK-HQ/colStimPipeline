formattedDates = getDates(-1);
set(0,'DefaultFigureWindowStyle','docked')
plotTrainingSummary('Y:\Pepper\Meta\training\', '32', formattedDates, 1)
function plotTrainingSummary(datapath, monkeyID, sessionID, saveFlag)
%% Extract conditions
condDifficulty={'(Easy)', '(Hard)'};
colors={'b','r'};
mkr={'v','^'};
counter=1;
nRows=2;nCols=2;
gap=.15;marginV=.1;marginH=.075;
nBlockPlot=10;

correctTrialsTotal=nan(nBlockPlot,4,numel(sessionID));
titlePage=0;
for sessID=1:numel(sessionID)
    %% Load
    [fullFileNames]=findREGEX([datapath 'M' monkeyID 'D' sessionID{sessID} 'R*TS.mat']);
    savefilename=['training/M' monkeyID 'summary'];

    if isempty(fullFileNames)
        continue
    else
        titlePage=titlePage+1;
        nBlocks=numel(fullFileNames);
        fprintf('%s found, loading %.0f blocks\n',sessionID{sessID}, nBlocks)
        
        % Init
        dataY=[];
        dataYCounts=[];
        dataX=[];
        dataXOrt=[];
        % horizontal as negative
        horizontalIdx=[];
        dataX=[];
            
        for blockNo=1:nBlocks
            load(fullFileNames{blockNo})

            disp(TS.Header.ConditionTable)
            
            %% Get relevant parameters for OD task (Size, SF, Ort, StimCon)
            % Size
            Sz=TS.Header.Conditions.GaborSize;
            uniqueSz=unique(Sz);
            nSz=numel(uniqueSz);

            % Spatial frequency
            SF=TS.Header.Conditions.GaborSF;
            uniqueSF=unique(SF);
            nSF=numel(uniqueSF);

            % Orientations
            Ort=TS.Header.Conditions.GaborOrt;
            uniqueOrt=unique(Ort);
            nUniqueOrt=numel(uniqueOrt);
            
            % Stim contrast
            if isfield(TS.Header.Conditions,'StimCon')
                StimCon=TS.Header.Conditions.StimCon;
                uniqueStimCon=unique(StimCon);
                nUniqueStimCon=numel(uniqueStimCon);
            else
                continue
            end
            

            %% Add stimulus page
            if titlePage==1
                figure('Name',['M' monkeyID 'S' sessionID{sessID} '-CondTable']); set(gcf,'Color','w')
                formattedStr={};
                % Add protocol text
                text1=TS.Header.ConditionTable([1:5]);
                    annotation('TextBox',[0.05,0.6,0.9,0.1], ...
                    'String',text1,'HorizontalAlignment','Center', ...
                    'FontWeight','bold', ...
                    'FontSize',12, ...
                    'LineStyle','None','Interpreter','none');
                
                % Add header
                text1=TS.Header.ConditionTable(7);
                    annotation('TextBox',[0.05,0.44,0.9,0.1], ...
                    'String',text1,'HorizontalAlignment','Center', ...
                    'FontWeight','bold', ...
                    'FontSize',12, ...
                    'LineStyle','None','Interpreter','none');

                % Add parameters
                condTableSize=size(TS.Header.ConditionTable,1);
                for i=8:condTableSize-1
                    inputStr = TS.Header.ConditionTable(i);
                    numSpaces = 4; % Number of additional spaces to add
                    formattedStr(i-7,1)= addSpaces(inputStr, numSpaces);
                end

                annotation('TextBox',[0.038,0.4,0.9,0.1], ...
                'String',formattedStr,'HorizontalAlignment','Center', ...
                'FontWeight','normal', ...
                'FontSize',14, ...
                'LineStyle','None');
                savePDF(savefilename, 'Pepper', saveFlag, 0, counter+1); offwarning
            end
            
            plotIDs=[1 2];
            % Determine manipulated parameter
            if isfield(TS.Header.Conditions,'StimCon')
                paramVal=StimCon;
                uniqueParamVal=unique(paramVal);
                titleCond={'Gabor stimuli, varied contrast'};
            else
                disp('Error')
            end

            %% Extract result
            trialCorrect=TS.Header.Outcomes.CountCondSuccess;
            trialComplete=TS.Header.Outcomes.CountCondTotalValid;
            trialCompletePrct=trialComplete ./sum(trialComplete);
            trialCorrectPrct=trialCorrect * 100 ./ trialComplete;
            trialTotal=TS.Header.Outcomes.CountBlockTotal;
            %Collect X and Y plotting data
            dataY=[dataY;trialCorrectPrct];
            dataYCounts=[dataYCounts;trialCorrect];
            dataX=[dataX;StimCon'];
            dataXOrt=[dataXOrt;Ort'];
        end
        % horizontal as negative
        horizontalIdx=find(dataXOrt==0);
        dataX(horizontalIdx)=dataX(horizontalIdx)*-1;
        % Plot
        figure()
        % Shading, guide lines
        shadePlot([30 100], [80 100], 'g', .5);
        shadePlot([-30 -100], [0 20], 'g', .5);
        yline(50,'LineStyle','--','color',[.3 .3 .3],'HandleVisibility','off');
        xline(0,'LineStyle','--','color',[.3 .3 .3],'HandleVisibility','off');
        % Scatter
        scatter(dataX(dataX<0),100-dataY(dataX<0),500,'Marker', 'v', ...
        'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'k','LineWidth',2); hold on
        scatter(dataX(dataX>0),dataY(dataX>0),500,'Marker', '^', ...
        'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'k','LineWidth',2); hold on
    
        hText=text(dataX, dataY, num2str(dataYCounts), 'Color', 'w', 'HorizontalAlignment', 'center',...
            'VerticalAlignment', 'middle', 'Margin', 100,  'FontSize',4);hold on
        hText=text(dataX(dataX<0), 100-dataY(dataX<0), num2str(dataYCounts(dataX<0)), 'Color', 'w', 'HorizontalAlignment', 'center',...
            'VerticalAlignment', 'middle', 'Margin', 100,  'FontSize',4);hold on

        %Axes
        xlim([-100 100]);
        ylim([0 100])
        suplabel(['M' monkeyID 'S' sessionID{sessID}],'t',[.1 .1 .85 .85]);
        upFontSize
        ax=gca;
        yyaxis left; ylabel('Correct H (%)'); yticks([0:10:100]);ax.YTickLabel = flipud(ax.YTickLabel);
        ax=gca;
        yyaxis right; ylabel('Correct V (%)'); ylim([0 100]); yticks([0:10:100]); set(gca,'ycolor','k') 
        xlabel('Gabor contrast (%)'); xticks([-100:25:100]); axis square
        axis square
        upFontSize(32,.01)
        legend('H-Gabor', 'V-Gabor', 'Location','southeast')
        
        %% Save
        if sessID<numel(sessionID)
            savePDF(savefilename, 'Pepper', saveFlag, counter, counter+1); offwarning
            counter=counter+1;
        elseif sessID==numel(sessionID)
            savePDF(savefilename, 'Pepper', saveFlag, counter, counter); offwarning
        end 
    end
end

% Init figure
figure('Name','Session average')
for i=1:4
    scatter(1:numel(sessionID),squeeze(nanmean(correctTrialsPrct(:,i,:))))
end
ylim([0 100])
xlim([0 10])

end

function output = addSpaces(inputStr, numSpaces)
    % Function to add spaces in a string after each character followed by a space

    % Define the pattern: any character followed by a space
    pattern = '(?<=\S) ';

    % Define the replacement: the same character followed by extra spaces
    replacement = [repmat(' ', 1, numSpaces) ' '];

    % Use regexprep to replace the pattern with the replacement
    output = regexprep(inputStr, pattern, replacement);
end
