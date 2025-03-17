formattedDates = getDates(-13);
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

        for blockNo=1:nBlocks
            load(fullFileNames{blockNo})

            
            disp(TS.Header.ConditionTable)
            
            %% Get relevant parameters for OD task
            % Spatial frequency
            SF=TS.Header.Conditions.GaborSF;
            uniqueSF=unique(SF);
            nSF=numel(uniqueSF);

            % Incorrect target contrast
            IncorTarCon=TS.Header.Conditions.IncorTarCon;
            uniqueIncorTarCon=unique(IncorTarCon);
            nIncorTarCon=numel(uniqueIncorTarCon);

            % Orientations
            Ort=TS.Header.Conditions.GaborOrt;
            uniqueOrt=unique(Ort);
            nUniqueOrt=numel(uniqueOrt);
            
            if isfield(TS.Header.Conditions,'BarContrast')
                barContrast=TS.Header.Conditions.BarContrast;
                uniqueBarContrast=unique(barContrast);
                nUniqueBar=numel(uniqueBarContrast);
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
            if blockNo==1
                % Init figure
                figure('Name',['M' monkeyID 'S' sessionID{sessID} '-Psychometric'])
                [hAx,~]=tight_subplot(nRows,nCols,[gap gap], [marginV marginV], [marginH marginH]);
            end
            
            plotIDs=[1 2];
            % Determine manipulated parameter
            if nIncorTarCon>1 % Easy = lower incorrect target contrast vs full contrast targets
                % target contrast < 100% == easy
                paramVal=IncorTarCon;
                uniqueParamVal=unique(paramVal);
                titleCond={'Dual contrast targets', '100% contrast targets'};
            elseif nSF>1 % Easy = bar+gabor vs gabor
                % spatial frequency < 1 == easy
                paramVal=SF;
                uniqueParamVal=unique(paramVal);
                titleCond={'Gabor+bar stimuli', 'Gabor stimuli'};
            elseif nUniqueBar>1 || isfield(TS.Header.Conditions,'BarContrast')
                % spatial frequency < 1 == easy
                paramVal=barContrast;
                uniqueParamVal=sort(unique(paramVal),'descend');
                if numel(uniqueParamVal)==1
                    uniqueParamVal=repmat(uniqueParamVal,1,2);
                end
                titleCond={'Gabor+bar stimuli', 'Gabor+bar stimuli'};            
                disp([fprintf('Block %.0f, bar contrast',blockNo) '(' num2str(unique(uniqueParamVal)) ')\n'])
            elseif uniqueSF==1 && nIncorTarCon==1 % Hard = bar+gabor only
                % spatial frequency < 1 == easy
                paramVal=[SF];
                uniqueParamVal=paramVal;
                titleCond={'Gabor+bar stimuli', 'Gabor stimuli'};
                plotIDs=2;
            elseif all(SF < 1) && nUniqueBar<2
                continue
            end

            %% Extract result
            trialCorrect=TS.Header.Outcomes.CountCondSuccess;
            trialComplete=TS.Header.Outcomes.CountCondTotalValid;
            trialCompletePrct=trialComplete ./sum(trialComplete);
            trialCorrectPrct=trialCorrect * 100 ./ trialComplete;
            trialTotal=TS.Header.Outcomes.CountBlockTotal;
            
            %% Plot counts
            for plotID=plotIDs % for H-stim, V-stim or combined
                if plotID==2 && numel(unique(uniqueParamVal))==1 
                    continue
                end
                axes(hAx(plotID)); hold on;
                for ortNo=1:nUniqueOrt
                    condNo=find(paramVal==uniqueParamVal(plotID));
                    scatter(blockNo, sum(trialCompletePrct(condNo)), 500, 'k', 'filled', 'o', 'LineWidth', 2.5, 'Marker', 'o'); hold on;
                end
            end
            
            %% Plot percent correct
            for plotID=plotIDs % for H-stim & V-stim combined
                if plotID==2 && numel(unique(uniqueParamVal))==1 
                    continue
                end
                axes(hAx(plotID+2)); hold on;
                condNo=find(paramVal==uniqueParamVal(plotID));
                scatter(blockNo, mean(trialCorrectPrct(condNo)), 500, 'k',  'filled', 'o', 'LineWidth', 2.5, 'Marker', 'o'); hold on;
                text(blockNo+.4, mean(trialCorrectPrct(condNo)), num2str(sum(trialCorrect(condNo))), 'HorizontalAlignment', 'center',...
                    'VerticalAlignment', 'middle', 'Margin', 100,  'FontSize',4)
                % Bar con
                if exist('barContrast')
                    text(blockNo,53.5, [num2str(unique(barContrast(condNo))) '%'], 'HorizontalAlignment', 'center',...
                        'VerticalAlignment', 'middle', 'Margin', 100,  'FontSize',10)
                end
                for ortNo=1:nUniqueOrt % for H-stim & V-stim individually
                    condNo=intersect(find(paramVal==uniqueParamVal(plotID)), find(Ort==uniqueOrt(ortNo)));
                    scatter(blockNo, mean(trialCorrectPrct(condNo)), 500, colors{ortNo}, 'filled', 'o', 'LineWidth', 2.5, 'Marker', mkr{ortNo}); hold on;
                    hText=text(blockNo, mean(trialCorrectPrct(condNo)), num2str(sum(trialCorrect(condNo))), 'Color', 'w', 'HorizontalAlignment', 'center',...
                        'VerticalAlignment', 'middle', 'Margin', 100,  'FontSize',4);
                    % Ensure the text appears on top
                    uistack(hText, 'top');
                    set(gca, 'Layer', 'top');
                    
                    % Save correct counts
                    %correctTrialsPrct(blockNo, ortNo+2*(plotID-1), sessID)=trialCorrectPrct(condNo);
                end
            end
        end

         for plotID=1:2
            axes(hAx(plotID)); hold on;
            title(['Proportion per condition ' condDifficulty{plotID}]); hold on;
            xlim([0 10]); ylim([0 1]);
            addSkippedTicks(0,10,1,'x')
            addSkippedTicks(0,1,.125,'y')
            if plotID==1
                ylabel('Proportion')
                legend({'Total trials'})
            end
            upFontSize
         end
        
        for plotID=1:2
            axes(hAx(plotID+2)); hold on;
            title([titleCond{plotID} ' ' condDifficulty{plotID}]); hold on;
            xlim([0 nBlockPlot]); ylim([50 100]);
            addSkippedTicks(0,nBlockPlot,1,'x')
            addSkippedTicks(50,100,5,'y')
            yChance=yline(50,'--'); yChance.Color=.5*ones(1,3); yChance.HandleVisibility='off';  yChance.LineWidth=1.5; yChance.Label='Chance';
            yCriterion=yline(80,'-'); yCriterion.Color=[0 99 50]/255; yCriterion.Alpha=1; yCriterion.LineWidth=1.5;yCriterion.Label='Criterion';yCriterion.LabelVerticalAlignment='bottom';
            if plotID==1
                xlabel('Block no.')
                ylabel('Correct %')
                h=legend({'Average','Horizontal', 'Vertical'},'Location','southwest');
            end
            upFontSize
        end
        suplabel(['M' monkeyID 'S' sessionID{sessID}],'t',[.1 .1 .85 .85]);
        upFontSize

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
