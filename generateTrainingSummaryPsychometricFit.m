formattedDates = getDates(-0);
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
            if isfield(TS.Header.Conditions,'StimCon') && numel(unique(TS.Header.Conditions.StimCon))>=4
                StimCon=TS.Header.Conditions.StimCon;
                uniqueStimCon=unique(StimCon);
                nUniqueStimCon=numel(uniqueStimCon);
            else
                continue
            end
            

            %% Add stimulus page
            if blockNo==1
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
            
            %fit here
            figure()
            horizontalIdx=find(dataXOrt==0);
            dataX(horizontalIdx)=dataX(horizontalIdx)*-1;
            dataY(horizontalIdx)=100-dataY(horizontalIdx);
            [threshold,betaFit,muFit,sigmaFit,param,x,y] = fitPsyDualCDFv3(dataX',dataY',sessionID{sessID},[0 0 0],'o',gca)
            %cosmetics
            xlim([-100 100]);
            ylim([0 100])
            xlabel('Gabor contrast (%)'); xticks([-100:25:100]); axis square
            ylabel('Vertical reports (%)'); yticks([0:10:100]);
            upFontSize(32,.01)
            legend('Training', 'Location','southeast')
            upFontSize
            
            % Add params table
            header={'Threshold','Beta','n','C50'};
            dat=round([threshold, betaFit*100 sigmaFit muFit],1);

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
                'Units','Normalized','Position',[.75 -.45 1 1]);
            
            dataXOrt=[];
            dataX=[];
            dataY=[];
        end

        %% Save
        if sessID<numel(sessionID)
            %savePDF(savefilename, 'Pepper', saveFlag, counter, counter+1); 
            counter=counter+1;
        elseif sessID==numel(sessionID)
            %savePDF(savefilename, 'Pepper', saveFlag, counter, counter); 
        end 
    end
end
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
