function plotPsychometricOptostimV2(dataDir,outputDir,taskName,plotFlags,verboseTitle,allSessions,unwantedSessions,mergeBothOrientations,mergeAcrossBlocks)
%% Plots psychometric curves, thresholds, biases across sessions per monkey
%analysis params
minimumConds=3;
gaborEccWanted=1.46;%1.46;%5;
gaborSizeWanted=1;%1;%1
gaborSFWanted=6;%6;%4
stimTimeWanted=300;%6;%4

%{
gaborSize<=gaborSizeWanted...
& gaborSF>=gaborSFWanted...
& abs(gaborEccentricity-gaborEccWanted)<=1
%}
fprintf('Analysis: Z=%g, SF=%g, E=%g \n',gaborSizeWanted,gaborSFWanted,gaborEccWanted)

cd(dataDir)
%% Plot raw psychometric curve per session
%get data filenames from data directory
dataFilenames = dir(strcat(dataDir,'*TS.mat'));
dataFilenames=struct2cell(dataFilenames);
dataFilenames=dataFilenames(1,:);

%define params
monkeyIDs={'28'};%{'19'};
uniqueStr='[-1.5,-1]';

%get all unique sessionID per monkey in folder
sessionIDs={};
wantedSessions=setdiff(allSessions,unwantedSessions);
count=1;
for sess=wantedSessions
    wantedSession=contains(dataFilenames,num2str(sess));
    if sum(wantedSession)>0
        sessionIDs{count}=num2str(sess);
        count=count+1;
    end
end


nMonkey=length(monkeyIDs);
nSessions=length(sessionIDs);
%subplot dims
nCol=4;
nRow=4;%ceil(nSessions/nCol)*3;
lastNsessions=min([3-1,nSessions]);
nSessionsExample=10;
%init storage containers
allThresholds=[];
allBias=[];
sessionsUsed={};
switch plotFlags(1)
    case {1}
        plotFlag='plot';
    case {0}
        plotFlag='';
end

%filename structure
%{'monkeyID'|'sessionID' |   'runID'
%    '19'   | '20210225' |'[1 2 3 4 5]'}
runIndexPos = 13;%strfind(sessionFilenames, 'R')
sessionFilenamesWanted={};
% For each monkey, get all sessions
pdfcreated=0;
plotID=1;
allDeltaOrientations=[];
allCW=[];
allCCW=[];
allnTrialsPerCond=[];
mergedSessionData=[];
sessionOutputData=[1;1;1;1];
sessionOutputDataUnsort=[];
sessionOutputData=[1;1;1;1];
for monkey=1:nMonkey
    %get monkey's data files
    monkeyFilenameIndexes=find(contains(dataFilenames,monkeyIDs{monkey}));
    monkeyFilenames=dataFilenames(monkeyFilenameIndexes);
    % For each session get all runs
    for session=1:nSessions
        sessionFilenamesWanted={};
        fprintf('Session: %s (%g/%g) \n',sessionIDs{session},session,nSessions)
        figure;hold on

        %pdf filename
        pdfFilename=[outputDir 'M' monkeyIDs{monkey} 'Z' num2str(gaborSizeWanted*100,'%05g') 'S' num2str(gaborSFWanted*100,'%05g')...
            'E' num2str(gaborEccWanted*100,'%05g') 'N' num2str(nSessions,'%03g')];

        %get session's data files
        sessionFilenameIndexes=find(contains(monkeyFilenames,sessionIDs{session}));
        sessionFilenames=monkeyFilenames(sessionFilenameIndexes);
        counter=1;

        for files=1:length(sessionFilenames)
            fName=sessionFilenames{files};
            fileRunNo=fName((runIndexPos+1):end-4);
            sessionFilenamesWanted{files}=fName;
        end
      
        
        nRunIDs=length(sessionFilenamesWanted);
        
        % Compile the CW/CCW and delta orientations across runs
        for run=1:nRunIDs
            %init
            %allDeltaOrientations=[];
            %allCW=[];
            %allCCW=[];
            %allnTrialsPerCond=[];
            %mergedSessionData=[];
            
            %load single session data from datafolder
            runData=load(sessionFilenamesWanted{run});


            
            
            
            
            if isfield(runData.TS.Header.Conditions,'GaborSize')
                %get spatial frequency wanted only (0.5 or 2 or 8), else skip
                %{
                gaborSize=runData.TS.Header.Conditions.GaborSize;%runData.TS.Header.Conditions.StimSize;
                gaborPosition=runData.TS.Header.Conditions.StimPos;
                gaborEccentricity=sqrt(gaborPosition(1)^2+gaborPosition(2)^2);
                gaborSF=runData.TS.Header.Conditions.GaborSF;
                stimTime=runData.TS.Header.Conditions.DurationTOn;
                %}
                %get CW and CCW delta orientations presented (usually symmetrical)
                %there are 4 conditions (2 orientations X 2 phases)
                switch taskName
                  case {'optostim'}
                    %xaxisMeasure=runData.TS.Header.Conditions.OrtCond;
                    %xaxisMeasureStr=xaxisMeasureStr;
                    
                    
                    
                    
                    % relevant DV
                    condData=runData.TS.Header.Conditions;
                    condIDs=1:numel(condData.TypeCond);
                    stimType=condData.TypeCond; %0 1 3
                    gaborOrt=condData.GaborOrt;
                    projectedImageConds=runData.TS.Header.Conditions.ProjImg;
                    visualControl=find(stimType==1);
                    optoControl=contains(projectedImageConds,'L000');
                    opto0=contains(projectedImageConds,'O00000');
                    opto90=contains(projectedImageConds,'O09000');
                    successData=runData.TS.Header.Outcomes.CountCondSuccess;
                    completedData=runData.TS.Header.Outcomes.CountCondTotalValid;
                    

                    
                  %% for H and V, plot con and incon for early/late/combined trials
                    [outcomeMat, outcomeTbl]=getTrialOutcomes(runData);                    

                    %
                    mkrCond={'o','^','v'};
                    
                    % colormap
                    cmap=fireice;
                    cmap=cmap([22 1+10 43-10],:);
                    nTrials=10;
                    
                     Hdata=outcomeMat(sign(outcomeMat(:,2))==1,:);
                     Vdata=outcomeMat(sign(outcomeMat(:,2))==-1,:);
                     Vdata(:,2)=-Vdata(:,2); %code as +
                     % --- groups for plotting ---
                     titles={'Horizontal gabor','Vertical gabor'};
                     datasets={Hdata,Vdata};
                     
                    % fig
                    nSubplotRow=2;
                    nSubplotCol=3;
                    figure;
                    [hAx,~]=tight_subplot(nSubplotRow,nSubplotCol,[.11 -.3]); maxSubplots=nSubplotRow*nSubplotCol; plotCounter=1;
                    for ds=1:size(datasets,2)
                        Data=datasets{ds};

                        %tiledlayout(1,3, 'Padding', 'loose', 'TileSpacing', 'loose');
                       
                        % --- grab data for 3 lines, baseline/con/incon ---
                        baselineData=Data(Data(:,4)==0,:);
                        congruentData=Data(Data(:,4)==1,:);
                        incongruentData=Data(Data(:,4)==-1,:);
                       
                        % --- split early and late ---
                        [ortBL1,outcomeBL1,ortBL2,outcomeBL2]=splitData(baselineData,2,'trial');
                        [ortCon1,outcomeCon1,ortCon2,outcomeCon2]=splitData(congruentData,2,'trial');
                        [ortIncon1,outcomeIncon1,ortIncon2,outcomeIncon2]=splitData(incongruentData,2,'trial');
                        % --- nosplit just extract columns ---
                        [ortBL,outcomeBL,~,~]=splitData(baselineData,1,'trial');
                        [ortCon,outcomeCon,~,~]=splitData(congruentData,1,'trial');
                        [ortIncon,outcomeIncon,~,~]=splitData(incongruentData,1,'trial');
                        
                        
                        % --- preprocess early data ---
                        [ortBL_x, correctBL_x, errorBL_x]=collapseData(ortBL1, outcomeBL1);
                        [ortCon_x, correctCon_x, errorCon_x]=collapseData(ortCon1, outcomeCon1);
                        [ortIncon_x, correctIncon_x, errorIncon_x]=collapseData(ortIncon1, outcomeIncon1);
                        dsEarly={ortBL_x,ortCon_x,ortIncon_x; correctBL_x,correctCon_x,correctIncon_x; errorBL_x,errorCon_x,errorIncon_x};
                        % --- fit early ---
                        axes(hAx(plotCounter));
                        if sum(dsEarly{1,1})<0
                            xMin=0;%-40;
                            xMax=40;%0;
                        else
                            xMin=0;
                            xMax=40;
                        end

                        for lineNo=1:size(dsEarly,2)
                            [thresholdGiac, shape, bias, lapse, likelihood, ver, threshold] = FitPsyMLConstant(...
                                dsEarly{1,lineNo},dsEarly{2,lineNo},dsEarly{3,lineNo},nTrials,sessionIDs{session},...
                                gaborSizeWanted,gaborEccWanted,gaborSFWanted,verboseTitle,plotFlag,[xMin,xMax],cmap(lineNo,:),mkrCond{lineNo});
                            if plotCounter==4
                                ylabel('% Correct')
                                xlabel('Gabor contrast')
                                title('1st half block')
                                legend({'Baseline','Con. optostim','Incon. optostim'},'Location','best')
                            else
                                title('1st half block')
                            end

                            axis square
                            upFontSize(14,.01)
                            addSkippedTicks(0, 35, 5, 'x')
                            addSkippedTicks(0, 100, 10, 'y')
                        end
                         plotCounter=plotCounter+1;%nexttile
                        % +++ add boxplot per con/incon/BL +++

                        
                        % --- preprocess late data ---
                        [ortBL_x, correctBL_x, errorBL_x]=collapseData(ortBL2, outcomeBL2);
                        [ortCon_x, correctCon_x, errorCon_x]=collapseData(ortCon2, outcomeCon2);
                        [ortIncon_x, correctIncon_x, errorIncon_x]=collapseData(ortIncon2, outcomeIncon2);
                        dsLate={ortBL_x,ortCon_x,ortIncon_x; correctBL_x,correctCon_x,correctIncon_x; errorBL_x,errorCon_x,errorIncon_x};
                        % --- fit late ---
                        axes(hAx(plotCounter));
                        for lineNo=1:size(dsLate,2)
                            [thresholdGiac, shape, bias, lapse, likelihood, ver, threshold] = FitPsyMLConstant(...
                                dsLate{1,lineNo},dsLate{2,lineNo},dsLate{3,lineNo},nTrials,sessionIDs{session},...
                                gaborSizeWanted,gaborEccWanted,gaborSFWanted,verboseTitle,plotFlag,[xMin,xMax],cmap(lineNo,:),mkrCond{lineNo});
                            title('2nd half block')
                            axis square
                            upFontSize(14,.01)
                            addSkippedTicks(0, 35, 5, 'x')
                            addSkippedTicks(0, 100, 10, 'y')
                        end
                        upFontSize(14,.01)
                        plotCounter=plotCounter+1;%nexttile
                        
                        
                        % --- preprocess whole block data ---
                        [ortBL_x, correctBL_x, errorBL_x]=collapseData(ortBL, outcomeBL);
                        [ortCon_x, correctCon_x, errorCon_x]=collapseData(ortCon, outcomeCon);
                        [ortIncon_x, correctIncon_x, errorIncon_x]=collapseData(ortIncon, outcomeIncon);
                        dsBlock={ortBL_x,ortCon_x,ortIncon_x; correctBL_x,correctCon_x,correctIncon_x; errorBL_x,errorCon_x,errorIncon_x};
                        % --- fit whole block ---
                        axes(hAx(plotCounter));
                        for lineNo=1:size(dsBlock,2)
                            [thresholdGiac, shape, bias, lapse, likelihood, ver, threshold] = FitPsyMLConstant(...
                                dsBlock{1,lineNo},dsBlock{2,lineNo},dsBlock{3,lineNo},nTrials,sessionIDs{session},...
                                gaborSizeWanted,gaborEccWanted,gaborSFWanted,verboseTitle,plotFlag,[xMin,xMax],cmap(lineNo,:),mkrCond{lineNo});
                            title('Whole block')
                            addSkippedTicks(0, 35, 5, 'x')
                            addSkippedTicks(0, 100, 10, 'y')
                        end
                        axis square
                        upFontSize(14,.01)

                        plotCounter=plotCounter+1;%nexttile
                        if plotCounter==maxSubplots
                            plotCounter=1;
                        end

                    end
                    
                        [~,h]=suplabel(num2str(allSessions),'t',[.08 .08 .84 .86]);
                        [~,h]=suplabel('Horizontal              Vertical','y',[.2 .08 .84 .86]);
                                                set(h,'FontSize',14)

                    export_fig(pdfFilename,'-pdf','-nocrop');
                   
                    
                    %{
                    
                    negative0=-1.*(stimType>0 & gaborOrt==0);
                    tmp=negative0'.*-1;
                    negative0(negative0==0)=1;
                    
                    
                    gaborCon=condData.StimCon.*negative0;
                    vertReportData=successData.*negative0' + completedData.*tmp;
                    vertReportPrctData=(vertReportData./completedData)'.*100;

                    % --- Sort them into blank/gabor 0 (noopto, 0opto, 90opto)/gabor 90 (noopto, 0opto, 90opto) ---
                    blankCond=find(stimType==0);
                    
                    baselineNoOpto=~isempty(visualControl)==1 && isempty(optoControl)==1;
                    
                    switch baselineNoOpto
                        case {1}
                            visCond0control=find(stimType==1 & gaborOrt==0); %[2:4]
                            visCond90control=find(stimType==1 & gaborOrt==90); %[5:7]
                        case {0}
                            visCond0control=find(optoControl==1 & gaborOrt==0);
                            visCond90control=find(optoControl==1 & gaborOrt==90);
                    end

                    vis0_opto0=find(stimType==3 & gaborOrt==0 & opto0);
                    vis0_opto90=find(stimType==3 & gaborOrt==0 & opto90);

                    vis90_opto0=find(stimType==3 & gaborOrt==90 & opto0);
                    vis90_opto90=find(stimType==3 & gaborOrt==90 & opto90);
                    
                    % --- compile into mat --- (change to con/incon optostim if need be)
                    conds0=[visCond0control;vis0_opto0;vis0_opto90];
                    conds90=[visCond90control;vis90_opto0;vis90_opto90];
                    
                    success0=successData(conds0); %horz
                    success90=successData(conds90); %vert
                    successPrct0=vertReportPrctData(conds0);
                    successPrct90=vertReportPrctData(conds90);
                    
                    %
                    mkrCond={'o','v','^'};
                    
                    % colormap
                    cmap=fireice;
                    cmap=cmap([22 1+10 43-10],:);
                    
                   %% Raw, whole session
                    for condSet=1:3
                        conds2extract=[conds0(condSet,:), conds90(condSet,:)]; %[-5 -12 -30] [5 12 30]
                        deltaOrientation=gaborCon(conds2extract)
                        percentHorReport=100-vertReportPrctData(conds2extract)
                        percentVertReport=vertReportPrctData(conds2extract)
                        nTrials=completedData(conds2extract);
                        
                        % ----- Giac -----
                        [thresholdGiac, shape, bias, lapse, likelihood, ver, threshold] = FitPsyMLConstant(...
                            deltaOrientation,percentVertReport,percentHorReport,nTrials,sessionIDs{session},...
                            gaborSize,gaborEccentricity,gaborSFWanted,verboseTitle,plotFlag,[-xMin,xMax],cmap(condSet,:),mkrCond{condSet});
                        ylabel('% Vertical reports')
                        xlabel('Gabor contrast')
                        legend({'Baseline','H stim (0-columns)','V-stim (90-columns)'})
                    end
                    
                    %% Early and late split
                    for condSet=1:3
                        
                        conds2extract=[conds0(condSet,:), conds90(condSet,:)]; %[-5 -12 -30] [5 12 30]
                        deltaOrientation=gaborCon(conds2extract)
                        percentHorReport=100-vertReportPrctData(conds2extract)
                        percentVertReport=vertReportPrctData(conds2extract)
                        nTrials=completedData(conds2extract);
                        
                        % ----- Giac -----
                        [thresholdGiac, shape, bias, lapse, likelihood, ver, threshold] = FitPsyMLConstant(...
                            deltaOrientation,percentVertReport,percentHorReport,nTrials,sessionIDs{session},...
                            gaborSize,gaborEccentricity,gaborSFWanted,verboseTitle,plotFlag,[-xMin,xMax],cmap(condSet,:),mkrCond{condSet});
                        ylabel('% Vertical reports')
                        xlabel('Gabor contrast')
                        legend({'Baseline','H stim (0-columns)','V-stim (90-columns)'})
                    end           
                    %}
                end
            end
        end
    end
end