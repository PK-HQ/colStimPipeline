function plotPsychometricOptostimV3(dataDir,outputDir,taskName,plotFlags,verboseTitle,allSessions,unwantedSessions,mergeBothOrientations,mergeAcrossBlocks)
%% Plots psychometric curves, thresholds, biases across sessions per monkey
%analysis params
minimumConds=3;
gaborEccentricityWanted=1.46;%1.46;%5;
gaborSizeWanted=1;%1;%1
gaborSFWanted=6;%6;%4
stimTimeWanted=300;%6;%4


%{
gaborSize<=gaborSizeWanted...
& gaborSF>=gaborSFWanted...
& abs(gaborEccentricity-gaborEccentricityWanted)<=1
%}
fprintf('Analysis: Z=%g, SF=%g, E=%g \n',gaborSizeWanted,gaborSFWanted,gaborEccentricityWanted)

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
            'E' num2str(gaborEccentricityWanted*100,'%05g') 'N' num2str(nSessions,'%03g')];

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


            
            
            
            
            if isfield(runData.TS.Header.Conditions,'GaborSize')%'GaborSize')
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
                    condsOptoControl=[visCond0control;visCond90control];
                    condsOpto0=[vis0_opto0;vis90_opto0];
                    condsOpto90=[vis0_opto90;vis90_opto90];

                    vertReportOptoControl=vertReportData(condsOptoControl); %control
                    vertReportOpto0=vertReportData(condsOpto0); %horz
                    vertReportOpto90=vertReportData(condsOpto90); %vert
                       
                    completedDataOptoControl=completedData(condsOptoControl); %horz
                    completedDataOpto0=completedData(condsOpto0); %horz
                    completedDataOpto90=completedData(condsOpto90); %vert
                    
                    
                    
                    fprintf('\nBlank optostim: %.0d vertical reports / %.0d completed trials\n',...
                        sum(vertReportOptoControl(:)), sum(completedDataOptoControl(:)))
                    fprintf('H optostim: %.0d vertical reports / %.0d completed trials\n',...
                        sum(vertReportOpto0(:)), sum(completedDataOpto0(:)))
                    fprintf('V optostim: %.0d vertical reports / %.0d completed trials\n\n',...
                        sum(vertReportOpto90(:)), sum(completedDataOpto90(:)))
                    %
                    mkrCond={'o','v','^'};
                    
                    % colormap
                    cmap=fireice;
                    cmap=cmap([22 1+10 43-10],:);
                    for condSet=1:3
                        
                        conds2extract=[conds0(condSet,:), conds90(condSet,:)]; %[-5 -12 -30] [5 12 30]
                        deltaOrientation=gaborCon(conds2extract)
                        percentHorReport=100-vertReportPrctData(conds2extract)
                        percentVertReport=vertReportPrctData(conds2extract)
                        
                        fitPsyDualCDF(...
                            deltaOrientation,percentVertReport,sessionIDs{session},...
                            cmap(condSet,:),mkrCond{condSet});
                        ylabel('% Vertical reports')
                        xlabel('Gabor orientation')
                        legend({'Baseline','H stim (0-columns)','V-stim (90-columns)'})
                        
                        
                    end
                end
            end
        end
    end
end