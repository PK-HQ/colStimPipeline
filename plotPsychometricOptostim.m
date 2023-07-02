function deltaSorted=plotPsychometricOptostim(dataDir,taskName,plotFlags,verboseTitle,allSessions,mergeBothOrientations,mergeAcrossBlocks)
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
%{
sessionIDs={};
wantedSessions=allSessions;
count=1;
for sess=wantedSessions
    wantedSession=contains(dataFilenames,num2str(sess));
    if sum(wantedSession)>0
        sessionIDs{count}=num2str(sess);
        count=count+1;
    end
end
%}

sessionIDs{1}=allSessions; %temp hack 

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
        pdfFilename=[dataDir 'M' monkeyIDs{monkey} 'Z' num2str(gaborSizeWanted*100,'%05g') 'S' num2str(gaborSFWanted*100,'%05g')...
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
                    
                    prct90opto=sum(vertReportOpto90(:))*100/sum(completedDataOpto90(:));
                    prct0opto=sum(vertReportOpto0(:))*100/sum(completedDataOpto0(:));
                    prctBaselineopto=sum(vertReportOptoControl(:))*100/sum(completedDataOptoControl(:));

                    fprintf('\nH-0 optostim: %f %% vertical reports\n',...
                        prct0opto)
                    fprintf('Blank optostim: %f %% vertical reports\n',...
                        prctBaselineopto)
                    fprintf('V-90 optostim: %f %% vertical reports\n\n',...
                        prct90opto)
                    %
                    mkrCond={'o','v','^'};
                    
                    % colormap
                    cmap=fireice;
                    cmap=cmap([22 1+10 43-10],:);
                    for condSet=1:3 %baseline, opto0, opto90
                        conds2extract=[conds0(condSet,:), conds90(condSet,:)]; %[-5 -12 -30] [5 12 30]
                        deltaOrientation=gaborCon(conds2extract);
                        percentHorReport=100-vertReportPrctData(conds2extract);
                        percentVertReport=vertReportPrctData(conds2extract);
                        nTrials=completedData(conds2extract);
                        
                        [~,yScaling] = fitPsyDualCDF(deltaOrientation,percentVertReport,sessionIDs{session},cmap(condSet,:),mkrCond{condSet});
                        delta(condSet)=yScaling;
                    end
                    
                    deltaSorted=[delta(2) delta(1) delta(3)];
                    % Add title
                    h=title({sessionIDs{session},...
                     ['${\Delta}$: ' num2str(deltaSorted(1),'%.2f') ', '  num2str(deltaSorted(2),'%.2f') ', ' num2str(deltaSorted(3),'%.2f')]},...
                     'FontWeight','normal','FontSize',14,'FontName','Arial','interpreter','latex');

                    pubfig(h,12,2,0.005);
                    
                    % axis labels
                    ylabel('% Vertical reports')
                    xlabel('Gabor orientation')
                    
                    % legend
                    legend({'Baseline','H stim (0\circ)','V-stim (90\circ)'})

                    upFontSize(14,0.01)
                    % ax ticks
                    addSkippedTicks(-60,60,5,'x')
                    addSkippedTicks(0,100,12.5,'y')
                end
            end
        end
    end
end
        
%{
        
        

                    allDeltaOrientations=[];
                    allCW=[];
                    allCCW=[];
                    allnTrialsPerCond=[];
                    mergedSessionData=[];
                    sessionOutputData=[1;1;1;1];
                    sessionOutputDataUnsort=[];

                    %add xy labels to the bottom leftmost plot
                    if plotID==1%max(intersect(1:nSessions,1:nCol:nRow*nCol))
                        xlabel(xaxisMeasureStr,'FontSize', 12);
                        ylabel('Percent CW (%)','FontSize', 12);
                    else
                        xlabel(' ');
                        ylabel(' ');
                    end
                    %threshold=randi([3 5],1,1);
                    %bias=randf(-0.5,0.5,1);


                        %save
                    if pdfcreated==0 && threshold>5 && threshold<30
                        export_fig(pdfFilename,'-pdf','-nocrop','-r600','-png','-c2,2,2,2');
                        %save threshold (alpha)
                        allThresholds(end+1)=threshold;
                        %save bias (behavioral bias, should be close to 0)
                        allBias(end+1)=bias;
                        sessionsUsed{end+1}=sessionIDs{session};
                        pdfcreated=1;
                    elseif pdfcreated==1 && threshold>5 && threshold<30
                        export_fig(pdfFilename,'-append','-pdf','-nocrop','-r600','-png','-c2,2,2,2');
                        %save threshold (alpha)
                        allThresholds(end+1)=threshold;
                        %save bias (behavioral bias, should be close to 0)
                        allBias(end+1)=bias;
                        sessionsUsed{end+1}=sessionIDs{session};
                    end
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    %get spatial frequency wanted only (0.5 or 2 or 8), else skip
                    gaborSize=runData.TS.Header.Conditions.GaborSize;%runData.TS.Header.Conditions.StimSize;
                    gaborPosition=runData.TS.Header.Conditions.StimPos;
                    gaborEccentricity=sqrt(gaborPosition(1)^2+gaborPosition(2)^2);
                    gaborSF=runData.TS.Header.Conditions.GaborSF;
                    stimTime=runData.TS.Header.Conditions.DurationTOn;
                    
                    
                    
                    
                  case {'coarseOD'}
                    xaxisMeasure=runData.TS.Header.Conditions.StimCon;%runData.TS.Header.Conditions.ContrastCond;
                    xaxisMeasure(length(xaxisMeasure)/2+1:end)=xaxisMeasure(length(xaxisMeasure)/2+1:end)*-1;
                    %xaxisMeasure=xaxisMeasure.*repmat([1 -1],1,numel(xaxisMeasure)/2);
                    %xaxisMeasureOrt=runData.TS.Header.Conditions.OrtCond;
                    xaxisMeasureStr='Contrast (%)';
                  case {'fineOD'}
                    xaxisMeasure=runData.TS.Header.Conditions.OrtCond;
                    xaxisMeasureStr=xaxisMeasureStr;
                end
                %print
                %("\n[%s-R%02g] SF:%gcpd, Sz:%g°, E:%.3g°, Ort:±%.3g\n",sessionIDs{session},run,spatialFrequency,gaborSize,...
                %    gaborEccentricity, unique(abs(xaxisMeasure)))
                if unique(gaborSize)<=gaborSizeWanted...
                        & unique(gaborSF)>=gaborSFWanted...
                        & abs(gaborEccentricity-gaborEccentricityWanted)<=1 ...
                        & unique(stimTime)==stimTimeWanted
                %if gaborSize==gaborSizeWanted & gaborSF==gaborSFWanted & gaborEccentricity>=gaborEccentricityWanted
                    %success is arranged by [O1,P1 O2,P2 O3,P3 O4,P4]
                    successData=runData.TS.Header.Outcomes.CountCondSuccess;
                    nTrialsPerCondData=runData.TS.Header.Outcomes.CountCondTotalValid;
                    nTrialsCompletedData=runData.TS.Header.Outcomes.CountBlockCompleted;
                    %collapse data into one CW & CCW value per DO, e.g. -3 (10,3) and -3 (12,1) into -3 (22,4)
                    nTrialsPerCond=[];
                    success=[];
                    uniqueXaxisMeasure=unique(xaxisMeasure);

                    for nthDeltaOrt=1:length(uniqueXaxisMeasure)
                        DO=uniqueXaxisMeasure(nthDeltaOrt);
                        DOidx=find(xaxisMeasure==DO);
                        nTrialsPerCond=[nTrialsPerCond,sum(nTrialsPerCondData(DOidx))];
                        success=[success,sum(successData(DOidx))];
                    end

                    CWtrials=uniqueXaxisMeasure<=0; %-ve or 0
                    CCWtrials=uniqueXaxisMeasure>0; %+ve

                    %get % CW report for CW and % CCW report for CCW
                    CW=[nTrialsPerCond(CWtrials)-success(CWtrials), success(CCWtrials)];
                    CCW=[success(CWtrials), nTrialsPerCond(CCWtrials)-success(CCWtrials)];

                    if CW+CCW~=nTrialsPerCond;
                        disp('Error in converting success to CW/CCW!')
                    end

                    %compile values across runs of a session for plotting later
                    allDeltaOrientations=[allDeltaOrientations,unique(xaxisMeasure)];
                    allCW=[allCW,CW];
                    allCCW=[allCCW,CCW];
                    allnTrialsPerCond=[allnTrialsPerCond,nTrialsPerCond];

                   %% --------------------
                    %put data together into a 2D matrix for sorting by delta orient.
                    sessionOutputDataUnsort=[allDeltaOrientations;allCW;allCCW;allnTrialsPerCond];

                    %sort by delta orient, preserving the CW/CCW of each delta orient
                    %value
                    sessionOutputData = sortrows(sessionOutputDataUnsort.',1).';

                    %% merge data for both CW/CCW orientations
                    absUniqueXaxisMeasure=unique(abs(sessionOutputData(1,:)));
                    
                    switch mergeBothOrientations
                        case {1}
                            for i=1:length(absUniqueXaxisMeasure)
                                %get columns with same delta orientation
                                duplicateDOIdx1=find(sessionOutputData(1,:)==absUniqueXaxisMeasure(i));
                                duplicateDOIdx2=find(sessionOutputData(1,:)==-absUniqueXaxisMeasure(i));

                                %get row 3 (numerator) and row 4 (denom) values for -ve k
                                negKcorrect=sum(sessionOutputData(3,duplicateDOIdx2));
                                negKtotal=sum(sessionOutputData(4,duplicateDOIdx2));

                                %get row 2 (numerator) and row 4 (denom) values for +ve k
                                posKcorrect=sum(sessionOutputData(2,duplicateDOIdx1));
                                posKtotal=sum(sessionOutputData(4,duplicateDOIdx1));

                                %sum num and denom
                                Kcorrect=negKcorrect+posKcorrect;
                                Ktotal=negKtotal+posKtotal;

                                %save as abs(k)
                                mergedSessionData(1,i)=abs(absUniqueXaxisMeasure(i));
                                mergedSessionData(2,i)=Kcorrect*100/Ktotal;
                                mergedSessionData(3,i)=100-(Kcorrect*100/Ktotal);
                                mergedSessionData(4,i)=Ktotal;
                            end
                            sessionOutputData=mergedSessionData;
                        case {0}
                            %convert corr/err counts to %
                            %sessionOutputData(2,:)=sessionOutputData(2,:)*100./sessionOutputData(4,:);
                            %sessionOutputData(3,:)=sessionOutputData(3,:)*100./sessionOutputData(4,:);
                    end
                    
                else
                    continue
                end
            else
                %dont
            end
            
        end
        
        
      %% Merge across blocks
        uniqueXaxisMeasure=unique((sessionOutputData(1,:)));

        switch mergeAcrossBlocks
            case {1}
                for i=1:length(uniqueXaxisMeasure)
                    %get columns with same delta orientation
                    duplicateDOIdx1=find(sessionOutputData(1,:)==uniqueXaxisMeasure(i));

                    %get row 2 (numerator) and row 4 (denom) values for +ve k
                    posKcorrect=sum(sessionOutputData(2,duplicateDOIdx1));
                    totalTrials=sum(sessionOutputData(4,duplicateDOIdx1));

                    %save as abs(k)
                    mergedSessionData(1,i)=(uniqueXaxisMeasure(i));
                    mergedSessionData(2,i)=posKcorrect*100/totalTrials;
                    mergedSessionData(3,i)=100-(posKcorrect*100/totalTrials);
                    mergedSessionData(4,i)=totalTrials;
                end
                sessionOutputData=mergedSessionData;
            case {0}
                %convert corr/err counts to %
                sessionOutputData(2,:)=sessionOutputData(2,:)*100./sessionOutputData(4,:);
                sessionOutputData(3,:)=sessionOutputData(3,:)*100./sessionOutputData(4,:);
        end
       %% add remove datapoint if >30% con, and <80% performance
        %fit and plot curve
        deltaOrientation = sessionOutputData(1,:);
        percentCWreport = sessionOutputData(2,:);
        percentCCWreport = sessionOutputData(3,:);
        nTrials = sessionOutputData(4,:);
        sessionStr= [sessionIDs{session}];% '-R' num2str(run)];
        
        
        
        if numel(unique(deltaOrientation))/2>=minimumConds %if has >=3 n contrasts
            % ----- Weibull -----
            % data
            % Load in your results (trial wise, not combined)
            results.response=[] %1 or 0, hit or miss
            results.intensity=[] %range=0:1

            % Strip out the trials where there was an invalid response (if there are any).
            
            goodTrials = ~isnan(results.response);
            results.response = results.response(goodTrials);
            results.intensity = results.intensity(goodTrials);
            
            % Starting parameters
            pInit.t = .1; %threshold
            pInit.b = 2; %slope

            [pBest,logLikelihoodBest] = fit('fitPsychometricFunction',pInit,{'b','t'},results,'Weibull')
            pBest 
            logLikelihoodBest 
            
            % log likelihood surface
            plot(log(pBest.t),pBest.b,'o','MarkerFaceColor','k');
            
            y = Weibull(pBest,x);
            figure(1)
            plot(log(x),100*y,'g-','LineWidth',2);
            threshold=pBest.t; %80%
            
            % ----- Giac -----
            [thresholdGiac, shape, bias, lapse, likelihood, ver, threshold] = FitPsyML(...
                deltaOrientation,percentCWreport,...
                percentCCWreport,nTrials,sessionStr,...
                gaborSize,gaborEccentricity,gaborSFWanted,verboseTitle,plotFlag,[-60,60]);
            
            allDeltaOrientations=[];
            allCW=[];
            allCCW=[];
            allnTrialsPerCond=[];
            mergedSessionData=[];
            sessionOutputData=[1;1;1;1];
            sessionOutputDataUnsort=[];

            %add xy labels to the bottom leftmost plot
            if plotID==1%max(intersect(1:nSessions,1:nCol:nRow*nCol))
                xlabel(xaxisMeasureStr,'FontSize', 12);
                ylabel('Percent CW (%)','FontSize', 12);
            else
                xlabel(' ');
                ylabel(' ');
            end
            %threshold=randi([3 5],1,1);
            %bias=randf(-0.5,0.5,1);


                %save
            if pdfcreated==0 && threshold>5 && threshold<30
                export_fig(pdfFilename,'-pdf','-nocrop','-r600','-png','-c2,2,2,2');
                %save threshold (alpha)
                allThresholds(end+1)=threshold;
                %save bias (behavioral bias, should be close to 0)
                allBias(end+1)=bias;
                sessionsUsed{end+1}=sessionIDs{session};
                pdfcreated=1;
            elseif pdfcreated==1 && threshold>5 && threshold<30
                export_fig(pdfFilename,'-append','-pdf','-nocrop','-r600','-png','-c2,2,2,2');
                %save threshold (alpha)
                allThresholds(end+1)=threshold;
                %save bias (behavioral bias, should be close to 0)
                allBias(end+1)=bias;
                sessionsUsed{end+1}=sessionIDs{session};
            end
            plotID=plotID+1;
            
            
        else
            plotID=plotID;
        end
        %[ax1,h1]=suplabel(['Contrast'],'Y');
        %set(h1,'FontSize',16)
        %set(h1,'FontWeight','bold')
        hold off;
        
    end
    

    
end
hold off;

%% Simulate 10 sessions with rand data (comment out if not used)
allThresholds=allThresholds;%randi([3 5],nSessionsExample,1); %allThresholds
allBias=allBias;%randf(-0.5,0.5,nSessionsExample)';
%nSessions=nSessionsExample;
filteredAllThresholds=allThresholds(abs(allBias)<10);
filteredAllBias=allBias(abs(allBias)<10);
%% Plot threshold across sessions
switch plotFlags(2)
    case {1}
        for monkey=1:nMonkey
            %get monkey's threshold
            sessionThreshold=mean(filteredAllThresholds,1);
            nSessions=numel(sessionThreshold);
            if nSessions < 10
                maxSessions=10;
            else
                maxSessions=nSessions;
            end
            %get mean of last 5 points
            if length(sessionThreshold)<3
                wantedSessionThresholds=sessionThreshold(:);
                meanLastSessions=mean(wantedSessionThresholds);
            else
                wantedSessionThresholds=sessionThreshold(end-lastNsessions:end);
                meanLastSessions=mean(wantedSessionThresholds);
            end
            
            %plot staircase/line plot of threshold across sessions (threshold x session #)
            figure
            
            h1=plot(1:length(sessionThreshold),sessionThreshold,'Color','k','LineStyle','-','markersize',35,'marker','.','markerfacecolor','w','markeredgecolor','k');hold on
            h1=plot(1:length(sessionThreshold),sessionThreshold,'Color','g','LineStyle','-','markersize',25,'marker','.','markerfacecolor','w','markeredgecolor','g');hold on

            %line(xlim,[meanLastSessions meanLastSessions],'Color','k','LineStyle','--');
            
            %ax limits
            set(gca,'xtick',1:numel(sessionsUsed));
            set(gca,'xticklabel',sessionsUsed,'fontsize',10)
            xtickangle(45)
            
            ylim([0, 30])
            xlim([0 maxSessions])
            
            %text
            xlabel('Session number')
            ylabel('Threshold contrast (%)')
            title(['Threshold across sessions, last ' num2str(length(wantedSessionThresholds)) ' sessions = ' num2str(mean(meanLastSessions),'%.1f') '\circ'])
            pubfig(h1,12,2,0.015);

            %save
            export_fig(pdfFilename,'-append','-pdf','-nocrop','-r600','-png','-c2,2,2,2');
        end
end

%% Plot bias across sessions
switch plotFlags(3)
    case {1}
    for monkey=1:nMonkey
        %get monkey's threshold
        sessionBiases=filteredAllBias(:);
        
        
        %get mean of last 5 points
        if nSessions<3
            meanLastSessions=mean(sessionBiases);
        else
            meanLastSessions=mean(sessionBiases(end-lastNsessions:end));
        end

        %plot staircase/line plot of bias across sessions (threshold x session #)
        figure
        %line(xlim,[meanLastSessions meanLastSessions],'Color','k','LineStyle','--');
        h1=plot(1:length(sessionBiases),sessionBiases,'Color','k','LineStyle','-','markersize',35,'marker','.','markerfacecolor','w','markeredgecolor','k');hold on
        h1=plot(1:length(sessionBiases),sessionBiases,'Color','r','LineStyle','-','markersize',25,'marker','.','markerfacecolor','w','markeredgecolor','r');hold on
        ylim([-10, 10])
        xlim([0 maxSessions])
        line(xlim,[0 0],'Color','k');

        %ax limits
        %ax limits
        set(gca,'xtick',1:numel(sessionsUsed));
        set(gca,'xticklabel',sessionsUsed,'fontsize',10)
        xtickangle(45)        

        
        %text
        xlabel('Session number')
        ylabel('Bias (%)')
        title(['Bias across sessions, last ' num2str(length(sessionBiases)) ' sessions = ' num2str(mean(meanLastSessions),'%.1f') '\circ'])
        pubfig(h1,12,2,0.015);
        

        %save
        export_fig(pdfFilename,'-append','-pdf','-nocrop','-r600','-png','-c2,2,2,2');
    end
end

end
%}
