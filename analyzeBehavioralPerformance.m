%% Main script to run all behavioral and imaging analyses
%% Figure behavior
set(0,'DefaultFigureWindowStyle','docked')
taskName='optostim';

%% Define dirs
if ispc
    pcID=getenv('COMPUTERNAME');
    switch pcID
        case {'SEIDEMANN2PANAL'}
            userStr='pktan';
        case {'LA-CPSD077020WD','LA-CPSA07019WD'}
            userStr='esexpt';
    end
    monkeyDir='Y:/Chip/PsychometricData/'; %'Y:/Apollo';
    boxDataDir=['C:/Users/' userStr '/Box/_Eyal/Columnar read-write/code/plots/Training/'];
    cd(monkeyDir)
elseif ismac
    cd '~/'
    monkeyDir='~/Box/_Eyal/Columnar read-write/data/'; %bc data is here not lab drive
    boxDataDir='~/Box/_Eyal/Columnar read-write/code/plots/Training/';
    cd(monkeyDir)
end
switch taskName
  case {'coarseOD'}
    taskStr='gaborOD/';%'roe';
  case {'fineOD'}
    taskStr='gabor_nostim/';
  case {'optostim'}
    taskStr='optostim/';

end

%% Define datastruct that contains info about runs
%open('Y:\users\PK\colStimPipeline\exptListBiasing.m')
run('Y:\users\PK\colStimPipeline\exptListBiasing.m')

%% Save filename
monkeyID=28;
gaborEcc=1.46;%1.46;%5;
gaborSize=1;%1;%1
gaborSF=6;%6;%4
nSessions=numel(2:length(datastruct));
pdfFilename=['Y:/Chip/Meta/plots/' 'M28' num2str(monkeyID) 'Z' num2str(gaborSize*100,'%05g') 'S' num2str(gaborSF*100,'%05g')...
    'E' num2str(gaborEcc*100,'%05g') 'N' num2str(nSessions,'%03g')];

%% Plot psychometric curve per biasing session
deltaSortedCont=[];
saveCounter=1;
for sessionID=2:length(datastruct)
    sessionRange=datastruct(sessionID).date;
    runNo=datastruct(sessionID).run;
    
    mergeBothOrientations=0;
    mergeAcrossBlocks=0;

    %% Plot psychometric curve
    plotCurve=1;
    plotThreshold=1;
    plotBias=1;
    verboseTitle=0;

    monkeyDir=['Y:/Chip/Chip' num2str(sessionRange) '/'];
    boxDataDir=['Y:/Chip/Chip' num2str(sessionRange) '/run' num2str(runNo) '/'];
    cd(monkeyDir)

    dataInputDir = [monkeyDir 'run' num2str(runNo) '/'];

    deltaSorted=plotPsychometricOptostimV4(dataInputDir,taskName,...
        [plotCurve,plotThreshold,plotBias],verboseTitle,sessionRange,mergeBothOrientations,mergeAcrossBlocks);
    deltaSortedCont(sessionID,:)=deltaSorted;
    if saveCounter==1
        export_fig(pdfFilename,'-pdf','-nocrop');
    else
        export_fig(pdfFilename,'-pdf','-nocrop','-append');
    end
    saveCounter=saveCounter+1;
end
deltaSortedContTruncated=deltaSortedCont(2:end,:);
for stimcond=1:3
    fprintf('Cond %.0f: %.2f (%.2f)\n',stimcond,mean(deltaSortedContTruncated(:,stimcond)), std(deltaSortedContTruncated(:,stimcond)))
end




%% hist

% hardcode        (1-seg significant PMC)                                (2-seg significant PMC)
h_opto=deltaSortedContTruncated(:,1);
bl_opto=deltaSortedContTruncated(:,2);
v_opto=deltaSortedContTruncated(:,3);

%% Run MWU test, calculate p-value that both dist different, 2-tailed
[p_hv,~]=ranksum(h_opto,v_opto,'tail','both')
[p_hbl,~]=ranksum(h_opto,bl_opto,'tail','both')
[p_vbl,~]=ranksum(v_opto,bl_opto,'tail','both')

%% Plot distributions and display medians and MWU p-value
figure
xlim([0 1])
ylim([0 10])

%h
histogram(h_opto,'FaceColor','b','FaceAlpha',.6,'BinWidth',0.05,'LineWidth',2); hold on
%bl
histogram(bl_opto,'FaceColor',[.85 .85 .85],'FaceAlpha',.6,'BinWidth',0.05,'LineWidth',2); hold on
%v
histogram(v_opto,'FaceColor','r','FaceAlpha',.6,'BinWidth',0.05,'LineWidth',2); hold on

xlim([0 1])
ylim([0 8])

%median lines
scatter(median(h_opto),6.8,150,'vk','filled','HandleVisibility','off')
scatter(median(h_opto),6.8,60,'v','filled','MarkerFaceColor','b')
scatter(median(bl_opto),6.8,150,'vk','filled','HandleVisibility','off')
scatter(median(bl_opto),6.8,60,'v','filled','MarkerFaceColor',[.85 .85 .85])
scatter(median(v_opto),6.8,150,'vk','filled','HandleVisibility','off')
scatter(median(v_opto),6.8,60,'vk','filled','MarkerFaceColor','r')


xlabel('\Delta')
ylabel('Count')
addSkippedTicks(0,1,.1,'x')

addSkippedTicks(0,8,1,'y')
hLegend=legend({'H-opto','Baseline','V-opto',...
    },'location','eastoutside');
hTitle=title('\Delta distributions for stimulation conditions','FontWeight','Normal');

upFontSize(16,0.01);
hLegend.FontSize=14;
text((median(h_opto)+median(v_opto))/2, 7.4, '***','FontSize',20)

text((median(h_opto)+median(bl_opto))/2, 7.8, 'n.s.','FontSize',20)

text((median(v_opto)+median(bl_opto))/2, 7.4, '***','FontSize',20)












%% Storage
%plotPsychometricOptostimV2(dataInputDir,plotOutputDir,taskName,...
%    [plotCurve,plotThreshold,plotBias],verboseTitle,sessionRange,unwantedSessions,mergeBothOrientations,mergeAcrossBlocks);

% -----  Gabor, contrast series  -----
%{
plotPsychometricCurveSessionContrast(dataInputDir,plotOutputDir,taskName,...
    [plotCurve,plotThreshold,plotBias],verboseTitle,sessionRange,unwantedSessions,mergeBothOrientations,mergeAcrossBlocks)
%}
%plotPsychometricCurve(dataInputDir,plotOutputDir,taskName,...
%    [plotCurve,plotThreshold,plotBias],verboseTitle,sessionRange,unwantedSessions,mergeBothOrientations,mergeAcrossBlocks)

%plotPsychometricCurveBlockContrast(matDataDirRoeContrast,plotOutputDirRoeContrast,taskName,...
%    [plotCurve,plotThreshold,plotBias],verboseTitle,sessionRange,unwantedSessions,mergeBothOrientations)

% -----  ROE, contrast series series  -----
% Sessionwise ROE
%plotPsychometricCurveSession(matDataDirRoe,plotOutputDirRoe,taskName,[plotCurve,plotThreshold,plotBias],verboseTitle,sessionRange,unwantedSessions)
% Blockwise ROE
%plotPsychometricCurveBlock(matDataDirRoe,plotOutputDirRoe,taskName,[plotCurve,plotThreshold,plotBias],verboseTitle,sessionRange,unwantedSessions)

% -----  Unused  -----
%plotPsychometricCurve(matDataDir,plotOutputDir,taskName,[plotCurve,plotThreshold,plotBias],verboseTitle)
%plotPsychometricCurve2(matDataDir,plotOutputDir,taskName,[plotCurve,plotThreshold,plotBias],verboseTitle)

%% Quantify shift and slope change of psychometric curve