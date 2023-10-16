function [betaSorted,muSorted,sigmaSorted]=analyzeBiasingBlock(dsCurrentSess,saveFlag,saveType)
%% Main script to run all behavioral biasing analyses
filenameStructCurrent=generateFilenames(dsCurrentSess);
if strcmp(saveType,'psychfit')
    pdfFilename=filenameStructCurrent.psychfitPDF;
elseif strcmp(saveType,'full')
    pdfFilename=filenameStructCurrent.neurometricPDF;
end
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

%% Save filename
monkeyID=28;
gaborEcc=1.46;%1.46;%5;
gaborSize=1;%1;%1
gaborSF=6;%6;%4
nSessions=1;

%% Plot psychometric curve per biasing session
plotCurve=1;
plotThreshold=1;
plotBias=1;
if isempty(filenameStructCurrent.baselinerun) % all opto and baseline conds in same block
    figure('Name','Psychometrics'); hold on
    runData=load(filenameStructCurrent.TS);
    [betaSorted,muSorted,sigmaSorted]=plotPsychometricOptostimV4(dsCurrentSess,runData,taskName,...
    [plotCurve,plotThreshold,plotBias]);
    xticks(-100:10:100)
    hold off
    xtickangle(gca,0)
    upFontSize(18,0.005)
    axis square
    switch saveFlag
      case {1}
        export_fig(pdfFilename,'-pdf','-nocrop','-append');
    end
else  % opto and baseline conds in separate block
    runData=load(filenameStructCurrent.TS);
    baselineData=load(filenameStructCurrent.baselinerun);
    figure('Name','Psychometrics'); hold on
    [betaBaseline,muBaseline,sigmaBaseline]=plotPsychometricOptostimV4(dsCurrentSess,baselineData,taskName,...
    [plotCurve,plotThreshold,plotBias]); hold on
    [betaOptostim,muOptostim,sigmaOptostim]=plotPsychometricOptostimV4(dsCurrentSess,runData,taskName,...
    [plotCurve,plotThreshold,plotBias]);
    betaSorted=[betaOptostim(2) betaBaseline(2) betaOptostim(1)];
    muSorted=[muOptostim(2) muBaseline(1) muOptostim(1)];
    sigmaSorted=[sigmaOptostim(2) sigmaBaseline(1) sigmaOptostim(1)];
    xticks(-100:10:100)
    hold off
    xtickangle(gca,0)
    upFontSize(18,0.005)
    axis square
    switch saveFlag
      case {1}
        export_fig(pdfFilename,'-pdf','-nocrop','-append');
    end
end


for stimcond=1:3
    fprintf('Cond %.0f: %.2f (%.2f)\n',stimcond,mean(betaSorted(:,stimcond)), std(betaSorted(:,stimcond)))
end

end