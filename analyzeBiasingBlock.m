function [betaSorted,muSorted,sigmaSorted,contrastsSorted,percentVerticalSorted]=analyzeBiasingBlock(dsCurrentSess,hAx,saveFlag,saveType)
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
        case {'CVIS-A64882','PSYC-A77304'}
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
if isempty(filenameStructCurrent.baselineTS) % all opto and baseline conds in same block
    runData=load(filenameStructCurrent.TS);
    [betaSorted,muSorted,sigmaSorted,contrastsSorted,percentVerticalSorted]=plotPsychometricOptostimV4(dsCurrentSess,runData,taskName,...
    [plotCurve,plotThreshold,plotBias],hAx); hold on;
    % Add title
    %{
    h=title({dsCurrentSess.date,...
         ['${\beta}$: ' num2str(betaSorted(2),'%.2f') ', '  num2str(betaSorted(1),'%.2f') ', ' num2str(betaSorted(3),'%.2f')]},...
         'FontWeight','normal','FontSize',14,'FontName','Arial','interpreter','latex');
    %}
    h=title({dsCurrentSess.date,...
         ['\beta: ' num2str(betaSorted(2),'%.2f') ', '  num2str(betaSorted(1),'%.2f') ', ' num2str(betaSorted(3),'%.2f')]},...
         'FontWeight','normal','FontSize',14,'FontName','Arial','interpreter','tex');
    axis square
    switch saveFlag
      case {1}
        export_fig(pdfFilename,'-pdf','-nocrop','-append','-transparent');
    end
else  % opto and baseline conds in separate block
    runData=load(filenameStructCurrent.TS);
    baselineData=load(filenameStructCurrent.baselineTS);
    [betaBaseline,muBaseline,sigmaBaseline,contrastsBaseline,percentVerticalBaseline]=plotPsychometricOptostimV4(dsCurrentSess,baselineData,taskName,...
    [plotCurve,plotThreshold,plotBias],hAx); hold on
    [betaOptostim,muOptostim,sigmaOptostim,contrastsOptostim,percentVerticalOptostim]=plotPsychometricOptostimV4(dsCurrentSess,runData,taskName,...
    [plotCurve,plotThreshold,plotBias],hAx); hold on;
    betaSorted=[betaOptostim(2) betaBaseline(2) betaOptostim(1)];
    muSorted=[muOptostim(2) muBaseline(1) muOptostim(1)];
    sigmaSorted=[sigmaOptostim(2) sigmaBaseline(1) sigmaOptostim(1)];
    
    % if number of contrasts dont match in baseline and optostim blocks
    condDifference=numel(contrastsOptostim(:,1))-numel(contrastsBaseline(:,1));
    if condDifference>0
        contrastsBaseline(end+condDifference,:)=NaN;
        percentVerticalBaseline(end+condDifference,:)=NaN;
    end
    
    contrastsSorted=[contrastsOptostim(:,2) contrastsBaseline(:,2) contrastsOptostim(:,1)];
    percentVerticalSorted=[percentVerticalOptostim(:,2) percentVerticalBaseline(:,2) percentVerticalOptostim(:,1)];

    % Add title
    %{
    h=title({[dsCurrentSess.date 'R' dsCurrentSess.run],...
        ['Mean columns: ' num2str(mean(dsCurrentSess.nColumns))],...
         ['${\beta}$: ' num2str(betaOptostim(2),'%.2f') ', '  num2str(betaBaseline(2),'%.2f') ', ' num2str(betaOptostim(1),'%.2f')]},...
         'FontWeight','normal','FontSize',14,'FontName','Arial','interpreter','latex');
     %}
    temporalpower=unique(runData.TS.Header.Conditions.ProjTTLPulseOn);
     h=title({[dsCurrentSess.date 'R' dsCurrentSess.run],...
        ['Columns: ' num2str(mean(dsCurrentSess.nColumns)) ' (' num2str(temporalpower) 'ms)'],...
         ['\beta: ' num2str(betaOptostim(2),'%.2f') ', '  num2str(betaBaseline(2),'%.2f') ', ' num2str(betaOptostim(1),'%.2f')]},...
         'FontWeight','normal','FontSize',14,'FontName','Arial','interpreter','tex');
    axis square
    upFontSize(24, .01)
    switch saveFlag
      case {1}
        export_fig(pdfFilename,'-pdf','-nocrop','-append','-transparent');
    end
end
for stimcond=1:3
    fprintf('Cond %.0f: %.2f (%.2f)\n',stimcond,mean(betaSorted(:,stimcond)), std(betaSorted(:,stimcond)))
end
hold off

end