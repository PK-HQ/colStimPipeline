function [betaSorted,muSorted,sigmaSorted,contrastsSorted,percentVerticalSorted]=analyzePsychometrics(dsCurrentSess,hAx,saveFlag,saveType)
%% Generates split behavioral plots
% Analysis type
taskName='optostim';

% Generate filename
filenameStructCurrent=generateFilenames(dsCurrentSess);

% Set save file name
if strcmp(saveType,'psychfit') ||  strcmp(saveType,'psychometrics') 
    pdfFilename=filenameStructCurrent.psychfitPDF;
elseif strcmp(saveType,'full')
    pdfFilename=filenameStructCurrent.neurometricPDF;
end

% Figure behavior
set(0,'DefaultFigureWindowStyle','docked')

%% Define dirs
setupEnvironment()
%{
if ispc
    pcID=getenv('COMPUTERNAME');
    switch pcID
        case {'CVIS-A64882','PSYC-A77304'} % Upper PC
            userStr='pktan';
        case {'LA-CPSD077020WD','LA-CPSA07019WD'} % Lower PC
            userStr='esexpt';
    end
    monkeyDir='Y:/Chip/PsychometricData/';
    boxDataDir=['C:/Users/' userStr '/Box/_Eyal/Columnar read-write/code/plots/Training/'];
    cd(monkeyDir)
elseif ismac % Laptop
    cd '~/'
    monkeyDir='~/Box/_Eyal/Columnar read-write/data/';
    boxDataDir='~/Box/_Eyal/Columnar read-write/code/plots/Training/';
    cd(monkeyDir)
end
%}

%% Plot psychometric curve per biasing session
plotCurve=1;
plotThreshold=1;
plotBias=1;
figure('Name',['Psychometrics (' dsCurrentSess.date 'R' dsCurrentSess.run ')']);
baselineSeparate=isempty(filenameStructCurrent.baselineTS);

% If opto and baseline conds are in the same block
switch baselineSeparate
    case 0
        runData=load(filenameStructCurrent.TS);
        [hAx,~]=tight_subplot(1,3);%,[.1 .1],[.1 .1],[.1 .1]);
        [betaSorted,muSorted,sigmaSorted,contrastsSorted,percentVerticalSorted]=plotPsychometricFit(dsCurrentSess,runData,taskName,...
        [plotCurve,plotThreshold,plotBias],hAx); hold on;
        axis square
        [~,h]=suplabel([dsCurrentSess.date 'R' dsCurrentSess.run],'t',[.08 .08 .84 .8]);
        upFontSize(24,0.01)
    case 1  % opto and baseline conds in separate block
        runData=load(filenameStructCurrent.TS);
        baselineData=load(filenameStructCurrent.baselineTS);
        [hAx,~]=tight_subplot(1,3);%,[.1 .1],[.1 .1],[.1 .1]);
        [betaBaseline,muBaseline,sigmaBaseline,contrastsBaseline,percentVerticalBaseline]=plotPsychometricFit(dsCurrentSess,baselineData,taskName,...
        [plotCurve,plotThreshold,plotBias],hAx); hold on
        [betaOptostim,muOptostim,sigmaOptostim,contrastsOptostim,percentVerticalOptostim]=plotPsychometricFit(dsCurrentSess,runData,taskName,...
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

        axis square
        [~,h]=suplabel([dsCurrentSess.date 'R' dsCurrentSess.run],'t',[.08 .08 .84 .8]);
        upFontSize(24,0.01)
end
hold off

% Save
switch saveFlag
  case {1}
    export_fig(pdfFilename,'-pdf','-nocrop','-append','-transparent');
end
end