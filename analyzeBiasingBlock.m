function [deltaSorted,muSorted,sigmaSorted]=analyzeBiasingBlock(dsCurrentSess,saveFlag)
%% Main script to run all behavioral biasing analyses
filenameStructCurrent=generateFilenames(dsCurrentSess);
pdfFilename=filenameStructCurrent.neurometricPDF;

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
dataInputDir = dsCurrentSess;
[deltaSorted,muSorted,sigmaSorted]=plotPsychometricOptostimV4(dataInputDir,taskName,...
    [plotCurve,plotThreshold,plotBias]);
upFontSize(18,0.005)
export_fig(pdfFilename,'-pdf','-nocrop','-append');

for stimcond=1:3
    fprintf('Cond %.0f: %.2f (%.2f)\n',stimcond,mean(deltaSorted(:,stimcond)), std(deltaSorted(:,stimcond)))
end

end