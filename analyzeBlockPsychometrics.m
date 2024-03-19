function [behavioralData]=...
    analyzeBlockPsychometrics(filenameStructCurrent, behavioralData, bitmapData, blockID, ...
    pdfFilename, reportType, saveFlag)
disp('Analyzing psychometrics...')

% Plot psychometric curve per biasing session
betaSorted=1;
muSorted=1;
sigmaSorted=1;
contrastsSorted=1;
percentVerticalSorted=1;
% Initialization
plotCurve = 1;
plotThreshold = 1;
plotBias = 1;
baselineSeparate = ~isempty(filenameStructCurrent.baselineTS);

% Prepare figure
figure('Name', ['Psychometrics (' filenameStructCurrent.date 'R' filenameStructCurrent.run ')']);

% Load run data
runData = load(filenameStructCurrent.TS);

% Setup subplot axes
[hAx, ~] = tight_subplot(1, 3); % Adjust subplot parameters as needed
gaborContrasts.V0=nan(3,6);gaborContrasts.V90=nan(3,6);
percentageCorrect.V0=nan(3,6);percentageCorrect.V90=nan(3,6);
% Conditional plotting based on baseline separation
if baselineSeparate
    % Baseline and opto conditions in separate blocks
    baselineData = load(filenameStructCurrent.baselineTS);
    [gaborContrasts, percentageCorrect] = extractPsychometricData(baselineData,gaborContrasts,percentageCorrect);
else
    % Opto and baseline conditions in the same block
    [betaBaseline, muBaseline, sigmaBaseline, contrastsBaseline, percentVerticalBaseline] = deal([], [], [], [], []); % Initialize for consistency
end

% Plot data for current session
[gaborContrasts, percentageCorrect] = extractPsychometricData(runData,gaborContrasts,percentageCorrect);

% Sort by baseline/congruent/incongruent for each subplot
[gaborContrasts,percentageCorrect]=sortByCongruence(gaborContrasts,percentageCorrect);


% Fit psychometric function and plot
[behavioralData] = fitPsychometricData(behavioralData, hAx, gaborContrasts.baseline, percentageCorrect.baseline, 'baseline', blockID);upFontSize(24, 0.01);
hold on;
[behavioralData] = fitPsychometricData(behavioralData, hAx, gaborContrasts.congruent, percentageCorrect.congruent, 'congruent', blockID);upFontSize(24, 0.01);
hold on;
[behavioralData] = fitPsychometricData(behavioralData, hAx, gaborContrasts.incongruent, percentageCorrect.incongruent, 'incongruent', blockID);upFontSize(24, 0.01);
hold on;
[behavioralData] = fitPsychometricData(behavioralData, hAx, gaborContrasts.congruent, percentageCorrect.congruent-percentageCorrect.incongruent, 'difference', blockID);upFontSize(24, 0.01);
hold on;
axes(hAx(1)); title('Vertical stimuli','FontWeight','normal'); ylabel('Correct (%)'); xlabel('Absolute gabor contrast (%)'); legend({'Baseline','Congruent','Incongruent','Difference'},'Location', 'east', 'FontSize',18);
axes(hAx(2)); title('Horizontal stimuli','FontWeight','normal');
axes(hAx(3)); title('Combined','FontWeight','normal');



% Supertitle
[~,h]=suplabel({'Behavioral biasing',...
    [filenameStructCurrent.date 'R' filenameStructCurrent.run ' (' num2str(bitmapData.nColumns(1,blockID)) ' & ' num2str(bitmapData.nColumns(2,blockID)) ' column @ ' num2str(mean(bitmapData.adjustedSPD_uW(1,:,blockID)),2) ' ' char(0177) ' ' num2str(std(bitmapData.adjustedSPD_uW(1,:,blockID)),2) ' ' char(181) 'W mm^{-2})']},'t',[.1 .1 .795 .78]);
h.FontSize=28;

% Saving 
if saveFlag
    if blockID==1 && strcmp(reportType,'psychometric')
        export_fig(pdfFilename,'-pdf','-nocrop');
    elseif (blockID>1 && strcmp(reportType,'psychometric')) ||  strcmp(reportType,'summary')
        export_fig(pdfFilename,'-pdf','-nocrop','-append');
    end
end
end

% Process and sort data for plotting
%[betaSorted, muSorted, sigmaSorted, contrastsSorted, percentVerticalSorted] = processAndSortData(baselineSeparate, betaBaseline, muBaseline, sigmaBaseline, contrastsBaseline, percentVerticalBaseline, betaCurrent, muCurrent, sigmaCurrent, contrastsCurrent, percentVerticalCurrent);

% Finalize plot
%finalizePlot(dsCurrentSess, betaSorted, muSorted, sigmaSorted, contrastsSorted, percentVerticalSorted);










%% === Processing data ===
function [betaSorted, muSorted, sigmaSorted, contrastsSorted, percentVerticalSorted] = processAndSortData(baselineSeparate, betaBaseline, muBaseline, sigmaBaseline, contrastsBaseline, percentVerticalBaseline, betaCurrent, muCurrent, sigmaCurrent, contrastsCurrent, percentVerticalCurrent)
% Process and sort data based on whether baseline conditions are separate

if baselineSeparate
    % Handle differences in the number of conditions
    condDifference = numel(contrastsCurrent(:, 1)) - numel(contrastsBaseline(:, 1));
    if condDifference > 0
        contrastsBaseline(end + condDifference, :) = NaN;
        percentVerticalBaseline(end + condDifference, :) = NaN;
    end
    % Sort data for separate baseline conditions
    betaSorted = [betaCurrent(2), betaBaseline(2), betaCurrent(1)];
    muSorted = [muCurrent(2), muBaseline(1), muCurrent(1)];
    sigmaSorted = [sigmaCurrent(2), sigmaBaseline(1), sigmaCurrent(1)];
    contrastsSorted = [contrastsCurrent(:, 2), contrastsBaseline(:, 2), contrastsCurrent(:, 1)];
    percentVerticalSorted = [percentVerticalCurrent(:, 2), percentVerticalBaseline(:, 2), percentVerticalCurrent(:, 1)];
else
    % Data sorting for combined conditions
    betaSorted = betaCurrent;
    muSorted = muCurrent;
    sigmaSorted = sigmaCurrent;
    contrastsSorted = contrastsCurrent;
    percentVerticalSorted = percentVerticalCurrent;
end
end

%% === Extract data ===
function [gaborContrasts, percentageCorrect] = extractPsychometricData(runData,gaborContrasts,percentageCorrect)
%% Plots psychometric curves, thresholds, biases across sessions per monkey
% Check if GaborSize field exists in runData
if isfield(runData.TS.Header.Conditions, 'GaborSize')
    % Extract condition and outcome data
    condData = runData.TS.Header.Conditions;
    
    visualOptoTrials= find(condData.TypeCond == 3);

    condType=condData.TypeCond(visualOptoTrials);
    gaborOrt = condData.GaborOrt(visualOptoTrials);
    projectedImageConds = condData.ProjImg(visualOptoTrials);
    successData = runData.TS.Header.Outcomes.CountCondSuccess(visualOptoTrials);
    completedData = runData.TS.Header.Outcomes.CountCondTotalValid(visualOptoTrials);

    % Process condition data
    optoControl = contains(projectedImageConds, 'L000');
    opto0 = contains(projectedImageConds, 'O00000');
    opto90 = contains(projectedImageConds, 'O09000');
    optoBL = contains(projectedImageConds, 'L000');
    gaborCon = condData.StimCon(visualOptoTrials);

    % Initialize variables for data extraction and plotting
    [gaborContrasts, percentageCorrect] = ...
        initializeVariables(gaborContrasts,percentageCorrect, condType, gaborOrt, opto0, opto90, optoBL, visualOptoTrials, optoControl, successData, completedData, gaborCon);
end
end

function [gaborContrasts, percentageCorrect] = initializeVariables(gaborContrasts,percentageCorrect, stimType, gaborOrt, opto0, opto90, optoBL, visualControl, optoControl, successData, completedData, gaborCon)
% Function to initialize variables for data extraction and plotting
gaborOrientations = -1 .* (stimType>0 & gaborOrt == 0);
tmp = gaborOrientations' .* -1;
gaborOrientations(gaborOrientations == 0) = 1;
gaborCon = gaborCon .* gaborOrientations;
vertReportData = successData .* gaborOrientations' + completedData .* tmp;
percentVertical = (vertReportData ./ completedData)' .* 100;

% V0
if sum(gaborOrientations<0 & optoBL==1)>0 % Visual-0 + Opto-BL
    nConds=length(gaborCon(gaborOrientations<0 & optoBL==1));
    gaborContrasts.V0(1,1:nConds)=gaborCon(gaborOrientations<0 & optoBL==1);
    percentageCorrect.V0(1,1:nConds)=100-percentVertical(gaborOrientations<0 & optoBL==1);
end
if sum(gaborOrientations<0 & opto0==1)>0 % V0 + Opto-0
    nConds=length(gaborCon(gaborOrientations<0 & opto0==1));
    gaborContrasts.V0(2,1:nConds)=gaborCon(gaborOrientations<0 & opto0==1);
    percentageCorrect.V0(2,1:nConds)=100-percentVertical(gaborOrientations<0 & opto0==1);
end
if sum(gaborOrientations<0 & opto90==1)>0
    nConds=length(gaborCon(gaborOrientations<0 & opto90==1));
    gaborContrasts.V0(3,1:nConds)=gaborCon(gaborOrientations<0 & opto90==1);
    percentageCorrect.V0(3,1:nConds)=100-percentVertical(gaborOrientations<0 & opto90==1);
end

% V90
if sum(gaborOrientations>0 & optoBL==1)>0
    nConds=length(gaborCon(gaborOrientations>0 & optoBL==1));
    gaborContrasts.V90(1,1:nConds)=gaborCon(gaborOrientations>0 & optoBL==1);
    percentageCorrect.V90(1,1:nConds)=percentVertical(gaborOrientations>0 & optoBL==1);
end
if sum(gaborOrientations>0 & opto0==1)>0
    nConds=length(gaborCon(gaborOrientations>0 & opto0==1));
    gaborContrasts.V90(2,1:nConds)=gaborCon(gaborOrientations>0 & opto0==1);
    percentageCorrect.V90(2,1:nConds)=percentVertical(gaborOrientations>0 & opto0==1);
end
if sum(gaborOrientations>0 & opto90==1)>0
    nConds=length(gaborCon(gaborOrientations>0 & opto90==1));
    gaborContrasts.V90(3,1:nConds)=gaborCon(gaborOrientations>0 & opto90==1);
    percentageCorrect.V90(3,1:nConds)=percentVertical(gaborOrientations>0 & opto90==1);
end

% Sorting conditions
%[conds0, conds90, condsOptoControl, condsOpto0, condsOpto90] = sortConditions(stimType, gaborOrt, opto0, opto90, visualControl, optoControl);
end

function [conds0, conds90, condsOptoControl, condsOpto0, condsOpto90] = sortConditions(stimType, gaborOrt, opto0, opto90, visualControl, optoControl)
% Sort conditions based on the specified criteria
baselineNoOpto = ~isempty(visualControl) == 1 && isempty(optoControl) == 1;
if baselineNoOpto
    visCond0control = find(stimType == 1 & gaborOrt == 0);
    visCond90control = find(stimType == 1 & gaborOrt == 90);
else
    visCond0control = find(optoControl == 1 & gaborOrt == 0);
    visCond90control = find(optoControl == 1 & gaborOrt == 90);
end
vis0_opto0 = find(stimType == 3 & gaborOrt == 0 & opto0);
vis0_opto90 = find(stimType == 3 & gaborOrt == 0 & opto90);
vis90_opto0 = find(stimType == 3 & gaborOrt == 90 & opto0);
vis90_opto90 = find(stimType == 3 & gaborOrt == 90 & opto90);
conds0 = [visCond0control; vis0_opto0; vis0_opto90];
conds90 = [visCond90control; vis90_opto0; vis90_opto90];
condsOptoControl = [visCond0control; visCond90control];
condsOpto0 = [vis0_opto0; vis90_opto0];
condsOpto90 = [vis0_opto90; vis90_opto90];
end

function [gaborContrasts,percentageCorrect]=sortByCongruence(gaborContrasts,percentageCorrect)
gaborContrasts.baseline=[-gaborContrasts.V0(1,:);...
                                               gaborContrasts.V90(1,:); ...
                                               mean([-gaborContrasts.V0(1,:); gaborContrasts.V90(1,:)])];
gaborContrasts.congruent=[-gaborContrasts.V0(2,:);...
                                                  gaborContrasts.V90(3,:); ...
                                                  mean([-gaborContrasts.V0(2,:);gaborContrasts.V90(3,:)])];
gaborContrasts.incongruent=[-gaborContrasts.V0(3,:);...
                                                    gaborContrasts.V90(2,:); ...
                                                    mean([-gaborContrasts.V0(3,:);gaborContrasts.V90(2,:)])];
percentageCorrect.baseline=[percentageCorrect.V0(1,:);...
                                               percentageCorrect.V90(1,:); ...
                                               mean([percentageCorrect.V0(1,:);percentageCorrect.V90(1,:)])];
percentageCorrect.congruent=[percentageCorrect.V0(2,:);...
                                                  percentageCorrect.V90(3,:); ...
                                                  mean([percentageCorrect.V0(2,:);percentageCorrect.V90(3,:)])];
percentageCorrect.incongruent=[percentageCorrect.V0(3,:);...
                                                    percentageCorrect.V90(2,:); ...
                                                    mean([percentageCorrect.V0(3,:);percentageCorrect.V90(2,:)])];
end
%% === Fit and plotting code ===
function behavioralData = fitPsychometricData(behavioralData, hAx, gaborContrastSlice, percentCorrectSlice, optoTypeStr, blockID)
% Fit psychometric function and plot results

switch optoTypeStr
    case 'baseline'
        lineScheme=1; % black with 'o'
    case 'congruent'
        lineScheme=2; % blue with '^'
    case 'incongruent'
        lineScheme=3; % red with 'v'
    case 'difference'
        lineScheme=4; % magenta with 'square'
end

sessionID = 'Session Date'; % Placeholder for session date
cmap=fireice; markerColor = [1 1 1; cmap([43-10 1+10],:); .52 0 .84]; % colors:  Black, blue, red
markerType = {'o-','^-','v-','s--'}; % marker type

% Initialize output variables
beta = []; mu = []; sigma = []; % Placeholders for fitting parameters
param=[];
% Loop through each condition set to fit psychometric function and plot
%line([0 0],[0 100],'Color',.6*[1 1 1],'LineStyle','--','LineWidth',.15,'HandleVisibility','off')
maxcontrast=round(max(gaborContrastSlice(:))/10) * 10;

for set = 1:3
    % Get line data
    x = gaborContrastSlice(set,:);
    y = percentCorrectSlice(set,:);

    % Fit psychometric function (placeholder logic)
    %[betaFit, muFit, sigmaFit, param] = fitPsyFunction(currentX, currentY, sessionID, markerColor, markerType, hAx, half);
    
    % Store fitting results
    %beta = [beta, betaFit]; mu = [mu, muFit]; sigma = [sigma, sigmaFit];
    
    % Plot results for current condition set
    plotFitCurve(x, y, sessionID, markerColor(lineScheme,:), markerType{lineScheme}, hAx, set);
    behavioralData.gaborContrast(lineScheme,set,:,blockID)=x;
    behavioralData.percentageCorrect(lineScheme,set,:,blockID)=y;
    behavioralData.auc(lineScheme,set,blockID)=trapz(x(~isnan(x)),y(~isnan(y))) *100 / trapz(x(~isnan(x)),repmat(100,1,numel(y(~isnan(y))))); % as percentage of possible max
    
    upFontSize(24,0.01); axis square; xlim([0 maxcontrast]);ylim([0 100]);
end

% Sort and organize output parameters for consistency
[betaSorted, muSorted, sigmaSorted, contrastsSorted, percentVerticalSorted] = sortParameters(beta, mu, sigma, gaborContrastSlice, percentCorrectSlice);
end

function [beta, mu, sigma, param] = fitPsyFunction(x, y, sessionID, markerColor, markerType, hAx, half)
% Placeholder function for fitting psychometric curve
% Implement fitting logic as needed, this is a simplified placeholder

% Example fitting results
beta = 0.5; mu = 0; sigma = 1; % Placeholder values
param = [mu, sigma, beta, 0, 0]; % Example parameter vector

% Example log likelihood (not calculated in this placeholder)
logLikelihood = -1;

end

function plotFitCurve(x, y, sessionStr, markerColor, markerType, hAx, set)
% Function to plot fitted curve with scatter plot overlay
% Implement plotting based on provided example

% Select axis for plotting
axes(hAx(set));
hold on;
line([-100 100],[50 50],'Color',.6*[1 1 1], 'LineStyle','--','LineWidth',.15,'HandleVisibility','off')

% Curve plotting
if markerColor==[1 1 1]
    lineColor='k';
else
    lineColor=markerColor;
end
plot(x, y, markerType, 'MarkerSize', 15, 'Color', lineColor, 'MarkerFaceColor', markerColor, 'MarkerEdgeColor', 'k', 'linewidth', 2.5);
end

function [betaSorted, muSorted, sigmaSorted, contrastsSorted, percentVerticalSorted] = sortParameters(beta, mu, sigma, contrasts, percentCorrect)
% Sort and organize fitting results for output
% Placeholder logic, adjust sorting as needed

betaSorted = beta; % Example sorted beta values
muSorted = mu; % Example sorted mu values
sigmaSorted = sigma; % Example sorted sigma values
contrastsSorted = contrasts; % Example sorted contrasts
percentVerticalSorted = percentCorrect; % Example sorted percent vertical reports

end

% === finalizePlot ===
function finalizePlot(dsCurrentSess, betaSorted, muSorted, sigmaSorted, contrastsSorted, percentVerticalSorted)
% Final adjustments to the plot
axis square;
[~, h] = suplabel([dsCurrentSess.date 'R' dsCurrentSess.run], 't', [.08 .08 .84 .8]);
hold off;
end
