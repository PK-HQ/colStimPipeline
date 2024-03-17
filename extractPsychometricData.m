function [gaborContrasts, percentageCorrect] = extractPsychometricData(runData, plotFlags)
%% Plots psychometric curves, thresholds, biases across sessions per monkey
% Determine plot flag
plotFlag = ''; if plotFlags(1) == 1; plotFlag = 'plot'; end

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
        initializeVariables(condType, gaborOrt, opto0, opto90, optoBL, visualOptoTrials, optoControl, successData, completedData, gaborCon);
end
end

function [gaborContrasts, percentageCorrect] = initializeVariables(stimType, gaborOrt, opto0, opto90, optoBL, visualControl, optoControl, successData, completedData, gaborCon)
% Function to initialize variables for data extraction and plotting
negative0 = -1 .* (stimType>0 & gaborOrt == 0);
tmp = negative0' .* -1;
negative0(negative0 == 0) = 1;
gaborCon = gaborCon .* negative0;
vertReportData = successData .* negative0' + completedData .* tmp;
percentVertical = (vertReportData ./ completedData)' .* 100;




gaborContrasts.congruent(1,:)=gaborCon(negative0<0);
percentageCorrect.congruent(1,:)=[100-percentVertical(negative0<0)];
gaborContrasts.incongruent(1,:)=gaborCon(negative0>0);
percentageCorrect.incongruent(1,:)=percentVertical(negative0>0);




% Sorting conditions
%[conds0, conds90, condsOptoControl, condsOpto0, condsOpto90] = sortConditions(stimType, gaborOrt, opto0, opto90, visualControl, optoControl);

% Split into vis 0 or vis 90 

% Split into opto 0 or opto 90

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