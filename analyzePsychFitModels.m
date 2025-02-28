%% Define chamber and model type
chambers={'R','L'};
monkeyName='Chip';
[mainPath, datastruct]=setupEnv('users/PK/colStimPipeline/exptListBiasingFull.m');

%% Initial fit to find median params
% Fitting
objFunc='MLE';
modelTypes={'base','beta','deltax'};
fieldName='AICc';
titleStr='/beta vs /Deltax';
modelTypeStrs={'/beta','/Deltax'};
constrainedParamStr='';
saveFlag=0;
for chamberID=[2 1]
    chamberWanted=chambers{chamberID};
    switch chamberWanted
        case 'L'
            load([mainPath 'Chip/Meta/summary/statisticsL.mat']);
        case 'R'
            load([mainPath 'Chip/Meta/summary/statisticsR.mat']); 
    end
    for modelID=[2 3]
        modelTypeStr=modelTypes{modelID};
        switch modelID
            case 1
                paramIdx=6;
            case 2
                paramIdx=6;
        end
        %% Fit psychometric data
        lines=1:3;
        plots=3;
        contrasts=1:size(behavioralData.gaborContrasts,2);
        blocks=1:size(behavioralData.gaborContrasts,3);
        xBlocks=squeeze(behavioralData.gaborContrasts(lines,contrasts,blocks));
        yBlocks=squeeze(behavioralData.percentageCorrect(lines,contrasts,blocks));
        constrainedParams=[nan nan nan];
        mdl=fitNakaRushtonMLE3(xBlocks, yBlocks, modelTypeStr); % Beta/exp/C50
        mdlStruct.([chamberWanted, modelTypeStr]) = mdl;

        %% Plot param distributions
        %{
        fieldNames = {'beta', 'exponent', 'c50', 'xDeltaBL', 'xDeltaOpto',  'AICc'};
        titleStrs = {'/beta', 'Exponent', 'C_{50}', '/Deltax-BL',  '/Deltax-Opto', 'AIC_{c}'};
        paramFitted=mdl.fittedParams(:,paramIdx,1:2);
        fieldName=fieldNames{paramIdx};
        titleStr=titleStrs{paramIdx};
        savefilename=['psychometrics/' chamberWanted '-chamber/' modelType '/' 'paramdist-' chamberWanted '-' modelType constrainedParamStr];
        plotParamDistributions(paramFitted,fieldName,titleStr, modelType, monkeyName, savefilename, saveFlag)
        %}
        
        %% Plot fits
        xFit=sort([-nlinspace(0, 100, 100, 'nonlinear') nlinspace(0, 100, 100, 'nonlinear')]);
        savefilename=['psychometrics/' chamberWanted '-chamber/' modelTypeStr '/' 'psychfit-' chamberWanted '-' modelTypeStr];
        plotNakaRushtonFit3(behavioralData, bitmapData, datastruct, analysisBlockID, mdl, mdl.fittedParams(:,:,1), ...
            xFit, monkeyName, saveFlag, savefilename)
    end
    
    %% Filter and compare models with non-parametric t-test (per chamber)
    
    savefilename=['psychometrics/' chamberWanted '-chamber/' 'paramdist-' chamberWanted '-' modelTypes{2} modelTypes{3}];
    % get power
    bitmapSPD.(chamberWanted)=squeeze(bitmapData.adjustedSPD_uW)';
    % get blocks with similar power (median +-3 MAD)
    [powerFiltered, powerFilteredMask] = powerFiltering(squeeze(bitmapData.adjustedSPD_uW)');
    AICCdeltax.(chamberWanted)=mdlStruct.([chamberWanted, 'deltax']).fittedParams(:,end).*powerFilteredMask;
    AICCbeta.(chamberWanted)=mdlStruct.([chamberWanted, 'beta']).fittedParams(:,end).*powerFilteredMask;
    plotParamDistributions(AICCbeta.(chamberWanted),AICCdeltax.(chamberWanted),fieldName,titleStr, modelTypeStrs, monkeyName, savefilename, saveFlag)
    CF;
end
save([mainPath '/Chip/Meta/psychometrics/fittingParams' chambers{1} chambers{2} '-' modelTypes{2} modelTypes{3} '.mat'],'mdlStruct','bitmapSPD','AICCdeltax','AICCbeta')

%% Filter and compare models with non-parametric t-test (combined chamber)
savefilename=['psychometrics/paramdist-' chambers{1} chambers{2}  '-' modelTypes{2} modelTypes{3} constrainedParamStr];
plotParamDistributions([AICCbeta.('L');AICCbeta.('R')],[AICCdeltax.('L');AICCdeltax.('R')],fieldName,titleStr, modelTypeStrs, monkeyName, savefilename, saveFlag)

ranksum(mdlStruct.Rbeta.fittedParams(:,1),mdlStruct.Lbeta.fittedParams(:,1))
ranksum(mdlStruct.Rdeltax.fittedParams(:,4),mdlStruct.Ldeltax.fittedParams(:,4))
