%% Define chamber and model type
chambers={'R','L'};
monkeyName='Chip';
[mainPath, datastruct]=setupEnv('users/PK/colStimPipeline/exptListBiasingFull.m');

%% Initial fit to find median params
% Fitting
objFunc='MLE';
modelTypes={'base','beta','deltax'};
constrainedParamStr='';
saveFlag=1;
for chamberID=[1 2]
    chamberWanted=chambers{chamberID};
    switch chamberWanted
        case 'L'
            load('Y:\Chip\Meta\summary\statisticsL.mat');
        case 'R'
            load('Y:\Chip\Meta\summary\statisticsR.mat'); 
    end
    for modelID=[2 3]
        modelType=modelTypes{modelID};
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
        mdl=fitNakaRushtonMLE3(xBlocks, yBlocks, modelType); % Beta/exp/C50
        mdlStruct.([chamberWanted, modelType]) = mdl;

        %% Plot param distributions
        %{
        fieldNames = {'beta', 'exponent', 'c50', 'xDeltaBL', 'xDeltaOpto',  'AICc'};
        titleStrs = {'\beta', 'Exponent', 'C_{50}', '\Deltax-BL',  '\Deltax-Opto', 'AIC_{c}'};
        paramFitted=mdl.fittedParams(:,paramIdx,1:2);
        fieldName=fieldNames{paramIdx};
        titleStr=titleStrs{paramIdx};
        savefilename=['psychometrics\' chamberWanted '-chamber/' modelType '/' 'paramdist-' chamberWanted '-' modelType constrainedParamStr];
        plotParamDistributions(paramFitted,fieldName,titleStr, modelType, monkeyName, savefilename, saveFlag)
        %}
        
        %% Plot fits
        xFit=sort([-nlinspace(0, 100, 100, 'nonlinear') nlinspace(0, 100, 100, 'nonlinear')]);
        savefilename=['psychometrics\' chamberWanted '-chamber/' modelType '/' 'psychfit-' chamberWanted '-' modelType];
        plotNakaRushtonFit3(behavioralData, bitmapData, datastruct, analysisBlockID, mdl, mdl.fittedParams(:,:,1), ...
            xFit, monkeyName, saveFlag, savefilename)
    end
    
    fieldName='AICc';
    titleStr='\beta vs \Deltax';
    savefilename=['psychometrics\' chamberWanted '-chamber/' 'paramdist-' chamberWanted '-' modelTypes{2} modelTypes{3}];
    modelType={'\beta','\Deltax'};
    plotParamDistributions(mdlStruct.beta,mdlStruct.deltax,fieldName,titleStr, modelType, monkeyName, savefilename, saveFlag)

end
    [pValL, ~, ~] = signrank(mdlStruct.Lbeta.fittedParams(:,end),mdlStruct.Ldeltax.fittedParams(:,end), 'tail', 'right')

    [pValR, ~, ~] = signrank(mdlStruct.Rbeta.fittedParams(:,end),mdlStruct.Rdeltax.fittedParams(:,end), 'tail', 'right')

