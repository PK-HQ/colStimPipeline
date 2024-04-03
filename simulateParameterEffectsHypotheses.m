function simulateParameterEffectsHypotheses(inputParams, hAx, modelComplexity)
    %% Initialize Parameters
    [baselineParams,modelParams]=initParamsStruct(inputParams);
    
    %% Parameters: Visual stimulus
    visV = nlinspace(0, 100, 100, 'nonlinear'); 
    visH = fliplr(-visV);
    contrasts = [visH visV];

    %% Update parameters with single values from struct
    paramNames = fieldnames(modelParams);
    for i = 1:length(paramNames)
        if length(modelParams.(paramNames{i})) == 1
            eval([paramNames{i} ' = ' num2str(modelParams.(paramNames{i})) ';']);
        end
    end

    %% Initialize figure
    hold on;

    %% Loop through parameters to find the one with a range and plot
    for i = 1:length(paramNames)
        paramName = paramNames{i};
        paramRange = modelParams.(paramName);

        % Check if the current parameter has a range
        if length(paramRange) > 1
            mainParamName=paramName;
            mainParamRange=paramRange;

            exclusionBuffer=2; exclusionBuffers=round(exclusionBuffer*1) + 2;
            colormapCon = (slanCM('sunburst', numel(paramRange)+exclusionBuffers)); 
            colormapCon=[colormapCon([exclusionBuffer+1:end-exclusionBuffer],:)];%plasma
            colormapIncon = (slanCM('freeze', numel(paramRange)+exclusionBuffers));  
            colormapIncon=[colormapIncon([exclusionBuffer+1:end-exclusionBuffer],:)];%plasma

            for j = 1:length(paramRange)
                % Update the varying parameter
                eval([paramName ' = paramRange(j);']);
                prctVertical = zeros(1, numel(contrasts));

                %% Predicted vertical reports
                prctVertical = calculateResponses(modelParams, j, contrasts, modelComplexity);
                prctVerticalBaseline = calculateResponses(baselineParams, j, contrasts, modelComplexity);
                % Now, split plotting based on the sign of contrasts
                negativeIndices = contrasts <= 0;
                positiveIndices = contrasts >= 0;
                
                % Plot on hAx(1) for negative contrasts
                if paramRange(j)>0 && paramRange(j)<.1
                    lineName=num2str(paramRange(j),1);
                elseif paramRange(j)>0 && paramRange(j)>.1
                    lineName=num2str(paramRange(j),'%.1f');
                end
                hold on;
                plot(contrasts(positiveIndices), prctVertical(positiveIndices), 'LineWidth', 2, 'Color', colormapCon(j, :), 'DisplayName', lineName);
                hold on;
                plot(-1*contrasts(negativeIndices), 100-prctVertical(negativeIndices), 'LineWidth', 2, 'Color', colormapIncon(j, :), 'HandleVisibility','off');
            end
            break; % Exit after plotting the parameter with a range
        end
    end
    plot(contrasts(positiveIndices), prctVerticalBaseline(positiveIndices), 'LineWidth', 2, 'Color', 'k', 'DisplayName', 'Baseline');
    %plot(-1*contrasts(negativeIndices), 100-prctVerticalBaseline(negativeIndices), 'LineWidth', 2, 'Color', 'k', 'DisplayName', 'Baseline');

    %% Finalize plot
    [paramNameSymbol, paramNameMini, paramNameFull, unitStr, hypothesisNameFull]=decodeParamName(mainParamName, modelParams.hypothesis);
    suplabel({hypothesisNameFull,sprintf(['%s: %.1f - %.1f' unitStr],...
        paramNameSymbol,min(mainParamRange),max(mainParamRange))},'t', [.1 .1 .77 .77]);

    for ax=1:length(hAx)
        if length(hAx)>1
            axes(hAx(ax));
        end
        xMax=40; yMax=100;
        xlim([0 xMax]); xticks([0:10:xMax])
        ylim([0 yMax]); yticks([0:10:yMax])

        if ax==1
            %title('Vertical gabor')
            ylabel('Correct reports (%)'); xlabel('Absolute gabor contrast');
            
            % Parameters and their values for the textbox
            [paramsTextA, paramsTextB] = generateParamsText(alphaCon, alphaIncon, visBL, betaCon, betaIncon, optoH, optoV, optoBL, c50, nPower, mainParamName);

            % Adding the textbox annotation outside the plot, under the x-axis
            text(1.35+.04,.8,paramsTextA,'Units','normalized', 'VerticalAlignment','top')
            text(1.35+.3,.8,paramsTextB,'Units','normalized', 'VerticalAlignment','top')
        elseif ax==2
            title('Horizontal gabor')
        end
        
        hold off;
        %% Guide lines
        line([-100 100],[50 50],'Color',.6*[1 1 1], 'LineStyle','--','LineWidth',.15,'HandleVisibility','off')
        line([0 0],[-100 100],'Color',.6*[1 1 1], 'LineStyle','--','LineWidth',.15,'HandleVisibility','off')
        
        legend('show', 'Location', 'southoutside','NumColumns',2,'FontSize',12);
        legendtitle(paramNameMini)
        
        upFontSize(24,0.005)
        axis square
    end
end

%% === Subfunctions ===
function [baselineParams, modelParams] = initParamsStruct(inputParams)
    % Initialize default parameter values and add hypothesis
    baselineParams = struct(...
        'alphaCon', 0.1, ...
        'alphaIncon', 0, ...
        'betaCon', 1, ...
        'betaIncon', 0.001, ...
        'visBL', 0, ...
        'optoH', 0, ...
        'optoV', 0, ...
        'optoBL', 0.00001, ...
        'c50', 1, ...
        'nPower', 2, ...
        'rMax', 1, ...
        'visualBias', 0, ...
        'integrationFactor', 0, ...
        'dominanceFactor', 0, ...
        'rtAdjustment', 0, ...
        'hypothesis', '');

    % Copy baselineParams to create modelParams and update with inputParams
    modelParams = baselineParams;
    modelParams.hypothesis = inputParams.hypothesis;  % Set hypothesis

    paramNames = fieldnames(baselineParams);
    for i = 1:length(paramNames)
        paramName = paramNames{i};
        if isfield(inputParams, paramName) && ~isempty(inputParams.(paramName))
            modelParams.(paramName) = inputParams.(paramName);
        end
    end
    
    % Adjust modelParams based on the specified hypothesis
    switch modelParams.hypothesis
        case 'snr'
            modelParams.c50 = modelParams.c50 * 0.8;
            modelParams.nPower = modelParams.nPower + 0.2;
            modelParams.rtAdjustment = -0.1;
        case 'decisionbias'
            modelParams.decisionShift = -0.05;
            modelParams.rtAdjustment = 0.05;
        case 'optopercept'
            modelParams.betaCon = modelParams.betaCon * 1.01;
            modelParams.rtAdjustment = -0.15;
        case 'plaid'
            modelParams.integrationFactor = 1;
            modelParams.rtAdjustment = -0.1;
        case 'winner'
            modelParams.dominanceFactor = 2.5;
            modelParams.rtAdjustment = -0.2;
    end
end


% Simulate columnar responses
function prctVertical = calculateResponses(paramsStruct, paramNum, contrasts, modelComplexity)
    prctVertical = zeros(1, numel(contrasts));
    
    paramsStructIter = updateModelParam(paramsStruct, paramNum);
    % Unpack parameters from paramsStruct for clarity
    alphaCon = paramsStructIter.alphaCon;
    alphaIncon = paramsStructIter.alphaIncon;
    betaCon = paramsStructIter.betaCon;
    betaIncon = paramsStructIter.betaIncon;
    visBL = paramsStructIter.visBL;
    optoH = paramsStructIter.optoH;
    optoV = paramsStructIter.optoV;
    optoBL = paramsStructIter.optoBL;
    c50 = paramsStructIter.c50;
    nPower = paramsStructIter.nPower;
    
    % Additional parameters
    visualBias = paramsStructIter.visualBias;
    integrationFactor = paramsStructIter.integrationFactor;
    dominanceFactor = paramsStructIter.dominanceFactor;
    
    for nContrast = 1:numel(contrasts)
        contrast = contrasts(nContrast);
        contrastH = max(0, -contrast);
        contrastV = max(0, contrast);

        switch modelComplexity
            case 'full'
                % Incorporating visualBias and integrationFactor into full model
                visualH = (alphaCon * contrastH + visualBias) + (alphaIncon * contrastV);
                visualV = (alphaIncon * contrastH) + (alphaCon * contrastV + visualBias);
                
                % Adjusting optoH and optoV based on integrationFactor and dominanceFactor
                optoH = betaCon * (optoH + optoBL + dominanceFactor) + betaIncon * (optoV + optoBL + integrationFactor);
                optoV = betaIncon * (optoH + optoBL + integrationFactor) + betaCon * (optoV + optoBL + dominanceFactor);

            case 'simple'
                % Adjusting simple model with dominanceFactor only for simplicity
                optoH = (betaCon * (optoH + dominanceFactor + optoBL));
                optoV = (betaCon * (optoV + dominanceFactor + optoBL));

                visualH = (alphaCon * contrastH);
                visualV = (alphaCon * contrastV);
        end

        % Calculate responses
        respH = visualH + optoH;
        respV = visualV + optoV;

        % Feed through Naka Rushton
        RespHmean = ComputeNakaRushton([c50, nPower], respH);
        RespVmean = ComputeNakaRushton([c50, nPower], respV);

        deltaResp = RespVmean - RespHmean;
        scaling = sqrt(RespHmean + RespVmean);
        respTrial = 2 * deltaResp / scaling;
        prctVertical(nContrast) = normcdf(respTrial, 0, 1) * 100;  % Ensure percent conversion
    end
end

% Update param for model iter
function paramsStructIter = updateModelParam(paramsStruct, paramNum)
    
    paramsStructIter=paramsStruct;
    
    % Get fieldnames of paramsStruct
    paramNames = fieldnames(paramsStruct);

    % Iterate through each field
    for i = 1:length(paramNames)
        % Get the field name
        fieldName = paramNames{i};

        % Get the value of the field
        fieldValue = paramsStruct.(fieldName);

        % Check if fieldValue has length greater than 1
        if length(fieldValue) > 1 && ~ischar(fieldValue)
            % Apply an index to fieldValue
            paramsStructIter.(fieldName)=paramsStruct.(fieldName)(paramNum);
        end
    end

end

% Param details
function [paramNameSymbol, paramNameMini, paramNameFull, unitStr, hypothesisNameFull]=decodeParamName(paramName, hypothesis)
unitStr=' a.u.';
switch paramName
    case 'alphaCon'
        paramNameSymbol='\alpha_{iso}';
        paramNameMini='Pref. column weight';
        paramNameFull=['Visual-weights of iso-tuned columns'];
    case 'alphaIncon'
        paramNameSymbol='\alpha_{ortho}';
        paramNameMini='Non-Pref. column weight';
        paramNameFull=['Visual-weights of ortho-tuned columns'];
    case 'betaCon'
        paramNameSymbol='\beta_{Con}';
        paramNameMini='Con. opto weight';
        paramNameFull=['Opto-weights of congruent optostim'];
    case 'betaIncon'
        paramNameSymbol='\beta_{Incon}';
        paramNameMini='Incon. opto weight';
        paramNameFull=['Opto-weights of incongruent optostim'];
    case 'optoH'
        paramNameSymbol='opto_{H}';
        paramNameMini='opto_{H} strength';
        paramNameFull=['Strength of opto_{H}'];
        unitStr='%%';
    case 'optoV'
        paramNameSymbol='opto_{V}';
        paramNameMini='opto_{V} strength';
        paramNameFull=['Strength of opto_{V}'];
        unitStr='%%';
    case 'optoBL'
        paramNameSymbol='opto_{BL}';
        paramNameMini='opto_{BL} strength';
        paramNameFull=['Strength of opto_{BL}'];
        unitStr='%%';
    case 'c50'
        paramNameSymbol='C_{50}';
        paramNameMini='N-R C_{50}';
        paramNameFull=['Naka-Rushton C_{50}'];
    case 'nPower'
        paramNameSymbol='n';
        paramNameMini='N-R n';
        paramNameFull=['Naka-Rushton exponent'];
end
switch hypothesis
    case ''
        hypothesisNameFull='Base model';
    case 'snr'
        hypothesisNameFull='\Delta SNR (no opto-percept)';
    case 'decisionbias'
        hypothesisNameFull='\Delta Decision bias(no opto-percept)';
    case 'optopercept'
        hypothesisNameFull='? (Single opto-percept)';
    case 'plaid'
        hypothesisNameFull='Integration (Plaid opto-percept)';
    case 'winner'
        hypothesisNameFull='Competition/Winner-takes-all (Single opto-percept)';
end
end

% Param text
function [paramsTextA, paramsTextB] = generateParamsText(alphaCon, alphaIncon, visBL, betaCon, betaIncon, optoH, optoV, optoBL, c50, nPower, mainParamName)
    % Adaptation for paramsTextA and paramsTextB to include mainParamName check
    
    % Ensure each variable for paramsTextA is processed
    alphaConText = formatVar(alphaCon, 'alphaCon', mainParamName);
    alphaInconText = formatVar(alphaIncon, 'alphaIncon', mainParamName);
    visBLText = formatVar(visBL, 'visBL', mainParamName);
    betaConText = formatVar(betaCon, 'betaCon', mainParamName);
    betaInconText = formatVar(betaIncon, 'betaIncon', mainParamName);
    optoHText = formatVar(optoH, 'optoH', mainParamName);
    optoVText = formatVar(optoV, 'optoV', mainParamName);
    optoBLText = formatVar(optoBL, 'optoBL', mainParamName);
    
    % Construct paramsTextA
    paramsTextA = sprintf(...
                 ['{\\bfVisual}\n',...
                 '\\alpha_{iso}: %s\n', ...
                  '\\alpha_{ortho}: %s\n', ...
                  'Vis_{BL}: %s\n\n', ...
                  '{\\bfOpto}\n',...
                  '\\beta_{Con}: %s\n', ...
                  '\\beta_{Incon}: %s\n', ...
                  'Opto_{H}: %s\n', ...
                  'Opto_{V}: %s\n', ...
                  'Opto_{BL}: %s'], ...
                  alphaConText, alphaInconText, visBLText,...
                  betaConText, betaInconText, optoHText, optoVText, optoBLText);

    % Ensure each variable for paramsTextB is processed
    c50Text = formatVar(c50, 'c50', mainParamName);
    nPowerText = formatVar(nPower, 'nPower', mainParamName);
    
    % Construct paramsTextB
    paramsTextB = sprintf(...
                 ['{\\bfN-R}\n',...
                  'C_{50}: %s\n', ...
                  'n: %s'], ...
                  c50Text, nPowerText);

    % Subfunction defined within the same file, after the main function
    function output = formatVar(var, varName, mainParamName)
        if strcmp(varName, mainParamName)
            output = 'varies';
        else
            % If the variable name does not match mainParamName, proceed with unique check
            uniqueVar = unique(var);
            if numel(uniqueVar) == 1
                % If there's exactly one unique value, convert it to string
                output = num2str(uniqueVar);
            else
                % If there are multiple unique values, return 'varies'
                output = 'varies';
            end
        end
    end
end

