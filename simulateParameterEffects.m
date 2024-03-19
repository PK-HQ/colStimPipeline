function simulateParameterEffects(paramsStruct, hAx, modelComplexity)
    %% Initialize Parameters
    visBL = 0;
    alphaCon = 0.1;
    alphaIncon = 0;
    optoH = 0;
    optoV = 0;
    optoBL = 0;
    betaCon = 0;
    betaIncon = 0;
    c50 = .5;
    nPower = 2;
    
    %% Parameters: Visual stimulus
    visV = nlinspace(0, 100, 100, 'nonlinear'); %-100:1:0;
    visH = fliplr(-visV);
    contrasts = [visH visV];

    %% Update parameters with single values from struct
    paramNames = fieldnames(paramsStruct);
    for i = 1:length(paramNames)
        if length(paramsStruct.(paramNames{i})) == 1
            eval([paramNames{i} ' = ' num2str(paramsStruct.(paramNames{i})) ';']);
        end
    end

    %% Initialize figure
    hold on;

    %% Loop through parameters to find the one with a range and plot
    for i = 1:length(paramNames)
        paramName = paramNames{i};
        paramRange = paramsStruct.(paramName);

        % Check if the current parameter has a range
        if length(paramRange) > 1
            mainParamName=paramName;
            mainParamRange=paramRange;

            exclusionBuffer=5; exclusionBuffers=exclusionBuffer*2 + 1;
            colormapCon = (slanCM('sunburst', numel(paramRange)+exclusionBuffers)); 
            colormapCon=[0 0 0; colormapCon([exclusionBuffer+1:end-exclusionBuffer],:)];%plasma
            colormapIncon = (slanCM('freeze', numel(paramRange)+exclusionBuffers));  
            colormapIncon=[0 0 0; colormapIncon([exclusionBuffer+1:end-exclusionBuffer],:)];%plasma

            for j = 1:length(paramRange)
                % Update the varying parameter
                eval([paramName ' = paramRange(j);']);
                prctVertical = zeros(1, numel(contrasts));

                %% Predicted vertical reports
                for nContrast = 1:numel(contrasts)
                    contrast = contrasts(nContrast);
                    contrastH = max(0, -contrast);
                    contrastV = max(0, contrast);

                    % Calculate responses
                    switch modelComplexity
                        case 'full'
                            % H-columns
                            visualH=(alphaCon * contrastH) + (alphaIncon * contrastV) + visBL;
                            optoH=(betaCon * (optoH + optoBL)) + (betaIncon * (optoV + optoBL));
                            respH =  visualH + optoH;
                            % V-columns
                            visualV=(alphaIncon * contrastH) + (alphaCon * contrastV) + visBL;
                            optoV=(betaIncon * (optoH + optoBL)) + (betaCon * (optoV + optoBL));
                            respV = visualV + optoV;
                        case 'simple'
                            % H-columns
                            visualH=(alphaCon * contrastH);
                            optoH=(betaCon * (optoH + optoBL));
                            respH =  visualH + optoH;
                            % V-columns
                            visualV=(alphaCon * contrastV);
                            optoV=(betaCon * (optoV + optoBL));
                            respV = visualV + optoV;
                    end
                    % Feed through Naka Rushton
                    RespHmean = ComputeNakaRushton([c50, nPower], respH);
                    RespVmean = ComputeNakaRushton([c50, nPower], respV);

                    deltaResp = RespVmean - RespHmean;
                    scaling = sqrt(RespHmean + RespVmean);
                    respTrial = 2 * deltaResp / scaling;
                    prctVertical(nContrast) = normcdf(respTrial);
                end

                % Plotting with the specified color for current parameter range
                
                % Now, split plotting based on the sign of contrasts
                negativeIndices = contrasts <= 0;
                positiveIndices = contrasts >= 0;
                
                % Plot on hAx(1) for negative contrasts
                axes(hAx(1));
                hold on;
                if paramRange(j)==0
                    lineName='Baseline';
                elseif paramRange(j)>0 && paramRange(j)<.1
                    lineName=num2str(paramRange(j),1);
                elseif paramRange(j)>0 && paramRange(j)>.1
                    lineName=num2str(paramRange(j),'%.1f');
                end
                plot(contrasts(positiveIndices), prctVertical(positiveIndices)*100, 'LineWidth', 2, 'Color', colormapCon(j, :), 'DisplayName', lineName);
                
                axes(hAx(2));
                hold on;
                plot(-1*contrasts(negativeIndices), 100-prctVertical(negativeIndices)*100, 'LineWidth', 2, 'Color', colormapIncon(j, :), 'DisplayName', lineName);
                
            end
            break; % Exit after plotting the parameter with a range
        end
    end

    %% Finalize plot
    [paramNameSymbol, paramNameMini, paramNameFull, unitStr]=decodeParamName(mainParamName);
    suplabel({paramNameFull,sprintf(['%s: %.1f - %.1f' unitStr],...
        paramNameSymbol,min(mainParamRange),max(mainParamRange))},'t', [.1 .1 .77 .77]);

    for ax=1:2
        axes(hAx(ax));
        xMax=40; yMax=100;
        xlim([0 xMax]); xticks([0:10:xMax])
        ylim([0 yMax]); yticks([0:10:yMax])

        if ax==1
            title('Vertical gabor')
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



%% Subfunctions
function [paramNameSymbol, paramNameMini, paramNameFull, unitStr]=decodeParamName(paramName)
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
end
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

