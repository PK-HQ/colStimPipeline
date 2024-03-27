function simulateParameterEffectsMdl(paramsStruct, hAx)
    models = {'plaid', 'winnerTakeAll', 'optoBiasWithPerception', 'optoBiasWithoutPerception', 'SNRChange'};
    
    for modelIdx = 1:length(models)
        figure('Name',models{modelIdx});
        paramsStruct.initparams = initializeModelParams(); % Initialize parameters
        paramsStruct.modelType = models{modelIdx};
        
        resultsVertical = []; % Store results for vertical Gabor simulation
        resultsHorizontal = []; % Store results for horizontal Gabor simulation
        
        % Define contrast ranges for plotting
        contrastsVertical = linspace(0, 100, 100);
        contrastsHorizontal = linspace(0, 100, 100); % Absolute values for -100:0
        
        for optoCond = ["baseline", "congruent", "incongruent"]
            % Get predicted params for each of the 5 models
            adjustedParams = adjustParamsForModel(paramsStruct);
            % Set optostim condition
            adjustedParams.condition = optoCond;

            % Vertical Gabor simulation
            prctVertical = simulateModel(adjustedParams, contrastsVertical);
            resultsVertical(end+1, :) = prctVertical; % Append results
            
            % Horizontal Gabor simulation (use absolute values for simulation)
            prctHorizontal = simulateModel(adjustedParams, -contrastsHorizontal); % Simulate with negative contrasts
            resultsHorizontal(end+1, :) = prctHorizontal; % Append results
        end
        
        % Plot results for vertical and horizontal Gabors, and their mean
        axes(hAx(1)); % First subplot for vertical Gabor
        plotResults(contrastsVertical, resultsVertical);
        
        axes(hAx(2)); % Second subplot for horizontal Gabor
        plotResults(contrastsHorizontal, resultsHorizontal); % Use adjusted contrasts
        
        axes(hAx(3)); % Third subplot for the mean
        resultsMean = (resultsVertical + resultsHorizontal) / 2;
        plotResults(contrastsVertical, resultsMean); % Use vertical contrasts for x-axis
        
        % Additional formatting for the figure
        suplabel(['Model: ' models{modelIdx}],'t',[.1 .1 .77 .77]);
        adjustSubplots(hAx);
    end
end

%% Subfunctions
% Initialize
function paramsStruct = initializeModelParams()
    % Initialize model parameters including visualTuning and optoTuning
    paramsStruct = struct(...
        'alphaCon', 0.1, ...
        'alphaIncon', 0.001, ...
        'optoH', 0, ...
        'optoV', 0, ...
        'optoBL', 0.0001, ...
        'betaCon', 1, ...
        'betaIncon', 0.01, ...
        'c50', 50, ...
        'nPower', 2, ...
        'rMax', 100, ...
        'visBL', 0, ...
        'modelType', '', ...
        'visualTuning', 0, ... % Default value, adjust as needed
        'optoTuning', 0); % Default value, adjust as needed
end

% Select params for given model 
function adjustedParams = adjustParamsForModel(modelParams)
    % Adjust parameters based on the specified model type
    adjustedParams = modelParams.initparams; % Start with the original parameters
    
    switch modelParams.modelType
        case 'plaid'
            % Plaid: Both visual and opto inputs contribute equally, enhancing each other.
            adjustedParams.alphaCon = adjustedParams.alphaCon + 0.05; % Slightly enhance contribution of congruent visual input.
            adjustedParams.betaCon = adjustedParams.betaCon + 0.05; % Enhance contribution of congruent opto input.
            
        case 'winnerTakeAll'
            % Winner-takes-all: Opto input dominates over visual input.
            adjustedParams.alphaCon = adjustedParams.alphaCon * 0.5; % Reduce effectiveness of visual input.
            adjustedParams.betaCon = adjustedParams.betaCon + 0.2; % Significantly enhance contribution of congruent opto input.
            
        case 'optoBiasWithPerception'
            % Opto-bias with perception: Opto input biases the response but also adds to the perception.
            adjustedParams.optoBL = adjustedParams.optoBL + 0.1; % Increase baseline opto contribution, implying a bias.
            adjustedParams.alphaCon = adjustedParams.alphaCon + 0.05; % Opto input slightly enhances visual input's effect.
            
        case 'optoBiasWithoutPerception'
            % Opto-bias without perception: Opto input biases the response without contributing to the actual perception.
            adjustedParams.optoBL = adjustedParams.optoBL + 0.15; % Increase baseline opto contribution more than with perception.
            % This model assumes that while there's an opto bias, it doesn't enhance visual input's effectiveness.
            
        case 'SNRChange'
            % SNR Change: Improvement in signal-to-noise ratio, making it easier to detect the signal.
            adjustedParams.c50 = adjustedParams.c50 - 0.1; % Decrease C50, making the system more sensitive at lower contrasts.
            adjustedParams.nPower = adjustedParams.nPower + 0.5; % Increase the steepness of the response curve.
    end
end

% Simulate psychometric results
function prctVertical = simulateModel(adjustedParams, contrasts)
    prctVertical = zeros(1, length(contrasts));
    for nContrast = 1:length(contrasts)
        contrast = contrasts(nContrast);
        
        % Simplify by considering contrast positive for visualization
        contrast = abs(contrast);
        
        % Determine effects based on condition
        switch adjustedParams.condition
            case 'baseline'
                % Baseline condition: No opto stimulation
                optoEffectH = 0;
                optoEffectV = 0;
            case 'congruent'
                % Congruent condition: Opto supports visual stimuli
                optoEffectH = adjustedParams.betaCon * adjustedParams.optoH;
                optoEffectV = adjustedParams.betaCon * adjustedParams.optoV;
            case 'incongruent'
                % Incongruent condition: Opto opposes visual stimuli
                optoEffectH = adjustedParams.betaIncon * adjustedParams.optoH;
                optoEffectV = adjustedParams.betaIncon * adjustedParams.optoV;
            otherwise
                optoEffectH = 0;
                optoEffectV = 0;
        end

        % Apply visual effect based on visualTuning (assuming 0 for vertical, 90 for horizontal)
        visualEffect = adjustedParams.alphaCon * contrast;

        % Compute responses considering the condition effects
        respH = visualEffect + optoEffectH; % For horizontal Gabor simulation
        respV = visualEffect + optoEffectV; % For vertical Gabor simulation

        % Feed through Naka Rushton
        RespHmean = ComputeNakaRushton([adjustedParams.c50, adjustedParams.nPower, adjustedParams.rMax], respH);
        RespVmean = ComputeNakaRushton([adjustedParams.c50, adjustedParams.nPower, adjustedParams.rMax], respV);

        % Calculate percent correct
        deltaResp = RespVmean - RespHmean;
        scaling = sqrt(RespHmean^2 + RespVmean^2); % Corrected to square the terms before sqrt
        respTrial = deltaResp / scaling;
        prctVertical(nContrast) = normcdf(respTrial, 0, 1) * 100; % Convert to percentage
    end
end

% Plot fit
function plotResults(contrasts, results)
    hold on;
    colors = {[0, 0, 0], [220, 20, 60]/255, [0, 0, 128]/255}; % Black, Crimson Red, Navy Blue
    
    % Plot each condition in the appropriate color
    for i = 1:size(results, 1)
        plot(contrasts, results(i, :), 'LineWidth', 2, 'Color', colors{i});
    end

    xlim([0 100]);
    ylim([0 100]);
    hold off;
end

function adjustSubplots(hAx)
    % Adjust subplot settings uniformly, if needed
    for ax = 1:length(hAx)
        axes(hAx(ax));
        if ax==1
            xlabel('Absolute Gabor Contrast (%)');
            ylabel('Percent Correct (%)');
            legend({'Baseline', 'Congruent', 'Incongruent'}, 'Location', 'southeast');
        end
        upFontSize(24,0.01)
        axis square
    end
end