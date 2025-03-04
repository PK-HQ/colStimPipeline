function yPredicted = normMdlSim(x, params, optoStr)
    % Visual + Opto-stim Psychometric Functions
    % Inputs:
    %   x - Stimulus contrast levels
    %   params - Structure containing model parameters
    % Outputs:
    %   yOpto - Correct percentages for congruent and incongruent optostim
    %   yBaseline - Correct percentages for baseline
    
    %% Parameters
    a = params(1); % Weight of iso-tuned visual input
    b = params(2); % Weight of ortho-tuned visual input
    l = params(3); % Weight of iso-tuned visual input norm
    w = params(4); % Weight of ortho-tuned visual input norm
    g0 = params(5); % Normalization constant
    n = params(6); % Spiking exponent
    rmx = params(7); % Max response
    e = params(8); % Relative weight of optostim effects on excitation
    switch optoStr
        case 'baseline'
            o = 0; % Zero opto-stim effective contrast
        case 'opto'
            o = params(9); % Opto-stim effective contrast
    end

    % stimulation setup
    nTrials = 10000; % Trials per level of stimulus contrast
    contrasts = x; % Stimulus contrast levels
    nContrasts = numel(x);

    %% Storage for results
    normresp_cH_vH_oV = zeros(1, nContrasts); 
    normresp_cV_vH_oV = zeros(1, nContrasts); 
    normresp_cH_vH_oH = zeros(1, nContrasts); 
    normresp_cV_vH_oH = zeros(1, nContrasts);
    normresp_cH_vV_oV = zeros(1, nContrasts); 
    normresp_cV_vV_oV = zeros(1, nContrasts); 
    normresp_cH_vV_oH = zeros(1, nContrasts); 
    normresp_cV_vV_oH = zeros(1, nContrasts);
    
    % Calculate results for the current parameter set
    for i = 1:nContrasts
        % Excitatory and normalization signals, naming: column, visual ort, opto ort
        resp_cH_vH_oV = contrasts(i) + b * (1 - e) * o;
        norm_cH_vH_oV = contrasts(i) + w * e * o;
        normresp_cH_vH_oV(i) = rmx * (resp_cH_vH_oV^n) / (norm_cH_vH_oV^n + g0^n);
    
        resp_cV_vH_oV = a * contrasts(i) + (1 - e) * o;
        norm_cV_vH_oV = l * contrasts(i) + e * o;
        normresp_cV_vH_oV(i) = rmx * (resp_cV_vH_oV^n) / (norm_cV_vH_oV^n + g0^n);
    
        resp_cH_vH_oH = contrasts(i) + (1 - e) * o;
        norm_cH_vH_oH = contrasts(i) + e * o;
        normresp_cH_vH_oH(i) = rmx * (resp_cH_vH_oH^n) / (norm_cH_vH_oH^n + g0^n);
    
        resp_cV_vH_oH = a * contrasts(i) + b * (1 - e) * o;
        norm_cV_vH_oH = l * contrasts(i) + w * e * o;
        normresp_cV_vH_oH(i) = rmx * (resp_cV_vH_oH^n) / (norm_cV_vH_oH^n + g0^n);
    
        resp_cH_vV_oV = a * contrasts(i) + b * (1 - e) * o;
        norm_cH_vV_oV = l * contrasts(i) + w * e * o;
        normresp_cH_vV_oV(i) = rmx * (resp_cH_vV_oV^n) / (norm_cH_vV_oV^n + g0^n);
    
        resp_cV_vV_oV = contrasts(i) + (1 - e) * o;
        norm_cV_vV_oV = contrasts(i) + e * o;
        normresp_cV_vV_oV(i) = rmx * (resp_cV_vV_oV^n) / (norm_cV_vV_oV^n + g0^n);
    
        resp_cH_vV_oH = a * contrasts(i) + (1 - e) * o;
        norm_cH_vV_oH = l * contrasts(i) + e * o;
        normresp_cH_vV_oH(i) = rmx * (resp_cH_vV_oH^n) / (norm_cH_vV_oH^n + g0^n);
    
        resp_cV_vV_oH = contrasts(i) + b * (1 - e) * o;
        norm_cV_vV_oH = contrasts(i) + w * e * o;
        normresp_cV_vV_oH(i) = rmx * (resp_cV_vV_oH^n) / (norm_cV_vV_oH^n + g0^n);
    end

    %% Simulate trials and compute performance
    trials = zeros(nTrials, 4);
    for trial = 1:nTrials
        % Generate random input and contrast index
        inp = randi([0 1], 1, 2, 'single'); % Binary inputs for trials
        contrastIdx = randi([1 nContrasts]); % Random contrast index
        trials(trial, 1:2) = inp;
        trials(trial, 3) = contrastIdx;
        visualOrt=inp(1); %0=0, 1=90
        optoOrt=inp(2); %0=0, 1=90
        % Determine response based on input condition
        if visualOrt == 0 && optoOrt == 0
            % H/V column + V visual + V opto
            respH = randn() + normresp_cH_vV_oV(contrastIdx);
            respV = randn() + normresp_cV_vV_oV(contrastIdx);
        elseif visualOrt == 0 && optoOrt == 1
            % H/V column + V visual + H opto
            respH = randn() + normresp_cH_vV_oH(contrastIdx);
            respV = randn() + normresp_cV_vV_oH(contrastIdx);
        elseif visualOrt == 1 && optoOrt == 0
            % H/V column + H visual + V opto
            respH = randn() + normresp_cH_vH_oV(contrastIdx);
            respV = randn() + normresp_cV_vH_oV(contrastIdx);
        elseif visualOrt == 1 && optoOrt == 1
            % H/V column + H visual + H opto
            respH = randn() + normresp_cH_vH_oH(contrastIdx);
            respV = randn() + normresp_cV_vH_oH(contrastIdx);
        end
    
        % Calculate likelihood ratio
        lr = sum(exp(-0.5 * ((respH - [normresp_cH_vH_oV, normresp_cH_vH_oH]).^2 + ...
                             (respV - [normresp_cV_vH_oV, normresp_cV_vH_oH]).^2))) / ...
             sum(exp(-0.5 * ((respH - [normresp_cH_vV_oV, normresp_cH_vV_oH]).^2 + ...
                             (respV - [normresp_cV_vV_oV, normresp_cV_vV_oH]).^2)));
    
        % Determine trial outcome
        trials(trial, 4) = double(lr > 1.0 || (lr == 1.0 && rand() > 0.5));
    end


    %% Initialize output as a structure
    output = struct(...
        'stimulusContrast', zeros(nContrasts, 1), ...
        'effectiveOptoContrast', zeros(nContrasts, 1), ...
        'congruentTrialCount', zeros(nContrasts, 1), ...
        'incongruentHVTrialCount', zeros(nContrasts, 1), ...
        'incongruentVHTrialCount', zeros(nContrasts, 1), ...
        'baselineTrialCount', zeros(nContrasts, 1), ...
        'correctCongruentTrialCount', zeros(nContrasts, 1), ...
        'correctIncongruentHVTrialCount', zeros(nContrasts, 1), ...
        'correctIncongruentVHTrialCount', zeros(nContrasts, 1), ...
        'correctBaselineTrialCount', zeros(nContrasts, 1));

    %% Process trials
    for trial = 1:nTrials
        % Extract variable definitions from trials matrix
        stimCon = trials(trial, 1); % Stimulus condition (congruent/incongruent)
        optoCon = trials(trial, 2); % Optostim condition (horizontal/vertical)
        contrastIdx = trials(trial, 3); % Index for the current contrast level
        lr = trials(trial, 4); % Horizontal-vertical incongruent trial count
    
        % Increment counters based on input conditions
        if stimCon == 1 && optoCon == 1
            % Congruent condition
            output.congruentTrialCount(contrastIdx) = output.congruentTrialCount(contrastIdx) + 1; % Count congruent trials
            if stimCon == lr
                output.correctCongruentTrialCount(contrastIdx) = output.correctCongruentTrialCount(contrastIdx) + 1; % Count correct congruent trials
            end
        elseif stimCon == 1 && optoCon == 0
            % Horizontal-vertical incongruent condition
            output.incongruentHVTrialCount(contrastIdx) = output.incongruentHVTrialCount(contrastIdx) + 1; % Count incongruent trials (horizontal-vertical)
            if stimCon == lr
                output.correctIncongruentHVTrialCount(contrastIdx) = output.correctIncongruentHVTrialCount(contrastIdx) + 1; % Count correct incongruent trials
            end
        elseif stimCon == 0 && optoCon == 1
            % Vertical-horizontal incongruent condition
            output.incongruentVHTrialCount(contrastIdx) = output.incongruentVHTrialCount(contrastIdx) + 1; % Count incongruent trials (vertical-horizontal)
            if stimCon == lr
                output.correctIncongruentVHTrialCount(contrastIdx) = output.correctIncongruentVHTrialCount(contrastIdx) + 1; % Count correct incongruent trials
            end
        elseif stimCon == 0 && optoCon == 0
            % Baseline condition
            output.baselineTrialCount(contrastIdx) = output.baselineTrialCount(contrastIdx) + 1; % Count baseline trials
            if stimCon == lr
                output.correctBaselineTrialCount(contrastIdx) = output.correctBaselineTrialCount(contrastIdx) + 1; % Count correct baseline trials
            end
        end
    end

    %% Compute probabilities
    probCorrectCongruent = (output.correctCongruentTrialCount + output.correctBaselineTrialCount) ./ ...
                           (output.congruentTrialCount + output.baselineTrialCount); % total correct con+bl / total count
    
    probCorrectIncongruent = (output.correctIncongruentHVTrialCount + output.correctIncongruentVHTrialCount) ./ ...
                             (output.incongruentHVTrialCount + output.incongruentVHTrialCount); % total correct incon / total count
    
    probCorrectBaseline = (output.correctCongruentTrialCount + output.correctIncongruentHVTrialCount + ...
                           output.correctIncongruentVHTrialCount + output.correctBaselineTrialCount) ./ ...
                          (output.congruentTrialCount + output.incongruentHVTrialCount + ...
                           output.incongruentVHTrialCount + output.baselineTrialCount);

    % Define indices
    zeroIdx = find(contrasts == 0); % Index for zero contrast
    inconIdx = find(contrasts < 0); % Indices for incongruent contrasts
    conIdx = find(contrasts > 0); % Indices for congruent contrasts

    % Baseline and optostim outputs
    switch optoStr
        case 'baseline'
            yPredicted = [100 - (probCorrectBaseline(inconIdx) * 100); ...
                          rmnan(mean([100 - (probCorrectBaseline(zeroIdx) * 100), probCorrectBaseline(zeroIdx) * 100]))'; ...
                          probCorrectBaseline(conIdx) * 100]';
            %fprintf('%.0f ',contrasts)
            %fprintf('\nBaseline Lengths: %.0f into %.0f\n', numel(zeroIdx), numel(rmnan(mean([100 - (probCorrectBaseline(zeroIdx) * 100), probCorrectBaseline(zeroIdx) * 100]))'));
        case 'opto'
            yPredicted = [100 - (probCorrectIncongruent(inconIdx) * 100); ...
                          rmnan(mean([100 - (probCorrectIncongruent(zeroIdx) * 100), probCorrectCongruent(zeroIdx) * 100]))'; ...
                          probCorrectCongruent(conIdx) * 100]';
            %fprintf('%.0f ',contrasts)
            %fprintf('\nOpto Lengths: %.0f into %.0f\n', numel(zeroIdx), numel(rmnan(mean([100 - (probCorrectIncongruent(zeroIdx) * 100), probCorrectCongruent(zeroIdx) * 100]))'));
    end

end