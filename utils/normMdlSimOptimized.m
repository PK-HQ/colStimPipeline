function yPredicted = normMdlSimOptimized(x, params, optoStr, options)
    % Simplified version of normMdlSimOptimized to match original v_o_psychometric_3.m
    % Uses uniform trial allocation and simple random number generation
    % 
    % Inputs:
    %   x - Stimulus contrast levels
    %   params - Model parameters [a, b, l, w, g0, n, rmx, e, o]
    %   optoStr - Mode ('baseline' or 'opto')
    %   options - (optional) Structure with customization options:
    %       .nTrials - Number of trials per contrast level (default: 10000)
    %       .verbose - Set to true to enable debug output (default: false)
    % 
    % Outputs:
    %   yPredicted - Predicted response values
    
    % Set up default options if none provided
    if nargin < 4
        options = struct();
    end
    
    % Set default options if not specified
    if ~isfield(options, 'verbose')
        options.verbose = false;
    end
    if ~isfield(options, 'nTrials')
        options.nTrials = 10000;
    end
    
    % Start timer if verbose
    if options.verbose
        startTime = tic;
    end
    
    %% Parameters - exactly like original script
    a = params(1);    % Weight of iso-tuned visual input (alpha)
    b = params(2);    % Weight of ortho-tuned visual input (beta)
    l = params(3);    % Weight of iso-tuned visual input norm (lambda)
    w = params(4);    % Weight of ortho-tuned visual input norm (omega)
    g0 = params(5);   % Normalization constant
    n = params(6);    % Spiking exponent
    rmx = params(7);  % Max response
    e = params(8);    % Relative weight of optostim effects (eta)
    
    % Set optostim effective contrast based on mode
    if strcmp(optoStr, 'baseline')
        o = 0;  % Zero opto-stim effective contrast for baseline
    else
        o = params(9);  % Opto-stim effective contrast
    end

    % Stimulation setup
    contrasts = x(:);      % Stimulus contrast levels as column vector
    nContrasts = numel(contrasts);
    nTrials = options.nTrials;  % Trials per contrast level
    
    % Simple consistent trials per contrast level like the original
    trialsPerContrast = nTrials * ones(nContrasts, 1);
    
    if options.verbose
        fprintf('Using %d trials per contrast level (uniform allocation)\n', nTrials);
    end

    %% Precompute all normalized responses - same equations as original
    % H/V column + H visual + V opto (vH_oV)
    resp_cH_vH_oV = contrasts + b * (1 - e) * o;
    norm_cH_vH_oV = contrasts + w * e * o;
    normresp_cH_vH_oV = rmx * (resp_cH_vH_oV.^n) ./ (norm_cH_vH_oV.^n + g0^n);
    
    resp_cV_vH_oV = a * contrasts + (1 - e) * o;
    norm_cV_vH_oV = l * contrasts + e * o;
    normresp_cV_vH_oV = rmx * (resp_cV_vH_oV.^n) ./ (norm_cV_vH_oV.^n + g0^n);
    
    % H/V column + H visual + H opto (vH_oH)
    resp_cH_vH_oH = contrasts + (1 - e) * o;
    norm_cH_vH_oH = contrasts + e * o;
    normresp_cH_vH_oH = rmx * (resp_cH_vH_oH.^n) ./ (norm_cH_vH_oH.^n + g0^n);
    
    resp_cV_vH_oH = a * contrasts + b * (1 - e) * o;
    norm_cV_vH_oH = l * contrasts + w * e * o;
    normresp_cV_vH_oH = rmx * (resp_cV_vH_oH.^n) ./ (norm_cV_vH_oH.^n + g0^n);
    
    % H/V column + V visual + V opto (vV_oV)
    resp_cH_vV_oV = a * contrasts + b * (1 - e) * o;
    norm_cH_vV_oV = l * contrasts + w * e * o;
    normresp_cH_vV_oV = rmx * (resp_cH_vV_oV.^n) ./ (norm_cH_vV_oV.^n + g0^n);
    
    resp_cV_vV_oV = contrasts + (1 - e) * o;
    norm_cV_vV_oV = contrasts + e * o;
    normresp_cV_vV_oV = rmx * (resp_cV_vV_oV.^n) ./ (norm_cV_vV_oV.^n + g0^n);
    
    % H/V column + V visual + H opto (vV_oH)
    resp_cH_vV_oH = a * contrasts + (1 - e) * o;
    norm_cH_vV_oH = l * contrasts + e * o;
    normresp_cH_vV_oH = rmx * (resp_cH_vV_oH.^n) ./ (norm_cH_vV_oH.^n + g0^n);
    
    resp_cV_vV_oH = contrasts + b * (1 - e) * o;
    norm_cV_vV_oH = contrasts + w * e * o;
    normresp_cV_vV_oH = rmx * (resp_cV_vV_oH.^n) ./ (norm_cV_vV_oH.^n + g0^n);
    
    %% Initialize counters for all conditions and contrast levels
    congruentTrialCount = zeros(nContrasts, 1);
    incongruentHVTrialCount = zeros(nContrasts, 1);
    incongruentVHTrialCount = zeros(nContrasts, 1);
    baselineTrialCount = zeros(nContrasts, 1);
    correctCongruentTrialCount = zeros(nContrasts, 1);
    correctIncongruentHVTrialCount = zeros(nContrasts, 1);
    correctIncongruentVHTrialCount = zeros(nContrasts, 1);
    correctBaselineTrialCount = zeros(nContrasts, 1);
    
    %% Process each contrast level serially - like original script
    for contrastLevel = 1:nContrasts
        % Predefine trial counts per contrast level
        nTrialsForThisContrast = trialsPerContrast(contrastLevel);
        
        for trial = 1:nTrialsForThisContrast
            % Generate random stimulus conditions (0=V, 1=H)
            visualOrt = randi([0 1]); % Visual orientation: 0=vertical, 1=horizontal
            optoOrt = randi([0 1]);   % Opto orientation: 0=vertical, 1=horizontal
            
            % Get precomputed responses for this contrast level
            if visualOrt == 1 && optoOrt == 0  % Visual=H, Opto=V
                respH = randn() + normresp_cH_vH_oV(contrastLevel);
                respV = randn() + normresp_cV_vH_oV(contrastLevel);
                
                incongruentHVTrialCount(contrastLevel) = incongruentHVTrialCount(contrastLevel) + 1;
                
            elseif visualOrt == 1 && optoOrt == 1  % Visual=H, Opto=H
                respH = randn() + normresp_cH_vH_oH(contrastLevel);
                respV = randn() + normresp_cV_vH_oH(contrastLevel);
                
                congruentTrialCount(contrastLevel) = congruentTrialCount(contrastLevel) + 1;
                
            elseif visualOrt == 0 && optoOrt == 0  % Visual=V, Opto=V
                respH = randn() + normresp_cH_vV_oV(contrastLevel);
                respV = randn() + normresp_cV_vV_oV(contrastLevel);
                
                baselineTrialCount(contrastLevel) = baselineTrialCount(contrastLevel) + 1;
                
            else  % visualOrt == 0 && optoOrt == 1, Visual=V, Opto=H
                respH = randn() + normresp_cH_vV_oH(contrastLevel);
                respV = randn() + normresp_cV_vV_oH(contrastLevel);
                
                incongruentVHTrialCount(contrastLevel) = incongruentVHTrialCount(contrastLevel) + 1;
            end
            
            % Calculate likelihood ratios using the original approach
            % Sum over all possible likelihoods for horizontal and vertical hypotheses
            
            % Horizontal hypothesis (visual=H)
            d1 = (respH - normresp_cH_vH_oV(contrastLevel))^2 + (respV - normresp_cV_vH_oV(contrastLevel))^2;
            d2 = (respH - normresp_cH_vH_oH(contrastLevel))^2 + (respV - normresp_cV_vH_oH(contrastLevel))^2;
            num = exp(-0.5 * d1) + exp(-0.5 * d2);
            
            % Vertical hypothesis (visual=V)
            d1 = (respH - normresp_cH_vV_oV(contrastLevel))^2 + (respV - normresp_cV_vV_oV(contrastLevel))^2;
            d2 = (respH - normresp_cH_vV_oH(contrastLevel))^2 + (respV - normresp_cV_vV_oH(contrastLevel))^2;
            den = exp(-0.5 * d1) + exp(-0.5 * d2);
            
            % Compute likelihood ratio and make decision
            lr = num / den;
            
            % Decision rule
            chosenHorz = false; % Default is "chosen vertical"
            if lr > 1.0
                chosenHorz = true; % Chose horizontal
            elseif lr == 1.0
                if rand() > 0.5
                    chosenHorz = true; % Randomized decision for tie
                end
            end
            
            % Track correct responses
            isCorrect = (visualOrt == 1 && chosenHorz) || (visualOrt == 0 && ~chosenHorz);
            
            % Count correct trials by condition
            if visualOrt == 1 && optoOrt == 1 && isCorrect  % Congruent (H-H)
                correctCongruentTrialCount(contrastLevel) = correctCongruentTrialCount(contrastLevel) + 1;
            elseif visualOrt == 1 && optoOrt == 0 && isCorrect  % Incongruent (H-V)
                correctIncongruentHVTrialCount(contrastLevel) = correctIncongruentHVTrialCount(contrastLevel) + 1;
            elseif visualOrt == 0 && optoOrt == 1 && isCorrect  % Incongruent (V-H)
                correctIncongruentVHTrialCount(contrastLevel) = correctIncongruentVHTrialCount(contrastLevel) + 1;
            elseif visualOrt == 0 && optoOrt == 0 && isCorrect  % Baseline (V-V)
                correctBaselineTrialCount(contrastLevel) = correctBaselineTrialCount(contrastLevel) + 1;
            end
        end
    end

    %% Compute probabilities with safety against division by zero
    epsilon = 1e-10; % Small value to avoid division by zero
    
    probCorrectCongruent = (correctCongruentTrialCount + correctBaselineTrialCount) ./ ...
                           (congruentTrialCount + baselineTrialCount + epsilon);
    
    probCorrectIncongruent = (correctIncongruentHVTrialCount + correctIncongruentVHTrialCount) ./ ...
                             (incongruentHVTrialCount + incongruentVHTrialCount + epsilon);
    
    probCorrectBaseline = (correctCongruentTrialCount + correctIncongruentHVTrialCount + ...
                           correctIncongruentVHTrialCount + correctBaselineTrialCount) ./ ...
                          (congruentTrialCount + incongruentHVTrialCount + ...
                           incongruentVHTrialCount + baselineTrialCount + epsilon);

    % Define indices for contrast levels
    zeroIdx = find(contrasts == 0); % Index for zero contrast
    inconIdx = find(contrasts < 0); % Indices for incongruent contrasts
    conIdx = find(contrasts > 0); % Indices for congruent contrasts

    % Use the same concatenation approach to ensure compatibility
    if strcmp(optoStr, 'baseline')
        % Format exactly as in original function to ensure matching output structure
        zeroPart = rmnan(mean([100 - (probCorrectBaseline(zeroIdx) * 100), probCorrectBaseline(zeroIdx) * 100]))';
        yPredicted = [100 - (probCorrectBaseline(inconIdx) * 100); ...
                      zeroPart; ...
                      probCorrectBaseline(conIdx) * 100]';
    else % 'opto'
        % Format exactly as in original function to ensure matching output structure
        zeroPart = rmnan(mean([100 - (probCorrectIncongruent(zeroIdx) * 100), probCorrectCongruent(zeroIdx) * 100]))';
        yPredicted = [100 - (probCorrectIncongruent(inconIdx) * 100); ...
                      zeroPart; ...
                      probCorrectCongruent(conIdx) * 100]';
    end
    
    % Report total execution time only if verbose
    if options.verbose
        executionTime = toc(startTime);
        fprintf('normMdlSimSimplified completed in %.3f seconds\n', executionTime);
    end
end

% Helper function to remove NaN values (simple implementation)
function cleaned = rmnan(x)
    cleaned = x;
    cleaned(isnan(cleaned)) = [];
end