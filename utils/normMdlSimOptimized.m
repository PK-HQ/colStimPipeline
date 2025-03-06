function yPredicted = normMdlSimOptimized(x, params, optoStr, options)
    % Fully customizable optimized version of normMdlSim
    % Modified to be more like the original v_o_psychometric_3.m:
    % 1. Uses even trial allocation (no smart allocation)
    % 2. Uses default random number generator
    % 3. Uses exponential lookup table for faster calculations
    % 
    % Inputs:
    %   x - Stimulus contrast levels
    %   params - Model parameters [a, b, l, w, g0, n, rmx, e, o]
    %   optoStr - Mode ('baseline' or 'opto')
    %   options - (optional) Structure with customization options:
    %       .useGPU - Set to true to use GPU (default: auto-detect)
    %       .useParallel - Set to true to use parallel processing (default: false)
    %       .verbose - Set to true to enable debug output (default: false)
    %       .nTrials - Number of trials per contrast level (default: 1000)
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
        options.nTrials = 1000;
    end
    
    % Start timer if verbose
    if options.verbose
        startTime = tic;
    end
    
    % Check for GPU availability if not explicitly set
    if ~isfield(options, 'useGPU')
        try
            options.useGPU = gpuDeviceCount > 0;
        catch
            options.useGPU = false;
        end
    end
    
    % Display GPU status
    if options.verbose
        if options.useGPU
            fprintf('Using GPU acceleration\n');
        else
            fprintf('Using CPU computation\n');
        end
    end
    
    % Check for parallel processing option
    if ~isfield(options, 'useParallel')
        options.useParallel = false;
    end
    
    % Display parallel status
    if options.verbose && options.useParallel
        fprintf('Using parallel processing\n');
    end
    
    %% Parameters
    a = params(1);    % Weight of iso-tuned visual input
    b = params(2);    % Weight of ortho-tuned visual input
    l = params(3);    % Weight of iso-tuned visual input norm
    w = params(4);    % Weight of ortho-tuned visual input norm
    g0 = params(5);   % Normalization constant
    n = params(6);    % Spiking exponent
    rmx = params(7);  % Max response
    e = params(8);    % Relative weight of optostim effects on excitation
    
    % Set optostim effective contrast based on mode
    if strcmp(optoStr, 'baseline')
        o = 0;  % Zero opto-stim effective contrast
    else
        o = params(9);  % Opto-stim effective contrast
    end

    % Stimulation setup
    contrasts = x(:);      % Stimulus contrast levels as column vector
    nContrasts = numel(contrasts);
    nTrials = options.nTrials;  % Trials per contrast level

    % Use uniform trial allocation like the original
    trialsPerContrast = nTrials * ones(nContrasts, 1);
    
    if options.verbose
        fprintf('Using uniform trial allocation (%d trials per contrast)\n', nTrials);
    end

    % Create exponential lookup table for fast exp calculation
    persistent expLookup expLookupRange
    if isempty(expLookup)
        expLookupRange = -20:0.01:0;  % Range for exp(x) inputs
        expLookup = exp(expLookupRange);
    end
    
    % Fast exponential function using lookup table
    fastExp = @(x) lookup_exp(x, expLookup, expLookupRange);

    %% Precompute all normalized responses - vectorized
    if options.useGPU
        % Use GPU for matrix calculations
        try
            d_contrasts = gpuArray(contrasts);
            d_a = gpuArray(a);
            d_b = gpuArray(b);
            d_l = gpuArray(l);
            d_w = gpuArray(w);
            d_g0 = gpuArray(g0);
            d_n = gpuArray(n);
            d_rmx = gpuArray(rmx);
            d_e = gpuArray(e);
            d_o = gpuArray(o);
            
            % H/V column + H visual + V opto (vH_oV)
            resp_cH_vH_oV = d_contrasts + d_b * (1 - d_e) * d_o;
            norm_cH_vH_oV = d_contrasts + d_w * d_e * d_o;
            normresp_cH_vH_oV = d_rmx * (resp_cH_vH_oV.^d_n) ./ (norm_cH_vH_oV.^d_n + d_g0^d_n);
            
            resp_cV_vH_oV = d_a * d_contrasts + (1 - d_e) * d_o;
            norm_cV_vH_oV = d_l * d_contrasts + d_e * d_o;
            normresp_cV_vH_oV = d_rmx * (resp_cV_vH_oV.^d_n) ./ (norm_cV_vH_oV.^d_n + d_g0^d_n);
            
            % H/V column + H visual + H opto (vH_oH)
            resp_cH_vH_oH = d_contrasts + (1 - d_e) * d_o;
            norm_cH_vH_oH = d_contrasts + d_e * d_o;
            normresp_cH_vH_oH = d_rmx * (resp_cH_vH_oH.^d_n) ./ (norm_cH_vH_oH.^d_n + d_g0^d_n);
            
            resp_cV_vH_oH = d_a * d_contrasts + d_b * (1 - d_e) * d_o;
            norm_cV_vH_oH = d_l * d_contrasts + d_w * d_e * d_o;
            normresp_cV_vH_oH = d_rmx * (resp_cV_vH_oH.^d_n) ./ (norm_cV_vH_oH.^d_n + d_g0^d_n);
            
            % H/V column + V visual + V opto (vV_oV)
            resp_cH_vV_oV = d_a * d_contrasts + d_b * (1 - d_e) * d_o;
            norm_cH_vV_oV = d_l * d_contrasts + d_w * d_e * d_o;
            normresp_cH_vV_oV = d_rmx * (resp_cH_vV_oV.^d_n) ./ (norm_cH_vV_oV.^d_n + d_g0^d_n);
            
            resp_cV_vV_oV = d_contrasts + (1 - d_e) * d_o;
            norm_cV_vV_oV = d_contrasts + d_e * d_o;
            normresp_cV_vV_oV = d_rmx * (resp_cV_vV_oV.^d_n) ./ (norm_cV_vV_oV.^d_n + d_g0^d_n);
            
            % H/V column + V visual + H opto (vV_oH)
            resp_cH_vV_oH = d_a * d_contrasts + (1 - d_e) * d_o;
            norm_cH_vV_oH = d_l * d_contrasts + d_e * d_o;
            normresp_cH_vV_oH = d_rmx * (resp_cH_vV_oH.^d_n) ./ (norm_cH_vV_oH.^d_n + d_g0^d_n);
            
            resp_cV_vV_oH = d_contrasts + d_b * (1 - d_e) * d_o;
            norm_cV_vV_oH = d_contrasts + d_w * d_e * d_o;
            normresp_cV_vV_oH = d_rmx * (resp_cV_vV_oH.^d_n) ./ (norm_cV_vV_oH.^d_n + d_g0^d_n);
            
            % Transfer back to CPU
            normresp_cH_vH_oV = gather(normresp_cH_vH_oV);
            normresp_cV_vH_oV = gather(normresp_cV_vH_oV);
            normresp_cH_vH_oH = gather(normresp_cH_vH_oH);
            normresp_cV_vH_oH = gather(normresp_cV_vH_oH);
            normresp_cH_vV_oV = gather(normresp_cH_vV_oV);
            normresp_cV_vV_oV = gather(normresp_cV_vV_oV);
            normresp_cH_vV_oH = gather(normresp_cH_vV_oH);
            normresp_cV_vV_oH = gather(normresp_cV_vV_oH);
            
            % Clear GPU variables to free memory
            clear d_contrasts d_a d_b d_l d_w d_g0 d_n d_rmx d_e d_o;
            clear resp_cH_vH_oV norm_cH_vH_oV resp_cV_vH_oV norm_cV_vH_oV;
            clear resp_cH_vH_oH norm_cH_vH_oH resp_cV_vH_oH norm_cV_vH_oH;
            clear resp_cH_vV_oV norm_cH_vV_oV resp_cV_vV_oV norm_cV_vV_oV;
            clear resp_cH_vV_oH norm_cH_vV_oH resp_cV_vV_oH norm_cV_vV_oH;
        catch
            if options.verbose
                fprintf('GPU processing failed, falling back to CPU\n');
            end
            options.useGPU = false;
            % Fall through to CPU calculation
        end
    end
    
    if ~options.useGPU
        % CPU calculation
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
    end

    %% Initialize counters for all contrast levels
    congruentTrialCount = zeros(nContrasts, 1);
    incongruentHVTrialCount = zeros(nContrasts, 1);
    incongruentVHTrialCount = zeros(nContrasts, 1);
    baselineTrialCount = zeros(nContrasts, 1);
    correctCongruentTrialCount = zeros(nContrasts, 1);
    correctIncongruentHVTrialCount = zeros(nContrasts, 1);
    correctIncongruentVHTrialCount = zeros(nContrasts, 1);
    correctBaselineTrialCount = zeros(nContrasts, 1);
    
    %% Process each contrast level (parallel or serial)
    if options.useParallel
        % Try to use parallel processing
        try
            % Check for parpool and start if needed with a short timeout
            if isempty(gcp('nocreate'))
                if options.verbose
                    fprintf('Starting parallel pool with 30 second timeout...\n');
                end
                parpool('local', 'Timeout', 30);
            end
            
            % Prepare for parallel processing
            parfor_congruentCounts = zeros(nContrasts, 1);
            parfor_incongruentHVCounts = zeros(nContrasts, 1);
            parfor_incongruentVHCounts = zeros(nContrasts, 1);
            parfor_baselineCounts = zeros(nContrasts, 1);
            parfor_correctCongruentCounts = zeros(nContrasts, 1);
            parfor_correctIncongruentHVCounts = zeros(nContrasts, 1);
            parfor_correctIncongruentVHCounts = zeros(nContrasts, 1);
            parfor_correctBaselineCounts = zeros(nContrasts, 1);
            
            % Process in parallel
            parfor contrastLevel = 1:nContrasts
                [c_count, ihv_count, ivh_count, b_count, ...
                 cc_count, cihv_count, civh_count, cb_count] = ...
                    processContrastLevel(contrastLevel, trialsPerContrast(contrastLevel), ...
                        normresp_cH_vH_oV, normresp_cV_vH_oV, ...
                        normresp_cH_vH_oH, normresp_cV_vH_oH, ...
                        normresp_cH_vV_oV, normresp_cV_vV_oV, ...
                        normresp_cH_vV_oH, normresp_cV_vV_oH, ...
                        fastExp); % Added fastExp parameter
                
                % Store results in parfor-compatible arrays
                parfor_congruentCounts(contrastLevel) = c_count;
                parfor_incongruentHVCounts(contrastLevel) = ihv_count;
                parfor_incongruentVHCounts(contrastLevel) = ivh_count;
                parfor_baselineCounts(contrastLevel) = b_count;
                parfor_correctCongruentCounts(contrastLevel) = cc_count;
                parfor_correctIncongruentHVCounts(contrastLevel) = cihv_count;
                parfor_correctIncongruentVHCounts(contrastLevel) = civh_count;
                parfor_correctBaselineCounts(contrastLevel) = cb_count;
            end
            
            % Merge results from parallel processing
            congruentTrialCount = parfor_congruentCounts;
            incongruentHVTrialCount = parfor_incongruentHVCounts;
            incongruentVHTrialCount = parfor_incongruentVHCounts;
            baselineTrialCount = parfor_baselineCounts;
            correctCongruentTrialCount = parfor_correctCongruentCounts;
            correctIncongruentHVTrialCount = parfor_correctIncongruentHVCounts;
            correctIncongruentVHTrialCount = parfor_correctIncongruentVHCounts;
            correctBaselineTrialCount = parfor_correctBaselineCounts;
            
            if options.verbose
                fprintf('Parallel processing completed successfully\n');
            end
        catch e
            % If parallel processing fails, fall back to serial
            if options.verbose
                fprintf('Parallel processing failed: %s\n', e.message);
                fprintf('Falling back to serial processing\n');
            end
            options.useParallel = false;
            % Fall through to serial processing
        end
    end
    
    % Serial processing if parallel is disabled or failed
    if ~options.useParallel
        if options.verbose
            fprintf('Using serial processing for %d contrast levels\n', nContrasts);
        end
        
        for contrastLevel = 1:nContrasts
            [congruentTrialCount(contrastLevel), ...
             incongruentHVTrialCount(contrastLevel), ...
             incongruentVHTrialCount(contrastLevel), ...
             baselineTrialCount(contrastLevel), ...
             correctCongruentTrialCount(contrastLevel), ...
             correctIncongruentHVTrialCount(contrastLevel), ...
             correctIncongruentVHTrialCount(contrastLevel), ...
             correctBaselineTrialCount(contrastLevel)] = ...
                processContrastLevel(contrastLevel, trialsPerContrast(contrastLevel), ...
                    normresp_cH_vH_oV, normresp_cV_vH_oV, ...
                    normresp_cH_vH_oH, normresp_cV_vH_oH, ...
                    normresp_cH_vV_oV, normresp_cV_vV_oV, ...
                    normresp_cH_vV_oH, normresp_cV_vV_oH, ...
                    fastExp); % Added fastExp parameter
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

    % Use the same concatenation approach as the original function to ensure compatibility
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
        fprintf('normMdlSim completed in %.3f seconds (GPU: %d, Parallel: %d)\n', ...
            executionTime, options.useGPU, options.useParallel);
    end
end

% Helper function for processing a single contrast level
function [c_count, ihv_count, ivh_count, b_count, ...
          cc_count, cihv_count, civh_count, cb_count] = ...
    processContrastLevel(contrastLevel, nTrials, ...
        normresp_cH_vH_oV, normresp_cV_vH_oV, ...
        normresp_cH_vH_oH, normresp_cV_vH_oH, ...
        normresp_cH_vV_oV, normresp_cV_vV_oV, ...
        normresp_cH_vV_oH, normresp_cV_vV_oH, ...
        fastExp) % Added fastExp parameter
    
    % Initialize local counters for this contrast level
    c_count = 0;   % congruentTrialCount
    ihv_count = 0; % incongruentHVTrialCount
    ivh_count = 0; % incongruentVHTrialCount
    b_count = 0;   % baselineTrialCount
    cc_count = 0;  % correctCongruentTrialCount
    cihv_count = 0; % correctIncongruentHVTrialCount
    civh_count = 0; % correctIncongruentVHTrialCount
    cb_count = 0;  % correctBaselineTrialCount
    
    % Access precomputed responses for this contrast
    resp_cH_vH_oV_i = normresp_cH_vH_oV(contrastLevel);
    resp_cV_vH_oV_i = normresp_cV_vH_oV(contrastLevel);
    resp_cH_vH_oH_i = normresp_cH_vH_oH(contrastLevel);
    resp_cV_vH_oH_i = normresp_cV_vH_oH(contrastLevel);
    resp_cH_vV_oV_i = normresp_cH_vV_oV(contrastLevel);
    resp_cV_vV_oV_i = normresp_cV_vV_oV(contrastLevel);
    resp_cH_vV_oH_i = normresp_cH_vV_oH(contrastLevel);
    resp_cV_vV_oH_i = normresp_cV_vV_oH(contrastLevel);
    
    % Process trials in efficient batched approach
    batchSize = min(1000, nTrials);
    numBatches = ceil(nTrials / batchSize);
    
    for batchIdx = 1:numBatches
        % Calculate batch size for this iteration
        startIdx = (batchIdx - 1) * batchSize + 1;
        endIdx = min(startIdx + batchSize - 1, nTrials);
        currentBatchSize = endIdx - startIdx + 1;
        
        % Generate random stimulus conditions for this batch
        visualOrt = randi([0 1], currentBatchSize, 1);  % 0=0°(V), 1=90°(H)
        optoOrt = randi([0 1], currentBatchSize, 1);    % 0=0°(V), 1=90°(H)
        
        % Generate random Gaussian noise
        noiseH = randn(currentBatchSize, 1);
        noiseV = randn(currentBatchSize, 1);
        
        % Preallocate response arrays
        respH = zeros(currentBatchSize, 1);
        respV = zeros(currentBatchSize, 1);
        
        % Compute responses using vectorized operations
        % Visual=V (0), Opto=V (0)
        mask_vV_oV = (visualOrt == 0) & (optoOrt == 0);
        if any(mask_vV_oV)
            respH(mask_vV_oV) = noiseH(mask_vV_oV) + resp_cH_vV_oV_i;
            respV(mask_vV_oV) = noiseV(mask_vV_oV) + resp_cV_vV_oV_i;
            b_count = b_count + sum(mask_vV_oV);
        end
        
        % Visual=V (0), Opto=H (1)
        mask_vV_oH = (visualOrt == 0) & (optoOrt == 1);
        if any(mask_vV_oH)
            respH(mask_vV_oH) = noiseH(mask_vV_oH) + resp_cH_vV_oH_i;
            respV(mask_vV_oH) = noiseV(mask_vV_oH) + resp_cV_vV_oH_i;
            ivh_count = ivh_count + sum(mask_vV_oH);
        end
        
        % Visual=H (1), Opto=V (0)
        mask_vH_oV = (visualOrt == 1) & (optoOrt == 0);
        if any(mask_vH_oV)
            respH(mask_vH_oV) = noiseH(mask_vH_oV) + resp_cH_vH_oV_i;
            respV(mask_vH_oV) = noiseV(mask_vH_oV) + resp_cV_vH_oV_i;
            ihv_count = ihv_count + sum(mask_vH_oV);
        end
        
        % Visual=H (1), Opto=H (1)
        mask_vH_oH = (visualOrt == 1) & (optoOrt == 1);
        if any(mask_vH_oH)
            respH(mask_vH_oH) = noiseH(mask_vH_oH) + resp_cH_vH_oH_i;
            respV(mask_vH_oH) = noiseV(mask_vH_oH) + resp_cV_vH_oH_i;
            c_count = c_count + sum(mask_vH_oH);
        end
        
        % For each trial in batch, calculate likelihood and decision
        for i = 1:currentBatchSize
            % Horizontal hypothesis (visual=H)
            d1 = (respH(i) - resp_cH_vH_oV_i)^2 + (respV(i) - resp_cV_vH_oV_i)^2;
            d2 = (respH(i) - resp_cH_vH_oH_i)^2 + (respV(i) - resp_cV_vH_oH_i)^2;
            % Use fastExp instead of standard exp
            lrNumerator = fastExp(-0.5 * d1) + fastExp(-0.5 * d2);
            
            % Vertical hypothesis (visual=V)
            d1 = (respH(i) - resp_cH_vV_oV_i)^2 + (respV(i) - resp_cV_vV_oV_i)^2;
            d2 = (respH(i) - resp_cH_vV_oH_i)^2 + (respV(i) - resp_cV_vV_oH_i)^2;
            % Use fastExp instead of standard exp
            lrDenominator = fastExp(-0.5 * d1) + fastExp(-0.5 * d2);
            
            % Decision
            lr = lrNumerator / lrDenominator;
            if lr > 1.0
                trialOutcome = 1;  % Decided Horizontal
            elseif lr == 1.0
                trialOutcome = (rand() > 0.5);  % Random decision on tie
            else
                trialOutcome = 0;  % Decided Vertical
            end
            
            % Update correct counters
            if visualOrt(i) == 0 && optoOrt(i) == 0 && trialOutcome == visualOrt(i)
                cb_count = cb_count + 1;
            elseif visualOrt(i) == 0 && optoOrt(i) == 1 && trialOutcome == visualOrt(i)
                civh_count = civh_count + 1;
            elseif visualOrt(i) == 1 && optoOrt(i) == 0 && trialOutcome == visualOrt(i)
                cihv_count = cihv_count + 1;
            elseif visualOrt(i) == 1 && optoOrt(i) == 1 && trialOutcome == visualOrt(i)
                cc_count = cc_count + 1;
            end
        end
    end
end

% Helper function for exponential lookup
function y = lookup_exp(x, expLookup, expRange)
    % Fast exponential approximation using table lookup
    
    if x <= expRange(1)
        y = 0;  % Underflow to zero for very small values
    elseif x >= expRange(end)
        y = 1;  % exp(0) = 1
    else
        % Find index in lookup table
        idx = round((x - expRange(1)) * 100) + 1;
        
        % Ensure index is within bounds
        idx = max(1, min(length(expLookup), idx));
        
        % Return value from lookup table
        y = expLookup(idx);
    end
end

% Helper function to remove NaN values
function y = rmnan(x)
    y = x;
    y(isnan(y)) = [];
end