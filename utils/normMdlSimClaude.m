function yPredicted = normMdlSimClaude(x, params, optoStr, options)
    % v_o_psychometric_func - Refactored from v_o_psychometric_3.m
    % Optimized for performance while maintaining statistical equivalence 
    %
    % Inputs:
    %   x - Stimulus contrast levels
    %   params - Model parameters [a, b, l, w, g0, n, rmx, e, o]
    %   optoStr - Mode ('baseline' or 'opto')
    %   options - (optional) Structure with customization options:
    %       .useGPU - Set to true to use GPU (default: auto-detect)
    %       .useParallel - Set to true to use parallel processing (default: false)
    %       .verbose - Set to true to enable debug output (default: false)
    %       .nTrials - Number of trials per contrast level (default: 100000)
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
        options.nTrials = 100000; % Default to high number of trials as in original
    end
    if ~isfield(options, 'useGPU')
        try
            options.useGPU = gpuDeviceCount > 0;
        catch
            options.useGPU = false;
        end
    end
    if ~isfield(options, 'useParallel')
        options.useParallel = false;
    end
    
    % Start timer if verbose
    if options.verbose
        startTime = tic;
        fprintf('Starting v_o_psychometric_func with %d trials per contrast level\n', options.nTrials);
    end
    
    %% Extract parameters
    a = params(1);    % alpha, visual excitatory cross-talk parameter
    b = params(2);    % beta, opto-stim excitatory cross-talk parameter
    l = params(3);    % lambda, visual normalization cross-talk parameter
    w = params(4);    % omega, opto-stim normalization cross-talk parameter
    g0 = params(5);   % normalization constant
    n = params(6);    % spiking exponent
    rmx = params(7);  % max response
    e = params(8);    % eta, relative weight of optostim on normalization vs. excitation
    
    % Set optostim effective contrast based on mode
    if strcmp(optoStr, 'baseline')
        o = 0;  % Zero opto-stim effective contrast
    else
        o = params(9);  % Opto-stim effective contrast
    end
    
    % Setup contrasts 
    contrasts = x(:);     % Ensure column vector
    nContrasts = length(contrasts);
    nTrials = options.nTrials;
    
    % Create exponential lookup table for fast exp calculation
    persistent expLookup expLookupRange
    if isempty(expLookup)
        expLookupRange = -20:0.01:0;  % Range for exp(x) inputs
        expLookup = exp(expLookupRange);
    end
    
    % Fast exponential function using lookup table
    fastExp = @(x) lookup_exp(x, expLookup, expLookupRange);
    
    %% Vectorized computation of neural responses for all contrasts
    if options.verbose
        if options.useGPU
            fprintf('Computing neural responses using GPU\n');
        else
            fprintf('Computing neural responses using CPU\n');
        end
    end
    
    % Use GPU if available and requested
    if options.useGPU
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
            
            % Calculate all mean responses on GPU
            
            % H (horizontal) column + V (vertical) visual + V (vertical) opto - notation: cH_vV_oV
            e_cH_vV_oV = d_a * d_contrasts + d_b * (1 - d_e) * d_o;
            g_cH_vV_oV = d_l * d_contrasts + d_w * d_e * d_o;
            u_cH_vV_oV = d_rmx * (e_cH_vV_oV.^d_n) ./ (g_cH_vV_oV.^d_n + d_g0^d_n);
            
            e_cV_vV_oV = d_contrasts + (1 - d_e) * d_o;
            g_cV_vV_oV = d_contrasts + d_e * d_o;
            u_cV_vV_oV = d_rmx * (e_cV_vV_oV.^d_n) ./ (g_cV_vV_oV.^d_n + d_g0^d_n);
            
            % H column + V visual + H opto
            e_cH_vV_oH = d_a * d_contrasts + (1 - d_e) * d_o;
            g_cH_vV_oH = d_l * d_contrasts + d_e * d_o;
            u_cH_vV_oH = d_rmx * (e_cH_vV_oH.^d_n) ./ (g_cH_vV_oH.^d_n + d_g0^d_n);
            
            e_cV_vV_oH = d_contrasts + d_b * (1 - d_e) * d_o;
            g_cV_vV_oH = d_contrasts + d_w * d_e * d_o;
            u_cV_vV_oH = d_rmx * (e_cV_vV_oH.^d_n) ./ (g_cV_vV_oH.^d_n + d_g0^d_n);
            
            % H column + H visual + V opto
            e_cH_vH_oV = d_contrasts + d_b * (1 - d_e) * d_o;
            g_cH_vH_oV = d_contrasts + d_w * d_e * d_o;
            u_cH_vH_oV = d_rmx * (e_cH_vH_oV.^d_n) ./ (g_cH_vH_oV.^d_n + d_g0^d_n);
            
            e_cV_vH_oV = d_a * d_contrasts + (1 - d_e) * d_o;
            g_cV_vH_oV = d_l * d_contrasts + d_e * d_o;
            u_cV_vH_oV = d_rmx * (e_cV_vH_oV.^d_n) ./ (g_cV_vH_oV.^d_n + d_g0^d_n);
            
            % H column + H visual + H opto
            e_cH_vH_oH = d_contrasts + (1 - d_e) * d_o;
            g_cH_vH_oH = d_contrasts + d_e * d_o;
            u_cH_vH_oH = d_rmx * (e_cH_vH_oH.^d_n) ./ (g_cH_vH_oH.^d_n + d_g0^d_n);
            
            e_cV_vH_oH = d_a * d_contrasts + d_b * (1 - d_e) * d_o;
            g_cV_vH_oH = d_l * d_contrasts + d_w * d_e * d_o;
            u_cV_vH_oH = d_rmx * (e_cV_vH_oH.^d_n) ./ (g_cV_vH_oH.^d_n + d_g0^d_n);
            
            % Gather results back to CPU
            u_cH_vV_oV = gather(u_cH_vV_oV);
            u_cV_vV_oV = gather(u_cV_vV_oV);
            u_cH_vV_oH = gather(u_cH_vV_oH);
            u_cV_vV_oH = gather(u_cV_vV_oH);
            u_cH_vH_oV = gather(u_cH_vH_oV);
            u_cV_vH_oV = gather(u_cV_vH_oV);
            u_cH_vH_oH = gather(u_cH_vH_oH);
            u_cV_vH_oH = gather(u_cV_vH_oH);
            
            % Clear GPU memory
            clear d_contrasts d_a d_b d_l d_w d_g0 d_n d_rmx d_e d_o;
            clear e_cH_vV_oV g_cH_vV_oV e_cV_vV_oV g_cV_vV_oV;
            clear e_cH_vV_oH g_cH_vV_oH e_cV_vV_oH g_cV_vV_oH;
            clear e_cH_vH_oV g_cH_vH_oV e_cV_vH_oV g_cV_vH_oV;
            clear e_cH_vH_oH g_cH_vH_oH e_cV_vH_oH g_cV_vH_oH;
            
        catch gpuErr
            if options.verbose
                fprintf('GPU processing failed: %s\nFalling back to CPU\n', gpuErr.message);
            end
            options.useGPU = false;
        end
    end
    
    % CPU computation if GPU not used or failed
    if ~options.useGPU
        % H (horizontal) column + V (vertical) visual + V (vertical) opto - notation: cH_vV_oV
        e_cH_vV_oV = a * contrasts + b * (1 - e) * o;
        g_cH_vV_oV = l * contrasts + w * e * o;
        u_cH_vV_oV = rmx * (e_cH_vV_oV.^n) ./ (g_cH_vV_oV.^n + g0^n);
        
        e_cV_vV_oV = contrasts + (1 - e) * o;
        g_cV_vV_oV = contrasts + e * o;
        u_cV_vV_oV = rmx * (e_cV_vV_oV.^n) ./ (g_cV_vV_oV.^n + g0^n);
        
        % H column + V visual + H opto
        e_cH_vV_oH = a * contrasts + (1 - e) * o;
        g_cH_vV_oH = l * contrasts + e * o;
        u_cH_vV_oH = rmx * (e_cH_vV_oH.^n) ./ (g_cH_vV_oH.^n + g0^n);
        
        e_cV_vV_oH = contrasts + b * (1 - e) * o;
        g_cV_vV_oH = contrasts + w * e * o;
        u_cV_vV_oH = rmx * (e_cV_vV_oH.^n) ./ (g_cV_vV_oH.^n + g0^n);
        
        % H column + H visual + V opto
        e_cH_vH_oV = contrasts + b * (1 - e) * o;
        g_cH_vH_oV = contrasts + w * e * o;
        u_cH_vH_oV = rmx * (e_cH_vH_oV.^n) ./ (g_cH_vH_oV.^n + g0^n);
        
        e_cV_vH_oV = a * contrasts + (1 - e) * o;
        g_cV_vH_oV = l * contrasts + e * o;
        u_cV_vH_oV = rmx * (e_cV_vH_oV.^n) ./ (g_cV_vH_oV.^n + g0^n);
        
        % H column + H visual + H opto
        e_cH_vH_oH = contrasts + (1 - e) * o;
        g_cH_vH_oH = contrasts + e * o;
        u_cH_vH_oH = rmx * (e_cH_vH_oH.^n) ./ (g_cH_vH_oH.^n + g0^n);
        
        e_cV_vH_oH = a * contrasts + b * (1 - e) * o;
        g_cV_vH_oH = l * contrasts + w * e * o;
        u_cV_vH_oH = rmx * (e_cV_vH_oH.^n) ./ (g_cV_vH_oH.^n + g0^n);
    end
    
    %% Initialize counters and run trials
    
    % Preallocate result matrices - one row per contrast
    output = zeros(nContrasts, 10);
    output(:,1) = contrasts;    % Contrast
    output(:,2) = o;            % Optostim contrast
    
    % Check if parallel processing should be used
    if options.useParallel
        if options.verbose
            fprintf('Using parallel processing for contrast levels\n');
        end
        
        try
            % Start parallel pool if not already running
            if isempty(gcp('nocreate'))
                parpool('local', 'IdleTimeout', 120);
            end
            
            % FIX: Create a parallel-safe array for collecting results
            parResults = zeros(nContrasts, 8);
            
            % Process in parallel for each contrast level
            parfor c_idx = 1:nContrasts
                % Store results in a parallel-safe way
                parResults(c_idx, :) = runTrialsForContrast(c_idx, nTrials, ...
                    u_cH_vV_oV(c_idx), u_cV_vV_oV(c_idx), ...
                    u_cH_vV_oH(c_idx), u_cV_vV_oH(c_idx), ...
                    u_cH_vH_oV(c_idx), u_cV_vH_oV(c_idx), ...
                    u_cH_vH_oH(c_idx), u_cV_vH_oH(c_idx), ...
                    fastExp);
            end
            
            % After parfor, assign to output matrix
            output(:, 3:10) = parResults;
            
        catch parErr
            if options.verbose
                fprintf('Parallel processing failed: %s\nFalling back to serial processing\n', parErr.message);
            end
            options.useParallel = false;
        end
    end
    
    % Serial processing if parallel is disabled or failed
    if ~options.useParallel
        if options.verbose
            fprintf('Using serial processing for contrast levels\n');
        end
        
        for c_idx = 1:nContrasts
            output(c_idx, 3:10) = runTrialsForContrast(c_idx, nTrials, ...
                u_cH_vV_oV(c_idx), u_cV_vV_oV(c_idx), ...
                u_cH_vV_oH(c_idx), u_cV_vV_oH(c_idx), ...
                u_cH_vH_oV(c_idx), u_cV_vH_oV(c_idx), ...
                u_cH_vH_oH(c_idx), u_cV_vH_oH(c_idx), ...
                fastExp);
        end
    end
    
    %% Calculate probabilities and format output
    
    % Calculate probabilities
    pc = (output(:,7) + output(:,10)) ./ (output(:,3) + output(:,6));
    pic = (output(:,8) + output(:,9)) ./ (output(:,4) + output(:,5));
    
    % Process baseline condition if needed
    if strcmp(optoStr, 'baseline')
        % Calculate all probabilities as in original
        pcntrl = (output(:,7) + output(:,8) + output(:,9) + output(:,10)) ./ ...
                 (output(:,3) + output(:,4) + output(:,5) + output(:,6));
        
        % Format output matching the original
        zeroIdx = find(contrasts == 0);
        inconIdx = find(contrasts < 0);
        conIdx = find(contrasts > 0);
        
        % Handle zero contrast (use mean for symmetry as in original)
        zeroPart = rmnan(mean([100 - (pcntrl(zeroIdx) * 100), pcntrl(zeroIdx) * 100]))';
        
        % Format exactly as in original function
        yPredicted = [100 - (pcntrl(inconIdx) * 100); ...
                      zeroPart; ...
                      pcntrl(conIdx) * 100]';
    else
        % For opto condition
        zeroIdx = find(contrasts == 0);
        inconIdx = find(contrasts < 0);
        conIdx = find(contrasts > 0);
        
        % Handle zero contrast (use mean for symmetry as in original)
        zeroPart = rmnan(mean([100 - (pic(zeroIdx) * 100), pc(zeroIdx) * 100]))';
        
        % Format output
        yPredicted = [100 - (pic(inconIdx) * 100); ...
                      zeroPart; ...
                      pc(conIdx) * 100]';
    end
    
    % Report timing if verbose
    if options.verbose
        executionTime = toc(startTime);
        fprintf('v_o_psychometric_func completed in %.3f seconds\n', executionTime);
    end
end

% Helper function to run trials for a specific contrast level
function result = runTrialsForContrast(c_idx, nTrials, ...
    u_cH_vV_oV, u_cV_vV_oV, u_cH_vV_oH, u_cV_vV_oH, ...
    u_cH_vH_oV, u_cV_vH_oV, u_cH_vH_oH, u_cV_vH_oH, fastExp)
    
    % Initialize results for this contrast
    % [nhh, nhv, nvh, nvv, chh, chv, cvh, cvv]
    result = zeros(1, 8);
    
    % Process trials in batches for better performance
    batchSize = min(10000, nTrials); % Balance memory and performance
    numBatches = ceil(nTrials / batchSize);
    
    for batch = 1:numBatches
        % Calculate current batch size
        currentBatchSize = min(batchSize, nTrials - (batch-1)*batchSize);
        
        % Generate random visual and optogenetic orientations (0=V, 1=H)
        visualOrt = randi([0 1], currentBatchSize, 1);
        optoOrt = randi([0 1], currentBatchSize, 1);
        
        % Generate Gaussian noise
        noiseH = randn(currentBatchSize, 1);
        noiseV = randn(currentBatchSize, 1);
        
        % Preallocate response arrays
        respH = zeros(currentBatchSize, 1);
        respV = zeros(currentBatchSize, 1);
        
        % Compute responses based on stimulus conditions (vectorized)
        mask_vV_oV = (visualOrt == 0) & (optoOrt == 0); % V visual, V opto
        mask_vV_oH = (visualOrt == 0) & (optoOrt == 1); % V visual, H opto
        mask_vH_oV = (visualOrt == 1) & (optoOrt == 0); % H visual, V opto
        mask_vH_oH = (visualOrt == 1) & (optoOrt == 1); % H visual, H opto
        
        % Count trials per condition
        n_vV_oV = sum(mask_vV_oV);
        n_vV_oH = sum(mask_vV_oH);
        n_vH_oV = sum(mask_vH_oV);
        n_vH_oH = sum(mask_vH_oH);
        
        % Add to result counts (nhh, nhv, nvh, nvv)
        result(1) = result(1) + n_vH_oH; % nhh (Horizontal visual, Horizontal opto)
        result(2) = result(2) + n_vH_oV; % nhv (Horizontal visual, Vertical opto)
        result(3) = result(3) + n_vV_oH; % nvh (Vertical visual, Horizontal opto)
        result(4) = result(4) + n_vV_oV; % nvv (Vertical visual, Vertical opto)
        
        % Calculate responses for each condition
        if n_vV_oV > 0
            respH(mask_vV_oV) = noiseH(mask_vV_oV) + u_cH_vV_oV;
            respV(mask_vV_oV) = noiseV(mask_vV_oV) + u_cV_vV_oV;
        end
        
        if n_vV_oH > 0
            respH(mask_vV_oH) = noiseH(mask_vV_oH) + u_cH_vV_oH;
            respV(mask_vV_oH) = noiseV(mask_vV_oH) + u_cV_vV_oH;
        end
        
        if n_vH_oV > 0
            respH(mask_vH_oV) = noiseH(mask_vH_oV) + u_cH_vH_oV;
            respV(mask_vH_oV) = noiseV(mask_vH_oV) + u_cV_vH_oV;
        end
        
        if n_vH_oH > 0
            respH(mask_vH_oH) = noiseH(mask_vH_oH) + u_cH_vH_oH;
            respV(mask_vH_oH) = noiseV(mask_vH_oH) + u_cV_vH_oH;
        end
        
        % Process each trial for decision
        for i = 1:currentBatchSize
            % Horizontal hypothesis (visual=H)
            d1 = (respH(i) - u_cH_vH_oV)^2 + (respV(i) - u_cV_vH_oV)^2;
            d2 = (respH(i) - u_cH_vH_oH)^2 + (respV(i) - u_cV_vH_oH)^2;
            lrNumerator = fastExp(-0.5 * d1) + fastExp(-0.5 * d2);
            
            % Vertical hypothesis (visual=V)
            d1 = (respH(i) - u_cH_vV_oV)^2 + (respV(i) - u_cV_vV_oV)^2;
            d2 = (respH(i) - u_cH_vV_oH)^2 + (respV(i) - u_cV_vV_oH)^2;
            lrDenominator = fastExp(-0.5 * d1) + fastExp(-0.5 * d2);
            
            % Decision
            lr = lrNumerator / lrDenominator;
            if lr > 1.0
                trialOutcome = 1; % Decided Horizontal
            elseif lr == 1.0
                trialOutcome = (rand() > 0.5); % Random on tie
            else
                trialOutcome = 0; % Decided Vertical
            end
            
            % Count correct trials
            if visualOrt(i) == trialOutcome
                if visualOrt(i) == 1 && optoOrt(i) == 1
                    result(5) = result(5) + 1; % chh
                elseif visualOrt(i) == 1 && optoOrt(i) == 0
                    result(6) = result(6) + 1; % chv
                elseif visualOrt(i) == 0 && optoOrt(i) == 1
                    result(7) = result(7) + 1; % cvh
                elseif visualOrt(i) == 0 && optoOrt(i) == 0
                    result(8) = result(8) + 1; % cvv
                end
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