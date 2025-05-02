function yPredicted = normMdlSimBillOptimized(x, params, optoStr, options)
    % Handle default options
    if nargin < 4
        options = struct();
    end
    if ~isfield(options, 'nTrials')
        options.nTrials = 100000;
    end
    if ~isfield(options, 'useGPU')
        options.useGPU = false;
    end
    
    % Function to handle array creation based on GPU preference
    function arr = createArray(data)
        if options.useGPU
            arr = gpuArray(data);
        else
            arr = data;
        end
    end
    
    % Extract parameters and move to GPU if needed
    a = createArray(params(1));    % alpha
    b = createArray(params(2));    % beta
    l = createArray(params(3));    % lambda
    w = createArray(params(4));    % omega
    g0 = createArray(params(5));   % normalization constant
    n = params(6);                 % spiking exponent - keep on CPU
    rmx = createArray(params(7));  % max response
    e = createArray(params(8));    % eta
    
    if strcmp(optoStr, 'baseline')
        o = createArray(0);
    else
        o = createArray(params(9));
    end
    
    % Setup contrasts
    c = createArray(x(:)');
    nlev = length(c);
    ntrl = options.nTrials;
    
    % Pre-compute all means in parallel for all contrast levels
    % Cast to complex for power operations
    % Numerator means
    ehhv = complex(c + b*(1-e)*o);
    ghhv = complex(c + w*e*o);
    uhhv = rmx * (ehhv.^n)./(ghhv.^n + g0^n);
    
    evhv = complex(a*c + (1-e)*o);
    gvhv = complex(l*c + e*o);
    uvhv = rmx * (evhv.^n)./(gvhv.^n + g0^n);
    
    ehhh = complex(c + (1-e)*o);
    ghhh = complex(c + e*o);
    uhhh = rmx * (ehhh.^n)./(ghhh.^n + g0^n);
    
    evhh = complex(a*c + b*(1-e)*o);
    gvhh = complex(l*c + w*e*o);
    uvhh = rmx * (evhh.^n)./(gvhh.^n + g0^n);
    
    % Denominator means
    ehvv = complex(a*c + b*(1-e)*o);
    ghvv = complex(l*c + w*e*o);
    uhvv = rmx * (ehvv.^n)./(ghvv.^n + g0^n);
    
    evvv = complex(c + (1-e)*o);
    gvvv = complex(c + e*o);
    uvvv = rmx * (evvv.^n)./(gvvv.^n + g0^n);
    
    ehvh = complex(a*c + (1-e)*o);
    ghvh = complex(l*c + e*o);
    uhvh = rmx * (ehvh.^n)./(ghvh.^n + g0^n);
    
    evvh = complex(c + b*(1-e)*o);
    gvvh = complex(c + w*e*o);
    uvvh = rmx * (evvh.^n)./(gvvh.^n + g0^n);
    
    % Take real parts after power operations
    uhhv = real(uhhv); uvhv = real(uvhv); uhhh = real(uhhh); uvhh = real(uvhh);
    uhvv = real(uhvv); uvvv = real(uvvv); uhvh = real(uhvh); uvvh = real(uvvh);
    
    % Pre-allocate trial arrays
    if options.useGPU
        inp1 = gpuArray.randi([0 1], 1, ntrl, 'logical');
        inp2 = gpuArray.randi([0 1], 1, ntrl, 'logical');
        ic = gpuArray.randi([1 nlev], 1, ntrl);
        noise = gpuArray.randn(1, ntrl);
    else
        inp1 = logical(randi([0 1], 1, ntrl));
        inp2 = logical(randi([0 1], 1, ntrl));
        ic = randi([1 nlev], 1, ntrl);
        noise = randn(1, ntrl);
    end
    
    % Use logical indexing for trial types
    type_vv = ~inp1 & ~inp2;
    type_vh = ~inp1 & inp2;
    type_hv = inp1 & ~inp2;
    type_hh = inp1 & inp2;
    
    % Pre-allocate response arrays
    if options.useGPU
        rh = gpuArray.zeros(1, ntrl);
        rv = gpuArray.zeros(1, ntrl);
        output = gpuArray.zeros(nlev, 10);
    else
        rh = zeros(1, ntrl);
        rv = zeros(1, ntrl);
        output = zeros(nlev, 10);
    end
    
    % Store basic parameters
    output(:,1) = c;
    output(:,2) = o;
    
    % Compute responses using logical indexing
    for i = 1:nlev
        mask = ic == i;
        if any(mask)
            rh(mask & type_vv) = noise(mask & type_vv) + uhvv(i);
            rv(mask & type_vv) = noise(mask & type_vv) + uvvv(i);
            
            rh(mask & type_vh) = noise(mask & type_vh) + uhvh(i);
            rv(mask & type_vh) = noise(mask & type_vh) + uvvh(i);
            
            rh(mask & type_hv) = noise(mask & type_hv) + uhhv(i);
            rv(mask & type_hv) = noise(mask & type_hv) + uvhv(i);
            
            rh(mask & type_hh) = noise(mask & type_hh) + uhhh(i);
            rv(mask & type_hh) = noise(mask & type_hh) + uvhh(i);
        end
    end
    
    % Process all trials
    for i = 1:nlev
        trial_mask = ic == i;
        if any(trial_mask)
            % Calculate likelihood ratios vectorized
            num = sum(exp(-0.5 * ((rh(trial_mask)-uhhv(i)).^2 + (rv(trial_mask)-uvhv(i)).^2)) + ...
                     exp(-0.5 * ((rh(trial_mask)-uhhh(i)).^2 + (rv(trial_mask)-uvhh(i)).^2)));
            den = sum(exp(-0.5 * ((rh(trial_mask)-uhvv(i)).^2 + (rv(trial_mask)-uvvv(i)).^2)) + ...
                     exp(-0.5 * ((rh(trial_mask)-uhvh(i)).^2 + (rv(trial_mask)-uvvh(i)).^2)));
            
            lr = num/den;
            hvrsp = lr > 1.0 | (lr == 1.0 & rand() > 0.5);
            
            % Count trials and correct responses
            mask_hh = trial_mask & type_hh;
            mask_hv = trial_mask & type_hv;
            mask_vh = trial_mask & type_vh;
            mask_vv = trial_mask & type_vv;
            
            output(i,3) = sum(mask_hh);
            output(i,4) = sum(mask_hv);
            output(i,5) = sum(mask_vh);
            output(i,6) = sum(mask_vv);
            
            output(i,7) = sum(mask_hh & (inp1 == hvrsp));
            output(i,8) = sum(mask_hv & (inp1 == hvrsp));
            output(i,9) = sum(mask_vh & (inp1 == hvrsp));
            output(i,10) = sum(mask_vv & (inp1 == hvrsp));
        end
    end
    
    % Move back to CPU if needed
    if options.useGPU
        output = gather(output);
        c = gather(c);
    end
    
    % Calculate probabilities
    pc = (output(:,7) + output(:,10)) ./ (output(:,3) + output(:,6));
    pic = (output(:,8) + output(:,9)) ./ (output(:,4) + output(:,5));
    
    % Format output
    if strcmp(optoStr, 'baseline')
        pcntrl = (output(:,7) + output(:,8) + output(:,9) + output(:,10)) ./ ...
                 (output(:,3) + output(:,4) + output(:,5) + output(:,6));
        
        zeroIdx = find(c == 0);
        inconIdx = find(c < 0);
        conIdx = find(c > 0);
        
        zeroPart = rmnan(mean([100 - (pcntrl(zeroIdx) * 100), pcntrl(zeroIdx) * 100]))';
        
        yPredicted = [100 - (pcntrl(inconIdx) * 100); ...
                      zeroPart; ...
                      pcntrl(conIdx) * 100]';
    else
        zeroIdx = find(c == 0);
        inconIdx = find(c < 0);
        conIdx = find(c > 0);
        
        zeroPart = rmnan(mean([100 - (pic(zeroIdx) * 100), pc(zeroIdx) * 100]))';
        
        yPredicted = [100 - (pic(inconIdx) * 100); ...
                      zeroPart; ...
                      pc(conIdx) * 100]';
    end
end

function y = rmnan(x)
    y = x;
    y(isnan(y)) = [];
end