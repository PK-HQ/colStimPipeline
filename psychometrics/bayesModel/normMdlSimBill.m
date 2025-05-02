function yPredicted = normMdlSimBill(x, params, optoStr, options)
    % v_o_psychometric_func_baseline - Minimal adaptation of v_o_psychometric_3.m
    % Provides identical input/output interface as normMdlSimOptimized while
    % preserving the original computation logic
    %
    % Inputs:
    %   x - Stimulus contrast levels
    %   params - Model parameters [a, b, l, w, g0, n, rmx, e, o]
    %   optoStr - Mode ('baseline' or 'opto')
    %   options - (optional) Structure with customization options:
    %       .nTrials - Number of trials per contrast level (default: 100000)
    % 
    % Outputs:
    %   yPredicted - Predicted response values
    
    % Set default options
    if nargin < 4
        options = struct();
    end
    
    if ~isfield(options, 'nTrials')
        options.nTrials = 100000; % Default from original script
    end
    
    % Extract parameters from input array
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
        o = 0;  % Zero opto-stim effective contrast for baseline
    else
        o = params(9);  % Opto-stim effective contrast for opto condition
    end
    
    % Setup contrasts - use the passed contrast levels
    c = x(:)';        % Ensure row vector (original used row vectors)
    nlev = length(c); % Number of contrast levels
    
    % Set number of trials
    ntrl = options.nTrials;
    
    % From here, the code is almost identical to the original v_o_psychometric_3.m
    % Only the variable initialization and final output formatting are changed
    
    % numerator mean storage
    uhhv = zeros(1,nlev);
    uvhv = zeros(1,nlev);
    uhhh = zeros(1,nlev);
    uvhh = zeros(1,nlev);
    % denominator mean storage
    uhvv = zeros(1,nlev);
    uvvv = zeros(1,nlev);
    uhvh = zeros(1,nlev);
    uvvh = zeros(1,nlev);
    
    for i = 1:nlev
        % numerator means  
        ehhv = c(i) + b*(1-e)*o;
        ghhv = c(i) + w*e*o;
        uhhv(i) = rmx * (ehhv^n)/(ghhv^n + g0^n);
        %
        evhv = a*c(i) + (1-e)*o;
        gvhv = l*c(i) + e*o;
        uvhv(i) = rmx * (evhv^n)/(gvhv^n + g0^n);
        %
        ehhh = c(i) + (1-e)*o;
        ghhh = c(i) + e*o;
        uhhh(i) = rmx * (ehhh^n)/(ghhh^n + g0^n);
        %
        evhh = a*c(i) + b*(1-e)*o;
        gvhh = l*c(i) + w*e*o;
        uvhh(i) = rmx * (evhh^n)/(gvhh^n + g0^n);
        % denominator means
        ehvv = a*c(i) + b*(1-e)*o;
        ghvv = l*c(i) + w*e*o;
        uhvv(i) = rmx * (ehvv^n)/(ghvv^n + g0^n);
        %
        evvv = c(i) + (1-e)*o;
        gvvv = c(i) + e*o;
        uvvv(i) = rmx * (evvv^n)/(gvvv^n + g0^n);
        %
        ehvh = a*c(i) + (1-e)*o;
        ghvh = l*c(i) + e*o;
        uhvh(i) = rmx * (ehvh^n)/(ghvh^n + g0^n);
        %
        evvh = c(i) + b*(1-e)*o;
        gvvh = c(i) + w*e*o;
        uvvh(i) = rmx * (evvh^n)/(gvvh^n + g0^n);
    end
    
    % run trials
    trls = zeros(ntrl,4); % inp(1),inp(2),ic,hvrsp
    for j = 1:ntrl
        inp = randi([0 1],1,2,"logical");
        ic = randi([1 nlev]);   % index to contrast
        trls(j,1:2) = inp;
        trls(j,3) = ic;
        if inp(1) == 0 && inp(2) == 0     % vert stm vert optstm
          rh = randn() + uhvv(ic);
          rv = randn() + uvvv(ic);
        elseif inp(1) == 0 && inp(2) == 1 % vert stm horz optstm
          rh = randn() + uhvh(ic);
          rv = randn() + uvvh(ic);
        elseif inp(1) == 1 && inp(2) == 0 % horz stm vert optstm
          rh = randn() + uhhv(ic);
          rv = randn() + uvhv(ic);      
        else                              % horz stm horz optstm
          rh = randn() + uhhh(ic);
          rv = randn() + uvhh(ic);  
        end
        num = sum(exp( -0.5*( (rh-uhhv).^2 + (rv-uvhv).^2 ) )...
            + exp( -0.5*( (rh-uhhh).^2 + (rv-uvhh).^2 ) ));
        den = sum(exp( -0.5*( (rh-uhvv).^2 + (rv-uvvv).^2 ) )...
            + exp( -0.5*( (rh-uhvh).^2 + (rv-uvvh).^2 ) ));
        lr = num/den;
        if lr > 1.0     % horz behav response
          trls(j,4) = 1;    
        elseif lr == 1  % random behav response
          if rand() > 0.5
             trls(j,4) = 1;
          end
        end             % else leave vert behav response
    end
    
    % process trials and save results
    % trls = inp(1),inp(2),ic,hvrsp
    % output = c,o,nhh,nhv,nvh,nvv,chh,chv,cvh,cvv
    output = zeros(nlev,10);
    for j = 1:ntrl
      ic = trls(j,3);
      output(ic,1) = c(ic);
      output(ic,2) = o;
      if trls(j,1) == 1 && trls(j,2) == 1
        output(ic,3) = output(ic,3) + 1;    % hh trial count
        if trls(j,1) == trls(j,4)
          output(ic,7) = output(ic,7) + 1;  % hh correct response count
        end
      elseif trls(j,1) == 1 && trls(j,2) == 0
        output(ic,4) = output(ic,4) + 1;    % hv trial count
        if trls(j,1) == trls(j,4)
          output(ic,8) = output(ic,8) + 1;  % hv correct response count
        end
      elseif trls(j,1) == 0 && trls(j,2) == 1
        output(ic,5) = output(ic,5) + 1;    % vh trial count
        if trls(j,1) == trls(j,4)
          output(ic,9) = output(ic,9) + 1;  % vh correct response count
        end
      elseif trls(j,1) == 0 && trls(j,2) == 0
        output(ic,6) = output(ic,6) + 1;    % vv trial count
        if trls(j,1) == trls(j,4)
          output(ic,10) = output(ic,10) + 1;% vv correct response count
        end
      end
    end
    
    % Calculate probabilities as in original
    pc = (output(:,7) + output(:,10)) ./ (output(:,3) + output(:,6));
    pic = (output(:,8) + output(:,9)) ./ (output(:,4) + output(:,5));
    
    % Format output matching normMdlSimOptimized interface
    if strcmp(optoStr, 'baseline')
        % Calculate all probabilities as in original
        pcntrl = (output(:,7) + output(:,8) + output(:,9) + output(:,10)) ./ ...
                 (output(:,3) + output(:,4) + output(:,5) + output(:,6));
        
        % Format output exactly like normMdlSimOptimized
        zeroIdx = find(c == 0);
        inconIdx = find(c < 0);
        conIdx = find(c > 0);
        
        % Handle zero contrast (use mean for symmetry as in original)
        zeroPart = rmnan(mean([100 - (pcntrl(zeroIdx) * 100), pcntrl(zeroIdx) * 100]))';
        
        % Final format
        yPredicted = [100 - (pcntrl(inconIdx) * 100); ...
                      zeroPart; ...
                      pcntrl(conIdx) * 100]';
    else
        % For opto condition
        zeroIdx = find(c == 0);
        inconIdx = find(c < 0);
        conIdx = find(c > 0);
        
        % Handle zero contrast
        zeroPart = rmnan(mean([100 - (pic(zeroIdx) * 100), pc(zeroIdx) * 100]))';
        
        % Final format
        yPredicted = [100 - (pic(inconIdx) * 100); ...
                      zeroPart; ...
                      pc(conIdx) * 100]';
    end
end

% Helper function to remove NaN values (copied from original)
function y = rmnan(x)
    y = x;
    y(isnan(y)) = [];
end