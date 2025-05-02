function yPredicted = v_o_psychometricFn(x, params, options)
% v_o_psychometric Simulates visual + opto-stim psychometric functions.
%
%   yPredicted = v_o_psychometric(x, params, options) 
%
%   Inputs:
%       x       - vector of contrast values.
%       params  - parameter vector: [a, b, l, w, g0, n, rmx, e, o].
%                 where a, b = excitatory parameters,
%                       l, w = normalization parameters,
%                       g0 = normalization constant,
%                       n = spiking exponent,
%                       rmx = maximum response,
%                       e = weight for opto-stim influence,
%                       o = effective opto-stim contrast.
%       options - struct with field 'nTrials' (number of trials to simulate).
%
%   Output:
%       yPredicted - struct with fields:
%                      pc      - percent correct for congruent conditions,
%                      pic     - percent correct for incongruent conditions,
%                      pcntrl  - percent correct for the control (no opto-stim) condition.
%                    (Each field is multiplied by 100.)
%
%   The function simulates responses for both a stimulation condition (with o from params)
%   and a control condition (o = 0) using random trial selection. The response on each trial
%   is determined by a likelihood ratio computed from noisy draws around the modeled means.

% Unpack parameters
a   = params(1);
b   = params(2);
l   = params(3);
w   = params(4);
g0  = params(5);
n   = params(6);
rmx = params(7);
e   = params(8);
o_stim = params(9);  % effective opto-stim contrast in the stimulation condition

% Number of trials and number of contrast levels
ntrl = options.nTrials;
nlev = length(x);

%% STIMULATION CONDITION (with o = o_stim)
% Preallocate mean-response arrays for numerator and denominator components
uhhv = zeros(1, nlev); uvhv = zeros(1, nlev);
uhhh = zeros(1, nlev); uvhh = zeros(1, nlev);
uhvv = zeros(1, nlev); uvvv = zeros(1, nlev);
uhvh = zeros(1, nlev); uvvh = zeros(1, nlev);

for i = 1:nlev
    c_val = x(i);
    % Numerator means
    ehhv = c_val + b*(1 - e)*o_stim;
    ghhv = c_val + w*e*o_stim;
    uhhv(i) = rmx * (ehhv^n) / (ghhv^n + g0^n);
    
    evhv = a * c_val + (1 - e)*o_stim;
    gvhv = l * c_val + e*o_stim;
    uvhv(i) = rmx * (evhv^n) / (gvhv^n + g0^n);
    
    ehhh = c_val + (1 - e)*o_stim;
    ghhh = c_val + e*o_stim;
    uhhh(i) = rmx * (ehhh^n) / (ghhh^n + g0^n);
    
    evhh = a * c_val + b*(1 - e)*o_stim;
    gvhh = l * c_val + w*e*o_stim;
    uvhh(i) = rmx * (evhh^n) / (gvhh^n + g0^n);
    
    % Denominator means
    ehvv = a * c_val + b*(1 - e)*o_stim;
    ghvv = l * c_val + w*e*o_stim;
    uhvv(i) = rmx * (ehvv^n) / (ghvv^n + g0^n);
    
    evvv = c_val + (1 - e)*o_stim;
    gvvv = c_val + e*o_stim;
    uvvv(i) = rmx * (evvv^n) / (gvvv^n + g0^n);
    
    ehvh = a * c_val + (1 - e)*o_stim;
    ghvh = l * c_val + e*o_stim;
    uhvh(i) = rmx * (ehvh^n) / (ghvh^n + g0^n);
    
    evvh = c_val + b*(1 - e)*o_stim;
    gvvh = c_val + w*e*o_stim;
    uvvh(i) = rmx * (evvh^n) / (gvvh^n + g0^n);
end

% Run trials for stimulation condition
trls = zeros(ntrl, 4);  % columns: [input1, input2, contrast index, behavioral response]
for j = 1:ntrl
    % Generate random binary stimulus conditions for two inputs
    inp = randi([0 1], 1, 2);
    ic  = randi([1, nlev]);  % randomly choose one contrast level
    trls(j, 1:2) = inp;
    trls(j, 3) = ic;
    
    % Depending on condition, assign response means
    if (inp(1) == 0 && inp(2) == 0)         % vertical stim, vertical opto
        rh = randn() + uhvv(ic);
        rv = randn() + uvvv(ic);
    elseif (inp(1) == 0 && inp(2) == 1)         % vertical stim, horizontal opto
        rh = randn() + uhvh(ic);
        rv = randn() + uvvh(ic);
    elseif (inp(1) == 1 && inp(2) == 0)         % horizontal stim, vertical opto
        rh = randn() + uhhv(ic);
        rv = randn() + uvhv(ic);
    else                                        % horizontal stim, horizontal opto
        rh = randn() + uhhh(ic);
        rv = randn() + uvhh(ic);
    end
    
    % Compute likelihood ratio using the modeled responses
    num = sum(exp(-0.5 * ((rh - uhhv).^2 + (rv - uvhv).^2)) + ...
              exp(-0.5 * ((rh - uhhh).^2 + (rv - uvhh).^2)));
    den = sum(exp(-0.5 * ((rh - uhvv).^2 + (rv - uvvv).^2)) + ...
              exp(-0.5 * ((rh - uhvh).^2 + (rv - uvvh).^2)));
    lr = num / den;
    
    % Determine behavioral response based on likelihood ratio
    if lr > 1.0
        trls(j, 4) = 1;
    elseif lr == 1
        if rand() > 0.5
            trls(j, 4) = 1;
        else
            trls(j, 4) = 0;
        end
    else
        trls(j, 4) = 0;
    end
end

% Process stimulation trials: accumulate counts for each contrast level.
outputStim = zeros(nlev, 10);  
% Columns: 1 = contrast, 2 = o, 3 = hh trial count, 4 = hv trial count,
% 5 = vh trial count, 6 = vv trial count, 7 = hh correct, 8 = hv correct,
% 9 = vh correct, 10 = vv correct.
for j = 1:ntrl
    ic = trls(j, 3);
    outputStim(ic, 1) = x(ic);
    outputStim(ic, 2) = o_stim;
    if (trls(j,1) == 1 && trls(j,2) == 1)  % hh condition
        outputStim(ic, 3) = outputStim(ic, 3) + 1;
        if trls(j,1) == trls(j,4)
            outputStim(ic, 7) = outputStim(ic, 7) + 1;
        end
    elseif (trls(j,1) == 1 && trls(j,2) == 0)  % hv condition
        outputStim(ic, 4) = outputStim(ic, 4) + 1;
        if trls(j,1) == trls(j,4)
            outputStim(ic, 8) = outputStim(ic, 8) + 1;
        end
    elseif (trls(j,1) == 0 && trls(j,2) == 1)  % vh condition
        outputStim(ic, 5) = outputStim(ic, 5) + 1;
        if trls(j,1) == trls(j,4)
            outputStim(ic, 9) = outputStim(ic, 9) + 1;
        end
    elseif (trls(j,1) == 0 && trls(j,2) == 0)  % vv condition
        outputStim(ic, 6) = outputStim(ic, 6) + 1;
        if trls(j,1) == trls(j,4)
            outputStim(ic, 10) = outputStim(ic, 10) + 1;
        end
    end
end

% Calculate psychometric functions for stimulation condition
pc  = (outputStim(:,7) + outputStim(:,10)) ./ (outputStim(:,3) + outputStim(:,6));
pic = (outputStim(:,8) + outputStim(:,9)) ./ (outputStim(:,4) + outputStim(:,5));

%% CONTROL CONDITION (with o = 0)
o_ctrl = 0;
% Preallocate arrays for control condition
uhhv = zeros(1, nlev); uvhv = zeros(1, nlev);
uhhh = zeros(1, nlev); uvhh = zeros(1, nlev);
uhvv = zeros(1, nlev); uvvv = zeros(1, nlev);
uhvh = zeros(1, nlev); uvvh = zeros(1, nlev);

for i = 1:nlev
    c_val = x(i);
    % Numerator means (control, so o = 0)
    ehhv = c_val + b*o_ctrl;
    ghhv = c_val + w*o_ctrl;
    uhhv(i) = rmx * (ehhv^n) / (ghhv^n + g0^n);
    
    evhv = a*c_val + o_ctrl;
    gvhv = l*c_val + o_ctrl;
    uvhv(i) = rmx * (evhv^n) / (gvhv^n + g0^n);
    
    ehhh = c_val + o_ctrl;
    ghhh = c_val + o_ctrl;
    uhhh(i) = rmx * (ehhh^n) / (ghhh^n + g0^n);
    
    evhh = a*c_val + b*o_ctrl;
    gvhh = l*c_val + w*o_ctrl;
    uvhh(i) = rmx * (evhh^n) / (gvhh^n + g0^n);
    
    % Denominator means
    ehvv = a*c_val + b*o_ctrl;
    ghvv = l*c_val + w*o_ctrl;
    uhvv(i) = rmx * (ehvv^n) / (ghvv^n + g0^n);
    
    evvv = c_val + o_ctrl;
    gvvv = c_val + o_ctrl;
    uvvv(i) = rmx * (evvv^n) / (gvvv^n + g0^n);
    
    ehvh = a*c_val + o_ctrl;
    ghvh = l*c_val + o_ctrl;
    uhvh(i) = rmx * (ehvh^n) / (ghvh^n + g0^n);
    
    evvh = c_val + b*o_ctrl;
    gvvh = c_val + w*o_ctrl;
    uvvh(i) = rmx * (evvh^n) / (gvvh^n + g0^n);
end

% Run trials for control condition
trls = zeros(ntrl, 4);
for j = 1:ntrl
    inp = randi([0 1], 1, 2);
    ic  = randi([1, nlev]);
    trls(j, 1:2) = inp;
    trls(j, 3) = ic;
    
    if (inp(1) == 0 && inp(2) == 0)
        rh = randn() + uhvv(ic);
        rv = randn() + uvvv(ic);
    elseif (inp(1) == 0 && inp(2) == 1)
        rh = randn() + uhvh(ic);
        rv = randn() + uvvh(ic);
    elseif (inp(1) == 1 && inp(2) == 0)
        rh = randn() + uhhv(ic);
        rv = randn() + uvhv(ic);
    else
        rh = randn() + uhhh(ic);
        rv = randn() + uvhh(ic);
    end
    
    num = sum(exp(-0.5 * ((rh - uhhv).^2 + (rv - uvhv).^2)) + ...
              exp(-0.5 * ((rh - uhhh).^2 + (rv - uvhh).^2)));
    den = sum(exp(-0.5 * ((rh - uhvv).^2 + (rv - uvvv).^2)) + ...
              exp(-0.5 * ((rh - uhvh).^2 + (rv - uvvh).^2)));
    lr = num / den;
    
    if lr > 1.0
        trls(j, 4) = 1;
    elseif lr == 1
        if rand() > 0.5
            trls(j, 4) = 1;
        else
            trls(j, 4) = 0;
        end
    else
        trls(j, 4) = 0;
    end
end

% Process control trials
outputCtrl = zeros(nlev, 10);
for j = 1:ntrl
    ic = trls(j, 3);
    outputCtrl(ic, 1) = x(ic);
    outputCtrl(ic, 2) = o_ctrl;
    if (trls(j,1) == 1 && trls(j,2) == 1)
        outputCtrl(ic, 3) = outputCtrl(ic, 3) + 1;
        if trls(j,1) == trls(j,4)
            outputCtrl(ic, 7) = outputCtrl(ic, 7) + 1;
        end
    elseif (trls(j,1) == 1 && trls(j,2) == 0)
        outputCtrl(ic, 4) = outputCtrl(ic, 4) + 1;
        if trls(j,1) == trls(j,4)
            outputCtrl(ic, 8) = outputCtrl(ic, 8) + 1;
        end
    elseif (trls(j,1) == 0 && trls(j,2) == 1)
        outputCtrl(ic, 5) = outputCtrl(ic, 5) + 1;
        if trls(j,1) == trls(j,4)
            outputCtrl(ic, 9) = outputCtrl(ic, 9) + 1;
        end
    elseif (trls(j,1) == 0 && trls(j,2) == 0)
        outputCtrl(ic, 6) = outputCtrl(ic, 6) + 1;
        if trls(j,1) == trls(j,4)
            outputCtrl(ic, 10) = outputCtrl(ic, 10) + 1;
        end
    end
end

% Calculate control psychometric measure
pcntrl = (outputCtrl(:,7) + outputCtrl(:,8) + ...
          outputCtrl(:,9) + outputCtrl(:,10)) ./ ...
          (outputCtrl(:,3) + outputCtrl(:,4) + outputCtrl(:,5) + outputCtrl(:,6));

%% Output structure (all measures multiplied by 100)
yPredicted.pc = pc * 100;
yPredicted.pic = pic * 100;
yPredicted.pcntrl = pcntrl * 100;

end
