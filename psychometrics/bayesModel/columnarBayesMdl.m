function yPredicted = columnarBayesMdl(x, params, options)
% columnarBayesMdl: Accelerated simulation of visual + opto-stim psychometric functions (from Bill's v_o_psychometric_3)
%
%   yPredicted = columnarBayesMdl(x, params, options)
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
%   This accelerated version vectorizes all heavy computations and uses GPU
%   acceleration when available (if Parallel Computing Toolbox is installed).
%   The logic, indexing, and calculations are preserved to within 1–2 decimal 
%   places of the original implementation.

%% Setup GPU (if available)
useGPU = (gpuDeviceCount > 0);
if useGPU
    % Convert x to a gpuArray.
    x = gpuArray(x);
end

% Unpack parameters
a      = params(1);
b      = params(2);
l      = params(3);
w      = params(4);
g0     = params(5);
n      = params(6);
rmx    = params(7);
e      = params(8);
o_stim = params(9);  % effective opto-stim contrast for stimulation condition

ntrl = options.nTrials;    % number of trials
nlev = length(x);          % number of contrast levels
g0_n = g0^n;               % precompute g0^n

%% --- STIMULATION CONDITION (with o = o_stim) ---
% Compute the modeled means for each contrast (vectorized):
% Numerator means
uhhv = rmx * ((x + b*(1-e)*o_stim).^n) ./ (((x + w*e*o_stim).^n) + g0_n);
uvhv = rmx * ((a*x + (1-e)*o_stim).^n) ./ (((l*x + e*o_stim).^n) + g0_n);
uhhh = rmx * ((x + (1-e)*o_stim).^n) ./ (((x + e*o_stim).^n) + g0_n);
uvhh = rmx * ((a*x + b*(1-e)*o_stim).^n) ./ (((l*x + w*e*o_stim).^n) + g0_n);
% Denominator means
uhvv = rmx * ((a*x + b*(1-e)*o_stim).^n) ./ (((l*x + w*e*o_stim).^n) + g0_n);
uvvv = rmx * ((x + (1-e)*o_stim).^n) ./ (((x + e*o_stim).^n) + g0_n);
uhvh = rmx * ((a*x + (1-e)*o_stim).^n) ./ (((l*x + e*o_stim).^n) + g0_n);
uvvh = rmx * ((x + b*(1-e)*o_stim).^n) ./ (((x + w*e*o_stim).^n) + g0_n);

% --- Simulate Trials (vectorized) ---
% Generate all random binary stimulus conditions and contrast indices at once.
inp = randi([0, 1], ntrl, 2);      % each row: [stimulus, opto]
ic  = randi([1, nlev], ntrl, 1);    % contrast level index for each trial

% Preallocate response arrays (rh and rv)
rh = zeros(ntrl,1, 'like', x);
rv = zeros(ntrl,1, 'like', x);

% Determine trial conditions using logical indices.
cond1 = (inp(:,1)==0 & inp(:,2)==0);   % (0,0): vertical stim, vertical opto -> use (uhvv, uvvv)
cond2 = (inp(:,1)==0 & inp(:,2)==1);   % (0,1): vertical stim, horizontal opto -> use (uhvh, uvvh)
cond3 = (inp(:,1)==1 & inp(:,2)==0);   % (1,0): horizontal stim, vertical opto -> use (uhhv, uvhv)
cond4 = (inp(:,1)==1 & inp(:,2)==1);   % (1,1): horizontal stim, horizontal opto -> use (uhhh, uvhh)

% For each condition, add a random draw (randn) to the corresponding precomputed mean.
if any(cond1)
    rh(cond1) = randn(sum(cond1),1, 'like', x) + uhvv(ic(cond1))';
    rv(cond1) = randn(sum(cond1),1, 'like', x) + uvvv(ic(cond1))';
end
if any(cond2)
    rh(cond2) = randn(sum(cond2),1, 'like', x) + uhvh(ic(cond2))';
    rv(cond2) = randn(sum(cond2),1, 'like', x) + uvvh(ic(cond2))';
end
if any(cond3)
    rh(cond3) = randn(sum(cond3),1, 'like', x) + uhhv(ic(cond3))';
    rv(cond3) = randn(sum(cond3),1, 'like', x) + uvhv(ic(cond3))';
end
if any(cond4)
    rh(cond4) = randn(sum(cond4),1, 'like', x) + uhhh(ic(cond4))';
    rv(cond4) = randn(sum(cond4),1, 'like', x) + uvhh(ic(cond4))';
end

% --- Vectorized Likelihood Ratio Computation ---
% For each trial, compute differences between its response and the full vector of means.
% (Implicit expansion creates an ntrl x nlev matrix for each term.)
num_term1 = exp(-0.5 * ( (rh - uhhv).^2 + (rv - uvhv).^2 ));
num_term2 = exp(-0.5 * ( (rh - uhhh).^2 + (rv - uvhh).^2 ));
num_val   = sum(num_term1 + num_term2, 2);

den_term1 = exp(-0.5 * ( (rh - uhvv).^2 + (rv - uvvv).^2 ));
den_term2 = exp(-0.5 * ( (rh - uhvh).^2 + (rv - uvvh).^2 ));
den_val   = sum(den_term1 + den_term2, 2);

lr = num_val ./ den_val;

% --- Determine Behavioral Response ---
% If likelihood ratio > 1, response = 1; if equal to 1, use a coin toss; else response = 0.
resp = zeros(ntrl, 1, 'like', x);
resp(lr > 1) = 1;
equal_idx = (lr == 1);
if any(equal_idx)
    coin = rand(sum(equal_idx), 1, 'like', x);
    resp(equal_idx) = coin > 0.5;
end

%% --- Accumulate Trial Counts for Stimulation ---
% For accumulation we gather the data to CPU if computed on GPU.
if useGPU
    inp = gather(inp);
    ic  = gather(ic);
    resp = gather(resp);
    x_cpu = gather(x);
else
    x_cpu = x;
end

% Preallocate output matrix (nlev x 10)
% Columns: 1 = contrast, 2 = o, 3 = hh count, 4 = hv count, 5 = vh count,
% 6 = vv count, 7 = hh correct, 8 = hv correct, 9 = vh correct, 10 = vv correct.
outputStim = zeros(nlev, 10);
outputStim(:,1) = x_cpu(:);
outputStim(:,2) = o_stim;

% Identify trial types based on the two binary inputs.
hh_idx = (inp(:,1)==1 & inp(:,2)==1);
hv_idx = (inp(:,1)==1 & inp(:,2)==0);
vh_idx = (inp(:,1)==0 & inp(:,2)==1);
vv_idx = (inp(:,1)==0 & inp(:,2)==0);

% Use accumarray to count trials and correct responses.
outputStim(:,3) = accumarray(ic(hh_idx), 1, [nlev,1], @sum, 0);
outputStim(:,7) = accumarray(ic(hh_idx & (resp==1)), 1, [nlev,1], @sum, 0);

outputStim(:,4) = accumarray(ic(hv_idx), 1, [nlev,1], @sum, 0);
outputStim(:,8) = accumarray(ic(hv_idx & (resp==1)), 1, [nlev,1], @sum, 0);

outputStim(:,5) = accumarray(ic(vh_idx), 1, [nlev,1], @sum, 0);
outputStim(:,9) = accumarray(ic(vh_idx & (resp==0)), 1, [nlev,1], @sum, 0);

outputStim(:,6) = accumarray(ic(vv_idx), 1, [nlev,1], @sum, 0);
outputStim(:,10)= accumarray(ic(vv_idx & (resp==0)), 1, [nlev,1], @sum, 0);

% Compute psychometric measures for stimulation condition.
pc  = (outputStim(:,7) + outputStim(:,10)) ./ (outputStim(:,3) + outputStim(:,6));
pic = (outputStim(:,8) + outputStim(:,9)) ./ (outputStim(:,4) + outputStim(:,5));

%% --- CONTROL CONDITION (with o = 0) ---
o_ctrl = 0;

% Compute control means (with o = 0; note many terms simplify)
uhhv_ctrl = rmx * (x.^n) ./ (x.^n + g0_n);
uvhv_ctrl = rmx * ((a*x).^n) ./ (((l*x).^n) + g0_n);
uhhh_ctrl = rmx * (x.^n) ./ (x.^n + g0_n);
uvhh_ctrl = rmx * ((a*x).^n) ./ (((l*x).^n) + g0_n);

% For denominator, note:
uhvv_ctrl = uvhh_ctrl; 
uvvv_ctrl = uhhv_ctrl;
uhvh_ctrl = uvhv_ctrl;
uvvh_ctrl = uhhv_ctrl;

% Generate control condition trials.
inp_ctrl = randi([0,1], ntrl, 2);
ic_ctrl  = randi([1, nlev], ntrl, 1);

rh_ctrl = zeros(ntrl,1, 'like', x);
rv_ctrl = zeros(ntrl,1, 'like', x);

cond1_ctrl = (inp_ctrl(:,1)==0 & inp_ctrl(:,2)==0);
cond2_ctrl = (inp_ctrl(:,1)==0 & inp_ctrl(:,2)==1);
cond3_ctrl = (inp_ctrl(:,1)==1 & inp_ctrl(:,2)==0);
cond4_ctrl = (inp_ctrl(:,1)==1 & inp_ctrl(:,2)==1);

if any(cond1_ctrl)
    rh_ctrl(cond1_ctrl) = randn(sum(cond1_ctrl),1, 'like', x) + uhvv_ctrl(ic_ctrl(cond1_ctrl))';
    rv_ctrl(cond1_ctrl) = randn(sum(cond1_ctrl),1, 'like', x) + uvvv_ctrl(ic_ctrl(cond1_ctrl))';
end
if any(cond2_ctrl)
    rh_ctrl(cond2_ctrl) = randn(sum(cond2_ctrl),1, 'like', x) + uhvh_ctrl(ic_ctrl(cond2_ctrl))';
    rv_ctrl(cond2_ctrl) = randn(sum(cond2_ctrl),1, 'like', x) + uvvh_ctrl(ic_ctrl(cond2_ctrl))';
end
if any(cond3_ctrl)
    rh_ctrl(cond3_ctrl) = randn(sum(cond3_ctrl),1, 'like', x) + uhhv_ctrl(ic_ctrl(cond3_ctrl))';
    rv_ctrl(cond3_ctrl) = randn(sum(cond3_ctrl),1, 'like', x) + uvhv_ctrl(ic_ctrl(cond3_ctrl))';
end
if any(cond4_ctrl)
    rh_ctrl(cond4_ctrl) = randn(sum(cond4_ctrl),1, 'like', x) + uhhh_ctrl(ic_ctrl(cond4_ctrl))';
    rv_ctrl(cond4_ctrl) = randn(sum(cond4_ctrl),1, 'like', x) + uvhh_ctrl(ic_ctrl(cond4_ctrl))';
end

% Likelihood ratio for control trials (vectorized):
num_term1_ctrl = exp(-0.5 * ( (rh_ctrl - uhhv_ctrl).^2 + (rv_ctrl - uvhv_ctrl).^2 ));
num_term2_ctrl = exp(-0.5 * ( (rh_ctrl - uhhh_ctrl).^2 + (rv_ctrl - uvhh_ctrl).^2 ));
num_val_ctrl   = sum(num_term1_ctrl + num_term2_ctrl, 2);

den_term1_ctrl = exp(-0.5 * ( (rh_ctrl - uhvv_ctrl).^2 + (rv_ctrl - uvvv_ctrl).^2 ));
den_term2_ctrl = exp(-0.5 * ( (rh_ctrl - uhvh_ctrl).^2 + (rv_ctrl - uvvh_ctrl).^2 ));
den_val_ctrl   = sum(den_term1_ctrl + den_term2_ctrl, 2);

lr_ctrl = num_val_ctrl ./ den_val_ctrl;

% Determine responses for control trials.
resp_ctrl = zeros(ntrl,1, 'like', x);
resp_ctrl(lr_ctrl > 1) = 1;
equal_idx_ctrl = (lr_ctrl == 1);
if any(equal_idx_ctrl)
    coin_ctrl = rand(sum(equal_idx_ctrl), 1, 'like', x);
    resp_ctrl(equal_idx_ctrl) = coin_ctrl > 0.5;
end

% Gather control trial data to CPU for accumulation.
if useGPU
    inp_ctrl = gather(inp_ctrl);
    ic_ctrl  = gather(ic_ctrl);
    resp_ctrl = gather(resp_ctrl);
end

% Accumulate control trial counts.
outputCtrl = zeros(nlev, 10);
outputCtrl(:,1) = x_cpu(:);
outputCtrl(:,2) = o_ctrl;

hh_idx_ctrl = (inp_ctrl(:,1)==1 & inp_ctrl(:,2)==1);
hv_idx_ctrl = (inp_ctrl(:,1)==1 & inp_ctrl(:,2)==0);
vh_idx_ctrl = (inp_ctrl(:,1)==0 & inp_ctrl(:,2)==1);
vv_idx_ctrl = (inp_ctrl(:,1)==0 & inp_ctrl(:,2)==0);

outputCtrl(:,3) = accumarray(ic_ctrl(hh_idx_ctrl), 1, [nlev,1], @sum, 0);
outputCtrl(:,7) = accumarray(ic_ctrl(hh_idx_ctrl & (resp_ctrl==1)), 1, [nlev,1], @sum, 0);

outputCtrl(:,4) = accumarray(ic_ctrl(hv_idx_ctrl), 1, [nlev,1], @sum, 0);
outputCtrl(:,8) = accumarray(ic_ctrl(hv_idx_ctrl & (resp_ctrl==1)), 1, [nlev,1], @sum, 0);

outputCtrl(:,5) = accumarray(ic_ctrl(vh_idx_ctrl), 1, [nlev,1], @sum, 0);
outputCtrl(:,9) = accumarray(ic_ctrl(vh_idx_ctrl & (resp_ctrl==0)), 1, [nlev,1], @sum, 0);

outputCtrl(:,6) = accumarray(ic_ctrl(vv_idx_ctrl), 1, [nlev,1], @sum, 0);
outputCtrl(:,10)= accumarray(ic_ctrl(vv_idx_ctrl & (resp_ctrl==0)), 1, [nlev,1], @sum, 0);

% Compute control psychometric measure.
pcntrl = (outputCtrl(:,7) + outputCtrl(:,8) + outputCtrl(:,9) + outputCtrl(:,10)) ./ ...
          (outputCtrl(:,3) + outputCtrl(:,4) + outputCtrl(:,5) + outputCtrl(:,6));

%% --- Final Output ---
yPredicted.pc      = pc * 100;
yPredicted.pic     = pic * 100;
yPredicted.pcntrl  = pcntrl * 100;

end
