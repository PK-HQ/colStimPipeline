function yPredicted = columnarBayesMdlGPT52(xRow, params, options)
% columnarBayesMdl: Accelerated simulation of visual + opto-stim psychometric functions
%
%   yPredicted = columnarBayesMdl(xRow, params, options)
%
%   Inputs:
%       xRow       - vector of contrast values.
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
%   This accelerated version vectorizes heavy computations, uses GPU acceleration 
%   when available, and further speeds up the likelihood ratio computation by looping
%   over the (typically small) number of contrast levels.

%% Options / reproducibility / GPU
if nargin < 3 || isempty(options), options = struct(); end

% Deterministic objective during fitting (if seed provided)
if isfield(options, 'seed') && ~isempty(options.seed)
    rng(options.seed, 'twister');
    if gpuDeviceCount > 0
        try, gpurng(options.seed, 'Threefry'); end %#ok<TRYNC>
    end
end

% Allow forcing GPU off for reproducibility
useGPU = (gpuDeviceCount > 0);
if isfield(options, 'useGPU')
    useGPU = logical(options.useGPU) && useGPU;
end

% Mode can skip half the work for speed
mode = 'both';
if isfield(options, 'mode') && ~isempty(options.mode)
    mode = options.mode;   % 'stimOnly' | 'ctrlOnly' | 'both'
end

% Preserve original xRow shape for outputs (prevents row/col broadcasting bugs)
xSize = size(xRow);
xRow = xRow(:)';          % always 1×nlev for internal math
if useGPU
    xRow = gpuArray(xRow);
end

%% Parameters
a      = params(1);
b      = params(2);
l      = params(3);
w      = params(4);
g0     = params(5);
n      = params(6);
rmx    = params(7);
e      = params(8);
o_stim = params(9);  % effective opto-stim contrast for stimulation condition

%ntrl = options.nTrials;    % number of trials
nlev = numel(xRow);          % number of contrast levels
g0_n = g0^n;               % precompute g0^n

%% --- STIMULATION CONDITION (with o = o_stim) ---
% Compute the modeled means for each contrast (vectorized):
uhhv = rmx * ((xRow + b*(1-e)*o_stim).^n) ./ (((xRow + w*e*o_stim).^n) + g0_n);
uvhv = rmx * ((a*xRow + (1-e)*o_stim).^n) ./ (((l*xRow + e*o_stim).^n) + g0_n);
uhhh = rmx * ((xRow + (1-e)*o_stim).^n) ./ (((xRow + e*o_stim).^n) + g0_n);
uvhh = rmx * ((a*xRow + b*(1-e)*o_stim).^n) ./ (((l*xRow + w*e*o_stim).^n) + g0_n);
% Denominator means
uhvv = rmx * ((a*xRow + b*(1-e)*o_stim).^n) ./ (((l*xRow + w*e*o_stim).^n) + g0_n);
uvvv = rmx * ((xRow + (1-e)*o_stim).^n) ./ (((xRow + e*o_stim).^n) + g0_n);
uhvh = rmx * ((a*xRow + (1-e)*o_stim).^n) ./ (((l*xRow + e*o_stim).^n) + g0_n);
uvvh = rmx * ((xRow + b*(1-e)*o_stim).^n) ./ (((xRow + w*e*o_stim).^n) + g0_n);

% --- Auto-derive nTrialsPerLevel from trialScale/minTrialsPerLevel if provided ---
% This is critical for plotting, where xRow can be a dense grid.
if (~isfield(options,'nTrialsPerLevel') || isempty(options.nTrialsPerLevel)) && ...
        (isfield(options,'trialScale') || isfield(options,'minTrialsPerLevel'))

    ts = 1;
    if isfield(options,'trialScale') && ~isempty(options.trialScale)
        ts = options.trialScale;
    end

    mmin = 0;
    if isfield(options,'minTrialsPerLevel') && ~isempty(options.minTrialsPerLevel)
        mmin = options.minTrialsPerLevel;
    end

    % Base experimental counts per contrast (your design): 40 at 0%, 20 otherwise
    xHost = gather(xRow);          % safe even if CPU
    baseCounts = 20 * ones(nlev,1);
    baseCounts(xHost(:) == 0) = 40;

    nPerAuto = round(ts * baseCounts);
    nPerAuto = max(nPerAuto, mmin);

    options.nTrialsPerLevel = nPerAuto;   % length nlev
end

% --- Simulate Trials (unbalanced per contrast level if provided) ---
if isfield(options, 'nTrialsPerLevel') && ~isempty(options.nTrialsPerLevel)
    nPer = options.nTrialsPerLevel(:);
    if numel(nPer) ~= nlev
        error('options.nTrialsPerLevel must have length equal to numel(x).');
    end
else
    % fallback: uniform
    if ~isfield(options,'nTrials') || isempty(options.nTrials)
        options.nTrials = 1000;
    end
    nPer = repmat(ceil(options.nTrials / nlev), nlev, 1);
end

ic = repelem((1:nlev)', nPer);
ntrl = numel(ic);

% Generate trial types. Optionally balance HH/HV/VH/VV within each contrast level.
balanceTrialTypes = true;
if isfield(options, 'balanceTrialTypes')
    balanceTrialTypes = logical(options.balanceTrialTypes);
end

if balanceTrialTypes
    inp = localBalancedInpPerLevel(nPer);
else
    inp = randi([0,1], ntrl, 2);
end

% --- Pre-generate trial noise with antithetic variates (variance reduction) ---
nHalf = ceil(ntrl/2);
z = randn(nHalf, 2, 'like', xRow);
Z = [z; -z];
Z = Z(1:ntrl, :);

zH = Z(:,1);
zV = Z(:,2);

rh = zeros(ntrl,1, 'like', xRow);
rv = zeros(ntrl,1, 'like', xRow);

muH = zeros(ntrl,1,'like',xRow);
muV = zeros(ntrl,1,'like',xRow);

% condition masks
cond1 = inp(:,1)==0 & inp(:,2)==0;  % VV
cond2 = inp(:,1)==0 & inp(:,2)==1;  % VH
cond3 = inp(:,1)==1 & inp(:,2)==0;  % HV
cond4 = inp(:,1)==1 & inp(:,2)==1;  % HH

% assign means (NO transposes; force column)
muH(cond1) = uhvv(ic(cond1)).';
muV(cond1) = uvvv(ic(cond1)).';

muH(cond2) = uhvh(ic(cond2)).';
muV(cond2) = uvvh(ic(cond2)).';

% HV
muH(cond3) = uhhv(ic(cond3)).';
muV(cond3) = uvhv(ic(cond3)).';

% HH
muH(cond4) = uhhh(ic(cond4)).';
muV(cond4) = uvhh(ic(cond4)).';

% simulate responses
rh = zH + muH;
rv = zV + muV;

% --- Likelihood Ratio Computation for Stimulation Condition ---
% Instead of forming an ntrl-by-nlev matrix via implicit expansion,
% loop over nlev (which is small) to accumulate the sums.
num_val = zeros(ntrl,1, 'like', rh);
den_val = zeros(ntrl,1, 'like', rh);
for i = 1:nlev
    term1 = exp(-0.5 * ( (rh - uhhv(i)).^2 + (rv - uvhv(i)).^2 ));
    term2 = exp(-0.5 * ( (rh - uhhh(i)).^2 + (rv - uvhh(i)).^2 ));
    num_val = num_val + term1 + term2;
    
    term3 = exp(-0.5 * ( (rh - uhvv(i)).^2 + (rv - uvvv(i)).^2 ));
    term4 = exp(-0.5 * ( (rh - uhvh(i)).^2 + (rv - uvvh(i)).^2 ));
    den_val = den_val + term3 + term4;
end
lr = num_val ./ den_val;

resp = zeros(ntrl, 1, 'like', xRow);
resp(lr > 1) = 1;
equal_idx = (lr == 1);
if any(equal_idx)
    coin = rand(sum(equal_idx), 1, 'like', xRow);
    resp(equal_idx) = coin > 0.5;
end

%% --- Accumulate Trial Counts for Stimulation ---
if useGPU
    inp = gather(inp);
    ic  = gather(ic);
    resp = gather(resp);
    x_cpu = gather(xRow);
else
    x_cpu = xRow;
end

outputStim = zeros(nlev, 10);
outputStim(:,1) = x_cpu(:);
outputStim(:,2) = o_stim;

hh_idx = (inp(:,1)==1 & inp(:,2)==1);
hv_idx = (inp(:,1)==1 & inp(:,2)==0);
vh_idx = (inp(:,1)==0 & inp(:,2)==1);
vv_idx = (inp(:,1)==0 & inp(:,2)==0);

outputStim(:,3) = accumarray(ic(hh_idx), 1, [nlev,1], @sum, 0);
outputStim(:,7) = accumarray(ic(hh_idx & (resp==1)), 1, [nlev,1], @sum, 0);
outputStim(:,4) = accumarray(ic(hv_idx), 1, [nlev,1], @sum, 0);
outputStim(:,8) = accumarray(ic(hv_idx & (resp==1)), 1, [nlev,1], @sum, 0);
outputStim(:,5) = accumarray(ic(vh_idx), 1, [nlev,1], @sum, 0);
outputStim(:,9) = accumarray(ic(vh_idx & (resp==0)), 1, [nlev,1], @sum, 0);
outputStim(:,6) = accumarray(ic(vv_idx), 1, [nlev,1], @sum, 0);
outputStim(:,10)= accumarray(ic(vv_idx & (resp==0)), 1, [nlev,1], @sum, 0);

den_pc  = (outputStim(:,3) + outputStim(:,6));
den_pic = (outputStim(:,4) + outputStim(:,5));

pc  = (outputStim(:,7) + outputStim(:,10)) ./ den_pc;
pic = (outputStim(:,8) + outputStim(:,9))  ./ den_pic;

pc(den_pc==0)   = 0.5;
pic(den_pic==0) = 0.5;

%% --- CONTROL CONDITION (with o = 0) ---
o_ctrl = 0;

uhhv_ctrl = rmx * (xRow.^n) ./ (xRow.^n + g0_n);
uvhv_ctrl = rmx * ((a*xRow).^n) ./ (((l*xRow).^n) + g0_n);
uhhh_ctrl = rmx * (xRow.^n) ./ (xRow.^n + g0_n);
uvhh_ctrl = rmx * ((a*xRow).^n) ./ (((l*xRow).^n) + g0_n);
% For denominator, note:
uhvv_ctrl = uvhh_ctrl; 
uvvv_ctrl = uhhv_ctrl;
uhvh_ctrl = uvhv_ctrl;
uvvh_ctrl = uhhv_ctrl;

% --- CONTROL: use the SAME per-contrast trial counts as stim (reduces jaggedness) ---

% Build ic_ctrl deterministically per contrast level using the same nTrialsPerLevel logic
if isfield(options, 'nTrialsPerLevel') && ~isempty(options.nTrialsPerLevel)
    nPer_ctrl = options.nTrialsPerLevel(:);
    if numel(nPer_ctrl) ~= nlev
        error('options.nTrialsPerLevel must have length equal to numel(xRow).');
    end
else
    % fallback: uniform
    if ~isfield(options,'nTrials') || isempty(options.nTrials)
        options.nTrials = 1000;
    end
    nPer_ctrl = repmat(ceil(options.nTrials / nlev), nlev, 1);
end

ic_ctrl = repelem((1:nlev)', nPer_ctrl);
ntrl = numel(ic_ctrl);

% Control has no opto; keep inp_ctrl(:,2)=0
inp_ctrl = zeros(ntrl, 2);
inp_ctrl(:,1) = randi([0,1], ntrl, 1);  % stimulus only


% --- Antithetic noise: pre-generate once, reuse for all conditions ---
nHalf = ceil(ntrl/2);
z = randn(nHalf, 2, 'like', xRow);
Z = [z; -z];
Z = Z(1:ntrl, :);

zH_ctrl = Z(:,1);
zV_ctrl = Z(:,2);

% --- Allocate means for each trial, then add noise once ---
muH_ctrl = zeros(ntrl, 1, 'like', xRow);
muV_ctrl = zeros(ntrl, 1, 'like', xRow);

cond1_ctrl = (inp_ctrl(:,1)==0 & inp_ctrl(:,2)==0);
cond2_ctrl = (inp_ctrl(:,1)==0 & inp_ctrl(:,2)==1);
cond3_ctrl = (inp_ctrl(:,1)==1 & inp_ctrl(:,2)==0);
cond4_ctrl = (inp_ctrl(:,1)==1 & inp_ctrl(:,2)==1);

if any(cond1_ctrl)
    muH_ctrl(cond1_ctrl) = reshape(uhvv_ctrl(ic_ctrl(cond1_ctrl)), [], 1);
    muV_ctrl(cond1_ctrl) = reshape(uvvv_ctrl(ic_ctrl(cond1_ctrl)), [], 1);
end
if any(cond2_ctrl)
    muH_ctrl(cond2_ctrl) = reshape(uhvh_ctrl(ic_ctrl(cond2_ctrl)), [], 1);
    muV_ctrl(cond2_ctrl) = reshape(uvvh_ctrl(ic_ctrl(cond2_ctrl)), [], 1);
end
if any(cond3_ctrl)
    muH_ctrl(cond3_ctrl) = reshape(uhhv_ctrl(ic_ctrl(cond3_ctrl)), [], 1);
    muV_ctrl(cond3_ctrl) = reshape(uvhv_ctrl(ic_ctrl(cond3_ctrl)), [], 1);
end
if any(cond4_ctrl)
    muH_ctrl(cond4_ctrl) = reshape(uhhh_ctrl(ic_ctrl(cond4_ctrl)), [], 1);
    muV_ctrl(cond4_ctrl) = reshape(uvhh_ctrl(ic_ctrl(cond4_ctrl)), [], 1);
end

% --- Simulate responses (single vectorized add) ---
rh_ctrl = zH_ctrl + muH_ctrl;
rv_ctrl = zV_ctrl + muV_ctrl;

% --- Likelihood Ratio Computation for Control Condition ---
num_val_ctrl = zeros(ntrl,1, 'like', rh_ctrl);
den_val_ctrl = zeros(ntrl,1, 'like', rh_ctrl);
for i = 1:nlev
    term1_ctrl = exp(-0.5 * ( (rh_ctrl - uhhv_ctrl(i)).^2 + (rv_ctrl - uvhv_ctrl(i)).^2 ));
    term2_ctrl = exp(-0.5 * ( (rh_ctrl - uhhh_ctrl(i)).^2 + (rv_ctrl - uvhh_ctrl(i)).^2 ));
    num_val_ctrl = num_val_ctrl + term1_ctrl + term2_ctrl;
    
    term3_ctrl = exp(-0.5 * ( (rh_ctrl - uhvv_ctrl(i)).^2 + (rv_ctrl - uvvv_ctrl(i)).^2 ));
    term4_ctrl = exp(-0.5 * ( (rh_ctrl - uhvh_ctrl(i)).^2 + (rv_ctrl - uvvh_ctrl(i)).^2 ));
    den_val_ctrl = den_val_ctrl + term3_ctrl + term4_ctrl;
end
lr_ctrl = num_val_ctrl ./ den_val_ctrl;

resp_ctrl = zeros(ntrl,1, 'like', xRow);
resp_ctrl(lr_ctrl > 1) = 1;
equal_idx_ctrl = (lr_ctrl == 1);
if any(equal_idx_ctrl)
    coin_ctrl = rand(sum(equal_idx_ctrl), 1, 'like', xRow);
    resp_ctrl(equal_idx_ctrl) = coin_ctrl > 0.5;
end

if useGPU
    inp_ctrl = gather(inp_ctrl);
    ic_ctrl  = gather(ic_ctrl);
    resp_ctrl = gather(resp_ctrl);
end

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

den_ctrl = (outputCtrl(:,3) + outputCtrl(:,4) + outputCtrl(:,5) + outputCtrl(:,6));
pcntrl = (outputCtrl(:,7) + outputCtrl(:,8) + outputCtrl(:,9) + outputCtrl(:,10)) ./ den_ctrl;
pcntrl(den_ctrl==0) = 0.5;

%% --- Final Output ---
% Reshape back to match the input x shape (prevents broadcasting bugs upstream)
yPredicted.pc     = reshape(pc * 100, xSize);
yPredicted.pic    = reshape(pic * 100, xSize);
yPredicted.pcntrl = reshape(pcntrl * 100, xSize);

end

function inp = localBalancedInpPerLevel(nPer)
% Returns N×2 inp with roughly equal counts of the 4 trial types within each contrast level.
% Trial types correspond to:
%   [0 0], [0 1], [1 0], [1 1]
    patterns = [0 0; 0 1; 1 0; 1 1];
    inp = zeros(sum(nPer), 2);
    idx = 1;
    for i = 1:numel(nPer)
        Ni = nPer(i);
        base = floor(Ni/4);
        rem  = Ni - 4*base;

        block = repmat(patterns, base, 1);

        if rem > 0
            ord = randperm(4, rem);
            block = [block; patterns(ord, :)]; %#ok<AGROW>
        end

        % Shuffle within this contrast level
        block = block(randperm(size(block,1)), :);

        inp(idx:idx+Ni-1, :) = block;
        idx = idx + Ni;
    end
end
