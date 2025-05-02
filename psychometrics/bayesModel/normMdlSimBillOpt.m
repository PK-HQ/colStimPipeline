function yPredicted = normMdlSimBillOpt(x, params, optoStr, options)
% normMdlSimBillOpt - Optimized version of normMdlSimBill using GPU vectorization,
% batch processing, and caching of invariant terms.
%
%   This function computes predicted response values (yPredicted) based on 
%   stimulus contrast levels (x) using a Monte-Carlo simulation.
%
%   Inputs:
%       x       - Stimulus contrast levels (vector)
%       params  - Model parameters [a, b, l, w, g0, n, rmx, e, o]
%       optoStr - Mode ('baseline' or 'opto')
%       options - (optional) structure with fields:
%                   .nTrials - Number of trials per contrast level (default: 100000)
%                   .useGPU  - Flag to use GPU vectorization (default: auto-detect)
%                   .verbose - Flag for verbose output (default: false)
%
%   Output:
%       yPredicted - Predicted response values (formatted as in the original)
%
%   This version allows minor numerical differences (within ~2 decimal places) 
%   from the original normMdlSimBill.
%

% ----- Set default options -----
if nargin < 4
    options = struct();
end
if ~isfield(options, 'nTrials')
    options.nTrials = 100000;
end
if ~isfield(options, 'verbose')
    options.verbose = false;
end
if ~isfield(options, 'useGPU')
    try
        options.useGPU = (gpuDeviceCount > 0);
    catch
        options.useGPU = false;
    end
end

if options.verbose
    fprintf('Running normMdlSimBillOpt with %d trials per contrast level. Use GPU: %d\n', options.nTrials, options.useGPU);
end

% ----- Extract parameters -----
a   = params(1);
b   = params(2);
l   = params(3);
w   = params(4);
g0  = params(5);
n   = params(6);
rmx = params(7);
e   = params(8);
if strcmp(optoStr, 'baseline')
    o = 0;
else
    o = params(9);
end

% Compute constant: g0^n
g0n = g0^n;

% ----- Define power function to avoid complex results on GPU -----
% For even exponents, use abs(x).^n. For odd, preserve sign.
if mod(n,2)==0
    power_fn = @(x) (abs(x)).^n;
else
    power_fn = @(x) (sign(x)).*(abs(x)).^n;
end

% ----- Precompute mean responses per contrast (Vectorized & Cached) -----
% Ensure x is a row vector for contrast levels.
c    = x(:)';  
nlev = numel(c);

% Numerator means for the four conditions:
uhhv = rmx * ( power_fn(c + b*(1-e)*o) )   ./ ( power_fn(c + w*e*o) + g0n ); % horizontal stim, vertical opto (visual==1, opto==0)
uvhv = rmx * ( power_fn(a*c + (1-e)*o) )    ./ ( power_fn(l*c + e*o)  + g0n ); % vertical component for above
uhhh = rmx * ( power_fn(c + (1-e)*o) )      ./ ( power_fn(c + e*o)    + g0n ); % horizontal stim, horizontal opto (visual==1, opto==1)
uvhh = rmx * ( power_fn(a*c + b*(1-e)*o) )   ./ ( power_fn(l*c + w*e*o) + g0n ); % vertical component for above

% Denominator means for the four conditions:
uhvv = rmx * ( power_fn(a*c + b*(1-e)*o) )   ./ ( power_fn(l*c + w*e*o) + g0n ); % vertical stim, vertical opto (visual==0, opto==0)
uvvv = rmx * ( power_fn(c + (1-e)*o) )       ./ ( power_fn(c + e*o)    + g0n ); % vertical component for above
uhvh = rmx * ( power_fn(a*c + (1-e)*o) )      ./ ( power_fn(l*c + e*o)    + g0n ); % vertical stim, horizontal opto (visual==0, opto==1)
uvvh = rmx * ( power_fn(c + b*(1-e)*o) )      ./ ( power_fn(c + w*e*o)    + g0n ); % vertical component for above

% Convert all mean arrays to column vectors for consistent indexing
uhhv = uhhv(:);
uvhv = uvhv(:);
uhhh = uhhh(:);
uvhh = uvhh(:);
uhvv = uhvv(:);
uvvv = uvvv(:);
uhvh = uhvh(:);
uvvh = uvvh(:);

% ----- Generate trial data (Using GPU if enabled) -----
ntrl = options.nTrials;
if options.useGPU
    ic     = gpuArray(randi([1, nlev], ntrl, 1));  % Contrast index for each trial
    visual = gpuArray(randi([0, 1], ntrl, 1));       % Visual stimulus (0=vertical, 1=horizontal)
    opto   = gpuArray(randi([0, 1], ntrl, 1));       % Opto stimulus (0=vertical, 1=horizontal)
    noiseH = gpuArray(randn(ntrl, 1));
    noiseV = gpuArray(randn(ntrl, 1));
else
    ic     = randi([1, nlev], ntrl, 1);
    visual = randi([0, 1], ntrl, 1);
    opto   = randi([0, 1], ntrl, 1);
    noiseH = randn(ntrl, 1);
    noiseV = randn(ntrl, 1);
end

% Preallocate response arrays (batch processing)
rh = zeros(ntrl, 1, 'like', noiseH);
rv = zeros(ntrl, 1, 'like', noiseV);

% ----- Vectorized Assignment Based on Conditions -----
% Condition 1: visual==0 and opto==0 (vertical, vertical)
idx = (visual == 0) & (opto == 0);
rh(idx) = noiseH(idx) + uhvv(ic(idx));
rv(idx) = noiseV(idx) + uvvv(ic(idx));

% Condition 2: visual==0 and opto==1 (vertical, horizontal)
idx = (visual == 0) & (opto == 1);
rh(idx) = noiseH(idx) + uhvh(ic(idx));
rv(idx) = noiseV(idx) + uvvh(ic(idx));

% Condition 3: visual==1 and opto==0 (horizontal, vertical)
idx = (visual == 1) & (opto == 0);
rh(idx) = noiseH(idx) + uhhv(ic(idx));
rv(idx) = noiseV(idx) + uvhv(ic(idx));

% Condition 4: visual==1 and opto==1 (horizontal, horizontal)
idx = (visual == 1) & (opto == 1);
rh(idx) = noiseH(idx) + uhhh(ic(idx));
rv(idx) = noiseV(idx) + uvhh(ic(idx));

% ----- Compute Likelihood Ratio (Vectorized with Cached Means) -----
% Retrieve the precomputed means for each trial using the contrast index.
uhhv_trial = uhhv(ic);
uvhv_trial = uvhv(ic);
uhhh_trial = uhhh(ic);
uvhh_trial = uvhh(ic);
uhvv_trial = uhvv(ic);
uvvv_trial = uvvv(ic);
uhvh_trial = uhvh(ic);
uvvh_trial = uvvh(ic);

% Compute numerator and denominator terms in one vectorized step.
termA = exp(-0.5 * ((rh - uhhv_trial).^2 + (rv - uvhv_trial).^2));
termB = exp(-0.5 * ((rh - uhhh_trial).^2 + (rv - uvhh_trial).^2));
num   = termA + termB;

termC = exp(-0.5 * ((rh - uhvv_trial).^2 + (rv - uvvv_trial).^2));
termD = exp(-0.5 * ((rh - uhvh_trial).^2 + (rv - uvvh_trial).^2));
den   = termC + termD;

lr = num ./ den;

% ----- Decision Rule: Determine Trial Outcome -----
% If lr > 1, choose 1 (horizontal); if lr < 1, choose 0 (vertical).
% For lr exactly equal to 1, make a random choice.
trialOutcome = zeros(ntrl, 1, 'like', visual);
trialOutcome(lr > 1) = 1;
tieIdx = (lr == 1);
if any(tieIdx)
    if options.useGPU
        trialOutcome(tieIdx) = gpuArray(double(rand(sum(tieIdx), 1) > 0.5));
    else
        trialOutcome(tieIdx) = double(rand(sum(tieIdx), 1) > 0.5);
    end
end

% ----- Gather GPU Results (if needed) -----
if options.useGPU
    ic           = gather(ic);
    visual       = gather(visual);
    opto         = gather(opto);
    trialOutcome = gather(trialOutcome);
end

% ----- Aggregate Trial Counts and Correct Responses (Vectorized) -----
% Build masks for the four stimulus conditions:
hh_mask = (visual == 1) & (opto == 1);  % horizontal, horizontal
hv_mask = (visual == 1) & (opto == 0);  % horizontal, vertical
vh_mask = (visual == 0) & (opto == 1);  % vertical, horizontal
vv_mask = (visual == 0) & (opto == 0);  % vertical, vertical

% Count trials per condition using accumarray.
hh_counts = accumarray(ic(hh_mask), 1, [nlev, 1], @sum, 0);
hv_counts = accumarray(ic(hv_mask), 1, [nlev, 1], @sum, 0);
vh_counts = accumarray(ic(vh_mask), 1, [nlev, 1], @sum, 0);
vv_counts = accumarray(ic(vv_mask), 1, [nlev, 1], @sum, 0);

% Count correct responses per condition (correct if trialOutcome equals visual).
hh_correct = accumarray(ic(hh_mask & (trialOutcome == 1)), 1, [nlev, 1], @sum, 0);
hv_correct = accumarray(ic(hv_mask & (trialOutcome == 1)), 1, [nlev, 1], @sum, 0);
vh_correct = accumarray(ic(vh_mask & (trialOutcome == 0)), 1, [nlev, 1], @sum, 0);
vv_correct = accumarray(ic(vv_mask & (trialOutcome == 0)), 1, [nlev, 1], @sum, 0);

% Build output matrix (columns: contrast, o, hh count, hv count, vh count, vv count,
% hh correct, hv correct, vh correct, vv correct)
output = zeros(nlev,10);
output(:,1) = c';
output(:,2) = o;
output(:,3) = hh_counts;
output(:,4) = hv_counts;
output(:,5) = vh_counts;
output(:,6) = vv_counts;
output(:,7) = hh_correct;
output(:,8) = hv_correct;
output(:,9) = vh_correct;
output(:,10)= vv_correct;

% ----- Compute Probabilities and Format yPredicted -----
% Compute probabilities pc and pic as in the original code.
if strcmp(optoStr, 'baseline')
    yPredicted = (output(:,7) + output(:,8) + output(:,9) + output(:,10)) ./ ...
              (output(:,3) + output(:,4) + output(:,5) + output(:,6) + eps);
elseif strcmp(optoStr, 'inconOpto')
    yPredicted = (output(:,8) + output(:,9)) ./ (output(:,4) + output(:,5) + eps);
elseif strcmp(optoStr, 'conOpto')
    yPredicted  = (output(:,7) + output(:,10)) ./ (output(:,3) + output(:,6) + eps);
end

%{
if strcmp(optoStr, 'baseline')
    % Control condition: average all correct responses.
    pcntrl = (output(:,7) + output(:,8) + output(:,9) + output(:,10)) ./ ...
              (output(:,3) + output(:,4) + output(:,5) + output(:,6) + eps);
    zeroIdx = find(c == 0);
    inconIdx = find(c < 0);
    conIdx = find(c > 0);
    if ~isempty(zeroIdx)
        zeroPart = nanmean([100 - (pcntrl(zeroIdx)*100), pcntrl(zeroIdx)*100], 2);
    else
        zeroPart = [];
    end
    yPredicted = [100 - (pcntrl(inconIdx)*100); zeroPart; pcntrl(conIdx)*100]';
else
    % Opto condition: use separate probabilities for incongruent responses.
    zeroIdx = find(c == 0);
    inconIdx = find(c < 0);
    conIdx = find(c > 0);
    if ~isempty(zeroIdx)
        zeroPart = nanmean([100 - (pic(zeroIdx)*100), pic(zeroIdx)*100], 2);
    else
        zeroPart = [];
    end
    yPredicted = [100 - (pic(inconIdx)*100); zeroPart; pic(conIdx)*100]';
end
%}
end
