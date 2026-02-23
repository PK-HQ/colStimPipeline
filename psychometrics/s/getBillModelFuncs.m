function [columnarBayesModel, objectiveFunction] = getBillModelFuncs(common)

    % --- Simulation control knobs (tweak these) ---
    baseOptions = struct();
    baseOptions.useGPU = false;              % reproducibility during fitting
    baseOptions.seed = 1;                    % fixed seed => deterministic objective
    baseOptions.trialScale = 50;             % scales experimental trial counts (40/20) -> 2000/1000 sims
    baseOptions.minTrialsPerLevel = 200;     % floor so small bins never get too noisy
    baseOptions.balanceTrialTypes = true;    % ensures no empty HH/HV/VH/VV bins

    epsilon = 1e-10;
    sumAll = @(x) sum(x(:));

    % Convenience wrapper: simulate using per-contrast trial counts + optional budget override
    simWithCounts = @(x, params, sumY, mode, seedOffset, budget) ...
        localSim(x, params, sumY, mode, seedOffset, localMergeOptions(baseOptions, budget));


    % Expose model handle (not used directly in the objective anymore, but keep API)
    opts = struct();
    opts.useGPU = false;
    opts.seed = 1;
    opts.trialScale = 300;          % e.g., 40->8000, 20->4000
    opts.minTrialsPerLevel = 200;
    opts.balanceTrialTypes = true;
    opts.mode = 'both';
    columnarBayesModel = @(x, params) columnarBayesMdlGPT52(x, params, opts);

    objectiveFunction = @(params, data) localNLL(params, data, simWithCounts, epsilon, sumAll);

end

% ====== Local helpers ======

function mdl = localSim(x, params, sumY, mode, seedOffset, baseOptions)
    opts = baseOptions;

    % Use experimental counts (unbalanced) but scale them up for smoother estimates
    nPer = max(opts.minTrialsPerLevel, round(opts.trialScale * sumY(:)));
    opts.nTrialsPerLevel = nPer;

    opts.mode = mode;

    % Different calls (baseline vs stim) get different fixed seeds, but still deterministic
    opts.seed = opts.seed + seedOffset;

    mdl = columnarBayesMdlGPT52(x, params, opts);
end

function nLL = localNLL(params, data, simWithCounts, epsilon, sumAll)
    % Budget override from data (used for staged fitting)
    budget = struct();
    if isfield(data, 'simBudget') && ~isempty(data.simBudget)
        budget = data.simBudget;
    end

    % --- Baseline: control-only sim (faster) ---
    mdlBase = simWithCounts(data.xBaseline, params, data.sumBaseline, 'ctrlOnly', 101, budget);
    pBase = mdlBase.pcntrl(:) / 100;
    pBase(~isfinite(pBase)) = 0.5;
    pBase = min(max(pBase, epsilon), 1-epsilon);

    % --- Opto: if x grids are identical, simulate once and reuse pc + pic ---
    sameGrid = isequal(data.xConOpto(:), data.xInconOpto(:)) && isequal(data.sumConOpto(:), data.sumInconOpto(:));

    if sameGrid
        mdlStim = simWithCounts(data.xConOpto, params, data.sumConOpto, 'stimOnly', 202, budget);

        pCon = mdlStim.pc(:) / 100;
        pInc = mdlStim.pic(:) / 100;
    else
        mdlCon = simWithCounts(data.xConOpto, params, data.sumConOpto, 'stimOnly', 202, budget);
        mdlInc = simWithCounts(data.xInconOpto, params, data.sumInconOpto, 'stimOnly', 303, budget);

        pCon = mdlCon.pc(:) / 100;
        pInc = mdlInc.pic(:) / 100;
    end

    pCon(~isfinite(pCon)) = 0.5;
    pInc(~isfinite(pInc)) = 0.5;

    pCon = min(max(pCon, epsilon), 1-epsilon);
    pInc = min(max(pInc, epsilon), 1-epsilon);

    % Ensure counts are column vectors (prevents implicit expansion bugs)
    sB = data.successBaseline(:);   nB = data.sumBaseline(:);
    sI = data.successInconOpto(:);  nI = data.sumInconOpto(:);
    sC = data.successConOpto(:);    nC = data.sumConOpto(:);

    nLL = 0;

    nLL = nLL - sumAll( sB .* log(pBase) + (nB - sB) .* log(1 - pBase) );
    nLL = nLL - sumAll( sI .* log(pInc)  + (nI - sI) .* log(1 - pInc) );
    nLL = nLL - sumAll( sC .* log(pCon)  + (nC - sC) .* log(1 - pCon) );

    % Tiny bound barrier to discourage exact-bound solutions (very small effect)
    if isfield(data,'useBoundBarrier') && data.useBoundBarrier
        soft = 1e-4;
        p = params(:)';
    
        lb = [0, 0, 0, 0, 30, 1, 70, 0, 0];
        ub = [0.3, 0.3, 0.3, 0.3, 100, 6, 100, 1, 30];
    
        boundPenalty = sum( -log((p - lb) + soft) - log((ub - p) + soft) );
        nLL = nLL + 1e-3 * boundPenalty;  % keep tiny
    end

end
function out = localMergeOptions(base, override)
% Merge two structs; fields in override replace base.
    out = base;
    if nargin < 2 || isempty(override), return; end
    f = fieldnames(override);
    for i = 1:numel(f)
        out.(f{i}) = override.(f{i});
    end
end

%{
function [columnarBayesModel, objectiveFunction] = getBillModelFuncs(common)
    options = struct('nTrials', 100);
    columnarBayesModel = @(x, params) columnarBayesMdlGPT52(x, params, options);
        
    sumAll = @(x) sum(x(:));
    epsilon = 1e-10;
    objectiveFunction = @(params, data) ...
        - sumAll(data.successBaseline .* log(max(columnarBayesModel(data.xBaseline, params).pcntrl/100, epsilon)) + ...
                  (data.sumBaseline - data.successBaseline) .* log(max(1 - columnarBayesModel(data.xBaseline, params).pcntrl/100, epsilon))) ...
        - sumAll(data.successInconOpto .* log(max(columnarBayesModel(data.xInconOpto, params).pic/100, epsilon)) + ...
                  (data.sumInconOpto - data.successInconOpto) .* log(max(1 - columnarBayesModel(data.xInconOpto, params).pic/100, epsilon))) ...
        - sumAll(data.successConOpto .* log(max(columnarBayesModel(data.xConOpto, params).pc/100, epsilon)) + ...
                  (data.sumConOpto - data.successConOpto) .* log(max(1 - columnarBayesModel(data.xConOpto, params).pc/100, epsilon)));
end
%}