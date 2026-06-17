function [baselineMode, isSeparate, baselineSourceValue] = ...
        detectBaselineModeByBlock(datastruct, analysisBlockID, blockIndices)
% Classify each experiment as combined-baseline or separate-baseline.
%
% The authoritative source is datastruct(...).baselineTS, matching
% analyzeBlockPsychometrics/generateFilenames. Empty or missing baselineTS
% means baseline trials are combined with the opto block.

    if nargin < 3 || isempty(blockIndices)
        blockIndices = 1:numel(analysisBlockID);
    end

    blockIndices = blockIndices(:);
    baselineMode = strings(numel(blockIndices), 1);
    isSeparate = false(numel(blockIndices), 1);
    baselineSourceValue = strings(numel(blockIndices), 1);

    for idx = 1:numel(blockIndices)
        blockIdx = blockIndices(idx);
        baselineMode(idx) = "combined";
        baselineSourceValue(idx) = "";

        if ~isfinite(blockIdx) || blockIdx ~= round(blockIdx) || ...
                blockIdx < 1 || blockIdx > numel(analysisBlockID)
            baselineMode(idx) = "unknown";
            continue;
        end

        datastructIdx = analysisBlockID(blockIdx);
        if ~isfinite(datastructIdx) || datastructIdx ~= round(datastructIdx) || ...
                datastructIdx < 1 || datastructIdx > numel(datastruct)
            baselineMode(idx) = "unknown";
            continue;
        end

        if ~isfield(datastruct, 'baselineTS')
            baselineMode(idx) = "combined";
            continue;
        end

        baselineValue = datastruct(datastructIdx).baselineTS;
        baselineSourceValue(idx) = stringifyBaselineValue(baselineValue);
        isSeparate(idx) = ~isBaselineFieldEmpty(baselineValue);
        if isSeparate(idx)
            baselineMode(idx) = "separate";
        else
            baselineMode(idx) = "combined";
        end
    end
end

function tf = isBaselineFieldEmpty(value)
    if isempty(value)
        tf = true;
    elseif isstring(value) || ischar(value)
        tf = all(strlength(string(value(:))) == 0);
    elseif iscell(value)
        tf = isempty(value) || all(cellfun(@isBaselineFieldEmpty, value));
    else
        tf = false;
    end
end

function textValue = stringifyBaselineValue(value)
    if isempty(value)
        textValue = "";
    elseif iscell(value)
        textValue = strjoin(string(value(:)), ";");
    else
        textValue = strjoin(string(value(:)), ";");
    end
end
