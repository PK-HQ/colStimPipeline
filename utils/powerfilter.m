function [selectedBlocks, selectedBlockIDs]=powerfilter(varargin)
if length(varargin)==4
    bitmapData=varargin{1};
    analysisBlockID=varargin{2};
    powerFlag=varargin{3};
    nColumnRange=varargin{4};
elseif length(varargin)==3
    bitmapData=varargin{1};
    analysisBlockID=varargin{2};
    powerFlag=varargin{3};
    nColumnRange=[1 40];
end

switch powerFlag
    case {'fixed'}
        % Extract the mean adjusted SPD across trials
        adjustedSPD_uW = squeeze(mean(bitmapData.adjustedSPD_uW, 2));
    
        % Extract by blockID, 5-40 nColumns
        nonPowerSeries=find(analysisBlockID<=105 & ... % avoid power series
            mean(bitmapData.nColumns,1)>=nColumnRange(1) & mean(bitmapData.nColumns,1)<=nColumnRange(2)); % within range of columns wanted
        
        % Calculate IQR and identify outliers
        Q1 = quantile(adjustedSPD_uW, 0.25);
        Q3 = quantile(adjustedSPD_uW, 0.75);
        IQR = Q3 - Q1;
        lowerBound = Q1 - 1.5 * IQR;
        upperBound = Q3 + 1.5 * IQR;
        selectedBlocks = intersect(nonPowerSeries, find(adjustedSPD_uW >= lowerBound & adjustedSPD_uW <= upperBound))';
end
selectedBlockIDs = analysisBlockID(selectedBlocks);
end