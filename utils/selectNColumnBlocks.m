function [desiredBlocks,blockColumns,blockEnergy] = selectNColumnBlocks(bitmapData, varargin)
    % selectNColumnBlocks selects blocks based on column count and energy threshold,
    % with optional filtering by specified cluster blocks.
    %
    % Inputs:
    %   bitmapData      - Struct containing 'nColumns' and 'energy' data.
    %   varargin        - Optional inputs:
    %                     1. clusterBlocks: Indices of cluster blocks to filter (default: all blocks).
    %                     2. columnsDesired: Desired number of columns for filtering (required).
    %                     3. columnSpread: Allowable spread around the desired number of columns (required).
    %                     4. energyThreshold: Maximum energy threshold for filtering (required).
    %
    % Outputs:
    %   desiredBlocks   - Indices of blocks meeting the criteria.

    % Parse inputs
    if numel(varargin) < 3
        error('You must specify at least columnsDesired, columnSpread, and energyThreshold.');
    end

    % Extract required and optional inputs
    clusterBlocks = [];
    if numel(varargin) == 3
        clusterBlocks = varargin{1};
        columnsDesired = varargin{2};
        columnSpread = varargin{3};
        energyThreshold = 999; %max
    elseif numel(varargin) == 4
        clusterBlocks = varargin{1};
        columnsDesired = varargin{2};
        columnSpread = varargin{3};
        energyThreshold = varargin{4};
    end

    % Compute mean columns and energy
    columns=bitmapData.nColumns';
    columnsMean = mean(columns, 2);
    bitmapEnergy = squeeze(bitmapData.energy(:,:,:))';
    bitmapEnergyMean = mean(bitmapEnergy, 2);

    % Filter blocks based on column and energy criteria
    blockIds = find(columnsMean >= columnsDesired - columnSpread & ...
                    columnsMean <= columnsDesired + columnSpread & ...
                    bitmapEnergyMean < energyThreshold);

    % If clusterBlocks is specified, filter blockIds further
    if ~isempty(clusterBlocks)
        desiredBlocks = intersect(clusterBlocks, blockIds);
        blockColumns=columns(desiredBlocks,:);
        blockEnergy=bitmapEnergy(desiredBlocks,:);
    else
        desiredBlocks = blockIds;
        bitmapEnergyMean=bitmapEnergyMean(desiredBlocks);
    end
end
