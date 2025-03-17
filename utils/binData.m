function [xAverage, yAverage, counts] = binData(x, y)
    % Define the bins
    edges = 0:5:100;
    binCenters = edges(1:end-1) + 2.5;
    
    % Initialize the size of xAverage
    xAverage = repmat(binCenters, size(x, 1), 1);
    yAverage = NaN(size(x, 1), length(binCenters));
    counts = zeros(size(x, 1), length(binCenters));
    
    % Iterate over each condition and each block
    for condition = 1:size(x, 1)
        % Aggregate all blocks for this condition
        conditionX = squeeze(x(condition, :, :));  % Reducing one dimension
        conditionY = squeeze(y(condition, :, :));  % Reducing one dimension
        
        % Initialize temporary storage for y values per bin across blocks
        tempY = cell(1, length(binCenters));
        for i = 1:length(binCenters)
            tempY{i} = [];
        end
        
        for block = 1:size(x, 3)
            % Extract this block's data
            xBlock = conditionX(:, block);
            yBlock = conditionY(:, block);
            
            % Bin the data based on x values
            [~, ~, binIdx] = histcounts(xBlock, edges);
            
            % Aggregate y values in each bin
            for j = 1:length(binCenters)
                validIdx = binIdx == j;
                tempY{j} = [tempY{j}; yBlock(validIdx)];  % Append to the list
            end
        end
        
        % Average y values in each bin and count the entries
        for j = 1:length(binCenters)
            if ~isempty(tempY{j})
                yAverage(condition, j) = nanmean(tempY{j});
                counts(condition, j) = sum(~isnan(tempY{j}));
            end
        end
    end
    
    [yAverage,idx]=rmnan(yAverage);
    
    for row=1:3
        cont(row,:)=xAverage(row,idx(row,:),:);
    end
    xAverage=cont;
end
