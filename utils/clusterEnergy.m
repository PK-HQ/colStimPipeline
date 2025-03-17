function [bins, binEdges, clusterIdx, validIndices] = clusterEnergy(bitmapEnergies, bitmapColumns, method, nClusters, analysisParams)
    % Validate inputs
    if ~ismatrix(bitmapEnergies) || size(bitmapEnergies, 1) ~= 2
        error('Etotal must be a 2xn array.');
    end
    
    % Compute the mean power of each block
    meanPowers = mean(bitmapEnergies, 1);
    
    % If nClusters is not provided, determine it unsupervised
    if nargin < 3 || isempty(nClusters)
        maxClusters = 10; % Maximum number of clusters to consider
        if strcmpi(method, 'bin')
            nClusters = 3; % Default for binning if not specified
        else
            silhouetteScores = zeros(1, maxClusters);
            for k = 2:maxClusters % Start from 2 clusters to maxClusters
                if strcmpi(method, 'kmeans')
                    [clusterLabels, ~] = kmeans(meanPowers', k, 'MaxIter',10000, 'Replicates', 1000, 'Start', 'plus');
                elseif strcmpi(method, 'ward')
                    Z = linkage(meanPowers', 'ward');
                    clusterLabels = cluster(Z, 'maxclust', k);
                end
                % Compute silhouette scores
                silhouetteScores(k) = mean(silhouette(meanPowers', clusterLabels));
            end
            % Find the number of clusters with the maximum silhouette score
            [~, nClusters] = max(silhouetteScores(2:end));
            nClusters = nClusters + 1; % Correct for the index offset caused by starting at 2 clusters
        end
    end
    
    % Perform the clustering
    if strcmpi(method, 'kmeans')
        [clusterIdx, ~] = kmeans(meanPowers', nClusters);
    elseif strcmpi(method, 'ward')
        Z = linkage(meanPowers', 'ward');
        clusterIdx = cluster(Z, 'maxclust', nClusters);
    elseif strcmpi(method, 'bin')
        % Generate bin edges that double sequentially
        minPower = min(meanPowers);
        maxPower = max(meanPowers);
        binEdges = minPower;
        while binEdges(end) < maxPower
            binEdges = [binEdges, binEdges(end) * 2];
        end
        if binEdges(end) < maxPower
            binEdges = [binEdges, maxPower];
        end
        clusterIdx = discretize(meanPowers, binEdges);
    end
    
    % Organize blocks into bins
    uniqueClusters = unique(clusterIdx);
    m = length(uniqueClusters);
    bins = cell(1, m);
    clusterMedians = zeros(1, m);
    clusterMAD = zeros(1, m);
    clusterSilhouettes = zeros(1, m);
    for i = 1:m
        bins{i} = bitmapEnergies(:, clusterIdx == uniqueClusters(i));
        clusterMedians(i) = median(meanPowers(clusterIdx == uniqueClusters(i)));
        clusterMAD(i) = mad(meanPowers(clusterIdx == uniqueClusters(i)), 1);
        if ~strcmpi(method, 'bin') % Silhouette does not apply directly to bin method
            clusterSilhouettes(i) = mean(silhouette(meanPowers', clusterIdx, 'Euclidean'));
        end
    end

    %% Plot datapoints
    figure;
    axesHandle = gca; % Get the current axes handle
    % Set the colormap to 'inferno' before plotting
    colormap(axesHandle, slanCM('inferno')); % Apply the colormap to the axes
    hold on; 
    % Add bin boundaries for the 'bin' method
    if strcmpi(method, 'bin')
        for j = 1:length(binEdges)
            line([binEdges(j) binEdges(j)], [0 50], 'Color', 'k', 'LineStyle', ':', 'HandleVisibility', 'off');
        end
    end    
    % Plotting the clustered data blocks
    for i = 1:m
        scatter(meanPowers(clusterIdx == uniqueClusters(i)), find(clusterIdx == uniqueClusters(i)), 100, 'filled');
    end
    title('Clustered optostim blocks');
    xlabel('Block power (mW)');
    ylabel('Block number');
    if strcmpi(method, 'bin')
        legendEntries = arrayfun(@(x) sprintf('Cluster %d (%.1f \\leq Energy < %.1f)', uniqueClusters(x), binEdges(x), binEdges(x+1)), 1:m, 'UniformOutput', false);
        legend(legendEntries, 'Interpreter', 'tex', 'FontName', 'Arial');
    else
        legend(arrayfun(@(x) sprintf('Cluster %d (median=%.1f ± %.3f)', uniqueClusters(x), clusterMedians(x), clusterMAD(x)), 1:m, 'UniformOutput', false));
    end
    addSkippedTicks(0,50,5,'y')
    hold off;
    upFontSize

    %% Filter for n-columns
    % Define columns per bitmap
    columns=mean(bitmapColumns,1);
    
    %% Filter for n-columns within specified criteria
    if ~isempty(analysisParams)
        % Find indices where conditions are NOT met
        invalidIndices = ~(ismember(clusterIdx, analysisParams.cluster) & ...
                          (columns >= (analysisParams.columnMean - analysisParams.columnRange)) & ...
                          (columns <= (analysisParams.columnMean + analysisParams.columnRange)));
    
        % Set values at invalid indices to NaN
        clusterIdx(invalidIndices) = NaN;
        validIndices=find(~isnan(clusterIdx));
    else
        validIndices=[];
    end
end
