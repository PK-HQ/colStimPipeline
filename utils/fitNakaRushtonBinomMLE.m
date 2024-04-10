function [fitParams, negLogLikelihood] = fitNakaRushtonBinomMLE(xBlocks, yBlocks)
    % Assuming x and y are 3D matrices: 3 opto-conditions x 6 contrasts x 16 blocks
    % Define bin edges according to the specified range and bin width
    binEdges = 0:5:100;
    
    % Pre-allocate arrays to store fit results and likelihoods
    fitParams = zeros(16, 3, size(xBlocks, 1)); % Blocks x Parameters x Conditions
    negLogLikelihood = zeros(16, size(xBlocks, 1)); % Blocks x Conditions
    
    % Loop through each block and condition
    for block = 1:size(xBlocks, 3)
        for cond = 1:size(xBlocks, 1)
            % Extract and bin data for the current block and condition
            %[binnedX, binnedY] = binData(squeeze(xBlocks(cond, :, block)), squeeze(yBlocks(cond, :, block)), binEdges);
            binnedX=squeeze(xBlocks(cond, :, block)); binnedX=binnedX(~isnan(binnedX));
            binnedY=squeeze(yBlocks(cond, :, block)); binnedY=binnedY(~isnan(binnedY));

            % Initial parameter guess
            [initRmax, initExp, initC50, initBeta] = guessInitialParams(binnedX, binnedY);
            
            % Modify the Naka-Rushton function as needed to use the parameters correctly
            nakarushton = @(params, xi) ((1-params(1)) .* xi.^params(2)) ./ (params(3).^params(2) + xi.^params(2)) + params(1);
            
            % Define the objective function for the MLE fit
            objectiveFunc = @(params) -sum(binopdf(round(binnedY'*100), 100, nakarushton(params, binnedX)));
            
            % Optimization settings
            options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off');
            lb = [0, 0, 0]; % Lower bounds: Beta, n, C50
            ub = [1, 10, 100]; % Upper bounds: Beta, n, C50
            
            % Use initial guesses based on condition
            if cond == 1
                initialGuess = [initBeta, initExp, initC50]; % Condition-specific adjustments
            else
                initialGuess = [initBeta, initExp, initC50]; % General case
            end
            
            % Perform the optimization
            [params, fval] = fmincon(objectiveFunc, initialGuess, [], [], [], [], lb, ub, [], options);
            
            % Store the fitted parameters and negative log likelihood
            fitParams(block, :, cond) = params;
            negLogLikelihood(block, cond) = -fval;
        end
    end
end

%% === Subfunctions ===
% Bin data
function [binnedX, binnedY] = binData(x, y, binEdges)
    % Ensure x and y are column vectors
    x = x(:);
    y = y(:);
    
    % Perform binning of x data
    [counts, edges, binIdx] = histcounts(x, binEdges);
    binnedX = edges(1:end-1) + diff(edges) / 2; % Calculate bin centers

    % Initialize binnedY with NaNs to handle bins without data
    binnedY = NaN(length(binnedX), 1);

    % Filter out zero indices from binIdx
    validIndices = binIdx > 0; % Exclude indices that do not correspond to a bin
    binIdxFiltered = binIdx(validIndices);
    yFiltered = y(validIndices);

    % Use accumarray with filtered indices and values
    accumY = accumarray(binIdxFiltered, yFiltered, [], @nanmean);
    
    % Assign calculated means to valid bins
    binnedY(1:length(accumY)) = accumY;
    
    % Remove bins with no data if necessary. However, keeping NaNs might be useful for plotting or analysis to retain bin structure.
    % Uncomment the next two lines to remove NaNs, though this might not be necessary or desired.
    % validBins = ~isnan(binnedY);
    % binnedX = binnedX(validBins); binnedY = binnedY(validBins);
end

% Initial guess
function [initRmax, initExp, initC50, initBeta] = guessInitialParams(x, yMeans)
    % guess NR: rmax (delta between mean responses)
    initRmax = max(yMeans) - min(yMeans);
    
    % guess NR: exp (intermediate between min and max exponent)
    initExp = mean([3]);
    
    % guess NR: c50 (midpoint between the delta of mean responses)
    initC50 = computeInitC50(yMeans, initRmax, x);
    
    % guess vertical offset: beta (lowest mean response)
    initBeta = min(yMeans);
end

function initC50 = computeInitC50(yMeans, initRmax, x)
    absoluteDelta = abs(yMeans - (initRmax / 2 + min(yMeans)));
    closest_index = find(absoluteDelta == min(absoluteDelta));

    % Adjust for cases with multiple equally close values
    if numel(closest_index) > 2
        closest_index = closest_index([1, 2]); % Take first two if more than two indices are equally close
    elseif numel(closest_index) == 2
        % Already ideal scenario, take them as they are
    else
        % For a single closest index, attempt to add the next index if it doesn't exceed bounds
        if closest_index < numel(yMeans)
            closest_index = [closest_index, closest_index + 1];
        % If the closest index is the last element, consider the previous one as well, if possible
        elseif closest_index > 1
            closest_index = [closest_index - 1, closest_index];
        end
    end

    % Compute the mean of the x values at the adjusted closest indices
    initC50 = mean(x(closest_index));
end