function [mse, ssimValues, ncc] = imageSimilarity(img1, img2)
    % Ensure both images are double for calculations
    img1 = double(img1);
    img2 = double(img2);
    
    % Check that the images have the same first two dimensions
    if ~isequal(size(img1, 1), size(img2, 1)) || ~isequal(size(img1, 2), size(img2, 2))
        error('The first two dimensions of both images must match for similarity comparison.');
    end
    
    % Get the number of slices in img1
    numSlices = size(img1, 3);
    
    % Preallocate arrays for the outputs
    mse = zeros(numSlices, 1);
    ssimValues = zeros(numSlices, 1);
    ncc = zeros(numSlices, 1);
    
    % Loop through each slice
    for i = 1:numSlices
        % Extract the i-th slice from img1 and compare with img2
        img1_slice = img1(:,:,i);
        
        % Create a mask for non-NaN areas in both images
        validMask = ~isnan(img1_slice) & ~isnan(img2);
        
        % Apply the mask to select only non-NaN values
        img1_valid = img1_slice(validMask);
        img2_valid = img2(validMask);
        
        % Calculate Mean Squared Error (MSE) for non-NaN areas
        mse(i) = mean((img1_valid - img2_valid).^2);
        
        % Calculate Structural Similarity Index (SSIM) for non-NaN areas
        ssimValues(i) = ssim(img1_valid, img2_valid);
        
        % Calculate Normalized Cross-Correlation (NCC) for non-NaN areas
        img1Mean = mean(img1_valid);
        img2Mean = mean(img2_valid);
        numerator = sum((img1_valid - img1Mean) .* (img2_valid - img2Mean));
        denominator = sqrt(sum((img1_valid - img1Mean).^2) * sum((img2_valid - img2Mean).^2));
        ncc(i) = numerator / denominator;
    end
    
    % Display the results for each slice
    for i = 1:numSlices
        fprintf('Slice %d - Mean Squared Error (MSE): %f\n', i, mse(i));
        fprintf('Slice %d - Structural Similarity Index (SSIM): %f\n', i, ssimValues(i));
        fprintf('Slice %d - Normalized Cross-Correlation (NCC): %f\n', i, ncc(i));
    end
end
