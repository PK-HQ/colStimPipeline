figure;
% Init figure
make_it_tight = true;
hmarg=.15;
subplot = @(m,n,p) subtightplot (m, n, p, [0.1 0.1], [hmarg hmarg], [0.1 0.1]);
if ~make_it_tight,  clear subplot;  end

figure('Name','Statics');
%Masks
load('Y:\Chip\Chip20240112\run0\M28D20240112R0OrientationP2.mat')
nanMaskChip=double(Mask);nanMaskChip(nanMaskChip==0)=nan;
load('Y:\Pepper\Pepper20250114\run5\M32D20250114R5Orientation.mat')
nanMaskPepper=double(Mask);nanMaskPepper(nanMaskPepper==0)=nan;

%GCaMPs
subplot(2,3,1)
imageData1=double(imread('Y:\Chip\Chip20240112\SL1EX480DM505EM520L100OD16_binned.tif')).*(nanMaskChip)./255;
imgsc(imageData1);axis square;colormap(fire);colorbar; 
title('Chip-GCAMP')
subplot(2,3,2)
imageData2=double(imread('Y:\Pepper\Pepper20250114\SL2EX480DM505EM520L100OD016_binned.tif')).*(nanMaskPepper)./255;
imgsc(imageData2);axis square;colormap(fire);colorbar;
title('Pepper-GCAMP')
subplot(2,3,3)
imageData3 = calcPercentWithOutlierHandling(imageData1, imageData2).*nanMaskPepper; disp(mean(imageData3(:),'omitnan'))
imgsc(imageData3);axis square;colormap(fire);colorbar;
colormap(gca, fireice); % Set colormap for this subplot only
title('% Pepper / Chip - GCaMP')
title('Pepper / Chip = 58% (GCaMP should increase by 74%)')

%ChRmines
subplot(2,3,4)
imageData1=double(imread('Y:\Chip\Chip20240112\SL1EX570DM600EM632L13OD0_binned.tif')).*(nanMaskChip)./255;
imgsc(imageData1);axis square;colormap(fire);colorbar;
title('Chip-ChRmine')
subplot(2,3,5)
imageData2 =double(imread('Y:\Pepper\Pepper20250114\SL2EX570DM600EM632L13OD0_binned.tif')).*(nanMaskPepper)./255;
imgsc(imageData2);axis square;colormap(fire);colorbar;
title('Pepper-ChRmine')
subplot(2,3,6)
imageData3 = calcPercentWithOutlierHandling(imageData1, imageData2).*nanMaskPepper; disp(mean(imageData3(:),'omitnan'))
imgsc(imageData3);axis square;colormap(fire);colorbar;
colormap(gca, fire); % Set colormap for this subplot only
title('Pepper / Chip = 135% (ChRmine should reduce by 25%)')

function percentDifference = calcPercent(img1, img2)
    % Calculate the percentage of img2 relative to img1 per pixel.
    % Inputs:
    %   img1 - First image (matrix of size MxN)
    %   img2 - Second image (matrix of size MxN)
    % Output:
    %   percentDifference - Percentage of img2 relative to img1 per pixel.

    % Ensure the inputs are of the same size
    if ~isequal(size(img1), size(img2))
        error('Images must be of the same size.');
    end

    % Convert images to double for accurate calculations
    img1 = double(img1);
    img2 = double(img2);

    % Initialize percentDifference
    percentDifference = zeros(size(img1));

    % Calculate percentage (img2 ./ img1) * 100 for valid pixels
    mask = img1 > 0; % Pixels where division is valid
    percentDifference(mask) = (img2(mask) ./ img1(mask)) * 100;

    % Handle cases where img1 = 0
    maskZero = img1 == 0;
    percentDifference(maskZero & img2 > 0) = Inf; % Division by zero results in Inf
    percentDifference(maskZero & img2 == 0) = 0;  % Both are zero, set to 0

    % Display the results
    fprintf('Percentage of img2 relative to img1 calculated for each pixel.\n');
end
function percentDifference = calcPercentWithOutlierHandling(img1, img2)
    % Calculate the percentage of img2 relative to img1 per pixel and handle outliers.
    % Inputs:
    %   img1 - First image (matrix of size MxN)
    %   img2 - Second image (matrix of size MxN)
    % Output:
    %   percentDifference - Percentage of img2 relative to img1 per pixel with outliers set to NaN.

    % Ensure the inputs are of the same size
    if ~isequal(size(img1), size(img2))
        error('Images must be of the same size.');
    end

    % Convert images to double for accurate calculations
    img1 = double(img1);
    img2 = double(img2);

    % Initialize percentDifference
    percentDifference = zeros(size(img1));

    % Calculate percentage (img2 ./ img1) * 100 for valid pixels
    mask = img1 > 0; % Pixels where division is valid
    percentDifference(mask) = (img2(mask) ./ img1(mask)) * 100;

    % Handle cases where img1 = 0
    maskZero = img1 == 0;
    percentDifference(maskZero & img2 > 0) = Inf; % Division by zero results in Inf
    percentDifference(maskZero & img2 == 0) = 0;  % Both are zero, set to 0

    % Handle outliers: Replace outliers with NaN
    % Define outliers as values > mean ± 3*std (can be adjusted)
    validValues = percentDifference(isfinite(percentDifference)); % Exclude Inf and NaN
    meanVal = mean(validValues, 'omitnan');
    stdVal = std(validValues, 'omitnan');
    lowerBound = meanVal - 3 * stdVal;
    upperBound = meanVal + 3 * stdVal;

    % Set outliers to NaN
    percentDifference(percentDifference < lowerBound | percentDifference > upperBound) = NaN;

    % Display the results
end
