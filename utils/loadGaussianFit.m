function ROIMaskgaussian = loadGaussianFit(filenameStruct)
    % Load Gaussian fit data, return empty if file does not exist
    if isfile(filenameStruct.gaussianFit)
        load(filenameStruct.gaussianFit); % .gaussian file
        ROIMaskgaussian = gaussianFit; % Assuming gaussianFit is the variable name in the file
    else
        ROIMaskgaussian = [];
    end
end