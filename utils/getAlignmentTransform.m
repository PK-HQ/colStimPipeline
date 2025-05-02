function alignmentTransform = getAlignmentTransform(dataStruct)
    % Let user pick the known (original) BMP file
    cd('V:\PK\ColSeries\')
    [bmpFile, bmpPath] = uigetfile({'*.bmp', 'BMP Files (*.bmp)'}, ...
        'Select the known/original BMP file');
    if bmpFile == 0
        error('No BMP file selected');
    end
    bmpOriginal = double(imread(fullfile(bmpPath, bmpFile)));
    
    % Let user pick the empirical (measured) response file
    cd('E:\')
    [responseFile, responsePath] = uigetfile({'*.bmp', 'BMP Files (*.bmp)'}, ...
        'Select the empirical/measured response file');
    if responseFile == 0
        error('No response file selected');
    end
    responseOriginal = double(imread(fullfile(responsePath, responseFile)));
    
    % Process response image
    thresholdBmpEst = 98;
    responseOriginalHE = adapthisteq(rescale(responseOriginal,0,1), ...
        'NumTiles', size(responseOriginal)./16, 'Range', 'full');
    responseOriginalFilt = double(responseOriginalHE >= ...
        prctile(responseOriginalHE(:), thresholdBmpEst));
    figure();imagesc(responseOriginalFilt)
  
    % Setup basic structures needed
    blockID=1;
    currentBlockStruct =  generateFilenames(dataStruct);
    behavioralData = struct(); % Not critical for cam2proj
    behavioralData.optoTS(blockID).Header.Conditions.ProjTTLPulseOn=0;
    behavioralData.optoTS(blockID).Header.Conditions.ProjTTLPulseOff=50;
    
    % Required fields for bitmapData
    bitmapData = struct();
    bitmapData.columnarbitmapCoreg(:,:,1,1) = responseOriginalFilt; % Your processed response image
    bitmapData.orts = zeros(1,1,1);  % Orientation (0 or 90 degrees)
    bitmapData.gaussianContourLevel =nan;  % Need this from you (1-200)
    bitmapData.gaussianContourLevelMax = nan;
    bitmapData.nColumnsWanted = 1;  % Need this from you
    bitmapData.projectorx = 1920;  % Standard projector resolution
    bitmapData.projectory = 1080;
    bitmapData.pixelsizemm = 0.0054;  % Need this from you - projector pixel size in mm
    bitmapData.gaussianCond(1)=nan;
    
    % Required fields for imagingData
    imagingData = struct();
    imagingData.pixels = [512 512]; % Camera resolution
    imagingData.pixelsizemm = 0.0161; % Example value, adjust as needed
    imagingData.mask(:,:, blockID)=ones(imagingData.pixels(1), imagingData.pixels(1));
    pdfFilename='test.pdf';
    % Generate bitmap, correct for projector properties and camera-projector alignment
    [bitmapData1] = convertForProjectorGPT(behavioralData, imagingData, bitmapData,...
        currentBlockStruct, 'cam2proj', blockID,...
        pdfFilename, 1, 0, 0);

    % Normalize estimated BMP
    bmpEstimated = rescale(double(bitmapData1.columnarbitmapTFprojspace));

    % set levels
    regParams.PyramidLevels=3;

    scale=.67;
    transX=285;
    transY=195;
    shear=0;
    rotate=0;
    % Create initial transformation with your parameters
    initTransform = affine2d([
        scale   shear         0;    % Scale & rotation
        rotate    scale         0;    % Scale & rotation
        transX  transY    1.0000   % Translation
    ]);

    % Then use it in your registration
    [optimizer,metric] = imregconfig('multimodal');
    optimizer.InitialRadius = 6.25*10^-6;
    optimizer.GrowthFactor = 1.05;
    transformType = 'affine';

    % Your registration call with the initial transform
    alignmentTransform = imregtform(bmpEstimated, bmpOriginal, transformType, ...
        optimizer, metric, 'DisplayOptimization', false, ...
        'PyramidLevels', regParams.PyramidLevels, ...
        'InitialTransformation', initTransform);

    % Display results
    displayResults(bmpOriginal, bmpEstimated, alignmentTransform);
    
    save([responsePath 'alignmentTransform.mat'],'alignmentTransform')
end

function displayResults(bmpOriginal, bmpEstimated, transform)
    % Displays original, estimated, and transformed images
    figure;
    [hA,~] = tight_subplot(2,3);
    
    axes(hA(1)); imagesc(bmpOriginal); title('bmpOriginal');
    axes(hA(2)); imagesc(bmpEstimated); title('bmpEstimated');
    
    % Pre-registration comparison
    axes(hA(3));
    imagesc(imfuse(bmpOriginal, bmpEstimated, 'ColorChannels', 'red-cyan'));
    title('Pre-coregistration');

    % Show transformed result
    bmpEstimatedTransformed = imwarp(bmpEstimated, transform, ...
        'OutputView', imref2d(size(bmpEstimated)));
    axes(hA(5));
    imagesc(imfuse(bmpOriginal, bmpEstimatedTransformed, ...
        'ColorChannels', 'red-cyan'));
    title('Final Result');
end

function val = getOr(struct, field, default)
    % Utility to get field from struct or return default
    if isfield(struct, field)
        val = struct.(field);
    else
        val = default;
    end
end

function bmpPath = constructBmpPath(folder, ort, gamma)
    % Constructs BMP file path with standardized format
    bmpPath = fullfile('X:', 'PK', 'ColSeries', folder, ...
        sprintf('O%05.fHE1000%sT09900.bmp', ort, gamma));
end

function [responseOriginal, responseOriginalFilt] = loadAndProcessResponse(...
    currentSession, referenceBMPfolder, threshold)
    % Loads and processes the response image
    responsePath = fullfile('F:', 'OIdata', currentSession, ...
        ['calib100_' referenceBMPfolder '.bmp']);
    responseOriginal = double(imread(responsePath));
    
    % Process response image
    responseOriginalHE = adapthisteq(rescale(responseOriginal,0,1), ...
        'NumTiles', size(responseOriginal)./16, 'Range', 'full');
    responseOriginalFilt = double(responseOriginalHE >= ...
        prctile(responseOriginalHE(:), threshold));
end

function dataStruct = createDataStruct(date, bmpPath, responsePath)
    % Creates standardized data structure for convertForProjector
    dataStruct.monkeyNo = '28';
    dataStruct.monkey = 'Chip';
    dataStruct.date = date;
    dataStruct.run = '0';
    dataStruct.site = '30';
    dataStruct.baselineTS = '0';
    dataStruct.modality = 'GCaMP';
    dataStruct.bmpOriginal = bmpPath;
    dataStruct.responseOriginal = responsePath;
    dataStruct.gaussianContourLevel = [1 1];
    dataStruct.gaussianContourLevelMax = 200;
end
