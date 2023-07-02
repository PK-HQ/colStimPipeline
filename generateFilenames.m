function filenameStruct=generateFilenames(dataStruct)
% def mainPath
if ispc
  mainPath='Y:/';
elseif contains(getenv('HOSTNAME'),'psy.utexas.edu')
  mainPath='/eslab/data/';
end

% set general calibration session (unused)
calibrationDate='20221104';%'20220904';
uniqueStr='S';
PCAsuffix='P2';%70%

%% Set dir, folder paths
dataPathOrt=[mainPath dataStruct.monkey '/' dataStruct.monkey dataStruct.date '/run0/']; %usually its run 0
dataPath=[mainPath dataStruct.monkey '/' dataStruct.monkey dataStruct.date '/run' dataStruct.run '/'];
imagePath=[mainPath dataStruct.monkey '/' dataStruct.monkey dataStruct.date '/'];
plotPath=imagePath;
calImagePath=[mainPath dataStruct.monkey '/' dataStruct.monkey calibrationDate 'Calibration/'];
folderfilePrefixStr=[dataPath 'M' dataStruct.monkeyNo 'D' dataStruct.date 'R' dataStruct.run];
filePrefixStr=['M' dataStruct.monkeyNo 'D' dataStruct.date 'R' dataStruct.run];

%% Define cross-session files
filenameStruct.plotPath=plotPath;
filenameStruct.calImg=[calImagePath 'SL00EX570DM600EM632L10_50um_aperture_reduced_binned.bmp'];
filenameStruct.calBmp=[calImagePath 'Calibrationseries/Cross0125Z00050L100.bmp'];

%% Define session image names
switch dataStruct.modality %isequal(dataStruct.modality,'GCaMP')
  
    case {'GCaMP','GEVI'}
        %% Trial structure file (TS)
        filenameStruct.TS=[folderfilePrefixStr 'TS.mat'];

        %% Functional responses (FFT/Intg)
        % Search for filenames filenames
        tmpRaw=dir([folderfilePrefixStr '*StabBin*.mat']);
        tmpFFT=dir([folderfilePrefixStr '*FFTAmp*PF0400.mat']);
        tmpIntg=dir([folderfilePrefixStr '*Intg*.mat']);
        % if raw binned exists         
        if ~isempty(tmpRaw)
            filenameStruct.Raw=[folderfilePrefixStr(1:end-numel(filePrefixStr)) tmpRaw.name];
        else
            filenameStruct.Raw='';
        end
        % if fft exists         
        if ~isempty(tmpFFT)
            filenameStruct.FFTAmp=[folderfilePrefixStr(1:end-numel(filePrefixStr)) tmpFFT.name];
        else
            filenameStruct.FFTAmp='';
        end
        % if integration exists
        if ~isempty(tmpIntg)
            filenameStruct.Intg=[folderfilePrefixStr(1:end-numel(filePrefixStr)) tmpIntg.name];
        else
            filenameStruct.Intg='';
        end

        %% Orientation map
        tmpOrientationP2=dir([folderfilePrefixStr '*OrientationP2*.mat']); %find with Orientation regex
        tmpOrientation=dir([folderfilePrefixStr '*Orientation.mat']); %find with Orientation regex

        %[folderfilePrefixStr 'Orientation' PCAsuffix '.mat'];
        if ~isempty(tmpOrientationP2) %to get orientation map for analysis of vis+opto block 
            filenameStruct.Orientation=[folderfilePrefixStr(1:end-numel(filePrefixStr)) tmpOrientationP2.name];
        elseif ~isempty(tmpOrientation) && isempty(tmpOrientationP2)
            filenameStruct.Orientation=[folderfilePrefixStr(1:end-numel(filePrefixStr)) tmpOrientation.name];%filenameStruct.Orientation=[dataPathOrt 'M' dataStruct.monkeyNo 'D' dataStruct.date 'R0Orientation.mat'];
        else 
            filenameStruct.Orientation='';
        end

        %% Static images
        % Unbinned for quality check
        tmpGreen=convertCharsToStrings(dir([imagePath '*EX540L_binned.tif']));
        filenameStruct.greenStatic=[imagePath tmpGreen.name];
        tmpGFP=convertCharsToStrings(dir([imagePath '*EX480*OD10_binned.tif']));
        filenameStruct.gfpStatic=[imagePath tmpGFP.name];
        tmpMCherry=convertCharsToStrings(dir([imagePath '*EX570*OD10_binned.tif']));
        filenameStruct.mcherryStatic=[imagePath tmpMCherry.name];

        % Binned green image (for coregistration)
        filenameStruct.greenImg=convertCharsToStrings([imagePath 'green_binned.bmp']);
        if ~isfile(filenameStruct.greenImg) && isfile([imagePath 'SL30EX540L_binned.tif'])
            greenImage=imread((convertCharsToStrings([imagePath 'SL30EX540L_binned.tif'])));
            imwrite(greenImage,filenameStruct.greenImg);
        end
        
       %% ROI timecourse
        tmpROITC=dir([folderfilePrefixStr '*ROITC.mat']); %find with ROITC regex
        filenameStruct.ROITC=[folderfilePrefixStr(1:end-numel(filePrefixStr)) tmpROITC.name];
        filenameStruct.vdaq=[dataPath 'Data_vdaqlog.mat']; % contains allegroimaging info
        
        
       %% Coregistration landmarks (manual) and transformation parameters (solution for coregistration)
        filenameStruct.landmarks=[imagePath 'landmarks.mat'];
        filenameStruct.transformParams=[imagePath 'transformParams.mat'];

        if isfield(dataStruct,'bmpOriginal') && isfield(dataStruct,'responseOriginal')
            filenameStruct.bmpOriginal=dataStruct.bmpOriginal;
            filenameStruct.responseOriginal=dataStruct.responseOriginal;
        end

       %% Save file names
        filenameStruct.colmapPDF=[filenameStruct.plotPath 'M' dataStruct.monkeyNo uniqueStr dataStruct(1).date PCAsuffix '.pdf'];
        filenameStruct.neurometricPDF=[mainPath 'Chip/Meta/summary/' 'M' dataStruct.monkeyNo uniqueStr dataStruct(1).date  'R' dataStruct.run 'summary.pdf'];%[filenameStruct.plotPath 'M' dataStruct.monkeyNo uniqueStr dataStruct(1).date 'NeurometricSummary.pdf'];
        filenameStruct.colmapMat=[filenameStruct.plotPath 'M' dataStruct.monkeyNo uniqueStr dataStruct(1).date PCAsuffix '.mat'];




    case {'VSD'}
        tmpFFT=dir([folderfilePrefixStr '*FFTAmp*0.mat']);
        filenameStruct.FFTAmp=[folderfilePrefixStr(1:end-numel(filePrefixStr)) tmpFFT.name];
        filenameStruct.FFTAmp
        %filenameStruct.FFTAmp=[basicStr 'DFFTAmpS016E115PF0500.mat'];
        filenameStruct.Orientation=[folderfilePrefixStr 'Orientation.mat'];
        filenameStruct.TS=[folderfilePrefixStr 'TS.mat'];
        filenameStruct.greenImg=convertCharsToStrings([imagePath 'SL30EX540L_binned_gainX5.tif']);
        filenameStruct.redImg=convertCharsToStrings([imagePath 'SL00Red_binned.bmp']);
        filenameStruct.initTransformParams=convertCharsToStrings([imagePath 'initTransformParams.mat']);
end
end