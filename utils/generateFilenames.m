function filenameStruct=generateFilenames(dataStruct)
%% --- Names ---

% pdfFilename=filenameStructSession.neurometricPDF;
% mapsFilename=filenameStructSession.colmapMat;

% def filenameStruct.mainPath
if ispc
  filenameStruct.mainPath='Y:/';
elseif contains(getenv('HOSTNAME'),'psy.utexas.edu')
  filenameStruct.mainPath='/eslab/data/';
end

% set general calibration session (unused)
calibrationDate='20221104';%'20220904';
uniqueStr='S';
PCAsuffix='P2';%70%

%% Set dir, folder paths
filenameStruct.monkey=dataStruct.monkey;
filenameStruct.monkeyNo=dataStruct.monkeyNo;
if isfield(dataStruct,'chamber')
    filenameStruct.chamber=dataStruct.chamber;
end
filenameStruct.date=dataStruct.date;
filenameStruct.run=dataStruct.run;
filenameStruct.modality=dataStruct.modality;

dataPathOrt=[filenameStruct.mainPath dataStruct.monkey '/' dataStruct.monkey dataStruct.date '/run' dataStruct.run '/']; %usually its run 0
dataPath=[filenameStruct.mainPath dataStruct.monkey '/' dataStruct.monkey dataStruct.date '/run' dataStruct.run '/'];
baselinePath=[filenameStruct.mainPath dataStruct.monkey '/' dataStruct.monkey dataStruct.date '/run' dataStruct.baselineTS '/'];

imagePath=[filenameStruct.mainPath dataStruct.monkey '/' dataStruct.monkey dataStruct.date '/'];
plotPath=imagePath;
calImagePath=[filenameStruct.mainPath dataStruct.monkey '/' dataStruct.monkey calibrationDate 'Calibration/'];
dataPrefixStr=[dataPath 'M' dataStruct.monkeyNo 'D' dataStruct.date 'R' dataStruct.run];
baselinePrefixStr=[baselinePath 'M' dataStruct.monkeyNo 'D' dataStruct.date 'R' dataStruct.baselineTS];
filePrefixStr=['M' dataStruct.monkeyNo 'D' dataStruct.date 'R' dataStruct.run];

%% Define cross-session files
filenameStruct.plotPath=plotPath;
filenameStruct.calImg=[calImagePath 'SL00EX570DM600EM632L10_50um_aperture_reduced_binned.bmp'];
filenameStruct.calBmp=[calImagePath 'Calibrationseries/Cross0125Z00050L100.bmp'];

%% Define session image names
switch dataStruct.modality %isequal(dataStruct.modality,'GCaMP')
  
    case {'GCaMP','GEVI'}
        %% Trial structure file (TS)
        filenameStruct.TS=[dataPrefixStr 'TS.mat'];
         
        %% Load baseline imaging
        filenameStruct.baselineTS=[baselinePrefixStr 'TS.mat'];
        baselineIntg=dir([baselinePrefixStr '*Intg*.mat']);
        if ~isempty(dataStruct.baselineTS)
            % if integration exists
            if ~isempty(baselineIntg)
                filenameStruct.baselineIntg=fullfile(baselineIntg(1).folder, baselineIntg(1).name);%[baselinePrefixStr(1:end-numel(filePrefixStr)) baselineIntg.name];
            else
                filenameStruct.baselineIntg='';
            end
        else
            filenameStruct.baselineTS=[];
            filenameStruct.baselineIntg='';
        end
        
        %% Functional responses (FFT/Intg)
        % Search for filenames filenames
        tmpRaw=dir([dataPrefixStr '*StabBin*.mat']);
        tmpFFT=dir([dataPrefixStr '*FFTAmp*PF0400.mat']);
        tmpIntg=dir([dataPrefixStr '*Intg*.mat']);
        tmpStab=dir([dataPrefixStr '*StabCfg.mat']);

        % if raw binned exists         
        if ~isempty(tmpRaw)
            filenameStruct.Stab=[dataPrefixStr(1:end-numel(filePrefixStr)) tmpStab.name];
        else
            filenameStruct.Stab='';
        end
        % if raw binned exists         
        if ~isempty(tmpRaw)
            filenameStruct.Raw=[dataPrefixStr(1:end-numel(filePrefixStr)) tmpRaw.name];
        else
            filenameStruct.Raw='';
        end
        % if fft exists         
        if ~isempty(tmpFFT)
            filenameStruct.FFTAmp=fullfile(tmpFFT(1).folder, tmpFFT(1).name);%[dataPrefixStr(1:end-numel(filePrefixStr)) tmpFFT.name];
        else
            filenameStruct.FFTAmp='';
        end
        %% Load opto/opto+baseline imaging
        if ~isempty(tmpIntg)
            %prefix=dataPrefixStr(1:end-numel(filePrefixStr));
            %suffix=tmpIntg.name;
            filenameStruct.Intg=fullfile(tmpIntg(1).folder, tmpIntg(1).name);%[prefix suffix];
        else
            filenameStruct.Intg='';
        end

        %% Orientation map
        tmpOrientationP2=dir([dataPathOrt '*OrientationP2*.mat']); %find with Orientation regex
        tmpOrientation=dir([dataPathOrt '*Orientation.mat']); %find with Orientation regex

        %[folderfilePrefixStr 'Orientation' PCAsuffix '.mat'];
        if ~isempty(tmpOrientationP2) %to get orientation map for analysis of vis+opto block 
            filenameStruct.Orientation=fullfile(tmpOrientationP2(1).folder, tmpOrientationP2(1).name);%[tmpOrientationP2.folder '\' tmpOrientationP2.name];
        elseif ~isempty(tmpOrientation) && isempty(tmpOrientationP2)
            filenameStruct.Orientation=fullfile(tmpOrientation(1).folder, tmpOrientation(1).name);%filenameStruct.Orientation=[dataPathOrt 'M' dataStruct.monkeyNo 'D' dataStruct.date 'R0Orientation.mat'];
        else 
            filenameStruct.Orientation='';
        end

        %% Static images
        filenameStruct.greenImgID=dataStruct.greenImgID;
        % Binned gain for quality check
        %{
        filenameStruct.greenStatic=getFirstFile(imagePath,'*EX540L_binned.tif');
        filenameStruct.gfpStatic=getFirstFile(imagePath,'*EX480*D*_binned_*.tif');
        filenameStruct.mcherryStatic=getFirstFile(imagePath,'*EX570*D*_binned*.tif');

        % Binned green image (for coregistration)
        filenameStruct.greenImg=convertCharsToStrings([imagePath 'green_binned.bmp']);
        if ~isfile(filenameStruct.greenImg) && isfile([imagePath 'SL30EX540L_binned.tif'])
            greenImage=imread((convertCharsToStrings([imagePath 'SL30EX540L_binned.tif'])));
            imwrite(greenImage,filenameStruct.greenImg);
        end
        %}

        %% Green images for coreg
        filenameStruct.greenServer=[filenameStruct.mainPath  dataStruct.monkey '/' dataStruct.monkey dataStruct.date '/green' num2str(dataStruct.greenImgID) '_binned.bmp'];
        filenameStruct.greenOI=['D:/' dataStruct.monkey dataStruct.date '/green' num2str(dataStruct.greenImgID) '_binned.bmp'];
        
       %% ROI timecourse
        tmpROITC=dir([dataPrefixStr '*ROITC.mat']); %find with ROITC regex
        filenameStruct.ROITC=[dataPrefixStr(1:end-numel(filePrefixStr)) tmpROITC.name];
        filenameStruct.vdaq=[dataPath 'Data_vdaqlog.mat']; % contains allegroimaging info
               
       %% Coregistration landmarks (manual) and transformation parameters (solution for coregistration)
        filenameStruct.landmarks=[imagePath 'landmarks.mat'];
        filenameStruct.transformParams=[imagePath 'transformParams.mat'];

        if isfield(dataStruct,'bmpOriginal') && isfield(dataStruct,'responseOriginal')
            filenameStruct.bmpOriginal=dataStruct.bmpOriginal;
            filenameStruct.responseOriginal=dataStruct.responseOriginal;
        end

       %% 2D Gaussian fit
        filenameStruct.gaussianFit=[imagePath 'gaussianFit.mat'];

        if ~isempty(dataStruct.gaussianResponse)  & ~isnan(dataStruct.gaussianResponse)
            % Find the last occurrence of '/'
            lastSlashIndex = find(dataStruct.gaussianResponse == '/', 3, 'last');
            % Extract the substring from the last '/' to the end
            extractedString = dataStruct.gaussianResponse(lastSlashIndex(1)+1:end);
            filenameStruct.gaussianResponse=[filenameStruct.mainPath dataStruct.monkey '/'  extractedString];
        else
            filenameStruct.gaussianResponse=[filenameStruct.mainPath dataStruct.monkey '/gaussianResponse.mat'];
        end
            
       %% Save file names
        filenameStruct.colmapPDF=[filenameStruct.plotPath 'M' dataStruct.monkeyNo uniqueStr dataStruct(1).date PCAsuffix '.pdf'];
        filenameStruct.psychneuroPDF=[filenameStruct.mainPath 'Chip/Meta/summary/' 'M' dataStruct.monkeyNo uniqueStr dataStruct(1).date  'R' dataStruct.run 'Summary.pdf'];%[filenameStruct.plotPath 'M' dataStruct.monkeyNo uniqueStr dataStruct(1).date 'NeurometricSummary.pdf'];
        filenameStruct.colmapMat=[filenameStruct.plotPath 'M' dataStruct.monkeyNo uniqueStr dataStruct(1).date PCAsuffix '.mat'];
        filenameStruct.psychfitPDF=[filenameStruct.mainPath 'Chip/Meta/summary/' 'M' dataStruct.monkeyNo 'Psychometrics.pdf'];%[filenameStruct.plotPath 'M' dataStruct.monkeyNo uniqueStr dataStruct(1).date 'NeurometricSummary.pdf'];
        filenameStruct.recruitmentPDF=[filenameStruct.mainPath 'Chip/Meta/recruitment/' 'M' dataStruct.monkeyNo uniqueStr dataStruct(1).date  'R' dataStruct.run 'Recruitment.pdf'];%[filenameStruct.plotPath 'M' dataStruct.monkeyNo uniqueStr dataStruct(1).date 'NeurometricSummary.pdf'];



    case {'VSD'}
        tmpFFT=dir([dataPrefixStr '*FFTAmp*0.mat']);
        filenameStruct.FFTAmp=[dataPrefixStr(1:end-numel(filePrefixStr)) tmpFFT.name];
        filenameStruct.FFTAmp
        %filenameStruct.FFTAmp=[basicStr 'DFFTAmpS016E115PF0500.mat'];
        filenameStruct.Orientation=[dataPrefixStr 'Orientation.mat'];
        filenameStruct.TS=[dataPrefixStr 'TS.mat'];
        filenameStruct.greenImg=convertCharsToStrings([imagePath 'SL30EX540L_binned_gainX5.tif']);
        filenameStruct.redImg=convertCharsToStrings([imagePath 'SL00Red_binned.bmp']);
        filenameStruct.initTransformParams=convertCharsToStrings([imagePath 'initTransformParams.mat']);
end
end