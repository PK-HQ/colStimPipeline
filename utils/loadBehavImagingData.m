function [behavioralData,imagingData,successFlag] = loadBehavImagingData(currentEntry, referenceEntry, alignmentEntry, behavioralData, imagingData, blockID)
% Init;
TS=[];
DataTrial=[];
RespCondPCA=[];
Ort=[];
Mask=[];
ROIMaskgaussian=[];
FFTCond=[];
alignmentTransform=[];
MapAmpOrt=[];
PCAExpl=[];
centerCoords=[];
% TS
successFlag = true;

try
    load(currentEntry.Stab, 'refcfg');
    % Stabilizer mat
    xMvmt=mean(refcfg.dx);
    yMvmt=mean(refcfg.dy);
    imagingData.Stab(blockID,:)=[xMvmt yMvmt];
end
try
    load(currentEntry.TS, 'TS');
    % TS 
    behavioralData.optoTS(blockID)=TS;
end
try
    load(currentEntry.baselineTS, 'TS');
    % TS 
    behavioralData.baselineTS(blockID)=TS;
end

 % Reference map
try
    load(referenceEntry.TS, 'TS');
    behavioralData.referenceTS(blockID)=TS;
end

load(referenceEntry.Orientation, 'Mask', 'RespCondPCA', 'Ort', 'MapAmpOrt', 'PCAExpl','nPCAComp');
% Gaussian Fit
if isfile(currentEntry.gaussianFit)
    load(currentEntry.gaussianFit,'ROIMaskgaussian','centerCoords'); % .gaussian file
end

if exist(currentEntry.gaussianResponse)==2 && isfile(currentEntry.gaussianResponse)
    load(currentEntry.gaussianResponse,'FFTCond'); % .gaussian file
else
    FFTCond=nan(128,128,3);
end
%% Proj2cam
if isfile(fullfile(alignmentEntry, 'alignmentTransform.mat'))
    load(fullfile(alignmentEntry, 'alignmentTransform.mat'), 'alignmentTransform');
end
%% refmap FOV to current FOV
load(currentEntry.transformParams, 'transformParams');
imagingData.transformRef2Block(blockID)=transformParams;

if ~isempty(currentEntry.Intg) % Load intg imaging data (either optostim or optostim+baseline, latter if both are in the same block)
        % Intg data
        load(currentEntry.Intg, 'DataTrial');      
        % Imaging block integrated data (300ms window, frames 4:10)
        if isempty(DataTrial)
            DataTrial=[];
            imagingData.optoIntg(:,:,:,blockID)= DataTrial;%padArray2(DataTrial, 400, 3, NaN, 'pad');%padArray(DataTrial, 400, 3, NaN); % pad ntrials in z-dimension to 400

        else
            disp(['IntgOpto file loaded:' currentEntry.Intg])
            imagingData.optoIntg(:,:,:,blockID)= DataTrial;%padArray2(DataTrial, 400, 3, NaN, 'pad');%padArray(DataTrial, 400, 3, NaN); % pad ntrials in z-dimension to 400
        end
else
    disp(['=== IntgOpto missing:' currentEntry.OptoTS])
end


if ~isempty(currentEntry.baselineIntg) % Load baseline data if optostim+baseline are not in the same block
    % Intg data
    load(currentEntry.baselineIntg, 'DataTrial');
    % Imaging block integrated data (300ms window, frames 4:10)
    if isempty(DataTrial)
        DataTrial=[];
        imagingData.baselineIntg(:,:,:,blockID)= DataTrial;%padArray(DataTrial, 400, 3, NaN); % pad ntrials in z-dimension to 400
    else
        disp(['IntgBL file loaded:' currentEntry.baselineIntg])
        imagingData.baselineIntg(:,:,:,blockID)= DataTrial;%padArray(DataTrial, 400, 3, NaN); % pad ntrials in z-dimension to 400
    end
else
    disp(['=== IntgBL missing:' currentEntry.baselineTS])
end


% Compile into struct
validRefConds=behavioralData.referenceTS(blockID).Header.Conditions.TypeCond(behavioralData.referenceTS(blockID).Header.Conditions.TypeCond>0);
validRefConds=validRefConds==2;

% Ort map
imagingData.ortpca(:,:,:,blockID)=RespCondPCA(:,:,validRefConds);
imagingData.orts(:,:,blockID)=unique(Ort);
imagingData.mask(:,:,blockID)=double(Mask);
imagingData.nanmask(:,:,blockID)=imagingData.mask(:,:,blockID); imagingData.nanmask(imagingData.nanmask==0)=NaN;
imagingData.ortampmap(:,:,:,blockID)=imagingData.mask(:,:,blockID).*MapAmpOrt;
imagingData.npca(1,blockID)=nPCAComp;
imagingData.pcaexpl(:,blockID)=PCAExpl(1:12);

% Gaussian fitting
imagingData.gaussresp(1:128,1:128,1:3,blockID)=padArray(FFTCond, 3, 3, NaN);
imagingData.gaussfit(:,:,blockID)=padStruct(ROIMaskgaussian, 200);
if isempty(centerCoords)
    imagingData.centerCoords(1,:,blockID)=[NaN NaN];
end

% Cam pixel sizes
imagingData.pixelsizemm(1,blockID)=behavioralData.referenceTS(blockID).Header.Imaging.SizePxl;
imagingData.pixels(1,blockID)=behavioralData.referenceTS(blockID).Header.Imaging.FrameWidth;
% Cam-proj transformation
imagingData.transformmatrix(:,blockID)=alignmentTransform;


%% Check if opto intg loaded 
if size(imagingData.optoIntg,3)>=blockID
    if isempty(imagingData.optoIntg(:,:,:,blockID))
        disp(['Check INTG: D' datastruct(analysisBlockID(blockID)).date 'R' datastruct(analysisBlockID(blockID)).run])
        %errorIntgBlocks{end+1}= ['D' datastruct(analysisBlockID(blockID)).date 'R' datastruct(analysisBlockID(blockID)).run];
    end
elseif size(imagingData.optoIntg,3)<blockID
    disp(['Check INTG: D' datastruct(analysisBlockID(blockID)).date 'R' datastruct(analysisBlockID(blockID)).run])
    %errorIntgBlocks{end+1}= ['D' datastruct(analysisBlockID(blockID)).date 'R' datastruct(analysisBlockID(blockID)).run];
end
end