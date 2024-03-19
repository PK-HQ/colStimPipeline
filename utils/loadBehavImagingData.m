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
    load(currentEntry.TS, 'TS');
catch ME
    % If an error occurs, print the error message and set successFlag to 0
    disp(['Error loading file: ', ME.message]);
    successFlag = false;
end
 % Reference map
try
    refTS=load(referenceEntry.TS, 'TS');
catch ME
    % If an error occurs, print the error message and set successFlag to 0
    disp(['Error loading file: ', ME.message]);
    successFlag = false;
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

if isfile(fullfile(alignmentEntry, 'alignmentTransform.mat'))
    load(fullfile(alignmentEntry, 'alignmentTransform.mat'), 'alignmentTransform');
end
%{
if ~isempty(currentEntry.Intg)
        % Intg data
        load(currentEntry.Intg, 'DataTrial');      
        % Imaging block integrated data (300ms window, frames 4:10)
        if isempty(DataTrial)
            DataTrial=[];
        end
        imagingData.intg(:,:,:,blockID)=padArray(DataTrial, 400, 3, NaN);
end
%}
% Compile into struct
validRefConds=refTS.TS.Header.Conditions.TypeCond(refTS.TS.Header.Conditions.TypeCond>0);
validRefConds=validRefConds==2;
% TS 
behavioralData.TS=TS;

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
imagingData.pixelsizemm(1,blockID)=refTS.TS.Header.Imaging.SizePxl;
imagingData.pixels(1,blockID)=refTS.TS.Header.Imaging.FrameWidth;
% Cam-proj transformation
imagingData.transformmatrix(:,blockID)=alignmentTransform;
end