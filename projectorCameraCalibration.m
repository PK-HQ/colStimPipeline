%% Correcting for projector vs camera offsets
%% ----- Get desired response: Project a columnar bmp to cortex, save the projected response ------
% Setup system with blank board facing camera. EM-Tdt600, full aperture
% Project col image with AllegroClient, and manually capture and save
% ('F8-preview')
set(0,'DefaultFigureWindowStyle','docked')

%% Get estimated bitmap
% set current session
currentSession='calibration20240105';
thresholdBmpEst=97;%99;
% set bmp folder (that you used to project)
referenceBMPfolder='20221104';%'20230802';%'20221104/20230804';
referenceBMPOrt=0;
referenceBMPGamma='G050'; %'G050/G120'

% ----- Load original bmp -----t
bmpOriginalFilename=['V:\PK\ColSeries\' referenceBMPfolder '\O' num2str(referenceBMPOrt,'%05.f') 'HE1000' referenceBMPGamma 'T09900.bmp'];%['V:\PK\ColSeries\' referenceBMPfolder '\O' num2str(referenceBMPOrt*100,'%05.f') 'HE0032' referenceBMPGamma 'S00001.bmp'];
%['V:\PK\ColSeries\' referenceBMPfolder '\O' num2str(referenceBMPOrt,'%05.f') 'HE1000' referenceBMPGamma 'T09900.bmp'];
bmpOriginal=double(imread(bmpOriginalFilename));

% ----- Load projected response, threshold it to get estimated bmp -----
responseOriginalfolder=['D:\OIdata\' currentSession '\'];%['D:\Chip' currentSession '\'];
responseOriginalFilename=[responseOriginalfolder 'calib100_' referenceBMPfolder '.bmp'];% 'calib50.bmp'];
dataSaveFolder=responseOriginalfolder;
responseOriginal=double(imread(responseOriginalFilename));
responseOriginalHE=adapthisteq(rescale(responseOriginal,0,1),'NumTiles',size(responseOriginal)./16,'Range','full');     %histEqImg=histeq(rescaled,1000);
responseOriginalFilt=double(responseOriginalHE>=prctile(responseOriginalHE(:),thresholdBmpEst));

% just datastruct for dummy input
entryNo=99;
dataStruct(entryNo).monkeyNo='28';
dataStruct(entryNo).monkey='Chip';
dataStruct(entryNo).date=referenceBMPfolder;
dataStruct(entryNo).run='0';
dataStruct(entryNo).site='30';
dataStruct(entryNo).baselineTS='0';
dataStruct(entryNo).modality='GCaMP';
dataStruct(entryNo).bmpOriginal=bmpOriginalFilename;
dataStruct(entryNo).responseOriginal=responseOriginalFilename;
dataStruct(entryNo).gaussianContourLevel=[1 1];
dataStruct(entryNo).gaussianContourLevelMax=200; 
%[bmpEstimated]=convertForProjector(dataStruct(entryNo),responseOriginalFilt,referenceBMPOrt,...
%1,1,1,'cam2proj',[],0);

[bmpEstimated,~,~]=convertForProjector(dataStruct(entryNo),dataStruct(entryNo),responseOriginalFilt,referenceBMPOrt,...
    1,1,1,'cam2proj',[],0,0);

% ----- Plot both original and estimated -----
figure
[hA,~]=tight_subplot(2,3);
axes(hA(1))
imgsc(bmpOriginal)
title('bmpOriginal')
axes(hA(2))
imgsc(bmpEstimated)
title('bmpEstimated')

%% Get transformation matrix
bmpEstimated=rescale(double(bmpEstimated)); %.*ROImask

axes(hA(3))
imagesc(imfuse(bmpOriginal,bmpEstimated,'ColorChannels','red-cyan'));
title(sprintf('Pre-coregistration, %.2f', coregStats.Similarity))

%optimizer
[optimizer,metric] = imregconfig('multimodal');

% Stage 1 = fast rigid (rotate translate)
%{
%optimizer.MaximumIterations = 500;
transformType='rigid';
transformCoarse = imregtform(bmpEstimated,bmpOriginal,transformType,optimizer,metric,'DisplayOptimization',false);
bmpEstimatedTransformedRigid = imwarp(bmpEstimated,transformCoarse,'OutputView',imref2d(size(bmpEstimated)));
axes(hA(4))
imagesc(imfuse(bmpOriginal,bmpEstimatedTransformedRigid,'ColorChannels','red-cyan'));
title(sprintf('Rigid'));
1;
%}
% Stage 2 = slow affine (rotate translate scale skew)
%optimizer.MaximumIterations = 1000;

transformType='similarity';
alignmentTransformSim = imregtform(bmpEstimated,bmpOriginal,transformType,optimizer,metric,'DisplayOptimization',false);
bmpEstimatedTransformedSim = imwarp(bmpEstimated,alignmentTransformSim,'OutputView',imref2d(size(bmpEstimated)));
axes(hA(5))
imagesc(imfuse(bmpOriginal,bmpEstimatedTransformedSim,'ColorChannels','red-cyan'));
title(sprintf('Similarity'));
1;

transformType='affine';
alignmentTransformAff = imregtform(bmpEstimated,bmpOriginal,transformType,optimizer,metric,'DisplayOptimization',false);
bmpEstimatedTransformed = imwarp(bmpEstimated,alignmentTransformAff,'OutputView',imref2d(size(bmpEstimated)));
axes(hA(6))
imagesc(imfuse(bmpOriginal,bmpEstimatedTransformed,'ColorChannels','red-cyan'));
title(sprintf('Affine'))
1;

%% Save
alignmentTransform=alignmentTransformSim;
save([dataSaveFolder 'alignmentTransform.mat'],'alignmentTransform')



