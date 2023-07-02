%% Correcting for projector vs camera offsets
%% ----- Get desired response: Project a columnar bmp to cortex, save the projected response ------
% Setup system with blank board facing camera. EM-Tdt600, full aperture
% Project col image with AllegroClient, and manually capture and save
% ('F8-preview')
set(0,'DefaultFigureWindowStyle','docked')

%% Get estimated bitmap
% set current session
currentSession='20230209';
thresholdBmpEst=99;
% set bmp folder (that you used to project)
referenceBMPfolder='20221104';
referenceBMPOrt=0;
referenceBMPGamma='G050';


% ----- Load original bmp -----t
bmpOriginalFilename=['X:\PK\ColSeries\' referenceBMPfolder '\O' num2str(referenceBMPOrt,'%05.f') 'HE1000' referenceBMPGamma 'T09900.bmp'];
bmpOriginal=double(imread(bmpOriginalFilename));

% ----- Load projected response, threshold it to get estimated bmp -----
responseOriginalfolder=['D:\Chip' currentSession '\'];
responseOriginalFilename=[responseOriginalfolder 'calib50.bmp'];%'cal100_binned.bmp']; %[responseOriginalfolder 'calibration_binned_gainX2.bmp'];
dataSaveFolder=responseOriginalfolder;
responseOriginal=double(imread(responseOriginalFilename));
responseOriginalHE=adapthisteq(rescale(responseOriginal,0,1),'NumTiles',size(responseOriginal)./16,'Range','full');     %histEqImg=histeq(rescaled,1000);
responseOriginalFilt=double(responseOriginalHE>prctile(responseOriginalHE(:),thresholdBmpEst));

% just datastruct for dummy input
entryNo=99;
dataStruct(entryNo).monkeyNo='28';
dataStruct(entryNo).monkey='Chip';
dataStruct(entryNo).date=referenceBMPfolder;
dataStruct(entryNo).run='0';
dataStruct(entryNo).site='30';
dataStruct(entryNo).modality='GCaMP';
dataStruct(entryNo).bmpOriginal=bmpOriginalFilename;
dataStruct(entryNo).responseOriginal=responseOriginalFilename;


[bmpEstimated]=convertForProjector(dataStruct(entryNo),responseOriginalFilt,referenceBMPOrt,...
1,1,1,'cam2proj',[],0);
% ----- Plot both original and estimated -----
figure
[hA,~]=tight_subplot(2,2);
axes(hA(1))
imgsc(bmpOriginal)
title('bmpOriginal')
axes(hA(2))
imgsc(bmpEstimated)
title('bmpEstimated')

%% Get transformation matrix
% ----- coregister estimated and original bitmap -----
[optimizer,metric] = imregconfig('multimodal');
optimizer.MaximumIterations = 1000;
regParams.PyramidLevels=3;
transformType='similarity';
alignmentTransform = imregtform(bmpEstimated,bmpOriginal,transformType,optimizer,metric,'DisplayOptimization',false,...
'PyramidLevels',regParams.PyramidLevels);
% ----- apply transform -----
bmpEstimatedTransformed = imwarp(bmpEstimated,alignmentTransform,'OutputView',imref2d(size(bmpEstimated)));

%{
transformType='similarity';
transformParams = imregtform(bmpEstimatedTransformed,bmpOriginal,transformType,optimizer,metric,'DisplayOptimization',false,...
'PyramidLevels',regParams.PyramidLevels);
% ----- apply transform -----
bmpEstimatedTransformed = imwarp(bmpEstimatedTransformed,transformParams,'OutputView',imref2d(size(bmpEstimated)));
%}
% ----- get transform params -----

transformCorrPost = corrcoef(bmpOriginal,bmpEstimatedTransformed);
transformCorrPre = corrcoef(bmpOriginal,bmpEstimated);
coregStats.Similarity = transformCorrPost(1,2);
coregStats.Translate = alignmentTransform.T(3,[1,2]);
coregStats.Rotation = atan2(alignmentTransform.T(2,1),alignmentTransform.T(1,1))*180/pi;
coregStats.Scale = abs(alignmentTransform.T(2,1)+1i*alignmentTransform.T(1,1));
sprintf('(X=%.3g, Y=%.3g, Rot=%.3g, Scale=%.3g)',coregStats.Translate(1),coregStats.Translate(2),coregStats.Rotation,coregStats.Scale)

subplot(2,2,3)
imagesc(imfuse(bmpOriginal,bmpEstimated,'ColorChannels','red-cyan')); axis equal
title(sprintf('Pre-coregistration, %.2f', coregStats.Similarity))

subplot(2,2,4)
imagesc(imfuse(bmpOriginal,bmpEstimatedTransformed,'ColorChannels','red-cyan')); axis equal
title(sprintf('Post-coregistration, %.2f', coregStats.Similarity))
1;


%% Apply to c
transformParamsOriginal=alignmentTransform;

[bmpEstimated]=convertForProjector(dataStruct(entryNo),responseOriginalFilt,referenceBMPOrt,...
1,1,1,'cam2proj',transformParamsOriginal,0);

% ----- Plot both original and estimated -----
figure
subplot(2,2,1)
imgsc(bmpOriginal)
title('bmpOriginal')
subplot(2,2,2)
imgsc(bmpEstimated)
title('bmpEstimated')

%% Get transformation matrix
% ----- coregister estimated and original bitmap -----
[optimizer,metric] = imregconfig('multimodal');
optimizer.MaximumIterations = 100;
regParams.PyramidLevels=3;
transformType='similarity';
transformParamsControl = imregtform(bmpEstimated,bmpOriginal,transformType,optimizer,metric,'DisplayOptimization',false,...
'PyramidLevels',regParams.PyramidLevels);

% ----- apply transform -----
bmpEstimatedTransformed = imwarp(bmpEstimated,transformParamsControl,'OutputView',imref2d(size(bmpEstimated)));

% ----- get transform params -----
transformCorrPost = corrcoef(bmpOriginal,bmpEstimatedTransformed);
transformCorrPre = corrcoef(bmpOriginal,bmpEstimated);
coregStats.Similarity = transformCorrPost(1,2);
coregStats.Translate = transformParamsControl.T(3,[1,2]);
coregStats.Rotation = atan2(transformParamsControl.T(2,1),transformParamsControl.T(1,1))*180/pi;
coregStats.Scale = abs(transformParamsControl.T(2,1)+1i*transformParamsControl.T(1,1));
sprintf('(X=%.3g, Y=%.3g, Rot=%.3g, Scale=%.3g)',coregStats.Translate(1),coregStats.Translate(2),coregStats.Rotation,coregStats.Scale)

subplot(2,2,3)
imagesc(imfuse(bmpOriginal,bmpEstimated,'ColorChannels','red-cyan')); axis equal
title(sprintf('[Should be perfect] Pre-coregistration, %.2f', coregStats.Similarity))

subplot(2,2,4)
imagesc(imfuse(bmpOriginal,bmpEstimatedTransformed,'ColorChannels','red-cyan')); axis equal
title(sprintf('[Sanity] Post-coregistration, %.2f', coregStats.Similarity))
1;


save([dataSaveFolder 'alignmentTransform.mat'],'alignmentTransform')

%% ----- Apply transform params to bitmap correction pipeline -----
%{
[projBitmapTRBB]=convertForProjector(dataStruct(4),pipelined,1,...
1,1,1,'calibrate');
[optimizer,metric] = imregconfig("multimodal");
imregtform(double(orig),pipelined,"similarity");
regParams.PyramidLevels=3;
optParams.MaximumIterations=500;

titleStats=sprintf('(X=%.2g, Y=%.2g, Rot=%.2g, Scale=%.2g)',coregStats.Translate(1),coregStats.Translate(2),coregStats.Rotation,coregStats.Scale);
transformCorrPre = corrcoef(double(orig),pipelined);
transformCorrPost = corrcoef(double(orig),bmpEstimatedTransformed);
bmpEstimatedTransformed = imwarp(pipelined,transformParams,'OutputView',imref2d(size(pipelined)));
%}





