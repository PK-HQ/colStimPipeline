function [movingImgMaskCoregistered,coregStats,transformParams]=coregisterGreenImages(ImgReference,ImgTarget,method,filenameStructSession,Mask)
coregMasked=1;
imgLength=size(ImgReference,1);
imgSize=[1 imgLength];

fixedImg=double(ImgTarget);
movingImg=double(ImgReference);

%% Critical step: filter for blood vessels
ImgTarget = FilterFermi2D(ImgTarget,0.01,inf,1);
ImgReference = FilterFermi2D(ImgReference,0.01,inf,1);

switch coregMasked
    case {0}
        fixedImgMask=fixedImg;
        movingImgMask=movingImg;
    case {1}
        %Mask image to focus on center
        [x,y]=find(Mask==~isnan(Mask));
        ROImask=NaN(size(Mask));
        yOffset=25;
        xOffset=25;
        ROImask(min(x)-xOffset:max(x)+xOffset,min(y)-yOffset:max(y)+yOffset)=1;
        movingImgMask=rescale(double(ImgReference).*ROImask);
        
        %manual mask
        %{
        xOffset=round(.025*imgLength); %0.1
        yOffset=round(0.05*imgLength); %~.05
        imgHalfLength=imgLength/2;
        boxHalfLength=round(.275*imgLength);

        centerROIx=[imgHalfLength-boxHalfLength:imgHalfLength+boxHalfLength]-xOffset;
        centerROIy=centerROIx-yOffset;
        Mask = zeros(size(ImgReference));
        Mask(centerROIy,centerROIx)=1;
        %}
        
end

%% Correct for gross misalignments
figure('name','Green image coregistration')
[hA,~]=tight_subplot(3,3,[.1 .1]);
% pre
axes(hA(1))
imgsc(fixedImg,'Target (fixed)')
axes(hA(2))
imgsc(movingImgMask,'Reference (moving)')
axes(hA(3))
imagesc(imfuse(fixedImg,movingImgMask,'ColorChannels','red-cyan')); axis square; title('Raw')

switch method
    case {'Manual'}
        % Load if there already exists a landmarks file, else open control
        % point selector
        %load('Y:\Chip\landmarks_fp_mp.mat');
        %fpReference=fp;
        %mpReference=mp;
        if isfile(filenameStructSession.landmarks)
            %[mp,fp] = cpselect((movingImgMask),(fixedImg),'Wait',true);

            load(filenameStructSession.landmarks);
            
            
        else
            [mp,fp] = cpselect((movingImgMask),(fixedImg),'Wait',true);
        end    
        % warp moving image with 1st transform estimate
        transformParams = fitgeotrans(mp,fp,'projective'); %or affine
        movingImgMaskCoregistered = imwarp(movingImgMask,transformParams,'OutputView',imref2d(size(fixedImg)));
        %save(filenameStructSession.landmarks,'fp','mp')
        
        %{
        regParams.PyramidLevels=3;
        optParams.MaximumIterations=300;
        %optimizer
        [optimizer,metric] = imregconfig('multimodal');
        optimizer.MaximumIterations = optParams.MaximumIterations;
        transformType='affine';
        transformParams = imregtform(double(movingImgCoreg1>0),double(fixedImg>0),transformType,optimizer,metric,'DisplayOptimization',false,...
        'PyramidLevels',regParams.PyramidLevels);
        movingImgCoreg2 = imwarp(movingImgCoreg1,transformParams,'OutputView',imref2d(size(fixedImg)));
        %}
        
    case {'Auto'}
        regParams.PyramidLevels=3;
        %optimizer
        [optimizer,metric] = imregconfig('multimodal');
        %optimizer.MaximumIterations = 500;
        optimizer.InitialRadius = 6.25*10^-4;
        optimizer.GrowthFactor = 1.05;
        
        % Stage 1 = fast rigid (rotate translate)
        transformType='translation';
        transformRigid = imregtform(movingImgMask,fixedImg,transformType,optimizer,metric,'DisplayOptimization',false,...
        'PyramidLevels',regParams.PyramidLevels);
        movingImgMaskRigid = imwarp(movingImgMask,transformRigid,'OutputView',imref2d(size(fixedImg)));
        
        axes(hA(4))
        imgsc(fixedImg,'Target (fixed, similarity)')
        axes(hA(5))
        imgsc(movingImgMask,'Reference (moving, similarity)')
        axes(hA(6))
        imagesc(imfuse(fixedImg,movingImgMaskRigid,'ColorChannels','red-cyan')); axis square; title('Similarity')

       
        % Stage 2 = slow affine (rotate translate scale skew)
        regParams.PyramidLevels=3;       
        optimizer.MaximumIterations = 1000;
        optimizer.InitialRadius = 6.25*10^-3;
        optimizer.GrowthFactor = 1.05;%1.05;        
        transformType='affine';
        transformParams = imregtform(movingImgMask,fixedImg,transformType,optimizer,metric,'DisplayOptimization',false,...
        'PyramidLevels',regParams.PyramidLevels,'InitialTransformation',transformRigid);
        movingImgMaskCoregistered = imwarp(movingImgMask,transformParams,'OutputView',imref2d(size(fixedImg)));
        
        axes(hA(7))
        imgsc(fixedImg,'Target (fixed, affine)')
        axes(hA(8))
        imgsc(movingImgMask,'Reference (moving, affine)')
        axes(hA(9))
        imagesc(imfuse(fixedImg,movingImgMaskCoregistered,'ColorChannels','red-cyan')); axis square; title('Affine')
        
end


%calc correlation, translate, rotate, scale
transformCorrPrecoreg = corrcoef(double(fixedImg),double(movingImg));
transformCorrPrecoregMasked = corrcoef(fixedImg,movingImgMask);
transformCorrPostCoreg = corrcoef(fixedImg,movingImgMaskCoregistered);
coregStats.SimilarityPrecoreg = transformCorrPrecoreg(1,2);
coregStats.SimilarityPrecoregMasked = transformCorrPrecoregMasked(1,2);

coregStats.Similarity = transformCorrPostCoreg(1,2);
coregStats.Translate = transformParams.T(3,[1,2]);
coregStats.Rotation = atan2(transformParams.T(2,1),transformParams.T(1,1))*180/pi;
coregStats.Scale = abs(transformParams.T(2,1)+1i*transformParams.T(1,1));

%{
figure
subplot(1,3,1)
imagesc(imfuse(ImgTarget,ImgReference,'ColorChannels','red-cyan')); axis square
addPix2MM(min(imgSize),max(imgSize),min(imgSize),max(imgSize),2,2,2);
title(sprintf('Pre-coregistration: Similarity r=%.2g',coregStats.SimilarityPrecoreg))

subplot(1,3,2)
imagesc(imfuse(fixedImg,movingImgCoreg,'ColorChannels','red-cyan')); axis square
addPix2MM(min(imgSize),max(imgSize),min(imgSize),max(imgSize),2,2,2);
title(sprintf('CP: Similarity r=%.2g',coregStats.SimilarityPrecoregMasked))

subplot(1,3,3)
imagesc(imfuse(fixedImg,movingImgCoreg2,'ColorChannels','red-cyan')); axis square
addPix2MM(min(imgSize),max(imgSize),min(imgSize),max(imgSize),2,2,2);
titleMethodAndSimilarity=sprintf('Auto: Similarity r=%.2g',method,coregStats.Similarity);
titleStats=sprintf('(X=%.2g, Y=%.2g, Rot=%.2g, Scale=%.2g)',coregStats.Translate(1),coregStats.Translate(2),coregStats.Rotation,coregStats.Scale);
title({titleMethodAndSimilarity,titleStats});     
%}


end