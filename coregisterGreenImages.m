function [movingImgMaskCoregistered,coregStats,transformParams,ImgReference,ImgTarget]=coregisterGreenImages(...
  currentBlockStruct,referenceBlockStruct,Mask)
%% [README] Loads green images of reference and current session, and 
%%% coregister with auto (same session) or manual mode (diff session or override)

%% load green images
[ImgReference,ImgTarget]=loadGreenImg(currentBlockStruct,referenceBlockStruct);

%% Coreg
% Get coregistration params for greenImageReference (to be transformed) to greenImageSession (anchor image)
sameSession=isequal(currentBlockStruct.date,referenceBlockStruct.date);
sameSession=1; %OVERRIDE
switch sameSession
    case {1}
        method='auto'; %manual/auto
    case {0}
        method='manual';
end

coregMask=1;
autoStages=2;

ImgTarget=double(ImgTarget);
movingImg=double(ImgReference);

%% Critical step: filter for blood vessels
ImgTarget = FilterFermi2D(ImgTarget,0.01,inf,1);
ImgReference = FilterFermi2D(ImgReference,0.01,inf,1);

switch coregMask
    case {0}
        movingImgMask=movingImg;
    case {1}
        %Mask image to focus on center
        [x,y]=find(Mask==~isnan(Mask));
        ROImask=zeros(size(Mask));
        yOffset=50;
        xOffset=50;
        movingImgMask=rescale(double(ImgReference)); %.*ROImask
        ImgTarget=rescale(double(ImgTarget));        
end

%% Correct for gross misalignments
figure('name','Vasculature coregistration')
[hA,~]=tight_subplot(2,2);%1+(1*autoStages),3,[.05 .05]);
% pre
axes(hA(1))
imgsc(ImgTarget,'Session green image (fixed)'); upFontSize(14,.015); colormap(gray)
axes(hA(2))
imgsc(movingImgMask,'Reference green image (moving)'); upFontSize(14,.015); colormap(gray)
axes(hA(3))
imagesc(imfuse(ImgTarget,movingImgMask,'ColorChannels','red-cyan')); axis square; title('Superimposed (Pre-coregistration)','FontWeight','Normal'); upFontSize(14,.015)

switch method
    case {'manual'}
        polyExponent=2;
        if isfile(referenceBlockStruct.transformParams)
            % load transformaiton file if exist
            load(referenceBlockStruct.transformParams);
        else
            % Load if there already exists a landmarks file, else open control point selector
            if isfile(referenceBlockStruct.landmarks)
                load(referenceBlockStruct.landmarks);
            else
                [mp,fp] = cpselect((movingImgMask),(ImgTarget),'Wait',true);
                save(referenceBlockStruct.landmarks,'fp','mp')
            end    
            % warp moving image with transform estimate
            transformParams = fitgeotrans(mp,fp,'polynomial',polyExponent); %fitgeotrans(mp,fp,'affine')%fitgeotrans(mp,fp,'polynomial',2) %polynomial/projective/affine
            % save transformation file
            save(referenceBlockStruct.transformParams,'transformParams')
        end
        titleStr=['Superimposed (Polynomial=' num2str(polyExponent) ')'];
        
    case {'auto'}
        switch autoStages
            case {2}
                if isfile(referenceBlockStruct.transformParams)
                    titleStr='Superimposed (similarity + affine)';
                    % load transformaiton file if exist
                    load(referenceBlockStruct.transformParams);
                else 
                    titleStr='Superimposed (similarity + affine)';

                    % set levels
                    regParams.PyramidLevels=3;

                    %optimizer
                    [optimizer,metric] = imregconfig('multimodal');

                    % Stage 1 = fast rigid (rotate translate)
                    %optimizer.MaximumIterations = 500;
                    optimizer.InitialRadius = 6.25*10^-4;
                    optimizer.GrowthFactor = 1.05;
                    transformType='rigid';
                    transformCoarse = imregtform(movingImgMask,ImgTarget,transformType,optimizer,metric,'DisplayOptimization',false,...
                    'PyramidLevels',regParams.PyramidLevels);
                    movingImgMaskCoarse = imwarp(movingImgMask,transformCoarse,'OutputView',imref2d(size(ImgTarget)));

                    % Stage 2 = slow affine (rotate translate scale skew)
                    %optimizer.MaximumIterations = 1000;
                    regParams.PyramidLevels=3;       
                    optimizer.InitialRadius = 6.25*10^-3;
                    optimizer.GrowthFactor = 1.05;       
                    transformType='affine';
                    transformParams = imregtform(movingImgMask,ImgTarget,transformType,optimizer,metric,'DisplayOptimization',false,...
                    'PyramidLevels',regParams.PyramidLevels,'InitialTransformation',transformCoarse);

                    % Previous function used
                    %[tformCoarse,movingImgMaskDemons]=imregdemons(movingImgMask,ImgTarget,'PyramidLevels',5,'AccumulatedFieldSmoothing',3);

                    % save transformation file
                    %save(referenceBlockStruct.transformParams,'transformParams')

                    %{
                    if isfile(filenameStructSession.transformParams)
                        % load transformaiton file if exist
                        load(filenameStructSession.transformParams);
                    else
                        regParams.PyramidLevels=3;

                        %optimizer
                        [optimizer,metric] = imregconfig('multimodal');

                        % Stage 1 = fast rigid (rotate translate)
                        %optimizer.MaximumIterations = 500;
                        optimizer.InitialRadius = 6.25*10^-4;
                        optimizer.GrowthFactor = 1.05;
                        transformType='rigid';
                        transformCoarse = imregtform(movingImgMask,ImgTarget,transformType,optimizer,metric,'DisplayOptimization',false,...
                        'PyramidLevels',regParams.PyramidLevels);
                        movingImgMaskCoarse = imwarp(movingImgMask,transformCoarse,'OutputView',imref2d(size(ImgTarget)));

                        % Stage 2 = slow affine (rotate translate scale skew)
                        %optimizer.MaximumIterations = 1000;
                        regParams.PyramidLevels=3;       
                        optimizer.InitialRadius = 6.25*10^-3;
                        optimizer.GrowthFactor = 1.05;       
                        transformType='affine';
                        transformParams = imregtform(movingImgMask,ImgTarget,transformType,optimizer,metric,'DisplayOptimization',false,...
                        'PyramidLevels',regParams.PyramidLevels,'InitialTransformation',transformCoarse);

                        % Previous function used
                        %[tformCoarse,movingImgMaskDemons]=imregdemons(movingImgMask,ImgTarget,'PyramidLevels',5,'AccumulatedFieldSmoothing',3);

                        % save transformation file
                        save(filenameStructSession.transformParams,'transformParams')
                    end
                      %}
                     upFontSize(14,.015)
                    1;
                end
            case {1}
                titleStr='Superimposed (affine)';
                if isfile(referenceBlockStruct.transformParams)
                    % load transformaiton file if exist
                    load(referenceBlockStruct.transformParams);
                else
                    % Stage 2 = slow affine (rotate translate scale skew)
                    regParams.PyramidLevels=3;                      
                    %optimizer
                    [optimizer,metric] = imregconfig('multimodal');
                    optimizer.MaximumIterations = 3000;
                    optimizer.InitialRadius = 6.25*10^-3;
                    optimizer.GrowthFactor = 1.05;%1.05;        
                    transformType='similarity';
                    transformParams = imregtform(movingImgMask,ImgTarget,transformType,optimizer,metric,'DisplayOptimization',false,...
                    'PyramidLevels',regParams.PyramidLevels);
                    
                    % save transformation file
                    save(referenceBlockStruct.transformParams,'transformParams')
                end
        end
end


movingImgMaskCoregistered = imwarp(movingImgMask,transformParams,'OutputView',imref2d(size(ImgTarget)));
axes(hA(4))
imagesc(imfuse(ImgTarget,movingImgMaskCoregistered,'ColorChannels','red-cyan')); axis square; title(titleStr); upFontSize(14,.015)

%calc correlation, translate, rotate, scale
transformCorrPrecoreg = corrcoef(double(ImgTarget),double(movingImg));
transformCorrPrecoregMasked = corrcoef(ImgTarget,movingImgMask);
transformCorrPostCoreg = corrcoef(ImgTarget,movingImgMaskCoregistered);
coregStats.SimilarityPrecoreg = transformCorrPrecoreg(1,2);
coregStats.SimilarityPrecoregMasked = transformCorrPrecoregMasked(1,2);

coregStats.Similarity = transformCorrPostCoreg(1,2);
%coregStats.Translate = transformParams.T(3,[1,2]);
%coregStats.Rotation = atan2(transformParams.T(2,1),transformParams.T(1,1))*180/pi;
%coregStats.Scale = abs(transformParams.T(2,1)+1i*transformParams.T(1,1));

end