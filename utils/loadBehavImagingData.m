function [behavioralData,imagingData,successFlag] = loadBehavImagingData(currentEntry, referenceEntry, alignmentEntry, behavioralData, imagingData, blockID, pipelineMode)

successFlag = true;

switch pipelineMode
    case 'expt'  % Load minimal components for bitmap generation
        % 1. Reference map components
        load(referenceEntry.Orientation, 'Mask', 'RespCondPCA', 'Ort', 'MapAmpOrt', 'PCAExpl', 'nPCAComp');
        
        % 2. Camera-projector alignment transform
        if isfile(fullfile(alignmentEntry, 'alignmentTransform.mat'))
            load(fullfile(alignmentEntry, 'alignmentTransform.mat'), 'alignmentTransform');
        else
            error('Missing alignment transform file');
        end
        
        % 3. Reference map to current FOV transform
        %load(currentEntry.transformParams, 'transformParams');
        %imagingData.transformRef2Block(blockID) = transformParams;
        
        % 4. Gaussian response for selecting column subsets
        load(referenceEntry.Orientation, 'Mask', 'RespCondPCA', 'Ort', 'MapAmpOrt', 'PCAExpl','nPCAComp');
        % Gaussian Fit
        if isfile(currentEntry.gaussianFit)
            load(currentEntry.gaussianFit,'ROIMaskgaussian','centerCoords');
        else
            ROIMaskgaussian=[];
            centerCoords=[];
        end

        if exist(currentEntry.gaussianResponse)==2 && isfile(currentEntry.gaussianResponse)
            load(currentEntry.gaussianResponse,'FFTCond');
        else
            FFTCond=nan(128,128,3);
        end

        
        % 5. Compile essential data into structs
        imagingData.ortpca(:,:,:,blockID) = RespCondPCA;
        imagingData.orts(:,:,blockID) = unique(Ort);
        imagingData.mask(:,:,blockID) = double(Mask);
        imagingData.nanmask(:,:,blockID) = imagingData.mask(:,:,blockID); 
        imagingData.nanmask(imagingData.nanmask==0) = NaN;
        imagingData.ortampmap(:,:,:,blockID) = imagingData.mask(:,:,blockID).*MapAmpOrt;
        imagingData.npca(1,blockID) = nPCAComp;
        imagingData.pcaexpl(:,blockID) = PCAExpl(1:12);
        imagingData.gaussresp(1:128,1:128,1:3,blockID)=padArray(FFTCond, 3, 3, NaN);
        imagingData.gaussfit(:,:,blockID)=padStruct(ROIMaskgaussian, 200);
        if isempty(centerCoords)
            imagingData.centerCoords(1,:,blockID)=[NaN NaN];
        end

        % Camera parameters
        imagingData.pixelsizemm(1,blockID) = 0.016; % Standard pixel size in mm
        imagingData.pixels(1,blockID) = 512;  % Standard image size
        
        % Camera-projector transformation
        imagingData.transformmatrix(:,blockID) = alignmentTransform;
        
    otherwise  % Original full functionality
        % Init
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
            load(currentEntry.gaussianFit,'ROIMaskgaussian','centerCoords');
        end

        if exist(currentEntry.gaussianResponse)==2 && isfile(currentEntry.gaussianResponse)
            load(currentEntry.gaussianResponse,'FFTCond');
        else
            FFTCond=nan(128,128,3);
        end

        % Proj2cam
        if isfile(fullfile(alignmentEntry, 'alignmentTransform.mat'))
            load(fullfile(alignmentEntry, 'alignmentTransform.mat'), 'alignmentTransform');
        end

        % refmap FOV to current FOV
        load(currentEntry.transformParams, 'transformParams');
        imagingData.transformRef2Block(blockID)=transformParams;

        if ~isempty(currentEntry.Intg)
            load(currentEntry.Intg, 'DataTrial');      
            if isempty(DataTrial)
                DataTrial=[];
                imagingData.optoIntg(:,:,:,blockID)= DataTrial;
            else
                disp(['IntgOpto file loaded:' currentEntry.Intg])
                imagingData.optoIntg(:,:,:,blockID)= DataTrial;
            end
        else
            disp(['=== IntgOpto missing'])
        end

        if ~isempty(currentEntry.baselineIntg)
            load(currentEntry.baselineIntg, 'DataTrial');
            if isempty(DataTrial)
                DataTrial=[];
                imagingData.baselineIntg(:,:,:,blockID)= DataTrial;
            else
                disp(['IntgBL file loaded:' currentEntry.baselineIntg])
                imagingData.baselineIntg(:,:,:,blockID)= DataTrial;
              %% Check if opto intg loaded 
                if size(imagingData.optoIntg,3)>=blockID
                    if isempty(imagingData.optoIntg(:,:,:,blockID))
                        disp(['Check INTG: Block ' num2str(blockID)])
                    end
                elseif size(imagingData.optoIntg,3)<blockID
                    disp(['Check INTG: Block ' num2str(blockID)])
                end
            end
        else
            disp(['=== IntgBL missing'])
        end

        % Compile into struct
        validRefConds=behavioralData.referenceTS(blockID).Header.Conditions.TypeCond(behavioralData.referenceTS(blockID).Header.Conditions.TypeCond>0);
        validRefConds=validRefConds==2;

        % Ort map
        imagingData.ortpca(:,:,:,blockID)=RespCondPCA(:,:,validRefConds);
        imagingData.orts(:,:,blockID)=unique(Ort);
        imagingData.mask(:,:,blockID)=double(Mask);
        imagingData.nanmask(:,:,blockID)=imagingData.mask(:,:,blockID); 
        imagingData.nanmask(imagingData.nanmask==0)=NaN;
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
end
end