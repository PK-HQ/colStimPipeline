%% Define analysis and current session to be analysed
pipelineMode='';% beta
analysisMode='summary'; %summary/psychometrics/beta/neurometric/stability/neurometricsfix
currentSessID=112;%for biasing expt

% Save, plot flags
saveFlag=0;
saveFlagBMP=0; 
plotFlag=0;
chambers={'R', 'L'};
%% Load dataStruct for the desired chamber
[mainPath, datastruct]=setupEnv('users/PK/colStimPipeline/exptListBiasingFull.m');
for chamberID=1%1:2
    nColumnsWanted=[]; chamberWanted=chambers{chamberID};
    analysisBlockID = organizeBlocks(datastruct, chamberWanted, nColumnsWanted);

    %% Pipeline: Bitmap generation
    switch pipelineMode
        case {'beta'}        
            %% Fast pipeline (coregisters within session ort map to optostim block green image)
            % Get correction for camera-projector alignment
            % imregtform of orignal and recovered bitmap, apply to correct for
            % camera-projector alignment
            open projectorCameraCalibration
            
            % Load projector-camera alignment transformation matrix
            load([alignmentSessionFolder 'alignmentTransform.mat'],'alignmentTransform')
    
            % Define reference (orientation map) session and current session
            % (for optostim experiment)
            dsReferenceSess=datastruct(referenceSessID); %moving reference ort map
            dsCurrentSess=datastruct(currentSessID); %fixed current session map, which ort map will be projected onto
            filenameStructCurrent=generateFilenames(dsCurrentSess);
            pdfFilename=filenameStructCurrent.neurometricPDF;
    
            % Get reference session orientation map
            [columnarBitmap,VERpca,columnarmapStats]=getColumnarBitmapV4(mainPath,dsReferenceSess,dsCurrentSess,bitmapParams, ...
              plotFlag, saveFlag, pdfFilename);
    
            % Coregister reference session orientation map to the cortical view of current session 
            [columnarBitmapCoregistered, columnarPCAsCoregistered]=coregisterBitmap2GreenImgV2(dsReferenceSess,dsCurrentSess, ...
              columnarBitmap,VERpca,plotFlag,saveFlag, pdfFilename);
            
            % create bitmap from orientation map, applying alignment matrix and
            % accounting for projector dimensions (i.e. a rectangle)
            orts=[0 90];%0:15:165;
            HE=1000;
            [projBitmapTRBB,bitmapsCamSpace,nBlobs,medianBlobAreas]=convertForProjector(dsReferenceSess,dsCurrentSess,columnarBitmapCoregistered,orts,...
                bitmapParams.gridSize,bitmapParams.gammaCorrFactor,bitmapParams.sensitivity,'cam2proj',alignmentTransform,plotFlag,saveFlagBMP,saveFlag, pdfFilename);
    end
    
    %% Run the desired analysis pipeline variant
    switch analysisMode
        case {'summary'}
            if ~exist('behavioralData','var')
                behavioralData=[];
                imagingData=[];
                bitmapData=[];
            end
            
            for blockID=1:numel(analysisBlockID)
                disp(['=== Block ' num2str(blockID)  '/' num2str(numel(analysisBlockID)) '==='])
                tic
                if isfield(behavioralData,'auc') && size(behavioralData.auc,3)>=blockID
                     disp('(Skipping completed block)')
                     continue
                end
                % === Load block data === (CHANGE SO THAT THIS WILL SAVE PAST DATA)
                 [currentBlockStruct,referenceBlockStruct,...
                    behavioralData, imagingData, bitmapData, successFlag]=loadBlockData(datastruct, analysisBlockID, behavioralData, imagingData, bitmapData, blockID);
                 if ~successFlag
                     continue
                 end
                % === Select pdf save file name ===
                 pdfFilename=currentBlockStruct.psychneuroPDF;
        
                  % ===  Get orientation map ===
                  bitmapData=getColumnarBitmapV4(currentBlockStruct, imagingData, bitmapData, blockID, ...
                    pdfFilename, plotFlag, saveFlag);
        
                  % === Transform columnar positions to current cortical view ===
                  bitmapData=coregisterBitmap2GreenImgV2(currentBlockStruct,referenceBlockStruct, ...
                    imagingData,bitmapData,blockID, ...
                    pdfFilename, plotFlag,saveFlag);
        
                  % === Generate bitmap, correct for projector properties and camera-projector alignment===
                   [bitmapData]=convertForProjectorGPT(behavioralData, imagingData, bitmapData,...
                        currentBlockStruct,'proj2cam', blockID, ...
                        pdfFilename, plotFlag,saveFlagBMP,saveFlag); % or proj2cam
        
                    % === Plot behavioral biasing results ===
                    reportType='summary';
                    behavioralData=analyzeBlockPsychometrics(currentBlockStruct, behavioralData, bitmapData, blockID,...
                        pdfFilename, reportType, saveFlag);
                    toc
                    
                    % === Save data ===
                    imagingData.intg=[];behavioralData.TS=[];
                    CF; % close fig
            end
            save([mainPath '/Chip/Meta/summary/statistics' chamberWanted '.mat'],'bitmapData','behavioralData','imagingData','analysisBlockID','datastruct', '-v7.3')
            sprintf('Saved!')
           %% Load data if needed
           % load('Y:\Chip\Meta\summary\statisticsHPC.mat')
           % load('Y:\Chip\Meta\summary\statisticsR.mat')
    
            %% Psychometrics
            % === Session summary ===
            %behavioralData=analyzeSessionPsychometrics(behavioralData, bitmapData, datastruct, analysisBlockID, chamberWanted, saveFlag);
            %export_fig('Y:\users\PK\Eyal\meetings\summary\summarypsychometrics1to40.pdf','-pdf','-nocrop');
    
            % === Power series ===
            %nColumns=20;  
            %plotPowerSeries(bitmapData, behavioralData, nColumns);
            %export_fig('Y:/Chip/Meta/powerSeries/powerSeries.pdf','-pdf','-nocrop');
    
            % === Min columns ===
            %plotMinColumns(bitmapData, behavioralData, analysisBlockID)
            %export_fig('Y:/Chip/Meta/minColumnSeries/minColumns.pdf','-pdf','-nocrop');
    
            %% Modelling
            %open npmodelParameterEval.m
    end
end