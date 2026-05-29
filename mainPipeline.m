    %% Define analysis and current session to be analysed
% === Select analysis to run ===
% Main modes:
% 1. expt = generate bitmaps for experiments
% 2. summary = generate expt + psychometric + neurometric fits from experiment
% 3. psycluster = fit psychometric curves only
% 4. neurometrics = fit neurometric curves only
% 5. psyphiscatter = plot psycluster data as histogram/scatter (requires psycluster to be run first)
% 6. stability = quantifies stability of vasculature image across sessions
% 7. neurometricsfix = neurometrics for columnar optostim, fixation-state
% 8. PRF = fit and plot PRF for single sessions
% 9. SIRF = fit and plot SIRF across sessions
% 10. psyphidist***

%% Change these for experiment runs
analysisMode='psycluster';%psyphidist
monkeyName='Pepper';%Pepper or Chip
currentSessID=81;%for biasing expt

% Saving and plotting flags
saveFlag=1;
saveFlagBMP=0;
plotFlag=1;
skipImaging=1;

%% Load dataStruct for the desired chamber
[mainPath, datastruct]=setupEnv(['users/PK/colStimPipeline/exptListBiasingFull' monkeyName '.m']);
chambers={'R', 'L'};
for chamberID=1
    nColumnsWanted=[]; chamberWanted=chambers{chamberID};
    analysisBlockID = organizeBlocks(datastruct, chamberWanted, nColumnsWanted);
    nBlockStr=num2str(numel(analysisBlockID));
    %% Run the desired analysis pipeline variant
    switch analysisMode
        case 'camproj'
            alignmentTransform = getAlignmentTransform(datastruct(currentSessID));
            %20250319O0, opto10%, optodicroic,
            %blackcard, 20hz
        case 'expt'
            % Setup for single chamber
            nColumnsWanted = [];
            chamberWanted = chambers{chamberID};

            % Initialize empty data structures
            behavioralData = [];
            imagingData = [];
            bitmapData = [];

            % Set current session as the only block to analyze
            blockID = 1;
            analysisBlockID = currentSessID;

            % Get current and reference session structs directly from datastruct
            currentBlockStruct = datastruct(currentSessID);
            referenceBlockStruct = datastruct(currentBlockStruct.referenceBlockNo);

            % Load single block data
            [currentBlockStruct, referenceBlockStruct,...
                behavioralData, imagingData, bitmapData, successFlag] = loadBlockData(datastruct, analysisBlockID, blockData,behavioralData, imagingData, bitmapData, blockID, analysisMode, skipImaging);

            if ~successFlag
                error('Failed to load block data');
            end

            % Set PDF filename for saving
            pdfFilename = currentBlockStruct.psychneuroPDF;

            % Get orientation map from reference session
            bitmapData = getColumnarBitmapV4(currentBlockStruct, imagingData, bitmapData, blockID,...
                pdfFilename, plotFlag, saveFlag);

            % Transform columnar positions to current cortical view
            bitmapData = coregisterBitmap2GreenImgV2(currentBlockStruct, referenceBlockStruct,...
                imagingData, bitmapData, blockID, analysisMode,...
                pdfFilename, plotFlag, saveFlag);

            % Generate bitmap, correct for projector properties and camera-projector alignment
            [bitmapData] = convertForProjectorGPT(behavioralData, imagingData, bitmapData,...
                currentBlockStruct, 'cam2proj', blockID,...
                pdfFilename, plotFlag, saveFlagBMP, saveFlag);
            
        case {'summary'}
            if ~exist('imagingData','var')
                blockData=[];
                behavioralData=[];
                imagingData=[];
                bitmapData=[];
            end
            
            for blockID=1:numel(analysisBlockID)
                disp(['=== Block ' num2str(blockID)  '/' nBlockStr ' (entry: ' num2str(analysisBlockID(blockID)) ')==='])
                tic
                if isfield(behavioralData,'auc') && size(behavioralData.auc,3)>=blockID
                     disp('(Skipping completed block)')
                     continue
                end
                % === Load block data ===
                 [currentBlockStruct,referenceBlockStruct,...
                    blockData, behavioralData, imagingData, bitmapData, successFlag]=loadBlockData(datastruct, analysisBlockID, blockData, behavioralData, imagingData, bitmapData, blockID, analysisMode, skipImaging);
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
                    analysisMode, pdfFilename, plotFlag,saveFlag);
        
                  % === Generate bitmap, correct for projector properties and camera-projector alignment===
                   [bitmapData]=convertForProjectorGPT2(behavioralData, imagingData, bitmapData,...
                        currentBlockStruct,'proj2cam', blockID, ...
                        pdfFilename, plotFlag,saveFlagBMP,saveFlag); % or proj2cam
        
                    % === Plot behavioral biasing results ===
                    reportType='summary';
                    behavioralData=analyzeBlockPsychometrics(currentBlockStruct, behavioralData, blockID,...
                        pdfFilename, reportType, saveFlag);
                    
                    %% Neurometrics
                    trialOutcomeType='averageCorrect';
                    %phiFull=analyzeProjection1D(currentBlockStruct, behavioralData, imagingData, bitmapData, trialOutcomeType, ...
                    %    blockID, gcf, 1, 'full', saveFlag); upFontSize(32,.01);
                    
                    % === Save data ===
                    %imagingData.optoIntg=[];imagingData.baselineIntg=[]; ...
                    %behavioralData.optoTS(blockID)=[]; behavioralData.baselineTS(blockID)=[]; behavioralData.referenceTS(blockID)=[]; % imagingData.gaussfit(:,:,blockID)=[];
                    CF; % close fig
                    behavioralData=clearFields(behavioralData, {'gaborContrasts', 'percentageCorrect','visualStim'});

            end
            behavioralData=clearFields(behavioralData, {'gaborContrasts', 'percentageCorrect','visualStim'});
            dataTag=chamberWanted;
            %save([mainPath '/' monkeyName '/Meta/summary/statistics' chamberWanted 'tag.mat'],'-v7.3','bitmapData','behavioralData','analysisBlockID','datastruct')
            save([mainPath '/' monkeyName '/Meta/summary/statistics' chamberWanted '-full' nBlockStr '.mat'],'blockData','bitmapData','behavioralData','imagingData','analysisBlockID','datastruct','dataTag')

        case {'psyclusterPre'}
            
            if ~exist('behavioralData','var')
                
                behavioralData=[];
                imagingData=[];
                bitmapData=[];
                load([mainPath '/' monkeyName '/Meta/summary/statistics' chamberWanted '-full' nBlockStr '.mat']);
            end
            
            for blockID=numel(analysisBlockID):-1:1
                disp(['=== Block ' num2str(blockID)  '/' nBlockStr '==='])
                tic
                if isfield(behavioralData,'auc') && size(behavioralData.auc,3)>=blockID
                     disp('(Skipping completed block)')
                     continue
                end
                % === Load block data ===
                 [currentBlockStruct,referenceBlockStruct,...
                    blockData, behavioralData, imagingData, bitmapData, successFlag]=loadBlockData(datastruct, analysisBlockID, blockData,behavioralData, imagingData, bitmapData, blockID, analysisMode, skipImaging);
                 if ~successFlag
                     continue
                 end
                pdfFilename=currentBlockStruct.psychneuroPDF;

                % === Plot behavioral biasing results ===
                reportType='summary';
                behavioralData=analyzeBlockPsychometrics(currentBlockStruct, behavioralData, blockID,...
                    pdfFilename, reportType, saveFlag);
                toc
                
                % === Save data ===
                %imagingData.optoIntg=[];imagingData.baselineIntg=[]; imagingData.gaussfit(:,:,blockID)=[]; behavioralData.optoTS(blockID)=[]; behavioralData.baselineTS(blockID)=[]; behavioralData.referenceTS(blockID)=[];
                CF; % close fig
            end
            behavioralData=clearFields(behavioralData, {'gaborContrasts', 'percentageCorrect','visualStim'});
            dataTag=chamberWanted;
            save([mainPath '/' monkeyName '/Meta/summary/statistics' chamberWanted '-psychometricsPre' nBlockStr '.mat'],'blockData','bitmapData','behavioralData','analysisBlockID','datastruct','dataTag')

        case {'psycluster'}
            %% Psychometrics: Cluster blocks by binning mean energy per block
            filterTag=true;
            % Load only if its not loaded
            if ~exist('dataTag')
                load([mainPath '/' monkeyName '/Meta/summary/statistics' chamberWanted '-psychometricsPre' nBlockStr '.mat']);
            elseif exist('dataTag')
                if ~strcmp(dataTag,chamberWanted)
                    load([mainPath '/' monkeyName '/Meta/summary/statistics' chamberWanted '-psychometricsPre' nBlockStr '.mat']);
                end
            end
            if filterTag==true
                analysisParams.columnMean=20; % mean
                analysisParams.columnRange=4; % stdev
            end
            analysisParams=[];
            nBlocks=numel(analysisBlockID);
            method='kmeans';
            [bins, binEdges, clusterIdx, validBlocks] = clusterEnergy(squeeze(bitmapData.meanPowerDensityWithinROI_mWmm2),...
                squeeze(bitmapData.nColumns), method, 1, analysisParams);
            %save([mainPath '/' monkeyName '/Meta/psychometrics/powercluster_' method monkeyName chamberWanted '_' analysisParams.columnMean '±' analysisParams.columnRange '.mat'],'bitmapData','behavioralData','analysisBlockID','datastruct','dataTag')

            %savePDF(['psychometrics/' chamberWanted '-chamber/' num2str(numel(bins)) 'clusters'], 'Chip', 1, 1, 1)
            monkeyName=datastruct(analysisBlockID(1)).monkey;
            
            objFunc='MLE';
            modelTypes={'bill','weibullfreeAll'};
            fieldName='AICc';
            constrainedParamStr='';
            plotLine=1;
                       
            validBlocks=analysisBlockID;
            mdlStruct=analyzePsychometricModels(monkeyName, chamberWanted, modelTypes, mainPath, ...
                behavioralData, bitmapData, datastruct, analysisBlockID, clusterIdx, plotFlag, plotLine, saveFlag);
            
            save([mainPath '/' monkeyName '/Meta/summary/statistics' chamberWanted '-final' nBlockStr '.mat'],'blockData','bitmapData','behavioralData','analysisBlockID','datastruct','dataTag','mdlStruct')

            %{
            %% 20 column power x biasing
            figure
            conY = mdlStruct.([chamberWanted, 'weibullfreeAll' , 'C1']).yConOpto;
            inconY = mdlStruct.([chamberWanted, 'weibullfreeAll' , 'C1']).yInconOpto;
            deltaY = nanmean(conY-inconY,2);
            bitmapEnergy=mean(squeeze(bitmapData.energy(:,:,:))',2);
            columns=mean(bitmapData.nColumns(:,:))';
            columnsDesired=20; columnSpread=2;
            blocksActual = find(columns >= columnsDesired-columnSpread &...
                columns <= columnsDesired+columnSpread &...
                mean(squeeze(bitmapData.orts)==[0;90])');
            
            blocksControl = find(columns >= columnsDesired-columnSpread &...
                columns <= columnsDesired+columnSpread &...
                mean(squeeze(bitmapData.orts)==[45;135])');
            
            yline(0,'LineStyle','--','color',[.5 .5 .5],'LineWidth',1.5,'HandleVisibility','off'); hold on
            % Plot actual data
            mdlScatter=fitSaturatingCurve(bitmapEnergy(blocksActual),deltaY(blocksActual),  'k', 1); hold on
            actualX=bitmapEnergy(blocksActual);
            actualY=deltaY(blocksActual);
            scatter(actualX,actualY,200,'Marker', 'square', ...
                'MarkerFaceColor', 'magenta', 'MarkerEdgeColor', 'k', 'LineWidth',2); hold on                
            % Plot control data
            
            controlX=bitmapEnergy(blocksControl);
            controlY=deltaY(blocksControl);
            if ~isempty(blocksControl)
                scatter(controlX+rand(1,5)'/10,controlY,200,'Marker', 'square', ...
                    'MarkerFaceColor', [1 1 1]*.7, 'MarkerEdgeColor', 'k', 'LineWidth',2); hold on
            end            

            save([mainPath '/' monkeyName '/Meta/psychometrics/psychfit' chamberWanted '-' modelTypes{2} '.mat'], 'mdlStruct', 'bitmapEnergy',...
                'behavioralData','analysisBlockID','datastruct','dataTag');
            
            title({'Power x Δcon-incon (M2-R)', sprintf('%.0f ± %.0f columns',columnsDesired,round(std(columns)))})
            axis square
            xlim([0 5])
            ylim([-10 20])
            addSkippedTicks(0,5,.5,'x')
            addSkippedTicks(-10,20,5,'y')
            %[p,h,stats] = ranksum(actualY(actualX>=3.5), [controlY(controlX>=3.5)' 2 -1 -3], 'tail','both')
            xlabel('Total energy (mW)')
            ylabel('Δ correct %')
            upFontSize(20,.02)
            export_fig(['Y:\users\PK\colStimPipeline\figures\powerxdeltaY-' monkeyName '-' chamberWanted],'-svg','-png','-nocrop','-r600');


            %% Per cluster: Minimum columns x beta
            nClusters=sort(unique(clusterIdx));
            mkrColors=slanCM('bold',20);mkrColors=mkrColors([2 6 4 8 10],:);
            figure('Name','Min. columns')
            yline(50,'LineStyle','--', 'LineWidth', 2.5); hold on
            for cluster=nClusters
                disp(['=== Cluster ' num2str(cluster) '===='])
                clusterBlocks=find(clusterIdx==cluster);
                nClusterBlocks=numel(clusterBlocks);
                % Get values
                nColumnsCluster=mean(bitmapData.nColumns(:,clusterBlocks)',2);
                bitmapEnergyCluster=mean(squeeze(bitmapData.energy(:,:,clusterBlocks))',2);
                betaCluster = mdlStruct.([chamberWanted, 'beta' , 'C' num2str(cluster)]).fittedParams(clusterBlocks, 1);
                % plot
                plot(nColumnsCluster,betaCluster, 'square', 'MarkerSize', 15, 'MarkerFaceColor', mkrColors(cluster,:), 'MarkerEdgeColor', 'k', 'linewidth', 2.5); hold on
                ylabel('\beta', 'FontName', 'Arial'); % Label the x-axis as 'AUC'
                xlabel('No. of columns (block average)'); % Label the y-axis as 'Mean nColumns'
                title('Effect of power and no. of columns stimulated','FontWeight','normal', 'Interpreter', 'none') % Title for the plot
                legend({'Chance','C1','C2','C3','C4','C5'},'Location','southwest', 'NumColumns', 2)
                xticks([0:5:40])
                yticks([0:10:100])
                xlim([0 40])
                ylim([30 100])
                upFontSize(24,.01)
            end
            hold off
            save([mainPath '/Chip/Meta/summary/statistics' chamberWanted 'full.mat'],'bitmapData','behavioralData','imagingData','analysisBlockID','datastruct','mdlStruct')
            %% Effect over sessions
            desiredPos=setFig([1000         558         560         420])

            columns=mean(bitmapData.nColumns',2);
            columnsDesired=20; columnSpread=4;
            bitmapEnergy=mean(squeeze(bitmapData.energy(:,:,:))',2);
            blocksDesired = find(columns >= columnsDesired-columnSpread & columns <= columnsDesired+columnSpread & bitmapEnergy<3.5);

            % Set xyz data
            xDataStr=vertcat(datastruct(analysisBlockID(blocksDesired)).date);
            xData=1:length(blocksDesired);
            yData=mdlStruct.RbetaC1.fittedParams(blocksDesired,1)';
            zData=bitmapEnergy(blocksDesired);
            
            % Normalize zData to the range of the colormap
            zMin = min(zData);
            zMax = max(zData);
            zNormalized = (zData - zMin) / (zMax - zMin);

            % Choose a colormap (e.g., 'jet', 'parula', etc.)
            colormap(viridis)
            colormapData = colormap(viridis);
            % Convert normalized zData to indices of the colormap
            numColors = size(colormapData, 1);
            colorIndices = round(1 + zNormalized * (numColors - 1)); 


            % Plot each point with the colormap for face color and black edge color
            scatter(xData, yData, 300, colormapData(colorIndices, :), 's', ...
                'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerFaceColor', 'flat');
            hold on;
            yline(50,'LineStyle','--','color',[.5 .5 .5],'LineWidth',2,'HandleVisibility','off')
            title({'Effect of biasing across sessions', sprintf('(%.0f ± %.0f columns)',columnsDesired,columnSpread)})
            xlim([0 16])
            ylim([0 100])
            addSkippedTicks(0,20,2,'x')
            addSkippedTicks(0,100,10,'y')
            xlabel('Session')
            ylabel('?_{opto}')
            axis square
            upFontSize(32,.015)
            % Add colorbar to the plot
            cb = colorbar(); cb.LineWidth=2; cb.FontSize=18
            cbL=ylabel(cb,'Power (mW)','FontSize',18,'Rotation',90);
            caxis([0 zMax]); % Set color axis range based on zData
            export_fig(['Y:\users\PK\posters\figures\2024\betaSessions' chamberWanted],'-svg','-png','-nocrop','-r600');
            %}
                                    
            %% Neurometrics
        case {'neurometrics'}
            clc; neuroStruct=[];
            nColumnsWanted=[]; chamberWanted=chambers{chamberID};
            
            % Build the full path once (more robust than manual concatenation)
            fileName = fullfile(mainPath, monkeyName, 'Meta', 'summary', ...
                                ['statistics' chamberWanted 'tag.mat']);
            
            % Load only when needed:
            %   – `dataTag` is missing,  OR
            %   – `dataTag` refers to a different chamber
            if ( ~exist('dataTag','var') || ~strcmp(dataTag, chamberWanted) ) && isfile(fileName)
                load(fileName, 'behavioralData', 'bitmapData', 'imagingData', 'dataTag');
            end
            
            analysisBlockID = organizeBlocks(datastruct, chamberWanted, nColumnsWanted); %RESET

            %[bins, clusterIdx] = clusterEnergy(squeeze(bitmapData.energy), 'bin', 2);          
            analysisParams=[];
            [bins, binEdges, clusterIdx] =  clusterEnergy(squeeze(bitmapData.energy), squeeze(bitmapData.nColumns), 'bin', 5, analysisParams);

            monkeyName='Chip';
            trialOutcomeType='average';  % 'average' or 'averageColumn'
            optostimMaskType='gaussian'; filterColumns=0; 
            ROIs={'full'};
            for roiID=1:numel(ROIs)
                %neuroStruct=analyzeProjection1D(mainPath, datastruct, behavioralData, imagingData, bitmapData, analysisBlockID, trialOutcomeType, ...
                %    gcf, 1, ROIs{roiID}, clusterIdx, filterColumns, optostimMaskType, saveFlag);
                neuroStruct=analyzeProjection1DRefactored(mainPath, datastruct, behavioralData, imagingData, bitmapData, analysisBlockID, trialOutcomeType, ...
                    gcf, 1, ROIs{roiID}, clusterIdx, filterColumns, optostimMaskType, saveFlag);
                %save([mainPath 'Chip/Meta/neurometric/neuroStruct' chamberWanted ROIs{roiID} optostimMaskType '.mat'], 'neuroStruct');
            end
        1;
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

        case {'psy'}
            chamberIDs=2; saveFlag=0; filterColumns=0;
            plotPsyRaw(datastruct, mainPath, chambers, chamberIDs, filterColumns, saveFlag)

        case {'psyphiscatter'}
                chamberIDs=1; saveFlag=0; filterColumns=20;
                plotPsyPhiCorrelation(datastruct, mainPath, chambers, chamberIDs, filterColumns, saveFlag)
            
        case {'psyphidist'}
                chamberIDs=chamberID; saveFlag=1; filterColumns=1;
                plotPsyPhiDist(datastruct, mainPath, monkeyName, chambers, chamberIDs, nBlockStr, filterColumns, saveFlag)

        case {'PRF'}
            fitPRFv2
        case {'SIRF'}
            fitSIRF
        case {'TS'}
            if ~exist('behavioralData','var')
                behavioralData=[];
                imagingData=[];
                bitmapData=[];
            end
            
            for blockID=1:numel(analysisBlockID)%1:numel(analysisBlockID)
                disp(['=== Block ' num2str(blockID)  '/' nBlockStr '==='])
                tic
                if isfield(behavioralData,'auc') && size(behavioralData.auc,3)>=blockID
                     disp('(Skipping completed block)')
                     continue
                end
                % === Load block data ===
                 [currentBlockStruct,referenceBlockStruct,...
                    behavioralData, imagingData, bitmapData, successFlag]=loadBlockData(datastruct, analysisBlockID, blockData,behavioralData, imagingData, bitmapData, blockID, analysisMode);

                % Save for table
                optoTrialsIdx=behavioralData.optoTS(blockID).Header.Conditions.TypeCond>0;
                param(blockID).ort=unique(behavioralData.optoTS(blockID).Header.Conditions.GaborOrt(optoTrialsIdx));
                param(blockID).sz=unique(behavioralData.optoTS(blockID).Header.Conditions.GaborSize(optoTrialsIdx));
                param(blockID).sf=unique(behavioralData.optoTS(blockID).Header.Conditions.GaborSF(optoTrialsIdx));
                param(blockID).phs=unique(behavioralData.optoTS(blockID).Header.Conditions.GaborPhase(optoTrialsIdx));
                param(blockID).contrast=unique(behavioralData.optoTS(blockID).Header.Conditions.StimCon(optoTrialsIdx));
                param(blockID).pos=unique(behavioralData.optoTS(blockID).Header.Conditions.StimPosCond(optoTrialsIdx));
            end
            1;
        case {'expt-deprecated'}        
            %% Fast pipeline (coregisters within session ort map to optostim block green image)
            % Get correction for camera-projector alignment
            % imregtform of orignal and recovered bitmap, apply to correct for
            % camera-projector alignment
            %open projectorCameraCalibration
            
            % Load projector-camera alignment transformation matrix
            alignmentSessionFolder='Y:\Chip\Chip20230815\' % test
            load([alignmentSessionFolder 'alignmentTransform.mat'],'alignmentTransform')
    
            % Define reference (orientation map) session and current session
            % (for optostim experiment)
            dsCurrentSess=datastruct(currentSessID); %fixed current session map, which ort map will be projected onto
            dsReferenceSess=datastruct(dsCurrentSess.referenceBlockNo); %moving reference ort map
            filenameStructCurrent=generateFilenames(dsCurrentSess);
            pdfFilename=[filenameStructCurrent.psychneuroPDF 'test']%filenameStructCurrent.psychneuroPDF;
    
            % Get reference session orientation map
            [columnarBitmap,VERpca,columnarmapStats]=getColumnarBitmapV4(mainPath,dsReferenceSess,dsCurrentSess,bitmapParams, ...
              plotFlag, saveFlag, pdfFilename);
    
            % Coregister reference session orientation map to the cortical view of current session 
            [columnarBitmapCoregistered, columnarPCAsCoregistered]=coregisterBitmap2GreenImgV2(dsReferenceSess,dsCurrentSess, ...
              columnarBitmap,VERpca,plotFlag,saveFlag, pdfFilename);
            
            % create bitmap from orientation map, applying alignment matrix and
            % accounting for projector dimensions (i.e. a rectangle)
            orts=[0 90]; %0:15:165;
            HE=1000;
            [projBitmapTRBB,bitmapsCamSpace,nBlobs,medianBlobAreas]=convertForProjector(dsReferenceSess,dsCurrentSess,columnarBitmapCoregistered,orts,...
                bitmapParams.gridSize,bitmapParams.gammaCorrFactor,bitmapParams.sensitivity,'cam2proj',alignmentTransform,plotFlag,saveFlagBMP,saveFlag, pdfFilename);
    end
end


