%% Define analysis and current session to be analysed
% === Select analysis to run ===
% Main modes:
% 1. expt = generate bitmaps for experiments
% 2. summary = generate expt + psychometric + neurometric fits from experiment
% 3. psycluster = fit psychometric curves only, 
% 4. neurometrics = fit neurometric curves only
% 5. psyphidist/psyphiscatter = plot psycluster data as histogram/scatter (requires psycluster to be run first)
% 6. stability = quantifies stability of vasculature image across sessions
% 7. neurometricsfix = neurometrics for columnar optostim, fixation-state

analysisMode='psycluster';
monkeyName=''; % Pepper or blank
referenceSessID=1;%for ort map
currentSessID=2;%for biasing expt

% Save, plot flags
saveFlag=0;
saveFlagBMP=0;
plotFlag=1;
chambers={'R', 'L'};
%% Load dataStruct for the desired chamber
[mainPath, datastruct]=setupEnv(['users/PK/colStimPipeline/exptListBiasingFull' monkeyName '.m']);

for chamberID=1%1:2%:2
    nColumnsWanted=[]; chamberWanted=chambers{chamberID};
    analysisBlockID = organizeBlocks(datastruct, chamberWanted, nColumnsWanted);
    
    %% Run the desired analysis pipeline variant
    switch analysisMode
        case {'expt'}        
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
            dsReferenceSess=datastruct(referenceSessID); %moving reference ort map
            dsCurrentSess=datastruct(currentSessID); %fixed current session map, which ort map will be projected onto
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

        case 'summaryExpt'
            % Process only right chamber for Pepper
            for chamberID = 1  % Only process R chamber
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
                referenceBlockStruct = datastruct(referenceSessID);

                % Load single block data
                [currentBlockStruct, referenceBlockStruct,...
                    behavioralData, imagingData, bitmapData, successFlag] = loadBlockData(datastruct, analysisBlockID, behavioralData, imagingData, bitmapData, blockID, analysisMode);

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
                    imagingData, bitmapData, blockID,...
                    pdfFilename, plotFlag, saveFlag);

                % Generate bitmap, correct for projector properties and camera-projector alignment
                [bitmapData] = convertForProjectorGPT(behavioralData, imagingData, bitmapData,...
                    currentBlockStruct, 'cam2proj', blockID,...
                    pdfFilename, plotFlag, saveFlagBMP, saveFlag);

                % Clear temporary data
                imagingData.optoIntg = [];
                imagingData.baselineIntg = [];
                if isfield(imagingData, 'gaussfit')
                    imagingData.gaussfit(:,:,blockID) = [];
                end
                CF; % close figures
            end
            
        case {'summary'}
            if ~exist('behavioralData','var')
                behavioralData=[];
                imagingData=[];
                bitmapData=[];
            end
            
            for blockID=1:numel(analysisBlockID)%1:numel(analysisBlockID)
                disp(['=== Block ' num2str(blockID)  '/' num2str(numel(analysisBlockID)) '==='])
                tic
                if isfield(behavioralData,'auc') && size(behavioralData.auc,3)>=blockID
                     disp('(Skipping completed block)')
                     continue
                end
                % === Load block data ===
                 [currentBlockStruct,referenceBlockStruct,...
                    behavioralData, imagingData, bitmapData, successFlag]=loadBlockData(datastruct, analysisBlockID, behavioralData, imagingData, bitmapData, blockID, analysisMode);
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
                    
                    %% Neurometrics
                    trialOutcomeType='averageCorrect';
                    phiFull=analyzeProjection1D(currentBlockStruct, behavioralData, imagingData, bitmapData, trialOutcomeType, ...
                        blockID, gcf, 1, 'full', saveFlag); upFontSize(32,.01);
                    
                    % === Save data ===
                    imagingData.optoIntg=[];imagingData.baselineIntg=[]; imagingData.gaussfit(:,:,blockID)=[]; %behavioralData.optoTS(blockID)=[]; behavioralData.referenceTS(blockID)=[];
                    CF; % close fig
            end
            behavioralData=clearFields(behavioralData, {'gaborContrasts', 'percentageCorrect'});
            save([mainPath '/Chip/Meta/summary/statistics' chamberWanted '.mat'],'bitmapData','behavioralData','imagingData','analysisBlockID','datastruct')

        case {'psycluster'}
            %% Psychometrics: Cluster blocks by binning mean energy per block
            filterTag=true;
            for c=2%2:-1:1%numel(chambers):-1:1
                chamberWanted=chambers{c};
                % Load only if its not loaded
                if ~exist('dataTag')
                    load([mainPath 'Chip/Meta/summary/statistics' chamberWanted 'tag.mat']);
                elseif exist('dataTag')
                    if ~strcmp(dataTag,chamberWanted)
                        load([mainPath 'Chip/Meta/summary/statistics' chamberWanted 'tag.mat']);
                    end
                end
                if filterTag==true & c==1
                    analysisParams.cluster=[1]; % mean
                    analysisParams.columnMean=20; % mean
                    analysisParams.columnRange=4; % stdev    
                elseif filterTag==true & c==2
                    analysisParams.cluster=[3 4]; % mean
                    analysisParams.columnMean=20; % mean
                    analysisParams.columnRange=4; % stdev
                end
                analysisParams=[];
                [bins, binEdges, clusterIdx, validBlocks] = clusterEnergy(squeeze(bitmapData.energy), squeeze(bitmapData.nColumns), 'bin', 2, analysisParams);
                %savePDF(['psychometrics/' chamberWanted '-chamber/' num2str(numel(bins)) 'clusters'], 'Chip', 1, 1, 1)
                monkeyName=datastruct(analysisBlockID(1)).monkey;
                %save([mainPath '/Chip/Meta/psychometrics/fittingParams' chamberWanted '-' modelTypes{2} modelTypes{3}  '-C' num2str(cluster) '.mat'], 'mdlStruct', 'bitmapEnergy', 'AICCdeltax', 'AICCbeta');
                
                monkeyName='Chip';
                objFunc='MLE';
                modelTypes={'base','beta','deltax','bill','weibull','weibullfreeL','weibullfreeLS','weibullfreeAll','weibullbeta'};
                fieldName='AICc';
                constrainedParamStr='';
                saveFlag=0;plotFlag=1;
                mdlStruct=analyzePsychometricModels(monkeyName, chamberWanted, modelTypes, mainPath, ...
                    behavioralData, bitmapData, datastruct, analysisBlockID, clusterIdx, plotFlag, saveFlag);
                
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
                
                %% 20 column power x biasing
                figure
                betas = mdlStruct.([chamberWanted, 'beta' , 'C1']).fittedParams(:, 1);
                bitmapEnergy=mean(squeeze(bitmapData.energy(:,:,:))',2);

                blocksDesired = find(columns >= columnsDesired-columnSpread & columns <= columnsDesired+columnSpread & bitmapEnergy<3.5);
                xData=bitmapEnergy(blocksDesired);
                yData=betas(blocksDesired);
                scatter(xData,yData,300,'Marker', 'square', ...
                    'MarkerFaceColor', mkrColors(5,:), 'MarkerEdgeColor', 'k', 'LineWidth',2); hold on                
                fitSaturatingCurve(xData,yData,  mkrColors(5,:), gca); hold on
                yline(50,'LineStyle','--','color',[.5 .5 .5],'LineWidth',2,'HandleVisibility','off')
                title({'Effect of power on ?', sprintf('(%.0f ± %.0f columns)',columnsDesired,columnSpread)})
                switch c
                    case 2
                        xlim([0 .6])
                        ylim([0 100])
                        addSkippedTicks(0,.6,.05,'x')
                        addSkippedTicks(0,100,10,'y')
                    case 1
                        xlim([0 7])
                        ylim([0 100])
                        addSkippedTicks(0,7,.5,'x')
                        addSkippedTicks(0,100,10,'y')
                end
                xlabel('Total power (mW)')
                ylabel('?_{opto}')
                upFontSize(32,.015)
                export_fig(['Y:\users\PK\posters\figures\2024\betaPower' chamberWanted],'-svg','-png','-nocrop','-r600');
            end
                                    
            %% Neurometrics
        case {'neurometrics'}
            clc; neuroStruct=[];
            for chamberID=2%numel(chambers):-1:1
                nColumnsWanted=[]; chamberWanted=chambers{chamberID};
                
                if ~exist('dataTag')
                    load([mainPath 'Chip/Meta/summary/statistics' chamberWanted '.mat'], 'behavioralData', 'bitmapData', 'imagingData');
                elseif exist('dataTag')
                    if ~strcmp(dataTag,chamberWanted)
                        load([mainPath 'Chip/Meta/summary/statistics' chamberWanted '.mat'], 'behavioralData', 'bitmapData', 'imagingData');
                    end
                end
                                
                analysisBlockID = organizeBlocks(datastruct, chamberWanted, nColumnsWanted); %RESET

                %[bins, clusterIdx] = clusterEnergy(squeeze(bitmapData.energy), 'bin', 2);          
                [bins, binEdges, clusterIdx] = clusterEnergy(squeeze(bitmapData.energy), 'bin', 2);

                monkeyName='Chip';
                saveFlag=1;
                trialOutcomeType='averageColumn';  % 'average' or 'averageColumn'
                optostimMaskType='gaussian'; filterColumns=0; 
                ROIs={'full'};
                for roiID=1:numel(ROIs)
                    %neuroStruct=analyzeProjection1D(mainPath, datastruct, behavioralData, imagingData, bitmapData, analysisBlockID, trialOutcomeType, ...
                    %    gcf, 1, ROIs{roiID}, clusterIdx, filterColumns, optostimMaskType, saveFlag);
                    neuroStruct=analyzeProjection1DRefactored(mainPath, datastruct, behavioralData, imagingData, bitmapData, analysisBlockID, trialOutcomeType, ...
                        gcf, 1, ROIs{roiID}, clusterIdx, filterColumns, optostimMaskType, saveFlag);
                    %save([mainPath 'Chip/Meta/neurometric/neuroStruct' chamberWanted ROIs{roiID} optostimMaskType '.mat'], 'neuroStruct');
                end
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
                chamberIDs=2; saveFlag=0; filterColumns=20;
                plotPsyPhiCorrelation(datastruct, mainPath, chambers, chamberIDs, filterColumns, saveFlag)
            
        case {'psyphidist'}
                chamberIDs=1:2; saveFlag=0; filterColumns=20;
                plotPsyPhiDist(datastruct, mainPath, chambers, chamberIDs, filterColumns, saveFlag)

        case {'TS'}
            if ~exist('behavioralData','var')
                behavioralData=[];
                imagingData=[];
                bitmapData=[];
            end
            
            for blockID=1:numel(analysisBlockID)%1:numel(analysisBlockID)
                disp(['=== Block ' num2str(blockID)  '/' num2str(numel(analysisBlockID)) '==='])
                tic
                if isfield(behavioralData,'auc') && size(behavioralData.auc,3)>=blockID
                     disp('(Skipping completed block)')
                     continue
                end
                % === Load block data ===
                 [currentBlockStruct,referenceBlockStruct,...
                    behavioralData, imagingData, bitmapData, successFlag]=loadBlockData(datastruct, analysisBlockID, behavioralData, imagingData, bitmapData, blockID, analysisMode);

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
    end
end
