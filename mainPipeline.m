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
% 11. demo-optostim = 1x3 demonstration figure for a single session (PCA diff + column targeting)

%% Change these for experiment runs
analysisMode='demo-optostim';
monkeyName='Chip';%Pepper or Chip
currentSessID=89; %81;%for biasing expt

% Saving and plotting flags
saveFlag=1;
saveFlagBMP=0;
plotFlag=1;
skipImaging=1;

% Targets for analysisMode = 'metatable': {animalName, chamberLetter, srcFilename}
metatableTargets = { ...
    'Pepper', 'R', 'statisticsR-final38.mat';...
    'Chip',   'L', 'statisticsL-final43.mat'; ...
    'Chip',   'R', 'statisticsR-final16.mat'; ...
    };

%% Load dataStruct for the desired chamber
[mainPath, datastruct]=setupEnv(['users/PK/colStimPipeline/exptListBiasingFull' monkeyName '.m']);
chambers={'R', 'L'};
for chamberID=1:2
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
            
        case {'demo-optostim'}
            %% Single-session 1x3 demonstration figure
            % Initialize data structures
            blockData      = [];
            behavioralData = [];
            imagingData    = [];
            bitmapData     = [];
            blockID        = 1;
            analysisBlockID = currentSessID;

            % Grab session structs directly from datastruct
            currentBlockStruct   = datastruct(currentSessID);
            referenceBlockStruct = datastruct(currentBlockStruct.referenceBlockNo);

            % Load single-block data (7-output form matches summary/psycluster)
            [currentBlockStruct, referenceBlockStruct, ...
                blockData, behavioralData, imagingData, bitmapData, successFlag] = ...
                loadBlockData(datastruct, analysisBlockID, blockData, ...
                behavioralData, imagingData, bitmapData, blockID, analysisMode, skipImaging);

            if ~successFlag
                error('demo-optostim: failed to load block data for session %d', currentSessID);
            end

            pdfFilename = currentBlockStruct.psychneuroPDF;

            % --- Suppress intermediate figures created by component functions ---
            existingFigs = findall(0, 'Type', 'figure');

            origVisible = get(groot, 'DefaultFigureVisible');
            set(groot, 'DefaultFigureVisible', 'off');
            try
                % Run the same imaging + bitmap-processing chain as 'summary'
                [bitmapData, columnarProducts] = getColumnarBitmapV4( ...
                    currentBlockStruct, imagingData, bitmapData, blockID, ...
                    pdfFilename, 0, 0);

                bitmapData = coregisterBitmap2GreenImgV2( ...
                    currentBlockStruct, referenceBlockStruct, ...
                    imagingData, bitmapData, blockID, analysisMode, pdfFilename, 0, 0);

                [bitmapData, demoProducts] = convertForProjectorGPT2( ...
                    behavioralData, imagingData, bitmapData, ...
                    currentBlockStruct, 'proj2cam', blockID, pdfFilename, 0, 0, 0);

            catch ME
                set(groot, 'DefaultFigureVisible', origVisible);
                rethrow(ME);
            end

            % Close only the figures created during processing
            newFigs = setdiff(findall(0, 'Type', 'figure'), existingFigs);
            if ~isempty(newFigs)
                delete(newFigs);
            end

            % Restore visibility before drawing the demo figure
            set(groot, 'DefaultFigureVisible', origVisible);

            % --- Build and optionally save the 1x3 demo figure ---
            [figHandle, demoPlotData, outputFiles] = plotDemoOptostim( ...
                currentBlockStruct, imagingData, bitmapData, ...
                columnarProducts, demoProducts, blockID, currentSessID, ...
                mainPath, monkeyName, chamberWanted, saveFlag);

        case {'metatable'}
            % Generate metatables for all configured animal/chamber targets.
            % Independent of chamberWanted/monkeyName/clusterIdx/modelTypes.
            fprintf('=== metatable mode: %d targets ===\n', size(metatableTargets, 1));
            metatableSummary = struct([]);
            for tIdx = 1:size(metatableTargets, 1)
                tAnimal  = metatableTargets{tIdx, 1};
                tChamber = metatableTargets{tIdx, 2};
                tSource  = metatableTargets{tIdx, 3};
                fprintf('\n--- Target %d/%d: %s %s (%s) ---\n', ...
                    tIdx, size(metatableTargets, 1), tAnimal, tChamber, tSource);
                try
                    tResult = runMetatableForTarget(mainPath, tAnimal, tChamber, tSource);
                catch ME
                    fprintf('ERROR for %s %s: %s\n', tAnimal, tChamber, ME.message);
                    tResult = struct( ...
                        'animal',         tAnimal, ...
                        'chamber',        tChamber, ...
                        'sourcePath',     '', ...
                        'nBlocks',        0, ...
                        'matOutputPath',  '', ...
                        'xlsxOutputPath', '', ...
                        'success',        false);
                end
                if isempty(metatableSummary)
                    metatableSummary = tResult;
                else
                    metatableSummary(end+1) = tResult; %#ok<AGROW>
                end
            end
            fprintf('\n=== metatable batch summary ===\n');
            for si = 1:numel(metatableSummary)
                s = metatableSummary(si);
                fprintf('  %s %s | blocks=%d | success=%d | %s\n', ...
                    s.animal, s.chamber, s.nBlocks, s.success, s.matOutputPath);
            end

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
            nBehaviorBlocks=size(behavioralData.gaborContrasts,3);
            nBitmapBlocks=size(bitmapData.nColumns,2);
            if nBehaviorBlocks ~= nBitmapBlocks
                error(['psycluster block mismatch: behavioralData has %d blocks, ' ...
                    'but bitmapData.nColumns has %d blocks.'], ...
                    nBehaviorBlocks, nBitmapBlocks);
            end
            if numel(analysisBlockID) < nBehaviorBlocks
                error(['psycluster block mismatch: analysisBlockID has %d entries, ' ...
                    'but the loaded psychometric data have %d blocks.'], ...
                    numel(analysisBlockID), nBehaviorBlocks);
            end
            nBlocks=nBehaviorBlocks;
            clusterMethod='orderedPowerEffect';
            % Fit all eligible sessions once so per-session deltaBias is
            % available for the ordered power-band segmentation below.
            clusterIdx=ones(nBlocks,1);
            validBlocks=analysisBlockID(1:nBlocks);
            %save([mainPath '/' monkeyName '/Meta/psychometrics/powercluster_' method monkeyName chamberWanted '_' analysisParams.columnMean '??' analysisParams.columnRange '.mat'],'bitmapData','behavioralData','analysisBlockID','datastruct','dataTag')

            %savePDF(['psychometrics/' chamberWanted '-chamber/' num2str(numel(bins)) 'clusters'], 'Chip', 1, 1, 1)
            monkeyName=datastruct(analysisBlockID(1)).monkey;
            
            objFunc='MLE';
            modelTypes={'bill','weibullSignedBX0'}; %'weibullfreeall'
            fieldName='AICc';
            constrainedParamStr='';
            plotLine=1;
                       
            validBlocks=analysisBlockID(1:nBlocks);
            experimentLabels = makeExperimentLabels( ...
                datastruct, analysisBlockID(1:nBlocks));
            excludeFromPsycluster = false(nBlocks, 1);
            if strcmpi(monkeyName, 'Chip') && strcmpi(chamberWanted, 'L')
                excludeFromPsycluster = ...
                    strcmpi(string(experimentLabels(:)), "20230811R4");
            end
            if nnz(excludeFromPsycluster) > 1
                error('mainPipeline:DuplicateExcludedExperiment', ...
                    'Multiple Chip-L blocks matched 20230811R4.');
            end
            if strcmpi(monkeyName, 'Chip') && strcmpi(chamberWanted, 'L') && ...
                    nnz(excludeFromPsycluster) == 0
                warning('mainPipeline:ExcludedExperimentNotFound', ...
                    'Chip-L exclusion target 20230811R4 was not found.');
            end
            clusterIdx(excludeFromPsycluster) = 0;
            excludedAnalysisRows = find(excludeFromPsycluster);
            if any(excludeFromPsycluster)
                fprintf(['Excluding from psycluster | experiment=20230811R4 | ' ...
                    'analysis row=%d | monkey=Chip | chamber=L\n'], ...
                    excludedAnalysisRows(1));
            end

            supportedMonkey = any(strcmpi(monkeyName, {'Chip', 'Pepper'}));
            supportedChamber = any(strcmpi(chamberWanted, {'L', 'R'}));
            enableDeltaPermutationStats = supportedMonkey && supportedChamber;
            nDeltaPermutations = 5000;
            deltaPlotOpts = struct( ...
                'showDeltaPermutationStats', enableDeltaPermutationStats, ...
                'nDeltaPermutations', nDeltaPermutations, ...
                'deltaPermutationBaseSeed', 99173, ...
                'saveDeltaPermutationExamples', false);
            if enableDeltaPermutationStats
                fprintf(['Delta permutation statistics enabled | monkey=%s | ' ...
                    'chamber=%s | permutations=%d\n'], ...
                    monkeyName, chamberWanted, nDeltaPermutations);
            end
            fitPassPlotOpts = deltaPlotOpts;
            fitPassPlotOpts.showDeltaPermutationStats = false;
            fitPassPlotOpts.saveDeltaPermutationExamples = false;
            fitPassPlotOpts.figureVisible = 'off';
            fitPassPlotOpts.closeAfterRender = true;
            fprintf(['Pre-clustering fit pass: hidden figures, no PDF, ' ...
                'no permutation statistics.\n']);
            mdlStruct=analyzePsychometricModels(monkeyName, chamberWanted, modelTypes, mainPath, ...
                behavioralData, bitmapData, datastruct, analysisBlockID, clusterIdx, plotFlag, plotLine, false, fitPassPlotOpts);

            if plotFlag
                aggregateOpts = struct('binWidth', 5, ...
                    'useMedianForPlot', false, 'minSessionsPerBin', 1);
                clusterOpts = struct( ...
                    'minSessionsPerCluster', 3, ...
                    'effectWeight', 0.75, ...
                    'monotonicPenalty', 2, ...
                    'kCandidates', [2 3], ...
                    'minRelativeImprovementFor3Clusters', 0.10);
                [baselineModeByBlock, ~, baselineSourceByBlock] = ...
                    detectBaselineModeByBlock(datastruct, analysisBlockID, ...
                    1:nBlocks);
                powerAuditExperimentLabels = makeExperimentLabels( ...
                    datastruct, analysisBlockID(1:nBlocks));
                powerMetricsAudit = auditPowerMetrics(bitmapData, 1:nBlocks, ...
                    powerAuditExperimentLabels, [], baselineModeByBlock);
                totalPowerByBlock = powerMetricsAudit.Ptotal_recomputed;
                storedTotalPowerByBlock = mean( ...
                    squeeze(bitmapData.totalPowerToOnPixelsWithinROI_mW)', ...
                    2, 'omitnan');
                missingTotalPower = ~isfinite(totalPowerByBlock);
                totalPowerByBlock(missingTotalPower) = ...
                    storedTotalPowerByBlock(missingTotalPower);
                powerDensityByBlock = powerMetricsAudit.PDROI_recomputed;
                fprintf(['Power clustering source: canonical recomputed Ptotal ' ...
                    '= PDDMD * AreaON * tDC\n']);
                fprintf('Power clustering units: mW\n');
                fprintf('Per-session PDROI values (mW mm^{-2}): %s\n', ...
                    formatNumericVectorForLog(powerDensityByBlock));
                fprintf('Per-session Ptotal values (mW): %s\n', ...
                    formatNumericVectorForLog(totalPowerByBlock));
                aggregateModelType = modelTypes{2};
                sourceMdlField = [chamberWanted, aggregateModelType, 'C1'];

                if isfield(mdlStruct, sourceMdlField)
                    clusterMdl = mdlStruct.(sourceMdlField);
                    clusterBlocks = clusterMdl.clusterBlocksIdx(:);
                    if any(ismember(clusterBlocks, excludedAnalysisRows))
                        error('mainPipeline:ExcludedExperimentEnteredModel', ...
                            'Excluded experiment 20230811R4 entered clusterMdl.');
                    end
                    clusterMdl.baselineMode = baselineModeByBlock(clusterBlocks);
                    clusterMdl.baselineModeSourceValue = ...
                        baselineSourceByBlock(clusterBlocks);
                    clusterMdl.baselineModeSourceField = ...
                        'datastruct(analysisBlockID(blockIdx)).baselineTS';
                    mdlStruct.(sourceMdlField) = clusterMdl;
                    sessionPower = totalPowerByBlock(clusterBlocks);
                    sessionBias = clusterMdl.deltaBias(:);
                    sessionMask = clusterMdl.deltaMask(:);
                    [experimentalByBlock, controlByBlock] = ...
                        classifyOptostimColumnTargets(bitmapData.orts, nBlocks);
                    experimentalSessionMask = experimentalByBlock(clusterBlocks);
                    controlSessionMask = controlByBlock(clusterBlocks);
                    unknownSessionMask = ...
                        ~(experimentalSessionMask | controlSessionMask);
                    if any(unknownSessionMask)
                        warning('%d sessions have unsupported bitmapData.orts pairs and will be excluded from power-cluster aggregates.', ...
                            sum(unknownSessionMask));
                    end
                    if ~any(experimentalSessionMask)
                        error('No 0/90-degree optostim sessions are available to define power-cluster boundaries.');
                    end

                    experimentalRows = find(experimentalSessionMask);
                    [experimentalCluster, clusterSummary, clusterDiagnostics] = ...
                        clusterOrderedPowerEffect( ...
                        sessionPower(experimentalSessionMask), ...
                        sessionBias(experimentalSessionMask), ...
                        sessionMask(experimentalSessionMask), clusterOpts);
                    powerEffectCluster = assignOrderedPowerBands( ...
                        sessionPower, clusterDiagnostics.boundaries);
                    if any(powerEffectCluster(experimentalSessionMask) ~= ...
                            experimentalCluster)
                        error('Fixed-boundary assignment does not match the experimental clustering.');
                    end
                    powerEffectCluster(controlSessionMask) = ...
                        clusterDiagnostics.chosenK;
                    powerEffectCluster(unknownSessionMask) = NaN;
                    clusterDiagnostics.controlAssignmentMethod = ...
                        'highestPowerBand';
                    powerEffectClusterByBlock = nan(nBlocks, 1);
                    powerEffectClusterByBlock(clusterBlocks) = powerEffectCluster;
                    powerMetricsAudit.powerClusterAssignment = ...
                        powerEffectClusterByBlock;
                    groupMasks = {experimentalSessionMask, controlSessionMask};
                    groupCodes = {'experimental', 'control'};
                    groupLabels = {'0^{\circ}/90^{\circ}', ...
                        '45^{\circ}/135^{\circ}'};
                    aggregatePsychometrics = struct([]);
                    for groupIdx = 1:numel(groupMasks)
                        groupCluster = powerEffectCluster;
                        groupCluster(~groupMasks{groupIdx}) = NaN;
                        groupAggregate = ...
                            aggregatePsychometricCountsByPowerCluster( ...
                            clusterMdl, groupCluster, aggregateOpts);
                        if isempty(groupAggregate)
                            continue;
                        end

                        for aggregateClusterIdx = 1:numel(groupAggregate)
                            clusterID = groupAggregate(aggregateClusterIdx).clusterID;
                            localBlocks = find(groupCluster == clusterID);
                            sharedClusterBlocks = find(powerEffectCluster == clusterID);
                            sourceBlocks = clusterBlocks(localBlocks);
                            clusterPowers = sessionPower(localBlocks);
                            sharedClusterPowers = sessionPower(sharedClusterBlocks);

                            groupAggregate(aggregateClusterIdx).clusterMethod = ...
                                clusterMethod;
                            groupAggregate(aggregateClusterIdx).powerRange = ...
                                [min(sharedClusterPowers), max(sharedClusterPowers)];
                            groupAggregate(aggregateClusterIdx).powerMean = ...
                                mean(clusterPowers);
                            groupAggregate(aggregateClusterIdx).powerUnits = 'mW';
                            groupAggregate(aggregateClusterIdx).powerClusterCount = ...
                                clusterDiagnostics.chosenK;
                            groupAggregate(aggregateClusterIdx).columnTargetGroup = ...
                                groupCodes{groupIdx};
                            groupAggregate(aggregateClusterIdx).columnTargetLabel = ...
                                groupLabels{groupIdx};
                            groupAggregate(aggregateClusterIdx).targetSortOrder = groupIdx;
                            groupAggregate(aggregateClusterIdx).mdlRowIndices = ...
                                localBlocks(:)';
                            groupAggregate(aggregateClusterIdx).sourceBlockIndices = ...
                                sourceBlocks(:)';
                            groupAggregate(aggregateClusterIdx).baselineModes = ...
                                baselineModeByBlock(sourceBlocks);
                            groupAggregate(aggregateClusterIdx).baselineModeSummary = ...
                                summarizeBaselineModeLabels( ...
                                baselineModeByBlock(sourceBlocks));
                            groupAggregate(aggregateClusterIdx).stimulationStats = ...
                                summarizeClusterStimulationParameters( ...
                                bitmapData, sourceBlocks);
                        end

                        if isempty(aggregatePsychometrics)
                            aggregatePsychometrics = groupAggregate;
                        else
                            aggregatePsychometrics = ...
                                [aggregatePsychometrics; groupAggregate]; %#ok<AGROW>
                        end
                    end

                    for summaryIdx = 1:numel(clusterSummary)
                        localBlocks = experimentalRows( ...
                            clusterSummary(summaryIdx).sessionIndices);
                        clusterSummary(summaryIdx).mdlRowIndices = localBlocks(:)';
                        clusterSummary(summaryIdx).blockIndices = ...
                            clusterBlocks(localBlocks)';
                        clusterSummary(summaryIdx).datastructIndices = ...
                            analysisBlockID(clusterBlocks(localBlocks));
                    end
                    [clusterSummary, clusterMdl] = ...
                        attachPowerClusterModelFields( ...
                        clusterSummary, clusterMdl, sessionBias, ...
                        experimentalSessionMask, powerEffectCluster, ...
                        powerEffectClusterByBlock, analysisBlockID, ...
                        clusterBlocks);
                    mdlStruct.(sourceMdlField) = clusterMdl;
                    modelSpecificMdlStructPath = ...
                        savePsychometricModelStruct( ...
                        mainPath, monkeyName, chamberWanted, ...
                        aggregateModelType, mdlStruct);
                    fprintf('Power-cluster model-specific mdlStruct saved: %s\n', ...
                        modelSpecificMdlStructPath);

                    excludedAggregateContribution = countExcludedAggregateContribution(...
                        aggregatePsychometrics, excludedAnalysisRows);
                    if excludedAggregateContribution > 0
                        error('mainPipeline:ExcludedExperimentEnteredAggregates', ...
                            'Excluded experiment 20230811R4 entered aggregate structures.');
                    end
                    if strcmpi(monkeyName, 'Chip') && strcmpi(chamberWanted, 'L') && ...
                            ~isempty(excludedAnalysisRows)
                        fprintf(['Psycluster exclusion validated | excluded=20230811R4 | ' ...
                            'individual pages=%d | aggregate contribution=%d\n'], ...
                            numel(clusterBlocks), excludedAggregateContribution);
                    end

                    reportOutputDir = fullfile('Y:\users\PK\colStimPipeline', 'outputs');
                    if ~exist(reportOutputDir, 'dir')
                        mkdir(reportOutputDir);
                    end
                    animalID = animalIDForPsychometricReport(monkeyName);
                    reportFilename = fullfile(reportOutputDir, ...
                        sprintf('M%d-psychfit-%s-%s.pdf', ...
                        animalID, chamberWanted, aggregateModelType));
                    fprintf('Final psychometric PDF target: %s\n', reportFilename);
                    deltaPermutationPlotOpts = deltaPlotOpts;
                    reportState = [];
                    if saveFlag
                        reportState = initializeReportPDFAssembly( ...
                            reportFilename, monkeyName);
                        if size(clusterMdl.fittedParams, 1) ~= numel(clusterBlocks)
                            error('mainPipeline:ClusterMdlRowMismatch', ...
                                ['clusterMdl has %d fitted rows, but ' ...
                                'clusterBlocks has %d entries.'], ...
                                size(clusterMdl.fittedParams, 1), ...
                                numel(clusterBlocks));
                        end
                        clusterLabelsForBlocks = ...
                            prepareClusterLabelsForIndividualPages( ...
                            datastruct, analysisBlockID, clusterBlocks, ...
                            powerEffectClusterByBlock);
                        renderExperimentLabels = makeExperimentLabels(...
                            datastruct, analysisBlockID(clusterBlocks));
                        if any(strcmpi(string(renderExperimentLabels(:)), "20230811R4"))
                            error('mainPipeline:ExcludedExperimentRendered', ...
                                'Excluded experiment 20230811R4 reached individual-page rendering.');
                        end
                        fprintf('Post-clustering individual render: nRenderBlocks=%d | numel(clusterLabelsForBlocks)=%d | plotAverageFlag=%d\n', ...
                            numel(clusterBlocks), ...
                            numel(clusterLabelsForBlocks), false);
                        individualXFit = sort(nlinspace(0, 100, 100, ...
                            'linear'));
                        fprintf(['Post-clustering saved render: %d pages with real cluster ' ...
                            'labels and permutation statistics.\n'], numel(clusterBlocks));
                        [clusterMdl, reportState] = plotNakaRushtonFit5(behavioralData, bitmapData, ...
                            datastruct, analysisBlockID, clusterMdl, ...
                            clusterMdl.fittedParams(:, :, 1), ...
                            individualXFit, monkeyName, clusterBlocks, ...
                            0, plotLine, saveFlag, 1, ...
                            aggregateModelType, reportFilename, ...
                            clusterLabelsForBlocks, reportState, deltaPermutationPlotOpts);
                    end

                    powerBiasFigure = plotDeltaBiasByPower( ...
                        sessionPower, sessionBias);
                    if saveFlag
                        reportState = stageAndCloseReportFigure( ...
                            powerBiasFigure, reportState);
                    end

                    diagnosticOpts = struct('controlMask', controlSessionMask);
                    clusterDiagnosticFigure = plotOrderedPowerEffectClusters( ...
                        sessionPower, sessionBias, powerEffectCluster, ...
                        clusterDiagnostics, diagnosticOpts);
                    if saveFlag
                        reportState = stageAndCloseReportFigure( ...
                            clusterDiagnosticFigure, reportState);
                    end

                    experimentNumber = analysisBlockID(clusterBlocks);
                    chronologyOpts = struct('controlMask', controlSessionMask);
                    chronologyFigure = plotDeltaBiasChronology( ...
                        experimentNumber, sessionBias, powerEffectCluster, ...
                        chronologyOpts);
                    if saveFlag
                        reportState = stageAndCloseReportFigure( ...
                            chronologyFigure, reportState);
                    end

                    aggregatePlotOpts = deltaPermutationPlotOpts;
                    aggregatePlotOpts.modelTypeStr = aggregateModelType;
                    [aggregatePsychometricFits, aggregateFigures] = ...
                        plotAggregatedPowerClusterPsychometrics(...
                        aggregatePsychometrics, aggregatePlotOpts);
                    if saveFlag
                        reportState = stageAndCloseReportFigures( ...
                            aggregateFigures, reportState);
                    end

                    distributionOpts = struct( ...
                        'experimentIDsByBlock', analysisBlockID(1:nBlocks), ...
                        'baselineModeByBlock', baselineModeByBlock, ...
                        'showSourceAuditText', false);
                    [meanDistributionFigures, meanDistributionData, ...
                        meanDistributionAudit, meanDistributionStatsAudit] = ...
                        plotPowerClusterMeanDistributions( ...
                        clusterMdl, aggregatePsychometricFits, ...
                        distributionOpts);
                    if saveFlag
                        reportState = stageAndCloseReportFigures( ...
                            meanDistributionFigures, reportState);
                    end

                    parameterDistributionOpts = distributionOpts;
                    parameterDistributionOpts.modelFieldName = sourceMdlField;
                    [parameterFigures, parameterDistributionData, ...
                        parameterDistributionAudit, ...
                        parameterDistributionStatsAudit] = ...
                        plotPowerClusterFitParameters( ...
                        clusterMdl, aggregatePsychometricFits, ...
                        parameterDistributionOpts);
                    expectedParameterClusterIDs = sort([aggregatePsychometricFits.clusterID]);
                    generatedParameterClusterIDs = getFigurePowerClusterIDs(parameterFigures);
                    assert(isequal(generatedParameterClusterIDs(:)', ...
                        expectedParameterClusterIDs(:)'), ...
                        'Parameter figure cluster IDs do not match aggregate fits.');
                    fprintf('Parameter figures expected: %d | generated: %d | cluster IDs: [%s]\n', ...
                        numel(expectedParameterClusterIDs), numel(parameterFigures), ...
                        strjoin(cellstr(string(generatedParameterClusterIDs(:)')), ' '));
                    if saveFlag
                        stagedParameterClusterIDs = generatedParameterClusterIDs;
                        reportState = stageAndCloseReportFigures( ...
                            parameterFigures, reportState);
                        fprintf('Staged parameter cluster IDs: [%s] | pages=%d\n', ...
                            strjoin(cellstr(string(stagedParameterClusterIDs(:)')), ' '), ...
                            numel(stagedParameterClusterIDs));
                        assert(isequal(stagedParameterClusterIDs(:)', ...
                            expectedParameterClusterIDs(:)'), ...
                            'Staged parameter cluster IDs do not match aggregate fits.');
                        finalPdfPath = reportState.outputFilename;
                        fprintf('Saving final psychometric PDF to:\n%s\n', finalPdfPath);
                        finalPdfPath = finalizeReportPDFAssembly(reportState);
                        assert(isfile(finalPdfPath), ...
                            'Final psychometric PDF was not created: %s', finalPdfPath);
                        pdfInfo = dir(finalPdfPath);
                        assert(pdfInfo.bytes > 0, ...
                            'Final psychometric PDF is empty: %s', finalPdfPath);
                        fprintf('Saved final psychometric PDF: %s | %d bytes\n', ...
                            finalPdfPath, pdfInfo.bytes);
                    end

                    distributionSourceAudit = [ ...
                        meanDistributionAudit; parameterDistributionAudit];
                    distributionStatsAudit = [ ...
                        meanDistributionStatsAudit; ...
                        parameterDistributionStatsAudit];
                    auditBaseName = sprintf( ...
                        'distributionSourceAudit_%s_%s_%s', ...
                        regexprep(monkeyName, '[^A-Za-z0-9_-]', '_'), ...
                        regexprep(chamberWanted, '[^A-Za-z0-9_-]', '_'), ...
                        regexprep(aggregateModelType, '[^A-Za-z0-9_-]', '_'));
                    statsAuditBaseName = sprintf( ...
                        'distributionStatsAudit_%s_%s_%s', ...
                        regexprep(monkeyName, '[^A-Za-z0-9_-]', '_'), ...
                        regexprep(chamberWanted, '[^A-Za-z0-9_-]', '_'), ...
                        regexprep(aggregateModelType, '[^A-Za-z0-9_-]', '_'));
                    distributionAuditPaths = struct();
                    distributionStatsAuditPaths = struct();
                    powerAuditPaths = struct();
                    deltaPointAuditPaths = struct();
                    if false && saveFlag
                        auditOutputBase = fullfile(mainPath, monkeyName, ...
                            'Meta', 'summary', auditBaseName);
                        distributionAuditPaths = saveDistributionSourceAudit( ...
                            distributionSourceAudit, auditOutputBase);
                        statsAuditOutputBase = fullfile(mainPath, monkeyName, ...
                            'Meta', 'summary', statsAuditBaseName);
                        distributionStatsAuditPaths = ...
                            saveDistributionStatsAudit( ...
                            distributionStatsAudit, statsAuditOutputBase);
                        powerAuditBaseName = sprintf( ...
                            'powerMetricsAudit_%s_%s_%s', ...
                            regexprep(monkeyName, '[^A-Za-z0-9_-]', '_'), ...
                            regexprep(chamberWanted, '[^A-Za-z0-9_-]', '_'), ...
                            regexprep(aggregateModelType, '[^A-Za-z0-9_-]', '_'));
                        powerAuditOutputBase = fullfile(mainPath, monkeyName, ...
                            'Meta', 'summary', powerAuditBaseName);
                        powerAuditPaths = savePowerMetricsAudit( ...
                            powerMetricsAudit, powerAuditOutputBase);
                        if isfield(clusterMdl, 'deltaPointAudit') && ...
                                istable(clusterMdl.deltaPointAudit)
                            deltaAuditBaseName = sprintf( ...
                                'deltaPointAudit_%s_%s_%s', ...
                                regexprep(monkeyName, '[^A-Za-z0-9_-]', '_'), ...
                                regexprep(chamberWanted, '[^A-Za-z0-9_-]', '_'), ...
                                regexprep(aggregateModelType, '[^A-Za-z0-9_-]', '_'));
                            deltaAuditOutputBase = fullfile(mainPath, monkeyName, ...
                                'Meta', 'summary', deltaAuditBaseName);
                            deltaPointAuditPaths = saveDeltaPointAudit( ...
                                clusterMdl.deltaPointAudit, deltaAuditOutputBase);
                        end
                    end
                    aggregateField = [chamberWanted, aggregateModelType, ...
                        'PowerClusterAggregate'];
                    mdlStruct.(aggregateField) = aggregatePsychometrics;
                    mdlStruct.([aggregateField, 'Fit']) = aggregatePsychometricFits;
                    mdlStruct.([aggregateField, 'MeanDistributions']) = ...
                        meanDistributionData;
                    mdlStruct.([aggregateField, 'ParameterDistributions']) = ...
                        parameterDistributionData;
                    mdlStruct.([aggregateField, 'SourceAudit']) = ...
                        distributionSourceAudit;
                    mdlStruct.([aggregateField, 'StatsAudit']) = ...
                        distributionStatsAudit;
                    mdlStruct.([aggregateField, 'PowerMetricsAudit']) = ...
                        powerMetricsAudit;
                    if isfield(clusterMdl, 'deltaPointAudit')
                        mdlStruct.([aggregateField, 'DeltaPointAudit']) = ...
                            clusterMdl.deltaPointAudit;
                    end
                    if saveFlag
                        mdlStruct.([aggregateField, 'SourceAuditPaths']) = ...
                            distributionAuditPaths;
                        mdlStruct.([aggregateField, 'StatsAuditPaths']) = ...
                            distributionStatsAuditPaths;
                        mdlStruct.([aggregateField, 'PowerAuditPaths']) = ...
                            powerAuditPaths;
                        mdlStruct.([aggregateField, 'DeltaPointAuditPaths']) = ...
                            deltaPointAuditPaths;
                    end
                    mdlStruct.([aggregateField, 'Summary']) = clusterSummary;
                    mdlStruct.([aggregateField, 'Diagnostics']) = clusterDiagnostics;
                    mdlStruct.([aggregateField, 'Labels']) = powerEffectCluster;
                    mdlStruct.([aggregateField, 'ExperimentalMask']) = ...
                        experimentalSessionMask;
                    mdlStruct.([aggregateField, 'ControlMask']) = ...
                        controlSessionMask;
                    mdlStruct.([aggregateField, 'LabelsByBlock']) = ...
                        powerEffectClusterByBlock;
                end
            end

            save([mainPath '/' monkeyName '/Meta/summary/statistics' chamberWanted '-final' nBlockStr '.mat'],'blockData','bitmapData','behavioralData','analysisBlockID','datastruct','dataTag','mdlStruct')
            %{
            %% 20 column power x biasing
            figure
            conY = mdlStruct.([chamberWanted, 'weibullfreeAll' , 'C1']).yConOpto;
            inconY = mdlStruct.([chamberWanted, 'weibullfreeAll' , 'C1']).yInconOpto;
            deltaY = nanmean(conY-inconY,2);
            bitmapEnergy=mean(squeeze(bitmapData.totalPowerToOnPixelsWithinROI_mW(:,:,:))',2);
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
            
            title({'Power x ??con-incon (M2-R)', sprintf('%.0f ?? %.0f columns',columnsDesired,round(std(columns)))})
            axis square
            xlim([0 5])
            ylim([-10 20])
            addSkippedTicks(0,5,.5,'x')
            addSkippedTicks(-10,20,5,'y')
            %[p,h,stats] = ranksum(actualY(actualX>=3.5), [controlY(controlX>=3.5)' 2 -1 -3], 'tail','both')
            xlabel('Total energy (mW)')
            ylabel('?? correct %')
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
                bitmapEnergyCluster=mean(squeeze(bitmapData.totalPowerToOnPixelsWithinROI_mW(:,:,clusterBlocks))',2);
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
            bitmapEnergy=mean(squeeze(bitmapData.totalPowerToOnPixelsWithinROI_mW(:,:,:))',2);
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
            title({'Effect of biasing across sessions', sprintf('(%.0f ?? %.0f columns)',columnsDesired,columnSpread)})
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
                                    
        %% Metatable
        case {'metatable'}
            load([mainPath '/' monkeyName '/Meta/summary/statistics' chamberWanted '-final' nBlockStr '.mat'],'blockData','bitmapData','behavioralData','analysisBlockID','datastruct','dataTag','mdlStruct')
            columnsDesired=20; columnSpread=4;
            MetaTable = buildBlockMetadata(behavioralData, bitmapData, columnsDesired, columnSpread, ...
                blockData, datastruct, analysisBlockID, mdlStruct);
        %% Neurometrics
        case {'neurometrics'}
            clc; neuroStruct=[];
            nColumnsWanted=[]; chamberWanted=chambers{chamberID};
            
            % Build the full path once (more robust than manual concatenation)
            fileName = fullfile(mainPath, monkeyName, 'Meta', 'summary', ...
                                ['statistics' chamberWanted 'tag.mat']);
            
            % Load only when needed:
            %   ??? `dataTag` is missing,  OR
            %   ??? `dataTag` refers to a different chamber
            if ( ~exist('dataTag','var') || ~strcmp(dataTag, chamberWanted) ) && isfile(fileName)
                load(fileName, 'behavioralData', 'bitmapData', 'imagingData', 'dataTag');
            end
            
            analysisBlockID = organizeBlocks(datastruct, chamberWanted, nColumnsWanted); %RESET

            %[bins, clusterIdx] = clusterEnergy(squeeze(bitmapData.totalPowerToOnPixelsWithinROI_mW), 'bin', 2);
            analysisParams=[];
            [bins, binEdges, clusterIdx] =  clusterEnergy(squeeze(bitmapData.totalPowerToOnPixelsWithinROI_mW), squeeze(bitmapData.nColumns), 'bin', 5, analysisParams);

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
                plotPsyPhiDist(datastruct, mainPath, monkeyName, chambers, chamberIDs, filterColumns, saveFlag)

        case {'psydeltahist'}
                plotMultiChamberDeltaBiasHistogram();

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

function label = formatNumericVectorForLog(values)
    values = values(:)';
    if isempty(values)
        label = '[]';
        return;
    end
    finiteMask = isfinite(values);
    parts = strings(1, numel(values));
    for idx = 1:numel(values)
        if finiteMask(idx)
            parts(idx) = sprintf('%0.6g', values(idx));
        else
            parts(idx) = "NaN";
        end
    end
    label = ['[', char(strjoin(parts, ', ')), ']'];
end

function labels = makeExperimentLabels(datastruct, blockIndices)
    blockIndices = blockIndices(:);
    labels = strings(size(blockIndices));
    for idx = 1:numel(blockIndices)
        dataIdx = blockIndices(idx);
        if ~isfinite(dataIdx) || dataIdx < 1 || dataIdx > numel(datastruct)
            labels(idx) = "";
            continue;
        end
        labels(idx) = string(datastruct(dataIdx).date) + "R" + ...
            string(datastruct(dataIdx).run);
    end
end


function nMatches = countExcludedAggregateContribution(aggregatePsychometrics, excludedAnalysisRows)
    nMatches = 0;
    if isempty(excludedAnalysisRows) || isempty(aggregatePsychometrics)
        return;
    end
    for aggregateIdx = 1:numel(aggregatePsychometrics)
        if isfield(aggregatePsychometrics(aggregateIdx), 'sourceBlockIndices')
            nMatches = nMatches + sum(ismember(...
                aggregatePsychometrics(aggregateIdx).sourceBlockIndices(:), ...
                excludedAnalysisRows(:)));
        end
    end
end
function clusterLabels = prepareClusterLabelsForIndividualPages( ...
        datastruct, analysisBlockID, clusterBlocks, ...
        powerEffectClusterByBlock)
    validClusterIDs = powerEffectClusterByBlock( ...
        isfinite(powerEffectClusterByBlock) & powerEffectClusterByBlock > 0);
    if isempty(validClusterIDs)
        error('mainPipeline:NoPowerClusterAssignments', ...
            'No finalized finite positive power-cluster assignments exist.');
    end

    nClusters = max(validClusterIDs);
    expectedCount = zeros(1, nClusters);
    renderedCount = zeros(1, nClusters);
    for clusterID = 1:nClusters
        expectedCount(clusterID) = sum( ...
            powerEffectClusterByBlock(clusterBlocks) == clusterID);
    end

    clusterLabels = cell(numel(clusterBlocks), 1);
    renderClusterIDs = nan(numel(clusterBlocks), 1);
    fprintf('Cluster mapping ready | render pages %d | assignments %d | clusters %d\n', ...
        numel(clusterBlocks), numel(clusterBlocks), nClusters);

    for rowIdx = 1:numel(clusterBlocks)
        globalBlockIdx = clusterBlocks(rowIdx);
        experimentID = experimentLabelFromBlock( ...
            datastruct, analysisBlockID, globalBlockIdx);
        clusterID = NaN;
        if globalBlockIdx >= 1 && ...
                globalBlockIdx <= numel(powerEffectClusterByBlock)
            clusterID = powerEffectClusterByBlock(globalBlockIdx);
        end

        if ~isfinite(clusterID) || clusterID < 1
            error('mainPipeline:IncludedBlockMissingCluster', ...
                ['Included block %d has no finalized power-cluster ' ...
                'assignment.'], globalBlockIdx);
        end

        clusterLabels{rowIdx} = sprintf( ...
            'Cluster: %d/%d', round(clusterID), nClusters);
        renderClusterIDs(rowIdx) = round(clusterID);
        renderedCount(round(clusterID)) = ...
            renderedCount(round(clusterID)) + 1;
        if rowIdx <= 3
            fprintf('render page %d | experiment %s | cluster %d/%d | source row %d\n', ...
                rowIdx, experimentID, round(clusterID), nClusters, globalBlockIdx);
        end
    end

    for clusterID = 1:nClusters
        prepared = sum(renderClusterIDs == clusterID);
        if prepared ~= expectedCount(clusterID) || ...
                renderedCount(clusterID) ~= expectedCount(clusterID)
            error('mainPipeline:ClusterLabelCountMismatch', ...
                ['Cluster %d individual-page labels (%d) do not match ' ...
                'finalized assignment count (%d).'], clusterID, ...
                prepared, expectedCount(clusterID));
        end
    end
    fprintf('Prepared %d post-clustering individual-page cluster labels.\n', ...
        numel(clusterLabels));
end

function label = experimentLabelFromBlock(datastruct, analysisBlockID, blockIdx)
    if blockIdx < 1 || blockIdx > numel(analysisBlockID)
        label = sprintf('block%d', blockIdx);
        return;
    end

    dataIdx = analysisBlockID(blockIdx);
    if dataIdx < 1 || dataIdx > numel(datastruct)
        label = sprintf('block%d', blockIdx);
        return;
    end

    label = char(string(datastruct(dataIdx).date) + "R" + ...
        string(datastruct(dataIdx).run));
end

function label = formatClusterValueForAudit(value)
    if isfinite(value)
        label = sprintf('%d', round(value));
    else
        label = 'NaN';
    end
end


function animalID = animalIDForPsychometricReport(monkeyName)
    switch lower(char(string(monkeyName)))
        case 'chip'
            animalID = 28;
        case 'pepper'
            animalID = 32;
        otherwise
            error('Unsupported monkey name for psychometric report filename: %s', ...
                char(string(monkeyName)));
    end
end
function [clusterSummary, clusterMdl] = attachPowerClusterModelFields( ...
        clusterSummary, clusterMdl, sessionBias, experimentalSessionMask, ...
        powerEffectCluster, powerEffectClusterByBlock, analysisBlockID, ...
        clusterBlocks)
    positiveEffectTest = ...
        'Wilcoxon signed-rank versus zero, right-tailed, raw p<0.05 and mean>0';
    powerEffectCluster = powerEffectCluster(:);
    powerEffectClusterByBlock = powerEffectClusterByBlock(:);
    sessionBias = sessionBias(:);
    experimentalSessionMask = experimentalSessionMask(:);
    clusterBlocks = clusterBlocks(:);

    nModelRows = size(clusterMdl.fittedParams, 1);
    assert(numel(powerEffectCluster) == nModelRows, ...
        'mainPipeline:PowerClusterRowMismatch', ...
        'powerEffectCluster must align to clusterMdl.fittedParams rows.');
    assert(numel(clusterBlocks) == nModelRows, ...
        'mainPipeline:ClusterBlockRowMismatch', ...
        'clusterBlocks must align to clusterMdl.fittedParams rows.');
    assert(numel(sessionBias) == nModelRows, ...
        'mainPipeline:SessionBiasRowMismatch', ...
        'sessionBias must align to clusterMdl.fittedParams rows.');
    assert(numel(experimentalSessionMask) == nModelRows, ...
        'mainPipeline:ExperimentalMaskRowMismatch', ...
        'experimentalSessionMask must align to clusterMdl.fittedParams rows.');

    finiteCluster = powerEffectCluster(isfinite(powerEffectCluster));
    assert(all(finiteCluster > 0 & finiteCluster == round(finiteCluster)), ...
        'mainPipeline:InvalidPowerClusterID', ...
        'Finite powerCluster assignments must be positive integers.');
    finiteBlockCluster = powerEffectClusterByBlock(isfinite(powerEffectClusterByBlock));
    assert(all(finiteBlockCluster > 0 & ...
        finiteBlockCluster == round(finiteBlockCluster)), ...
        'mainPipeline:InvalidPowerClusterByBlockID', ...
        'Finite powerClusterByBlock assignments must be positive integers.');

    for summaryIdx = 1:numel(clusterSummary)
        clusterID = clusterSummary(summaryIdx).clusterID;
        clusterValues = sessionBias( ...
            experimentalSessionMask & powerEffectCluster == clusterID);
        [clusterP, clusterSignificant, reason] = ...
            classifyPositivePowerClusterEffect(clusterValues);
        clusterSummary(summaryIdx).deltaBiasOneSidedP = clusterP;
        clusterSummary(summaryIdx).positiveEffectSignificant = ...
            double(clusterSignificant);
        clusterSummary(summaryIdx).positiveEffectTest = positiveEffectTest;
        clusterSummary(summaryIdx).positiveEffectReason = reason;
    end

    powerCluster = powerEffectCluster(:)';
    powerClusterSig = nan(size(powerCluster));
    powerClusterP = nan(size(powerCluster));
    powerClusterMeanDeltaBias = nan(size(powerCluster));
    powerClusterByBlock = powerEffectClusterByBlock(:)';
    powerClusterByBlockSig = nan(size(powerClusterByBlock));
    powerClusterByBlockP = nan(size(powerClusterByBlock));
    powerClusterByBlockMeanDeltaBias = nan(size(powerClusterByBlock));

    for summaryIdx = 1:numel(clusterSummary)
        clusterID = clusterSummary(summaryIdx).clusterID;
        rowMask = powerCluster == clusterID;
        blockMask = powerClusterByBlock == clusterID;
        clusterSig = clusterSummary(summaryIdx).positiveEffectSignificant;
        clusterP = clusterSummary(summaryIdx).deltaBiasOneSidedP;
        clusterMean = clusterSummary(summaryIdx).deltaBiasMean;
        powerClusterSig(rowMask) = clusterSig;
        powerClusterP(rowMask) = clusterP;
        powerClusterMeanDeltaBias(rowMask) = clusterMean;
        powerClusterByBlockSig(blockMask) = clusterSig;
        powerClusterByBlockP(blockMask) = clusterP;
        powerClusterByBlockMeanDeltaBias(blockMask) = clusterMean;
    end

    clusterMdl.powerCluster = powerCluster;
    clusterMdl.powerClusterSig = powerClusterSig;
    clusterMdl.powerClusterP = powerClusterP;
    clusterMdl.powerClusterMeanDeltaBias = powerClusterMeanDeltaBias;
    clusterMdl.powerClusterExperimentID = analysisBlockID(clusterBlocks)';
    clusterMdl.powerClusterTest = positiveEffectTest;
    clusterMdl.powerClusterByBlock = powerClusterByBlock;
    clusterMdl.powerClusterByBlockSig = powerClusterByBlockSig;
    clusterMdl.powerClusterByBlockP = powerClusterByBlockP;
    clusterMdl.powerClusterByBlockMeanDeltaBias = ...
        powerClusterByBlockMeanDeltaBias;

    assert(numel(clusterMdl.powerCluster) == nModelRows, ...
        'mainPipeline:PowerClusterLengthMismatch', ...
        'powerCluster must align to fittedParams rows.');
    assert(numel(clusterMdl.powerClusterSig) == nModelRows, ...
        'mainPipeline:PowerClusterSigLengthMismatch', ...
        'powerClusterSig must align to fittedParams rows.');
    assert(numel(clusterMdl.powerClusterExperimentID) == nModelRows, ...
        'mainPipeline:PowerClusterExperimentIDLengthMismatch', ...
        'powerClusterExperimentID must align to fittedParams rows.');

    finiteSig = clusterMdl.powerClusterSig(isfinite(clusterMdl.powerClusterSig));
    assert(all(finiteSig == 0 | finiteSig == 1), ...
        'mainPipeline:InvalidPowerClusterSig', ...
        'Finite powerClusterSig values must be exactly 0 or 1.');
    unassignedRows = ~isfinite(clusterMdl.powerCluster);
    assert(all(isnan(clusterMdl.powerClusterSig(unassignedRows))) && ...
        all(isnan(clusterMdl.powerClusterP(unassignedRows))) && ...
        all(isnan(clusterMdl.powerClusterMeanDeltaBias(unassignedRows))), ...
        'mainPipeline:UnassignedPowerClusterMetadata', ...
        'Unassigned model rows must not receive cluster statistics.');

    for rowIdx = 1:nModelRows
        clusterID = clusterMdl.powerCluster(rowIdx);
        if ~isfinite(clusterID)
            continue;
        end
        summaryIdx = find([clusterSummary.clusterID] == clusterID, 1);
        assert(~isempty(summaryIdx), ...
            'mainPipeline:PowerClusterMissingSummary', ...
            'Power-cluster row has no matching clusterSummary row.');
        assert(isequaln(clusterMdl.powerClusterP(rowIdx), ...
            clusterSummary(summaryIdx).deltaBiasOneSidedP), ...
            'mainPipeline:PowerClusterPMismatch', ...
            'powerClusterP does not match clusterSummary.');
        assert(isequaln(clusterMdl.powerClusterMeanDeltaBias(rowIdx), ...
            clusterSummary(summaryIdx).deltaBiasMean), ...
            'mainPipeline:PowerClusterMeanMismatch', ...
            'powerClusterMeanDeltaBias does not match clusterSummary.');
        assert(isequaln(clusterMdl.powerClusterSig(rowIdx), ...
            clusterSummary(summaryIdx).positiveEffectSignificant), ...
            'mainPipeline:PowerClusterSigMismatch', ...
            'powerClusterSig does not match clusterSummary.');
    end

    fprintf('Power-cluster model fields saved | experiments=%d | clusters=%d\n', ...
        nModelRows, numel(clusterSummary));
    for summaryIdx = 1:numel(clusterSummary)
        fprintf('C%d: n=%d | mean DeltaBias=%0.6g | p(right)=%0.6g | significant=%d\n', ...
            clusterSummary(summaryIdx).clusterID, ...
            clusterSummary(summaryIdx).nSessions, ...
            clusterSummary(summaryIdx).deltaBiasMean, ...
            clusterSummary(summaryIdx).deltaBiasOneSidedP, ...
            clusterSummary(summaryIdx).positiveEffectSignificant);
    end
end

function [pValue, significant, reason, clusterMean] = ...
        classifyPositivePowerClusterEffect(clusterValues)
    finiteValues = clusterValues(isfinite(clusterValues));
    clusterMean = mean(clusterValues, 'omitnan');
    pValue = NaN;
    significant = false;
    if numel(finiteValues) < 2
        reason = 'insufficient finite values';
        return;
    end
    if numel(unique(finiteValues)) < 2
        reason = 'nonvarying finite values';
        return;
    end
    if exist('signrank', 'file') ~= 2
        reason = 'signrank unavailable';
        return;
    end
    try
        pValue = signrank(finiteValues, 0, 'tail', 'right');
    catch ME
        reason = ME.message;
        return;
    end
    if ~isfinite(pValue)
        reason = 'signrank returned nonfinite p-value';
        return;
    end
    significant = clusterMean > 0 && pValue < 0.05;
    if significant
        reason = 'significant positive effect';
    elseif clusterMean <= 0
        reason = 'mean deltaBias is not positive';
    else
        reason = 'raw one-sided p >= 0.05';
    end
end

function mdlStructPath = savePsychometricModelStruct( ...
        mainPath, monkeyName, chamberWanted, modelTypeStr, mdlStruct)
    columnsDesired = 20;
    mdlStructDir = fullfile(mainPath, monkeyName, 'Meta', ...
        'psychometrics', [chamberWanted '-chamber'], modelTypeStr);
    if ~exist(mdlStructDir, 'dir')
        mkdir(mdlStructDir);
    end
    mdlStructPath = fullfile(mdlStructDir, ...
        ['mdlStruct' chamberWanted num2str(columnsDesired) '.mat']);
    save(mdlStructPath, 'mdlStruct');
end

function saveDeltaBiasPermutationInspectionAudit(monkeyName, chamberWanted, modelType, clusterMdl, aggregateFits)
    outputDir = fullfile('Y:\users\PK\colStimPipeline', 'outputs', 'deltaBiasPermutation');
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    tag = sprintf('%s_%s_%s', monkeyName, chamberWanted, modelType);

    experimentContrasts = table();
    experimentSummary = table();
    clusterContrasts = table();
    clusterSummary = table();

    if isfield(clusterMdl, 'deltaBiasPermutationExperimentContrasts')
        experimentContrasts = clusterMdl.deltaBiasPermutationExperimentContrasts;
    end
    if isfield(clusterMdl, 'deltaBiasPermutationExperimentSummary')
        experimentSummary = clusterMdl.deltaBiasPermutationExperimentSummary;
    end
    for idx = 1:numel(aggregateFits)
        if isfield(aggregateFits(idx), 'merged') && ...
                isfield(aggregateFits(idx).merged, 'deltaBiasPermutationClusterContrasts')
            clusterContrasts = [clusterContrasts; ...
                aggregateFits(idx).merged.deltaBiasPermutationClusterContrasts]; %#ok<AGROW>
        end
        if isfield(aggregateFits(idx), 'merged') && ...
                isfield(aggregateFits(idx).merged, 'deltaBiasPermutationClusterSummary')
            clusterSummary = [clusterSummary; ...
                aggregateFits(idx).merged.deltaBiasPermutationClusterSummary]; %#ok<AGROW>
        end
    end

    matPath = fullfile(outputDir, sprintf('deltaBiasPermutationAudit_%s.mat', tag));
    expContrastPath = fullfile(outputDir, sprintf('deltaBiasPermutationExperimentContrasts_%s.csv', tag));
    expSummaryPath = fullfile(outputDir, sprintf('deltaBiasPermutationExperimentSummary_%s.csv', tag));
    clusterContrastPath = fullfile(outputDir, sprintf('deltaBiasPermutationClusterContrasts_%s.csv', tag));
    clusterSummaryPath = fullfile(outputDir, sprintf('deltaBiasPermutationClusterSummary_%s.csv', tag));

    save(matPath, 'experimentContrasts', 'experimentSummary', ...
        'clusterContrasts', 'clusterSummary');
    writetable(experimentContrasts, expContrastPath);
    writetable(experimentSummary, expSummaryPath);
    writetable(clusterContrasts, clusterContrastPath);
    writetable(clusterSummary, clusterSummaryPath);

    fprintf('Saved deltaBias permutation audit outputs:\n');
    fprintf('  %s\n', matPath);
    fprintf('  %s\n', expContrastPath);
    fprintf('  %s\n', expSummaryPath);
    fprintf('  %s\n', clusterContrastPath);
    fprintf('  %s\n', clusterSummaryPath);
end

function reportState = stageAndCloseReportFigures(figureHandles, reportState)
    for figureIdx = 1:numel(figureHandles)
        reportState = stageAndCloseReportFigure( ...
            figureHandles(figureIdx), reportState);
    end
end

function clusterIDs = getFigurePowerClusterIDs(figureHandles)
    clusterIDs = nan(size(figureHandles));
    for figureIdx = 1:numel(figureHandles)
        figHandle = figureHandles(figureIdx);
        if isempty(figHandle) || ~isgraphics(figHandle)
            error('mainPipeline:InvalidParameterFigureHandle', ...
                'Parameter figure %d is not a valid graphics handle.', figureIdx);
        end
        if ~isappdata(figHandle, 'PowerClusterID')
            error('mainPipeline:MissingParameterFigureClusterID', ...
                'Parameter figure %d is missing PowerClusterID appdata.', figureIdx);
        end
        clusterIDs(figureIdx) = getappdata(figHandle, 'PowerClusterID');
    end
end
function reportState = stageAndCloseReportFigure(figHandle, reportState)
    if isempty(figHandle) || ~isgraphics(figHandle)
        error(['Plot producer did not return a live graphics handle. ' ...
            'class=%s, size=%s'], class(figHandle), ...
            mat2str(size(figHandle)));
    end

    reportState = stageReportPDFPage(reportState, figHandle);
    close(figHandle);
end

