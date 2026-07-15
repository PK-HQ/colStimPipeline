function neurometricOptostim = runNeurometricOptostimAnalysis( ...
    mainPath, datastruct, currentSessID, monkeyName, chamberWanted, saveFlag, plotFlag, skipImaging)
%RUNNEUROMETRICOPTOSTIMANALYSIS Condition-average optostim neurometrics.
%
% This single-session analysis reuses the demo-optostim loading and
% coregistration chain, then quantifies nContrasts x 6 condition-average
% response maps with one shared broad Gaussian and continuous 0/90-degree
% column-domain templates.

if nargin < 8 || isempty(skipImaging)
    skipImaging = 0;
end
if nargin < 7 || isempty(plotFlag)
    plotFlag = 1;
end
if nargin < 6 || isempty(saveFlag)
    saveFlag = 0;
end

currentBlockStruct = datastruct(currentSessID);
if ~strcmp(chamberWanted, currentBlockStruct.chamber)
    error('runNeurometricOptostim:WrongChamber', ...
        'Session %d is chamber %s, not requested chamber %s.', ...
        currentSessID, currentBlockStruct.chamber, chamberWanted);
end

referenceBlockStruct = datastruct(currentBlockStruct.referenceBlockNo);
blockData = [];
behavioralData = [];
imagingData = [];
bitmapData = [];
blockID = 1;
analysisBlockID = currentSessID;
analysisMode = 'neurometric-optostim';

[currentBlockStruct, referenceBlockStruct, blockData, ...
    behavioralData, imagingData, bitmapData, successFlag] = ...
    loadBlockData(datastruct, analysisBlockID, blockData, behavioralData, ...
    imagingData, bitmapData, blockID, analysisMode, skipImaging);

if ~successFlag
    error('runNeurometricOptostim:LoadFailed', ...
        'Failed to load block data for session %d.', currentSessID);
end

filenameStruct = currentBlockStruct;
pdfFilename = currentBlockStruct.psychneuroPDF;

existingFigs = findall(0, 'Type', 'figure');
origVisible = get(groot, 'DefaultFigureVisible');
set(groot, 'DefaultFigureVisible', 'off');
try
    [bitmapData, columnarProducts] = getColumnarBitmapV4( ...
        currentBlockStruct, imagingData, bitmapData, blockID, ...
        pdfFilename, 0, 0);

    bitmapData = coregisterBitmap2GreenImgV2( ...
        currentBlockStruct, referenceBlockStruct, imagingData, bitmapData, ...
        blockID, analysisMode, pdfFilename, 0, 0);

    [bitmapData, demoProducts] = convertForProjectorGPT2( ...
        behavioralData, imagingData, bitmapData, currentBlockStruct, ...
        'proj2cam', blockID, pdfFilename, 0, 0, 0);
catch ME
    set(groot, 'DefaultFigureVisible', origVisible);
    rethrow(ME);
end
newFigs = setdiff(findall(0, 'Type', 'figure'), existingFigs);
if ~isempty(newFigs)
    delete(newFigs);
end
set(groot, 'DefaultFigureVisible', origVisible);

options = defaultNeurometricOptostimOptions();
options.plotFlag = plotFlag;
options.saveFlag = saveFlag;
repoRoot = fileparts(mfilename('fullpath'));
options.outputDir = fullfile(repoRoot, 'outputs', 'neurometric-optostim');
options.sessionLabel = sprintf('%s%sR%s', currentBlockStruct.monkey, ...
    currentBlockStruct.date, currentBlockStruct.run);

dataAudit = auditNeurometricOptostimData(currentBlockStruct, filenameStruct, ...
    behavioralData, imagingData, bitmapData, columnarProducts, demoProducts, blockID);

[responseMaps, conditionIndices, conditionDefinitions, stimContrastValues, ...
    responseSources, conditionMatrixAudit] = buildNeurometricConditionMatrix( ...
    behavioralData, filenameStruct, blockID);

[templates, analysisMask, templateAudit] = buildNeurometricTemplates( ...
    imagingData, bitmapData, columnarProducts, demoProducts, blockID, options);

[templates.G_shared, sharedFit, gaussianQC] = fitSharedGaussianFootprint( ...
    responseMaps, templates.G_reference, analysisMask, ...
    dataAudit.pixelsPerMM, options);
templates.G_shared = templates.G_shared ./ max(templates.G_shared(analysisMask));

[templates.W_COL0, templates.W_COL90, templateQC] = ...
    buildColumnReadoutWeights(templates, analysisMask, options);

[filteredMaps, metrics, ratioQC] = quantifyNeurometricOptostimMaps( ...
    responseMaps, templates, analysisMask, dataAudit.pixelsPerMM, options);

zeroContrastQC = calculateZeroContrastQC(responseMaps, stimContrastValues);

neurometricOptostim = struct();
neurometricOptostim.metadata = struct( ...
    'analysisMode', analysisMode, ...
    'created', datestr(now, 30), ...
    'monkeyName', monkeyName, ...
    'chamber', chamberWanted, ...
    'currentSessID', currentSessID, ...
    'blockID', blockID, ...
    'sessionLabel', options.sessionLabel, ...
    'spatialFilterCyclesPerMM', options.columnarBandCyclesPerMM, ...
    'lowpassCutoffCyclesPerMM', options.broadHighCutoffCyclesPerMM, ...
    'signConvention', 'A_COL_signed negative=0deg dominance, positive=90deg dominance');
neurometricOptostim.dataAudit = mergeStructs(dataAudit, conditionMatrixAudit);
neurometricOptostim.sourceFiles = filenameStruct;
neurometricOptostim.conditionDefinitions = conditionDefinitions;
neurometricOptostim.stimContrastValues = stimContrastValues;
neurometricOptostim.conditionIndices = conditionIndices;
neurometricOptostim.responseSources = responseSources;
neurometricOptostim.responseMaps = responseMaps;
neurometricOptostim.filteredMaps = filteredMaps;
neurometricOptostim.analysisMask = analysisMask;
neurometricOptostim.templates = templates;
neurometricOptostim.templates.templateQC = mergeStructs(templateAudit, templateQC);
neurometricOptostim.metrics = metrics;
neurometricOptostim.gaussianQC = mergeStructs(sharedFit, gaussianQC);
neurometricOptostim.ratioQC = ratioQC;
neurometricOptostim.zeroContrastQC = zeroContrastQC;
neurometricOptostim.options = options;

printNeurometricOptostimAudit(neurometricOptostim);

if plotFlag || saveFlag
    outputFiles = plotNeurometricOptostimQC(neurometricOptostim, options);
    neurometricOptostim.outputFiles = outputFiles;
end

if saveFlag
    if ~exist(options.outputDir, 'dir')
        mkdir(options.outputDir);
    end
    outFile = fullfile(options.outputDir, ...
        sprintf('%s_neurometricOptostim.mat', options.sessionLabel));
    save(outFile, 'neurometricOptostim', '-v7.3');
    neurometricOptostim.outputMatFile = outFile;
end
end

function options = defaultNeurometricOptostimOptions()
options = struct();
options.columnarBandCyclesPerMM = [0.8 3];
options.broadHighCutoffCyclesPerMM = 0.8;
options.maskErodeRadiusPx = 8;
options.centerBoundMM = 0.5;
options.widthFractionBound = 0.25;
options.rotationBoundDeg = 15;
options.flexCenterBoundMM = 0.5;
options.flexWidthFractionBound = 0.25;
options.flexRotationBoundDeg = 15;
options.fitFlexibleGaussian = true;
options.minTemplatePixels = 500;
options.templateOverlapWarning = 0.5;
options.signalWeightPercentile = 70;
options.dcStabilityMultiplier = 3;
options.epsilon = eps;
end

function dataAudit = auditNeurometricOptostimData(currentBlockStruct, filenameStruct, ...
    behavioralData, imagingData, bitmapData, columnarProducts, demoProducts, blockID)

dataAudit = struct();
dataAudit.session = currentBlockStruct;
dataAudit.PCAFieldNames = intersect(fieldnames(imagingData), ...
    {'ortpca'; 'ortampmap'; 'orts'; 'pcaexpl'; 'npca'; 'mask'; 'nanmask'});
dataAudit.bitmapFieldNames = intersect(fieldnames(bitmapData), ...
    {'columnarbitmap'; 'columnarbitmapCoreg'; 'columnarbitmapTFcamspace'; ...
    'orts'; 'transformParams'; 'transformClass'; 'gaussianContourLevel'});
dataAudit.demoProductFieldNames = fieldnames(demoProducts);
dataAudit.columnarProductFieldNames = fieldnames(columnarProducts);
dataAudit.ortpcaSize = sizeOrEmpty(imagingData, 'ortpca');
dataAudit.orts = squeeze(imagingData.orts(:,:,blockID));
dataAudit.coregisteredBitmapSize = sizeOrEmpty(bitmapData, 'columnarbitmapCoreg');
dataAudit.targetedColumnsCamspaceSize = sizeOrEmpty(demoProducts, 'targetedColumnsCamspace');
dataAudit.gaussianMaskCamspaceSize = sizeOrEmpty(demoProducts, 'gaussianMaskCamspace');
dataAudit.coordinateSystem = struct( ...
    'imagingData_ortpca', 'reference-session camera before coregistration', ...
    'bitmapData_columnarbitmapCoreg', 'current-session camera', ...
    'demoProducts_camspace', 'current-session camera', ...
    'responseMaps', 'current-session camera / integrated-response frame');
dataAudit.pixelsPerMM = 1 ./ imagingData.pixelsizemm(blockID);
dataAudit.pixelSizeMM = imagingData.pixelsizemm(blockID);
dataAudit.files = filenameStruct;
dataAudit.TSConditionFields = fieldnames(behavioralData.optoTS(blockID).Header.Conditions);
if isfield(behavioralData, 'baselineTS') && ...
        numel(behavioralData.baselineTS) >= blockID && ...
        ~isempty(behavioralData.baselineTS(blockID).Header)
    dataAudit.baselineMode = 'separate baseline file';
else
    dataAudit.baselineMode = 'combined baseline/opto file';
end
dataAudit.responseMapStatus = ...
    ['condition averages from getUsableTrials after blank subtraction; ' ...
     'trial-level maps retained in getUsableTrials output but not quantified here'];
dataAudit.temporalIntegrationStatus = ...
    ['integrated response maps loaded as DataCond from filenameStruct.Intg/' ...
     'baselineIntg when available; no local trial-level regrouping is performed'];
end

function [responseMaps, conditionIndices, conditionDefinitions, stimContrastValues, responseSources, audit] = ...
    buildNeurometricConditionMatrix(behavioralData, filenameStruct, blockID)

conditionDefinitions = { ...
    'Visual 0 + no optostimulation'; ...
    'Visual 0 + optostimulation of 0 columns'; ...
    'Visual 0 + optostimulation of 90 columns'; ...
    'Visual 90 + no optostimulation'; ...
    'Visual 90 + optostimulation of 0 columns'; ...
    'Visual 90 + optostimulation of 90 columns'};
conditionFields = {'V0', 'V0O0', 'V0O90', 'V90', 'V90O0', 'V90O90'};

responseSources = buildResponseSources(behavioralData, filenameStruct, blockID);
stimContrastValues = inferSourceContrastValues(responseSources(1));
nContrasts = numel(stimContrastValues);
if nContrasts == 0
    error('runNeurometricOptostim:NoStimContrasts', ...
        'No stimulus contrasts could be inferred from getUsableTrials condIDs.');
end

conditionIndices = struct();
for ii = 1:numel(conditionFields)
    conditionIndices.(conditionFields{ii}) = nan(nContrasts, 1);
end

firstMap = responseSources(1).DataCond(:,:,1);
responseMaps = nan(size(firstMap,1), size(firstMap,2), nContrasts, 6);
conditionMappingTable = initializeConditionMappingRows(nContrasts, 6);

for contrastIndex = 1:nContrasts
    contrastValue = stimContrastValues(contrastIndex);
    for conditionIndex = 1:6
        fieldName = conditionFields{conditionIndex};
        sourceNumber = sourceNumberForCondition(responseSources, fieldName);
        source = responseSources(sourceNumber);
        [dataCondIndex, matchedSourceContrast, nSourceContrastMatches] = ...
            selectConditionIndexForContrast(source, fieldName, contrastValue);

        if isnan(dataCondIndex)
            error('runNeurometricOptostim:MissingConditionMap', ...
                ['No getUsableTrials condition map found for contrast=%g, ' ...
                 'output column=%d, condition=%s.'], ...
                contrastValue, conditionIndex, fieldName);
        end

        responseMaps(:,:,contrastIndex,conditionIndex) = ...
            double(source.DataCond(:,:,dataCondIndex));
        conditionIndices.(fieldName)(contrastIndex) = dataCondIndex;
        conditionMappingTable((contrastIndex-1)*6 + conditionIndex) = ...
            makeConditionMappingRow(contrastValue, conditionIndex, ...
            conditionDefinitions{conditionIndex}, sourceNumber, source, ...
            fieldName, dataCondIndex, matchedSourceContrast, ...
            nSourceContrastMatches);
    end
end

zeroIdx = find(stimContrastValues == 0, 1);
if ~isempty(zeroIdx)
    zeroNoOptoDifference = responseMaps(:,:,zeroIdx,1) - responseMaps(:,:,zeroIdx,4);
else
    zeroNoOptoDifference = [];
end

audit = struct();
audit.nResponseSources = numel(responseSources);
audit.responseSourceSummary = summarizeResponseSources(responseSources);
audit.nStimContrasts = nContrasts;
audit.stimContrastValues = stimContrastValues;
audit.conditionIndexFields = conditionFields;
audit.conditionIndices = conditionIndices;
audit.conditionMappingTable = conditionMappingTable;
audit.conditionMappingTableColumns = { ...
    'contrast', 'outputColumn', 'outputCondition', 'sourceNumber', ...
    'sourceName', 'conditionField', 'dataCondIndex', 'trialCount', ...
    'matchedSourceContrast', 'nSourceContrastMatches'};
audit.zeroContrastNoOptoMeanAbsDifference = meanFinite(abs(zeroNoOptoDifference(:)));
audit.conditionIndexSource = ...
    'getUsableTrials condIDs matched by source StimCon and converted to raw DataCond indices';
audit.baselineSemantics = determineBaselineSemantics(responseSources);
end

function responseSources = buildResponseSources(behavioralData, filenameStruct, blockID)
mainSource = loadResponseSource('main', behavioralData.optoTS(blockID), ...
    filenameStruct.TS, filenameStruct.Intg);
responseSources = mainSource;

hasBaselineTS = isfield(behavioralData, 'baselineTS') && ...
    numel(behavioralData.baselineTS) >= blockID && ...
    isfield(behavioralData.baselineTS(blockID), 'Header') && ...
    ~isempty(behavioralData.baselineTS(blockID).Header);

if hasBaselineTS
    if ~isfield(filenameStruct, 'baselineIntg') || isempty(filenameStruct.baselineIntg)
        error('runNeurometricOptostim:MissingBaselineIntegratedResponseFile', ...
            ['A separate baselineTS is present (%s), but generateFilenames ' ...
             'did not find a paired baseline integrated-response file.'], ...
            filenameStruct.baselineTS);
    end
    responseSources(2) = loadResponseSource('separate baseline', ...
        behavioralData.baselineTS(blockID), filenameStruct.baselineTS, ...
        filenameStruct.baselineIntg);
end
end

function source = loadResponseSource(sourceName, TS, tsFilename, intgFilename)
if isempty(tsFilename) || exist(tsFilename, 'file') ~= 2
    error('runNeurometricOptostim:MissingTSFile', ...
        'Missing TS file for source %s: %s', sourceName, tsFilename);
end
if isempty(intgFilename) || exist(intgFilename, 'file') ~= 2
    error('runNeurometricOptostim:MissingIntegratedResponseFile', ...
        'Missing integrated-response file for source %s: %s', ...
end

availableVars = who('-file', intgFilename);
if ~any(strcmp(availableVars, 'DataCond'))
    error('runNeurometricOptostim:MissingIntegratedResponseVariable', ...
        ['Integrated-response file for source %s does not contain DataCond: %s. ' ...
         'Variables found: %s'], ...
        sourceName, intgFilename, formatAvailableVariables(availableVars));
end
loaded = load(intgFilename, 'DataCond');
DataCond = double(loaded.DataCond);
[condIDs, groupingAudit] = inferSourceConditionIDs(TS);
maxConditionIndex = maxConditionIndexFromCondIDs(condIDs) + numel(condIDs.blankConds);
if size(DataCond,3) < maxConditionIndex
    error('runNeurometricOptostim:DataCondConditionMismatch', ...
        ['Source %s DataCond has %d condition slices, but getUsableTrials ' ...
         'requires raw DataCond index %d.'], ...
        sourceName, size(DataCond,3), maxConditionIndex);
end

source = struct();
source.sourceName = sourceName;
source.TS = TS;
source.DataCond = DataCond;
source.condIDs = condIDs;
source.groupingAudit = groupingAudit;
source.tsFilename = tsFilename;
source.intgFilename = intgFilename;
source.loadedVariableName = 'DataCond';
source.availableVariables = {availableVars{:}};
source.dataCondConstruction = 'loaded directly from integrated-response MAT file';
source.isTemporallyIntegrated = ~isempty(strfind(intgFilename, 'Intg'));
source.baselineSubtractionStatus = ...
    'not modified here; DataCond is used as stored by the integrated-response file';
end

function varText = formatAvailableVariables(availableVars)
if isempty(availableVars)
    varText = '<none>';
else
    varText = strjoin(availableVars(:)', ', ');
end
end

function contrasts = inferSourceContrastValues(source)
candidateFields = {'V0O0', 'V0O90', 'V90O0', 'V90O90', 'V0', 'V90'};
contrasts = [];
for ii = 1:numel(candidateFields)
    fieldName = candidateFields{ii};
    if isfield(source.condIDs, fieldName) && ~isempty(source.condIDs.(fieldName))
        contrasts = sourceContrastsForConditionIDs(source, source.condIDs.(fieldName));
        if ~isempty(contrasts)
            break;
        end
    end
end
contrasts = unique(contrasts, 'stable');
end

function [condIDs, groupingAudit] = inferSourceConditionIDs(TS)
nCompletedTrials = countCompletedImagingTrials(TS);
dummyImages = zeros(1, 1, max(nCompletedTrials, 1));
[groupingTS, groupingAudit] = canonicalizeGroupingTSProjImg(TS);
[~, ~, condIDs] = getUsableTrials(groupingTS, dummyImages);
end

function [groupingTS, audit] = canonicalizeGroupingTSProjImg(TS)
groupingTS = TS;
audit = struct('nO045Canonicalized', 0, ...
    'nO135Canonicalized', 0, ...
    'canonicalizedAny', false, ...
    'status', 'no O045/O135 labels found');
if ~isfield(groupingTS.Header.Conditions, 'ProjImg')
    return;
end
projImg = groupingTS.Header.Conditions.ProjImg;
if iscell(projImg)
    for ii = 1:numel(projImg)
        if ischar(projImg{ii}) || isstring(projImg{ii})
            [projImg{ii}, nO045, nO135] = canonicalizeOneProjImgLabel(char(projImg{ii}));
            audit.nO045Canonicalized = audit.nO045Canonicalized + nO045;
            audit.nO135Canonicalized = audit.nO135Canonicalized + nO135;
        end
    end
elseif ischar(projImg) || isstring(projImg)
    [projImg, audit.nO045Canonicalized, audit.nO135Canonicalized] = ...
        canonicalizeOneProjImgLabel(char(projImg));
end
groupingTS.Header.Conditions.ProjImg = projImg;
audit.canonicalizedAny = audit.nO045Canonicalized > 0 || audit.nO135Canonicalized > 0;
if audit.canonicalizedAny
    audit.status = sprintf('canonicalized O045->O000 count=%d; O135->O090 count=%d', ...
        audit.nO045Canonicalized, audit.nO135Canonicalized);
end
end

function [labelOut, nO045, nO135] = canonicalizeOneProjImgLabel(labelIn)
nO045 = countSubstring(labelIn, 'O045');
nO135 = countSubstring(labelIn, 'O135');
labelOut = strrep(strrep(labelIn, 'O045', 'O000'), 'O135', 'O090');
end

function n = countSubstring(textValue, pattern)
n = numel(strfind(textValue, pattern));
end

function nCompletedTrials = countCompletedImagingTrials(TS)
trialOutcomeID = double([TS.Trial.Outcome]);
trialOIID = double([TS.Trial.FlagOIBLK]);
trialOIID(trialOIID == 0) = NaN;
trialOutcomeID = trialOutcomeID .* trialOIID;
nCompletedTrials = sum(trialOutcomeID == 10 | trialOutcomeID == 11 | ...
    trialOutcomeID == -10 | trialOutcomeID == -11);
end

function contrasts = sourceContrastsForConditionIDs(source, conditionIDs)
dataCondIndices = getUsableCondIDsToDataCondIndices(source, conditionIDs);
if isempty(dataCondIndices)
    contrasts = [];
    return;
end
contrasts = double(source.TS.Header.Conditions.StimCon(dataCondIndices));
contrasts = contrasts(:);
end

function dataCondIndices = getUsableCondIDsToDataCondIndices(source, conditionIDs)
% getUsableTrials subtracts nBlanks after constructing blank-removed images.average;
% integrated DataCond is stored raw TS-condition-indexed, as in plotDemoOptostim2.
conditionIDs = conditionIDs(:);
conditionIDs(isnan(conditionIDs)) = [];
dataCondIndices = conditionIDs + numel(source.condIDs.blankConds);
end

function maxConditionIndex = maxConditionIndexFromCondIDs(condIDs)
fields = fieldnames(condIDs);
maxConditionIndex = 0;
for ii = 1:numel(fields)
    values = condIDs.(fields{ii});
    if ~isempty(values)
        maxConditionIndex = max(maxConditionIndex, max(values(:)));
    end
end
end

function rows = initializeConditionMappingRows(nContrasts, nConditions)
emptyRow = struct( ...
    'contrast', NaN, ...
    'outputColumn', NaN, ...
    'outputCondition', '', ...
    'sourceNumber', NaN, ...
    'sourceName', '', ...
    'conditionField', '', ...
    'dataCondIndex', NaN, ...
    'trialCount', NaN, ...
    'matchedSourceContrast', NaN, ...
    'nSourceContrastMatches', NaN);
rows = repmat(emptyRow, nContrasts .* nConditions, 1);
end

function sourceNumber = sourceNumberForCondition(responseSources, fieldName)
if isBaselineConditionField(fieldName) && numel(responseSources) >= 2
    sourceNumber = 2;
else
    sourceNumber = 1;
end
end

function tf = isBaselineConditionField(fieldName)
tf = strcmp(fieldName, 'V0') || strcmp(fieldName, 'V90');
end

function [dataCondIndex, matchedSourceContrast, nMatches] = ...
    selectConditionIndexForContrast(source, fieldName, requestedContrast)
if ~isfield(source.condIDs, fieldName)
    dataCondIndex = NaN;
    matchedSourceContrast = NaN;
    nMatches = 0;
    return;
end
candidateCondIDs = source.condIDs.(fieldName);
candidateDataCondIndices = getUsableCondIDsToDataCondIndices(source, candidateCondIDs);
candidateContrasts = double(source.TS.Header.Conditions.StimCon(candidateDataCondIndices));
matchIdx = find(abs(candidateContrasts(:) - requestedContrast) <= eps(max(1, abs(requestedContrast))) .* 16);
nMatches = numel(matchIdx);
if nMatches == 0
    error('runNeurometricOptostim:MissingConditionContrast', ...
        ['No %s condition in source %s matches requested contrast %g. ' ...
         'Candidate DataCond indices: %s; candidate StimCon values: %s.'], ...
        fieldName, source.sourceName, requestedContrast, ...
        mat2str(candidateDataCondIndices(:)'), mat2str(candidateContrasts(:)'));
elseif nMatches > 1
    error('runNeurometricOptostim:AmbiguousConditionContrast', ...
        ['Multiple %s conditions in source %s match requested contrast %g. ' ...
         'Candidate DataCond indices: %s; candidate StimCon values: %s.'], ...
        fieldName, source.sourceName, requestedContrast, ...
        mat2str(candidateDataCondIndices(:)'), mat2str(candidateContrasts(:)'));
end
dataCondIndex = candidateDataCondIndices(matchIdx);
matchedSourceContrast = candidateContrasts(matchIdx);
end

function row = makeConditionMappingRow(contrastValue, outputColumn, outputCondition, ...
    sourceNumber, source, conditionField, dataCondIndex, matchedSourceContrast, ...
    nSourceContrastMatches)
row = initializeConditionMappingRows(1, 1);
row.contrast = contrastValue;
row.outputColumn = outputColumn;
row.outputCondition = outputCondition;
row.sourceNumber = sourceNumber;
row.sourceName = source.sourceName;
row.conditionField = conditionField;
row.dataCondIndex = dataCondIndex;
row.trialCount = countUsableTrialsForCondition(source.TS, dataCondIndex);
row.matchedSourceContrast = matchedSourceContrast;
row.nSourceContrastMatches = nSourceContrastMatches;
end

function nTrials = countUsableTrialsForCondition(TS, conditionIndex)
trialCondID = double([TS.Trial.CurrCond]);
trialOutcomeID = double([TS.Trial.Outcome]);
trialOIID = double([TS.Trial.FlagOIBLK]);
trialOIID(trialOIID == 0) = NaN;
trialOutcomeID = trialOutcomeID .* trialOIID;
usableTrialIdx = find(trialOutcomeID == 10 | trialOutcomeID == 11 | ...
    trialOutcomeID == -10 | trialOutcomeID == -11);
nTrials = sum(trialCondID(usableTrialIdx) == conditionIndex);
end

function summary = summarizeResponseSources(responseSources)
summary = struct([]);
for ii = 1:numel(responseSources)
    summary(ii).sourceNumber = ii;
    summary(ii).sourceName = responseSources(ii).sourceName;
    summary(ii).tsFilename = responseSources(ii).tsFilename;
    summary(ii).intgFilename = responseSources(ii).intgFilename;
    summary(ii).loadedVariableName = responseSources(ii).loadedVariableName;
    summary(ii).dataCondSize = size(responseSources(ii).DataCond);
    summary(ii).isTemporallyIntegrated = responseSources(ii).isTemporallyIntegrated;
    summary(ii).baselineSubtractionStatus = responseSources(ii).baselineSubtractionStatus;
    summary(ii).groupingAudit = responseSources(ii).groupingAudit;
end
end

function semantics = determineBaselineSemantics(responseSources)
if numel(responseSources) < 2
    semantics = 'No separate baseline source; no-opto cells are mapped from the main source.';
    return;
end
if ~isempty(responseSources(2).condIDs.V0) || ~isempty(responseSources(2).condIDs.V90)
    semantics = ['Separate baseline source contains measured no-opto ' ...
        'conditions from getUsableTrials and is used for columns 1 and 4.'];
else
    semantics = ['Separate baseline source is present, but no measured ' ...
        'V0/V90 no-opto conditions were returned by getUsableTrials.'];
end
end

function firstCond = findFirstFiniteCondition(conditionIndices)
names = fieldnames(conditionIndices);
firstCond = [];
for ii = 1:numel(names)
    values = conditionIndices.(names{ii});
    if ~isempty(values)
        firstCond = values(1);
        return;
    end
end
error('runNeurometricOptostim:NoConditions', ...
    'No nonblank condition indices were found.');
end

function [templates, analysisMask, templateAudit] = buildNeurometricTemplates( ...
    imagingData, bitmapData, columnarProducts, demoProducts, blockID, options)

currentCameraSize = size(bitmapData.columnarbitmapCoreg(:,:,1,blockID));
currentCameraRef = imref2d(currentCameraSize);
pcaOrts = squeeze(imagingData.orts(1,:,blockID));
pcaIdx0 = findOrientationIndex(pcaOrts, 0, 'imagingData.orts');
pcaIdx90 = findOrientationIndex(pcaOrts, 90, 'imagingData.orts');
coregOrts = squeeze(bitmapData.orts(1,:,blockID));
coregIdx0 = findOrientationIndex(coregOrts, 0, 'bitmapData.orts');
coregIdx90 = findOrientationIndex(coregOrts, 90, 'bitmapData.orts');

pca0Reference = double(imagingData.ortpca(:,:,pcaIdx0,blockID));
pca90Reference = double(imagingData.ortpca(:,:,pcaIdx90,blockID));
pca0Current = warpToCurrentCamera(pca0Reference, bitmapData, blockID, currentCameraRef);
pca90Current = warpToCurrentCamera(pca90Reference, bitmapData, blockID, currentCameraRef);
allPcaReference = double(imagingData.ortpca(:,:,:,blockID));
allPcaCurrent = nan([currentCameraSize size(allPcaReference,3)]);
for ii = 1:size(allPcaReference,3)
    allPcaCurrent(:,:,ii) = warpToCurrentCamera(allPcaReference(:,:,ii), ...
        bitmapData, blockID, currentCameraRef);
end

roiReference = isfinite(imagingData.nanmask(:,:,blockID)) & ...
    imagingData.nanmask(:,:,blockID) ~= 0;
roiMask = warpBinaryToCurrentCamera(roiReference, bitmapData, blockID, currentCameraRef);
validTransformMask = warpBinaryToCurrentCamera(true(size(roiReference)), ...
    bitmapData, blockID, currentCameraRef);
analysisMask = roiMask & validTransformMask;
analysisMask = imerode(analysisMask, strel('disk', options.maskErodeRadiusPx, 0));
if nnz(analysisMask) < options.minTemplatePixels
    warning('runNeurometricOptostim:SmallAnalysisMask', ...
        'Analysis mask has only %d pixels.', nnz(analysisMask));
end

allOrientationActivity = meanFiniteDim3(allPcaCurrent);
allOrientationActivity(~analysisMask) = NaN;
broadAllOrientation = FilterFermi2D(fillMissingWithMedian(allOrientationActivity, analysisMask), ...
    0, options.broadHighCutoffCyclesPerMM, imagingData.pixelsizemm(blockID));
broadAllOrientation(~analysisMask) = NaN;

referenceGaussian = fitReferenceGaussian(broadAllOrientation, analysisMask, ...
    imagingData.pixelsizemm(blockID));
referenceGaussian.map = referenceGaussian.map ./ max(referenceGaussian.map(analysisMask));

signedDifference = pca90Current - pca0Current;
signInfo = calibratePcaSign(pca0Current, pca90Current, signedDifference, analysisMask);
signedDifference = signInfo.signMultiplier .* signedDifference;
preference0 = max(-signedDifference, 0);
preference90 = max(signedDifference, 0);
preference0(~analysisMask) = 0;
preference90(~analysisMask) = 0;

templates = struct();
templates.G_reference = referenceGaussian.map;
templates.W_COL0 = [];
templates.W_COL90 = [];
templates.preference0Current = preference0;
templates.preference90Current = preference90;
templates.pca0Current = pca0Current;
templates.pca90Current = pca90Current;
templates.broadAllOrientation = broadAllOrientation;
templates.binaryProjectorMask0 = demoProducts.targetedColumnsCamspace(:,:,coregIdx0) > 0;
templates.binaryProjectorMask90 = demoProducts.targetedColumnsCamspace(:,:,coregIdx90) > 0;
templates.coregColumnMap0 = bitmapData.columnarbitmapCoreg(:,:,coregIdx0,blockID);
templates.coregColumnMap90 = bitmapData.columnarbitmapCoreg(:,:,coregIdx90,blockID);
templates.signCalibration = signInfo;

templateAudit = struct();
templateAudit.pcaIdx0 = pcaIdx0;
templateAudit.pcaIdx90 = pcaIdx90;
templateAudit.coregIdx0 = coregIdx0;
templateAudit.coregIdx90 = coregIdx90;
templateAudit.referenceGaussian = rmfield(referenceGaussian, 'map');
templateAudit.analysisMaskPixels = nnz(analysisMask);
templateAudit.analysisMaskErodeRadiusPx = options.maskErodeRadiusPx;
templateAudit.preferenceSource = ...
    'imagingData.ortpca 90-0 difference warped to current-session camera';
templateAudit.projectorMasksUsedAsAmplitudeTemplates = false;
end

function [W0, W90, templateQC] = buildColumnReadoutWeights(templates, analysisMask, options)
W0 = templates.G_shared .* templates.preference0Current;
W90 = templates.G_shared .* templates.preference90Current;
if any(W0(analysisMask) < 0) || any(W90(analysisMask) < 0)
    error('runNeurometricOptostim:SignedColumnWeight', ...
        'Column-domain weights must be nonnegative before normalization.');
end
W0(~analysisMask) = 0;
W90(~analysisMask) = 0;
if sum(W0(analysisMask)) <= 0 || sum(W90(analysisMask)) <= 0
    error('runNeurometricOptostim:EmptyColumnWeights', ...
        'Column-domain template weights are empty.');
end
W0 = W0 ./ sum(W0(analysisMask));
W90 = W90 ./ sum(W90(analysisMask));

w0 = W0(analysisMask);
w90 = W90(analysisMask);
templateQC = struct();
templateQC.spatialCorrelation = corr(w0(:), w90(:), 'rows', 'complete');
templateQC.softOverlapCoefficient = sum(min(w0(:), w90(:)));
templateQC.effectiveAreaCOL0Pixels = 1 ./ sum(w0(:).^2);
templateQC.effectiveAreaCOL90Pixels = 1 ./ sum(w90(:).^2);
templateQC.minCOL0 = min(w0(:));
templateQC.maxCOL0 = max(w0(:));
templateQC.sumCOL0 = sum(w0(:));
templateQC.nonzeroPixelCountCOL0 = nnz(w0(:) > 0);
templateQC.minCOL90 = min(w90(:));
templateQC.maxCOL90 = max(w90(:));
templateQC.sumCOL90 = sum(w90(:));
templateQC.nonzeroPixelCountCOL90 = nnz(w90(:) > 0);
[c0x, c0y] = weightedCenter(W0, analysisMask);
[c90x, c90y] = weightedCenter(W90, analysisMask);
templateQC.weightedCenterCOL0 = [c0x c0y];
templateQC.weightedCenterCOL90 = [c90x c90y];
templateQC.weightedCenterDistancePx = hypot(c90x - c0x, c90y - c0y);
templateQC.binaryMaskCorrelationCOL0 = safeCorr(w0(:), ...
    double(templates.binaryProjectorMask0(analysisMask)));
templateQC.binaryMaskCorrelationCOL90 = safeCorr(w90(:), ...
    double(templates.binaryProjectorMask90(analysisMask)));
templateQC.reducedSeparability = templateQC.softOverlapCoefficient > ...
    options.templateOverlapWarning;
if templateQC.reducedSeparability
    warning('runNeurometricOptostim:ReducedColumnSeparability', ...
        'Column templates overlap substantially (soft overlap %.3f).', ...
        templateQC.softOverlapCoefficient);
end
end

function [Gshared, sharedFit, gaussianQC] = fitSharedGaussianFootprint( ...
    responseMaps, Greference, analysisMask, pixelsPerMM, options)

referenceParams = gaussianParamsFromMap(Greference, analysisMask);
weights = signalWeights(responseMaps, analysisMask, options.signalWeightPercentile);
fitMaps = reshape(responseMaps, size(responseMaps,1), size(responseMaps,2), []);
for ii = 1:numel(weights)
    if weights(ii) <= 0
        fitMaps(:,:,ii) = NaN;
    end
end

lb = referenceParams;
ub = referenceParams;
pixelBound = options.centerBoundMM .* pixelsPerMM;
lb(1:2) = referenceParams(1:2) - pixelBound;
ub(1:2) = referenceParams(1:2) + pixelBound;
lb(3:4) = referenceParams(3:4) .* (1 - options.widthFractionBound);
ub(3:4) = referenceParams(3:4) .* (1 + options.widthFractionBound);
lb(5) = referenceParams(5) - deg2radLocal(options.rotationBoundDeg);
ub(5) = referenceParams(5) + deg2radLocal(options.rotationBoundDeg);

objective = @(p) sharedGaussianObjective(p, fitMaps, analysisMask, weights);
fitOptions = optimset('Display', 'off', 'MaxIter', 250, 'TolX', 1e-4, 'TolFun', 1e-4);
param0 = referenceParams(1:5);
lb = lb(1:5);
ub = ub(1:5);
boundedObjective = @(u) objective(lb + (ub - lb) ./ (1 + exp(-u)));
u0 = log((param0 - lb) ./ max(ub - param0, eps));
u0(~isfinite(u0)) = 0;
[uFit, fval, exitflag] = fminsearch(boundedObjective, u0, fitOptions);
pFit = lb + (ub - lb) ./ (1 + exp(-uFit));
Gshared = gaussianMapFromParams(pFit, size(Greference));
Gshared = Gshared ./ max(Gshared(analysisMask));

[A, B, R2, RMSE] = fitAmplitudeBackground(responseMaps, Gshared, analysisMask);
sharedFit = struct();
sharedFit.referenceParams = referenceParams(1:5);
sharedFit.sharedParams = pFit;
sharedFit.sharedBoundsLower = lb;
sharedFit.sharedBoundsUpper = ub;
sharedFit.sharedParameterNames = {'centerX', 'centerY', 'sigmaMajor', 'sigmaMinor', 'rotation'};
sharedFit.sharedBoundTolerance = 1e-3;
sharedFit.sharedReachedLowerBound = abs(pFit - lb) <= sharedFit.sharedBoundTolerance;
sharedFit.sharedReachedUpperBound = abs(pFit - ub) <= sharedFit.sharedBoundTolerance;
sharedFit.sharedReachedAnyBound = sharedFit.sharedReachedLowerBound | sharedFit.sharedReachedUpperBound;
sharedFit.objectiveValue = fval;
sharedFit.exitflag = exitflag;
sharedFit.weightRule = sprintf('Cells weighted by RMS above %.1f percentile; zero weights excluded from shape fit.', ...
    options.signalWeightPercentile);
sharedFit.weights = reshape(weights, size(responseMaps,3), size(responseMaps,4));

gaussianQC = struct();
gaussianQC.A_DC_shared = A;
gaussianQC.DC_background = B;
gaussianQC.DC_fit_R2 = R2;
gaussianQC.DC_fit_RMSE = RMSE;
gaussianQC.G_shared = Gshared;
end

function [filteredMaps, metrics, ratioQC] = quantifyNeurometricOptostimMaps( ...
    responseMaps, templates, analysisMask, pixelsPerMM, options)

pixelSizeMM = 1 ./ pixelsPerMM;
nContrasts = size(responseMaps, 3);
nConditions = size(responseMaps, 4);
maps3D = reshape(responseMaps, size(responseMaps,1), size(responseMaps,2), []);
broad3D = FilterFermi3D(fillMissingStack(maps3D, analysisMask), ...
    0, options.broadHighCutoffCyclesPerMM, pixelSizeMM);
columnar3D = FilterFermi3D(fillMissingStack(maps3D, analysisMask), ...
    options.columnarBandCyclesPerMM(1), options.columnarBandCyclesPerMM(2), pixelSizeMM);
broadMaps = reshape(broad3D, size(responseMaps));
columnarMaps = reshape(columnar3D, size(responseMaps));

filteredMaps = struct('broad', broadMaps, 'columnar', columnarMaps);
metrics = initializeMetricMatrices(nContrasts, nConditions);
ratioQC = struct();
ratioQC.ratioInvalidReason = repmat({''}, nContrasts, nConditions);

noiseSigma = estimateBackgroundNoise(broadMaps, analysisMask);
threshold = options.dcStabilityMultiplier .* noiseSigma;
ratioQC.dcStabilityThreshold = threshold;
ratioQC.noiseSigma = noiseSigma;

for ci = 1:nContrasts
    for ki = 1:nConditions
        result = quantifyConditionMap(responseMaps(:,:,ci,ki), ...
            broadMaps(:,:,ci,ki), columnarMaps(:,:,ci,ki), ...
            templates, analysisMask);
        names = fieldnames(result);
        for ni = 1:numel(names)
            if isfield(metrics, names{ni})
                metrics.(names{ni})(ci,ki) = result.(names{ni});
            end
        end
        stable = isfinite(result.A_DC_shared) && abs(result.A_DC_shared) > threshold;
        metrics.ratioIsStable(ci,ki) = stable;
        if stable
            denom = abs(result.A_DC_shared);
            metrics.COLsigned_to_DC(ci,ki) = result.A_COL_signed ./ denom;
            metrics.COLmagnitude_to_DC(ci,ki) = result.A_COL_magnitude ./ denom;
            metrics.COLenergy_to_DC(ci,ki) = result.A_COL_energy ./ denom;
            metrics.COL0abs_to_DC(ci,ki) = abs(result.M_COL0) ./ denom;
            metrics.COL90abs_to_DC(ci,ki) = abs(result.M_COL90) ./ denom;
        else
            ratioQC.ratioInvalidReason{ci,ki} = ...
                'abs(A_DC_shared) <= stability threshold';
        end
        metrics.columnDominance(ci,ki) = ...
            (result.M_COL90 - result.M_COL0) ./ ...
            (abs(result.M_COL90) + abs(result.M_COL0) + options.epsilon);
    end
end

if options.fitFlexibleGaussian
    flexible = fitFlexibleGaussianDiagnostics(broadMaps, templates.G_shared, ...
        analysisMask, pixelsPerMM, options);
    fNames = fieldnames(flexible);
    for fi = 1:numel(fNames)
        metrics.(fNames{fi}) = flexible.(fNames{fi});
    end
end
end

function result = quantifyConditionMap(responseMap, broadMap, columnarMap, templates, analysisMask)
[A, B, R2, RMSE] = fitAmplitudeBackground(broadMap, templates.G_shared, analysisMask);
Gnorm = templates.G_shared;
Gnorm(~analysisMask) = 0;
Gnorm = Gnorm ./ sum(Gnorm(analysisMask));
M0 = sum(templates.W_COL0(analysisMask) .* columnarMap(analysisMask));
M90 = sum(templates.W_COL90(analysisMask) .* columnarMap(analysisMask));
result = struct();
result.A_DC_shared = A;
result.DC_background = B;
result.DC_fit_R2 = R2;
result.DC_fit_RMSE = RMSE;
result.M_COL0 = M0;
result.M_COL90 = M90;
result.A_COL_signed = 0.5 .* (M90 - M0);
result.A_COL_magnitude = abs(result.A_COL_signed);
result.A_COL_common = 0.5 .* (M90 + M0);
result.A_COL_energy = sqrt(sum(Gnorm(analysisMask) .* columnarMap(analysisMask).^2));
result.DC_component_map = A .* templates.G_shared;
result.responseMean = meanFinite(responseMap(analysisMask));
end

function metrics = initializeMetricMatrices(nContrasts, nConditions)
metricNames = {'A_DC_shared', 'A_DC_flexible', 'DC_background', ...
    'DC_fit_R2', 'DC_fit_RMSE', 'M_COL0', 'M_COL90', 'A_COL_signed', ...
    'A_COL_magnitude', 'A_COL_common', 'A_COL_energy', ...
    'COL0abs_to_DC', 'COL90abs_to_DC', 'COLsigned_to_DC', ...
    'COLmagnitude_to_DC', 'COLenergy_to_DC', 'columnDominance', ...
    'centerShiftMM', 'majorAxisScale', 'minorAxisScale', ...
    'rotationDifferenceDeg', 'flexibleFitR2', 'deltaR2', 'responseMean'};
for ii = 1:numel(metricNames)
    metrics.(metricNames{ii}) = nan(nContrasts, nConditions);
end
metrics.ratioIsStable = false(nContrasts, nConditions);
end

function flexible = fitFlexibleGaussianDiagnostics(broadMaps, Gshared, ...
    analysisMask, pixelsPerMM, options)
nContrasts = size(broadMaps, 3);
nConditions = size(broadMaps, 4);
flexible = struct();
flexible.A_DC_flexible = nan(nContrasts, nConditions);
flexible.centerShiftMM = nan(nContrasts, nConditions);
flexible.majorAxisScale = nan(nContrasts, nConditions);
flexible.minorAxisScale = nan(nContrasts, nConditions);
flexible.rotationDifferenceDeg = nan(nContrasts, nConditions);
flexible.flexibleFitR2 = nan(nContrasts, nConditions);
flexible.deltaR2 = nan(nContrasts, nConditions);
sharedParams = gaussianParamsFromMap(Gshared, analysisMask);
for ci = 1:nContrasts
    for ki = 1:nConditions
        y = broadMaps(:,:,ci,ki);
        if rmsFinite(y(analysisMask)) <= 0
            continue;
        end
        p = fitSingleGaussianShape(y, sharedParams(1:5), analysisMask, pixelsPerMM, options);
        G = gaussianMapFromParams(p, size(Gshared));
        G = G ./ max(G(analysisMask));
        [A, ~, R2flex] = fitAmplitudeBackground(y, G, analysisMask);
        [~, ~, R2shared] = fitAmplitudeBackground(y, Gshared, analysisMask);
        flexible.A_DC_flexible(ci,ki) = A;
        flexible.centerShiftMM(ci,ki) = hypot(p(1)-sharedParams(1), ...
            p(2)-sharedParams(2)) ./ pixelsPerMM;
        flexible.majorAxisScale(ci,ki) = p(3) ./ sharedParams(3);
        flexible.minorAxisScale(ci,ki) = p(4) ./ sharedParams(4);
        flexible.rotationDifferenceDeg(ci,ki) = rad2degLocal(p(5)-sharedParams(5));
        flexible.flexibleFitR2(ci,ki) = R2flex;
        flexible.deltaR2(ci,ki) = R2flex - R2shared;
    end
end
end

function pFit = fitSingleGaussianShape(map, p0, analysisMask, pixelsPerMM, options)
pixelBound = options.flexCenterBoundMM .* pixelsPerMM;
lb = p0;
ub = p0;
lb(1:2) = p0(1:2) - pixelBound;
ub(1:2) = p0(1:2) + pixelBound;
lb(3:4) = p0(3:4) .* (1 - options.flexWidthFractionBound);
ub(3:4) = p0(3:4) .* (1 + options.flexWidthFractionBound);
lb(5) = p0(5) - deg2radLocal(options.flexRotationBoundDeg);
ub(5) = p0(5) + deg2radLocal(options.flexRotationBoundDeg);
objective = @(u) singleGaussianObjective(lb + (ub-lb)./(1+exp(-u)), map, analysisMask);
u0 = zeros(size(p0));
pFitU = fminsearch(objective, u0, optimset('Display', 'off', 'MaxIter', 100));
pFit = lb + (ub-lb)./(1+exp(-pFitU));
end

function err = singleGaussianObjective(p, map, analysisMask)
G = gaussianMapFromParams(p, size(map));
G = G ./ max(G(analysisMask));
[~, ~, ~, RMSE] = fitAmplitudeBackground(map, G, analysisMask);
err = RMSE;
end

function [A, B, R2, RMSE] = fitAmplitudeBackground(maps, G, analysisMask)
if ndims(maps) == 2
    y = maps(analysisMask);
    X = [G(analysisMask), ones(nnz(analysisMask),1)];
    valid = all(isfinite(X), 2) & isfinite(y);
    beta = X(valid,:) \ y(valid);
    yhat = X(valid,:) * beta;
    residual = y(valid) - yhat;
    A = beta(1);
    B = beta(2);
    R2 = 1 - sum(residual.^2) ./ max(sum((y(valid)-mean(y(valid))).^2), eps);
    RMSE = sqrt(mean(residual.^2));
else
    nContrasts = size(maps, 3);
    nConditions = size(maps, 4);
    A = nan(nContrasts, nConditions);
    B = nan(nContrasts, nConditions);
    R2 = nan(nContrasts, nConditions);
    RMSE = nan(nContrasts, nConditions);
    for ci = 1:nContrasts
        for ki = 1:nConditions
            [A(ci,ki), B(ci,ki), R2(ci,ki), RMSE(ci,ki)] = ...
                fitAmplitudeBackground(maps(:,:,ci,ki), G, analysisMask);
        end
    end
end
end

function referenceGaussian = fitReferenceGaussian(map, analysisMask, pixelSizeMM)
params = gaussianParamsFromMap(abs(map), analysisMask);
objective = @(p) singleGaussianObjective(p, abs(map), analysisMask);
pFit = fminsearch(objective, params(1:5), optimset('Display', 'off', 'MaxIter', 300));
G = gaussianMapFromParams(pFit, size(map));
[A, B, R2, RMSE] = fitAmplitudeBackground(abs(map), G, analysisMask);
referenceGaussian = struct('map', G, 'centerX', pFit(1), 'centerY', pFit(2), ...
    'sigmaMajorPx', pFit(3), 'sigmaMinorPx', pFit(4), ...
    'sigmaMajorMM', pFit(3).*pixelSizeMM, 'sigmaMinorMM', pFit(4).*pixelSizeMM, ...
    'rotationDeg', rad2degLocal(pFit(5)), ...
    'aspectRatio', pFit(3)./pFit(4), 'amplitude', A, ...
    'background', B, 'R2', R2, 'RMSE', RMSE);
end

function params = gaussianParamsFromMap(map, analysisMask)
[yy, xx] = ndgrid(1:size(map,1), 1:size(map,2));
w = abs(map);
w(~analysisMask | ~isfinite(w)) = 0;
if sum(w(:)) <= 0
    w = double(analysisMask);
end
w = w ./ sum(w(:));
cx = sum(xx(:).*w(:));
cy = sum(yy(:).*w(:));
x0 = xx - cx;
y0 = yy - cy;
Cxx = sum((x0(:).^2).*w(:));
Cyy = sum((y0(:).^2).*w(:));
Cxy = sum((x0(:).*y0(:)).*w(:));
[V,D] = eig([Cxx Cxy; Cxy Cyy]);
[evals, order] = sort(diag(D), 'descend');
V = V(:,order);
sigmaMajor = sqrt(max(evals(1), 1));
sigmaMinor = sqrt(max(evals(2), 1));
theta = atan2(V(2,1), V(1,1));
params = [cx cy sigmaMajor sigmaMinor theta];
end

function G = gaussianMapFromParams(p, mapSize)
[yy, xx] = ndgrid(1:mapSize(1), 1:mapSize(2));
x = xx - p(1);
y = yy - p(2);
ct = cos(p(5));
st = sin(p(5));
xp = x .* ct + y .* st;
yp = -x .* st + y .* ct;
G = exp(-0.5 .* ((xp ./ max(p(3), eps)).^2 + (yp ./ max(p(4), eps)).^2));
end

function err = sharedGaussianObjective(p, responseMaps, analysisMask, weights)
G = gaussianMapFromParams(p, [size(responseMaps,1) size(responseMaps,2)]);
G = G ./ max(G(analysisMask));
err = 0;
for ii = 1:numel(weights)
    if weights(ii) <= 0
        continue;
    end
    map = responseMaps(:,:,ii);
    [~, ~, ~, RMSE] = fitAmplitudeBackground(map, G, analysisMask);
    err = err + weights(ii) .* RMSE;
end
end

function weights = signalWeights(responseMaps, analysisMask, pct)
n = size(responseMaps,3) * size(responseMaps,4);
weights = nan(n,1);
for ii = 1:n
    map = responseMaps(:,:,ii);
    weights(ii) = rmsFinite(map(analysisMask));
end
cutoff = prctile(weights(isfinite(weights)), pct);
weights(weights < cutoff | ~isfinite(weights)) = 0;
if all(weights == 0)
    weights = ones(n,1);
end
weights = weights ./ max(weights);
end

function noiseSigma = estimateBackgroundNoise(broadMaps, analysisMask)
values = [];
for ii = 1:size(broadMaps,3)
    for jj = 1:size(broadMaps,4)
        map = broadMaps(:,:,ii,jj);
        edgeValues = map(analysisMask);
        values = [values; edgeValues(:)]; %#ok<AGROW>
    end
end
noiseSigma = 1.4826 .* medianFinite(abs(values - medianFinite(values)));
if ~isfinite(noiseSigma) || noiseSigma <= 0
    noiseSigma = stdFinite(values);
end
end

function signInfo = calibratePcaSign(pca0Current, pca90Current, signedDifference, analysisMask)
score0 = meanFinite(signedDifference(analysisMask) .* pca0Current(analysisMask));
score90 = meanFinite(signedDifference(analysisMask) .* pca90Current(analysisMask));
if score90 < score0
    signMultiplier = -1;
else
    signMultiplier = 1;
end
signInfo = struct('signMultiplier', signMultiplier, ...
    'score0BeforeSign', score0, 'score90BeforeSign', score90, ...
    'rule', 'Signed PCA difference is oriented so 90deg mapping response exceeds 0deg mapping response.');
end

function transformed = warpToCurrentCamera(image, bitmapData, blockID, currentCameraRef)
transformed = imwarp(image, bitmapData.transformParams{blockID}, ...
    'OutputView', currentCameraRef);
end

function transformed = warpBinaryToCurrentCamera(mask, bitmapData, blockID, currentCameraRef)
transformed = imwarp(double(mask), bitmapData.transformParams{blockID}, ...
    'nearest', 'OutputView', currentCameraRef) > 0.5;
end

function idx = findOrientationIndex(orts, target, sourceName)
ortVec = mod(double(orts(:)), 180);
target = mod(target, 180);
[delta, idx] = min(abs(ortVec - target));
if isempty(idx) || delta > 1e-6
    error('runNeurometricOptostim:MissingOrientation', ...
        'Could not find %.1f degree orientation in %s. Values: %s', ...
        target, sourceName, mat2str(ortVec(:)'));
end
end

function stack = fillMissingStack(stack, analysisMask)
for ii = 1:size(stack,3)
    stack(:,:,ii) = fillMissingWithMedian(stack(:,:,ii), analysisMask);
end
end

function map = fillMissingWithMedian(map, analysisMask)
fillValue = medianFinite(map(analysisMask));
if ~isfinite(fillValue)
    fillValue = 0;
end
map(~isfinite(map)) = fillValue;
end

function zeroContrastQC = calculateZeroContrastQC(responseMaps, stimContrastValues)
zeroContrastQC = struct();
zeroIdx = find(stimContrastValues == 0, 1);
if isempty(zeroIdx)
    zeroContrastQC.hasZeroContrast = false;
    return;
end
zeroContrastQC.hasZeroContrast = true;
zeroContrastQC.noOptoColumnsRetainedSeparately = true;
zeroContrastQC.noOptoAverage = 0.5 .* ...
    (responseMaps(:,:,zeroIdx,1) + responseMaps(:,:,zeroIdx,4));
zeroContrastQC.noOptoMeanAbsDifference = meanFinite(abs( ...
    responseMaps(:,:,zeroIdx,1) - responseMaps(:,:,zeroIdx,4)));
end

function printNeurometricOptostimAudit(result)
fprintf('\n=== neurometric-optostim implementation audit ===\n');
fprintf('Exact W_COL0/W_COL90 construction lines:\n');
fprintf('  W0 = templates.G_shared .* templates.preference0Current;\n');
fprintf('  W90 = templates.G_shared .* templates.preference90Current;\n');
fprintf('  preference0 = max(-signedDifference, 0);\n');
fprintf('  preference90 = max(signedDifference, 0);\n');
fprintf('  W0 = W0 ./ sum(W0(analysisMask));\n');
fprintf('  W90 = W90 ./ sum(W90(analysisMask));\n');
fprintf('Verification: templates are separate nonnegative soft-domain weights; no signed map is normalized by signed sum.\n');

qc = result.templates.templateQC;
fprintf('Template stats:\n');
fprintf('  W_COL0  min=%g max=%g sum=%g nonzero=%d weightedAreaPx=%g overlap=%g\n', ...
    qc.minCOL0, qc.maxCOL0, qc.sumCOL0, qc.nonzeroPixelCountCOL0, ...
    qc.effectiveAreaCOL0Pixels, qc.softOverlapCoefficient);
fprintf('  W_COL90 min=%g max=%g sum=%g nonzero=%d weightedAreaPx=%g overlap=%g\n', ...
    qc.minCOL90, qc.maxCOL90, qc.sumCOL90, qc.nonzeroPixelCountCOL90, ...
    qc.effectiveAreaCOL90Pixels, qc.softOverlapCoefficient);
fprintf('  spatialCorrelation=%g weightedCenterDistancePx=%g\n', ...
    qc.spatialCorrelation, qc.weightedCenterDistancePx);

pxPerMM = result.dataAudit.pixelsPerMM;
pxSizeMM = result.dataAudit.pixelSizeMM;
broadCutoff = result.options.broadHighCutoffCyclesPerMM;
band = result.options.columnarBandCyclesPerMM;
fprintf('Filter units for FilterFermi2D/3D:\n');
fprintf('  pixelsPerMM=%g; SizePxl=%g mm/pixel\n', pxPerMM, pxSizeMM);
fprintf('  broad LowCutOff=0 cycles/mm; broad HighCutOff=%g cycles/mm\n', broadCutoff);
fprintf('  columnar LowCutOff=%g cycles/mm; columnar HighCutOff=%g cycles/mm\n', band(1), band(2));
fprintf('  wavelength at %.3g cycles/mm = %.6g mm = %.6g pixels\n', ...
    band(1), 1 ./ band(1), pxPerMM ./ band(1));
fprintf('  wavelength at %.3g cycles/mm = %.6g mm = %.6g pixels\n', ...
    band(2), 1 ./ band(2), pxPerMM ./ band(2));

fprintf('Response sources:\n');
for ii = 1:numel(result.dataAudit.responseSourceSummary)
    src = result.dataAudit.responseSourceSummary(ii);
    fprintf('  source %d (%s): TS=%s\n', ...
        src.sourceNumber, src.sourceName, src.tsFilename);
    fprintf('    Intg=%s\n', src.intgFilename);
    fprintf('    variable=%s DataCondSize=%s temporallyIntegrated=%d\n', ...
        src.loadedVariableName, mat2str(src.dataCondSize), ...
        src.isTemporallyIntegrated);
    fprintf('    baselineSubtractionStatus=%s\n', src.baselineSubtractionStatus);
end
fprintf('Baseline semantics: %s\n', result.dataAudit.baselineSemantics);
fprintf('Response array dimensions: %s (rows x cols x nContrasts x 6 conditions)\n', ...
    mat2str(size(result.responseMaps)));

fprintf('Condition mapping table:\n');
fprintf(['  requestedContrast\toutCol\tsourceNo\tsourceName\tfield\t' ...
    'dataCond\ttrialCount\tmatchedSourceContrast\tnMatches\tcondition\n']);
rows = result.dataAudit.conditionMappingTable;
for ii = 1:numel(rows)
    fprintf('  %g\t%d\t%d\t%s\t%s\t%d\t%d\t%g\t%d\t%s\n', ...
        rows(ii).contrast, rows(ii).outputColumn, rows(ii).sourceNumber, ...
        rows(ii).sourceName, rows(ii).conditionField, ...
        rows(ii).dataCondIndex, rows(ii).trialCount, ...
        rows(ii).matchedSourceContrast, rows(ii).nSourceContrastMatches, ...
        rows(ii).outputCondition);
end

signInfo = result.templates.signCalibration;
fprintf('Mapping-derived sign calibration:\n');
fprintf('  signMultiplier=%g score0BeforeSign=%g score90BeforeSign=%g rule=%s\n', ...
    signInfo.signMultiplier, signInfo.score0BeforeSign, ...
    signInfo.score90BeforeSign, signInfo.rule);
[~, highContrastIndex] = max(result.stimContrastValues);
fprintf('  highContrast=%g visual-only A_COL_signed: col1(V0,noOpto)=%g col4(V90,noOpto)=%g\n', ...
    result.stimContrastValues(highContrastIndex), ...
    result.metrics.A_COL_signed(highContrastIndex,1), ...
    result.metrics.A_COL_signed(highContrastIndex,4));

fprintf('Shared Gaussian bounds:\n');
for ii = 1:numel(result.gaussianQC.sharedParameterNames)
    fprintf('  %s=%g lower=%g upper=%g reachedLower=%d reachedUpper=%d\n', ...
        result.gaussianQC.sharedParameterNames{ii}, ...
        result.gaussianQC.sharedParams(ii), ...
        result.gaussianQC.sharedBoundsLower(ii), ...
        result.gaussianQC.sharedBoundsUpper(ii), ...
        result.gaussianQC.sharedReachedLowerBound(ii), ...
        result.gaussianQC.sharedReachedUpperBound(ii));
end
fprintf('  anySharedGaussianParameterAtBound=%d\n', any(result.gaussianQC.sharedReachedAnyBound));

zeroIdx = find(result.stimContrastValues == 0, 1);
if ~isempty(zeroIdx)
    ratioFields = {'COLsigned_to_DC', 'COLmagnitude_to_DC', 'COLenergy_to_DC', ...
        'COL0abs_to_DC', 'COL90abs_to_DC'};
    unstableFinite = false;
    for ff = 1:numel(ratioFields)
        values = result.metrics.(ratioFields{ff})(zeroIdx,:);
        unstableFinite = unstableFinite || any(isfinite(values(~result.metrics.ratioIsStable(zeroIdx,:))));
    end
    fprintf('Zero-contrast ratio stability:\n');
    fprintf('  A_DC_shared: %s\n', mat2str(result.metrics.A_DC_shared(zeroIdx,:), 6));
    fprintf('  ratioIsStable: %s\n', mat2str(result.metrics.ratioIsStable(zeroIdx,:)));
    fprintf('  unstableFiniteRatioPresent=%d (expected 0; unstable ratios remain NaN)\n', unstableFinite);
else
    fprintf('Zero-contrast ratio stability: no exact zero contrast row found.\n');
end
fprintf('=== end neurometric-optostim implementation audit ===\n\n');
end

function outputFiles = plotNeurometricOptostimQC(result, options)
outputFiles = struct();
if options.saveFlag && ~exist(options.outputDir, 'dir')
    mkdir(options.outputDir);
end
fig1 = figure('Name', 'Neurometric optostim templates');
maps = {result.templates.broadAllOrientation, result.templates.G_reference, ...
    result.templates.G_shared, result.templates.preference0Current, ...
    result.templates.preference90Current, result.templates.W_COL0, ...
    result.templates.W_COL90, result.templates.W_COL90 - result.templates.W_COL0};
titles = {'Broad all-ort', 'Reference G', 'Shared G', 'Pref 0', ...
    'Pref 90', 'W COL0', 'W COL90', 'W90-W0'};
for ii = 1:numel(maps)
    subplot(3,3,ii);
    imagesc(maps{ii}); axis image off; colorbar; title(titles{ii});
end
subplot(3,3,9);
imagesc(result.analysisMask); axis image off; title('Analysis mask');
if options.saveFlag
    outputFiles.templateFigure = saveFigure(fig1, options, 'templates');
end

mapTypes = {'responseMaps', 'broad', 'columnar'};
for mt = 1:numel(mapTypes)
    if strcmp(mapTypes{mt}, 'responseMaps')
        data = result.responseMaps;
    else
        data = result.filteredMaps.(mapTypes{mt});
    end
    fig = figure('Name', ['Neurometric optostim ' mapTypes{mt}]);
    plotConditionMatrixMaps(data, result.stimContrastValues, ...
        result.conditionDefinitions, mapTypes{mt});
    if options.saveFlag
        outputFiles.([mapTypes{mt} 'Figure']) = saveFigure(fig, options, mapTypes{mt});
    end
end

metricNames = {'A_DC_shared', 'M_COL0', 'M_COL90', 'A_COL_signed', ...
    'A_COL_magnitude', 'A_COL_energy', 'COLsigned_to_DC', ...
    'COLmagnitude_to_DC', 'columnDominance'};
fig3 = figure('Name', 'Neurometric optostim scalar metrics');
for ii = 1:numel(metricNames)
    subplot(3,3,ii);
    values = result.metrics.(metricNames{ii});
    if contains(metricNames{ii}, '_to_DC')
        values(~result.metrics.ratioIsStable) = NaN;
    end
    imagesc(values); colorbar; title(metricNames{ii}, 'Interpreter', 'none');
    xlabel('Condition'); ylabel('Contrast');
    set(gca, 'XTick', 1:6, 'YTick', 1:numel(result.stimContrastValues));
end
if options.saveFlag
    outputFiles.scalarMetricFigure = saveFigure(fig3, options, 'scalar_metrics');
end

if isfield(result.metrics, 'flexibleFitR2')
    flexNames = {'centerShiftMM', 'majorAxisScale', 'minorAxisScale', ...
        'rotationDifferenceDeg', 'deltaR2'};
    fig4 = figure('Name', 'Neurometric optostim Gaussian sensitivity');
    for ii = 1:numel(flexNames)
        subplot(2,3,ii);
        imagesc(result.metrics.(flexNames{ii})); colorbar;
        title(flexNames{ii}, 'Interpreter', 'none');
    end
    if options.saveFlag
        outputFiles.gaussianSensitivityFigure = ...
            saveFigure(fig4, options, 'gaussian_sensitivity');
    end
end
end

function plotConditionMatrixMaps(data, contrasts, conditionDefinitions, mapType)
nContrasts = size(data,3);
nConditions = size(data,4);
finiteValues = data(isfinite(data));
clim = prctile(finiteValues, [2 98]);
for ci = 1:nContrasts
    for ki = 1:nConditions
        subplot(nContrasts, nConditions, (ci-1)*nConditions + ki);
        imagesc(data(:,:,ci,ki)); axis image off;
        if all(isfinite(clim)) && clim(1) < clim(2)
            caxis(clim);
        end
        if ci == 1
            title(conditionDefinitions{ki}, 'Interpreter', 'none');
        end
        if ki == 1
            ylabel(sprintf('Con %.3g', contrasts(ci)));
        end
    end
end
colormap(parula);
suplabel(mapType, 't');
end

function pathOut = saveFigure(fig, options, suffix)
pathOut = fullfile(options.outputDir, ...
    sprintf('%s_%s.png', options.sessionLabel, suffix));
saveas(fig, pathOut);
end

function s = sizeOrEmpty(parent, fieldName)
if isfield(parent, fieldName)
    s = size(parent.(fieldName));
else
    s = [];
end
end

function out = mergeStructs(a, b)
out = a;
names = fieldnames(b);
for ii = 1:numel(names)
    out.(names{ii}) = b.(names{ii});
end
end

function r = safeCorr(a, b)
valid = isfinite(a) & isfinite(b);
if nnz(valid) < 3
    r = NaN;
else
    r = corr(a(valid), b(valid));
end
end

function [cx, cy] = weightedCenter(W, mask)
[yy, xx] = ndgrid(1:size(W,1), 1:size(W,2));
w = W;
w(~mask) = 0;
w = w ./ sum(w(:));
cx = sum(xx(:).*w(:));
cy = sum(yy(:).*w(:));
end

function y = rmsFinite(x)
x = x(isfinite(x));
if isempty(x)
    y = NaN;
else
    y = sqrt(mean(x(:).^2));
end
end

function y = meanFinite(x)
x = x(isfinite(x));
if isempty(x)
    y = NaN;
else
    y = mean(x(:));
end
end

function y = meanFiniteDim3(x)
y = mean(x, 3, 'omitnan');
end

function y = medianFinite(x)
x = x(isfinite(x));
if isempty(x)
    y = NaN;
else
    y = median(x(:));
end
end

function y = stdFinite(x)
x = x(isfinite(x));
if numel(x) < 2
    y = NaN;
else
    y = std(x(:));
end
end
function y = deg2radLocal(x)
y = x .* pi ./ 180;
end

function y = rad2degLocal(x)
y = x .* 180 ./ pi;
end
