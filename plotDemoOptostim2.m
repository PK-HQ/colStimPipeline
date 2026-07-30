function [figHandle, demoPlotData, outputFiles] = plotDemoOptostim2( ...
    currentBlockStruct, imagingData, bitmapData, ...
    columnarProducts, demoProducts, blockID, currentSessID, ...
    mainPath, monkeyName, chamberWanted, saveFlag)
%PLOTDEMOOPTOSTIM2 Four-panel visual/opto demonstration figure.
%
% A: PCA-denoised visual response, 90 deg - 0 deg.
% B: Highest-contrast visual-only response, 90 deg - 0 deg, loaded from
%    the baseline run specified by currentBlockStruct.baselineTS.
% C: Zero-contrast opto response, OptoStim 90 deg - OptoStim 0 deg.
% D: Fitted model of Panel C.
%
% The validated calculation/registration path is retained verbatim in the
% +demoLegacy package. This wrapper adds Panel B and redraws the output.

fprintf('Running plotDemoOptostim2 version 2026-07-30-four-panel-visual-v2\n');

originalVisibility = get(groot, 'DefaultFigureVisible');
set(groot, 'DefaultFigureVisible', 'off');
legacyFig = [];
try
    [legacyFig, demoPlotData, outputFiles] = ...
        demoLegacy.plotDemoOptostim2( ...
            currentBlockStruct, imagingData, bitmapData, ...
            columnarProducts, demoProducts, blockID, currentSessID, ...
            mainPath, monkeyName, chamberWanted, saveFlag);
catch ME
    set(groot, 'DefaultFigureVisible', originalVisibility);
    if ~isempty(legacyFig) && ishghandle(legacyFig)
        delete(legacyFig);
    end
    rethrow(ME);
end
set(groot, 'DefaultFigureVisible', originalVisibility);
if ~isempty(legacyFig) && ishghandle(legacyFig)
    delete(legacyFig);
end

[visualTrialStructureFile, visualIntegratedResponseFile, visualSourceInfo] = ...
    resolveBaselineVisualFiles( ...
        currentBlockStruct, ...
        demoPlotData.trialStructureFile, ...
        demoPlotData.integratedResponseFile);

[visualDifferenceRaw, visualInfo] = loadHighestContrastVisualDifference( ...
    visualTrialStructureFile, visualIntegratedResponseFile);
visualInfo.source = visualSourceInfo;

referenceSize = size(demoPlotData.pcaDifference);
if ~isequal(size(visualDifferenceRaw), referenceSize)
    error('plotDemoOptostim:VisualResponseSizeMismatch', ...
        'Visual response size [%s] does not match map size [%s].', ...
        num2str(size(visualDifferenceRaw)), num2str(referenceSize));
end

bandpassSF = [0.8 3];
if isfield(demoPlotData, 'bandpassSFcyclesPerMM') && ...
        numel(demoPlotData.bandpassSFcyclesPerMM) == 2
    bandpassSF = double(demoPlotData.bandpassSFcyclesPerMM(:)');
end

% Reuse the selected activity transform from the validated optostim path so
% Panel B is displayed in the same Panel-A camera coordinates as Panels C/D.
activityTform = affine2d(double(demoPlotData.activityTransformMatrix));
visualDifferenceBandpassUnaligned = bandpassMap( ...
    visualDifferenceRaw, bandpassSF, visualInfo.imagingSizePxl);
visualDifferenceBandpass = transformMap( ...
    visualDifferenceBandpassUnaligned, activityTform, referenceSize);
visualActivity0CoregRaw = transformMap( ...
    visualInfo.response0, activityTform, referenceSize);
visualActivity90CoregRaw = transformMap( ...
    visualInfo.response90, activityTform, referenceSize);
visualDifferenceCoregRaw = ...
    visualActivity90CoregRaw - visualActivity0CoregRaw;

optoDifferenceBandpass = ...
    double(demoPlotData.zeroContrastOptoDifferenceBandpass);
modelDifferenceBandpass = ...
    double(demoPlotData.modelActivityDifference90Minus0Bandpass);
comparisonMask = logical(demoPlotData.gaussianMaskLargestCurrent);
[~, r2AB, nAB] = mapCorrelation( ...
    demoPlotData.pcaDifference, visualDifferenceBandpass, comparisonMask);
[~, r2AC, nAC] = mapCorrelation( ...
    demoPlotData.pcaDifference, optoDifferenceBandpass, comparisonMask);
[~, r2CD, nCD] = mapCorrelation( ...
    optoDifferenceBandpass, modelDifferenceBandpass, comparisonMask);

demoPlotData.visualTrialStructureFile = visualTrialStructureFile;
demoPlotData.visualIntegratedResponseFile = visualIntegratedResponseFile;
demoPlotData.visualSourceInfo = visualSourceInfo;
demoPlotData.highestContrastVisualConditionInfo = visualInfo;
demoPlotData.highestContrastVisualDifferenceUnaligned = visualDifferenceRaw;
demoPlotData.highestContrastVisualDifferenceBandpassUnaligned = ...
    visualDifferenceBandpassUnaligned;
demoPlotData.highestContrastVisualActivity0CoregRaw = ...
    visualActivity0CoregRaw;
demoPlotData.highestContrastVisualActivity90CoregRaw = ...
    visualActivity90CoregRaw;
demoPlotData.highestContrastVisualDifference90Minus0CoregRaw = ...
    visualDifferenceCoregRaw;
demoPlotData.highestContrastVisualDifference90Minus0CoregBandpass = ...
    visualDifferenceBandpass;
demoPlotData.fourPanelCorrelations = struct( ...
    'rSquaredAB_pcaVisual', r2AB, 'nAB_pcaVisual', nAB, ...
    'rSquaredAC_pcaOpto', r2AC, 'nAC_pcaOpto', nAC, ...
    'rSquaredCD_optoModel', r2CD, 'nCD_optoModel', nCD);

fprintf('\n--- Four-panel demo additions ---\n');
fprintf('Visual baseline run:              %s\n', ...
    visualSourceInfo.baselineRun);
fprintf('Visual TS file:                   %s\n', ...
    visualTrialStructureFile);
fprintf('Visual integrated-response file:  %s\n', ...
    visualIntegratedResponseFile);
fprintf('Highest visual contrast:          %.6g\n', ...
    visualInfo.highestCommonContrast);
fprintf('Visual 0-deg condition(s):        [%s]\n', ...
    num2str(visualInfo.conditionIndices0));
fprintf('Visual 90-deg condition(s):       [%s]\n', ...
    num2str(visualInfo.conditionIndices90));
fprintf('Panel A vs B Pearson R^2:         %.5f (n=%d)\n', r2AB, nAB);
fprintf('Panel A vs C Pearson R^2:         %.5f (n=%d)\n', r2AC, nAC);
fprintf('Panel C vs D Pearson R^2:         %.5f (n=%d)\n', r2CD, nCD);
fprintf('Target colors:                    cyan=0 deg, red=90 deg\n');
fprintf('Target blob outlines:             none\n');
fprintf('----------------------------------\n\n');

figHandle = figure( ...
    'Name', sprintf('demo-optostim_sess%d', currentSessID), ...
    'Color', 'white', ...
    'Units', 'inches', ...
    'Position', [0.10 0.55 36.0 7.4]);

xMM = demoPlotData.xMM;
yMM = demoPlotData.yMM;
sharedLimits = demoPlotData.pcaColorLimits;
roiMask = logical(demoPlotData.roiMaskCurrent);
gaussMask = logical(demoPlotData.gaussianMaskLargestCurrent);
target0 = logical(demoPlotData.targetDesigned0Current);
target90 = logical(demoPlotData.targetDesigned90Current);
pixelSizeMM = mean(abs([median(diff(xMM(:))); median(diff(yMM(:)))]));
gaussOrientation = demoPlotData.largestGaussianOrientation;

maps = { ...
    demoPlotData.pcaDifference, ...
    fillInvalid(visualDifferenceBandpass), ...
    fillInvalid(optoDifferenceBandpass), ...
    fillInvalid(modelDifferenceBandpass)};
colorbarLabels = { ...
    ['PCA response: 90' char(176) ' - 0' char(176)], ...
    ['Bandpassed visual response: 90' char(176) ' - 0' char(176)], ...
    ['Bandpassed response: OptoStim 90' char(176) ...
     ' - OptoStim 0' char(176)], ...
    ['Bandpassed model response: OptoStim 90' char(176) ...
     ' - OptoStim 0' char(176)]};
titles = { ...
    ['PCA-denoised response: 90' char(176) ' - 0' char(176)], ...
    sprintf(['Highest-contrast visual (run %s): 90%c - 0%c ' ...
             '(contrast %.4g; %.1f-%.1f cycles/mm)'], ...
        visualSourceInfo.baselineRun, char(176), char(176), ...
        visualInfo.highestCommonContrast, bandpassSF(1), bandpassSF(2)), ...
    sprintf(['0%% contrast: OptoStim 90%c - OptoStim 0%c ' ...
             '(%.1f-%.1f cycles/mm)'], ...
        char(176), char(176), bandpassSF(1), bandpassSF(2)), ...
    sprintf(['Model: OptoStim 90%c - OptoStim 0%c ' ...
             '(%.1f-%.1f cycles/mm)'], ...
        char(176), char(176), bandpassSF(1), bandpassSF(2))};

axesHandles = gobjects(1,4);
panelLabels = gobjects(1,4);
for panelIndex = 1:4
    axesHandles(panelIndex) = subplot(1,4,panelIndex);
    drawPanel(axesHandles(panelIndex), maps{panelIndex}, ...
        xMM, yMM, sharedLimits, target0, target90, ...
        roiMask, gaussMask, pixelSizeMM, ...
        colorbarLabels{panelIndex}, titles{panelIndex}, gaussOrientation);
    panelLabels(panelIndex) = text(axesHandles(panelIndex), ...
        -0.13, 1.08, char(64 + panelIndex), ...
        'Units', 'normalized', 'FontName', 'Arial', ...
        'FontSize', 18, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
        'Clipping', 'off');
end

if exist('upFontSize', 'file') == 2
    for panelIndex = 1:4
        axes(axesHandles(panelIndex)); %#ok<LAXES>
        upFontSize(18, 0.01);
    end
end
set(axesHandles, 'LineWidth', 2, 'TickDir', 'out', ...
    'TickLength', [0.01 0.01], 'FontName', 'Arial');
set(panelLabels, 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Arial');
panelX = [0.035 0.285 0.535 0.785];
for panelIndex = 1:4
    set(axesHandles(panelIndex), ...
        'Position', [panelX(panelIndex) 0.265 0.145 0.585]);
end

stimulationSummary = demoPlotData.stimulationSummary;
summaryText = sprintf([ ...
    'Columns: %.0f    Pixels ON: %.0f    Total power: %.3f mW    ' ...
    'ROI power density: %.3f mW/mm^2    ' ...
    'R^2_{AB} (PCA-visual): %.2f    ' ...
    'R^2_{AC} (PCA-opto): %.2f    ' ...
    'R^2_{CD} (opto-model): %.2f'], ...
    stimulationSummary.meanColumns, stimulationSummary.meanPixelsON, ...
    stimulationSummary.meanTotalPower_mW, ...
    stimulationSummary.meanPowerDensityWithinROI_mWmm2, ...
    r2AB, r2AC, r2CD);
annotation(figHandle, 'textbox', [0.015 0.015 0.970 0.060], ...
    'String', summaryText, 'Interpreter', 'tex', ...
    'FitBoxToText', 'off', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'FontName', 'Arial', ...
    'FontSize', 10.5, 'LineStyle', 'none', 'EdgeColor', 'none', ...
    'BackgroundColor', 'none', 'Margin', 1);
drawnow;

% The legacy call already created the established output paths and MAT
% contract. Overwrite only the figure files with the new four-panel figure,
% then append the visual maps to the existing mapData structure.
if saveFlag == 1
    set(figHandle, 'PaperPositionMode', 'auto');
    print(figHandle, outputFiles.svg, '-dsvg');
    print(figHandle, outputFiles.png, '-dpng', '-r300');
    savedData = load(outputFiles.mat, 'mapData', 'stimulationSummary');
    mapData = savedData.mapData;
    stimulationSummary = savedData.stimulationSummary;
    mapData.visualBaselineRun = visualSourceInfo.baselineRun;
    mapData.visualTrialStructureFile = visualTrialStructureFile;
    mapData.visualIntegratedResponseFile = visualIntegratedResponseFile;
    mapData.visualHighestContrast = visualInfo.highestCommonContrast;
    mapData.visualConditionIndices0 = visualInfo.conditionIndices0;
    mapData.visualConditionIndices90 = visualInfo.conditionIndices90;
    mapData.visualActivity0CoregRaw = visualActivity0CoregRaw;
    mapData.visualActivity90CoregRaw = visualActivity90CoregRaw;
    mapData.visualActivityDifference90Minus0CoregRaw = ...
        visualDifferenceCoregRaw;
    mapData.visualActivityDifference90Minus0CoregBandpass = ...
        visualDifferenceBandpass;
    stimulationSummary.visualBaselineRun = visualSourceInfo.baselineRun;
    stimulationSummary.highestVisualContrast = ...
        visualInfo.highestCommonContrast;
    stimulationSummary.fourPanelCorrelations = ...
        demoPlotData.fourPanelCorrelations;
    save(outputFiles.mat, 'mapData', 'stimulationSummary', '-v7.3');
    fprintf('Replaced SVG/PNG with four-panel figure and updated MAT maps.\n');
end
end


function [tsFile, responseFile, info] = resolveBaselineVisualFiles( ...
        currentBlockStruct, optoTsFile, optoResponseFile)
% Resolve the visual-only TS/DataCond pair from the baseline run recorded in
% the experiment metadata. For 20230906R2, baselineTS='1', so this resolves
% the visual source to 20230906R1 while leaving the optostim source at R2.

requiredFields = {'monkeyNo', 'date', 'run', 'baselineTS'};
for fieldIndex = 1:numel(requiredFields)
    fieldName = requiredFields{fieldIndex};
    if ~isfield(currentBlockStruct, fieldName) || ...
            isempty(currentBlockStruct.(fieldName))
        error('plotDemoOptostim:MissingBaselineMetadata', ...
            'currentBlockStruct.%s is required to resolve the visual run.', ...
            fieldName);
    end
end

monkeyNo = scalarText(currentBlockStruct.monkeyNo, 'monkeyNo');
sessionDate = scalarText(currentBlockStruct.date, 'date');
currentRun = scalarText(currentBlockStruct.run, 'run');
baselineRun = scalarText(currentBlockStruct.baselineTS, 'baselineTS');

optoRunDirectory = fileparts(optoTsFile);
sessionDirectory = fileparts(optoRunDirectory);
baselineDirectory = fullfile(sessionDirectory, ['run' baselineRun]);
if exist(baselineDirectory, 'dir') ~= 7
    error('plotDemoOptostim:MissingBaselineRunDirectory', ...
        'Baseline visual run directory does not exist: %s', ...
        baselineDirectory);
end

baselinePrefix = sprintf('M%sD%sR%s', ...
    monkeyNo, sessionDate, baselineRun);
currentPrefix = sprintf('M%sD%sR%s', ...
    monkeyNo, sessionDate, currentRun);

tsExact = fullfile(baselineDirectory, [baselinePrefix 'TS.mat']);
tsFile = resolveSingleBaselineFile( ...
    tsExact, fullfile(baselineDirectory, [baselinePrefix '*TS.mat']), ...
    'trial-structure');

[~, optoResponseBase, optoResponseExtension] = fileparts(optoResponseFile);
if strncmp(optoResponseBase, currentPrefix, numel(currentPrefix))
    responseSuffix = optoResponseBase(numel(currentPrefix)+1:end);
else
    responseSuffix = 'StabIntgS004E023';
end
responseExact = fullfile(baselineDirectory, ...
    [baselinePrefix responseSuffix optoResponseExtension]);
responseFile = resolveSingleBaselineFile( ...
    responseExact, ...
    fullfile(baselineDirectory, [baselinePrefix 'StabIntg*.mat']), ...
    'integrated-response');

info = struct();
info.currentRun = currentRun;
info.baselineRun = baselineRun;
info.baselineDirectory = baselineDirectory;
info.trialStructureFile = tsFile;
info.integratedResponseFile = responseFile;
end


function filePath = resolveSingleBaselineFile(exactPath, fallbackPattern, label)
if exist(exactPath, 'file') == 2
    filePath = exactPath;
    return;
end

matches = dir(fallbackPattern);
matches = matches(~[matches.isdir]);
if isempty(matches)
    error('plotDemoOptostim:MissingBaselineVisualFile', ...
        ['Could not find the baseline visual %s file. Tried:\n%s\n' ...
         'Fallback pattern:\n%s'], ...
        label, exactPath, fallbackPattern);
end
if numel(matches) > 1
    names = {matches.name};
    error('plotDemoOptostim:AmbiguousBaselineVisualFile', ...
        ['Found multiple baseline visual %s files matching:\n%s\n' ...
         'Matches: %s'], ...
        label, fallbackPattern, strjoin(names, ', '));
end
filePath = fullfile(matches(1).folder, matches(1).name);
end


function value = scalarText(rawValue, fieldName)
if iscell(rawValue)
    if numel(rawValue) ~= 1
        error('plotDemoOptostim:InvalidBaselineMetadata', ...
            'currentBlockStruct.%s must be scalar.', fieldName);
    end
    rawValue = rawValue{1};
end
if isnumeric(rawValue)
    if ~isscalar(rawValue) || ~isfinite(rawValue)
        error('plotDemoOptostim:InvalidBaselineMetadata', ...
            'currentBlockStruct.%s must be a finite scalar.', fieldName);
    end
    value = num2str(rawValue);
elseif ischar(rawValue)
    value = strtrim(rawValue);
elseif isstring(rawValue) && isscalar(rawValue)
    value = strtrim(char(rawValue));
else
    error('plotDemoOptostim:InvalidBaselineMetadata', ...
        'Unsupported currentBlockStruct.%s value.', fieldName);
end
if isempty(value)
    error('plotDemoOptostim:InvalidBaselineMetadata', ...
        'currentBlockStruct.%s is empty.', fieldName);
end
end


function [differenceMap, info] = ...
        loadHighestContrastVisualDifference(tsFile, responseFile)
% Match the visual-only conventions used by the neurometric analysis:
% ProjImg contains 'Dot', GaborOrt selects 0/90, and TypeCond is 3 when
% that field is available. The files passed here come from baselineTS.

loadedTS = load(tsFile, 'TS');
loadedResponse = load(responseFile, 'DataCond');
TS = loadedTS.TS;
DataCond = double(loadedResponse.DataCond);
conditions = TS.Header.Conditions;
contrast = double(conditions.StimCon(:));
orientation = double(conditions.GaborOrt(:));
nConditions = numel(contrast);
if numel(orientation) ~= nConditions || size(DataCond,3) ~= nConditions
    error('plotDemoOptostim:VisualConditionCountMismatch', ...
        'StimCon, GaborOrt, and DataCond must have matching condition counts.');
end
images = normalizeImages(conditions.ProjImg, nConditions);

isVisual = false(nConditions,1);
for conditionIndex = 1:nConditions
    isVisual(conditionIndex) = contains(lower(images{conditionIndex}), 'dot');
end
if isfield(conditions, 'TypeCond') && ...
        numel(conditions.TypeCond) == nConditions
    isVisual = isVisual & double(conditions.TypeCond(:)) == 3;
end
is0 = isVisual & orientationMatch(orientation, 0);
is90 = isVisual & orientationMatch(orientation, 90);
commonContrasts = intersect( ...
    unique(contrast(is0 & isfinite(contrast))), ...
    unique(contrast(is90 & isfinite(contrast))));
if isempty(commonContrasts)
    error('plotDemoOptostim:MissingMatchedVisualContrast', ...
        ['No visual-only contrast is shared by 0- and 90-degree ' ...
         'conditions in baseline files:\nTS: %s\nDataCond: %s'], ...
        tsFile, responseFile);
end
highestContrast = max(commonContrasts);
tolerance = max(1e-10, abs(highestContrast) * 1e-10);
indices0 = find(is0 & abs(contrast - highestContrast) <= tolerance);
indices90 = find(is90 & abs(contrast - highestContrast) <= tolerance);
response0 = mean(DataCond(:,:,indices0), 3, 'omitnan');
response90 = mean(DataCond(:,:,indices90), 3, 'omitnan');
differenceMap = response90 - response0;

info = struct();
info.highestCommonContrast = highestContrast;
info.conditionIndices0 = indices0(:)';
info.conditionIndices90 = indices90(:)';
info.response0 = response0;
info.response90 = response90;
info.imagingSizePxl = double(TS.Header.Imaging.SizePxl);
info.trialStructureFile = tsFile;
info.integratedResponseFile = responseFile;
end


function images = normalizeImages(rawImages, nConditions)
if iscell(rawImages)
    images = cell(numel(rawImages),1);
    for imageIndex = 1:numel(rawImages)
        if isempty(rawImages{imageIndex})
            images{imageIndex} = '';
        else
            images{imageIndex} = char(string(rawImages{imageIndex}));
        end
    end
elseif isstring(rawImages)
    images = cellstr(rawImages(:));
elseif ischar(rawImages) && size(rawImages,1) == nConditions
    images = cellstr(rawImages);
else
    error('plotDemoOptostim:InvalidProjImg', ...
        'Unsupported TS.Header.Conditions.ProjImg format.');
end
if numel(images) ~= nConditions
    error('plotDemoOptostim:ProjImgCountMismatch', ...
        'ProjImg has %d entries; expected %d.', numel(images), nConditions);
end
images = images(:);
end


function matches = orientationMatch(orientations, target)
orientations = mod(double(orientations), 180);
difference = abs(orientations - mod(target,180));
matches = min(difference, 180 - difference) <= 1e-6;
end


function outputMap = bandpassMap(inputMap, bandpassSF, imagingSizePxl)
validMask = isfinite(inputMap);
filterInput = double(inputMap);
filterInput(~validMask) = 0;
outputMap = real(FilterFermi3D( ...
    filterInput, bandpassSF(1), bandpassSF(2), imagingSizePxl));
outputMap(~validMask) = NaN;
end


function outputMap = transformMap(inputMap, tform, outputSize)
validMask = isfinite(inputMap);
inputFilled = double(inputMap);
inputFilled(~validMask) = 0;
outputRef = imref2d(outputSize);
outputMap = imwarp(inputFilled, tform, 'linear', ...
    'OutputView', outputRef, 'FillValues', 0);
outputValidity = imwarp(double(validMask), tform, 'linear', ...
    'OutputView', outputRef, 'FillValues', 0);
outputMap(outputValidity < 0.5) = NaN;
end


function outputMap = fillInvalid(inputMap)
outputMap = inputMap;
outputMap(~isfinite(outputMap)) = 0;
end


function [r, rSquared, nPixels] = mapCorrelation(map1, map2, mask)
validMask = logical(mask) & isfinite(map1) & isfinite(map2);
nPixels = nnz(validMask);
if nPixels < 3
    r = NaN;
    rSquared = NaN;
    return;
end
values1 = double(map1(validMask));
values2 = double(map2(validMask));
correlationMatrix = corrcoef(values1, values2);
r = correlationMatrix(1,2);
rSquared = r.^2;
end


function drawPanel(ax, imageData, xMM, yMM, colorLimits, ...
        target0, target90, roiMask, gaussMask, pixelSizeMM, ...
        colorbarLabel, titleText, gaussOrientation)

imageHandle = imagesc(ax, xMM, yMM, imageData);
imageHandle.AlphaData = isfinite(imageData);
imageHandle.HandleVisibility = 'off';
set(ax, 'Color', [0.82 0.82 0.82], 'CLim', colorLimits);
colormap(ax, gray(257));
hold(ax, 'on');

[targetHandle0, targetHandle90] = ...
    plotTargets(ax, target0, target90, xMM, yMM);
roiHandle = plotMask(ax, roiMask, pixelSizeMM, '-', 2.10, 'ROI');
gaussHandle = plotMask(ax, gaussMask, pixelSizeMM, '--', 0.85, ...
    sprintf('Gaussian footprint (%d%c)', gaussOrientation, char(176)));

axis(ax, 'image');
set(ax, 'YDir', 'reverse', ...
    'XLim', [min(xMM) max(xMM)], 'YLim', [min(yMM) max(yMM)], ...
    'Box', 'off');
xlabel(ax, 'Cortical distance (mm)');
ylabel(ax, 'Cortical distance (mm)');
colorbarHandle = colorbar(ax);
colorbarHandle.Label.String = colorbarLabel;
colorbarHandle.Label.FontName = 'Arial';
colorbarHandle.Label.FontSize = 10;
title(ax, titleText, 'FontWeight', 'normal', ...
    'FontName', 'Arial', 'FontSize', 12);
legendHandle = legend(ax, ...
    [targetHandle0 targetHandle90 roiHandle gaussHandle], ...
    'Location', 'southoutside', 'Orientation', 'horizontal', ...
    'FontName', 'Arial', 'FontSize', 8, 'Box', 'off');
try
    legendHandle.NumColumns = 2;
catch
end
end


function [handle0, handle90] = plotTargets(ax, target0, target90, xMM, yMM)
% Cyan/red fills only: no black or white target contours.
color0 = [0 1 1];
color90 = [1 0 0];
overlapColor = [1 0 1];
target0Only = logical(target0) & ~logical(target90);
target90Only = logical(target90) & ~logical(target0);
overlapMask = logical(target0) & logical(target90);
[nRows,nColumns] = size(target0Only);
rgbImage = zeros(nRows,nColumns,3);
alphaImage = zeros(nRows,nColumns);
for channelIndex = 1:3
    plane = rgbImage(:,:,channelIndex);
    plane(target0Only) = color0(channelIndex);
    plane(target90Only) = color90(channelIndex);
    plane(overlapMask) = overlapColor(channelIndex);
    rgbImage(:,:,channelIndex) = plane;
end
alphaImage(target0Only | target90Only) = 0.78;
alphaImage(overlapMask) = 0.90;
overlayHandle = image(ax, xMM, yMM, rgbImage);
overlayHandle.AlphaData = alphaImage;
overlayHandle.HandleVisibility = 'off';
handle0 = patch(ax, NaN, NaN, color0, ...
    'FaceAlpha', 0.78, 'EdgeColor', 'none', ...
    'DisplayName', ['0' char(176) ' designed targets']);
handle90 = patch(ax, NaN, NaN, color90, ...
    'FaceAlpha', 0.78, 'EdgeColor', 'none', ...
    'DisplayName', ['90' char(176) ' designed targets']);
end


function legendHandle = plotMask(ax, mask, pixelSizeMM, ...
        lineStyle, lineWidth, displayName)
legendHandle = plot(ax, NaN, NaN, lineStyle, ...
    'Color', [0 0 0], 'LineWidth', lineWidth, ...
    'DisplayName', displayName);
boundaries = bwboundaries(logical(mask));
for boundaryIndex = 1:numel(boundaries)
    boundary = boundaries{boundaryIndex};
    plot(ax, (boundary(:,2)-1)*pixelSizeMM, ...
        (boundary(:,1)-1)*pixelSizeMM, lineStyle, ...
        'Color', [0 0 0], 'LineWidth', lineWidth, ...
        'HandleVisibility', 'off');
end
end
