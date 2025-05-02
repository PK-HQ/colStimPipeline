%% Setup
clear;

%% Inputs
monkeyStr = {'Chip','Pepper'}; % Example for multiple monkeys
monkeyIDstr = {'28','32'}; % Example for multiple monkeys
sessionStr = {'20240131', '20250318'}; %'20250114'}; % Example for multiple sessions

pdfSuffix = 'test';

LumRuns = [7, 9]; % 0]; % Corresponding to sessions
SiteID = [1, 2]; % Corresponding to sessions
SiteIDStr = sessionStr; % Corresponding to sessions

%% Conditions
OptoLuminance = [1e-12, 5, 10, 20, 30, 50, 70, 100];
OptoMaxPower = 62.6; % mW, RM6, with ND0
OptolightGuideDi = 9; % mm
OptoProjectorArea = pi * (OptolightGuideDi / 2)^2; % mm^2
OptoLightPD = (OptoLuminance / 100) * OptoMaxPower / OptoProjectorArea; % mW/mm^2

%% Cosmetic
colormapRB = fireice;
midColormapRB = length(colormapRB) / 2;
colormapRB = colormapRB(midColormapRB - 15:midColormapRB + 16, :); % Midpoint adjustments
midColormapRB = length(colormapRB) / 2;

colormapR = [0 0 0; colormapRB(midColormapRB + 1:end, :)]; % Add zero = black
colormapB = [colormapRB(1:midColormapRB, :); 0 0 0]; % Add zero = black
colormapRB = [colormapB; 0 0 0; colormapR];

% Plot settings
set(0, 'DefaultLineLineWidth', 2);
set(0, 'DefaultFigureWindowStyle', 'docked');
set(0, 'DefaultAxesFontSize', 12);

%% Load Data
numSessions = numel(sessionStr);
DataFFT = cell(numSessions, 1);
DataTS = DataFFT;

for iSession = 1:numSessions
    fprintf('Processing session %s (Monkey %s, Site %d)...\n', ...
        sessionStr{iSession}, monkeyStr{iSession}, SiteID(iSession));
    
    DataPath = ['Y:/' monkeyStr{iSession} '/' monkeyStr{iSession} sessionStr{iSession} '/'];
    tRunNum = LumRuns(iSession);
    tRunDir = [DataPath 'run' num2str(tRunNum) '/'];
    
    % Load data
    Run = load([tRunDir 'M' monkeyIDstr{iSession} 'D' sessionStr{iSession} 'R' num2str(tRunNum) 'StabBin008.mat']);
    DataFFT{iSession} = ...
        load([tRunDir 'M' monkeyIDstr{iSession} 'D' sessionStr{iSession} 'R' num2str(tRunNum) 'StabDFFTAmpS004E023PF0400.mat']);
    DataTS{iSession} = ...
        load([tRunDir 'M' monkeyIDstr{iSession} 'D' sessionStr{iSession} 'R' num2str(tRunNum) 'TS.mat']);
end

%% Initialize Response Arrays
RespFFT = zeros(size(DataFFT{1}.DataCond, 1), ...
                size(DataFFT{1}.DataCond, 2), ...
                numSessions, numel(OptoLuminance));

for iSession = 1:numSessions
    for iLED = 1:numel(OptoLuminance)
        tData = DataFFT{iSession}.DataCond;
        RespFFT(:, :, iSession, iLED) = tData(:, :, iLED + 1) - tData(:, :, 1); % Blank-subtracted response
    end
end


%% Responses
figure;
[ha, ~] = tight_subplot(numel(SiteID) + 1, numel(OptoLuminance) + 1, .01, [.05, -.05], [.01, .07]);
numRow = numel(SiteID) + 1;
numCol = numel(OptoLuminance) + 1;
Clim = [min(RespFFT(:)), max(RespFFT(:))]*.8;

axes(ha(1));
axis off;

for iSite = 1:numel(monkeyStr)
    axes(ha(iSite * numCol + 1));
    txt = sprintf('%s', monkeyStr{iSite});
    text(0.92, 0.5, txt, 'Units', 'Normalized', ...
        'FontSize', 14, 'FontWeight', 'normal', 'HorizontalAlignment', 'center');
    axis off;
end

for iLED = 1:numel(OptoLuminance)
    axes(ha(iLED + 1));
    txt = sprintf('%.2f mW/mm^{2}\n(%g%%)', OptoLightPD(iLED), OptoLuminance(iLED));
    text(0.5, 0.01, txt, 'Units', 'Normalized', ...
        'FontSize', 14, 'FontWeight', 'normal', 'HorizontalAlignment', 'center');
    axis off;
end

for iSession = 1:numSessions
    for iLED = 1:numel(OptoLuminance)
        axes(ha(iSession * numCol + iLED + 1));
        imagesc(RespFFT(:, :, iSession, iLED), Clim);
        axis image;
        set(gca, 'xticklabel', [], 'yticklabel', []);
    end
    % Add colorbar
    origSiz = get(gca, 'Position');
    colorbar;
    set(gca, 'Position', origSiz);
end

% LABELS
[axX, hX] = suplabel(['Half-toned luminance'], 't');
[axY, hY] = suplabel(['Site ID'], 'y');
set(hX, 'FontWeight', 'bold');
set(hY, 'FontWeight', 'bold');

titlePos = get(hX, 'position');
titlePos(2) = 1.15;
set(hX, 'position', titlePos);

set(gcf, 'color', 'w');
colormap(colormapR); 

%% Responses with ROI Overlay
figure;
Dim = size(RespFFT, 1);
ROISize = 115; % ROI size
CenterResp = zeros(numSessions, numel(OptoLuminance));

for iSession = 1:numSessions
    for iLED = 1:numel(OptoLuminance)
        % Determine ROI using Gaussian filtering
        filtSigma = 30; % Gaussian filter
        filtWidth = filtSigma * 6; 
        imageFilter = fspecial('gaussian', filtWidth, filtSigma);

        x = RespFFT(:, :, iSession, iLED);
        dataFilteredX = nanconv(x, imageFilter, 'nanout');
        
        % Define the central 256x256 region (indices 129 to 384)
        centerRowStart = 129;
        centerRowEnd   = 384;
        centerColStart = 129;
        centerColEnd   = 384;
        centerRegion   = dataFilteredX(centerRowStart:centerRowEnd, centerColStart:centerColEnd);
        
        % Find the peak within the central region
        [maxVal, linearIdx] = max(centerRegion(:));
        [peakRow, peakCol] = ind2sub(size(centerRegion), linearIdx);
        
        % Convert the peak location back to the original image coordinates
        xMax = peakRow + centerRowStart - 1;
        yMax = peakCol + centerColStart - 1;
        
        % Define the ROI based on ROISize
        ROI.X = xMax - floor(ROISize / 2) + 1 : xMax + ceil(ROISize / 2);
        ROI.Y = yMax - floor(ROISize / 2) + 1 : yMax + ceil(ROISize / 2);


        %{
        x = RespFFT(:, :, iSession, iLED);
        dataFilteredX = nanconv(x, imageFilter, 'nanout');
        [xMax, yMax] = find(dataFilteredX == max(max(dataFilteredX)));
        ROI.X = xMax - floor(ROISize / 2) + 1 : xMax + ceil(ROISize / 2);
        ROI.Y = yMax - floor(ROISize / 2) + 1 : yMax + ceil(ROISize / 2);
        %}

        CenterResp(iSession, iLED) = mean(mean(RespFFT(ROI.X, ROI.Y, iSession, iLED), 1), 2) * 100;
    end
    
    subplot(1, numSessions, iSession);
    imagesc(RespFFT(:, :, iSession, end), Clim);
    colorbar;
    hold on;
    rectangle('Position', [min(ROI.Y), min(ROI.X), ROISize, ROISize], 'EdgeColor', 'k');
    axis image;
    set(gca, 'xticklabel', [], 'yticklabel', []);
    txt = sprintf('%s', monkeyStr{iSession});
    text(.5, 1.03, txt, 'Units', 'Normalized', ...
        'FontSize', 14, 'FontWeight', 'normal', 'HorizontalAlignment', 'center');
end

[axX, hX] = suplabel([{'ROI across sites', '(all responses for 100% luminance)'}], 't');
set(hX, 'FontWeight', 'bold');

set(gcf, 'color', 'w');
colormap(colormapR);

%% Naka-Rushton Fit
clc
% Modified Naka-Rushton with a baseline parameter:
% Param(1): baseline, Param(2): Rmax (amplitude above baseline),
% Param(3): C50, Param(4): exponent (n)
NakaRushton = @(Param, x) Param(1) + Param(2) .* ((x.^Param(4)) ./ ((x.^Param(4)) + (Param(3).^Param(4))));
NormCenterResp = CenterResp;

% Normalize responses
for iSession = 1:numSessions
    tResp = CenterResp(iSession, :);
    tResp = tResp - min(tResp);
    tResp = tResp ./ max(tResp);
    NormCenterResp(iSession, :) = tResp;
end

% Define preferred colormap
preferredColormap = slanCM('plasma',numSessions); % Choose a colormap (e.g., jet, parula, hot, cool)
Col = preferredColormap;%(round(linspace(1, size(preferredColormap, 1), numSessions)), :); % Assign colors based on sessions

figure;
hold on;

fitParams = zeros(numSessions, 6);
for iSession = 1:numSessions
    tLED = OptoLightPD;
    % Define initial guesses:
    % - baseline_guess: minimum response (nonzero if measured)
    % - amplitude_guess: range of response (max - min)
    % - C50: initial guess (e.g., 1)
    % - exponent: initial guess (e.g., 5)
    baseline_guess = min(CenterResp(iSession, :));
    amplitude_guess = max(CenterResp(iSession, :)) - baseline_guess;
    Param0 = [baseline_guess, amplitude_guess, 1, 3]; % Baseline, amplitude, C50, Exponent
    
    % Optionally add bounds to ensure physical realism:
    lb = [-Inf, 0, eps, 0];
    ub = [Inf, Inf, Inf, Inf];
    
    % Fit the model
    FitParam = lsqcurvefit(NakaRushton, Param0, tLED, CenterResp(iSession, :), lb, ub);
    
    % Calculate metrics like P50 relative to the baseline
    P_50 = NakaRushton(FitParam, 0.5);
    P_max = NakaRushton(FitParam, 1);
    
    fitParams(iSession, :) = [FitParam(1), FitParam(2), FitParam(3), FitParam(4), P_50, P_max];
    
    % Plot the fitted curve
    x = linspace(0, max(tLED), 100);
    plot(x, NakaRushton(FitParam, x), 'Color', Col(iSession, :), 'LineWidth', 2);
    scatter(tLED, CenterResp(iSession, :), 75, Col(iSession, :), 'filled', 'MarkerEdgeColor', 'k');
    addSkippedTicks(0,1,.1,'x');
end

% Axis customization
axis square;
xlabel('Power density (mW/mm^2)');
ylabel('ROI amplitude (\DeltaF/F%)');
title('Power-response function');
legend([monkeyStr], 'Location', 'northwest');
set(gcf, 'color', 'w');
upFontSize(24, .015);
startPos = [0.5, 3];
for iSession=1:numSessions
    % Add custom table with parameters
    if iSession == 1
        % Customize the starting position and spacing
        headers = {'Base', 'R_{max}', 'C_{50}', 'n', 'P_{50}', 'P_{max}'};
        xSpacing = 0.11; % Horizontal spacing between columns
        ySpacing = 1; % Vertical spacing between rows
    else
        headers = {'', '', '', '', '', ''}; % Blank headers for subsequent sessions
        startPos = startPos - [0, iSession*.5]; % Adjust starting position for subsequent sessions
    end
    createCustomTable2(gca, '', headers, round(fitParams(iSession, :), 1, 'decimal'), startPos, xSpacing, ySpacing);
end


%% Helper Function for Table
function createCustomTable2(ax, modelTypeStr, headers, fitParams, startPos, xSpacing, ySpacing)
    for i = 1:numel(headers)
        text(startPos(1) + (i - 1) * xSpacing, ...
             startPos(2), ...
             headers{i}, ...
             'Interpreter', 'tex', ...
             'HorizontalAlignment', 'center', ...
             'FontWeight', 'bold', ...
             'Parent', ax); % Headers
        text(startPos(1) + (i - 1) * xSpacing, ...
             startPos(2) - ySpacing, ...
             sprintf('%.1f', fitParams(i)), ...
             'HorizontalAlignment', 'center', ...
             'Parent', ax); % Values
    end
end
