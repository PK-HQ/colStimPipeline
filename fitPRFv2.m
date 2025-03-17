%% Setup
clear;

%% Inputs
monkeyStr = {'Chip', 'Chip','Pepper'}; % Example for multiple monkeys
monkeyIDstr = {'28', '28','32'}; % Example for multiple monkeys
sessionStr = {'20240112', '20240131', '20250114'}; % Example for multiple sessions

pdfSuffix = 'test';

LumRuns = [10,7,0]; % Corresponding to sessions
SiteID = [1, 1, 1]; % Corresponding to sessions
SiteIDStr = sessionStr; % Corresponding to sessions

%% Conditions
OptoLuminance = [5, 10, 20, 30, 50, 70, 100];
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
Clim = [min(RespFFT(:)), max(RespFFT(:))];

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
        [xMax, yMax] = find(dataFilteredX == max(max(dataFilteredX)));
        ROI.X = xMax - floor(ROISize / 2) + 1 : xMax + ceil(ROISize / 2);
        ROI.Y = yMax - floor(ROISize / 2) + 1 : yMax + ceil(ROISize / 2);

        CenterResp(iSession, iLED) = mean(mean(RespFFT(ROI.X, ROI.Y, iSession, iLED), 1), 2) * 100;
    end
    
    subplot(1, numSessions, iSession);
    imagesc(RespFFT(:, :, iSession, end), [min(RespFFT(:)), max(RespFFT(:))]);
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
Param0 = [0, 0, 1]; % Initial guesses for RMax, C50, and n
NakaRushton = @(Param, x) Param(1) .* ((x .^ Param(3)) ./ ((x .^ Param(3)) + (Param(2) .^ Param(3))));

FitParam = zeros(numSessions, numel(Param0));
NormCenterResp = CenterResp;

% Normalize responses
for iSession = 1:numSessions
    tResp = CenterResp(iSession, :);
    tResp = tResp - min(tResp);
    tResp = tResp ./ max(tResp);
    NormCenterResp(iSession, :) = tResp;
end

% Define preferred colormap
preferredColormap = slanCM('guppy',3); % Choose a colormap (e.g., jet, parula, hot, cool)
Col = preferredColormap(round(linspace(1, size(preferredColormap, 1), numSessions)), :); % Assign colors based on sessions

figure;
hold on;

fitParams = zeros(numSessions, 5);

for iSession = 1:numSessions
    tLED = OptoLightPD;

    % Fit the Naka-Rushton function
    FitParam = lsqcurvefit(NakaRushton, Param0, tLED, CenterResp(iSession, :));
    
    % Calculate P_{50} and P_{max}
    P_50 = NakaRushton(FitParam, 0.5); % Response at 50% power
    P_max = NakaRushton(FitParam, 1); % Response at 100% power

    % Collect parameters for display
    fitParams(iSession, :) = [FitParam(1), FitParam(2), FitParam(3), P_50, P_max];

    % Plot the fit
    x = linspace(0, max(tLED), 100);
    plot(x, NakaRushton(FitParam, x), ...
        'Color', Col(iSession, :), 'LineStyle', '-', 'LineWidth', 2); % Use colormap for line color

    % Scatter the data points
    scatter(tLED, CenterResp(iSession, :), 75, Col(iSession, :), 'filled', 'MarkerEdgeColor', 'k','HandleVisibility','off'); % Use colormap for scatter color


end
% Axis customization
axis square;
xlabel('Stimulation power density (mW/mm^2)');
ylabel('Average ROI response amplitude (\DeltaF/F%)');
title('Power-response function (Naka-Rushton fit)');
legend([monkeyStr], 'Location', 'northwest');
set(gcf, 'color', 'w');
upFontSize
for iSession=1:numSessions
    % Add custom table with parameters
    if iSession == 1
        % Customize the starting position and spacing
        headers = {'R_{max}', 'C_{50}', 'n', 'P_{50}', 'P_{max}'};
        startPos = [0.5, 3]; % Starting position in data coordinates
        xSpacing = 0.1; % Horizontal spacing between columns
        ySpacing = .5; % Vertical spacing between rows
    else
        headers = {'', '', '', '', ''}; % Blank headers for subsequent sessions
        startPos = [0.5, 3.5] - [0, iSession / 2]; % Adjust starting position for subsequent sessions
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
