%% Setup


%% Inputs
monkeyStr = 'Chip';
monkeyIDstr = '28';
sessionStr = {'20230812'}; % Example for multiple sessions
pdfSuffix = 'test';

LumRuns = [5, 3]; % Corresponding to sessions
SiteID = [1, 1]; % Corresponding to sessions
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
    fprintf('Processing session %s (Site %d)...\n', sessionStr{iSession}, SiteID(iSession));
    
    DataPath = ['Y:/' monkeyStr '/' monkeyStr sessionStr{iSession} '/'];
    tRunNum = LumRuns(iSession);
    tRunDir = [DataPath 'run' num2str(tRunNum) '/'];
    
    % Load data
    Run = load([tRunDir 'M' monkeyIDstr 'D' sessionStr{iSession} 'R' num2str(tRunNum) 'StabBin008.mat']);
    DataFFT{iSession} = ...
        load([tRunDir 'M' monkeyIDstr 'D' sessionStr{iSession} 'R' num2str(tRunNum) 'StabDFFTAmpS004E023PF0400.mat']);
    DataTS{iSession} = ...
        load([tRunDir 'M' monkeyIDstr 'D' sessionStr{iSession} 'R' num2str(tRunNum) 'TS.mat']);
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

%% Naka-Rushton Fit
clc

Param0 = [0, ... % RMax
          0, ... % C50 (Asymptotic maximum response amplitude)
          1, ... % n
          ];
NakaRushton = @(Param, x) Param(1) .* ((x .^ Param(3)) ./ ((x .^ Param(3)) + (Param(2) .^ Param(3))));

CenterResp = zeros(numSessions, numel(OptoLuminance));
ROISize = 115; % ROI size

for iSession = 1:numSessions
    for iLED = 1:numel(OptoLuminance)
        % Determine ROI based on max response
        filtSigma = 30; % Gaussian filter
        filtWidth = filtSigma * 6; 
        imageFilter = fspecial('gaussian', filtWidth, filtSigma);

        x = RespFFT(:, :, iSession, iLED);
        dataFilteredX = nanconv(x, imageFilter, 'nanout');
        [xMax, yMax] = find(dataFilteredX == max(max(dataFilteredX)));
        ROI.X = xMax - floor(ROISize / 2) + 1 : xMax + ceil(ROISize / 2);
        ROI.Y = yMax - floor(ROISize / 2) + 1 : yMax + ceil(ROISize / 2);

        % Calculate center response
        CenterResp(iSession, iLED) = mean(mean(RespFFT(ROI.X, ROI.Y, iSession, iLED), 1), 2)*100;
    end
end

FitParam = zeros(numSessions, numel(Param0));
NormCenterResp = CenterResp;

for iSession = 1:numSessions
    tResp = CenterResp(iSession, :);
    tResp = tResp - min(tResp);
    tResp = tResp ./ max(tResp);
    NormCenterResp(iSession, :) = tResp;
end

%% Plot Power-Response Function with Naka-Rushton Parameters
annot = {};
Col = colormap(colormapR);
Col = Col(round(linspace(1, size(Col, 1), numSessions)), :);

figure;
hold on;
ax = gca; % Get current axes

% Headers for Naka-Rushton parameters
headers1 = {'R_{max}', 'C_{50}', 'n', 'P_{50}', 'P_{max}'};
headers2 = {'', '', '', '', ''};
for iSession = 1:numSessions
    tLED = OptoLightPD;

    % Fit the Naka-Rushton function
    FitParam = lsqcurvefit(NakaRushton, Param0, tLED, CenterResp(iSession, :));
    
    % Calculate P_{50} and P_{max}
    P_50 = NakaRushton(FitParam, .5); % Response at 50% power
    P_max = NakaRushton(FitParam, 1); % Response at 100% power

    % Collect parameters for display
    fitParams(iSession,:) = [FitParam(1), FitParam(2), FitParam(3), P_50, P_max];

    % Prepare annotations
    txt = sprintf('Session %s: R_m_a_x = %.2f, C_5_0 = %.2f, n = %.2f', ...
                  sessionStr{iSession}, FitParam(1), FitParam(2), FitParam(3));
    annot{end+1} = txt;

    % Plot the fit
    x = linspace(0, max(tLED));
    plot(x, NakaRushton(FitParam, x), ...
        'Color', Col(iSession, :), 'LineStyle', '-', 'LineWidth', 4);

    % Scatter the data points
    scatter(tLED, CenterResp(iSession, :), 75, Col(iSession, :), 'filled', 'HandleVisibility','off');
    upFontSize
end

% Customize the starting position and spacing
startPos = [.5, 5]; % Starting position in data coordinates
xSpacing =.1; % Horizontal spacing between columns
ySpacing = 1; % Vertical spacing between rows
createCustomTable2(ax, '', headers1, fitParams(1,:), startPos, xSpacing, ySpacing); hold on;
createCustomTable2(ax, '', headers2, fitParams(2,:), startPos-[0 1], xSpacing, ySpacing); hold on;

axis square;
xlabel('Stimulation power density (mW/mm^2)');
ylabel('Average ROI response amplitude (\DeltaF/F%)');
title('Power-response function (Naka-Rushton fit)');
legend(sessionStr, 'Location', 'northwest');
set(gcf, 'color', 'w');

% Save the plot
%export_fig(pdfFilename, '-append', '-pdf', '-nocrop');

function createCustomTable2(ax, modelTypeStr, headers, fitParams, startPos, xSpacing, ySpacing)
    % Create a custom table on axes using text objects with subscripts.
    %
    % Parameters:
    %   ax - The axes handle where the table will be drawn.
    %   headers - A cell array of header names with subscripts.
    %   fitParams - A 1xN numeric array containing the values.
    %   startPos - Starting position of the table ([x, y]) in normalized coordinates.
    %   xSpacing - Horizontal spacing between columns.
    %   ySpacing - Vertical spacing between rows.

    % Plot headers
    for i = 1:numel(headers)
        text(startPos(1) + (i - 1) * xSpacing, ...
             startPos(2), ...
             headers{i}, ...
             'Interpreter', 'tex', ...
             'HorizontalAlignment', 'center', ...
             'FontWeight', 'bold', ...
             'Parent', ax); % Use the specified axes
    end

    % Plot values
    for i = 1:numel(fitParams)
        formattedValue = sprintf('%.2f', fitParams(i)); % Format with 2 decimal places
        text(startPos(1) + (i - 1) * xSpacing, ...
             startPos(2) - ySpacing, ...
             formattedValue, ...
             'HorizontalAlignment', 'center', ...
             'Parent', ax); % Use the specified axes
    end
end
