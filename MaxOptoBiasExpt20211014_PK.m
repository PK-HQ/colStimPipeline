clear;clc
close all
directory='Y:/';
cd(directory)

%% Cosmetic parameters
% Colour variables stolen from Matt
grey       = [.5 .5 .5 .5];
black      = [0 0 0];
blackSolid = black(1:3)+[ .3 .3 .3];
blue       = [.0 .3 .7 .3];
blueSolid  = blue(1:3)+[.3 .3 .3];
green      = [0 .7 .2 .5];
gold       = [.8 .6 0];
red        = [.8 0.1 .1 .5];
burnt      = [.9 .6 0 .5];
purple     = [.4 0 .6 .3];
pink       = [.9 .5 .65 .3];
sea        = [.3 .8 .7 .3];
lightBlue  = blue(1:3)+[.3 .3 .3];
lightRed   = red(1:3)+[.2 .3 .3];
lightGold  = gold(1:3)+[.1 .3 .3];
lightGreen = green(1:3) + [.3 .3 .3];

colours = {red blue green gold purple sea burnt pink lightRed lightBlue ...
  lightGreen lightGold red blue green gold purple sea burnt pink lightRed ...
  lightBlue lightGreen lightGold};

ImperialRed   = [230  57  70]./255;
TuftsBlue     = [ 48 131 220]./255;
OceanGreen    = [115 186 155]./255;
CarrotOrange  = [244 157  55]./255;
AfricanViolet = [180 126 179]./255;
YellowGreen   = [159 211  86]./255;
SafetyYellow  = [235 211  49]./255;
DarkSlateBlue = [ 82  72 156]./255;
Rhythm        = [119 125 167]./255;
WarmBlack     = [  3  71  72]./255;

Col = {...
  ImperialRed...
  TuftsBlue...
  SafetyYellow...
  CarrotOrange...
  OceanGreen...
  AfricanViolet...
  YellowGreen...
  DarkSlateBlue...
  Rhythm...
  WarmBlack...
  };


set(0,'DefaultLineLineWidth',3)
set(0,'DefaultFigureWindowStyle','default')
set(0,'DefaultAxesFontSize',12)

%% Parameters
LEDPercents = [5;...
               10;...
               25;...
               50;...
               75;...
               100];

%Filename params
monke = 'Max';
monkeID = '27';
session = '20211014';

%filename identifiers
runPvR = 1;
runLgGrt = 2;
runVsGss = 3;

%Stimuli ARs
vsGssARs = [-1.5, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 1.5];

%% Power/response function run

%   fprintf('Site %g!/n',Sites(iSite))

%datapath = ['/data/' monke '\' monke session '/run' num2str(runPvR) '\'];
datapathSess = [monke '\' monke session '\'];
datapathRun = [monke '\' monke session '\run' num2str(runPvR) '\'];

fftPvR = ...
  load([grabFilename(datapathRun, '*FFTAmp*.mat')]);
tsPvR  = ...
  load([grabFilename(datapathRun, '*TS.mat')]);

respPvR = fftPvR.DataCond;

respPvR = respPvR(:, :, 2:end) - respPvR(:, :, 1);

px2mm = tsPvR.TS.Header.Imaging.SizePxl;
if (isfield(fftPvR, 'BinSize'))
  px2mm = px2mm * fftPvR.BinSize;
end

% Filter stuff
for i = 1:size(respPvR,3)
  respPvR(:,:,i) = FilterFermi2DNoFig(respPvR(:,:,i), 0, 1, px2mm);
end

% RespFFT(RespFFT < 0) = 0;
%figure;imagesc(respPvR(:, :, end));axis image
% test = RespFFT(:, :, end);
% figure;histogram(test(:))

%% Get maximum optostim condition

% Get blank-subtracted 100% uniform opto response
OptoResponse = respPvR(:,:,end);
cprintf(sea(1:3),...
  'Code check: OptoResponse size is %d X %d X %d; should be 512 X 512 X 1./n',...
  size(OptoResponse,1),size(OptoResponse,2),size(OptoResponse,3))

% Responses across delta-AR and orientation
figure('name','max-uniform-optostim')

set(gcf,'position',get(0,'screensize'))

imagesc(OptoResponse)
colorbar
axis square
set(gca,'xticklabel',[],'yticklabel',[])
title('max uniform optostim response')

figure('name','max-response-2d')
set(gcf,'position',get(0,'screensize'))
plot(1:size(OptoResponse,1),OptoResponse(size(OptoResponse,1)/2-1,:));
xlabel('horizontal space')
ylabel('{/Delta}F/F')
title('2D slice of max response')

OptoResponse = OptoResponse-min(OptoResponse(:));
OptoResponse = OptoResponse./max(OptoResponse(:));

% Pick desired max amplitude
MaxAmp = 0.1;
Mask = OptoResponse >= MaxAmp;

%% ROIs

% Want to bin the image into some number of regions to calculate
% nonlinearity separately in

ROISize = 4;
Dim = size(respPvR,1);
nBin = Dim / ROISize;
assert(mod(Dim, nBin) == 0);
ROIs = zeros(nBin, ROISize);

for iBin = 1:nBin
  
  tROI = 1:ROISize;
  tROI = tROI + ((iBin - 1) * ROISize);
  ROIs(iBin, :) = tROI;
  
end

MaskBinned = zeros(nBin, nBin);

for xBin = 1:nBin
  
  for yBin = 1:nBin
    
    xROI = ROIs(xBin, :);
    yROI = ROIs(yBin, :);
    
    tBinMaxCond = mean(mean(OptoResponse(yROI, xROI), 1), 2);
    
    if (tBinMaxCond >= MaxAmp)
      
      MaskBinned(yBin, xBin) = 1;
      
    else
      
      MaskBinned(yBin, xBin) = 0;
      
    end
    
  end
  
end

MaskBinned = logical(MaskBinned);

CenterResp = zeros(nBin, nBin, numel(LEDPercents));

RespMask = respPvR .* Mask;

for iLED = 1:numel(LEDPercents)
  
  for xBin = 1:nBin
    
    for yBin = 1:nBin
    
    xROI = ROIs(xBin, :);
    yROI = ROIs(yBin, :);
    if (MaskBinned(yBin, xBin))
      
      CenterResp(yBin, xBin, iLED) = ...
        mean(mean(RespMask(yROI, xROI, iLED), 1), 2);
      
    else
      
      CenterResp(yBin, xBin, iLED) = 0;
      
    end
    
    end
  
  end
  
end

%% Fitting

Param0 = [1, ...    % RMax
          0, ...    % P50
          1, ...    % n
          0];       % DC addition
lb =     [0, ...    % RMax
          0, ...    % P50
          1, ...    % n
          0];       % DC addition
ub =     [2, ...   % RMax
          100, ...  % P50
          10, ...   % n
          1];       % DC addition
        
NakaRushton = @(Param, x) Param(1) .* ((x .^ Param(3)) ./ ((x .^ Param(3)) ...
                          + (Param(2) .^ Param(3)))) + Param(4);
PcLED = LEDPercents ./ 100;
optOptm = optimset('LargeScale', 'On', 'Display', 'Off',...
                   'TolFun', 1e-8, 'TolX', 1e-8,...
                   'MaxFunEvals', 1e6, 'MaxIter', 1e5);
                        
% Fit to the entire image first, hold n relative stable after
MeanResp = squeeze(mean(mean(respPvR, 1), 2));
MeanFit = (lsqcurvefit(NakaRushton, Param0, PcLED, MeanResp, lb, ub, optOptm));

% Display whole-image fit
figure('Name', 'Whole-Field Naka-Rushton')
hold on
scatter(PcLED, MeanResp, 150, Col{1}, 'LineWidth',2);
plot(linspace(0, max(PcLED)), NakaRushton(MeanFit, linspace(0, max(PcLED))))

% Prompt user for linear fit or not
if (~exist('bLinear', 'var'))
  answerLinear = questdlg('Use a linear fit?', 'Fit', 'Yes', 'No', 'No');
end

switch answerLinear
  case 'Yes'
    bLinear = true;
  case 'No'
    bLinear = false;
end

if (bLinear)
  
  P0 = [0, ...  % Slope
        0];     % Intercept
  lb = [0, ...
        0];
  ub = [100, ...
        0];
      
  FitParams = zeros(nBin, nBin, numel(P0));
  NormCenterResp = CenterResp - min(CenterResp(:));
  NormCenterResp = NormCenterResp ./ max(NormCenterResp(:));
  objectiveFunc = @(param, x) param(1) .* x + param(2);

  disp('Fitting linear parameters over space...')
  for xBin = 1:nBin

    for yBin = 1:nBin

      FitParams(yBin, xBin, :) = lsqcurvefit(objectiveFunc, P0, ...
         PcLED, squeeze(NormCenterResp(yBin, xBin, :)), lb, ub, optOptm);

    end

  end
  disp('Done!')
  
  RepFitParams = zeros(Dim, Dim, size(FitParams, 3));

  for i = 1:size(FitParams, 3)

    RepFitParams(:, :, i) = kron(FitParams(:, :, i), ones(ROISize));

  end
  
  MaskBinned = true(size(MaskBinned));
  
end

if (~bLinear)

  Param0(3) = MeanFit(3);
  lb(3)     = MeanFit(3) - 0.5;
  ub(3)     = MeanFit(3) + 0.5;

  FitParams = zeros(nBin, nBin, numel(Param0));
  NormCenterResp = CenterResp - min(CenterResp(:));
  NormCenterResp = NormCenterResp ./ max(NormCenterResp(:));

  disp('Fitting Naka-Rushton parameters over space...')
  for xBin = 1:nBin

    for yBin = 1:nBin

      FitParams(yBin, xBin, :) = lsqcurvefit(NakaRushton, Param0, PcLED, ...
        squeeze(CenterResp(yBin, xBin, :)), lb, ub, optOptm);

    end

  end
  disp('Done!')

  RepFitParams = zeros(Dim, Dim, size(FitParams, 3));

  for i = 1:size(FitParams, 3)

    RepFitParams(:, :, i) = kron(FitParams(:, :, i), ones(ROISize));

  end

  % RepFitParams = RepFitParams .* Mask;

  RepMaskBinned = kron(MaskBinned, ones(ROISize));

end

% Display fit variables over space

if (bLinear)
  params = {'slope', 'intercept'};
elseif (~bLinear)
  params = {'R_m_a_x', 'P_5_0', 'n'};
end

clims = zeros(3, 2);
for i = 1:numel(params)
  tX = FitParams(:, :, i);
  tXMasked = FitParams(:, :, i) .* MaskBinned;
  clims(i, :) = [min(tX(:)) max(tXMasked(:))];
end
  
figure('name', 'params-over-space')

for i = 1:numel(params)
  subplot(1, numel(params), i)
  imagesc(FitParams(:, :, i) .* MaskBinned, clims(1, :))
  axis image
  set(gca, 'xticklabel', [], 'yticklabel', [])
  title(sprintf([params{i} ' over space']))
  colorbar
end

%%

% W = sqrt(FilterFermi2DNoFig(OptoResponse,0,1.0,.1224));
% W = sqrt(OptoResponse);
W = OptoResponse;
W(Mask) = MaxAmp ./ W(Mask);
W(~Mask) = 0;

figure('name','optostim correction')

set(gcf,'position',get(0,'screensize'))

imagesc(W)
colorbar
axis square
set(gca,'xticklabel',[],'yticklabel',[])
title('optostim correction weights')

Proj_PixelSize = 10.8/1000;  % mm per px
fLensProj      = 50;
fBottomLens    = 35;
ProjPx2mm      = Proj_PixelSize/(fBottomLens/fLensProj); % mm per px
CamPx2mm       = tsPvR.TS.Header.Imaging.SizePxl; % mm per px
if (isfield(fftPvR, 'BinSize'))
  CamPx2mm = CamPx2mm * fftPvR.BinSize;
end


% Projector BMP dimensions
PSizeX = 608;  % pixel
PSizeY = 684;  % pixel
PHalfY = PSizeY/2;

% Factor to rescale camera images by to match projector pixel resolution
Cam2Proj = CamPx2mm/ProjPx2mm; 

%%
% Need to find common center between projector stimuli and camera image
ExptDir = [monke '\' monke session '\'];

% Use small dot projection to find projection center
DotImage = imread([ExptDir 'centering.tif']);

%%
figure('name','projector-centering')

set(gcf,'position',get(0,'screensize'))

imagesc(DotImage)
colormap gray
axis square
set(gca,'xticklabel',[],'yticklabel',[])
title('small dot stimulus projection')

cprintf(green(1:3),'Click dot center./n')
ProjCtr = ginput(1);

% %%
% DotImage = zeros(64,64);
% ProjCtr = [size(DotImage,1)/2 size(DotImage,2)/2];

% The non-run camera images have dimensions of 2048 X 2048 px
% Change ProjCtr to offset from center
ProjCtr(1) = ProjCtr(1)-size(DotImage,1)/2; % neg if left of center
ProjCtr(2) = size(DotImage,2)/2-ProjCtr(2); % neg if below center
% 
% % ProjCtr = [0,0];
% 
ProjCtr = round(ProjCtr*(size(OptoResponse,1)/2048)*Cam2Proj);
cprintf(pink(1:3),...
  'Projector offset is %d, %d px in projector pixel resolution/n',...
  ProjCtr(1),ProjCtr(2))

% End weighting matrix (W) has to be offset and 608 X 342 px @ 0.0108 mm/px
%   Input: 512 X 512 px @ 0.0153 mm/px

%   Find the section of W actually covered by projector
%   Remake stimuli this time multiplied by W

%   Scale W to have same pixel resolution (px2mm) as projector
WRescale = imresize(W,Cam2Proj,'method','bilinear'); 

% W is taller and narrower than projector coverage (P); find indices
% of W X relative to P and P Y relative to W assuming centred projection

% Indices corresponding to W's centred placement within P across X
PIndWX = ceil(PSizeX/2)-floor(size(WRescale,2)/2)+1: ...
         ceil(PSizeX/2)+ceil(size(WRescale,2)/2);
% Indices corresponding to P's centred placement within W across Y
WIndPY = ceil(size(WRescale,1)/2)-floor(PHalfY/2)+1: ...
         ceil(size(WRescale,1)/2)+ceil(PHalfY/2);

% Apply offset
PIndWX = round(PIndWX+ProjCtr(1));
WIndPY = round(WIndPY-ProjCtr(2)); % because up-down are switched for matrices

TestMat = zeros(PHalfY,PSizeX);
cprintf(sea(1:3),...
  ['Code check: P subset size is %d X %d and W subset size is %d X %d;/n  '...
  'the two should be the same size./n'],...
  size(TestMat(:,PIndWX),2),size(TestMat(:,PIndWX),1),...
  size(WRescale(WIndPY,:),2),size(WRescale(WIndPY,:),1))

%% Preparing to create projector BMPs 

% Things necessary to create projector BMPs
%   * Desired response
%   * Nonlinearity
%   * Opsin profile
%   * SIRF for deconvolution

% Load visual runs: large gratings
datapathRun = [datapathSess 'run' num2str(runLgGrt) '\'];
fftLgGrt = ...
  load([datapathRun 'M' monkeID 'D' session 'R' num2str(runLgGrt) ...
  'StabBin008FFTAmpS004E023PF0400.mat']);
tsLgGrt  = ...
  load([datapathRun 'M' monkeID 'D' session 'R' num2str(runLgGrt) 'TS.mat']);
nBlank = tsLgGrt.TS.Header.NumBlankCond;
respLgGrt = fftLgGrt.DataCond;
respLgGrt = respLgGrt(:,:,nBlank+1:end)-mean(respLgGrt(:,:,1:nBlank),3);
px2mm = tsLgGrt.TS.Header.Imaging.SizePxl;
if (isfield(fftLgGrt, 'BinSize'))
  px2mm = px2mm * fftLgGrt.BinSize;
end
% Filter stuff
for i = 1:size(respLgGrt,3)
  respLgGrt(:,:,i) = FilterFermi2DNoFig(respLgGrt(:,:,i), 0, 1, px2mm);
end

meanRespLgGrt = mean(respLgGrt,3);

% Load visual runs: Gaussians
datapathRun = [datapathSess 'run' num2str(runVsGss) '\'];
fftVsGss = ...
  load([datapathRun 'M' monkeID 'D' session 'R' num2str(runVsGss) ...
  'StabBin008FFTAmpS004E023PF0400.mat']);
tsVsGss  = ...
  load([datapathRun 'M' monkeID 'D' session 'R' num2str(runVsGss) 'TS.mat']);
nBlank = tsVsGss.TS.Header.NumBlankCond;
respVsGss = fftVsGss.DataCond;
respVsGss = respVsGss(:,:,nBlank+1:end)-mean(respVsGss(:,:,1:nBlank),3);
px2mm = tsVsGss.TS.Header.Imaging.SizePxl;
if (isfield(fftVsGss, 'BinSize'))
  px2mm = px2mm * fftVsGss.BinSize;
end
% Filter stuff
for i = 1:size(respVsGss,3)
  respVsGss(:,:,i) = FilterFermi2DNoFig(respVsGss(:,:,i), 0, 1, px2mm);
end

% Gaussian runs need to be sorted in order of dAR
% Since it's always circular first then alternating hori/vert in increasing
% dAR it should be easy to do just hardcoded
indsH = [2:2:size(respVsGss,3)];
indsV = [3:2:size(respVsGss,3)];
indCirc = ceil(size(respVsGss,3)/2);
sortInds = [flip(indsH) indCirc indsV];

respVsGss = respVsGss(:,:,sortInds);

% Fit Gaussians

% Fit to mean
respMean = mean(respVsGss, 3);

expression = meanRespLgGrt;
expression = expression - min(expression(:));
expression = expression ./ max(expression(:));

% Define coordinate space
coord.X = 1:size(respMean, 1);
coord.Y = 1:size(respMean, 2);
coord.expression = expression;
% Define initial 2D Gaussian parameters
PI  = [0,      ...                      % amplitude
       0, ...                           % center X coordinate
       0, ...                           % center Y coordinate
       0,                           ... % DC component
       45,                          ... % orientation
       size(respMean, 1) / 10,      ... % sigma 1
       size(respMean, 1) / 10];         % sigma 2
% Define lower bound for parameter fitting
PLB = [0,   ...
       1,   ...
       1,   ...
       0,  ...
       0, ...
       0,   ...
       0];
% Define upper bound for parameter fitting
PUB = [1,             ...
       size(respMean, 1), ...
       size(respMean, 2), ...
       1,             ...
       90,            ...
       size(respMean, 1),...
       size(respMean, 1)];    
% Define optimization options to pass to fitting function
optOptm = optimset('LargeScale', 'On',...
                   'TolFun', 1e-8, 'TolX', 1e-8,...
                   'MaxFunEvals', 1e6, 'MaxIter', 1e5);
fitParamsMean = lsqcurvefit('FuncWoNGaussian2DExpression', ...
                            PI, coord, respMean, PLB, PUB, optOptm);

% Get fitted Gaussian center coords and orientation
% fitPeak = 1;
fitCenX = fitParamsMean(2);
fitCenY = fitParamsMean(3);
fitDC   = fitParamsMean(4);
fitOri  = fitParamsMean(5);

% Hold these parameters constant across subsequent fittings
% [PI(1), PLB(1), PUB(1)] = deal(fitPeak);
[PI(2), PLB(2), PUB(2)] = deal(fitCenX);
[PI(3), PLB(3), PUB(3)] = deal(fitCenY);
[PI(4), PLB(4), PUB(4)] = deal(fitDC);
[PI(5), PLB(5), PUB(5)] = deal(fitOri);

paramVsGss = zeros(numel(fitParamsMean), size(respVsGss, 3));

for i = 1:size(respVsGss, 3)
  
  paramVsGss(:, i) = ...
    lsqcurvefit('FuncWoNGaussian2DArea', ...
    PI, coord, respVsGss(:, :, i), PLB, PUB, optOptm);
  
end

dARs = vsGssARs;

% Visual response fits vs projected response fits
param = paramVsGss;
% Normalize sigmas by across-condition means
param(6, :) = param(6, :) ./ mean(param(6, :), 2);
param(7, :) = param(7, :) ./ mean(param(7, :), 2);

sigMa = find(fitParamsMean == min([fitParamsMean(6), fitParamsMean(7)]));
sigMi = 13 - sigMa;

fitdARsNorm = zeros(size(paramVsGss, 2), 1);

for i = 1:size(param, 2)
  
  sig1 = param(sigMa, i);
  sig2 = param(sigMi, i);
  
  dar = sig1 / sig2;
  
  if (sig2 > sig1)
    
    dar = 1 / dar;
    dar = -(dar - 1);
    
  else
    
    dar = dar - 1;
    
  end
  
  fitdARsNorm(i) = dar;
  
end

% dar = -1.*dar;

[funcOpto, S] = polyfit(dARs, fitdARsNorm', 1);
x1 = linspace(min(dARs), max(dARs));
yOpto = polyval(funcOpto, x1);

r2Opto = 1 - (S.normr / norm(fitdARsNorm - mean(fitdARsNorm))) ^ 2;

% Load original VSD fits
% load('/eslab/public/users/shun/code/MaxVSDAnly20210830Params.mat')

% yOrig = polyval(arFunc, x1);


figure('Name', 'Visual vs Response ARs')
hold on
plot(x1, yOpto, 'color', Col{1})
scatter(dARs, fitdARsNorm, 150, Col{1}, 'filled', 'LineWidth', 2, ...
  'handlevisibility', 'off')
% plot(x1, yOrig, 'color', Col{2})
% scatter(dARs, corticalARs, 150, Col{2}, 'filled', 'LineWidth', 2, ...
%   'handlevisibility', 'off')
axis square
xlabel('Visual Stimulus dAR')
ylabel('Norm. Response Fit dAR')
title('dAR of Evoked Response Fits')
% legend(sprintf('Opto-evoked (R^2 = %g)', round(r2Opto, 4)), ...
%        sprintf('Visual-evoked (R^2 = %g)', round(r2, 4)), ...
%        'location', 'southeast')


% Create Gaussians based on fit lines to responses
desired = 0.5;
sigMajMu = fitParamsMean(sigMa);
sigMinMu = fitParamsMean(sigMi);
% AR of mean-condition fit
ARMu = sigMajMu/sigMinMu;

vertAR = polyval(funcOpto,desired)+1;
vertAR = vertAR*ARMu;
vFitSigMaj = sigMajMu*sqrt(vertAR);
vFitSigMin = sigMinMu/sqrt(vertAR);
vParams = fitParamsMean;
vParams(sigMa) = vFitSigMaj;
vParams(sigMi) = vFitSigMin;

% Get horizontal
horiAR = 1/(polyval(funcOpto,-desired)+1);
horiAR = horiAR*ARMu;
hFitSigMaj = sigMajMu/sqrt(horiAR);
hFitSigMin = sigMinMu*sqrt(horiAR);
hParams = fitParamsMean;
hParams(sigMa) = hFitSigMaj;
hParams(sigMi) = hFitSigMin;

figure
subplot(1,2,1)
  imagesc(FuncWoNGaussian2D(hParams, coord))
  axis image
  set(gca,'xticklabel',[],'yticklabel',[])
subplot(1,2,2)
  imagesc(FuncWoNGaussian2D(vParams, coord))
  axis image
  set(gca,'xticklabel',[],'yticklabel',[])

%% Create desired responses
% Get desired center coords

if (~exist('GssCtr','var'))
  figure('name','gaussian-centering')

  set(gcf,'position',get(0,'screensize'))

  imagesc(OptoResponse)
  colormap gray
  axis square
  set(gca,'xticklabel',[],'yticklabel',[])
  title('max uniform optostim')

  cprintf(green(1:3),'Select center coords for Gaussians./n')
  GssCtr = ginput(1);
end
  

% %   Gaussians of same two shapes but varying amplitudes
% ampScale = [1.00];
% sd = size(RespFFT,1) / 6; % size 1 Gss sigma = 1/6 of the imaging fov
% sz = 0.5; % 6 sigmas of the desired Gss is sz the imaging fov across
% ars = [1.1, 1.3, 1.5];
% orts = [0 90];
% coord.X = 1:size(RespFFT,2);
% coord.Y = 1:size(RespFFT,1);
% desiredResps = ...
%   zeros(size(RespFFT,1),size(RespFFT,2),numel(ampScale),numel(ars),numel(orts));
% for iAmp = 1:numel(ampScale)
%   for iAR = 1:numel(ars)
%     for iOrt = 1:numel(orts)
%       tParam = [ampScale(iAmp),...
%                 GssCtr(1),...
%                 GssCtr(2),...
%                 0,...
%                 orts(iOrt),...
%                 sd*sz/sqrt(ars(iAR)),...
%                 sd*sz*sqrt(ars(iAR))];              
%       desiredResps(:,:,iAmp,iAR,iOrt) = FuncWoNGaussian2D(tParam,coord);
%     end
%   end
% end

desiredParams = [hParams;vParams];

% Load fit params
desiredResps = zeros(size(respPvR, 1), ...
                     size(respPvR, 2), ...
                     size(desiredParams, 1),1,1);
coord.X = 1:size(respPvR,2);
coord.Y = 1:size(respPvR,1);
for i = 1:size(desiredParams, 1)
  tParam = desiredParams(i, :);
  tParam(1) = 1; % Amplitude
  tParam(2) = GssCtr(1);
  tParam(3) = GssCtr(2);
  tParam(4) = 0; % Base
  tParam(6:7) = 0.75 .* tParam(6:7); 
  desiredResps(:, :, i) = FuncWoNGaussian2D(tParam, coord);
end
ampScale = 1;
ars = unique(abs(vsGssARs));
orts = [0, 90];

% % Inverse Naka-Rushton nonlinearity
% INR = @(Param, R) ...
%   nthroot((R .* (Param(2) .^ Param(3))) ./ (Param(1) - R), ...
%   repmat(Param(3), size(R, 1), size(R, 2)));

% Construct spatial impulse response function
x = (-floor(size(respPvR, 2)/2)+1:ceil(size(respPvR,2)/2)).*CamPx2mm;
y = (-floor(size(respPvR, 1)/2)+1:ceil(size(respPvR,1)/2)).*CamPx2mm;
[X,Y] = meshgrid(x,y);
tau = .4; % cycles per mm
sirf = exp(-sqrt(X.^2+Y.^2)./tau);

volume = trapz(y, trapz(x, sirf, 2), 1);
sirf = sirf ./ volume;

%%

% % desiredResps = zeros(Dim, Dim);
% ampScale = 1;
% ars = 1;
% orts = 1;
% sz = 1;
% % 
% % desiredResps(end-Dim/2+1:end-Dim/4, Dim/4+1:Dim/2) = 1;
% % % figure;imagesc(desiredResps);axis image
% % 
% desiredResps = ones(Dim, Dim);
% % % desiredResps(end/2+1:end, :) = 0;

%% Create projector BMPs to recreate desired response

% desiredResps = ones(size(desiredResps, 1), size(desiredResps, 2));

if (~exist('pathname','var'))
  pathname = uigetdir('.','Select a Directory to Save BMP Files');
  drawnow;
end
if pathname==false
  disp('No directory is selected!');
  return;
end

bIndivNorm = false;

dar = [-desired, desired];
BMPs = zeros(size(WRescale,1),size(WRescale,2), numel(dar));

% figure;
% blah = 1;
disp('Creating BMPs...')
for iAR = 1:numel(dar)
  % Desired response
  tDesired = desiredResps(:,:,iAR);

  % Invert nonlinearity
  if bLinear
    tINR = (tDesired - RepFitParams(:, :, 2)) ./ RepFitParams(:, :, 1);
  else
    tINR = INR(RepFitParams, tDesired) .* RepMaskBinned;
    tINR = imgaussfilt(tINR, 3);
  end
%     fprintf('INR: Min is %0.8f; Max is %0.8f/n',...
%       min(abs(tINR(:))),max(abs(tINR(:))));   

  % Weight by opsin profile
%       tWeighted = tINR .* W;
  tWeighted = tINR;
%     fprintf('Weighted: Min is %0.8f; Max is %0.8f/n',...
%       min(abs(tWeighted(:))),max(abs(tWeighted(:))));

  % Deconvolve by SIRF
  %      tX = {tWeighted};
  %        tX = deconvlucy(tX,sirf);
  %      tX = tX{2};
%   tDesired = desiredResps(:,:,iAR);

%     if (bScaleBMPs)
%       scale = ampScale(tAmp);
%     elseif (~bScaleBMPs)
%       scale = 1;
%     end

  uniques = unique(tWeighted(:));
  tWeighted(tWeighted == inf) = 0; %uniques(end - 1);

  tDecon = deconvlucy(tWeighted,sirf);
  if (bIndivNorm)
    tDecon = tDecon/max(abs(tDecon(:)));
  end
%     fprintf('Decon: Min is %0.8f; Max is %0.8f/n',...
%       min(abs(tDecon(:))),max(abs(tDecon(:))));

  % Convert to projector resolution
  tPattern = imresize(tDecon,Cam2Proj);
  tPattern = tPattern./(max(tPattern(:))/max(tDecon(:)));
%     fprintf('tPattern: Min is %0.8f; Max is %0.8f/n',...
%       min(abs(tPattern(:))),max(abs(tPattern(:))));

%       tPattern = tPattern - mean(mean(tPattern(~RepMaskBinned), 1), 2);
%       tPattern = tPattern - min(tPattern(:));

  % Find mean response in top right
  meanBaseline = mean(mean(tPattern(1:100, end-100:end), 1), 2);

  tPattern = tPattern - meanBaseline;
  tPattern = tPattern ./ max(tPattern(:));

  BMPs(:,:,iAR) = tPattern;

%       sprintf('pattern mode = %g/n', round(mode(tPattern(:)), 10))
%       subplot(2, 3, blah);histogram(tPattern(:))
%       blah = blah + 1;
end

if (~bIndivNorm)
  % Normalize within BMPs
  BMPs = BMPs/max(abs(BMPs(:)));
end

for iAR = 1:numel(dar)
  if (dar(iAR) <= 0)
    tOrt = 0;
  elseif (dar(iAR) > 0)
    tOrt = 90;
  end
  
  % Create BMP
  tPattern = BMPs(:,:,iAR);
  tBMP = zeros(PHalfY,PSizeX);
  tBMP(:,PIndWX) = tPattern(WIndPY,:);
  tBMPProj = kron(tBMP,ones(2,1)); % double every row
  tBMPProj = floydHalftone(255.*tBMPProj); % floyd halftoning
  tBMPProj = flipud(255.*tBMPProj);

  tBMP = floydHalftone(255.*tBMP);
  tBMP = 255.*tBMP;

%       if (bProjectorFlip)
%         tBMP = rot90(tBMP, 2);
%       end

  filename = sprintf('ProjGssS%03dAR%05dO%05d.bmp',...
                     round(1*1e2), ...
                     round(abs(dar(iAR))*1e4), ...
                     tOrt*1e2);
  imwrite(tBMPProj,fullfile(pathname,filename),'bmp');

  filename = sprintf('DisplayBMPS%03dAR%05dO%05d.bmp',...
                     round(1*1e2), ...
                     round(abs(dar(iAR))*1e4), ...
                     tOrt*1e2);
  imwrite(tBMP,fullfile(pathname,filename),'bmp');

end

disp('Done!')

%% Functions

function out = INR(Param, R, Mask)

x = (R .* (Param(:, :, 2) .^ Param(:, :, 3))) ./ (Param(:, :, 1) - R);
% x = x .* Mask;
x(x < 0) = 0;

out = nthroot(x, Param(:, :, 3));

end


