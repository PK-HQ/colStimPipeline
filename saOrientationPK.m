function hgfMapOrt = saOrientation(fnTS,TS,DAVersion)

%% function hgfMapOrt = saOrientation(fnTS,TS,DAVersion)
%
% Calculate and display orientation map
%
% Input:
%   - fnTS is the TS.mat filename, including pathname and extention (.mat).
%   - TS is the trial structure.
%   - DAVersion is the DA version.
%
% Output:
%   - hgfMapOrt is the handle of figure for orientation map
%
%
% YC at ES lab
% Created on Apr. 20, 2012
% Last modified on May 9, 2022

%% Timer starts
TimerStart = now;
disp('Calculate orientation map!');

%% Check inputs and outputs
[PathName,fnRoot] = fileparts(fnTS);
fnRoot = regexprep(fnRoot,'TS','');
fnData = fullfile(PathName,[fnRoot,'Orientation.mat']);
fnFig = fullfile(PathName,[fnRoot,'Orientation']);

hgfMapOrt = [];

[FileName,PathName] = ...
  uigetfile('*FFT*.mat','Select a FFT file',PathName);
drawnow;pause(0.1);

if all(~FileName)
  return;
end

if ~contains(FileName,'FFT')|| ...
   contains(FileName,'Bin')|| ...
   ~contains(FileName,fnRoot)
  beep;
  disp('Wrong FFT file!');
  return;
end

fnFFT = fullfile(PathName,FileName);

%% Parameters
DefaultValues;

iCondBlank = find(TS.Header.Conditions.TypeCond==0);
iCondStim = find(TS.Header.Conditions.TypeCond>0);
nCondStim = length(iCondStim);

Ort = TS.Header.Conditions.GGOrtCond(iCondStim);
nOrt = nCondStim;

FermiLowCutOff = 0.8;  % cyc/mm
FermiHighCutOff = 3;  % cyc/mm

SizePxl = TS.Header.Imaging.SizePxl;
SigmaRMSPxl = 0.25/SizePxl;  % pixel

LabelStim = cell(1,nCondStim);
for i = 1:nCondStim
  LabelStim{i} = sprintf('%d^o',Ort(i));
end

%% Create a ps/pdf file
hgfFig = figure;
annotation('TextBox',[0.05,0.7,0.9,0.1], ...
           'String','Orientation map', ...
           'HorizontalAlignment','Center', ...
           'FontWeight','bold', ...
           'FontSize',18, ...
           'LineStyle','None');
annotation('TextBox',[0.05,0.5,0.9,0.1], ...
           'String',fnRoot, ...
           'HorizontalAlignment','Center', ...
           'FontWeight','bold', ...
           'FontSize',16, ...
           'LineStyle','None');
annotation('TextBox',[0.05,0.3,0.9,0.1], ...
           'String','', ...
           'HorizontalAlignment','Center', ...
           'FontWeight','bold', ...
           'FontSize',16, ...
           'LineStyle','None');
annotation('TextBox',[0.05,0.2,0.9,0.1], ...
           'String',date, ...
           'HorizontalAlignment','Center', ...
           'FontWeight','bold', ...
           'FontSize',14, ...
           'LineStyle','None');

%% Load files
load(fnFFT,'DataTrial');

if ~isreal(DataTrial)
  DataTrial = abs(DataTrial);
end

[Height,Width,nTrial] = size(DataTrial);

%% Low-pass filtration
DataTrial = FilterFermi2D(DataTrial,0,FermiHighCutOff,SizePxl);

%% Average
RespCond = zeros(Height,Width,nCondStim);
DPCond = zeros(Height,Width,nCondStim);
for i = 1:nCondStim
  RespCond(:,:,i) = ...
    mean(DataTrial(:,:,cell2mat(TS.Header.Index.iValidBLKCond ...
                                (iCondStim(i)))),3);
  DPCond(:,:,i) = ...
    CalculateDPrime( ...
      DataTrial(:,:,cell2mat(TS.Header.Index.iValidBLKCond ...
                             (iCondStim(i)))), ...
      DataTrial(:,:,cell2mat(TS.Header.Index.iValidBLKCond ...
                             (iCondBlank))),3);
end
RespCond = ...
  RespCond- ...
  repmat(mean(DataTrial(:,:,cell2mat(TS.Header.Index.iValidBLKCond ...
                                     (iCondBlank))),3),[1,1,nCondStim]);

DisplayMap(RespCond,LabelStim, ...
           sprintf(['Response map with low-pass filtration ', ...
                    '(%1.1f cyc/mm)'], ...
                   FermiHighCutOff));
DisplayMap(DPCond,LabelStim, ...
           sprintf('d'' map with low-pass filtration (%1.1f cyc/mm)', ...
                   FermiHighCutOff));

%% High-pass filtration
[Spct1D,SF1D] = ...
  CalculateFFTAmp1D(RespCond-repmat(mean(RespCond,3), ...
                                    [1,1,nCondStim]),SizePxl);
Spct1D = mean(Spct1D,3);
Spct1D = Spct1D/max(Spct1D(2:end));
nSF1D = length(SF1D);

RespCondFilt = FilterFermi2D(RespCond,FermiLowCutOff,inf,SizePxl);

%% Remove the mean response map
RespCondFilt = RespCondFilt-repmat(mean(RespCondFilt,3),[1,1,nCondStim]);

hgfFig(end+1) = ...
  DisplayMap(RespCondFilt,LabelStim, ...
             sprintf(['d'' map with band-pass filtration ', ...
                      '(%1.1f~%1.1f cyc/mm)'], ...
                     FermiLowCutOff,FermiHighCutOff));
colormap gray;

%% PCA for orientation data
nPCAComp = 2;
[PCACoef,PCAScore,~,~,PCAExpl] = ...
  pca(reshape(RespCondFilt,[Height*Width,nOrt]));

PCAComp = reshape(PCAScore,[Height,Width,nOrt]);

RespCondPCA = ...
  reshape(PCAScore(:,1:nPCAComp)*PCACoef(:,1:nPCAComp)', ...
          [Height,Width,nOrt]);

CorrResp = zeros(1,nOrt);
for i = 1:nOrt
  tCorr = corrcoef(RespCondPCA(:,:,i),RespCondFilt(:,:,i));
  CorrResp(i) = tCorr(1,2);
end

LabelPCA = cell(1,nCondStim);
for i = 1:nCondStim
  LabelPCA{i} = sprintf('PCA%d---%0.1f%%',i,PCAExpl(i));
end
hgfFig(end+1) = ...
  DisplayMap(PCAComp,LabelPCA,'PCA');
colormap gray;

%% Mask the area with poor signals
RespRMS = sqrt(mean(RespCondFilt.^2,3));
% RespRMS = ...
%   sqrt((RespRMS.^2+imag(hilbert(imag(hilbert(RespRMS))'))'.^2)/2);
tCood.X = -round(3*SigmaRMSPxl):round(3*SigmaRMSPxl);
tCood.Y = tCood.X;
GssRMS = ...
  FuncWoNGaussian2D([1,0,0,0,0,SigmaRMSPxl,SigmaRMSPxl],tCood);
GssRMS = GssRMS/sum(GssRMS(:));
RespRMS = CalculateConv2D(RespRMS,GssRMS);

MethodMask = 2;  % define valid area: 0-no mask, 1-response, 2-d', 3-RMS, 4-all
ThsdResp = max(RespCond(:))/3;  % define valid area based on response
ThsdDP = 6;  % define valid area based on d'
ThsdRMS = max(RespRMS(:))/3;  % define valid area based on RMS

Mask = mean(DPCond,3)>ThsdDP;
Crop = [1,1,Width,Height];  % [X0,Y0,Width,Height]

CLimResp = [0,max(RespCond(:))];
CLimDP = [0,10];
CLimRMS = [0,max(RespRMS(:))];

hgfFig(end+1) = figure;
orient landscape;
annotation('TextBox',[0,0.9,1,0.1], ...
           'String', ...
           sprintf('%s---Select a mask',fnRoot), ...
           'HorizontalAlignment','Center', ...
           'FontWeight','bold', ...
           'FontSize',12, ...
           'LineStyle','None');
axes('Position',[0.1,0.45,0.36,0.48]);
imagesc(mean(RespCond,3)*100,CLimResp*100);
colorbar;
axis image off;
hgoRect1 = ...
  rectangle('Position',Crop, ...
            'EdgeColor','r', ...
            'FaceColor','None', ...
            'Curvature',[0,0], ...
            'LineWidth',2);
title('1. Response (%)', ...
      'FontSize',12, ...
      'FontWeight','Bold');
axes('Position',[0.1,0,0.36,0.48]);
imagesc(mean(DPCond,3),CLimDP);
colorbar;
axis image off;
hgoRect2 = ...
  rectangle('Position',Crop, ...
            'EdgeColor','r', ...
            'FaceColor','None', ...
            'Curvature',[0,0], ...
            'LineWidth',2);
title('2. d''', ...
      'FontSize',12, ...
      'FontWeight','Bold');
axes('Position',[0.6,0.45,0.36,0.48]);
imagesc(RespRMS*100,CLimRMS*100);
colorbar;
axis image off;
hgoRect3 = ...
  rectangle('Position',Crop, ...
            'EdgeColor','r', ...
            'FaceColor','None', ...
            'Curvature',[0,0], ...
            'LineWidth',2);
title('3. RMS (%)', ...
      'FontSize',12, ...
      'FontWeight','Bold');
hgaMask = axes('Position',[0.6,0,0.36,0.48]);

while true
  LabelMask = ...
    {'No mask', ...
     sprintf('Mask (Resp>%0.2f%%)',ThsdResp*100), ...
     sprintf('Mask (d''>%g)',ThsdDP), ...
     sprintf('Mask (RMS>%0.2f%%)',ThsdRMS*100), ...
     sprintf('Mask (Resp>%0.2f%%,d''>%g,RMS>%0.2f%%)', ...
             ThsdResp*100,ThsdDP,ThsdRMS*100)};
  set(hgoRect1,'Position',Crop);
  set(hgoRect2,'Position',Crop);
  set(hgoRect3,'Position',Crop);
  axes(hgaMask);
  imagesc(Mask,[0,1]);
  colorbar;
  axis image off;
  title({LabelMask{MethodMask+1}, ...
         sprintf('ROI = [%d,%d,%d,%d]',Crop)}, ...
        'FontSize',12, ...
        'FontWeight','Bold');

  Answer = ...
    inputdlg({['Mask method:', ...
               '0-no mask, 1-response, 2-d'' (default)), 3-RMS, 4-all:'], ...
              'Threshold for response (%, 1/3 max)', ...
              'Threshold for d''', ...
              'Threshold for RMS (%, 1/3 max)', ...
              'Select ROI (pixel, x0,y0,width,height):'}, ...
             'Select a mask. Cancle when done',1, ...
             {sprintf('%g',MethodMask), ...
              sprintf('%g',ThsdResp*100), ...
              sprintf('%g',ThsdDP), ...
              sprintf('%g',ThsdRMS*100), ...
              sprintf('%d,%d,%d,%d',Crop)});
  drawnow;pause(0.1);

  if isempty(Answer)
    break;
  end

  MethodMask = str2double(Answer{1});
  ThsdResp = str2double(Answer{2})/100;
  ThsdDP = str2double(Answer{3});
  ThsdRMS = str2double(Answer{4})/100;
  Crop = str2num(Answer{5});

  if isempty(MethodMask)||~isscalar(MethodMask)|| ...
     MethodMask<0||MethodMask>4|| ...
     isempty(ThsdResp)||~isscalar(ThsdResp)|| ...
     isempty(ThsdDP)||~isscalar(ThsdDP) || ...
     isempty(ThsdRMS)||~isscalar(ThsdRMS) || ...
     isempty(Crop)||length(Crop(:))~=4|| ...
     Crop(1)<1||Crop(1)>Width|| ...
     Crop(2)<1||Crop(2)>Height|| ...
     Crop(3)<1||Crop(3)>Width|| ...
     Crop(4)<1||Crop(4)>Height|| ...
     (Crop(1)+Crop(3)-1)>Width|| ...
     (Crop(2)+Crop(4)-1)>Height
    beep;
    warndlg(['Some required parameters are wrong!', ...
             ' Please re-input.'],'Warning!','modal');
    drawnow;pause(0.1);
    uiwait;
  else
    MaskCrop = false(Height,Width);
    MaskCrop(Crop(2)-1+(1:Crop(4)),Crop(1)-1+(1:Crop(3))) = true;
    switch MethodMask
      case 1  % response
        Mask = MaskCrop&mean(RespCond,3)>ThsdResp;
      case 2  % d'
        Mask = MaskCrop&mean(DPCond,3)>ThsdDP;
      case 3  % RMS
        Mask = MaskCrop&RespRMS>ThsdRMS;
      case 4  % all
        Mask = ...
          MaskCrop&mean(RespCond,3)>ThsdResp& ...
          mean(DPCond,3)>ThsdDP&RespRMS>ThsdRMS;
      otherwise
        Mask = MaskCrop;
    end
  end
end

%% Calculate cross correlation
CorrOrt = zeros(nCondStim,nCondStim);
for i = 1:nCondStim
  tRespCondFilt1 = RespCondFilt(:,:,i);
  for j = 1:nCondStim
    tRespCondFilt2 = RespCondFilt(:,:,j);
    tCorr = corrcoef(tRespCondFilt1(Mask),tRespCondFilt2(Mask));
    CorrOrt(i,j) = tCorr(1,2);
  end
end

[Ort1,Ort2] = meshgrid(mod(Ort,180),mod(Ort,180));
OrtDiff = 90-abs(mod(Ort1-Ort2,180)-90);
OrtDiff1D = unique(OrtDiff(:))';
CorrOrt1D = zeros(1,length(OrtDiff1D));  % convert to 1D
for k = 1:length(OrtDiff1D)
  CorrOrt1D(k) = mean(CorrOrt(OrtDiff==OrtDiff1D(k)));
end

hgfFig(end+1) = figure;
annotation('TextBox',[0,0.9,1,0.1], ...
           'String',[fnRoot,'---Spectrum and correlation'], ...
           'HorizontalAlignment','Center', ...
           'FontWeight','bold', ...
           'FontSize',12, ...
           'LineStyle','None');
axes('Position',[0.1,0.3,0.36,0.48]);
plot(SF1D,Spct1D,'-k','LineWidth',2);
axis square;
axis([0,3,0,1]);
xlabel('SF (cyc/mm)', ...
       'FontSize',12, ...
       'FontWeight','Bold');
ylabel('Normalized amplitude', ...
       'FontSize',12, ...
       'FontWeight','Bold');
title('1D spectrum', ...
      'FontSize',12, ...
      'FontWeight','Bold');
axes('Position',[0.6,0.3,0.36,0.48]);
plot(OrtDiff1D,CorrOrt1D,'-k','LineWidth',2);
axis square;
axis([0,90,-1,1]);
line(xlim,[1,1]*0,'Color','k','LineStyle','--');
set(gca,'XTick',0:45:90);
xlabel('Orientation difference (deg)', ...
       'FontSize',12, ...
       'FontWeight','Bold');
ylabel('Cross-correlation', ...
       'FontSize',12, ...
       'FontWeight','Bold');
title('Correlation', ...
      'FontSize',12, ...
      'FontWeight','Bold');
box on;

% Calculate orientation map
[MapOrt,AmpOrt,TCOrt] = CalculateMapTCPROrt(RespCondFilt,mod(Ort,180),Mask);

OrtUniq = unique(mod(Ort,180));
nOrtUniq = length(OrtUniq);
TCOrtUniq = zeros(nOrtUniq,nOrtUniq);
for i = 1:nOrtUniq
  ti = (mod(Ort,180)==OrtUniq(i));
  for j = 1:nOrtUniq
    tj = (mod(Ort,180)==OrtUniq(j));
    TCOrtUniq(i,j) = mean(mean(TCOrt(ti,tj),1),2);
  end
end

OrtComb = circshift(OrtUniq,[0,floor(nOrtUniq/2)]);
OrtComb(OrtComb>=90) = OrtComb(OrtComb>=90)-180;
OrtComb = [OrtComb,180+OrtComb(1)];
for i = 1:nOrtUniq
  TCOrtUniq(i,:) = circshift(TCOrtUniq(i,:),[0,floor(nOrtUniq/2)+1-i]);
end

DCMean = mean(RespCond(repmat(Mask,[1,1,nCondStim])));  % DC component

TCOrtUniq = [TCOrtUniq,TCOrtUniq(:,1)]+DCMean;
TCOrtUniq = TCOrtUniq/max(TCOrtUniq(:));
TCOrtComb = mean(TCOrtUniq,1);
TCOrtComb = TCOrtComb/max(TCOrtComb(:));

PI = [0.1,0,0.8,30];
PLB = [0,-1e-5,0,1];
PUB = [1,1e-5,1,90];
PFGssOrt = ...
  lsqcurvefit('FuncWoNGaussianCirc1D',PI, ...
              OrtComb,TCOrtComb, ...
              PLB,PUB,sOptOptm);
TCOrtCombFit = FuncWoNGaussianCirc1D(PFGssOrt,OrtComb);

t = corrcoef(TCOrtCombFit,TCOrtComb);
RSQTCOrtComb = t(1,2)^2;

ColorMap = [0.5,0.5,0.5;hsv(18)];
ColorAmp = repmat(0:0.01:1,[3,1])';

tMapOrt = repmat(MapOrt,[1,1,3]);

% MapAmpOrt is colored version of MapOrt
MapAmpOrt = NaN(Height,Width,3);
for j = 1:size(ColorMap,1)
  % tiMapAmpOrt = logical of tMapOrt, if it falls between 2 orientations (binSize=10, -10:180; -10 = null)
  tiMapAmpOrt = (tMapOrt>=(j-2)*10 & tMapOrt<=(j-1)*10);
  %RGB color for this band (:,:,3) 
  tColor = repmat(reshape(ColorMap(j,:),[1,1,3]),[Height,Width,1]);
  %Assign color to each pixel
  MapAmpOrt(tiMapAmpOrt) = tColor(tiMapAmpOrt);
end

%tAmpOrt = amp map, scaled by max
tAmpOrt = AmpOrt/max(AmpOrt(:))*3;
tAmpOrt(tAmpOrt>1) = 1; % why clip like that? 
tAmpOrt(MapOrt==-10) = 1; % why assign like that?
%MapAmpOrt = MapAmpOrt.*Amp
MapAmpOrt = MapAmpOrt.*repmat(tAmpOrt,[1,1,3]);

hgfFig(end+1) = figure;
hgfMapOrt = hgfFig(end);
annotation('TextBox',[0,0.9,1,0.1], ...
           'String',[fnRoot,'---Orientation map and tuning curve'], ...
           'HorizontalAlignment','Center', ...
           'FontWeight','bold', ...
           'FontSize',12, ...
           'LineStyle','None');
hOrtMap=axes('Position',[0.2,0.5,0.3,0.4]);
imagesc(MapAmpOrt);
axis image;
set(gca,'XAxisLocation','Top','XTick',[],'YTick',[]);
xlabel(sprintf('%0.1fmm',Crop(3)*TS.Header.Imaging.SizePxl), ...
       'FontSize',10, ...
       'FontWeight','Bold');
ylabel(sprintf('%0.1fmm',Crop(4)*TS.Header.Imaging.SizePxl), ...
       'FontSize',10, ...
       'FontWeight','Bold');
annotation('Line',[0.5-(0.3/(Crop(3)*TS.Header.Imaging.SizePxl)),0.5],[0.52,0.52], ...
           'LineWidth',4);
annotation('TextBox',[0.3,0.52,0.2,0.06], ...
           'String','1 mm', ...
           'HorizontalAlignment','Right', ...
           'FontWeight','bold', ...
           'FontSize',12, ...
           'LineStyle','None');
axes('Position',[0.55,0.5,0.01,0.3]);
imagesc(permute(flipdim(ColorMap,1),[1,3,2]));
set(gca,'YAxisLocation','Right', ...
        'XTick',[], ...
        'YTick',[1,9.5,18], ...
        'YTickLabel',{'180','90','0'});
ylabel('Orientation ^o', ...
       'FontSize',10, ...
       'FontWeight','Bold');
box on;
axes('Position',[0.7,0.5,0.01,0.3]);
imagesc(permute(flipdim(ColorAmp,1),[1,3,2]));
set(gca,'YAxisLocation','Right', ...
        'XTick',[], ...
        'YTick',[1,101], ...
        'YTickLabel',{sprintf('%0.3f%%',max(AmpOrt(:))/3*100),'0'});
ylabel('Amplitude', ...
       'FontSize',10, ...
       'FontWeight','Bold');
box on;
axes('Position',[0.2,0.1,0.24,0.32]);
hold on;
plot(OrtComb,TCOrtComb,'bo');
plot(OrtComb,TCOrtCombFit,'-r', ...
     'LineWidth',2);
hold off;
axis square tight;
axis([-90,90,0,1.1]);
set(gca,'XTick',-90:45:90)
xlabel('Orientation difference (deg)', ...
       'FontSize',12, ...
       'FontWeight','Bold');
ylabel('Normalized amp', ...
       'FontSize',12, ...
       'FontWeight','Bold');
text(0,0.2,sprintf('\\sigma=%4.1f^o',PFGssOrt(4)), ...
     'FontSize',10, ...
     'FontWeight','Bold');
box on;
saveas(hgfFig(end),[fnFig,'.fig'],'fig');
saveas(hgfFig(end),[fnFig,'.jpg'],'jpeg');
system(['chgrp eslab ',fnFig,'.fig']);
system(['chgrp eslab ',fnFig,'.jpg']);

%% Convert ps to pdf
SavePDF(hgfFig,[fnFig,'.pdf']);
system(['chgrp eslab ',fnFig,'.pdf']);

%% Save files
save(fnData, ...
     'iCondBlank','iCondStim','nCondStim','DAVersion', ...
     'Ort','nOrt','Height','Width', ...
     'FermiLowCutOff','FermiHighCutOff','SizePxl','LabelStim', ...
     'RespCond','DPCond','RespCondFilt','RespRMS', ...
     'SF1D','Spct1D','nSF1D','Crop', ...
     'MethodMask','ThsdResp','ThsdDP','Mask', ...
     'CorrOrt','OrtDiff','OrtDiff1D','CorrOrt1D', ...
     'MapOrt','AmpOrt','TCOrt','', ...
     'OrtComb','TCOrtComb','DCMean','', ...
     'OrtUniq','nOrtUniq','TCOrtUniq','', ...
     'PFGssOrt','TCOrtCombFit','RSQTCOrtComb','', ...
     'MapAmpOrt','ColorMap','ColorAmp','', ...
     'nPCAComp','PCACoef','PCAScore','PCAExpl', ...
     'PCAComp','RespCondPCA','CorrResp','', ...
     '','','','', ...
     '','','','', ...
     '');
system(['chgrp eslab ',fnData]);

%% Timer ends
TimerEnd = now;
disp(['Session started at ',datestr(TimerStart)]);
disp(['Session ended at ',datestr(TimerEnd)]);


