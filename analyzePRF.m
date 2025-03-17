%% Setup
%close all
clear

%% Inputs
monkeyStr = 'Chip';
monkeyIDstr = '28';
sessionStr = '20230815';
pdfSuffix='test';

LumRuns = [5];%
SiteID = [1];%

%% 
%conditions
OptoLuminance = [5;10;...
               20;...
               30;...
               50;...
               70;...
               100];

OptoMaxPower = 62.6; %mW, RM4, with ND10
OptolightGuideDi=9; %mm
OptoProjectorArea = pi*(OptolightGuideDi/2)^2; % mm^2
OptoLightPD=(OptoLuminance/100)*OptoMaxPower/OptoProjectorArea;% mW/mm^2
SiteIDStr = {'#1'};

%% Cosmetic
%color
colormapRB=fireice;
midColormapRB=length(colormapRB)/2;
colormapRB=colormapRB(midColormapRB-15:midColormapRB+16,:); %this is correct as midpoint is a whole number, fireice is 64 x 3 not 65 x 3 (without black midpoint)
midColormapRB=length(colormapRB)/2;

colormapR=[0 0 0; colormapRB(midColormapRB+1:end,:)]; %add zero = black
colormapB=[colormapRB(1:midColormapRB,:); 0 0 0]; %add zero = black
colormapRB=[colormapB; 0 0 0; colormapR];

%plot
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultAxesFontSize',12)

%% Load data

DataFFT = cell(numel(SiteID),1);
DataTS = DataFFT;

for iSite = 1:numel(SiteID)
%   fprintf('Site %g!\n',Sites(iSite))
  DataPath = ['Y:/' monkeyStr '/' monkeyStr sessionStr '/'];
  tRunNum = LumRuns(iSite);
  tRunDir = [DataPath 'run' num2str(tRunNum) '/'];
    
  Run = load([tRunDir 'M' monkeyIDstr 'D' sessionStr 'R' num2str(tRunNum) 'StabBin008.mat']);
  DataFFT{iSite} = ...
    load([tRunDir 'M' monkeyIDstr 'D' sessionStr 'R' num2str(tRunNum) 'StabDFFTAmpS004E023PF0400.mat']);
    %'StabBin008FFTAmpS001E025PF0400.mat']);
  DataTS{iSite}  = ...
    load([tRunDir 'M' monkeyIDstr 'D' sessionStr 'R' num2str(tRunNum) 'TS.mat']);
end

RespFFT = ...
  zeros(size(DataFFT{1,1}.DataCond,1),...
        size(DataFFT{1,1}.DataCond,2),...
        numel(SiteID),numel(OptoLuminance));
      
for iSite = 1:numel(SiteID)
  for iLED = 1:numel(OptoLuminance)
    tData = DataFFT{iSite}.DataCond;
    %blank subtracted response
    RespFFT(:,:,iSite,iLED) = tData(:,:,iLED+1)-tData(:,:,1);
  end
end

% RespFFT = flip(RespFFT,4);

%% Define desired datapat and params (selects dir, monkey, session)
pcID=getenv('COMPUTERNAME');
switch pcID
  case {'SEIDEMANN2PANAL'} %SEA
    imagingDataDir='Y:';
    bitmapDataDir='X:';
    monkeyStr='Chip'; %'Roger';
    visDatapath=[imagingDataDir '\' monkeyStr '\' monkeyStr sessionStr '\expressionMonitoring'];
    bitmapDatapath=[bitmapDataDir '\users\PK\bitmap\dutycycleHalftoned\'];
    plotpath='C:\Users\pktan\Box\_Eyal\Columnar read-write\code\plots\PRF\';
  case {'LA-CPSA10646WD'} %ARC
    imagingDataDir='Y:';
    bitmapDataDir='Z:';
    monkeyStr='Chip';
    visDatapath=[imagingDataDir '\' monkeyStr '\' monkeyStr sessionStr '\expressionMonitoring'];
    bitmapDatapath=[bitmapDataDir '\users\PK\bitmap\dutycycleHalftoned\'];
    plotpath='C:\Users\esexpt\Box\_Eyal\Columnar read-write\code\plots\PRF\';
  case {'LA-CPSD077020WD'}
    imagingDataDir='Y:';
    bitmapDataDir='X:';
    monkeyStr='Chip'; %'Roger';
    visDatapath=[imagingDataDir '\' monkeyStr '\' monkeyStr sessionStr '\expressionMonitoring'];
    bitmapDatapath=[bitmapDataDir '\users\PK\bitmap\dutycycleHalftoned\'];
    plotpath='C:\Users\esexpt\Box\_Eyal\Columnar read-write\code\plots\PRF\';
  case {'PSYC-A77304'}
    imagingDataDir='Y:';
    bitmapDataDir='Y:';
    monkeyStr='Chip'; %'Roger';
    visDatapath=[imagingDataDir '\' monkeyStr '\' monkeyStr sessionStr '\expressionMonitoring'];
    bitmapDatapath=[bitmapDataDir '\users\PK\bitmap\dutycycleHalftoned\'];
    plotpath='Y:\Chip\Meta\neurometric\PRF\';
end

pdfFilename=[plotpath monkeyStr sessionStr '-PRF' pdfSuffix];

%% Responses
colormap(colormapR);
figure
[ha,~] = tight_subplot(numel(SiteID)+1,numel(OptoLuminance)+1,.01,[.05,-.05],[.01,.07]);
numRow = numel(SiteID)+1;
numCol = numel(OptoLuminance)+1;
Clim = [min(RespFFT(:)) max(RespFFT(:))];

axes(ha(1))
axis off

for iSite = 1:numel(SiteID)
  axes(ha(iSite*numCol+1))
  txt = sprintf('#%g',SiteID(iSite));
  text(0.92,0.5,txt,'Units','Normalized',...
    'FontSize',14,'FontWeight','normal','HorizontalAlignment','center')
  axis off
end

for iLED = 1:numel(OptoLuminance)
  axes(ha(iLED+1))
  txt = sprintf('%.2f mW/mm^{2}\n(%g%%)',OptoLightPD(iLED),OptoLuminance(iLED));
  %txt = sprintf('%g%% LED\n%0.2f mW/mm^2',LEDPercents(iLED),PowerDensities(iLED));

  text(0.5,0.01,txt,'Units','Normalized',...
    'FontSize',14,'FontWeight','normal','HorizontalAlignment','center')
  axis off
end

for iSite = 1:numel(SiteID)
  for iLED = 1:numel(OptoLuminance)
    axes(ha(iSite*numCol+iLED+1))
    imagesc(RespFFT(:,:,iSite,iLED),Clim)
    axis image
    set(gca,'xticklabel',[],'yticklabel',[])
  end
  % Colorbar
  origSiz = get(gca,'Position');
  colorbar
  set(gca,'Position',origSiz)
end


% LABELS
[axX,hX]=suplabel(['Half-toned luminance'],'t');
[axY,hY]=suplabel(['Site ID'],'y');
set(hX,'FontWeight','bold');
set(hY,'FontWeight','bold');

titlePos = get(hX,'position');
titlePos(2) = 1.05;
set( hX , 'position' , titlePos);

set(gcf,'color','w');
colormap(colormapR);

export_fig(pdfFilename,'-pdf','-nocrop');

%% Responses with ROI overlay
colormap(colormapR);
Dim = size(RespFFT,1);
ROISize = 115;%20;
figure;

CenterResp = zeros(numel(SiteID),numel(OptoLuminance));

for iSite = 1:numel(SiteID)
  for iLED = 1:numel(OptoLuminance)
    if unique(SiteID)==99
        ROI.X=ceil(5*Dim/10)-floor(ROISize/2)+1:ceil(5*Dim/10)+ceil(ROISize/2);
        ROI.Y=ceil((10-5.2)*Dim/10)-floor(ROISize/2)+1:ceil((10-5.2)*Dim/10)+ceil(ROISize/2);
        %ROI.Y=ROI.X; %center
    else
        %smoothing
        filtSigma = 30; %4 or 8
        filtWidth = filtSigma*6; 
        imageFilter=fspecial('gaussian',filtWidth,filtSigma);

        %find max pixel
        x=RespFFT(:,:,iSite,iLED);
        dataFilteredX = nanconv(x,imageFilter, 'nanout');
        [xMax,yMax]=find(dataFilteredX==max(max(dataFilteredX)));
        ROI.X=xMax-floor(ROISize/2)+1:xMax+ceil(ROISize/2);
        ROI.Y=yMax-floor(ROISize/2)+1:yMax+ceil(ROISize/2);
    end
    CenterResp(iSite,iLED) = mean(mean(RespFFT(ROI.X,ROI.Y,iSite,iLED),1),2);
  end
  
  %plot
  subplot(1,numel(SiteID),iSite)
  imagesc(RespFFT(:,:,iSite,end),Clim);colorbar
  hold on
  rectangle('Position',[min(ROI.Y) min(ROI.X) ROISize ROISize],'EdgeColor','k')
  axis image
  set(gca,'xticklabel',[],'yticklabel',[])
  
  %TEXT
  txt = sprintf('#%g',SiteID(iSite));
  text(.5,1.03,txt,'Units','Normalized',...
    'FontSize',14,'FontWeight','normal','HorizontalAlignment','center')
end

% LABELS
[axX,hX]=suplabel([{'ROI across sites','(all responses for 100% half-toned luminance)'}],'t');
%[axY,hY]=suplabel(['Blue light LED (%)'],'y');
set(hX,'FontWeight','bold');
%set(hY,'FontWeight','bold');

set(gcf,'color','w');
colormap(colormapR);

export_fig(pdfFilename,'-append','-pdf','-nocrop');

%% Fitting
clc
% NakaRushton = fittype( @(r,k,n,x) r.*((x.^n)./((x.^n)+(k.^n))));
% NakaRushton = fittype('r.*((x.^n)./((x.^n)+(k^n)))');
Param0 = [0,... % RMax
          0,... % C50 (Asymptotic maximum response amplitude)
          1,... % n
          ];
NakaRushton = ...
  @(Param,x) Param(1).*((x.^Param(3))./((x.^Param(3))+(Param(2).^Param(3))));

NakaRushN1 = ...
  @(Param,x) Param(1).*((x)./((x)+(Param(2))));

% Fit = lsqcurvefit(NakaRushton,Param0,LEDPercents,CenterResp(1,:)');
% 
% figure
% scatter(LEDPercents,CenterResp(1,:))
% hold on
% plot(LEDPercents,NakaRushton(Fit,LEDPercents))

FitParam = zeros(numel(SiteID),numel(Param0));
NormCenterResp = CenterResp;
for iSite = 1:numel(SiteID)
  tResp = CenterResp(iSite,:);
  tResp = tResp-min(tResp);
  tResp = tResp./max(tResp);
  NormCenterResp(iSite,:) = tResp;
end

annot = {};
Col=colormap(colormapR);
Col=Col(round(linspace(1,size(Col,1),4)),:);
figure
hold on
for iSite = 1:numel(SiteID)
  
  tLED = OptoLightPD;%(OptoMaxPower.*(OptoLuminance./100))./OptoProjectorArea;
    
  %scatter(tLED,CenterResp(iSite,:),75,Col(iSite,:),'filled');

  %scatter(tLED,CenterResp(iSite,:),75,Col(iSite,:),'filled');
  
  
%   if Sites(iSite)~=49
    FitParam(iSite,:) = ...
      lsqcurvefit(NakaRushton,Param0,tLED,CenterResp(iSite,:)');
%   else
%     tData = CenterResp(iSite,:)';
%     tData(1:4) = [];
%     tPD = PowerDensities;
%     tPD(1:4) = [];
%     FitParam(iSite,1:2) = ...
%       lsqcurvefit(NakaRushN1,[0 0],tPD,tData);
%     FitParam(iSite,3) = 1;
%   end
  txt = sprintf('R_m_a_x = %0.2f, C_5_0 = %0.2f, n = %0.2f',...
    FitParam(iSite,1),FitParam(iSite,2),FitParam(iSite,3));
  annot{end+1} = txt;
  x = linspace(0,max(tLED));

  plot(x,NakaRushton(FitParam(iSite,:),x)*100,...
    'Color',Col(iSite,:),'LineStyle','-','LineWidth',4)
end
axis square
xlabel(sprintf('Stimulation power density (mW/mm^2)'))
xlim([0 ceil(max(OptoLightPD))])
%xticks(0:0.2:max((OptoMaxPower.*(OptoLuminance./100))./OptoProjectorArea))

% ylim([0 1])
ylabel(sprintf('Average ROI response amplitude (\\DeltaF/F%%)'))
title('Power-response function (Naka-Rushton fit)')
hleg=legend(SiteIDStr,'Location','northwest')
htitle = get(hleg,'Title');
set(htitle,'String','Site ID')


for iSite = 1:numel(SiteID)
  
  scatter(tLED,CenterResp(iSite,:)*100,75,Col(iSite,:),'filled','HandleVisibility','off');
end


set(gcf,'color','w');




export_fig(pdfFilename,'-append','-pdf','-nocrop');
% annotation('textbox',[.35 .7 .2 .2],'string',annot,'fitboxtotext',...
%   'on','fontsize',15,'linewidth',0.5)

