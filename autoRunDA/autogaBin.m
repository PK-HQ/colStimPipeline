function [Message,TS] = autogaBin(fnTS,TS,DAVersion,mode,dostab,binParams)

%% function [Message,TS] = gaBin(fnTS,TS,DAVersion)
%
% Bin BLK data over space for each trial
%
% STEPS:
%   - Bin BLK data over space for each trial
%   - Find and remove outliers
%   - Normalize binned data
%   - Average data for each condition
%
% Input:
%   - fnTS is the TS.mat filename, including pathname and extention (.mat).
%   - TS is the trial structure.
%   - DAVersion is the DA version.
%
% Output:
%   - Message is a character string for processing information.
%   - TS is the trial structure after removing outliers
%
% Parameters:
%   - BinSize must be a positive integer.
%   - CrtrRmOl is the criterion of outlier in unit of standard deviation
%   - iFramesNmlTrial is the frames for normalization in trial
%
% Note:
%   - Update TS.Header.Index, including iValidTrial and iValidTrialCond
%     (YC, Nov. 5, 2010)
%   - Change BinTrial and BinCond to DataTrial and DataCond
%     (YC, Apr. 23, 2013)
%   - Add time course for individual trials in ~Avg.mat
%     (YC, Sep.10, 2013)
%   - Check illuminantion saturation in blkBin.m
%     (YC, Jan. 25, 2014)
%   - Find outliers by removing mean over all conditions rather than 
%     over each condition in "Bin BLK"
%     (YC, Feb. 9, 2014)
%   - Make normalization unique by trial
%     (YC, Sep. 17, 2014)
%   - Update maximum pixel value (saturation value) to work with PCO data
%     (SC, Dec. 14, 2016)
%   - Downsizing with pre-filter option added (set to default)
%   - Software stabilization option added
%     (SC, Dec. 14, 2016)
%
%
% YC at ES lab
% Created on Apr. 17, 2008
% Last modified on May 14, 2019

%% Timer starts
TimerStart = now;
disp('Downsizing BLK files ... busy');

%% Check inputs and outputs
[PathName,fnRoot] = fileparts(fnTS);
fnRoot = regexprep(fnRoot,'TS','');

%% Bin or Downsample
%{
mode = ...
  questdlg('Do you want to BIN or DOWNSAMPLE with Gaussian pre-filter?', ...
	'Image Stabilization', ...
	'Bin','Down','Down');
%}

% camera frame offset
if (isfield(TS.Header.Imaging,'CameraName'))
  FrameOffset = 1;
else
  FrameOffset = 3;
end  

if (isfield(TS.Header.Conditions,'DurationPreTOn'))
  DelayPreStim       = TS.Header.Conditions.DurationPreTOn / TS.Header.Imaging.FrameDuration;
else
  DelayPreStim       = 0;
end
iFramesNmlTrialStart = floor(DelayPreStim + 1);
iFramesNmlTrialEnd   = DelayPreStim + FrameOffset + ...
                       (TS.Header.Delay.PreStimulus+20)/ TS.Header.Imaging.FrameDuration;
iFramesNmlTrialEnd   = ceil(iFramesNmlTrialEnd);                     

switch (mode)

  case 'Bin'
    %{
    Answer = ...
      inputdlg({'Bin size', ...
                'Criterion for outlier (*SD)', ...
                'Remove saturated trials', ...
                'Index of frames for normalization by trial'}, ...
               'Parameters',1, ...
               {'8','5','0', ...
                sprintf('%d:%d',iFramesNmlTrialStart,iFramesNmlTrialEnd)});
    drawnow;pause(0.1);
    %}
    if isempty(binParams)
      TS = [];
      Message = 'Canceled!';
      disp(Message);
      return;
    end

    BinSize         = str2double(binParams{1});
    CrtrRmOl        = str2double(binParams{2});
    lRmSaturate     = str2double(binParams{3});
    iFramesNmlTrial = str2num(binParams{4});
    TargetRes       = 1; % dummy
    
  case 'Down'
    %{
    Answer = ...
      inputdlg({'Target resolution', ...
                'Criterion for outlier (*SD)', ...
                'Remove saturated trials', ...
                'Index of frames for normalization by trial'}, ...
               'Parameters',1, ...
               {'128','5','0', ...
                sprintf('%d:%d',iFramesNmlTrialStart,iFramesNmlTrialEnd)});
    drawnow;pause(0.1);
    %}

    if isempty(binParams)
      TS = [];
      Message = 'Canceled!';
      disp(Message);
      return;
    end

    TargetRes       = str2double(binParams{1});
    CrtrRmOl        = str2double(binParams{2});
    lRmSaturate     = str2double(binParams{3});
    iFramesNmlTrial = str2num(binParams{4});
    BinSize         = round(TS.Header.Imaging.Resolution/TargetRes); % dummy
    
end

if isempty(BinSize)||~isscalar(BinSize)||(BinSize<=0)|| ...
   isempty(CrtrRmOl)||~isscalar(CrtrRmOl)||(CrtrRmOl<=0)|| ...
   isempty(lRmSaturate)||~isscalar(lRmSaturate)||(lRmSaturate<0)|| ...
   isempty(iFramesNmlTrial)||any(iFramesNmlTrial<=0)|| ...
   any(iFramesNmlTrial>TS.Header.Imaging.BLKHeader.NFramesPerStim) || ...
   isempty(TargetRes)||~isscalar(TargetRes)||(TargetRes<=0)||(TargetRes>TS.Header.Imaging.Resolution)
  beep;
  TS = [];
  Message = 'Wrong parameters, please check!';
  disp(Message);
  return;
end

%% Init stabilization 
% Ask for stabilization
%{
dostab = questdlg('Would you like a run image stabilization?', ...
	'Image Stabilization', ...
	'Yes','No','Yes');
%}
dostab = strcmp(dostab,'Yes');

stabcfg   = [];
if dostab
  
  % load existing file if exists
  fnStabCfg = fullfile(PathName,sprintf('%sStabCfg.mat',fnRoot));
  if (exist(fnStabCfg,'file'))
    loadstabcfg = questdlg('Image stabilization configuration file found. Re-use?', ...
    'Image Stabilization', ...
    'Yes','No','Yes');  
    if (strcmp(loadstabcfg,'Yes'))
      stabcfg = load(fnStabCfg);
      stabcfg.src = fnStabCfg;
      stabcfg.ipu = WideField.ipu.ImageStabilizer(stabcfg.ipu);
    end
  end
  
  % user config
  if (isempty(stabcfg))
    
    if (isfield(TS.Header.Imaging,'CameraName'))
      divlinesX = '[]';
    else
      divlinesX = '[128 129 256 257 384 385]';
    end
    
    
    binParams = {'1.5',num2str(iFramesNmlTrial(ceil(numel(iFramesNmlTrial)/2))), ...
                '[138 138 375 375]',divlinesX,'[]'};

    if isempty(binParams)
      Message = 'Canceled!';
      disp(Message);
      return;
    end

    spatialHP = str2double(binParams{1});
    iFrameRef = str2double(binParams{2});
    regwin      = str2num(binParams{3});
    divlinesX   = str2num(binParams{4});
    iFramesStab = str2num(binParams{5});

    if isempty(spatialHP)||~isscalar(spatialHP)||(spatialHP<0)|| ...
       isempty(iFrameRef)||~isscalar(iFrameRef)||(iFrameRef<0)||(iFrameRef>TS.Header.Imaging.nFrame)|| ...
       isempty(regwin)||regwin(1)<1||regwin(2)<1||regwin(3)>TS.Header.Imaging.FrameWidth||regwin(4)>TS.Header.Imaging.FrameHeight|| ...
       any(divlinesX<1) || any(divlinesX>TS.Header.Imaging.FrameWidth)|| ...
       any(iFramesStab<1) || any(iFramesStab>TS.Header.Imaging.nFrame)
      beep;
      return;
    end

    stabcfg.cfg = cell(TS.Header.Index.nBLKTrial,1);
    stabcfg.ipu = WideField.ipu.ImageStabilizer( ...
      'spatialHP',     spatialHP, ...
      'refframe',      iFrameRef, ...
      'regwinX',       regwin([1 3]), ...
      'regwinY',       regwin([2 4]), ...
      'divlinesX',     divlinesX, ...
      'regressFrames', iFramesStab, ...
      'SizePxl',       TS.Header.Imaging.SizePxl, ...
      'FrameRate',     TS.Header.Imaging.FrameRate ...
      );
    stabcfg.ipu
  end
  
  fnBin = fullfile(PathName,sprintf('%sStab%s%03d.mat',fnRoot,mode,BinSize));
                            
else

  fnBin = fullfile(PathName,sprintf('%s%s%03d.mat',fnRoot,mode,BinSize));
  
end

if ~FileExistsOverwrite(fnBin)
  TS = [];
  Message = 'Bin file has already existed!';
  disp(Message);
  return;
end

%% Show wait bar
hWaitBar = ...
  waitbar(0,{'Downsizing BLK files. Please wait ...', ...
             'To cancel, close this window'}, ...
          'Name','0% done');
drawnow;pause(0.1);
Count = 0;
tic;

%% Parameters
BinHeight = floor(TS.Header.Imaging.BLKHeader.FrameHeight/BinSize);
BinWidth = floor(TS.Header.Imaging.BLKHeader.FrameWidth/BinSize);
nFrames = TS.Header.Imaging.BLKHeader.NFramesPerStim;
nBLKTrial = TS.Header.Index.nBLKTrial;

BinSizeRmOl = floor([BinHeight,BinWidth]/2);
BinSizeDisp = ceil([BinHeight,BinWidth]/10);

%% Bin
DataTrial = zeros(BinHeight,BinWidth,nFrames,nBLKTrial);
AvgTrial = ...
  zeros(TS.Header.Imaging.BLKHeader.FrameHeight, ...
        TS.Header.Imaging.BLKHeader.FrameWidth,nBLKTrial);
MaxBLK = zeros(1,nBLKTrial);

StabRefFrame = zeros(TS.Header.Imaging.BLKHeader.FrameHeight, ...
                TS.Header.Imaging.BLKHeader.FrameWidth, ...
                TS.Header.Index.nBLKTrial);          
         
tic;
for j = 1:nBLKTrial
  fprintf('Downsizing BLK file for trial %03d/%03d ... ', ...
          j,nBLKTrial);
  tTrial = ...
    blkGetData(fullfile(PathName, ...
                        TS.Trial(TS.Header.Index.iBLKTrial(j)). ...
                        fnBLK.name), ...
               TS.Header.Imaging.BLKHeader);
  MaxBLK(j) = max(tTrial(:));
  
  % stabilize
  if (dostab)
    % re-use previously calculated motion vector
    if (~isempty(stabcfg.cfg{j}))
      tTrial = stabcfg.ipu.imwarp(tTrial, stabcfg.cfg{j}.tforms, stabcfg.cfg{j});
    else
      if (iFrameRef > 0)
        StabRefFrame(:,:,j) = tTrial(:,:,iFrameRef);
      end
      [tTrial, stabcfg.cfg{j}] = stabcfg.ipu.procimage(tTrial);
    end
  end  
  
  switch (mode)
    case 'Bin'
      DataTrial(:,:,:,j) = BinND(tTrial,[BinSize,BinSize]);
    case 'Down'
      DataTrial(:,:,:,j) = WideField.tools.Downsize2D(tTrial,BinSize);
  end
  
  AvgTrial(:,:,j) = mean(tTrial,3);
  
  fprintf('%f seconds elapsed\n', toc);
  
  % Show wait bar
  Count = Count+1;
  Ratio = Count/(nBLKTrial+1);
  drawnow;pause(0.1);
  if ishandle(hWaitBar)
    waitbar(Ratio,hWaitBar);
    set(hWaitBar,'Name',sprintf('%d%% done, %d sec remaining', ...
                                floor(Ratio*100), ...
                                round(toc*(1/Ratio-1))));
    drawnow;pause(0.1);
  else
    if ~strcmp(questdlg('Are you sure to cancel?','Confirm to cancel'), ...
               'Yes')
      hWaitBar = ...
        waitbar(Ratio,{'Binning BLK files. Please wait ...  ', ...
                       'To cancel, close this window'}, ...
                'Name',sprintf('%d%% done, %d sec remaining', ...
                               floor(Ratio*100), ...
                               round(toc*(1/Ratio-1))));
      drawnow;pause(0.1);
    else
      beep;
      TS = [];
      Message = 'Canceled!';
      disp(Message);
      return;
    end
  end
end


% Spencer's update: Dec 14, 2016
% work out the saturation level from the data values
maxI = max(AvgTrial(:));
if (maxI <= 255)
  saturationLevel = 255;
elseif (maxI <= 4095)
  saturationLevel = 4095;
elseif (maxI <= 65535)
  saturationLevel = 65535;
else
  % can't really tell if software binning is employed: guess
  saturationLevel = 2^ceil(log2(maxI));
end

iBLKSaturation = find(MaxBLK>=saturationLevel);
iBLKSaturation = intersect(iBLKSaturation,TS.Header.Index.iValidBLK);
nBLKSaturation = length(iBLKSaturation);
TS.Header.Index.iBLKSaturation = iBLKSaturation;
TS.Header.Index.nBLKSaturation = nBLKSaturation;

%% Show wait bar
drawnow;pause(0.1);
if ishandle(hWaitBar)
  waitbar(Ratio,hWaitBar,{'Removing outliers. Please wait ...', ...
                          'To cancel, close this window'});
  set(hWaitBar,'Name','Removing outliers ...');
  drawnow;pause(0.1);
else
  if ~strcmp(questdlg('Are you sure to cancel?','Confirm to cancel'), ...
             'Yes')
    hWaitBar = ...
      waitbar(Ratio,{'Removing outliers. Please wait ...', ...
                     'To cancel, close this window'}, ...
              'Name','Removing outliers ...');
    drawnow;pause(0.1);
  else
    beep;
    Message = 'Canceled!';
    disp(Message);
    return;
  end
end

%% Find outliers
DataTrialRmOl = zeros(2,2,nFrames,nBLKTrial);
for j = 1:nBLKTrial
  DataTrialRmOl(:,:,:,j) = BinND(DataTrial(:,:,:,j),BinSizeRmOl);
end

DataTrialRmOl = ...  % normalize by each trial
  DataTrialRmOl./repmat(mean(DataTrialRmOl,3),[1,1,nFrames,1]);

DataTrialRmOl = ...  % remove mean over all conditions
  DataTrialRmOl- ...
  repmat(mean(DataTrialRmOl(:,:,:,TS.Header.Index.iValidBLK),4), ...
         [1,1,1,TS.Header.Index.nBLKTrial]);

tSD = std(DataTrialRmOl(:,:,:,TS.Header.Index.iValidBLK),[],4);

iBLKOutlier = ...
  find(squeeze(any(any(any(abs(DataTrialRmOl)>repmat(CrtrRmOl*tSD, ...
    [1,1,1,nBLKTrial]),1),2),3)))';
iBLKOutlier = intersect(iBLKOutlier,TS.Header.Index.iValidBLK);
nBLKOutlier = length(iBLKOutlier);

tnBLKOutlier = ...
  sum(any(abs(DataTrialRmOl(:,:,:,TS.Header.Index.iValidBLK))> ...
          repmat(CrtrRmOl*tSD,[1,1,1,TS.Header.Index.nValidBLK]),3),4);
tiValidBLK = TS.Header.Index.iValidBLK;

TS.Header.Index.iBLKOutlier = iBLKOutlier;
TS.Header.Index.nBLKOutlier = nBLKOutlier;

%% Remove outliers and illuminantion saturations
if lRmSaturate
  iBLKRemove = union(iBLKOutlier,iBLKSaturation);
else
  iBLKRemove = iBLKOutlier;
end
nBLKRemove = length(iBLKRemove);
TS.Header.Index.iBLKRemove = iBLKRemove;
TS.Header.Index.nBLKRemove = nBLKRemove;

TS.Header.Index.iValidTrial = ...
  setdiff(TS.Header.Index.iValidTrial, ...
          TS.Header.Index.iBLKTrial(iBLKRemove));
TS.Header.Index.nValidTrial = length(TS.Header.Index.iValidTrial);
TS.Header.Index.iValidBLKTrial = ...
  setdiff(TS.Header.Index.iValidBLKTrial, ...
          TS.Header.Index.iBLKTrial(iBLKRemove));
TS.Header.Index.nValidBLKTrial = length(TS.Header.Index.iValidBLKTrial);
TS.Header.Index.iValidBLK = ...
  setdiff(TS.Header.Index.iValidBLK,iBLKRemove);
TS.Header.Index.nValidBLK = length(TS.Header.Index.iValidBLK);

for i = 1:TS.Header.NumValidCond
  TS.Header.Index.iValidTrialCond{i} = ...
    setdiff(TS.Header.Index.iValidTrialCond{i}, ...
            TS.Header.Index.iBLKTrial(iBLKRemove));
  TS.Header.Index.nValidTrialCond(i) = ...
    length(TS.Header.Index.iValidTrialCond{i});
  TS.Header.Index.iValidBLKTrialCond{i} = ...
    setdiff(TS.Header.Index.iValidBLKTrialCond{i}, ...
            TS.Header.Index.iBLKTrial(iBLKRemove));
  TS.Header.Index.nValidBLKTrialCond(i) = ...
    length(TS.Header.Index.iValidBLKTrialCond{i});
  TS.Header.Index.iValidBLKCond{i} = ...
    setdiff(TS.Header.Index.iValidBLKCond{i},iBLKRemove);
  TS.Header.Index.nValidBLKCond(i) = ...
    length(TS.Header.Index.iValidBLKCond{i});
end

TS.Header.Misc.CrtrRmOl = CrtrRmOl;

%% Figure for outliers and saturations
hgfFig = figure;
annotation('TextBox',[0.05,0.9,0.9,0.1], ...
           'String', ...
           sprintf('%s---Find outliers (%d)', ...
                   fnRoot,nBLKOutlier), ...
           'HorizontalAlignment','Center', ...
           'FontWeight','bold', ...
           'FontSize',12, ...
           'LineStyle','None');
for i = 1:2
  for j = 1:2
    subplot(2,2,(i-1)*2+j);
    hold on;
    plot(squeeze(DataTrialRmOl(i,j,:,tiValidBLK)));
    plot(squeeze(CrtrRmOl*tSD(i,j,:)),'-r','LineWidth',2);
    plot(squeeze(-CrtrRmOl*tSD(i,j,:)),'-r','LineWidth',2);
    hold off;
    axis([1,nFrames,[-1,1]*(CrtrRmOl+1)*max(tSD(:))]);
    xlabel('Frame #', ...
           'FontWeight','bold', ...
           'FontSize',12);
    ylabel('\DeltaF/F', ...
           'FontWeight','bold', ...
           'FontSize',12);
    title(sprintf('%d outlier(s)',tnBLKOutlier(i,j)), ...
          'FontWeight','bold', ...
          'FontSize',12);
  end
end
fnFig = fullfile(PathName,[fnRoot,'Outlier']);
saveas(hgfFig,[fnFig,'.fig'],'fig');
saveas(hgfFig,[fnFig,'.jpg'],'jpeg');
system(['chgrp eslab ',fnFig,'.fig']);
system(['chgrp eslab ',fnFig,'.jpg']);

hgfFig = figure;
annotation('TextBox',[0.05,0.9,0.9,0.1], ...
           'String', ...
           sprintf('%s---Find saturations (%d)', ...
                   fnRoot,nBLKSaturation), ...
           'HorizontalAlignment','Center', ...
           'FontWeight','bold', ...
           'FontSize',12, ...
           'LineStyle','None');
subplot(2,1,1);
plot(MaxBLK,'*');
axis tight;
xlabel('Trial #', ...
       'FontWeight','bold', ...
       'FontSize',12);
ylabel('Maximum Grey level', ...
       'FontWeight','bold', ...
       'FontSize',12);

subplot(2,1,2);

% Spencer's update: Dec 14, 2016
% hist(MaxBLK,0:64:4096);
histx = linspace(0,saturationLevel,65);
hist(MaxBLK,histx);

axis tight;
xlabel('Maximum grey level', ...
       'FontWeight','bold', ...
       'FontSize',12);
ylabel('Count of trials', ...
       'FontWeight','bold', ...
       'FontSize',12);
fnFig = fullfile(PathName,[fnRoot,'Saturation']);
saveas(hgfFig,[fnFig,'.fig'],'fig');
saveas(hgfFig,[fnFig,'.jpg'],'jpeg');
system(['chgrp eslab ',fnFig,'.fig']);
system(['chgrp eslab ',fnFig,'.jpg']);

%% Show wait bar
drawnow;pause(0.1);
if ishandle(hWaitBar)
  waitbar(Ratio,hWaitBar, ...
          {'Averaging and normalizing. Please wait ...', ...
           'To cancel, close this window'});
  set(hWaitBar,'Name','Averaging and normalizing ...');
  drawnow;pause(0.1);
else
  if ~strcmp(questdlg('Are you sure to cancel?','Confirm to cancel'), ...
             'Yes')
    hWaitBar = ...
      waitbar(Ratio,{'Averaging and normalizing. Please wait ...', ...
                     'To cancel, close this window'}, ...
              'Name','Averaging and normalizing ...');
    drawnow;pause(0.1);
  else
    beep;
    Message = 'Canceled!';
    disp(Message);
    return;
  end
end

%% Average and normalize
AvgAll = nanmean(AvgTrial(:,:,TS.Header.Index.iValidBLK),3);
AvgAllBin = BinND(AvgAll,[BinSize,BinSize]);

% Spencer: saves normalization image
NmlTrial = zeros(size(DataTrial,1),size(DataTrial,2),size(DataTrial,4));

for j = 1:nBLKTrial
  tNml = mean(DataTrial(:,:,iFramesNmlTrial,j),3);
  NmlTrial(:,:,j) = squeeze(tNml);
  tNml = repmat(tNml,[1,1,nFrames]);
  DataTrial(:,:,:,j) = DataTrial(:,:,:,j)./tNml;
end

DataCond = zeros(BinHeight,BinWidth,nFrames,TS.Header.NumValidCond);
for i = 1:TS.Header.NumValidCond
  DataCond(:,:,:,i) = ...
    mean(DataTrial(:,:,:,TS.Header.Index.iValidBLKCond{i}),4);
end
AvgCondDisp = BinND(DataCond,BinSizeDisp);

DataTrialCenter = ...
  squeeze(mean(mean(DataTrial(floor(BinHeight/4)+(1:(BinHeight/2)), ...
                              floor(BinWidth/4)+(1:(BinWidth/2)), ...
                              :,:),1),2));
% DataTrialCenter = ...
%   DataTrialCenter./repmat(mean(DataTrialCenter,1),[nFrames,1])-1;

%% Show wait bar
drawnow;pause(0.1);
if ishandle(hWaitBar)
  waitbar(Ratio,hWaitBar,'Saving files. Please wait ...');
  set(hWaitBar,'Name','Saving files ...');
  drawnow;pause(0.1);
else
  if ~strcmp(questdlg('Are you sure to cancel?','Confirm to cancel'), ...
             'Yes')
    hWaitBar = ...
      waitbar(Ratio,'Saving files. Please wait ...', ...
              'Name','Saving files ...');
    drawnow;pause(0.1);
  else
    beep;
    Message = 'Canceled!';
    disp(Message);
    return;
  end
end

%% Save file
disp('Saving "Downsized BLK" file ...');

if (dostab && ~isfield(stabcfg,'src'))
  
  % extract X and Y motion vectors for convenience
  for bb = 1:numel(stabcfg.cfg)
    [stabcfg.cfg{bb}.dx, stabcfg.cfg{bb}.dy] = ...
        WideField.ipu.ImageStabilizer.getMotionVec(stabcfg.cfg{bb}.tforms);
  end
  stabstruct = struct(stabcfg.ipu);
    
  % calculate movement across all reference frames
  stabcfg.ipu.refframe = round(size(StabRefFrame,3)/2);
  [tforms, Xs] = stabcfg.ipu.imregtform(StabRefFrame);
  stabcfg.refcfg.tforms     = tforms;
  stabcfg.refcfg.lumProfile = Xs.lumProfile;
  stabcfg.refcfg.lumModel   = Xs.lumModel;
  [stabcfg.refcfg.dx, stabcfg.refcfg.dy] = WideField.ipu.ImageStabilizer.getMotionVec(tforms);
  
  stabcfg.ipu = stabstruct;
  save(fnStabCfg, '-struct', 'stabcfg');
  system(['chgrp eslab ',fnStabCfg]);
end

if numel(DataTrial) > 2^27
  save(fnBin, '-v7.3', ...
       'DataTrial','DataCond','BinSize','TargetRes','AvgAllBin', ...
       'iFramesNmlTrial','DAVersion','NmlTrial');
else
  save(fnBin, ...
       'DataTrial','DataCond','BinSize','TargetRes','AvgAllBin', ...
       'iFramesNmlTrial','DAVersion','NmlTrial');
end
save(fnTS,'TS','DAVersion');
save(fullfile(PathName,[fnRoot,'Avg.mat']), ...
     'AvgAll','AvgCondDisp', ...
     'CrtrRmOl','DataTrialRmOl','DataTrialCenter','DAVersion');
system(['chgrp eslab ',fnBin]);
system(['chgrp eslab ',fullfile(PathName,[fnRoot,'Avg.mat'])]);

%% Close wait bar
close(hWaitBar);

%% Timer ends
TimerEnd = now;
disp('Downsizing BLK files ... done!');
disp(['Session started at ',datestr(TimerStart)]);
disp(['Session ended at ',datestr(TimerEnd)]);
Message = 'Successful!';
disp(Message);


