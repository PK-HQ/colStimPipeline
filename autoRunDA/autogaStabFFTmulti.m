function Message = autogaStabFFTmulti(fnTS,TS,DAVersion)

%% function Message = gaFFT(fnTS,TS,DAVersion)
%
% FFT BLK data over time for each trial
%
% STEPS:
%   - Normalize BLK data
%   - FFT data over time for each trial
%   - Take data at peak frequency
%
% Input:
%   - fnTS is the TS.mat filename, including pathname and extention (.mat).
%   - TS is the trial structure.
%   - DAVersion is the DA version.
%
% Output:
%   - Message is a character string for processing information.
%
% Parameters:
%   - iFramesFFT is the frames for fft.
%   - PeakFreq is the peak frequency in Hz.
%   - lShowFig indicates if show FFT spectrum
%   - iFramesNmlTrial is the frames for normalization in trial
%
% Note:
%   - Add average for each condition (YC, Mar. 21, 2012)
%   - Convert FFT amplitude to peak-to-trough (YC, Mar. 21, 2012)
%   - Change FFTTrial and FFTCond to DataTrial and DataCond
%     (YC, Apr. 23, 2013)
%   - Add FFT amplitude
%     (YC, 02/4/2014)
%   - Make normalization unique by trial
%     (YC, Sep. 17, 2014)
%   - Automatic choose appropriate frame indices depending on camera
%     Add options to stabilize images
%   - Software stabilization option added
%     (SC, Apr. 3, 2017)
%   - Adapt automatica frame selection to work with Retinotopy data
%     Use common script and minimize software stabilization saved output.
%     (SC, Apr. 13, 2017)
%   - Adadpt dF/F normalizing frames for Allegro timing.
%     (SC, Mar. 27, 2019)
%
% YC/SC at ES lab
% Created on Apr. 14, 2008
% Last modified on Apr. 13, 2017

%% Timer starts
TimerStart = now;
disp('FFTing BLK files ... busy');


%% Check inputs and outputs
[PathName,fnRoot] = fileparts(fnTS);
fnRoot = regexprep(fnRoot,'TS','');

% estimate stimulus frequency
if (isfield(TS.Header.Conditions,'DurationOn'))
  StimPeriod = TS.Header.Conditions.DurationOn + ...
               TS.Header.Conditions.DurationOff;
  StimRep    = TS.Header.Conditions.FrameNum;  
elseif (isfield(TS.Header.Conditions,'ClipDuration'))
  StimPeriod = TS.Header.Conditions.ClipDuration;
  StimRep    = TS.Header.Conditions.ClipRepeat;  
else
  StimPeriod = 200;
  StimRep    = 5;  
end

StimFreq = 1000 / StimPeriod(1);

% camera frame offset
if (isfield(TS.Header.Imaging,'CameraName'))
  if (TS.Header.Version.PCLVersion(1) == '4')  %Allegro
    NmlFrameOffset = 3;
    FFTFrameOffset = 1;
  else
    NmlFrameOffset = 1;
    FFTFrameOffset = 1;
  end
else
  NmlFrameOffset = 3;
  FFTFrameOffset = 3;
end  

iFramesNmlTrialEnd = NmlFrameOffset + ceil((TS.Header.Delay.PreStimulus+20)/ TS.Header.Imaging.FrameDuration);
iFramesFFTStart    = iFramesNmlTrialEnd + FFTFrameOffset;
iFramesFFTEnd      = 1000 * StimRep(1) / StimFreq(1);
iFramesFFTEnd      = iFramesFFTEnd / TS.Header.Imaging.FrameDuration ;
iFramesFFTEnd      = round(iFramesFFTEnd + iFramesFFTStart - 1);

iFramesFFT = 14:140;
PeakFreq = 33; 
lShowFig = 0; 
iFramesNmlTrial = 1:13;

%% Additional Processing

% Ask for stabilization
stabcfg = autoaskAboutStabilization(fnTS,TS,iFramesNmlTrial);
dostab  = ~isempty(stabcfg);

if (dostab)
  stabstr = 'Stab';
else
  stabstr = '';
end


% Ask for detrend
dodetrend = 'Yes';%questdlg('Detrend FFT?','Detrend','Yes','No','Yes');
dodetrend = strcmp(dodetrend,'Yes');

if (dodetrend)
  dtrendstr = 'D';
else
  dtrendstr = '';
end


fnFFT = sprintf('%s%s%sFFTS%03dE%03dPF%04d.mat', ...
                fnRoot, ...
                stabstr, dtrendstr, ...
                iFramesFFT(1), iFramesFFT(end), ...
                round(PeakFreq*100))
              
fnFFT = fullfile(PathName, fnFFT);

                          
if ~FileExistsOverwrite(fnFFT)
  Message = 'FFT file has already existed!';
  disp(Message);
  return;
end

%% Load average
load(fullfile(PathName,[fnRoot,'Avg']),'AvgAll');

nFrames = length(iFramesFFT);
AvgAll = repmat(AvgAll,[1,1,nFrames]);

%% Parameter for figure
if lShowFig
  Freq = ...
    ((1:nFrames)-1)/nFrames/TS.Header.Imaging.BLKHeader.FrameDuration*1000;
  Width = size(AvgAll,2);
  Height = size(AvgAll,1);
  iWidth = round(Width/4):round(Width/4*3);
  iHeight = round(Height/4):round(Height/4*3);
end

%% FFT
DataTrial = ...
  complex(zeros(TS.Header.Imaging.BLKHeader.FrameHeight, ...
                TS.Header.Imaging.BLKHeader.FrameWidth, ...
                TS.Header.Index.nBLKTrial));

StabRefFrame = zeros(TS.Header.Imaging.BLKHeader.FrameHeight, ...
                TS.Header.Imaging.BLKHeader.FrameWidth, ...
                TS.Header.Index.nBLKTrial);          

NmlTrial = zeros(size(DataTrial,1),size(DataTrial,2),size(DataTrial,4));

% For fft lienar detrend
FFTtrend = fft(iFramesFFT - mean(iFramesFFT));

tic;
% figure

%% Do single trial FFT
% Get stim freq per condition
trialCond=extractfield(TS.Trial,'CurrCond');
StimFreqCond = 1000 ./ StimPeriod;

%% For each cond, get all trials and do FFT @ correct stimFreqCond for stim and blank conds







tmp=TS.Trial(TS.Header.Index.iBLKTrial(:));
stimConds=unique(double(extractfield(tmp,'CurrCond')));
OItrials=double(extractfield(tmp,'TrialNum'));
nStimConds=numel(unique(stimConds));
for condNo=stimConds
    tic
    fprintf('Loading BLK files for cond %03d/%03d ... \n', ...
          condNo,nStimConds);
      
    % find stimcond == cond
    condTrialIdx=find(stimConds==condNo);

    allTrialsOfCond=OItrials(condTrialIdx);
    
    for i=1:numel(allTrialsOfCond)
        %extract single trial, 512x512x140 frames
        tCondTrials(:,:,:,i) = ...
           blkGetData(fullfile(PathName, ...
                            TS.Trial(allTrialsOfCond(i)). ...
                            fnBLK.name), ...
                   TS.Header.Imaging.BLKHeader);
    end

    tCondAvgs(:,:,:,condNo)=mean(tCondTrials,4);
    tCondTrials=[];
    toc
end

for condNo=1:nStimConds
    tic
  fprintf('FFT-ing BLK file for trial cond %03d/%03d ... \n', ...
          condNo,nStimConds);
      
  % extract single cond, 512x512x140 frames
  tCondAvg = squeeze(tCondAvgs(:,:,:,condNo));
           
  % stabilize
  if (dostab)
    % re-use previously calculated motion vector
    if (~isempty(stabcfg.cfg{condNo}))
      tCondAvg = stabcfg.ipu.imwarp(tCondAvg, stabcfg.cfg{condNo}.tforms, stabcfg.cfg{condNo});
    else
      [tCondAvg, stabcfg.cfg{condNo}] = stabcfg.ipu.procimage(tCondAvg);
      StabRefFrame(:,:,condNo) = mean(tCondAvg,3);
    end
  end
            
  NmlTrial(:,:,condNo) = mean(tCondAvg(:,:,iFramesNmlTrial),3);

  tFFTCond = ...
    fft(tCondAvg(:,:,iFramesFFT)./ ...
        repmat(NmlTrial(:,:,condNo), [1,1,nFrames]),[],3);

  % FFT detrend
  if (dodetrend)
    tFFTCond = reshape(tFFTCond,[],numel(iFramesFFT));
    B         = FFTtrend(2:end)' \ tFFTCond(:,2:end)';
    BTrend    = B' * FFTtrend(2:end);
    tFFTCond(:,2:end) = tFFTCond(:,2:end) - BTrend;
    tFFTCond = reshape(tFFTCond, size(NmlTrial,1), size(NmlTrial,2), numel(iFramesFFT));
  end


  
  
  
  
   % here set the PeakFreq, iPeakFreq, FFT2Amp for this trial
   PeakFreqCond(condNo)=StimFreqCond(condNo);
   iPeakFreqCond(condNo) = ...
     round(PeakFreqCond(condNo)/1000*nFrames*TS.Header.Imaging.BLKHeader.FrameDuration)+1;
   if iPeakFreqCond(condNo)==1
     FTT2AmpCond = 1/nFrames;  % convert to peak to trough, DC component
   else
     FTT2AmpCond = 4/nFrames;  % convert to peak to trough, Non-DC component
   end

  DataCond(:,:,condNo) = tFFTCond(:,:,iPeakFreqCond(condNo))*FTT2AmpCond;
  
  fprintf('%f seconds elapsed\n', toc);
end



% % Stabilize across trials: NOT FULLY WORKING YET
% [StabNmlTrial, stabcfg.refcfg] = stabcfg.ipu.procimage(NmlTrial);
% StabDataTrial = stabcfg.ipu.imwarp(DataTrial .* NmlTrial, stabcfg.refcfg.tforms, stabcfg.refcfg);
% 
% DataTrial = StabDataTrial ./ StabNmlTrial;
% NmlTrial  = StabNmlTrial;


%% Average
%{
DataCond = ...
  complex(zeros(TS.Header.Imaging.BLKHeader.FrameHeight, ...
                TS.Header.Imaging.BLKHeader.FrameWidth, ...
                TS.Header.NumValidCond));
for i = 1:TS.Header.NumValidCond
  DataCond(:,:,i) = ...
    mean(DataTrial(:,:,TS.Header.Index.iValidBLKCond{i}),3);
end
%}

%% Subtract mean

if TS.Header.NumBlankCond
  MeanBlank = mean(DataTrial(:,:, ...
          cell2mat(TS.Header.Index.iValidBLKCond( ...
            1:TS.Header.NumBlankCond))),3);
  DataTrial = ...
    DataTrial- ...
        repmat(MeanBlank, ...
               [1,1,TS.Header.Index.nBLKTrial]);
  DataCond = zeros(size(DataCond));
  for i = 1:TS.Header.NumValidCond
    DataCond(:,:,i) = ...
      mean(DataTrial(:,:,TS.Header.Index.iValidBLKCond{i}),3);
  end
end


%% Save file
disp('Saving "FFT BLK" file ...');

if (dostab && ~isfield(stabcfg,'src'))  
  gaSaveStab(stabcfg, StabRefFrame, stabcfg.fnStabCfg);  
end

if numel(DataTrial) > 2^26
  save(fnFFT, '-v7.3', ...
       'DataTrial','DataCond', 'NmlTrial', 'iFramesFFT', 'MeanBlank', ...
       'PeakFreqCond','iPeakFreqCond','FTT2AmpCond', ...
       'iFramesNmlTrial','DAVersion');
else
  save(fnFFT, ...
       'DataTrial','DataCond','NmlTrial', 'iFramesFFT', 'MeanBlank', ...
       'PeakFreqCond','iPeakFreqCond','FTT2AmpCond', ...
       'iFramesNmlTrial','DAVersion');
end

   
system(['chgrp eslab ',fnFFT]);



%% Calculate FFT amplitude
if TS.Header.NumBlankCond
  DataTrial = abs(DataTrial);
  DataCond = zeros(size(DataCond));
  for i = 1:TS.Header.NumValidCond
    DataCond(:,:,i) = ...
      mean(DataTrial(:,:,TS.Header.Index.iValidBLKCond{i}),3);
  end
  save(regexprep(fnFFT,'FFT','FFTAmp'), ...
       'DataTrial','DataCond','NmlTrial','iFramesFFT', ...
       'PeakFreqCond','iPeakFreqCond','FTT2AmpCond', ...
       'iFramesNmlTrial','DAVersion');
  system(['chgrp eslab ',regexprep(fnFFT,'FFT','FFTAmp')]);
else
  warndlg('No blank condition, so no FFT (Amp) amplitude file!');
end


%% Timer ends
TimerEnd = now;
disp('FFTing BLK files ... done!');
disp(['Session started at ',datestr(TimerStart)]);
disp(['Session ended at ',datestr(TimerEnd)]);
Message = 'Successful!';
disp(Message);


