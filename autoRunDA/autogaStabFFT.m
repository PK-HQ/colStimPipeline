function Message = autogaStabFFT(fnTS,TS,DAVersion)

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

Answer = ...
  inputdlg({'Index of frames for FFT', ...
            'Peak frequency (Hz)', ...
            'Show power spectrum for each trial', ...
            'Index of frames for normalization by trial'}, ...
           'Parameters',1, ...
           {sprintf('%d:%d',iFramesFFTStart,iFramesFFTEnd), ...
            num2str(StimFreq),'0',...
            sprintf('1:%d',iFramesNmlTrialEnd)});
% drawnow;pause(0.1);

if isempty(Answer)
  Message = 'Canceled!';
  disp(Message);
  return;
end

iFramesFFT = str2num(Answer{1});
PeakFreq = str2double(Answer{2});
lShowFig = str2double(Answer{3});
iFramesNmlTrial = str2num(Answer{4});

if isempty(iFramesFFT)||any(iFramesFFT<=0)|| ...
   any(iFramesFFT>TS.Header.Imaging.BLKHeader.NFramesPerStim)|| ...
   isempty(PeakFreq)||~isscalar(PeakFreq)||(PeakFreq<0)|| ...
   isempty(lShowFig)||~isscalar(lShowFig)||(lShowFig<0)|| ...
   isempty(iFramesNmlTrial)||any(iFramesNmlTrial<=0)|| ...
   any(iFramesNmlTrial>TS.Header.Imaging.BLKHeader.NFramesPerStim)
  beep;
  Message = 'Wrong parameters, please check!';
  disp(Message);
  return;
end

%% Additional Processing

% Ask for stabilization
stabcfg = askAboutStabilization(fnTS,TS,iFramesNmlTrial);
dostab  = ~isempty(stabcfg);

if (dostab)
  stabstr = 'Stab';
else
  stabstr = '';
end


% Ask for detrend
dodetrend = questdlg('Detrend FFT?','Detrend','Yes','No','Yes');
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

%% Show wait bar
hWaitBar = ...
  waitbar(0,{'FFTing BLK files. Please wait ...  ', ...
             'To cancel, close this window'}, ...
          'Name','0% done');
% drawnow;pause(0.1);
Count = 0;
tic;

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
iPeakFreq = ...
  round(PeakFreq/1000*nFrames*TS.Header.Imaging.BLKHeader.FrameDuration)+1;
if iPeakFreq==1
  FTT2Amp = 1/nFrames;  % convert to peak to trough, DC component
else
  FTT2Amp = 4/nFrames;  % convert to peak to trough, Non-DC component
end
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
for j = 1:TS.Header.Index.nBLKTrial
  fprintf('FFTing BLK file for trial %03d/%03d ... ', ...
          j,TS.Header.Index.nBLKTrial);
  tTrial = ...
    blkGetData(fullfile(PathName, ...
                        TS.Trial(TS.Header.Index.iBLKTrial(j)). ...
                        fnBLK.name), ...
               TS.Header.Imaging.BLKHeader);

  % stabilize
  if (dostab)
    % re-use previously calculated motion vector
    if (~isempty(stabcfg.cfg{j}))
      tTrial = stabcfg.ipu.imwarp(tTrial, stabcfg.cfg{j}.tforms, stabcfg.cfg{j});
    else
      [tTrial, stabcfg.cfg{j}] = stabcfg.ipu.procimage(tTrial);
      StabRefFrame(:,:,j) = mean(tTrial,3);
    end
  end
            
  NmlTrial(:,:,j) = mean(tTrial(:,:,iFramesNmlTrial),3);

  tFFTTrial = ...
    fft(tTrial(:,:,iFramesFFT)./ ...
        repmat(NmlTrial(:,:,j), [1,1,nFrames]),[],3);

  % FFT detrend
  if (dodetrend)
    tFFTTrial = reshape(tFFTTrial,[],numel(iFramesFFT));
    B         = FFTtrend(2:end)' \ tFFTTrial(:,2:end)';
    BTrend    = B' * FFTtrend(2:end);
    tFFTTrial(:,2:end) = tFFTTrial(:,2:end) - BTrend;
    tFFTTrial = reshape(tFFTTrial, size(NmlTrial,1), size(NmlTrial,2), numel(iFramesFFT));
  end
  
  DataTrial(:,:,j) = tFFTTrial(:,:,iPeakFreq)*FTT2Amp;
  % Show figure
  if lShowFig
    if j==1
      hgfSpectrum = figure;
      set(hgfSpectrum,'NumberTitle','off', ...
                      'Name','FFT spectrum');
    else
      set(0,'CurrentFigure',hgfSpectrum);
    end
    tFFTTrialSpectrum = ...
      squeeze(abs(mean(mean(tFFTTrial(iHeight,iWidth,:),1),2)))*FTT2Amp;
    tFFTTrialSpectrum(1) = NaN;
    plot(Freq,tFFTTrialSpectrum, ...
         'Color','b', ...
         'LineWidth',2);
    axis tight;
    xlim([0,PeakFreq*4]);
    line([1,1]*Freq(iPeakFreq),ylim, ...
         'Color','r', ...
         'LineStyle','--', ...
         'LineWidth',2);
    xlabel('Freq (Hz)', ...
           'FontWeight','bold', ...
           'FontSize',12);
    ylabel('FFT amp', ...
           'FontWeight','bold', ...
           'FontSize',12);
    title('Spectum on center', ...
          'FontWeight','bold', ...
          'FontSize',12);
  end
  
  fprintf('%f seconds elapsed\n', toc);

  % Show wait bar
%   drawnow;pause(0.1);
  if ishandle(hWaitBar)
    Count = Count+1;
    Ratio = Count/(TS.Header.Index.nBLKTrial+1);
    waitbar(Ratio,hWaitBar);
    set(hWaitBar,'Name',sprintf('%d%% done, %d sec remaining', ...
                                floor(Ratio*100), ...
                                round(toc*(1/Ratio-1))));
%     drawnow;pause(0.1);
  else
    if ~strcmp(questdlg('Are you sure to cancel?','Confirm to cancel'), ...
               'Yes')
      hWaitBar = ...
        waitbar(Ratio,{'FFTing BLK files. Please wait ...  ', ...
                       'To cancel, close this window'}, ...
                'Name',sprintf('%d%% done, %d sec remaining', ...
                               floor(Ratio*100), ...
                               round(toc*(1/Ratio-1))));
%       drawnow;pause(0.1);
    else
      beep;
      Message = 'Canceled!';
      disp(Message);
      return;
    end
  end
end


% % Stabilize across trials: NOT FULLY WORKING YET
% [StabNmlTrial, stabcfg.refcfg] = stabcfg.ipu.procimage(NmlTrial);
% StabDataTrial = stabcfg.ipu.imwarp(DataTrial .* NmlTrial, stabcfg.refcfg.tforms, stabcfg.refcfg);
% 
% DataTrial = StabDataTrial ./ StabNmlTrial;
% NmlTrial  = StabNmlTrial;


%% close figure
if lShowFig
  close(hgfSpectrum);
end

%% Show wait bar
% drawnow;pause(0.1);
if ishandle(hWaitBar)
  waitbar(Ratio,hWaitBar, ...
          {'Averaging and normalizing. Please wait ...', ...
           'To cancel, close this window'});
  set(hWaitBar,'Name','Averaging for each condition ...');
%   drawnow;pause(0.1);
else
  if ~strcmp(questdlg('Are you sure to cancel?','Confirm to cancel'), ...
            'Yes')
    hWaitBar = ...
      waitbar(Ratio,{'Averaging and normalizing. Please wait ...', ...
                     'To cancel, close this window'}, ...
              'Name','Averaging and normalizing ...');
%     drawnow;pause(0.1);
  else
    beep;
    Message = 'Canceled!';
    disp(Message);
    return;
  end
end

%% Average
DataCond = ...
  complex(zeros(TS.Header.Imaging.BLKHeader.FrameHeight, ...
                TS.Header.Imaging.BLKHeader.FrameWidth, ...
                TS.Header.NumValidCond));
for i = 1:TS.Header.NumValidCond
  DataCond(:,:,i) = ...
    mean(DataTrial(:,:,TS.Header.Index.iValidBLKCond{i}),3);
end

%% Show wait bar
% drawnow;pause(0.1);
if ishandle(hWaitBar)
  waitbar(Ratio,hWaitBar,'Saving file. Please wait ...');
  set(hWaitBar,'Name','Saving file ...');
%   drawnow;pause(0.1);
else
  if ~strcmp(questdlg('Are you sure to cancel?','Confirm to cancel'), ...
             'Yes')
    hWaitBar = ...
      waitbar(Ratio,'Saving file. Please wait ...', ...
              'Name','Saving file ...');
%     drawnow;pause(0.1);
  else
    beep;
    Message = 'Canceled!';
    disp(Message);
    return;
  end
end

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
       'PeakFreq','iPeakFreq','FTT2Amp', ...
       'iFramesNmlTrial','DAVersion');
else
  save(fnFFT, ...
       'DataTrial','DataCond','NmlTrial', 'iFramesFFT', 'MeanBlank', ...
       'PeakFreq','iPeakFreq','FTT2Amp', ...
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
       'PeakFreq','iPeakFreq','FTT2Amp', ...
       'iFramesNmlTrial','DAVersion');
  system(['chgrp eslab ',regexprep(fnFFT,'FFT','FFTAmp')]);
else
  warndlg('No blank condition, so no FFT (Amp) amplitude file!');
end

%% Close wait bar
close(hWaitBar);

%% Timer ends
TimerEnd = now;
disp('FFTing BLK files ... done!');
disp(['Session started at ',datestr(TimerStart)]);
disp(['Session ended at ',datestr(TimerEnd)]);
Message = 'Successful!';
disp(Message);


