function [TS, ErrorCode] = autotsMatchBLK(TS, fnTS)

ErrorCode = 0;

[PathName,fnRoot] = fileparts(fnTS);
fnRoot = regexprep(fnRoot,'TS','');

%% Synchonize and verify the BLK files
TS.Header.Index.iBLKTrial = find([TS.Trial.FlagOIBLK]==true);
TS.Header.Index.nBLKTrial = length(TS.Header.Index.iBLKTrial);

fnBLKAll = dir(fullfile(PathName,'*.BLK'));
TS.Header.Imaging.nBLKFile = length(fnBLKAll);

BLKCond = zeros(1,TS.Header.Imaging.nBLKFile);
BLKSerial = zeros(1,TS.Header.Imaging.nBLKFile);
for i = 1:TS.Header.Imaging.nBLKFile
  ti = strfind(fnBLKAll(i).name,'_');
  ti2 = strfind(fnBLKAll(i).name(1:ti(1)),'C');
  % Find condition C##_ in the BLK filename
  BLKCond(i) = sscanf(fnBLKAll(i).name((ti2(end)+1):(ti(1)-1)),'%d');
  % Find serial _E##B###.BLK in the BLK filename
  BLKSerial(i) = ...
    sscanf(fnBLKAll(i).name((ti(end)+1):end),'E%*dB%d.BLK')'+1;
end

[t,ti] = sort(BLKSerial);
BLKCond = BLKCond(ti);
if isfield(TS.Trial,'OIStimID')
  TrialCond = [TS.Trial(TS.Header.Index.iBLKTrial).OIStimID];
else
  warning('OIStimID field not found. Using CurrCond instead.');
  TrialCond = [TS.Trial(TS.Header.Index.iBLKTrial).CurrCond];
  TrialCond = TrialCond - 1;
end
  

if TS.Header.Imaging.nBLKFile~=TS.Header.Index.nBLKTrial || ...
    any(TrialCond~=BLKCond)
  Message = ...
    sprintf(['BLK files (%d) do NOT match the log file (%d)! ', ...
             'Please check the message in workspace!'], ...
            TS.Header.Imaging.nBLKFile,TS.Header.Index.nBLKTrial);
  disp(Message);
  disp('');
  disp('');
  disp('List of unmatched BLK files!');
  disp('Note: trial condition number is OIStimID in log file');
  disp('');
  disp('Trial#  TrialCond#  BLKCond#  Unmatch      BLKFilename      ');
  for i = 1:min(TS.Header.Imaging.nBLKFile,TS.Header.Index.nBLKTrial)
    fprintf('%4d      %4d      %4d      %4d     %s\n', ...
            TS.Header.Index.iBLKTrial(i), ...
            TrialCond(i),BLKCond(i),TrialCond(i)~=BLKCond(i), ...
            fnBLKAll(ti(i)).name);
  end
  ErrorCode = 6;  % inconsistent value(s) in BLK files
  return;
end

%% Store BLK filenames in TS
% if min(BLKSerial)>1
%   BLKSerial = BLKSerial-min(BLKSerial)+1;
% end
%
% for i = 1:TS.Header.Imaging.nBLKFile
%   TS.Trial(TS.Header.Index.iBLKTrial(BLKSerial(i))).fnBLK = fnBLKAll(i);
% end

for i = 1:numel(TS.Header.Index.iBLKTrial)
  TS.Trial(TS.Header.Index.iBLKTrial(i)).fnBLK = fnBLKAll(ti(i));
end

%% Find valid BLK trials
TS.Header.Index.iValidBLKTrial = ...
  intersect(TS.Header.Index.iValidTrial,TS.Header.Index.iBLKTrial);
TS.Header.Index.nValidBLKTrial = length(TS.Header.Index.iValidBLKTrial);

TS.Header.Index.iValidBLKTrialCond = cell(1,TS.Header.NumValidCond);
TS.Header.Index.nValidBLKTrialCond = zeros(1,TS.Header.NumValidCond);
for i = 1:TS.Header.NumValidCond
  TS.Header.Index.iValidBLKTrialCond{i} = ...
    TS.Header.Index.iValidBLKTrial ...
    ([TS.Trial(TS.Header.Index.iValidBLKTrial).CurrCond]==i);
  TS.Header.Index.nValidBLKTrialCond(i) = ...
    length(TS.Header.Index.iValidBLKTrialCond{i});
end

[t,TS.Header.Index.iValidBLK,ti] = ...
  intersect(TS.Header.Index.iBLKTrial,TS.Header.Index.iValidBLKTrial);
TS.Header.Index.nValidBLK = length(TS.Header.Index.iValidBLK);
TS.Header.Index.iValidBLK = ...
  reshape(TS.Header.Index.iValidBLK,[1,TS.Header.Index.nValidBLK]);

TS.Header.Index.iValidBLKCond = cell(1,TS.Header.NumValidCond);
TS.Header.Index.nValidBLKCond = zeros(1,TS.Header.NumValidCond);
for i = 1:TS.Header.NumValidCond
  [t,TS.Header.Index.iValidBLKCond{i},ti] = ...
    intersect(TS.Header.Index.iBLKTrial, ...
              TS.Header.Index.iValidBLKTrialCond{i});  % intersect vector
  TS.Header.Index.nValidBLKCond(i) = ...
    length(TS.Header.Index.iValidBLKCond{i});
  TS.Header.Index.iValidBLKCond{i} = ...
    reshape(TS.Header.Index.iValidBLKCond{i}, ...
            [1,TS.Header.Index.nValidBLKCond(i)]);  % intersect
end

%% Read BLK header
TS.Header.Imaging.BLKHeader = ...
  blkGetHeader(fullfile(PathName,fnBLKAll(end).name));

%% Read custom user info (Spencer)
nfo = textscan(TS.Header.Imaging.BLKHeader.User,'%s','Delimiter','=;');
nfo = nfo{1};
if numel(nfo)>1
  for ii = 1:2:numel(nfo),
    num = str2double(nfo{ii+1});
    if (~isempty(num) && ~isnan(num))
      nfo{ii+1} = num;
    end
    TS.Header.Imaging.(nfo{ii}) = nfo{ii+1};
  end
end

%% Add some parameters (i.e size) in TS.Header.Imaging
if ~isfield(TS.Header.Imaging,'Resolution')
  if TS.Header.Imaging.BLKHeader.FrameWidth>504 || ...
     TS.Header.Imaging.BLKHeader.FrameHeight>504
   if TS.Header.Imaging.BLKHeader.FrameWidth<=1024 && ...
     TS.Header.Imaging.BLKHeader.FrameHeight<=1024
      TS.Header.Imaging.Resolution = 1024;
   else
      TS.Header.Imaging.Resolution = 2048;
   end
  else
    TS.Header.Imaging.Resolution = 512;
  end
end

TS.Header.Imaging.FrameRate = ...
  round(1000/TS.Header.Imaging.BLKHeader.FrameDuration);
TS.Header.Imaging.FrameWidth = ...
  TS.Header.Imaging.BLKHeader.FrameWidth;
TS.Header.Imaging.FrameHeight = ...
  TS.Header.Imaging.BLKHeader.FrameHeight;
TS.Header.Imaging.nFrame = ...
  TS.Header.Imaging.BLKHeader.NFramesPerStim;
TS.Header.Imaging.FrameDuration = ...
  TS.Header.Imaging.BLKHeader.FrameDuration;
TS.Header.Imaging.Duration = ...
  TS.Header.Imaging.BLKHeader.FrameDuration* ...
  TS.Header.Imaging.BLKHeader.NFramesPerStim;

if ~isfield(TS.Header.Imaging,'SizeCameraSensor')
  switch TS.Header.Imaging.BLKHeader.ActiveSystemID
    case 12
      TS.Header.Imaging.SizeCameraSensor = 0.014*1024;  % mm
    case 103
      TS.Header.Imaging.SizeCameraSensor = 0.012*1024;  % mm
    case 255 % sCMOS
      TS.Header.Imaging.SizeCameraSensor = 0.0065*2048; % mm
    case 254 % MV1
      TS.Header.Imaging.SizeCameraSensor = 0.0106*1024; % mm
    otherwise
      error('Wrong camera system!');
  end
end

% Bugfix 2020.7.25 [Spencer]: v422 saved OILens values as strings
% Moved here by Yuzhi, 2022/08/02
if ~isnumeric(TS.Header.ExptInfo.OILensTop)
  TS.Header.ExptInfo.OILensTop = sscanf(TS.Header.ExptInfo.OILensTop,'%dmm');
  TS.Header.ExptInfo.OILensBottom = sscanf(TS.Header.ExptInfo.OILensBottom,'%dmm');
end

% Correction factor for the focallength, Yuzhi, 2022/08/02
TS.Header.Imaging.SizePxlCorrectionFactor = 1.05;

TS.Header.Imaging.SizePxl = ...
  TS.Header.Imaging.SizeCameraSensor* ...
  TS.Header.ExptInfo.OILensBottom/TS.Header.ExptInfo.OILensTop/ ...
  TS.Header.Imaging.Resolution*TS.Header.Imaging.SizePxlCorrectionFactor;
TS.Header.Imaging.Width = ...
  TS.Header.Imaging.SizePxl*TS.Header.Imaging.FrameWidth;
TS.Header.Imaging.Height = ...
  TS.Header.Imaging.SizePxl*TS.Header.Imaging.FrameHeight;

switch TS.Header.ExptInfo.OIType
  case TS.Header.DEF.IMAGING.NONE
    TS.Header.Imaging.VSDI{1} = 'No imaging';
  case TS.Header.DEF.IMAGING.DYE
    TS.Header.Imaging.VSDI{1} = 'Voltage sensitive dye imaging';
  case TS.Header.DEF.IMAGING.INTRINSIC
    TS.Header.Imaging.VSDI{1} = 'Intrinsic imaging';
  case TS.Header.DEF.IMAGING.GCAMP
    TS.Header.Imaging.VSDI{1} = 'GCaMP imaging';
  otherwise
    TS.Header.Imaging.VSDI{1} = 'Unknown';
end
TS.Header.Imaging.VSDI{2} = ...
  sprintf('Resolution: %d x %d pixel at %dHz', ...
          TS.Header.Imaging.Resolution, ...
          TS.Header.Imaging.Resolution, ...
          TS.Header.Imaging.FrameRate);
TS.Header.Imaging.VSDI{3} = ...
  sprintf('Frame: %d x %d pixel with %d frames', ...
          TS.Header.Imaging.FrameWidth, ...
          TS.Header.Imaging.FrameHeight, ...
          TS.Header.Imaging.nFrame);
TS.Header.Imaging.VSDI{4} = ...
  sprintf('Region: %0.2fmm x %0.2fmm in %dms', ...
          TS.Header.Imaging.Width,TS.Header.Imaging.Height, ...
          round(TS.Header.Imaging.Duration));
        
TS.Header.Imaging.VSDI{5} = ...
  sprintf('Lenses: top %dmm, bottom %dmm', ...
          TS.Header.ExptInfo.OILensTop,TS.Header.ExptInfo.OILensBottom);

TS.Header.Imaging.VSDI{6} = ...
  sprintf('Pixel size: %0.3fmm', ...
          TS.Header.Imaging.SizePxl);
        
%% Save condition and outcome tables in a TXT file
fnTables = fullfile(PathName,[fnRoot,'Tables.txt']);
fid = fopen(fnTables,'wt');

fprintf(fid,'Protocol Name: %s\r\n',TS.Header.ProtocolName);
fprintf(fid,'TS(DA) version: %s\r\n',TS.Header.Version.TSVersion);
fprintf(fid,'PCL version: %s\r\n',TS.Header.Version.PCLVersion);

fprintf(fid,'\r\n');
fprintf(fid,'\r\n');

for k = 1:length(TS.Header.ConditionTable)
  fprintf(fid,'%s\r\n',TS.Header.ConditionTable{k});
end

fprintf(fid,'\r\n');
fprintf(fid,'\r\n');

for k = 1:length(TS.Header.OutcomeTable)
  fprintf(fid,'%s\r\n',TS.Header.OutcomeTable{k});
end

fprintf(fid,'\r\n');
fprintf(fid,'\r\n');

fprintf(fid,'Comments: %s\r\n',TS.Header.Comments);

fclose(fid);
system(['chgrp eslab ',fnTables]);

