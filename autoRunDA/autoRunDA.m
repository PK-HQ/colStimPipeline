%% Define wanted sessions
% def sessions
sessions={'20230502'};%,'20230426'};

% def mainPath
if ispc
  mainPath='Y:/';
elseif contains(getenv('HOSTNAME'),'psy.utexas.edu')
  mainPath='/eslab/data/';
end

%% Apply match blockk for TS files
% --- Gather TS files to be processed (i.e. with no binning done) --- %
%{
TSfilesMatch={};
files=[];
% for session, get TS files
for sessionNo=1:length(sessions)
    files = dir([mainPath 'Chip/Chip' sessions{sessionNo} '/run*/M*TS.mat']);
    %for files, check if already stabbinned. If not, stabbin
    for fileNo=1:length(files)
        % get stab-bin file
        stabBinFile=dir([files(fileNo).folder '/*StabBin008.mat']);
        % check if exist. If doesn't, add to TSfiles for stabbin
        if isempty(stabBinFile)
            TSfilesMatch{end+1,1} = fullfile(files(fileNo).folder,files(fileNo).name);
        end
    end
end

% --- Apply gaStabBin for each TSfile --- %
if ~isempty(TSfilesMatch)
    for TSfileNo=1%:length(TSfilesMatch)
        % load TS
        fnTS=TSfilesMatch{TSfileNo};
        load(fnTS)
        % match blk
        [Message,TS] = autotsMatchBLK(TS,fnTS);
    end
end

%% Apply stabilization and binning for TS files
% --- Gather TS files to be processed (i.e. with no binning done) --- %
TSfilesBin={};
files=[];
% for session, get TS files
for sessionNo=1:length(sessions)
    files = dir([mainPath 'Chip/Chip' sessions{sessionNo} '/run*/M*TS.mat']);
    %for files, check if already stabbinned. If not, stabbin
    for fileNo=1:length(files)
        % get stab-bin file
        stabBinFile=dir([files(fileNo).folder '/*StabBin008.mat']);
        % check if exist. If doesn't, add to TSfiles for stabbin
        if isempty(stabBinFile)
            TSfilesBin{end+1,1} = fullfile(files(fileNo).folder,files(fileNo).name);
        end
    end
end

% --- Apply gaStabBin for each TSfile --- %
if ~isempty(TSfilesBin)
    for TSfileNo=1%:length(TSfilesBin)
        % load TS
        fnTS=TSfilesBin{TSfileNo};
        load(fnTS);

        % load the vdaqlog file, get framerate and define frames to normalize
        [fnTSpath,~] = fileparts(fnTS);
        load([fnTSpath '/Data_vdaqlog.mat'],'VDaqSettings');
        framerateHz=VDaqSettings.datalog.framerate;
        if round(framerateHz)==100 %100Hz
            nmlFrames=num2str(13); %first 13 frames, or 130ms
        elseif round(framerateHz)==20 %20Hz
            nmlFrames=num2str(3); %first 3 frames, or 150ms
        end

        % params for binning
        mode='Bin';
        dostab='Yes';
        binParams={'8', '5', '0', nmlFrames}; %BinSize, CrtrRmOl, rmSaturatedTrials, nmlFrames

        % stabbin
        [Message,TS] = autogaBin(fnTS,TS,DAVersion,mode,dostab,binParams);
    end
end
%}
%% Apply stabilization and FFT for TS files
% --- Gather TS files to be processed (i.e. with no binning done) --- %
% get ts files
TSfilesFFT={};
files=[];
% for session, get TS files
for sessionNo=1:length(sessions)
    files = dir([mainPath 'Chip/Chip' sessions{sessionNo} '/run*/M*TS.mat']);
    %for files, check if already stabbinned. If not, stabbin
    for fileNo=1:length(files)
        % get stab-bin file
        stabFFTFile=dir([files(fileNo).folder '/*FFT*.mat']);
        % check if exist. If doesn't, add to TSfiles for stabbin
        if isempty(stabFFTFile)
            TSfilesFFT{end+1,1} = fullfile(files(fileNo).folder,files(fileNo).name);
        end
    end
end

% --- Apply multi-frequency FFT --- %
for TSfileNo=1%:length(TSfilesFFT)
    % load TS
    fnTS=TSfilesFFT{TSfileNo};
    load(fnTS);
    Message = autogaStabFFTmultitrials(fnTS,TS,DAVersion);
end

%% Apply ROI TC (WIP)
% autogaStabFFTmulti