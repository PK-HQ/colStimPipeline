function bmp=cloneBitmaps(currentBlockStruct, cloneLoadFlag)
% Check if TS file exists first
if ~isfile(currentBlockStruct.TS)
    % If TS file doesn't exist, return NaN matrix
    bmp = nan(1080, 1920, 2);
    return;
end

% Rest of the original function...
bmp = nan(1080, 1920, 2);
bmpPath = currentBlockStruct.mainPath;
savePath = currentBlockStruct.ROITC;

% Load the file containing the TS variable
loadedData = load(currentBlockStruct.TS);
TS = loadedData.TS;

% Extract the unique .bmp filenames from TS.Header.Conditions.ProjImg
bmpFilenamesRaw = unique(TS.Header.Conditions.ProjImg(cellfun('isempty', strfind(TS.Header.Conditions.ProjImg, 'Dot'))));
bmpFilenamesRaw=bmpFilenamesRaw(end-1:end);
bmpFilenamesRaw=strrep(bmpFilenamesRaw,'\','/');

% Ensure the sessionPath directory exists
if ~exist(savePath, 'dir')
    mkdir(savePath);
end

switch cloneLoadFlag
    case {'clone'}
        resultArray = cellfun(@(x) extractAfter(x, 'PK'), bmpFilenamesRaw, 'UniformOutput', false);
        bmpFilenameTmp = findMatchingFiles([bmpPath 'users/PK/colStimPipeline'],resultArray,currentBlockStruct.date);
        bmpFilenames = {bmpFilenameTmp.filePaths}';
        % Look for these filenames in directoryPath and copy them to sessionPath
        for i = 1:length(bmpFilenames)
            srcFile = bmpFilenames{i}; % was {}
            if exist(srcFile, 'file')
                [~, name, ext] = fileparts(bmpFilenames{i});
                destFile = fullfile(savePath, [name ext]);
                copyfile(srcFile, destFile);
                fprintf('Cloning... Copied %s to %s\n',[name ext],savePath);
                bmp(:,:,i)=double(imread(destFile));
            else
                fprintf('Cloning... File %s not found in %s\n', bmpFilenames{i}, bmpPath);
                bmp(:,:,i)=nan(1080,1920);
            end
        end
    case {'load'}
        % Extracts .bmp filenames (sans path)
        for i=1:length(bmpFilenamesRaw)
            [filepath,filename,fileext]=fileparts(bmpFilenamesRaw{i});

            fileNames{i}=[filename fileext];
        end
        
        bmpStruct = findMatchingFiles(currentBlockStruct.ROITC,fileNames, currentBlockStruct.date);
        bmpFilenames={bmpStruct.filePaths};
        % Look for these filenames in directoryPath and loads
        for i = 1:length(bmpFilenames)
            [~, name, ext] = fileparts(bmpFilenames{i});
            destFile = bmpFilenames{i}; % fullfile(savePath, [name ext])
            if exist(destFile, 'file')
                bmp(:,:,i)=double(imread(destFile));
                fprintf('Bitmap exists %\n',[name ext]);

            else
                bmp(:,:,i)=nan(1080,1920);
                disp('BITMAP MISSING!')
            end
        end
end
end

function results = findMatchingFiles(bmpPath, filePaths, currentDate, varargin)
% FINDMATCHINGFILES  Locate bitmap frames across several roots, with fall-back logic
% and detailed status reporting.
%
%   results = findMatchingFiles(bmpPath, filePaths, currentDate)
%
%   OPTIONAL NAME–VALUE PAIRS
%   -------------------------
%     'ReturnClosest'   (true/false)  – when exact file is missing, choose the
%                        best wildcard match (default = true)
%     'MaxEditDistance' (integer)     – maximum Levenshtein distance allowed
%                        for a “closest” match (default = 5)
%     'LogFile'         (char)        – path to a .txt file; missing files are
%                        appended here (default = '', i.e. no logging)
%
%   OUTPUT
%   ------
%   results  – struct array (one per requested file) with fields:
%              .requestedName   original name from filePaths{i}
%              .filePaths       full path if found; '' otherwise
%              .status          'found' | 'closest' | 'missing' | 'ambiguous'
%              .candidates      cell array of wildcard hits (may be empty)

% -------------------------------------------------------------------------
% Parse inputs
% -------------------------------------------------------------------------
p = inputParser;
addParameter(p,'ReturnClosest',true,@islogical);
addParameter(p,'MaxEditDistance',5,@(x)validateattributes(x,{'numeric'}, ...
    {'scalar','integer','nonnegative'}));
addParameter(p,'LogFile','',@ischar);
parse(p,varargin{:});

opt = p.Results;
writeLog = ~isempty(opt.LogFile);

% -------------------------------------------------------------------------
% Prepare
% -------------------------------------------------------------------------
nRequ = numel(filePaths);
results = repmat(struct('requestedName','','filePaths','','status','', ...
                        'candidates',{{}}), nRequ, 1);

% Minimal anonymous helper for fast suffix stripping
stripSuffix = @(fname) regexprep(fname,'([CS]\d{1,5})$','');  % trims trailing Cxx, Sxxxxx

% Pre-compute edit distance table only when needed
computeED = @(a,b) strEditDistance(a,b);  %#ok<NASGU> (helper below)

% -------------------------------------------------------------------------
% Loop through requested files
% -------------------------------------------------------------------------
for i = 1:nRequ
    reqRelPath          = filePaths{i};
    [~, reqName, reqExt]= fileparts(reqRelPath);
    baseName            = stripSuffix(reqName);          % generalised “prefix”
    filePattern         = [baseName, '*', reqExt];       % wildcard
    results(i).requestedName = [reqName, reqExt];

    % Build the search roots (you can add more here if needed)
    origDir   = fileparts(fullfile(bmpPath, reqRelPath));
    searchDirs= {fullfile('Y:/users/PK/colStimPipeline/ColSeries/', currentDate), fullfile('V:/PK/ColSeries', currentDate), origDir};

    hits      = [];  % cat struct of dir() hits
    for k = 1:numel(searchDirs)
        if isfolder(searchDirs{k})
            tmp = dir(fullfile(searchDirs{k}, filePattern));
            hits = [hits; tmp]; %#ok<AGROW>
        end
    end

    % Filter to .bmp only
    hits = hits(endsWith({hits.name}, '.bmp'));

    % ---------------------------------------------------------------------
    % Decision tree
    % ---------------------------------------------------------------------
    if any(strcmpi({hits.name}, [reqName, reqExt]))       % 1) exact hit
        idx                       = find(strcmpi({hits.name}, [reqName, reqExt]),1);
        results(i).filePaths      = fullfile(hits(idx).folder, hits(idx).name);
        results(i).status         = 'found';
        results(i).candidates     = {hits.name};
        continue
    end

    % 2) No exact match – look for “closest” if allowed
    if opt.ReturnClosest && ~isempty(hits)
        % Compute edit distances vs requested filename
        distances = arrayfun(@(h) strEditDistance(h.name, [reqName, reqExt]), hits);
        [minDist, idx] = min(distances);

        if minDist <= opt.MaxEditDistance                      % Acceptable distance?
            results(i).filePaths  = fullfile(hits(idx).folder, hits(idx).name);
            results(i).status     = 'closest';
        else
            results(i).status     = 'ambiguous';
        end
        results(i).candidates = arrayfun(@(h)h.name, hits, 'uni', false);
    else                                                     % 3) Missing
        results(i).status     = 'missing';
        results(i).candidates = {};
    end

    % Optional logging of miss / ambiguity
    if writeLog && ~strcmp(results(i).status,'found')
        fid = fopen(opt.LogFile,'a');
        if fid>0
            fprintf(fid,'%s\t%s\t%s\n', datestr(now,'yyyy-mm-dd HH:MM:SS'), ...
                    results(i).status, results(i).requestedName);
            fclose(fid);
        end
    end
end
end
% ===== helper ============================================================
function d = strEditDistance(s1,s2)
% Simple Levenshtein distance (case-sensitive). Dynamic programming.
    m = strlength(s1); n = strlength(s2);
    v = 0:n;
    for i = 1:m
        last = i;  v(1) = i;
        for j = 1:n
            new = min([v(j)+1, last+1, v(j)+(s1(i)~=s2(j))]);
            last = v(j+1); v(j+1) = new;
        end
    end
    d = v(end);
end

%{
function matchingFiles = findMatchingFiles(bmpPath, filePaths, currentDate)
    % Initialize output
    matchingFiles = {};

    % Loop over each file path
    for i = 1:length(filePaths)
        % Extract directory and filename from the file path
        [fileDirData, fileName, fileExt] = fileparts(fullfile(bmpPath, filePaths{i}));

        % Initialize a flag to indicate if a match has been found
        foundMatch = false;

        % Define the directories to search
        searchDirs = {fileDirData,['V:/PK/ColSeries/' currentDate], ['Y:/users/PK/colStimPipeline/ColSeries/' currentDate]}; % Search in either the server Y:\ or local experiment computer T:\

        % Iterate over search directories
        for k = 1:length(searchDirs)
            % Create a pattern to match the filename without the potential Cxx suffix
            filePattern = [fileName(1:end-8), '*', fileExt];  % e.g., O09000HE0032G040S00001*

            % Get list of matching files in the current directory
            files = dir(fullfile(searchDirs{k}, filePattern));

            % Filter out non-bmp files
            files = files(endsWith({files.name}, '.bmp'));

            % Check for direct match first
            directMatch = [];
            for j = 1:length(files)
                if strcmp(files(j).name, [fileName, fileExt])
                    directMatch = files(j).name;
                    break;
                end
            end

            % If a direct match is found, add it to the matchingFiles
            if ~isempty(directMatch)
                matchingFiles{end+1} = fullfile(searchDirs{k}, directMatch);
                foundMatch = true;
                break; % Stop searching once a match is found
            end
        end

        % If no direct match was found in any directory, add logic for closest match if needed
        if ~foundMatch
            % Repeat search logic for closest match (if applicable)
            % Similar code block can be added here if no matches are found
        end
    end
end

%}
